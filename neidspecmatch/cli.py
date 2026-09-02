"""Installed command-line entry points for NEIDSpecMatch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import pandas as pd

import neidspecmatch
from .calibration import PopulationRule


def fit_parser():
    parser = argparse.ArgumentParser(
        prog="neidspecmatch-fit",
        description="Fit stellar parameters from a NEID L2 spectrum.",
    )
    parser.add_argument("targetfile", type=Path)
    parser.add_argument("targetname")
    parser.add_argument("output", type=Path)
    parser.add_argument("--orders", type=int, nargs="+", default=[102])
    parser.add_argument("--library-path", type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--fits-path", type=Path)
    parser.add_argument("--library-id")
    parser.add_argument("--dataset-id")
    parser.add_argument(
        "--validation-summary", type=Path,
        help="Current 0.2 cross-validation summary matching the order/library",
    )
    parser.add_argument("--validation-population")
    parser.add_argument(
        "--validation-population-teff", type=float,
        help="Independent Teff used to establish bounded-population membership",
    )
    parser.add_argument("--feh-calibration", type=Path,
                        help="Validated JSON artifact; raw [Fe/H] is the safe default")
    parser.add_argument("--calibration-population")
    parser.add_argument(
        "--calibration-population-teff", type=float,
        help="Independent Teff used for bounded calibration membership",
    )
    parser.add_argument("--max-vsini", type=float, default=30.0)
    parser.add_argument("--rv-source", choices=("drp", "supplied", "custom_ccf"),
                        default="drp")
    parser.add_argument("--absrv", type=float, help="km/s; requires --rv-source supplied")
    parser.add_argument("--ccf-mask", type=Path,
                        help="requires --rv-source custom_ccf")
    parser.add_argument(
        "--ccf-mask-medium",
        choices=("air", "vacuum"),
        help=(
            "input wavelength medium for an unregistered custom CCF mask; "
            "the exact bundled ESPRESSO M3 mask has a SHA-256-bound policy"
        ),
    )
    parser.add_argument("--vsini", type=float)
    parser.add_argument("--vsini-window", type=float, default=0.0)
    parser.add_argument(
        "--input-is-deblazed", action="store_true",
        help="HDU1 already contains corrected flux; bypass normal L2 blaze correction",
    )
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--save-plot-data", action="store_true")
    parser.add_argument("--save-legacy-pickle", action="store_true",
                        help="also save trusted-files-only pickle; JSON is default")
    parser.add_argument("--allow-mixed-drp", action="store_true",
                        help="allow but mark a DRP mismatch as unvalidated")
    parser.add_argument(
        "--allow-unmanifested-library", action="store_true",
        help="allow but mark an absent/invalid library manifest as unvalidated",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-v", "--verbose", action="count", default=0)
    return parser


def fit_main(argv=None):
    parser = fit_parser()
    args = parser.parse_args(argv)
    if args.vsini is None and args.vsini_window != 0.0:
        parser.error("--vsini-window requires --vsini")
    if (args.rv_source == "supplied") != (args.absrv is not None):
        parser.error("--absrv and '--rv-source supplied' must be used together")
    if (args.rv_source == "custom_ccf") != (args.ccf_mask is not None):
        parser.error("--ccf-mask and '--rv-source custom_ccf' must be used together")
    if args.ccf_mask_medium is not None and args.rv_source != "custom_ccf":
        parser.error("--ccf-mask-medium requires '--rv-source custom_ccf'")
    library_path = Path(neidspecmatch.resolve_library_path(args.library_path))
    catalog = args.catalog or library_path / neidspecmatch.DEFAULT_LIBRARY_CATALOG
    fits_path = args.fits_path or library_path / "FITS"
    library_id = args.library_id or library_path.name
    try:
        neidspecmatch.validate_library(
            library_path, require_manifest=True, deep=True,
            expected_library_id=library_id, catalog_name=catalog.name,
            validate_fits_schema=True,
        )
    except (FileNotFoundError, ValueError):
        if not args.allow_unmanifested_library:
            raise
    args.output.mkdir(parents=True, exist_ok=True)
    neidspecmatch.run_specmatch_for_orders(
        targetfile=str(args.targetfile), targetname=args.targetname,
        outputdirectory=str(args.output), path_df_lib=str(catalog),
        path_df_lib_fits=str(fits_path), library_id=library_id,
        orders=args.orders, maxvsini=args.max_vsini,
        calibrate_feh=args.feh_calibration is not None,
        feh_calibration=args.feh_calibration,
        calibration_population=args.calibration_population,
        calibration_population_teff=args.calibration_population_teff,
        deblazed=args.input_is_deblazed, mode="HR", plot=args.plot,
        save_plot_data=args.save_plot_data,
        save_legacy_pickle=args.save_legacy_pickle,
        vsini=args.vsini, vsini_window=args.vsini_window,
        absrv=args.absrv, rv_source=args.rv_source,
        ccf_mask_path=args.ccf_mask, ccf_mask_medium=args.ccf_mask_medium,
        allow_mixed_drp=args.allow_mixed_drp,
        dataset_id=args.dataset_id, verbose=args.verbose,
        random_seed=args.seed,
        validation_summary=args.validation_summary,
        validation_population=args.validation_population,
        validation_population_teff=args.validation_population_teff,
        allow_unmanifested_library=args.allow_unmanifested_library,
    )
    return 0


def parse_population(value):
    parts = value.split(":")
    if len(parts) > 3 or not parts[0]:
        raise argparse.ArgumentTypeError("population must be NAME[:TEFF_MIN[:TEFF_MAX]]")
    parts += [""] * (3 - len(parts))
    try:
        return PopulationRule(
            parts[0], float(parts[1]) if parts[1] else None,
            float(parts[2]) if parts[2] else None,
        )
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def crossval_parser():
    parser = argparse.ArgumentParser(
        prog="neidspecmatch-crossval",
        description="Run spectral LOO and optional [Fe/H] calibration cross-fitting.",
    )
    parser.add_argument("--orders", type=int, nargs="+")
    parser.add_argument("--summarize-existing", type=Path, nargs="+", metavar="CSV")
    parser.add_argument("--library-path", type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--fits-path", type=Path)
    parser.add_argument("--library-id")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fit-feh-calibration", action="store_true")
    parser.add_argument("--population", action="append", type=parse_population)
    parser.add_argument("--min-calibration-stars", type=int, default=8)
    parser.add_argument("--max-vsini", type=float, default=30.0)
    parser.add_argument("--plot-summary", action="store_true")
    parser.add_argument("--plot-each-target", action="store_true")
    parser.add_argument("--allow-mixed-drp", action="store_true")
    parser.add_argument("--allow-unmanifested-library", action="store_true")
    parser.add_argument(
        "--no-resume", action="store_true",
        help="ignore an existing compatible per-fold checkpoint",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-v", "--verbose", action="count", default=0)
    return parser


def _order_from_path(path):
    match = re.search(r"(?:^|[_-])o(\d+)(?=$|[_-])", path.stem)
    if match is None:
        raise ValueError(f"Cannot infer an order from {path}; use a filename containing oNN.")
    return int(match.group(1))


def crossval_main(argv=None):
    parser = crossval_parser()
    args = parser.parse_args(argv)
    if bool(args.orders) == bool(args.summarize_existing):
        parser.error("provide exactly one of --orders or --summarize-existing")
    populations = args.population or [PopulationRule()]
    args.output.mkdir(parents=True, exist_ok=True)
    summary_paths = []
    if args.summarize_existing:
        if args.fit_feh_calibration:
            parser.error(
                "--fit-feh-calibration cannot be used with archived CSVs: "
                "successful fold/checkpoint and strict validation evidence are required"
            )
        if not args.library_id:
            parser.error("--library-id is required when summarizing existing results")
        for filename in args.summarize_existing:
            order = _order_from_path(filename)
            products = neidspecmatch.write_crossvalidation_products(
                pd.read_csv(filename), order=order,
                outputdir=args.output / f"o{order}_crossval",
                library_id=args.library_id, cv_source=str(filename.resolve()),
                fit_feh_calibration=args.fit_feh_calibration,
                population_rules=populations,
                min_calibration_stars=args.min_calibration_stars,
            )
            summary_paths.append(products["summary"])
    else:
        library_path = Path(neidspecmatch.resolve_library_path(args.library_path))
        catalog = args.catalog or library_path / neidspecmatch.DEFAULT_LIBRARY_CATALOG
        fits_path = args.fits_path or library_path / "FITS"
        library_id = args.library_id or library_path.name
        try:
            neidspecmatch.validate_library(
                library_path, require_manifest=True, deep=True,
                expected_library_id=library_id, catalog_name=catalog.name,
                validate_fits_schema=True,
            )
        except (FileNotFoundError, ValueError):
            if not args.allow_unmanifested_library:
                raise
        for order in args.orders:
            frame = neidspecmatch.run_crossvalidation_for_orders(
                order=order, df_lib=catalog, outputdir=args.output,
                path_df_lib_fits=str(fits_path),
                plot_results=args.plot_summary,
                plot_each_target=args.plot_each_target,
                fit_feh_calibration=args.fit_feh_calibration,
                population_rules=populations, library_id=library_id,
                verbose=args.verbose,
                min_calibration_stars=args.min_calibration_stars,
                random_seed=args.seed, allow_mixed_drp=args.allow_mixed_drp,
                resume=not args.no_resume,
                allow_unmanifested_library=args.allow_unmanifested_library,
                maxvsini=args.max_vsini,
            )
            summary_paths.append(frame.attrs["products"]["summary"])
    aggregate = pd.concat(
        [pd.read_csv(filename) for filename in summary_paths], ignore_index=True
    ).sort_values(["order", "population"])
    aggregate_path = args.output / "crossvalidation_summary.csv"
    aggregate.to_csv(aggregate_path, index=False)
    print(aggregate_path)
    return 0


def library_parser():
    parser = argparse.ArgumentParser(
        prog="neidspecmatch-library",
        description="Download or validate the empirical reference library.",
    )
    parser.add_argument("--library-path", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument(
        "--create-manifest", action="store_true",
        help="create/replace a deterministic manifest for an existing library",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--deep", action="store_true")
    parser.add_argument("--fits-schema", action="store_true")
    parser.add_argument("--library-id")
    parser.add_argument(
        "--catalog",
        help="catalog basename (defaults to manifest value when validating)",
    )
    parser.add_argument("--expected-fits", type=int)
    parser.add_argument("--source-url")
    return parser


def library_main(argv=None):
    parser = library_parser()
    args = parser.parse_args(argv)
    if args.download and args.create_manifest:
        parser.error("--download and --create-manifest are mutually exclusive")
    if args.download:
        result = neidspecmatch.get_library(
            overwrite=args.overwrite, library_path=args.library_path
        )
    elif args.create_manifest:
        if not args.library_id:
            parser.error("--create-manifest requires --library-id")
        manifest_path = neidspecmatch.build_library_manifest(
            args.library_path, library_id=args.library_id,
            catalog_name=(args.catalog or neidspecmatch.DEFAULT_LIBRARY_CATALOG),
            expected_fits=args.expected_fits,
            source_url=args.source_url,
        )
        result = neidspecmatch.validate_library(
            args.library_path, expected_fits=args.expected_fits,
            require_manifest=True, deep=True,
            validate_fits_schema=args.fits_schema,
            expected_library_id=args.library_id,
            catalog_name=(args.catalog or neidspecmatch.DEFAULT_LIBRARY_CATALOG),
        )
        result["manifest_path"] = str(manifest_path)
    else:
        result = neidspecmatch.validate_library(
            args.library_path, expected_fits=args.expected_fits,
            require_manifest=True, deep=args.deep,
            validate_fits_schema=args.fits_schema,
            expected_library_id=args.library_id,
            catalog_name=args.catalog,
        )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
