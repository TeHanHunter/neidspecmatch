#!/usr/bin/env python3
"""Create the publication-facing NEIDSpecMatch cross-validation figure.

The layout follows the original research-note figure: each reference star is
shown at its catalog location and connected to the leave-one-out recovered
location.  One row is shown for each of the four release-validation orders.

This script deliberately fails closed.  It accepts only complete raw results,
checkpoint receipts, and summary rows that the installed NEIDSpecMatch code can
reconstruct and validate against the supplied DRP-1.5 library manifest.  It
must not be used as a generic plotter for historical or exploratory CSV files.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import string

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import neidspec
import neidspecmatch.neidspecmatch as core
import neidspecmatch.utils as package_utils
from neidspecmatch.version import __version__


ORDERS = (55, 101, 102, 103)
COOL_LIMIT_K = 4500.0
EXPECTED_DRP_SERIES = "1.5"
EXPECTED_REFERENCE_COUNT = 78
POPULATION_BOUNDS = {
    "all": (None, None),
    "cool": (None, COOL_LIMIT_K),
    "hot": (COOL_LIMIT_K, None),
}
INTENDED_POPULATION = {
    55: "hot",
    101: "cool",
    102: "cool",
    103: "cool",
}
COLORS = {"cool": "#0072B2", "hot": "#D55E00"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _regular_file(path: Path, label: str) -> Path:
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be a regular, non-symlink file: {path}")
    return path.resolve()


def _strict_count_mapping(value: object, label: str) -> dict[str, int]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a mapping")
    output: dict[str, int] = {}
    for key, count in value.items():
        if isinstance(count, bool):
            raise ValueError(f"{label} contains a Boolean count")
        try:
            number = float(count)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{label} contains a non-numeric count") from exc
        if not np.isfinite(number) or not number.is_integer() or number < 0:
            raise ValueError(f"{label} contains an invalid count")
        name = str(key).strip()
        if not name or name in output:
            raise ValueError(f"{label} contains an invalid or duplicate key")
        output[name] = int(number)
    return dict(sorted(output.items()))


def _strict_integer(value: object, label: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be an integer")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be an integer") from exc
    if not np.isfinite(number) or not number.is_integer():
        raise ValueError(f"{label} must be an integer")
    return int(number)


def _safe_manifest_record(record: object, label: str) -> tuple[str, int, str]:
    if not isinstance(record, dict):
        raise ValueError(f"{label} must be a mapping")
    path = str(record.get("path", ""))
    if (
        not path
        or path.startswith("/")
        or "\\" in path
        or any(part in {"", ".", ".."} for part in path.split("/"))
    ):
        raise ValueError(f"{label} has an unsafe path")
    try:
        size = _strict_integer(record["size_bytes"], f"{label} byte size")
    except KeyError as exc:
        raise ValueError(f"{label} has an invalid byte size") from exc
    digest = str(record.get("sha256", ""))
    if size <= 0 or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError(f"{label} has invalid integrity metadata")
    return path, size, digest


def load_library_evidence(manifest_path: Path) -> dict[str, object]:
    """Load the release library receipt and derive its comparison state."""
    manifest_path = _regular_file(manifest_path, "library manifest")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("Library manifest is not readable JSON") from exc
    if not isinstance(manifest, dict):
        raise ValueError("Library manifest must be a JSON object")

    library_id = str(manifest.get("library_id", "")).strip()
    catalog = str(manifest.get("catalog", "")).strip()
    try:
        fits_count = _strict_integer(
            manifest["fits_count"], "library manifest FITS count"
        )
    except KeyError as exc:
        raise ValueError("Library manifest has an invalid FITS count") from exc
    if (
        _strict_integer(
            manifest.get("schema_version", -1), "library manifest schema version"
        )
        != package_utils.LIBRARY_MANIFEST_SCHEMA
        or not library_id
        or not catalog
        or fits_count != EXPECTED_REFERENCE_COUNT
    ):
        raise ValueError(
            "Library manifest does not describe the expected 78-star schema"
        )

    records = manifest.get("files")
    if not isinstance(records, list):
        raise ValueError("Library manifest has no file records")
    file_records: dict[str, dict[str, object]] = {}
    for index, record in enumerate(records):
        path, size, digest = _safe_manifest_record(
            record, f"library file record {index}"
        )
        if path in file_records:
            raise ValueError("Library manifest contains duplicate file paths")
        file_records[path] = {"size_bytes": size, "sha256": digest}
    fits_records = {
        path: record
        for path, record in file_records.items()
        if path.startswith("FITS/") and path.lower().endswith(".fits")
    }
    if (
        len(fits_records) != fits_count
        or set(file_records) != {catalog, *fits_records}
    ):
        raise ValueError(
            "Library manifest file list must contain one catalog and 78 FITS files"
        )
    catalog_sha256 = str(file_records[catalog]["sha256"])

    archive_products = manifest.get("archive_products")
    if not isinstance(archive_products, list) or len(archive_products) != fits_count:
        raise ValueError("Library manifest lacks one archive receipt per FITS file")
    version_by_path: dict[str, str] = {}
    for product in archive_products:
        if not isinstance(product, dict):
            raise ValueError("Library archive receipt must be a mapping")
        filename = str(product.get("l2filename", ""))
        relative = f"FITS/{filename}"
        version = str(product.get("swversion", "")).strip()
        if (
            relative not in fits_records
            or relative in version_by_path
            or product.get("flagged") not in (0, "0", False)
            or product.get("rejected") not in (0, "0", False)
            or core._drp_minor_series(version) != EXPECTED_DRP_SERIES
        ):
            raise ValueError(
                "Archive receipts do not uniquely bind unflagged DRP-1.5 products"
            )
        version_by_path[relative] = version
    if set(version_by_path) != set(fits_records):
        raise ValueError("Archive receipts do not cover the manifest FITS files")
    drp_counts = dict(sorted(
        pd.Series(list(version_by_path.values())).value_counts().astype(int).items()
    ))

    dq_records = manifest.get("reference_dq_records")
    if not isinstance(dq_records, list) or len(dq_records) != fits_count:
        raise ValueError("Library manifest lacks one DQ receipt per FITS file")
    dq_by_path: dict[str, str] = {}
    for record in dq_records:
        if not isinstance(record, dict):
            raise ValueError("Library DQ receipt must be a mapping")
        relative = str(record.get("path", ""))
        status = str(record.get("dq_status", "")).strip().lower()
        if relative not in fits_records or relative in dq_by_path or not status:
            raise ValueError("Library DQ receipts do not uniquely cover the FITS files")
        dq_by_path[relative] = status
    if set(dq_by_path) != set(fits_records):
        raise ValueError("Library DQ receipts do not cover the manifest FITS files")
    quality_counts = dict(sorted(
        pd.Series(list(dq_by_path.values())).value_counts().astype(int).items()
    ))

    source_archive = manifest.get("source_archive")
    if not isinstance(source_archive, dict):
        raise ValueError("Library manifest has no source-archive receipt")
    if (
        source_archive.get("all_unflagged_unrejected") is not True
        or str(source_archive.get("required_drp_major_minor", ""))
        != EXPECTED_DRP_SERIES
        or _strict_count_mapping(
            source_archive.get("swversion_counts"), "source DRP counts"
        ) != drp_counts
        or _strict_count_mapping(
            source_archive.get("fits_dq_status_counts"), "source DQ counts"
        ) != quality_counts
    ):
        raise ValueError("Library source-archive summary is inconsistent")

    return {
        "path": manifest_path,
        "sha256": sha256(manifest_path),
        "library_id": library_id,
        "catalog": catalog,
        "catalog_sha256": catalog_sha256,
        "fits_count": fits_count,
        "fits_records": fits_records,
        "version_by_path": version_by_path,
        "dq_by_path": dq_by_path,
        "drp_counts": drp_counts,
        "quality_counts": quality_counts,
        "drp_series": EXPECTED_DRP_SERIES,
        "drp_label": f"DRP {EXPECTED_DRP_SERIES}.x",
    }


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)


def _read_json_mapping(path: Path, label: str) -> dict[str, object]:
    path = _regular_file(path, label)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not readable JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def _same_optional_bound(observed: object, expected: float | None) -> bool:
    number = pd.to_numeric(pd.Series([observed]), errors="coerce").iloc[0]
    if expected is None:
        return bool(pd.isna(number))
    return bool(np.isfinite(number) and float(number) == float(expected))


def _validate_population_rows(
    summary: pd.DataFrame, results: pd.DataFrame, order: int
) -> None:
    if (
        "population" not in summary
        or len(summary) != len(POPULATION_BOUNDS)
        or summary["population"].duplicated().any()
        or set(summary["population"].astype(str)) != set(POPULATION_BOUNDS)
    ):
        raise ValueError(
            f"Order {order} summary must contain exactly all/cool/hot rows"
        )
    for name, (expected_min, expected_max) in POPULATION_BOUNDS.items():
        row = summary.loc[summary["population"].astype(str) == name].iloc[0]
        if not (
            _same_optional_bound(row.get("teff_min"), expected_min)
            and _same_optional_bound(row.get("teff_max"), expected_max)
        ):
            raise ValueError(
                f"Order {order} {name!r} population does not use the "
                f"adopted {COOL_LIMIT_K:.0f} K boundary"
            )

    truth = pd.to_numeric(results["teff_true"], errors="coerce")
    if (truth == COOL_LIMIT_K).any():
        raise ValueError(
            f"Order {order} has a reference exactly on the inclusive population "
            "boundary; cool/hot membership would overlap"
        )
    expected_counts = {
        "all": len(results),
        "cool": int((truth <= COOL_LIMIT_K).sum()),
        "hot": int((truth >= COOL_LIMIT_K).sum()),
    }
    observed_counts = {
        str(row.population): int(row.n_stars)
        for row in summary[["population", "n_stars"]].itertuples(index=False)
    }
    if observed_counts != expected_counts:
        raise ValueError(
            f"Order {order} population counts do not match teff_true"
        )


def _validate_checkpoint_identities(
    checkpoint: dict[str, object], order: int, library: dict[str, object]
) -> tuple[tuple[object, ...], ...]:
    key = checkpoint.get("key")
    if not isinstance(key, dict):
        raise ValueError(f"Order {order} checkpoint has no key")
    identities = key.get("spectrum_identities_in_order")
    if not isinstance(identities, list) or len(identities) != EXPECTED_REFERENCE_COUNT:
        raise ValueError(f"Order {order} checkpoint has incomplete identities")

    fits_records = library["fits_records"]
    versions = library["version_by_path"]
    dq_statuses = library["dq_by_path"]
    assert isinstance(fits_records, dict)
    assert isinstance(versions, dict)
    assert isinstance(dq_statuses, dict)
    observed_paths: set[str] = set()
    observed_objects: set[str] = set()
    projection: list[tuple[object, ...]] = []
    for identity in identities:
        if not isinstance(identity, dict):
            raise ValueError(f"Order {order} checkpoint identity is not a mapping")
        relative = f"FITS/{identity.get('basename', '')}"
        object_id = str(identity.get("object_id", "")).strip()
        loaded_range = identity.get("loaded_order_range")
        if (
            relative not in fits_records
            or relative in observed_paths
            or not object_id
            or object_id in observed_objects
            or not isinstance(loaded_range, list)
            or len(loaded_range) != 2
        ):
            raise ValueError(
                f"Order {order} identities do not uniquely cover the library"
            )
        start, stop = (
            _strict_integer(
                value, f"order {order} identity loaded-order range"
            )
            for value in loaded_range
        )
        if not start <= order < stop:
            raise ValueError(
                f"Order {order} falls outside a checkpoint loaded-order range"
            )
        record = fits_records[relative]
        assert isinstance(record, dict)
        if (
            str(identity.get("file_sha256")) != str(record["sha256"])
            or str(identity.get("drp_version")) != str(versions[relative])
            or str(identity.get("dq_status")).lower() != str(dq_statuses[relative])
        ):
            raise ValueError(
                f"Order {order} checkpoint identity disagrees with the manifest"
            )
        observed_paths.add(relative)
        observed_objects.add(object_id)
        projection.append((
            object_id,
            relative,
            str(identity["file_sha256"]),
            str(identity["drp_version"]),
            str(identity["dq_status"]).lower(),
        ))
    if observed_paths != set(fits_records):
        raise ValueError(f"Order {order} checkpoint omits library spectra")
    return tuple(projection)


def _validate_estimator_configuration(
    configuration: dict[str, object], order: int, library: dict[str, object]
) -> str:
    spectral_grid = configuration.get("spectral_grid")
    library_state = configuration.get("library")
    blaze = configuration.get("blaze")
    vsini = configuration.get("vsini")
    if not all(isinstance(value, dict) for value in (
        spectral_grid, library_state, blaze, vsini
    )):
        raise ValueError(f"Order {order} estimator configuration is incomplete")
    assert isinstance(spectral_grid, dict)
    assert isinstance(library_state, dict)
    assert isinstance(blaze, dict)
    assert isinstance(vsini, dict)
    reference_blaze = blaze.get("reference_library_state")
    if not isinstance(reference_blaze, dict):
        raise ValueError(f"Order {order} has no reference blaze state")
    records = reference_blaze.get("records")
    if (
        _strict_integer(
            spectral_grid.get("order", -1),
            f"order {order} estimator spectral-grid order",
        ) != order
        or str(library_state.get("library_id")) != str(library["library_id"])
        or str(library_state.get("manifest_sha256")) != str(library["sha256"])
        or str(blaze.get("target_blaze_source")).lower() != "l2"
        or reference_blaze.get("sources") != ["l2"]
        or _strict_integer(
            reference_blaze.get("reference_count", -1),
            f"order {order} estimator reference count",
        )
        != EXPECTED_REFERENCE_COUNT
        or not isinstance(records, list)
        or len(records) != EXPECTED_REFERENCE_COUNT
        or any(
            not isinstance(item, dict)
            or str(item.get("blaze_source", "")).lower() != "l2"
            for item in records
        )
        or vsini.get("mode") != "free"
    ):
        raise ValueError(
            f"Order {order} estimator is not the release DRP-1.5/L2/free-vsini state"
        )
    common = copy.deepcopy(configuration)
    common.pop("spectral_grid", None)
    return json.dumps(common, sort_keys=True, separators=(",", ":"), allow_nan=False)


def load_order(
    validation_dir: Path, order: int, library: dict[str, object]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Load and strictly validate one order's complete CV evidence bundle."""
    order_dir = Path(validation_dir) / f"o{order}_crossval"
    results_path = _regular_file(
        order_dir / f"crossvalidation_results_o{order}.csv",
        f"order-{order} raw results",
    )
    summary_path = _regular_file(
        order_dir / f"crossvalidation_summary_o{order}.csv",
        f"order-{order} summary",
    )
    checkpoint_path = _regular_file(
        order_dir / f"crossvalidation_checkpoint_o{order}.json",
        f"order-{order} checkpoint",
    )
    try:
        results = pd.read_csv(results_path)
        summary = pd.read_csv(summary_path)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise ValueError(f"Order {order} CSV evidence is unreadable") from exc
    required = {
        "teff", "feh", "logg", "vsini", "d_teff", "d_feh", "d_logg",
        "teff_true", "feh_true", "logg_true", "targetname",
    }
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"Order {order} results lack columns: {sorted(missing)}")
    numeric_required = sorted(required - {"targetname"})
    numeric = results[numeric_required].apply(pd.to_numeric, errors="coerce")
    if (
        len(results) != EXPECTED_REFERENCE_COUNT
        or numeric.isna().any().any()
        or not np.isfinite(numeric.to_numpy(dtype=float)).all()
        or results["targetname"].astype(str).duplicated().any()
    ):
        raise ValueError(f"Order {order} is not a complete finite 78-fold result")
    _validate_population_rows(summary, results, order)

    checkpoint = _read_json_mapping(checkpoint_path, f"order-{order} checkpoint")
    key = checkpoint.get("key")
    if not isinstance(key, dict):
        raise ValueError(f"Order {order} checkpoint has no key")
    if str(key.get("catalog_sha256")) != str(library["catalog_sha256"]):
        raise ValueError(f"Order {order} checkpoint catalog disagrees with manifest")
    try:
        configuration = json.loads(str(key["estimator_config_json"]))
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"Order {order} checkpoint estimator configuration is invalid"
        ) from exc
    if not isinstance(configuration, dict):
        raise ValueError(f"Order {order} estimator configuration is not a mapping")
    common_configuration = _validate_estimator_configuration(
        configuration, order, library
    )
    identity_projection = _validate_checkpoint_identities(
        checkpoint, order, library
    )

    population_status: dict[str, dict[str, object]] = {}
    for population, bounds in POPULATION_BOUNDS.items():
        membership_teff = None
        if population == "cool":
            membership_teff = COOL_LIMIT_K - 1.0
        elif population == "hot":
            membership_teff = COOL_LIMIT_K + 1.0
        status = core._validation_product_status(
            summary_path,
            order=order,
            library_id=library["library_id"],
            library_manifest={"manifest_sha256": library["sha256"]},
            drp_compatibility={
                "publication_validated": True,
                "versions": library["drp_counts"],
            },
            blaze_source="l2",
            validation_population=(None if population == "all" else population),
            validation_population_teff=membership_teff,
            calibration_metadata={"status": "not_applied"},
            reference_quality_counts=library["quality_counts"],
            reference_drp_versions=library["drp_counts"],
            estimator_config=configuration,
        )
        if not status.get("matched"):
            raise ValueError(
                f"Order {order} {population} evidence rejected: "
                f"{status.get('reason', 'unknown_reason')}"
            )
        population_details = status.get("population")
        if not isinstance(population_details, dict) or (
            str(population_details.get("name")) != population
            or not _same_optional_bound(population_details.get("teff_min"), bounds[0])
            or not _same_optional_bound(population_details.get("teff_max"), bounds[1])
        ):
            raise ValueError(
                f"Order {order} {population} validation selected the wrong row"
            )
        population_status[population] = status

    try:
        random_seed = _strict_integer(
            key["random_seed"], f"order {order} checkpoint seed"
        )
    except KeyError as exc:
        raise ValueError(f"Order {order} checkpoint seed is invalid") from exc
    return results, summary, {
        "results_path": results_path,
        "summary_path": summary_path,
        "checkpoint_path": checkpoint_path,
        "configuration": configuration,
        "common_configuration": common_configuration,
        "identity_projection": identity_projection,
        "random_seed": random_seed,
        "population_status": population_status,
    }


def load_validated_bundle(
    validation_dir: Path, manifest_path: Path
) -> tuple[
    dict[int, pd.DataFrame],
    dict[int, pd.DataFrame],
    dict[int, dict[str, object]],
    dict[str, object],
]:
    """Validate all four orders and their shared analysis state before plotting."""
    validation_dir = Path(validation_dir)
    if validation_dir.is_symlink() or not validation_dir.is_dir():
        raise ValueError(
            f"Validation directory must be a real directory: {validation_dir}"
        )
    library = load_library_evidence(manifest_path)
    data: dict[int, pd.DataFrame] = {}
    summaries: dict[int, pd.DataFrame] = {}
    audits: dict[int, dict[str, object]] = {}
    for order in ORDERS:
        data[order], summaries[order], audits[order] = load_order(
            validation_dir, order, library
        )

    first = ORDERS[0]
    shared_configuration = audits[first]["common_configuration"]
    shared_identities = audits[first]["identity_projection"]
    shared_seed = audits[first]["random_seed"]
    truth_columns = ["targetname", "teff_true", "feh_true", "logg_true"]
    baseline_truth = data[first][truth_columns].reset_index(drop=True)
    for order in ORDERS[1:]:
        if (
            audits[order]["common_configuration"] != shared_configuration
            or audits[order]["identity_projection"] != shared_identities
            or audits[order]["random_seed"] != shared_seed
        ):
            raise ValueError(
                f"Order {order} was not produced from the same library/configuration/seed"
            )
        comparison = data[order][truth_columns].reset_index(drop=True)
        if not comparison["targetname"].astype(str).equals(
            baseline_truth["targetname"].astype(str)
        ) or not np.allclose(
            comparison[["teff_true", "feh_true", "logg_true"]].to_numpy(float),
            baseline_truth[["teff_true", "feh_true", "logg_true"]].to_numpy(float),
            rtol=2e-12,
            atol=2e-12,
        ):
            raise ValueError(f"Order {order} uses different reference truth values")
    return data, summaries, audits, library


def add_displacements(
    ax: plt.Axes,
    data: pd.DataFrame,
    y_true: str,
    y_fit: str,
) -> None:
    populations = np.where(data["teff_true"] < COOL_LIMIT_K, "cool", "hot")
    for population in ("cool", "hot"):
        subset = data.loc[populations == population]
        color = COLORS[population]
        for row in subset.itertuples(index=False):
            ax.plot(
                [row.teff_true, row.teff],
                [getattr(row, y_true), getattr(row, y_fit)],
                color=color,
                alpha=0.30,
                linewidth=0.65,
                zorder=1,
            )
        ax.scatter(
            subset["teff_true"], subset[y_true], s=22, facecolor=color,
            edgecolor="white", linewidth=0.35, alpha=0.92, zorder=3,
        )
        ax.scatter(
            subset["teff"], subset[y_fit], s=19, marker="x", color=color,
            linewidth=0.85, alpha=0.85, zorder=4,
        )


def metrics_text(summary: pd.DataFrame, population: str) -> str:
    row = summary.loc[summary["population"] == population].iloc[0]
    label = "all" if population == "all" else f"{population}, intended"
    return (
        f"{label} (N={int(row.n_stars)}): "
        rf"RMSE $T_{{\rm eff}}$={row.teff_predictive_rmse:.0f} K, "
        rf"[Fe/H]={row.feh_raw_predictive_rmse:.3f} dex, "
        rf"$\log g$={row.logg_predictive_rmse:.3f} dex; "
        f"support-limited={int(row.source_support_limited_folds)}"
    )


def make_figure(
    order_data: dict[int, pd.DataFrame],
    order_summaries: dict[int, pd.DataFrame],
    *,
    drp_label: str,
) -> plt.Figure:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    n_rows = len(ORDERS)
    fig, axes = plt.subplots(
        n_rows,
        3,
        figsize=(11.2, 3.15 * n_rows + 0.7),
        constrained_layout=False,
    )

    all_teff = np.concatenate([order_data[o][["teff_true", "teff"]].to_numpy().ravel() for o in ORDERS])
    teff_limits = (all_teff.min() - 90, all_teff.max() + 90)

    for row_index, order in enumerate(ORDERS):
        data = order_data[order]
        summary = order_summaries[order]
        ax_feh, ax_logg, ax_resid = axes[row_index]

        add_displacements(ax_feh, data, "feh_true", "feh")
        add_displacements(ax_logg, data, "logg_true", "logg")

        populations = np.where(data["teff_true"] < COOL_LIMIT_K, "cool", "hot")
        for population in ("cool", "hot"):
            subset = data.loc[populations == population]
            ax_resid.scatter(
                subset["feh_true"], subset["d_feh"], s=24,
                facecolor=COLORS[population], edgecolor="white", linewidth=0.35,
                alpha=0.90, zorder=3,
            )
        fit = np.polyfit(data["feh_true"], data["d_feh"], deg=1)
        x_fit = np.linspace(data["feh_true"].min(), data["feh_true"].max(), 100)
        ax_resid.axhline(0, color="0.35", linewidth=0.8, linestyle="--", zorder=1)
        ax_resid.plot(x_fit, np.polyval(fit, x_fit), color="black", linewidth=1.2, zorder=2)
        ax_resid.text(
            0.04, 0.87,
            rf"raw residual fit: $\Delta$[Fe/H] = {fit[0]:+.3f}[Fe/H] {fit[1]:+.3f}",
            transform=ax_resid.transAxes, va="top", fontsize=8.2,
        )

        ax_feh.set_ylabel("[Fe/H] (dex)")
        ax_logg.set_ylabel(r"$\log g$ (dex)")
        ax_resid.set_ylabel(r"Recovered $-$ catalog [Fe/H] (dex)")
        ax_feh.set_xlim(teff_limits)
        ax_logg.set_xlim(teff_limits)
        ax_feh.axvline(COOL_LIMIT_K, color="0.65", linestyle=":", linewidth=0.8)
        ax_logg.axvline(COOL_LIMIT_K, color="0.65", linestyle=":", linewidth=0.8)

        intended = INTENDED_POPULATION[order]
        row_text = (
            f"Order {order}  |  {metrics_text(summary, 'all')}\n"
            f"{metrics_text(summary, intended)}"
        )
        ax_feh.text(
            0.0,
            1.18,
            row_text,
            transform=ax_feh.transAxes,
            ha="left",
            va="top",
            fontsize=9.0,
        )

        for ax in axes[row_index]:
            style_axis(ax)

    for row in axes:
        for ax in row[:2]:
            ax.set_xlabel(r"$T_{\rm eff}$ (K)")
        row[2].set_xlabel("Catalog [Fe/H] (dex)")

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["cool"],
               markeredgecolor="white", markersize=6,
               label=rf"Cool ($T_{{\rm eff}}<{COOL_LIMIT_K:.0f}$ K)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["hot"],
               markeredgecolor="white", markersize=6,
               label=rf"Hot ($T_{{\rm eff}}\geq{COOL_LIMIT_K:.0f}$ K)"),
        Line2D([0], [0], marker="o", color="0.25", markerfacecolor="0.25",
               markersize=5, linestyle="none", label="Catalog"),
        Line2D([0], [0], marker="x", color="0.25", markersize=5,
               linestyle="none", label="Recovered"),
    ]
    fig.legend(handles=legend, loc="upper right", bbox_to_anchor=(0.982, 0.995),
               ncol=4, frameon=False, fontsize=8.6, handletextpad=0.4, columnspacing=1.2)

    for label, ax in zip(string.ascii_lowercase, axes.ravel()):
        ax.text(0.015, 0.97, f"({label})", transform=ax.transAxes,
                ha="left", va="top", fontweight="bold", fontsize=9.5)

    fig.suptitle(
        f"NEIDSpecMatch leave-one-out cross-validation — {drp_label} reference library",
        x=0.055, y=0.995, ha="left", fontsize=11.5,
    )
    fig.subplots_adjust(
        left=0.075,
        right=0.982,
        bottom=0.055,
        top=0.93,
        hspace=0.62,
        wspace=0.28,
    )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, required=True)
    parser.add_argument(
        "--library-manifest",
        type=Path,
        required=True,
        help="library_manifest.json used by every plotted cross-validation run",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    validation_dir = args.validation_dir.resolve()
    output_dir = (args.output_dir or validation_dir / "figures").resolve()

    # Complete the entire provenance preflight before creating or overwriting
    # any release-facing output.
    data, summaries, audits, library = load_validated_bundle(
        validation_dir, args.library_manifest
    )
    figure = make_figure(data, summaries, drp_label=str(library["drp_label"]))

    output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = output_dir / "crossvalidation_summary_combined.csv"
    pd.concat(
        [summaries[order] for order in ORDERS], ignore_index=True
    ).to_csv(combined_path, index=False)

    png_path = output_dir / "crossvalidation_research_note_drp15.png"
    pdf_path = output_dir / "crossvalidation_research_note_drp15.pdf"
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    inputs = {
        "library_manifest.json": str(library["sha256"]),
    }
    validation_receipts: dict[str, object] = {}
    for order in ORDERS:
        audit = audits[order]
        order_dir = f"o{order}_crossval"
        for name in ("results_path", "summary_path", "checkpoint_path"):
            path = audit[name]
            assert isinstance(path, Path)
            inputs[f"{order_dir}/{path.name}"] = sha256(path)
        statuses = audit["population_status"]
        assert isinstance(statuses, dict)
        validation_receipts[str(order)] = {
            name: {
                "reason": status["reason"],
                "evidence_validation": status["evidence_validation"],
            }
            for name, status in statuses.items()
        }
    outputs = {
        path.name: sha256(path) for path in (png_path, pdf_path, combined_path)
    }
    provenance = {
        "figure_description": (
            "Strictly reconstructed raw leave-one-out validation; no [Fe/H] "
            "detrending applied. RMSE is a predictive error metric, not an "
            "automatic Gaussian one-sigma uncertainty."
        ),
        "validation_status": "matched_current_raw_checkpoint_summary_evidence",
        "library_id": library["library_id"],
        "library_manifest_sha256": library["sha256"],
        "reference_count": library["fits_count"],
        "drp_series": library["drp_series"],
        "drp_version_counts": library["drp_counts"],
        "reference_dq_status_counts": library["quality_counts"],
        "cool_hot_boundary_k": COOL_LIMIT_K,
        "population_rules": {
            name: {"teff_min": bounds[0], "teff_max": bounds[1]}
            for name, bounds in POPULATION_BOUNDS.items()
        },
        "references_exactly_on_boundary": 0,
        "orders": list(ORDERS),
        "intended_population_by_order": {
            str(order): population
            for order, population in INTENDED_POPULATION.items()
        },
        "neidspecmatch_version": __version__,
        "neidspec_version": getattr(neidspec, "__version__", "unknown"),
        "neidspecmatch_source_sha256": core._pipeline_source_fingerprint(),
        "neidspec_source_sha256": core._neidspec_source_fingerprint(),
        "plot_script_sha256": sha256(Path(__file__).resolve()),
        "strict_validation_receipts": validation_receipts,
        "inputs_sha256": inputs,
        "outputs_sha256": outputs,
    }
    provenance_path = output_dir / "crossvalidation_research_note_drp15.provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(png_path)
    print(pdf_path)
    print(combined_path)
    print(provenance_path)


if __name__ == "__main__":
    main()
