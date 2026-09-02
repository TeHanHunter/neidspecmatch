import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scripts.plot_crossvalidation_research_note as figure


class CrossvalidationFigureTests(unittest.TestCase):
    def _manifest_payload(self, drp_version="v1.5.3"):
        fits_paths = [f"FITS/reference-{index:02d}.fits" for index in range(78)]
        files = [{
            "path": "library.csv",
            "size_bytes": 100,
            "sha256": "a" * 64,
        }]
        files.extend({
            "path": path,
            "size_bytes": 1000 + index,
            "sha256": f"{index + 1:064x}",
        } for index, path in enumerate(fits_paths))
        return {
            "schema_version": 2,
            "library_id": "test-drp15-library",
            "catalog": "library.csv",
            "fits_count": 78,
            "files": files,
            "archive_products": [{
                "l2filename": Path(path).name,
                "swversion": drp_version,
                "flagged": "0",
                "rejected": "0",
            } for path in fits_paths],
            "reference_dq_records": [{
                "path": path,
                "dq_status": "pass",
            } for path in fits_paths],
            "source_archive": {
                "all_unflagged_unrejected": True,
                "required_drp_major_minor": "1.5",
                "swversion_counts": {drp_version: 78},
                "fits_dq_status_counts": {"pass": 78},
            },
        }

    def _write_manifest(self, root, drp_version="v1.5.3"):
        path = Path(root) / "library_manifest.json"
        path.write_text(
            json.dumps(self._manifest_payload(drp_version), sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return path

    def _results(self):
        cool = np.linspace(3300.0, 4400.0, 46)
        hot = np.linspace(4700.0, 6200.0, 32)
        truth = np.concatenate([cool, hot])
        residual = np.linspace(-25.0, 25.0, len(truth))
        feh_true = np.linspace(-0.8, 0.5, len(truth))
        logg_true = np.linspace(4.3, 5.0, len(truth))
        return pd.DataFrame({
            "teff": truth + residual,
            "feh": feh_true + residual / 1000.0,
            "logg": logg_true + residual / 1000.0,
            "vsini": np.full(len(truth), 2.0),
            "d_teff": residual,
            "d_feh": residual / 1000.0,
            "d_logg": residual / 1000.0,
            "teff_true": truth,
            "feh_true": feh_true,
            "logg_true": logg_true,
            "targetname": [f"star-{index}" for index in range(len(truth))],
        })

    def _summary(self):
        return pd.DataFrame({
            "population": ["all", "cool", "hot"],
            "teff_min": [np.nan, np.nan, 4500.0],
            "teff_max": [np.nan, 4500.0, np.nan],
            "n_stars": [78, 46, 32],
            "teff_predictive_rmse": [100.0, 70.0, 130.0],
            "feh_raw_predictive_rmse": [0.12, 0.15, 0.08],
            "logg_predictive_rmse": [0.06, 0.04, 0.09],
            "source_support_limited_folds": [2, 1, 1],
        })

    def _write_bound_order(self, root, order, library, results):
        order_dir = Path(root) / f"o{order}_crossval"
        order_dir.mkdir()
        raw_path = order_dir / f"crossvalidation_results_o{order}.csv"
        checkpoint_path = order_dir / f"crossvalidation_checkpoint_o{order}.json"
        results.to_csv(raw_path, index=False)

        class Reference:
            def __init__(self, name, rv):
                self.object = name
                self.rv = rv
                self.rv_source = "library_catalog"
                self.blaze_source_used = "l2"

        references = [
            Reference(name, float(index))
            for index, name in enumerate(results["targetname"])
        ]
        wavelength = np.linspace(5000.0 + order, 5010.0 + order, 20)
        estimator_config = figure.core._build_estimator_config(
            order=order,
            wavelength=wavelength,
            maxvsini=30.0,
            library_id=library["library_id"],
            library_manifest={"manifest_sha256": library["sha256"]},
            deblazed=False,
            blaze_source="l2",
            reference_rv_state=figure.core._reference_rv_state(references),
            reference_blaze_state=figure.core._reference_blaze_state(references),
        )
        estimator_json = figure.core._canonical_estimator_config_json(
            estimator_config
        )

        support_library = pd.DataFrame({
            "OBJECT_ID": [f"component-{index}" for index in range(8)],
            "Teff": np.arange(8, dtype=float) * 200.0 + 3000.0,
            "[Fe/H]": np.linspace(-0.8, 0.6, 8),
            "log(g)": np.linspace(4.4, 5.1, 8),
            "vsini": np.arange(8, dtype=float),
            "rss": np.arange(8, dtype=float) + 1.0,
        })
        support = figure.core._build_support_diagnostics(
            support_library,
            support_library.iloc[1:6].reset_index(drop=True),
            np.full(5, 0.2),
        )
        fits_paths = list(library["fits_records"])
        identities = []
        for index, (name, relative) in enumerate(
            zip(results["targetname"], fits_paths)
        ):
            identities.append({
                "object_id": name,
                "basename": Path(relative).name,
                "file_sha256": library["fits_records"][relative]["sha256"],
                "loaded_state_sha256": f"{index + 1000:064x}",
                "drp_version": library["version_by_path"][relative],
                "dq_status": library["dq_by_path"][relative],
                "loaded_order_range": [order, order + 1],
            })
        value_columns = [
            "teff", "feh", "logg", "vsini", "d_teff", "d_feh", "d_logg",
        ]
        folds = [{
            "index": index,
            "targetname": results.loc[index, "targetname"],
            "target_identity": identities[index],
            "values": [
                float(value) for value in results.loc[index, value_columns]
            ],
            "optimizer": {
                "success": True,
                "reported_success": True,
                "finite_objective": True,
                "message": "converged",
                "nfev": 1,
                "nit": 1,
                "objective_log_likelihood": -1.0,
                "composite_nonunique_optimum": False,
                "pairwise_all_references_success": True,
                "pairwise_failed_reference_count": 0,
            },
            "support": support,
        } for index in range(len(results))]
        checkpoint = {
            "key": {
                "checkpoint_schema_version": 4,
                "result_schema_version": figure.core.RESULT_SCHEMA_VERSION,
                "neidspecmatch_version": figure.__version__,
                "neidspec_version": getattr(
                    figure.neidspec, "__version__", "unknown"
                ),
                "neidspec_source_sha256": (
                    figure.core._neidspec_source_fingerprint()
                ),
                "pipeline_source_sha256": (
                    figure.core._pipeline_source_fingerprint()
                ),
                "selection_metric": "unweighted_residual_sum_of_squares",
                "order": order,
                "library_id": library["library_id"],
                "library_manifest_sha256": library["sha256"],
                "catalog_sha256": library["catalog_sha256"],
                "drp_versions": library["drp_counts"],
                "blaze_source": "l2",
                "estimator_config_schema_version": (
                    figure.core.ESTIMATOR_CONFIG_SCHEMA_VERSION
                ),
                "estimator_config_sha256": (
                    figure.core._hash_estimator_config(estimator_config)
                ),
                "estimator_config_json": estimator_json,
                "random_seed": 0,
                "spectrum_identities_in_order": identities,
                "reference_dq_status_counts": library["quality_counts"],
            },
            "status": "complete",
            "folds": folds,
            "result_csv": raw_path.name,
            "result_csv_sha256": figure.sha256(raw_path),
            "result_csv_hash_scope": "file_bytes",
            "optimizer": {
                "all_folds_success": True,
                "failed_fold_indices": [],
            },
        }
        checkpoint_path.write_text(
            json.dumps(checkpoint, sort_keys=True) + "\n", encoding="utf-8"
        )
        figure.core.write_crossvalidation_products(
            results,
            order=order,
            outputdir=order_dir,
            library_id=library["library_id"],
            cv_source=str(raw_path),
            population_rules=[
                figure.core.PopulationRule(name="all"),
                figure.core.PopulationRule(
                    name="cool", teff_max=figure.COOL_LIMIT_K
                ),
                figure.core.PopulationRule(
                    name="hot", teff_min=figure.COOL_LIMIT_K
                ),
            ],
            source_pipeline_version=figure.__version__,
            source_drp_compatibility_status="consistent",
            source_library_manifest_sha256=library["sha256"],
            source_neidspec_version=getattr(
                figure.neidspec, "__version__", "unknown"
            ),
            source_blaze_source="l2",
            source_drp_versions=json.dumps(library["drp_counts"], sort_keys=True),
            source_pipeline_sha256=figure.core._pipeline_source_fingerprint(),
            source_neidspec_source_sha256=(
                figure.core._neidspec_source_fingerprint()
            ),
            source_cv_results_sha256=figure.sha256(raw_path),
            source_cv_results_hash_scope="file_bytes",
            source_checkpoint=checkpoint_path.name,
            source_checkpoint_sha256=figure.sha256(checkpoint_path),
            source_checkpoint_hash_scope="file_bytes",
            source_reference_dq_status_counts=json.dumps(
                library["quality_counts"], sort_keys=True
            ),
            source_all_folds_optimizer_success=True,
            source_failed_optimizer_folds=0,
            source_estimator_config=estimator_config,
            source_fold_support=[support for _ in range(len(results))],
        )

    def test_release_constants_and_four_order_layout(self):
        self.assertEqual(figure.ORDERS, (55, 101, 102, 103))
        self.assertEqual(figure.COOL_LIMIT_K, 4500.0)
        data = {order: self._results() for order in figure.ORDERS}
        summaries = {order: self._summary() for order in figure.ORDERS}
        rendered = figure.make_figure(
            data, summaries, drp_label="DRP 1.5.x"
        )
        try:
            self.assertEqual(len(rendered.axes), 12)
            self.assertIn("DRP 1.5.x", rendered._suptitle.get_text())
            labels = [text.get_text() for text in rendered.legends[0].texts]
            self.assertTrue(any("4500" in label for label in labels))
            row_text = "\n".join(
                text.get_text() for axis in rendered.axes for text in axis.texts
            )
            for order in figure.ORDERS:
                self.assertIn(f"Order {order}", row_text)
        finally:
            plt.close(rendered)

    def test_manifest_preflight_accepts_only_complete_drp15_receipt(self):
        with tempfile.TemporaryDirectory() as temp_name:
            path = self._write_manifest(temp_name)
            evidence = figure.load_library_evidence(path)
            self.assertEqual(evidence["drp_counts"], {"v1.5.3": 78})
            self.assertEqual(evidence["quality_counts"], {"pass": 78})

            path = self._write_manifest(temp_name, drp_version="v1.4.0")
            with self.assertRaisesRegex(ValueError, "DRP-1.5"):
                figure.load_library_evidence(path)

    def test_order_loader_fails_closed_when_core_rejects_evidence(self):
        with tempfile.TemporaryDirectory() as temp_name:
            root = Path(temp_name)
            manifest_path = self._write_manifest(root)
            library = figure.load_library_evidence(manifest_path)
            order_dir = root / "o55_crossval"
            order_dir.mkdir()
            results = self._results()
            summary = self._summary()
            results.to_csv(order_dir / "crossvalidation_results_o55.csv", index=False)
            summary.to_csv(order_dir / "crossvalidation_summary_o55.csv", index=False)
            checkpoint = {
                "key": {
                    "catalog_sha256": library["catalog_sha256"],
                    "estimator_config_json": "{}",
                    "random_seed": 0,
                }
            }
            (order_dir / "crossvalidation_checkpoint_o55.json").write_text(
                json.dumps(checkpoint) + "\n", encoding="utf-8"
            )
            with (
                mock.patch.object(
                    figure,
                    "_validate_estimator_configuration",
                    return_value="shared-configuration",
                ),
                mock.patch.object(
                    figure,
                    "_validate_checkpoint_identities",
                    return_value=(("shared",),),
                ),
                mock.patch.object(
                    figure.core,
                    "_validation_product_status",
                    return_value={"matched": False, "reason": "stale_source"},
                ),
            ):
                with self.assertRaisesRegex(
                    ValueError, "evidence rejected: stale_source"
                ):
                    figure.load_order(root, 55, library)

    def test_full_four_order_preflight_accepts_current_bound_evidence(self):
        with tempfile.TemporaryDirectory() as temp_name:
            root = Path(temp_name)
            manifest_path = self._write_manifest(root)
            library = figure.load_library_evidence(manifest_path)
            results = self._results()
            for order in figure.ORDERS:
                self._write_bound_order(root, order, library, results)

            data, summaries, audits, accepted_library = (
                figure.load_validated_bundle(root, manifest_path)
            )
            self.assertEqual(tuple(data), figure.ORDERS)
            self.assertEqual(tuple(summaries), figure.ORDERS)
            self.assertEqual(tuple(audits), figure.ORDERS)
            self.assertEqual(
                accepted_library["sha256"], figure.sha256(manifest_path)
            )
            for order in figure.ORDERS:
                for status in audits[order]["population_status"].values():
                    self.assertTrue(status["matched"], status)

            output_dir = root / "figures"
            argv = [
                "plot_crossvalidation_research_note.py",
                "--validation-dir", str(root),
                "--library-manifest", str(manifest_path),
                "--output-dir", str(output_dir),
            ]
            with (
                mock.patch("sys.argv", argv),
                mock.patch("builtins.print"),
            ):
                figure.main()
            expected_outputs = {
                "crossvalidation_research_note_drp15.png",
                "crossvalidation_research_note_drp15.pdf",
                "crossvalidation_summary_combined.csv",
                "crossvalidation_research_note_drp15.provenance.json",
            }
            self.assertEqual(
                {path.name for path in output_dir.iterdir()}, expected_outputs
            )
            provenance = json.loads(
                (output_dir / (
                    "crossvalidation_research_note_drp15.provenance.json"
                )).read_text(encoding="utf-8")
            )
            self.assertEqual(provenance["orders"], list(figure.ORDERS))
            self.assertEqual(
                provenance["cool_hot_boundary_k"], figure.COOL_LIMIT_K
            )
            self.assertEqual(
                provenance["drp_version_counts"], {"v1.5.3": 78}
            )
            for order in figure.ORDERS:
                self.assertIn(
                    f"o{order}_crossval/crossvalidation_checkpoint_o{order}.json",
                    provenance["inputs_sha256"],
                )
                self.assertEqual(
                    set(provenance["strict_validation_receipts"][str(order)]),
                    {"all", "cool", "hot"},
                )

            raw_path = (
                root / "o102_crossval" / "crossvalidation_results_o102.csv"
            )
            corrupted = pd.read_csv(raw_path)
            corrupted.loc[0, "teff"] += 1.0
            corrupted.to_csv(raw_path, index=False)
            with self.assertRaisesRegex(ValueError, "evidence rejected"):
                figure.load_validated_bundle(root, manifest_path)

            stale_output_dir = root / "stale-figures"
            stale_argv = [
                "plot_crossvalidation_research_note.py",
                "--validation-dir", str(root),
                "--library-manifest", str(manifest_path),
                "--output-dir", str(stale_output_dir),
            ]
            with mock.patch("sys.argv", stale_argv):
                with self.assertRaisesRegex(ValueError, "evidence rejected"):
                    figure.main()
            self.assertFalse(stale_output_dir.exists())


if __name__ == "__main__":
    unittest.main()
