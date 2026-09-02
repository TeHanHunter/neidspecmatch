import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

import neidspec
import neidspecmatch.neidspecmatch as core
from neidspecmatch.version import __version__


class ValidationProvenanceTests(unittest.TestCase):
    def _crossvalidation_frame(self):
        index = np.arange(8, dtype=float)
        return pd.DataFrame({
            "teff": 3400.0 + 40.0 * index,
            "feh": -0.3 + 0.05 * index,
            "logg": 4.7 + 0.01 * index,
            "vsini": np.full(8, 2.0),
            "d_teff": np.linspace(-20.0, 20.0, 8),
            "d_feh": np.linspace(-0.04, 0.04, 8),
            "d_logg": np.linspace(-0.03, 0.03, 8),
            "teff_true": 3400.0 + 40.0 * index - np.linspace(-20.0, 20.0, 8),
            "feh_true": -0.3 + 0.05 * index - np.linspace(-0.04, 0.04, 8),
            "logg_true": 4.7 + 0.01 * index - np.linspace(-0.03, 0.03, 8),
            "targetname": [f"star-{value}" for value in range(8)],
        })

    def _bound_validation_product(
            self, root, *, vsini=None, vsini_window=None):
        """Write internally consistent raw/checkpoint/summary test evidence."""
        root = Path(root)
        raw_path = root / "crossvalidation_results_o102.csv"
        checkpoint_path = root / "crossvalidation_checkpoint_o102.json"
        frame = self._crossvalidation_frame()
        frame.to_csv(raw_path, index=False)
        manifest_sha = "a" * 64
        quality_counts = {"pass": len(frame)}
        wavelength = np.linspace(8500.0, 8510.0, 20)

        class Reference:
            def __init__(self, name, rv):
                self.object = name
                self.rv = rv
                self.rv_source = "library_catalog"
                self.blaze_source_used = "l2"

        references = [
            Reference(name, float(index))
            for index, name in enumerate(frame["targetname"])
        ]
        reference_rv_state = core._reference_rv_state(references)
        reference_blaze_state = core._reference_blaze_state(references)
        estimator_config = core._build_estimator_config(
            order=102, wavelength=wavelength, maxvsini=30.0,
            vsini=vsini, vsini_window=vsini_window,
            library_id="library-v1",
            library_manifest={"manifest_sha256": manifest_sha},
            deblazed=False, blaze_source="l2",
            reference_rv_state=reference_rv_state,
            reference_blaze_state=reference_blaze_state,
        )
        support_library = pd.DataFrame({
            "OBJECT_ID": [f"reference-{index}" for index in range(8)],
            "Teff": np.arange(8, dtype=float) * 200.0 + 3000.0,
            "[Fe/H]": np.linspace(-0.8, 0.6, 8),
            "log(g)": np.linspace(4.4, 5.1, 8),
            "vsini": np.arange(8, dtype=float),
            "rss": np.arange(8, dtype=float) + 1.0,
        })
        support = core._build_support_diagnostics(
            support_library,
            support_library.iloc[1:6].reset_index(drop=True),
            np.full(5, 0.2),
        )
        identities = [{
            "object_id": name,
            "basename": f"{name}.fits",
            "file_sha256": f"{index:064x}",
            "loaded_state_sha256": f"{index + 100:064x}",
            "drp_version": "1.3.0",
            "dq_status": "pass",
            "loaded_order_range": [102, 103],
        } for index, name in enumerate(frame["targetname"])]
        value_columns = [
            "teff", "feh", "logg", "vsini", "d_teff", "d_feh", "d_logg",
        ]
        folds = [{
            "index": index,
            "targetname": frame.loc[index, "targetname"],
            "target_identity": identities[index],
            "values": [float(value) for value in frame.loc[index, value_columns]],
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
        } for index in range(len(frame))]
        estimator_json = core._canonical_estimator_config_json(estimator_config)
        checkpoint = {
            "key": {
                "checkpoint_schema_version": 4,
                "result_schema_version": core.RESULT_SCHEMA_VERSION,
                "neidspecmatch_version": __version__,
                "neidspec_version": getattr(neidspec, "__version__", "unknown"),
                "neidspec_source_sha256": core._neidspec_source_fingerprint(),
                "pipeline_source_sha256": core._pipeline_source_fingerprint(),
                "selection_metric": "unweighted_residual_sum_of_squares",
                "order": 102,
                "library_id": "library-v1",
                "library_manifest_sha256": manifest_sha,
                "catalog_sha256": "b" * 64,
                "drp_versions": {"1.3.0": len(frame)},
                "blaze_source": "l2",
                "estimator_config_schema_version": (
                    core.ESTIMATOR_CONFIG_SCHEMA_VERSION
                ),
                "estimator_config_sha256": core._hash_estimator_config(
                    estimator_config
                ),
                "estimator_config_json": estimator_json,
                "random_seed": 0,
                "spectrum_identities_in_order": identities,
                "reference_dq_status_counts": quality_counts,
            },
            "status": "complete",
            "folds": folds,
            "result_csv": raw_path.name,
            "result_csv_sha256": core.sha256_file(raw_path),
            "result_csv_hash_scope": "file_bytes",
            "optimizer": {
                "all_folds_success": True,
                "failed_fold_indices": [],
            },
        }
        checkpoint_path.write_text(
            json.dumps(checkpoint, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        products = core.write_crossvalidation_products(
            frame, order=102, outputdir=root,
            library_id="library-v1", cv_source=str(raw_path),
            source_pipeline_version=__version__,
            source_drp_compatibility_status="consistent",
            source_library_manifest_sha256=manifest_sha,
            source_neidspec_version=getattr(neidspec, "__version__", "unknown"),
            source_blaze_source="l2",
            source_drp_versions=json.dumps({"1.3.0": len(frame)}),
            source_pipeline_sha256=core._pipeline_source_fingerprint(),
            source_neidspec_source_sha256=core._neidspec_source_fingerprint(),
            source_cv_results_sha256=core.sha256_file(raw_path),
            source_cv_results_hash_scope="file_bytes",
            source_checkpoint=checkpoint_path.name,
            source_checkpoint_sha256=core.sha256_file(checkpoint_path),
            source_checkpoint_hash_scope="file_bytes",
            source_reference_dq_status_counts=json.dumps(
                quality_counts, sort_keys=True
            ),
            source_all_folds_optimizer_success=True,
            source_failed_optimizer_folds=0,
            source_estimator_config=estimator_config,
            source_fold_support=[support for _ in range(len(frame))],
        )
        context = {
            "order": 102,
            "library_id": "library-v1",
            "library_manifest": {"manifest_sha256": manifest_sha},
            "drp_compatibility": {"versions": {"1.3.0": len(frame) + 1}},
            "blaze_source": "l2",
            "reference_quality_counts": quality_counts,
            "calibration_metadata": {"status": "not_applied"},
            "estimator_config": estimator_config,
        }
        return Path(products["summary"]), checkpoint_path, context

    def test_validation_gate_rejects_nonfree_and_malformed_vsini_modes(self):
        reason = (
            "fixed_or_bounded_vsini_requires_independent_"
            "broadlined_validation"
        )
        for mode, value, window in (
                ("fixed", 30.0, 0.0),
                ("bounded", 20.0, 3.0)):
            with (
                    self.subTest(mode=mode),
                    tempfile.TemporaryDirectory() as temp_name):
                summary_path, _, context = self._bound_validation_product(
                    temp_name, vsini=value, vsini_window=window,
                )
                status = core._validation_product_status(
                    summary_path, **context
                )
                self.assertFalse(status["matched"], status)
                self.assertEqual(status["reason"], reason)
                self.assertEqual(
                    status["summary_sha256"], core.sha256_file(summary_path)
                )

        with tempfile.TemporaryDirectory() as temp_name:
            summary_path, _, context = self._bound_validation_product(temp_name)
            free_status = core._validation_product_status(
                summary_path, **context
            )
            self.assertTrue(free_status["matched"], free_status)
            for malformed_vsini in (None, {}, {"mode": "future-mode"}):
                with self.subTest(malformed_vsini=malformed_vsini):
                    malformed_context = copy.deepcopy(context)
                    if malformed_vsini is None:
                        malformed_context["estimator_config"].pop("vsini")
                    else:
                        malformed_context["estimator_config"]["vsini"] = (
                            malformed_vsini
                        )
                    status = core._validation_product_status(
                        summary_path, **malformed_context
                    )
                    self.assertFalse(status["matched"], status)
                    self.assertEqual(
                        status["reason"],
                        "runtime_estimator_vsini_mode_missing_or_unrecognized",
                    )
                    self.assertEqual(
                        status["summary_sha256"],
                        core.sha256_file(summary_path),
                    )

    def test_validation_rejects_incomplete_checkpoint_even_with_matching_hash(self):
        with tempfile.TemporaryDirectory() as temp_name:
            summary_path, checkpoint_path, context = (
                self._bound_validation_product(temp_name)
            )
            checkpoint_path.write_text(
                json.dumps({"status": "complete"}) + "\n", encoding="utf-8"
            )
            summary = pd.read_csv(summary_path)
            summary.loc[0, "source_checkpoint_sha256"] = core.sha256_file(
                checkpoint_path
            )
            summary.to_csv(summary_path, index=False)
            status = core._validation_product_status(summary_path, **context)
            self.assertFalse(status["matched"], status)

            summary_path, checkpoint_path, context = (
                self._bound_validation_product(temp_name)
            )
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            checkpoint["folds"][0] = None
            checkpoint_path.write_text(
                json.dumps(checkpoint) + "\n", encoding="utf-8"
            )
            summary = pd.read_csv(summary_path)
            summary.loc[0, "source_checkpoint_sha256"] = core.sha256_file(
                checkpoint_path
            )
            summary.to_csv(summary_path, index=False)
            malformed = core._validation_product_status(
                summary_path, **context
            )
            self.assertFalse(malformed["matched"], malformed)

    def test_validation_rejects_summary_metric_not_derived_from_raw_sidecar(self):
        with tempfile.TemporaryDirectory() as temp_name:
            summary_path, _, context = self._bound_validation_product(temp_name)
            summary = pd.read_csv(summary_path)
            summary.loc[0, "teff_predictive_rmse"] *= 100.0
            summary.to_csv(summary_path, index=False)
            status = core._validation_product_status(summary_path, **context)
            self.assertFalse(status["matched"], status)

    def test_malformed_support_diagnostics_cannot_validate_parameters(self):
        class Spectrum:
            object = "target"
            drp_version = "1.3.0"
            instrument = "NEID"
            observation_type = "SCI"
            data_level = 2
            observing_mode = "HR"
            dq_status = "pass"
            dq_assessments = {}
            provenance = {"blaze_source": "l2"}

        validation_status = {
            "matched": True,
            "reason": "matched_current_product",
            "crossvalidation_error_metrics": {
                "residual_bias": {
                    "teff_k": 0.0, "feh_raw_dex": 0.0,
                    "feh_calibrated_dex": None, "logg_dex": 0.0,
                },
                "sample_scatter": {
                    "teff_k": 50.0, "feh_raw_dex": 0.1,
                    "feh_calibrated_dex": None, "logg_dex": 0.1,
                },
                "predictive_rmse": {
                    "teff_k": 50.0, "feh_raw_dex": 0.1,
                    "feh_calibrated_dex": None, "logg_dex": 0.1,
                },
            },
        }
        with mock.patch.object(
                core, "_validation_product_status",
                return_value=validation_status):
            provenance = core.build_result_provenance(
                target=Spectrum(), references=[],
                rv_metadata={
                    "target": {
                        "object": "target", "rv_kms": 0.0,
                        "rv_source": "supplied",
                    },
                    "references": [],
                },
                order=102, wavelength=np.linspace(8500.0, 8510.0, 20),
                deblazed=False, vsini=None, vsini_window=None, random_seed=0,
                feh_calibration_metadata={"status": "not_applied"},
                drp_compatibility={
                    "publication_validated": True,
                    "versions": {"1.3.0": 1},
                },
                library_id="library-v1",
                library_manifest={
                    "manifest_present": True, "recognized": True,
                    "deep_verified": True, "reference_pool_verified": True,
                    "library_id": "library-v1",
                },
                composite_optimizer_diagnostics={
                    "success": True, "nonunique_optimum": False,
                },
                pairwise_optimizer_diagnostics={
                    "all_references_success": True,
                },
                estimator_config={"schema_version": 1},
                support_diagnostics={"parameters": {}, "warnings": []},
            )
        for parameter in ("teff", "feh", "logg"):
            self.assertFalse(
                provenance["parameter_validation"][parameter][
                    "publication_validated"
                ],
                (parameter, provenance["overall_validation"]),
            )

    def test_matching_product_binds_sidecars_dq_and_drp_set(self):
        with tempfile.TemporaryDirectory() as temp_name:
            root = Path(temp_name)
            raw_path = root / "crossvalidation_results_o102.csv"
            summary_path, _, context = self._bound_validation_product(root)
            products = {"summary": str(summary_path)}
            manifest_sha = context["library_manifest"]["manifest_sha256"]
            quality_counts = context["reference_quality_counts"]
            estimator_config = context["estimator_config"]
            runtime_drp = context["drp_compatibility"]
            status = core._validation_product_status(
                products["summary"], order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertTrue(status["matched"], status)
            patch_compatible = core._validation_product_status(
                products["summary"], order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility={
                    "versions": {"1.3.0": 8, "1.3.1": 1},
                },
                reference_drp_versions={"1.3.0": 8},
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertTrue(patch_compatible["matched"], patch_compatible)
            self.assertNotIn("uncertainties_1sigma", status)
            metrics = status["crossvalidation_error_metrics"]
            self.assertEqual(metrics["schema_version"], 2)
            self.assertEqual(
                metrics["sample_scatter"]["statistic"],
                "sample_standard_deviation",
            )
            self.assertEqual(metrics["sample_scatter"]["ddof"], 1)
            self.assertEqual(
                metrics["predictive_rmse"]["statistic"],
                "root_mean_square_error",
            )
            self.assertIn(
                "none_is_automatically_a_gaussian_1sigma_uncertainty",
                metrics["uncertainty_semantics"],
            )
            summary_frame = pd.read_csv(products["summary"])
            self.assertEqual(
                summary_frame.loc[0, "uncertainty_statistic"],
                "not_automatically_defined",
            )
            self.assertTrue(pd.isna(
                summary_frame.loc[0, "uncertainty_ddof"]
            ))
            self.assertAlmostEqual(
                summary_frame.loc[0, "teff_sample_scatter"],
                summary_frame.loc[0, "teff_sigma"],
            )
            self.assertAlmostEqual(
                summary_frame.loc[0, "teff_predictive_rmse"],
                summary_frame.loc[0, "teff_rmse"],
            )

            summary_path = Path(products["summary"])
            original_summary = summary_path.read_bytes()
            wrong_definition = pd.read_csv(summary_path)
            wrong_definition.loc[0, "residual_scatter_ddof"] = 0
            wrong_definition.to_csv(summary_path, index=False)
            rejected_definition = core._validation_product_status(
                summary_path, order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(rejected_definition["matched"])
            self.assertIn(
                "error_metric_definition", rejected_definition["reason"]
            )
            summary_path.write_bytes(original_summary)

            mismatched = core._validation_product_status(
                products["summary"], order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts={"pass": 78},
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(mismatched["matched"])
            self.assertIn("dq_composition", mismatched["reason"])

            incomplete_frame = pd.read_csv(summary_path)
            incomplete_frame.loc[0, "teff_n_valid"] = 2
            incomplete_frame.loc[0, "teff_finite_coverage"] = 0.25
            incomplete_frame.to_csv(summary_path, index=False)
            incomplete = core._validation_product_status(
                summary_path, order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(incomplete["matched"])
            self.assertIn("teff_residuals_incomplete", incomplete["reason"])
            summary_path.write_bytes(original_summary)

            escaping_frame = pd.read_csv(summary_path)
            escaping_frame.loc[0, "cv_source"] = "../crossvalidation_results_o102.csv"
            escaping_frame.to_csv(summary_path, index=False)
            escaping = core._validation_product_status(
                summary_path, order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(escaping["matched"])
            self.assertIn("sidecar_hash_mismatch", escaping["reason"])
            summary_path.write_bytes(original_summary)

            mutation_paths = (
                ("spectral_grid", "order"),
                ("spectral_grid", "wavelength_float64_le_sha256"),
                ("library", "manifest_sha256"),
                ("rv", "application_policy"),
                ("rv", "reference_library_state"),
                ("blaze", "target_blaze_source"),
                ("blaze", "reference_library_state"),
                ("mask", "policy"),
                ("pairwise", "top_k"),
                ("pairwise", "continuum_degree"),
                ("pairwise", "limb_darkening_u1"),
                ("pairwise", "selection_metric"),
                ("composite", "loss"),
                ("vsini", "mode"),
                ("vsini", "max_kms"),
                ("vsini", "window_kms"),
            )
            for section, field in mutation_paths:
                changed = copy.deepcopy(estimator_config)
                current = changed[section][field]
                changed[section][field] = (
                    current + 1 if isinstance(current, (int, float))
                    else f"{current}-changed"
                )
                rejected = core._validation_product_status(
                    summary_path, order=102, library_id="library-v1",
                    library_manifest={"manifest_sha256": manifest_sha},
                    drp_compatibility=runtime_drp,
                    blaze_source="l2",
                    reference_quality_counts=quality_counts,
                    calibration_metadata={"status": "not_applied"},
                    estimator_config=changed,
                )
                self.assertFalse(
                    rejected["matched"], (section, field, rejected)
                )

            fixed_config = core._build_estimator_config(
                order=102, wavelength=np.linspace(8500.0, 8510.0, 20),
                maxvsini=30.0, vsini=30.0, vsini_window=0.0,
                library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                deblazed=False, blaze_source="l2",
                reference_rv_state=estimator_config["rv"][
                    "reference_library_state"
                ],
                reference_blaze_state=estimator_config["blaze"][
                    "reference_library_state"
                ],
            )
            fixed_rejected = core._validation_product_status(
                summary_path, order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=fixed_config,
            )
            self.assertFalse(fixed_rejected["matched"])

            legacy = pd.read_csv(summary_path).drop(columns=[
                "source_estimator_config_schema_version",
                "source_estimator_config_sha256",
                "source_estimator_config_json",
            ])
            legacy_rejected = core._validation_product_status(
                legacy, order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(legacy_rejected["matched"])
            self.assertIn("missing_fields", legacy_rejected["reason"])

            raw_path.write_text("corrupt\n", encoding="utf-8")
            corrupted = core._validation_product_status(
                products["summary"], order=102, library_id="library-v1",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility=runtime_drp,
                blaze_source="l2", reference_quality_counts=quality_counts,
                calibration_metadata={"status": "not_applied"},
                estimator_config=estimator_config,
            )
            self.assertFalse(corrupted["matched"])
            self.assertIn("sidecar_hash_mismatch", corrupted["reason"])

    def test_warning_target_is_never_publication_validated(self):
        class Spectrum:
            object = "target"
            drp_version = "1.3.0"
            instrument = "NEID"
            observation_type = "SCI"
            data_level = 2
            observing_mode = "HR"
            dq_status = "unvalidated"
            dq_assessments = {
                "DQLEVEL1": {"status": "warning"},
                "DQLEVEL2": {"status": "pass"},
            }
            provenance = {"blaze_source": "l2"}
            blaze_source_used = "legacy_response"

        target = Spectrum()
        provenance = core.build_result_provenance(
            target=target, references=[],
            rv_metadata={
                "target": {"object": "target", "rv_kms": 0.0, "rv_source": "drp"},
                "references": [],
            },
            order=102, wavelength=np.linspace(8500.0, 8510.0, 20),
            deblazed=False, vsini=None, vsini_window=None, random_seed=0,
            feh_calibration_metadata={"status": "not_applied"},
            drp_compatibility={"publication_validated": True, "versions": {"1.3.0": 1}},
            library_id="library-v1",
            library_manifest={
                "manifest_present": True, "recognized": True,
                "deep_verified": True, "reference_pool_verified": True,
                "library_id": "library-v1",
            },
            composite_optimizer_diagnostics={
                "success": True, "nonunique_optimum": False,
            },
            pairwise_optimizer_diagnostics={"all_references_success": True},
        )
        self.assertIn(
            "target_dq_status_not_pass",
            provenance["overall_validation"]["reasons"],
        )
        self.assertIn(
            "non_l2_blaze_treatment_unvalidated",
            provenance["overall_validation"]["reasons"],
        )
        teff_validation = provenance["parameter_validation"]["teff"]
        self.assertNotIn("uncertainty_1sigma_k", teff_validation)
        self.assertIn("crossvalidation_residual_bias_k", teff_validation)
        self.assertIn("crossvalidation_sample_scatter_k", teff_validation)
        self.assertIn("crossvalidation_predictive_rmse_k", teff_validation)


if __name__ == "__main__":
    unittest.main()
