import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import neidspec
import neidspecmatch.neidspecmatch as core

from neidspecmatch.calibration import (
    CalibrationError,
    FeHCalibrationSet,
    PopulationRule,
    UnvalidatedCalibrationError,
    build_calibration_products,
    raw_crossvalidation_summary,
)
from neidspecmatch.neidspecmatch import (
    FitLinCombSpec,
    LPFunctionLinComb,
    detrend_feh,
)


class CalibrationTests(unittest.TestCase):
    def artifact_evidence(self):
        return {
            "source_library_manifest_sha256": "1" * 64,
            "source_pipeline_sha256": "2" * 64,
            "source_all_folds_optimizer_success": True,
            "source_failed_optimizer_folds": 0,
            "source_neidspec_version": "0.2.0",
            "source_neidspec_sha256": "3" * 64,
            "source_result_schema_version": core.RESULT_SCHEMA_VERSION,
            "source_blaze_source": "neid_l2_blaze",
            "source_drp_version_set": ("1.3.0",),
        }

    def synthetic_crossval(self):
        true = np.linspace(-0.5, 0.5, 16)
        noise = np.array([
            -0.01, 0.01, -0.005, 0.005, 0.0, -0.008, 0.008, 0.003,
            -0.003, 0.006, -0.006, 0.004, -0.004, 0.002, -0.002, 0.0,
        ])
        recovered = 0.7 * true - 0.02 + noise
        return pd.DataFrame({
            "feh": recovered,
            "feh_true": true,
            "d_feh": recovered - true,
            "teff_true": np.linspace(3100, 5600, len(true)),
            "d_teff": np.linspace(-20, 20, len(true)),
            "d_logg": np.linspace(-0.05, 0.05, len(true)),
        })

    def test_stale_calibration_cannot_be_selected(self):
        with self.assertRaises(UnvalidatedCalibrationError):
            detrend_feh(0.1)
        with self.assertRaises(CalibrationError):
            PopulationRule("bad", teff_min=np.nan)

    def test_cv_summary_separates_bias_scatter_and_predictive_rmse(self):
        residuals = np.array([1.0, 2.0, 3.0, 4.0])
        frame = pd.DataFrame({
            "d_teff": residuals,
            "d_feh": residuals / 10.0,
            "d_logg": residuals / 100.0,
        })
        summary = raw_crossvalidation_summary(
            frame, order=102, library_id="test-library",
            cv_source="crossval.csv",
        )
        self.assertAlmostEqual(summary["teff_residual_bias"], 2.5)
        self.assertAlmostEqual(
            summary["teff_sample_scatter"], np.std(residuals, ddof=1)
        )
        self.assertAlmostEqual(
            summary["teff_predictive_rmse"],
            np.sqrt(np.mean(residuals ** 2)),
        )
        self.assertNotEqual(
            summary["teff_sample_scatter"],
            summary["teff_predictive_rmse"],
        )
        # The old columns remain exact, explicitly documented aliases.
        self.assertEqual(summary["teff_bias"], summary["teff_residual_bias"])
        self.assertEqual(summary["teff_sigma"], summary["teff_sample_scatter"])
        self.assertEqual(summary["teff_rmse"], summary["teff_predictive_rmse"])

    def test_raw_metallicity_is_default(self):
        w = np.linspace(0, 1, 5)
        target = np.ones(5)
        refs = np.ones((5, 5))
        lpf = LPFunctionLinComb(w, target, np.ones(5), refs, refs)
        fit = FitLinCombSpec(
            lpf,
            teffs=np.full(5, 4000.0),
            fehs=np.array([-0.2, -0.1, 0.0, 0.1, 0.2]),
            loggs=np.full(5, 4.7),
            vsinis=np.full(5, 2.0),
            calibrate_feh=False,
            verbose=0,
        )
        weights = np.full(5, 0.2)
        fit.calculate_stellar_parameters(weights)
        self.assertAlmostEqual(fit.feh, fit.feh_raw)
        self.assertEqual(fit.feh_calibration_metadata["status"], "not_applied")

    def test_crossfit_artifact_and_metrics(self):
        frame = self.synthetic_crossval()
        artifact, summary, validation = build_calibration_products(
            frame,
            order=102,
            library_id="test-library",
            cv_source="archived-crossval.csv",
            **self.artifact_evidence(),
        )
        calibration = artifact.calibrations[0]
        self.assertFalse(calibration.validated)
        self.assertEqual(calibration.validation_method,
                         "spectral_leave_one_out_plus_calibration_leave_one_out")
        self.assertIn("correlated", calibration.validation_limitation)
        self.assertLess(calibration.calibrated_sigma, calibration.raw_sigma)
        self.assertEqual(len(validation), len(frame))
        self.assertIn("d_feh_calibrated_crossfit", validation)
        self.assertEqual(summary.loc[0, "calibration_status"], "unvalidated")
        self.assertIn("not strict nested", calibration.validation_reason)
        self.assertFalse(artifact.publication_validation_eligible)
        with self.assertRaises(UnvalidatedCalibrationError):
            artifact.select(102, library_id="test-library")

        with tempfile.TemporaryDirectory() as temp:
            filename = Path(temp) / "calibration.json"
            artifact.save(filename)
            loaded = FeHCalibrationSet.load(filename)
            self.assertEqual(loaded.canonical_sha256, artifact.canonical_sha256)
            with self.assertRaises(UnvalidatedCalibrationError):
                loaded.select(102)
            payload = json.loads(filename.read_text())
            self.assertEqual(payload["schema_version"], 2)
            self.assertEqual(payload["artifact_sha256"], artifact.canonical_sha256)
            payload["calibrations"][0]["calibrated_sigma"] *= 2
            filename.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(CalibrationError):
                FeHCalibrationSet.load(filename)

    def test_order_library_and_population_mismatches_fail(self):
        frame = self.synthetic_crossval()
        artifact, _, _ = build_calibration_products(
            frame,
            order=102,
            library_id="test-library",
            cv_source="crossval.csv",
            populations=[PopulationRule("cool", teff_max=4500)],
            **self.artifact_evidence(),
        )
        with self.assertRaises(CalibrationError):
            artifact.select(101, population="cool", library_id="test-library")
        with self.assertRaises(CalibrationError):
            artifact.select(102, population="cool", library_id="other-library")
        # Bounded populations cannot be chosen automatically using recovered Teff.
        with self.assertRaises(CalibrationError):
            artifact.select(102, teff=3800, library_id="test-library")
        with self.assertRaises(CalibrationError):
            artifact.select(
                102, teff=np.nan, population="cool",
                library_id="test-library",
            )
        with self.assertRaises(UnvalidatedCalibrationError):
            artifact.select(
                102, teff=3800, population="cool", library_id="test-library"
            )
        self.assertEqual(artifact.calibrations[0].population.name, "cool")

    def test_editable_validation_booleans_cannot_bypass_scientific_gates(self):
        artifact, _, _ = build_calibration_products(
            self.synthetic_crossval(), order=102, library_id="test-library",
            cv_source="crossval.csv",
            **self.artifact_evidence(),
        )
        payload = artifact.to_dict()
        payload["calibrations"][0]["n_stars"] = 1
        canonical = {key: value for key, value in payload.items()
                     if key != "artifact_sha256"}
        import hashlib
        payload["artifact_sha256"] = hashlib.sha256(json.dumps(
            canonical, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("utf-8")).hexdigest()
        with self.assertRaises(CalibrationError):
            FeHCalibrationSet.from_mapping(payload)

    def test_improving_artifact_without_successful_cv_evidence_is_unusable(self):
        artifact, _, _ = build_calibration_products(
            self.synthetic_crossval(), order=102, library_id="test-library",
            cv_source="crossval.csv",
        )
        self.assertFalse(artifact.source_evidence_validated)
        with self.assertRaises(UnvalidatedCalibrationError):
            artifact.select(102)

    def test_even_self_asserted_validated_leaky_artifact_cannot_be_applied(self):
        artifact, _, _ = build_calibration_products(
            self.synthetic_crossval(), order=102, library_id="test-library",
            cv_source="crossval.csv", **self.artifact_evidence(),
        )
        with self.assertRaisesRegex(
            CalibrationError, "validated flag disagrees"
        ):
            replace(
                artifact.calibrations[0],
                validated=True,
                validation_reason="passed",
            )

    def test_calibration_runtime_context_is_exact_and_rejects_mixed_drp(self):
        manifest_sha = "a" * 64
        artifact, _, _ = build_calibration_products(
            self.synthetic_crossval(), order=102, library_id="test-library",
            cv_source="crossval.csv",
            source_library_manifest_sha256=manifest_sha,
            source_pipeline_sha256=core._pipeline_source_fingerprint(),
            source_all_folds_optimizer_success=True,
            source_failed_optimizer_folds=0,
            source_neidspec_version=getattr(neidspec, "__version__", "unknown"),
            source_neidspec_sha256=core._neidspec_source_fingerprint(),
            source_result_schema_version=core.RESULT_SCHEMA_VERSION,
            source_blaze_source="l2",
            source_drp_version_set=("1.3.0",),
        )
        returned = core._validate_calibration_runtime_context(
            artifact, library_id="test-library",
            library_manifest={"manifest_sha256": manifest_sha},
            drp_compatibility={
                "publication_validated": True,
                "versions": {"1.3.0": 79},
            },
            blaze_source="l2",
        )
        self.assertIs(returned, artifact)
        with self.assertRaises(UnvalidatedCalibrationError):
            core._validate_calibration_runtime_context(
                artifact, library_id="test-library",
                library_manifest={"manifest_sha256": "b" * 64},
                drp_compatibility={
                    "publication_validated": True,
                    "versions": {"1.3.0": 79},
                },
                blaze_source="l2",
            )
        with self.assertRaises(UnvalidatedCalibrationError):
            core._validate_calibration_runtime_context(
                artifact, library_id="test-library",
                library_manifest={"manifest_sha256": manifest_sha},
                drp_compatibility={
                    "publication_validated": False,
                    "versions": {"1.3.0": 78, "1.5.3": 1},
                },
                blaze_source="l2",
            )

    def test_near_singular_calibration_is_recorded_but_cannot_be_applied(self):
        true = np.linspace(-0.5, 0.5, 12)
        recovered = 0.001 * true
        frame = pd.DataFrame({
            "feh": recovered,
            "feh_true": true,
            "d_feh": recovered - true,
            "teff_true": np.linspace(3200, 5000, len(true)),
        })
        artifact, summary, _ = build_calibration_products(
            frame,
            order=93,
            library_id="test-library",
            cv_source="crossval.csv",
        )
        self.assertEqual(summary.loc[0, "calibration_status"], "unvalidated")
        with self.assertRaises(UnvalidatedCalibrationError):
            artifact.select(93)

    def test_finite_but_worsening_calibration_fails_performance_gate(self):
        true = np.linspace(-0.5, 0.5, 16)
        noise = np.array([0.03, -0.03] * 8)
        recovered = true + noise
        frame = pd.DataFrame({
            "feh": recovered,
            "feh_true": true,
            "d_feh": recovered - true,
            "teff_true": np.linspace(3200, 5200, len(true)),
        })
        artifact, summary, _ = build_calibration_products(
            frame,
            order=102,
            library_id="test-library",
            cv_source="crossval.csv",
        )
        self.assertGreater(
            summary.loc[0, "feh_calibrated_sigma"],
            summary.loc[0, "feh_raw_sigma"],
        )
        self.assertFalse(summary.loc[0, "feh_calibration_performance_passed"])
        self.assertIn("performance gate failed", summary.loc[0, "calibration_validation_reason"])
        with self.assertRaises(UnvalidatedCalibrationError):
            artifact.select(102)


if __name__ == "__main__":
    unittest.main()
