import tempfile
import json
import unittest
from contextlib import redirect_stdout
import io
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

import neidspecmatch.neidspecmatch as core


class _Spectrum:
    def __init__(self, name):
        self.object = name
        self.drp_version = "1.3.0"
        self.rv = 0.0
        self.rv_source = "library_catalog"
        self.blaze_source_used = "l2"


class _SpectrumList:
    def __init__(self, names):
        self.splist = [_Spectrum(name) for name in names]


def _valid_support_diagnostics():
    library = pd.DataFrame({
        "OBJECT_ID": [f"reference-{index}" for index in range(6)],
        "Teff": np.arange(6, dtype=float) * 200.0 + 3000.0,
        "[Fe/H]": np.linspace(-0.8, 0.6, 6),
        "log(g)": np.linspace(4.4, 5.0, 6),
        "vsini": np.arange(6, dtype=float),
        "rss": np.arange(6, dtype=float) + 1.0,
    })
    return core._build_support_diagnostics(
        library, library.iloc[:5].reset_index(drop=True), np.full(5, 0.2)
    )


class _FitResult:
    optimizer_diagnostics = {
        "success": True,
        "message": "converged",
        "nfev": 100,
        "nit": 10,
        "objective_log_likelihood": -1.0,
        "nonunique_optimum": False,
    }
    pairwise_optimizer_diagnostics = {
        "all_references_success": True,
        "failed_reference_count": 0,
    }
    support_diagnostics = _valid_support_diagnostics()


class CrossvalidationTests(unittest.TestCase):
    def test_checkpoint_fold_records_are_identity_and_finite_bound(self):
        identities = [{"object_id": "a"}, {"object_id": "b"}]
        valid = [{
            "index": 0, "target_identity": identities[0],
            "values": [1.0] * 7, "optimizer": {"success": True},
            "support": _FitResult.support_diagnostics,
        }]
        core._validate_checkpoint_folds(valid, identities)
        with self.assertRaises(ValueError):
            core._validate_checkpoint_folds(valid * 2, identities)
        invalid = [dict(valid[0], values=[1.0] * 6 + [np.nan])]
        with self.assertRaises(ValueError):
            core._validate_checkpoint_folds(invalid, identities)
        mismatched = [dict(valid[0], target_identity=identities[1])]
        with self.assertRaises(ValueError):
            core._validate_checkpoint_folds(mismatched, identities)

    def test_seven_column_result_construction_and_quiet_products(self):
        names = [f"star-{index}" for index in range(4)]
        library = pd.DataFrame({
            "OBJECT_ID": names,
            "Teff": [3500, 3700, 3900, 4100],
            "e_Teff": [50] * 4,
            "[Fe/H]": [-0.2, -0.1, 0.0, 0.1],
            "e_[Fe/H]": [0.05] * 4,
            "log(g)": [4.8] * 4,
            "e_log(g)": [0.05] * 4,
        })
        calls = []

        def fake_run(target, refs, ww, v, df_library, df_target, **kwargs):
            calls.append(target.object)
            row = df_target.iloc[0]
            # The real API returns nine values; CV consumes exactly the first seven.
            return (
                row["Teff"] + 10.0,
                row["[Fe/H]"] + 0.01,
                row["log(g)"] + 0.02,
                2.0,
                10.0,
                0.01,
                0.02,
                object(),
                _FitResult(),
            )

        with tempfile.TemporaryDirectory() as temp, \
             mock.patch.object(core, "run_specmatch", side_effect=fake_run):
            capture = io.StringIO()
            with redirect_stdout(capture):
                result = core.run_crossvalidation_for_orders(
                    order=102,
                    df_lib=library,
                    HLS=_SpectrumList(names),
                    outputdir=temp,
                    library_id="test-library",
                    verbose=0,
                    allow_unmanifested_library=True,
                )
            self.assertEqual(capture.getvalue(), "")
            self.assertEqual(calls, names)
            self.assertEqual(result.shape[0], 4)
            self.assertEqual(
                list(result.columns),
                [
                    "teff", "feh", "logg", "vsini", "d_teff", "d_feh",
                    "d_logg", "teff_true", "feh_true", "logg_true", "targetname",
                ],
            )
            summary = Path(result.attrs["products"]["summary"])
            self.assertTrue(summary.is_file())
            summary_frame = pd.read_csv(summary)
            self.assertAlmostEqual(summary_frame.loc[0, "feh_raw_sigma"], 0.0)

    def test_atomic_checkpoint_resumes_only_missing_folds(self):
        names = [f"star-{index}" for index in range(4)]
        library = pd.DataFrame({
            "OBJECT_ID": names,
            "Teff": [3500, 3600, 3700, 3800],
            "e_Teff": [50] * 4,
            "[Fe/H]": [-0.2, -0.1, 0.0, 0.1],
            "e_[Fe/H]": [0.05] * 4,
            "log(g)": [4.8] * 4,
            "e_log(g)": [0.05] * 4,
        })

        def result_for(target):
            row = library.loc[library.OBJECT_ID == target.object].iloc[0]
            return (
                row.Teff + 1, row["[Fe/H]"] + 0.01, row["log(g)"] + 0.01,
                2.0, 1.0, 0.01, 0.01, object(), _FitResult(),
            )

        first_calls = []

        def interrupted(target, *args, **kwargs):
            first_calls.append(target.object)
            if len(first_calls) == 2:
                raise RuntimeError("interrupted")
            return result_for(target)

        with tempfile.TemporaryDirectory() as temp:
            spectra = _SpectrumList(names)
            with mock.patch.object(core, "run_specmatch", side_effect=interrupted):
                with self.assertRaises(RuntimeError):
                    core.run_crossvalidation_for_orders(
                        102, df_lib=library, HLS=spectra, outputdir=temp,
                        library_id="test-library", verbose=0,
                        allow_unmanifested_library=True,
                    )
            checkpoint_path = Path(temp) / "o102_crossval" / (
                "crossvalidation_checkpoint_o102.json"
            )
            partial_checkpoint = json.loads(checkpoint_path.read_text())
            key = partial_checkpoint["key"]
            self.assertEqual(key["checkpoint_schema_version"], 4)
            self.assertEqual(
                key["estimator_config_schema_version"],
                core.ESTIMATOR_CONFIG_SCHEMA_VERSION,
            )
            self.assertEqual(
                key["estimator_config_sha256"],
                core._hash_estimator_config(json.loads(
                    key["estimator_config_json"]
                )),
            )
            checkpoint_config = json.loads(key["estimator_config_json"])
            self.assertEqual(
                checkpoint_config["rv"]["reference_library_state"][
                    "reference_count"
                ],
                len(names),
            )
            self.assertEqual(
                checkpoint_config["blaze"]["reference_library_state"][
                    "sources"
                ],
                ["l2"],
            )
            with mock.patch.object(core, "run_specmatch", side_effect=result_for):
                with self.assertRaisesRegex(ValueError, "checkpoint is incompatible"):
                    core.run_crossvalidation_for_orders(
                        102, df_lib=library, HLS=spectra, outputdir=temp,
                        library_id="test-library", verbose=0,
                        allow_unmanifested_library=True, maxvsini=31.0,
                    )
            spectra.splist[0].rv = 1.0
            with mock.patch.object(core, "run_specmatch", side_effect=result_for):
                with self.assertRaisesRegex(ValueError, "checkpoint is incompatible"):
                    core.run_crossvalidation_for_orders(
                        102, df_lib=library, HLS=spectra, outputdir=temp,
                        library_id="test-library", verbose=0,
                        allow_unmanifested_library=True,
                    )
            spectra.splist[0].rv = 0.0
            resumed_calls = []

            def resumed(target, *args, **kwargs):
                resumed_calls.append(target.object)
                return result_for(target)

            with mock.patch.object(core, "run_specmatch", side_effect=resumed):
                frame = core.run_crossvalidation_for_orders(
                    102, df_lib=library, HLS=spectra, outputdir=temp,
                    library_id="test-library", verbose=0,
                    allow_unmanifested_library=True,
                )
            self.assertEqual(first_calls, names[:2])
            self.assertEqual(resumed_calls, names[1:])
            self.assertEqual(len(frame), 4)
            checkpoint = json.loads(Path(frame.attrs["checkpoint"]).read_text())
            self.assertEqual(checkpoint["status"], "complete")
            self.assertEqual(len(checkpoint["folds"]), 4)
            self.assertTrue(checkpoint["optimizer"]["all_folds_success"])

    def test_support_counts_are_population_specific(self):
        frame = pd.DataFrame({
            "teff": np.arange(8) + 3500.0,
            "feh": np.zeros(8),
            "logg": np.full(8, 4.8),
            "vsini": np.ones(8),
            "d_teff": np.ones(8),
            "d_feh": np.linspace(-0.1, 0.1, 8),
            "d_logg": np.linspace(-0.05, 0.05, 8),
            "teff_true": [3300, 3500, 3700, 3900, 4700, 4900, 5100, 5300],
        })
        normal_library = pd.DataFrame({
            "OBJECT_ID": [f"normal-{index}" for index in range(7)],
            "Teff": np.arange(7, dtype=float) * 200.0 + 3000.0,
            "[Fe/H]": np.linspace(-0.6, 0.6, 7),
            "log(g)": np.linspace(4.4, 5.0, 7),
            "vsini": np.arange(7, dtype=float),
            "rss": np.arange(7, dtype=float) + 1.0,
        })
        edge_library = normal_library.copy()
        edge_library["OBJECT_ID"] = [f"edge-{index}" for index in range(7)]
        for column in ("Teff", "[Fe/H]", "log(g)"):
            edge_library.loc[1, column] = edge_library.loc[0, column]

        def support_receipt(*, edge, vertex):
            library = edge_library if edge else normal_library
            selected = (
                library.iloc[:5] if edge else library.iloc[1:6]
            ).reset_index(drop=True)
            if vertex:
                weights = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
            elif edge:
                weights = np.array([0.5, 0.5, 0.0, 0.0, 0.0])
            else:
                weights = np.full(5, 0.2)
            return core._build_support_diagnostics(
                library, selected, weights
            )

        supports = [
            support_receipt(
                edge=index in {0, 4, 5},
                vertex=index in {0, 1, 4, 6},
            )
            for index in range(8)
        ]
        with tempfile.TemporaryDirectory() as temp:
            raw = Path(temp) / "raw.csv"
            frame.to_csv(raw, index=False)
            paths = core.write_crossvalidation_products(
                frame, order=102, outputdir=temp,
                library_id="library-v1", cv_source=str(raw),
                population_rules=[
                    core.PopulationRule("all"),
                    core.PopulationRule("cool", teff_max=4500),
                ],
                source_fold_support=supports,
            )
            summary = pd.read_csv(paths["summary"]).set_index("population")
        self.assertEqual(summary.loc["all", "source_library_edge_folds"], 3)
        self.assertEqual(summary.loc["cool", "source_library_edge_folds"], 1)
        self.assertEqual(summary.loc["all", "source_simplex_vertex_folds"], 4)
        self.assertEqual(summary.loc["cool", "source_simplex_vertex_folds"], 2)


if __name__ == "__main__":
    unittest.main()
