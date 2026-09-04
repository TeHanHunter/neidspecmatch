import unittest
import json
import tempfile
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

import neidspecmatch.neidspecmatch as core
from neidspecmatch.neidspecmatch import get_data_ready


class _Spectrum:
    def __init__(self, name, rv, source):
        self.object = name
        self.rv = rv
        self.rv_source = source
        self.redshift_calls = []
        self.provenance = {"rv_source": source}
        self.drp_version = "1.3.0"
        self.blaze_source = "l2"
        self.blaze_source_used = "l2"
        self.deblaze_calls = []
        self.deblaze_hook = None

    def deblaze(self, force=False):
        self.deblaze_calls.append(bool(force))
        self.blaze_source_used = self.blaze_source
        if self.deblaze_hook is not None:
            self.deblaze_hook(bool(force))
        return None

    def redshift(self, rv):
        self.redshift_calls.append(float(rv))

    def rvabs_for_orders(self, *args, **kwargs):
        raise AssertionError("NEIDSpecMatch must not rerun the generic mask CCF")

    def resample_order(self, wavelength, **kwargs):
        return np.ones_like(wavelength), np.ones_like(wavelength) * 0.01


class _SpectrumList:
    def __init__(self, spectra):
        self.splist = spectra


class RVProvenanceTests(unittest.TestCase):
    def test_saved_result_preserves_nested_custom_ccf_provenance(self):
        unsupported = object()

        class Spectrum:
            object = "target"
            drp_version = "v1.5.2"
            instrument = "NEID"
            observation_type = "Sci"
            data_level = 2
            observing_mode = "HR"
            dq_status = "pass"
            dq_assessments = {}
            blaze_source_used = "l2"
            rv = 4.5
            rv_source = "explicit custom CCF"
            provenance = {
                "blaze_source": "l2",
                "rv_source": rv_source,
                "rv_mask_path": "/private/observer/masks/ESPRESSO_M3.fits",
                "rv_internal_object": unsupported,
                "rv_order_indices": [55, 56, 91],
                "rv_order_values_kms": np.array([4.4, 4.5, 4.6]),
                "rv_ccf_diagnostics": [
                    {
                        "order_index": np.int64(55),
                        "fitted_rv_kms": np.float64(4.4),
                        "fit_success": np.bool_(True),
                    }
                ],
                "rv_mask_wavelength_policy": {
                    "schema_version": 1,
                    "mask_sha256": "a" * 64,
                    "input_medium": "air",
                    "ccf_medium": "vacuum",
                    "air_to_vacuum_applied": True,
                    "policy_source": "neidspec_sha256_registry",
                },
            }

        provenance = core.build_result_provenance(
            target=Spectrum(), references=[],
            rv_metadata={
                "target": {
                    "object": "target", "rv_kms": 4.5,
                    "rv_source": "explicit custom CCF",
                },
                "references": [],
            },
            order=102, wavelength=np.linspace(8500.0, 8510.0, 20),
            deblazed=False, vsini=None, vsini_window=None, random_seed=0,
            feh_calibration_metadata={"status": "not_applied"},
            drp_compatibility={
                "publication_validated": True,
                "versions": {"v1.5.2": 1},
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
            pairwise_optimizer_diagnostics={"all_references_success": True},
        )
        result = {"provenance": provenance}

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            path.write_text(
                json.dumps(result, allow_nan=False), encoding="utf-8"
            )
            serialized = path.read_text(encoding="utf-8")
            saved = json.loads(serialized)["provenance"]["target"][
                "rv_diagnostics"
            ]

        self.assertEqual(saved["rv_order_indices"], [55, 56, 91])
        self.assertEqual(saved["rv_order_values_kms"], [4.4, 4.5, 4.6])
        self.assertEqual(saved["rv_ccf_diagnostics"][0]["order_index"], 55)
        self.assertTrue(saved["rv_ccf_diagnostics"][0]["fit_success"])
        self.assertEqual(
            saved["rv_mask_wavelength_policy"]["input_medium"], "air"
        )
        self.assertEqual(
            saved["rv_mask_wavelength_policy"]["ccf_medium"], "vacuum"
        )
        self.assertEqual(saved["rv_mask_basename"], "ESPRESSO_M3.fits")
        self.assertNotIn("rv_mask_path", saved)
        self.assertNotIn("rv_internal_object", saved)
        self.assertNotIn("/private/observer", serialized)

    def test_custom_ccf_mask_medium_is_forwarded_to_neidspec(self):
        captured = {}

        class StopAfterConstruction(Exception):
            pass

        def fake_spectrum(filename, **kwargs):
            captured.update(kwargs)
            raise StopAfterConstruction

        with mock.patch.object(
            core.neidspec, "NEIDSpectrum", side_effect=fake_spectrum
        ):
            with self.assertRaises(StopAfterConstruction):
                core.run_specmatch_for_orders(
                    "target.fits",
                    "target",
                    HLS=_SpectrumList([]),
                    orders=[102],
                    rv_source="custom_ccf",
                    ccf_mask_path="custom-mask.txt",
                    ccf_mask_medium="air",
                    verbose=0,
                )
        self.assertEqual(captured["ccf_mask_path"], "custom-mask.txt")
        self.assertEqual(captured["ccf_mask_medium"], "air")

    def test_ccf_mask_medium_is_rejected_for_drp_rv(self):
        with self.assertRaisesRegex(ValueError, "only valid with"):
            core.run_specmatch_for_orders(
                "target.fits",
                "target",
                HLS=_SpectrumList([]),
                orders=[102],
                rv_source="drp",
                ccf_mask_medium="air",
                verbose=0,
            )

    def test_ccf_mask_path_is_rejected_for_drp_rv(self):
        with self.assertRaisesRegex(ValueError, "ccf_mask_path is only valid"):
            core.run_specmatch_for_orders(
                "target.fits",
                "target",
                HLS=_SpectrumList([]),
                orders=[102],
                rv_source="drp",
                ccf_mask_path="ignored-mask.txt",
                verbose=0,
            )

    def test_default_preserves_neidspec_rvs_without_ccf(self):
        target = _Spectrum("target", 88.4, "NEID DRP HDU12 CCFJDSUM")
        references = [
            _Spectrum("ref-a", -3.2, "NEID DRP HDU12 CCFJDSUM"),
            _Spectrum("ref-b", 14.5, "NEID DRP HDU12 CCFJDSUM"),
        ]
        result = get_data_ready(
            target,
            _SpectrumList(references),
            np.linspace(5000, 5001, 8),
            np.linspace(-10, 10, 5),
            order=55,
            verbose=0,
            return_rv_metadata=True,
        )
        metadata = result[-1]
        self.assertEqual(target.redshift_calls, [])
        self.assertEqual([item.redshift_calls for item in references], [[], []])
        self.assertEqual(metadata["target"]["rv_kms"], 88.4)
        self.assertEqual(
            metadata["target"]["rv_source"], "NEID DRP HDU12 CCFJDSUM"
        )

    def test_explicit_target_rv_is_applied_once_and_recorded(self):
        target = _Spectrum("target", 1.0, "old")
        reference = _Spectrum("ref", 2.0, "library DRP")
        result = get_data_ready(
            target,
            _SpectrumList([reference]),
            np.linspace(5000, 5001, 8),
            np.linspace(-10, 10, 5),
            order=55,
            absrv=12.25,
            verbose=0,
            return_rv_metadata=True,
        )
        self.assertEqual(target.redshift_calls, [12.25])
        self.assertEqual(target.rv, 12.25)
        self.assertEqual(target.rv_source, "supplied to NEIDSpecMatch")
        self.assertEqual(result[-1]["target"]["rv_source"], "supplied to NEIDSpecMatch")

    def test_run_applies_all_rv_overrides_before_verification_and_ranking(self):
        target = _Spectrum("target", 1.0, "target DRP")
        references = [
            _Spectrum("ref-{}".format(index), float(index), "library DRP")
            for index in range(5)
        ]
        target._neidspecmatch_loaded_state_sha256 = (
            core._spectrum_state_fingerprint(target)
        )
        for reference in references:
            reference._neidspecmatch_loaded_state_sha256 = (
                core._spectrum_state_fingerprint(reference)
            )

        target_rv = 12.25
        reference_rvs = np.array([-4.0, -2.0, 0.0, 2.0, 4.0])
        wavelength = np.linspace(8500.0, 8501.0, 20)
        comparison_mask = np.zeros(wavelength.size, dtype=bool)
        comparison_mask[2:18] = True
        sequence = []
        target.deblaze_hook = lambda force: sequence.append(
            "blaze" if force else "unexpected_cached_blaze"
        )

        diagnostics = {
            "success": True,
            "linear_rank": 6,
            "objective_rss": 1.0,
            "valid_pixel_count": int(np.count_nonzero(comparison_mask)),
            "valid_pixel_fraction": float(np.mean(comparison_mask)),
        }
        pairwise = pd.DataFrame({
            "index": np.arange(5),
            "OBJECT_ID": [item.object for item in references],
            "rss": np.arange(1.0, 6.0),
            "poly_params": [np.array([1.0, 0, 0, 0, 0, 0])] * 5,
            "vsini": np.arange(5.0),
            "optimizer_diagnostics": [diagnostics.copy() for _ in references],
        })
        library = pd.DataFrame({
            "OBJECT_ID": [item.object for item in references],
            "Teff": np.linspace(3200.0, 4000.0, 5),
            "[Fe/H]": np.linspace(-0.2, 0.2, 5),
            "log(g)": np.linspace(4.7, 5.1, 5),
        })

        def assert_current_rvs(stage):
            sequence.append(stage)
            self.assertEqual(target.rv, target_rv)
            self.assertTrue(np.array_equal(
                [item.rv for item in references], reference_rvs
            ))
            self.assertEqual(target.redshift_calls, [target_rv])
            self.assertEqual(
                [item.redshift_calls for item in references],
                [[value] for value in reference_rvs],
            )
            for spectrum in [target, *references]:
                self.assertEqual(
                    spectrum._neidspecmatch_loaded_state_sha256,
                    core._spectrum_state_fingerprint(spectrum),
                )

        def fake_verify(actual_target, actual_references, library_path):
            self.assertIs(actual_target, target)
            self.assertEqual(list(actual_references), references)
            assert_current_rvs("verify")
            return {"reference_pool_verified": True}

        def fake_pairwise(ww, actual_target, actual_references, **kwargs):
            self.assertTrue(kwargs["return_comparison_mask"])
            assert_current_rvs("rank")
            return (
                pairwise.copy(), pairwise.copy(),
                _SpectrumList(references), comparison_mask.copy(),
            )

        def fake_ready(actual_target, actual_references, ww, velocities,
                       **kwargs):
            self.assertNotIn("absrv", kwargs)
            self.assertNotIn("reference_rvs", kwargs)
            assert_current_rvs("stage2")
            target_flux = np.ones_like(ww)
            target_error = np.full_like(ww, 0.01)
            reference_fluxes = np.vstack([
                target_flux + index * 1e-3 for index in range(5)
            ])
            reference_errors = np.full_like(reference_fluxes, 0.02)
            metadata = {
                "target": {
                    "object": target.object,
                    "rv_kms": target.rv,
                    "rv_source": target.rv_source,
                },
                "references": [
                    {
                        "object": item.object,
                        "rv_kms": item.rv,
                        "rv_source": item.rv_source,
                    }
                    for item in references
                ],
            }
            return (
                target_flux, target_error, reference_fluxes,
                reference_errors, metadata,
            )

        class FakeFit:
            def __init__(self, objective, teffs, fehs, loggs, vsinis, **kwargs):
                self.lpf = objective
                self.teffs = np.asarray(teffs)
                self.fehs = np.asarray(fehs)
                self.loggs = np.asarray(loggs)
                self.vsinis = np.asarray(vsinis)
                self.feh_calibration_metadata = {"status": "raw"}

            def minimize_PyDE(self, **kwargs):
                self.min_pv = np.full(4, 0.2)
                self.optimizer_diagnostics = {
                    "success": True,
                    "nonunique_optimum": False,
                }

            def calculate_stellar_parameters(self, weights):
                self.teff = float(weights @ self.teffs)
                self.feh_raw = float(weights @ self.fehs)
                self.feh = self.feh_raw
                self.logg = float(weights @ self.loggs)
                self.vsini = float(weights @ self.vsinis)

        captured_provenance = {}

        def fake_provenance(**kwargs):
            captured_provenance.update(kwargs)
            return {}

        with tempfile.TemporaryDirectory() as temp:
            library_root = Path(temp) / "library"
            library_root.mkdir()
            manifest_path = library_root / "library_manifest.json"
            manifest_path.write_text("{}\n", encoding="utf-8")
            receipt = {
                "recognized": True,
                "deep_verified": True,
                "manifest_sha256": core.sha256_file(manifest_path),
            }
            output = Path(temp) / "results"
            with mock.patch.object(
                    core, "_verify_reference_pool_against_manifest",
                    side_effect=fake_verify), \
                 mock.patch.object(
                    core, "chi2spectraPolyLoop", side_effect=fake_pairwise), \
                 mock.patch.object(
                    core, "get_data_ready", side_effect=fake_ready), \
                 mock.patch.object(core, "FitLinCombSpec", FakeFit), \
                 mock.patch.object(
                    core, "build_result_provenance",
                    side_effect=fake_provenance):
                result = core.run_specmatch(
                    target, references, wavelength,
                    np.linspace(-10.0, 10.0, 5), library,
                    savefolder=output, order=102, verbose=0,
                    absrv=target_rv, reference_rvs=reference_rvs,
                    library_id="test-library", library_path=library_root,
                    library_manifest=receipt,
                    drp_compatibility={"publication_validated": True},
                )

        self.assertEqual(sequence, ["verify", "blaze", "rank", "stage2"])
        self.assertEqual(target.deblaze_calls, [True])
        self.assertEqual(
            [item.deblaze_calls for item in references], [[True]] * 5
        )
        self.assertTrue(np.array_equal(result[-1].lpf.valid_mask, comparison_mask))
        self.assertEqual(target.redshift_calls, [target_rv])
        self.assertEqual(
            [item.redshift_calls for item in references],
            [[value] for value in reference_rvs],
        )
        estimator = captured_provenance["estimator_config"]
        self.assertEqual(
            estimator["rv"]["reference_library_state"],
            core._reference_rv_state(references),
        )
        self.assertEqual(
            estimator["blaze"]["reference_library_state"],
            core._reference_blaze_state(references),
        )

    def test_library_loader_uses_catalog_rv_only_for_invalid_drp_rv(self):
        with tempfile.TemporaryDirectory() as temp:
            fits_dir = Path(temp)
            for basename in ("good.fits", "blank.fits"):
                (fits_dir / basename).write_bytes(b"")
            catalog = pd.DataFrame({
                "OBJECT_ID": ["good-star", "fallback-star"],
                "basenames": ["good.fits", "blank.fits"],
                "rv": [1.2, -0.7434],
            })
            calls = []

            class FakeLoaded:
                def __init__(self, filename, targetname, rv_source, rv=None,
                             verbose=False, **kwargs):
                    calls.append((Path(filename).name, targetname, rv_source, rv))
                    if Path(filename).name == "blank.fits" and rv_source == "drp":
                        raise ValueError("blank.fits has an invalid DRP RV; provide an explicit rv")
                    self.object = "archive-header-name"
                    self.rv = 3.0 if rv is None else rv
                    self.rv_source = rv_source
                    self.provenance = {"rv_source": rv_source}

            with mock.patch.object(core.neidspec, "NEIDSpectrum", FakeLoaded), \
                 mock.patch.object(
                     core.neidspec, "NEIDSpecList",
                     side_effect=lambda splist: type("List", (), {"splist": splist})(),
                 ):
                loaded = core.load_reference_library(catalog, temp, verbose=0)
            self.assertEqual(len(loaded.splist), 2)
            self.assertEqual(loaded.splist[0].rv_source, "drp")
            self.assertEqual(loaded.splist[1].rv, -0.7434)
            self.assertEqual(loaded.splist[1].rv_source, "library_catalog")
            self.assertEqual(loaded.splist[0].object, "good-star")
            self.assertEqual(loaded.splist[1].object, "fallback-star")
            self.assertEqual(
                loaded.splist[0].fits_object, "archive-header-name"
            )
            self.assertEqual(
                loaded.splist[0].provenance["library_identity_source"],
                "catalog_OBJECT_ID",
            )
            self.assertEqual(
                calls,
                [
                    ("good.fits", "good-star", "drp", None),
                    ("blank.fits", "fallback-star", "drp", None),
                    ("blank.fits", "fallback-star", "supplied", -0.7434),
                ],
            )

    def test_drp_guard_rejects_mixed_and_unknown_by_default(self):
        target = _Spectrum("target", 1.0, "drp")
        reference = _Spectrum("ref", 2.0, "drp")
        status = core.validate_drp_compatibility(target, [reference])
        self.assertTrue(status["publication_validated"])
        reference.drp_version = "1.5.3"
        with self.assertRaises(ValueError):
            core.validate_drp_compatibility(target, [reference])
        override = core.validate_drp_compatibility(
            target, [reference], allow_mixed_drp=True
        )
        self.assertFalse(override["publication_validated"])
        target.drp_version = "unknown"
        reference.drp_version = "unknown"
        with self.assertRaises(ValueError):
            core.validate_drp_compatibility(target, [reference])

    def test_drp_guard_accepts_patch_versions_within_one_minor_series(self):
        target = _Spectrum("target", 1.0, "drp")
        reference = _Spectrum("ref", 2.0, "drp")
        target.drp_version = "v1.5.2"
        reference.drp_version = "v1.5.3"
        status = core.validate_drp_compatibility(target, [reference])
        self.assertTrue(status["publication_validated"])
        self.assertTrue(status["mixed_patch_versions"])
        self.assertFalse(status["mixed_drp"])
        self.assertEqual(status["minor_series"], {"1.5": 2})

    def test_drp_guard_rejects_malformed_version(self):
        target = _Spectrum("target", 1.0, "drp")
        reference = _Spectrum("ref", 2.0, "drp")
        target.drp_version = "v1.5"
        reference.drp_version = "v1.5.3"
        with self.assertRaises(ValueError):
            core.validate_drp_compatibility(target, [reference])

    def test_library_loader_closes_prior_spectra_on_later_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            Path(temp, "one.fits").write_bytes(b"")
            catalog = pd.DataFrame({
                "OBJECT_ID": ["one", "missing"],
                "basenames": ["one.fits", "missing.fits"],
                "rv": [0.0, 0.0],
            })
            opened = []

            class Loaded:
                def __init__(self, filename, targetname, **kwargs):
                    self.object = targetname
                    self.provenance = {}
                    self.closed = False
                    opened.append(self)

                def close(self):
                    self.closed = True

            with mock.patch.object(core.neidspec, "NEIDSpectrum", Loaded):
                with self.assertRaises(FileNotFoundError):
                    core.load_reference_library(catalog, temp, verbose=0)
            self.assertTrue(opened[0].closed)


if __name__ == "__main__":
    unittest.main()
