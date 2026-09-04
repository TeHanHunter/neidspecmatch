import unittest
from unittest import mock

import numpy as np
import pandas as pd
import scipy.optimize

import neidspecmatch.neidspecmatch as core
from neidspecmatch import rotbroad_help


class _Spectrum:
    def __init__(self, name):
        self.object = name

    def resample_order(self, wavelength, **kwargs):
        flux = np.ones_like(wavelength, dtype=float)
        return flux, np.full_like(flux, 0.01)


class VsiniTests(unittest.TestCase):
    def test_constraint_is_immutable_across_references(self):
        seen = []

        def fake_fit(*args, **kwargs):
            seen.append((
                kwargs["vsini"], kwargs["vsini_window"],
                kwargs["comparison_mask"].copy(),
            ))
            index = len(seen)
            return (
                float(index), 10.0 + index, np.zeros(6),
                {"success": True, "objective_rss": float(index)},
            )

        refs = [_Spectrum("ref-a"), _Spectrum("ref-b"), _Spectrum("ref-c")]
        fake_list = mock.Mock()
        fake_list.splist = refs
        with mock.patch.object(core, "chi2spectraPolyVsini", side_effect=fake_fit), \
             mock.patch.object(core.neidspec, "NEIDSpecList", return_value=fake_list):
            result, _, _, comparison_mask = core.chi2spectraPolyLoop(
                np.arange(8500.0, 8510.0, 0.01),
                _Spectrum("target"),
                refs,
                plot_chi=False,
                verbose=0,
                vsini=10.0,
                vsini_window=0.5,
                return_comparison_mask=True,
            )
        self.assertEqual([(item[0], item[1]) for item in seen], [(10.0, 0.5)] * 3)
        self.assertTrue(all(
            np.array_equal(item[2], comparison_mask) for item in seen
        ))
        self.assertGreater(np.count_nonzero(comparison_mask), 0)
        self.assertEqual(list(result["vsini"]), [11.0, 12.0, 13.0])

    def test_composite_honors_explicit_pairwise_mask_exactly(self):
        wavelength = np.arange(20.0)
        target = np.linspace(0.9, 1.1, wavelength.size)
        references = np.vstack([
            target + offset for offset in (-0.02, -0.01, 0.0, 0.01, 0.02)
        ])
        comparison_mask = np.zeros(wavelength.size, dtype=bool)
        comparison_mask[3:17:2] = True
        objective = core.LPFunctionLinComb(
            wavelength, target, np.full(wavelength.size, 0.01),
            references, np.full_like(references, 0.02),
            valid_mask=comparison_mask,
        )
        self.assertTrue(np.array_equal(objective.valid_mask, comparison_mask))

        unsupported = comparison_mask.copy()
        unsupported[0] = True
        references_with_gap = references.copy()
        references_with_gap[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "unsupported"):
            core.LPFunctionLinComb(
                wavelength, target, np.full(wavelength.size, 0.01),
                references_with_gap, np.full_like(references, 0.02),
                valid_mask=unsupported,
            )

    def test_composite_model_ignores_nan_in_zero_weight_reference(self):
        wavelength = np.arange(10.0)
        target = np.ones(wavelength.size)
        references = np.ones((5, wavelength.size))
        references[1, 2] = np.nan
        objective = core.LPFunctionLinComb(
            wavelength, target, np.full(wavelength.size, 0.01),
            references, np.full_like(references, 0.02),
        )

        inactive_gap_model = objective.compute_model([1.0, 0.0, 0.0, 0.0])
        self.assertTrue(np.isfinite(inactive_gap_model[2]))

        active_gap_model = objective.compute_model([0.0, 1.0, 0.0, 0.0])
        self.assertTrue(np.isnan(active_gap_model[2]))

    def test_exact_vsini_and_optimizer_smoke(self):
        wavelength = np.linspace(8564.5, 8565.5, 31)
        flux = 1.0 + 0.01 * (wavelength - wavelength.mean())
        error = np.full_like(flux, 0.01)
        objective = core.Chi2FunctionVsiniPolynomial(
            wavelength, flux, error, flux, error,
            maxvsini=30.0, vsini=7.25, vsini_window=0.0,
        )
        fit = core.FitTargetRefStarVsiniPolynomial(objective)
        with mock.patch.object(
                core.rotbroad_help, "broaden",
                side_effect=lambda w, f, vsini, u1: f):
            fit.minimize_AMOEBA(verbose=0)
        self.assertAlmostEqual(fit.min_pv[0], 7.25, places=12)
        self.assertTrue(np.isfinite(objective(fit.min_pv)))
        self.assertTrue(fit.optimizer_diagnostics["vsini_fixed_by_user"])
        self.assertFalse(fit.optimizer_diagnostics["vsini_boundary_hit"])

    def test_estimator_configuration_binds_vsini_and_exact_grid(self):
        wavelength = np.arange(8564.5, 8685.7, 0.01)

        class Reference:
            def __init__(self, object_id, rv, rv_source="library_catalog",
                         blaze_source="l2"):
                self.object = object_id
                self.rv = rv
                self.rv_source = rv_source
                self.blaze_source_used = blaze_source

        references = [
            Reference("reference-a", -1.25),
            Reference("reference-b", 2.5),
        ]
        rv_state = core._reference_rv_state(references)
        blaze_state = core._reference_blaze_state(references)
        common = dict(
            order=102,
            wavelength=wavelength,
            maxvsini=30.0,
            library_id="library-v1",
            library_manifest={"manifest_sha256": "a" * 64},
            deblazed=False,
            blaze_source="l2",
            reference_rv_state=rv_state,
            reference_blaze_state=blaze_state,
        )
        free = core._build_estimator_config(**common)
        fixed = core._build_estimator_config(
            **common, vsini=30.0, vsini_window=0.0
        )
        bounded = core._build_estimator_config(
            **common, vsini=20.0, vsini_window=3.0
        )
        self.assertEqual(free["vsini"]["mode"], "free")
        self.assertEqual(fixed["vsini"]["effective_bounds_kms"], [30.0, 30.0])
        self.assertEqual(bounded["vsini"]["effective_bounds_kms"], [17.0, 23.0])
        hashes = {
            core._hash_estimator_config(value)
            for value in (free, fixed, bounded)
        }
        self.assertEqual(len(hashes), 3)
        shifted_grid = core._build_estimator_config(
            **{**common, "wavelength": wavelength + 1e-8}
        )
        self.assertNotEqual(
            free["spectral_grid"]["wavelength_float64_le_sha256"],
            shifted_grid["spectral_grid"]["wavelength_float64_le_sha256"],
        )
        self.assertNotEqual(
            core._hash_estimator_config(free),
            core._hash_estimator_config(shifted_grid),
        )
        self.assertEqual(
            rv_state, core._reference_rv_state(list(reversed(references)))
        )
        self.assertEqual(
            blaze_state,
            core._reference_blaze_state(list(reversed(references))),
        )
        changed_rv = [
            Reference("reference-a", -1.0),
            Reference("reference-b", 2.5),
        ]
        changed_rv_config = core._build_estimator_config(
            **{**common, "reference_rv_state": core._reference_rv_state(
                changed_rv
            )}
        )
        self.assertNotEqual(
            core._hash_estimator_config(free),
            core._hash_estimator_config(changed_rv_config),
        )
        changed_blaze = [
            Reference("reference-a", -1.25, blaze_source="legacy_response"),
            Reference("reference-b", 2.5),
        ]
        changed_blaze_config = core._build_estimator_config(
            **{
                **common,
                "reference_blaze_state": core._reference_blaze_state(
                    changed_blaze
                ),
            }
        )
        self.assertNotEqual(
            core._hash_estimator_config(free),
            core._hash_estimator_config(changed_blaze_config),
        )
        duplicates = [Reference("same", 0.0), Reference("same", 1.0)]
        with self.assertRaises(ValueError):
            core._reference_rv_state(duplicates)
        with self.assertRaises(ValueError):
            core._reference_blaze_state(duplicates)

    def test_support_diagnostics_distinguish_vertex_from_library_edge(self):
        library = pd.DataFrame({
            "OBJECT_ID": [f"star-{index}" for index in range(7)],
            "Teff": [3000, 3200, 3400, 3600, 3800, 4000, 4200],
            "[Fe/H]": [-0.8, -0.3, 0.0, 0.2, 0.4, 0.6, 0.8],
            "log(g)": [4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0],
            "vsini": np.arange(7, dtype=float),
            "rss": np.arange(7, dtype=float) + 1.0,
        })
        selected = library.iloc[1:6].reset_index(drop=True)
        blend = core._build_support_diagnostics(
            library, selected, [0.5, 0.5, 0.0, 0.0, 0.0]
        )
        self.assertEqual(blend["simplex"]["active_count"], 2)
        self.assertAlmostEqual(blend["simplex"]["effective_component_count"], 2.0)
        self.assertFalse(blend["simplex"]["is_vertex"])
        self.assertFalse(blend["any_full_library_edge"])
        self.assertTrue(blend["publication_support_validated"])

        interior_vertex = core._build_support_diagnostics(
            library, selected, [0.0, 1.0, 0.0, 0.0, 0.0]
        )
        self.assertTrue(interior_vertex["simplex"]["is_vertex"])
        self.assertFalse(interior_vertex["any_full_library_edge"])
        self.assertTrue(interior_vertex["publication_support_validated"])

        edge_selected = library.iloc[:5].reset_index(drop=True)
        lower_vertex = core._build_support_diagnostics(
            library, edge_selected, [1.0, 0.0, 0.0, 0.0, 0.0]
        )
        self.assertTrue(lower_vertex["vertex_at_full_library_edge"])
        self.assertTrue(
            lower_vertex["parameters"]["teff"]["at_full_library_lower_edge"]
        )
        self.assertFalse(lower_vertex["publication_support_validated"])

    def test_support_edge_can_have_multiple_active_references(self):
        library = pd.DataFrame({
            "OBJECT_ID": ["low-a", "low-b", "mid-a", "mid-b", "high", "top"],
            "Teff": [3000, 3000, 3400, 3600, 4000, 4300],
            "[Fe/H]": [-0.5, 0.2, -0.1, 0.1, 0.4, 0.7],
            "log(g)": [4.6, 4.8, 4.7, 4.9, 5.0, 5.1],
            "vsini": [1, 2, 3, 4, 5, 6],
            "rss": [1, 2, 3, 4, 5, 6],
        })
        support = core._build_support_diagnostics(
            library, library.iloc[:5].copy(), [0.4, 0.6, 0.0, 0.0, 0.0]
        )
        self.assertEqual(support["simplex"]["active_count"], 2)
        self.assertFalse(support["simplex"]["is_vertex"])
        self.assertTrue(support["parameters"]["teff"]["support_limited"])
        self.assertFalse(support["parameters"]["feh_raw"]["support_limited"])
        self.assertAlmostEqual(
            support["parameters"]["feh_raw"]["raw_value"], -0.08
        )

    def test_exact_linear_combination_matches_scipy_and_repeats(self):
        wavelength = np.linspace(0, 1, 20)
        refs = np.array([
            1 + scale * wavelength for scale in (0.0, 0.1, -0.1, 0.2, -0.2)
        ])
        target = 0.2 * refs[0] + 0.2 * refs[1] + 0.2 * refs[2] + 0.2 * refs[3] + 0.2 * refs[4]
        lpf = core.LPFunctionLinComb(
            wavelength, target, np.full(20, 0.01), refs, np.ones_like(refs)
        )
        fit = core.FitLinCombSpec(lpf, verbose=0)
        fit.minimize_PyDE(
            npop=20, de_iter=12, mcmc=False, random_seed=42,
            tolerance=1e-5, verbose=0,
        )
        weights = np.asarray(lpf.get_pv_all(fit.min_pv))
        self.assertTrue(np.all(weights >= 0))
        self.assertAlmostEqual(np.sum(weights), 1.0, places=10)
        self.assertTrue(np.isfinite(fit.min_pv_lnval))
        self.assertEqual(
            fit.optimizer_diagnostics["name"],
            "enumerated_active_set_weighted_least_squares",
        )
        reference = scipy.optimize.minimize(
            lambda pv: -lpf(pv), np.full(4, 0.15), method="SLSQP",
            bounds=[(0.0, 1.0)] * 4,
            constraints={"type": "ineq", "fun": lambda pv: 1.0 - np.sum(pv)},
            options={"ftol": 1e-12, "maxiter": 500},
        )
        self.assertTrue(reference.success)
        self.assertLessEqual(-fit.min_pv_lnval, reference.fun + 1e-8)
        repeated = core.FitLinCombSpec(lpf, verbose=0)
        repeated.minimize_PyDE(mcmc=False, random_seed=999, verbose=0)
        self.assertTrue(np.array_equal(repeated.min_pv, fit.min_pv))

    def test_active_set_matches_analytic_simplex_projection(self):
        references = np.eye(5)
        target = np.array([-0.2, 0.2, 0.3, 0.1, 0.6])
        lpf = core.LPFunctionLinComb(
            np.arange(5.0), target, np.ones(5), references,
            np.ones_like(references),
        )
        fit = core.FitLinCombSpec(lpf, verbose=0)
        fit.minimize_PyDE(mcmc=False, verbose=0)
        weights = np.asarray(lpf.get_pv_all(fit.min_pv))
        self.assertTrue(np.allclose(
            weights, [0.0, 0.15, 0.25, 0.05, 0.55], atol=1e-12
        ))
        self.assertAlmostEqual(
            fit.optimizer_diagnostics["inactive_dual_minimum"], 0.5,
            places=10,
        )

    def test_duplicate_references_are_reported_nonidentifiable(self):
        wavelength = np.linspace(0.0, 1.0, 30)
        base = 1.0 + 0.1 * wavelength
        references = np.vstack([
            base, base, 1.0 - 0.1 * wavelength,
            1.0 + 0.2 * wavelength, 1.0 - 0.2 * wavelength,
        ])
        lpf = core.LPFunctionLinComb(
            wavelength, base, np.full(30, 0.01), references,
            np.ones_like(references),
        )
        fit = core.FitLinCombSpec(
            lpf, teffs=[3200, 4800, 3500, 3700, 3900], verbose=0
        )
        fit.minimize_PyDE(mcmc=False, verbose=0)
        diagnostics = fit.optimizer_diagnostics
        self.assertTrue(diagnostics["success"])
        self.assertTrue(diagnostics["nonunique_optimum"])
        self.assertGreater(diagnostics["distinct_tied_weight_vectors"], 1)
        self.assertGreater(
            diagnostics["tied_stellar_parameter_spreads"]["teff_k"], 0.0
        )

    def test_active_set_random_stress_matches_slsqp(self):
        rng = np.random.default_rng(90210)
        for _ in range(40):
            references = 1.0 + rng.normal(0.0, 0.03, (5, 35))
            true_weights = rng.dirichlet(np.full(5, 0.7))
            target = true_weights @ references + rng.normal(0.0, 0.002, 35)
            error = np.exp(rng.normal(-4.0, 0.3, 35))
            lpf = core.LPFunctionLinComb(
                np.arange(35.0), target, error, references,
                np.ones_like(references),
            )
            fit = core.FitLinCombSpec(lpf, verbose=0)
            fit.minimize_PyDE(mcmc=False, verbose=0)
            weights = np.asarray(lpf.get_pv_all(fit.min_pv))
            reference = scipy.optimize.minimize(
                lambda candidate: np.sum(
                    ((target - candidate @ references) / error) ** 2
                ),
                np.full(5, 0.2), method="SLSQP",
                bounds=[(0.0, 1.0)] * 5,
                constraints={"type": "eq", "fun": lambda value: value.sum() - 1.0},
                options={"ftol": 1e-12, "maxiter": 1000},
            )
            self.assertTrue(reference.success)
            ours = np.sum(((target - weights @ references) / error) ** 2)
            self.assertLessEqual(ours, reference.fun + 1e-7)

    def test_pairwise_profile_uses_one_frozen_mask_at_all_vsini(self):
        wavelength = np.arange(8500.0, 8510.0, 0.01)
        rng = np.random.default_rng(4)
        reference = np.ones_like(wavelength)
        target = reference + rng.normal(0.0, 0.01, wavelength.size)
        error = np.full_like(wavelength, 0.01)
        objective = core.PairwiseRSSFunctionVsiniPolynomial(
            wavelength, target, error, reference, error, maxvsini=30.0,
        )
        counts = []
        for value in (0.0, 15.0, 30.0):
            _, rss, diagnostics = objective.solve_coefficients(value)
            self.assertTrue(np.isfinite(rss))
            counts.append(diagnostics["valid_pixel_count"])
        self.assertEqual(counts, [objective.comparison_pixel_count] * 3)

    def test_scaled_chebyshev_is_stable_on_physical_wavelengths(self):
        wavelength = np.linspace(8446.3, 8565.3, 151)
        coordinate = core._scaled_chebyshev_coordinate(wavelength)
        coefficients = np.array([1.0, 0.03, -0.01, 0.004, 0.0, -0.002])
        reference = np.ones_like(wavelength)
        target = np.polynomial.chebyshev.chebval(coordinate, coefficients)
        objective = core.PairwiseRSSFunctionVsiniPolynomial(
            wavelength, target, np.ones_like(target), reference,
            np.ones_like(reference), maxvsini=30.0, vsini=4.0,
            vsini_window=0.0,
        )
        fit = core.FitTargetRefStarVsiniPolynomial(objective)
        with mock.patch.object(
                core.rotbroad_help, "broaden",
                side_effect=lambda w, f, vsini, u1: f), \
             mock.patch.object(
                core.rotbroad_help, "broaden_variance",
                side_effect=lambda w, variance, vsini, u1: variance):
            fit.minimize_AMOEBA(verbose=0)
            model = objective.compute_model(fit.min_pv)
        self.assertTrue(np.all(np.isfinite(fit.min_pv)))
        self.assertLess(np.max(np.abs(model - target)), 1e-5)
        self.assertLess(np.max(np.abs(fit.min_pv[1:])), 2.0)

    def test_stage_two_preparation_matches_pairwise_model(self):
        wavelength = np.linspace(8564.5, 8685.7, 41)
        coefficients = np.array([1.0, 0.02, -0.01, 0.0, 0.0, 0.0])
        flux = np.linspace(0.9, 1.1, wavelength.size)
        error = np.full_like(flux, 0.02)

        class Spectrum:
            def __init__(self, name, values):
                self.object = name
                self.rv = 0.0
                self.rv_source = "drp"
                self.provenance = {"rv_source": "drp"}
                self.values = values

            def deblaze(self):
                return None

            def resample_order(self, w, **kwargs):
                self.kwargs = kwargs
                return self.values.copy(), error.copy()

        target_spectrum = Spectrum("target", flux)
        reference = Spectrum("reference", flux)
        reference_list = type("List", (), {"splist": [reference]})()
        objective = core.PairwiseRSSFunctionVsiniPolynomial(
            wavelength, flux, error, flux, error, maxvsini=30.0,
            vsini=6.0, vsini_window=0.0,
        )
        parameters = np.r_[6.0, coefficients]
        with mock.patch.object(
                core.rotbroad_help, "broaden",
                side_effect=lambda w, f, vsini, u1: f), \
             mock.patch.object(
                core.rotbroad_help, "broaden_variance",
                side_effect=lambda w, variance, vsini, u1: variance):
            expected = objective.compute_model(parameters)
            prepared = core.get_data_ready(
                target_spectrum, reference_list, wavelength,
                np.linspace(-1, 1, 3), polyvals=[coefficients],
                vsinis=[6.0], order=102, verbose=0,
            )
        self.assertTrue(np.allclose(prepared[2][0], expected))
        continuum = np.polynomial.chebyshev.chebval(
            core._scaled_chebyshev_coordinate(wavelength), coefficients
        )
        self.assertTrue(np.allclose(prepared[3][0], error * np.abs(continuum)))
        self.assertNotIn("p", reference.kwargs)
        self.assertNotIn("vsini", reference.kwargs)

    def test_pairwise_metric_is_exact_unweighted_rss(self):
        wavelength = np.arange(8500.0, 8510.0, 0.05)
        target = np.linspace(1.0, 1.08, wavelength.size)
        reference = np.ones(wavelength.size)
        objective = core.PairwiseRSSFunctionVsiniPolynomial(
            wavelength, target, np.full(wavelength.size, 1000.0), reference,
            np.full(wavelength.size, 0.001), maxvsini=30.0, vsini=3.0,
            vsini_window=0.0,
        )
        parameters = np.r_[3.0, [1.0, 0, 0, 0, 0, 0]]
        with mock.patch.object(
            core.rotbroad_help, "broaden",
            side_effect=lambda w, f, vsini, u1: f):
            measured = objective(parameters)
        mask = objective.comparison_mask
        self.assertAlmostEqual(
            measured, np.sum((target[mask] - reference[mask]) ** 2)
        )

    def test_dynamic_rotational_kernel_support_and_flux_conservation(self):
        wavelength = np.arange(8564.5, 8685.7, 0.01)
        dv = (
            np.median(np.diff(wavelength)) / np.median(wavelength)
            * 2.99792458e5
        )
        constant = np.full(wavelength.size, 3.25)
        for vsini in (25.0, 50.0, 125.0):
            half_width = int(np.ceil(vsini / dv)) + 1
            velocity, kernel = rotbroad_help.rot(
                2 * half_width + 1, dv, vsini, u1=0.3
            )
            self.assertLess(velocity[0], -vsini)
            self.assertGreater(velocity[-1], vsini)
            self.assertTrue(np.all(kernel[np.abs(velocity) > vsini] == 0.0))
            self.assertAlmostEqual(float(kernel.sum()), 1.0, places=14)
            broadened = rotbroad_help.broaden(
                wavelength, constant, vsini, u1=0.3
            )
            finite = np.isfinite(broadened)
            self.assertGreater(np.count_nonzero(finite), wavelength.size // 2)
            self.assertTrue(np.allclose(
                broadened[finite], constant[finite], atol=1e-12
            ))

    def test_zero_rotation_is_exact_identity_and_grid_is_validated(self):
        wavelength = np.arange(5173.8, 5217.4, 0.01)
        flux = np.sin(wavelength)
        broadened = rotbroad_help.broaden(wavelength, flux, 0.0)
        self.assertTrue(np.array_equal(broadened, flux))
        invalid = wavelength.copy()
        invalid[10] = invalid[9]
        with self.assertRaises(ValueError):
            rotbroad_help.broaden(invalid, flux, 25.0)

    def test_pairwise_profile_finds_zero_and_upper_endpoint(self):
        wavelength = np.arange(8564.5, 8569.5, 0.01)
        reference = 1.0 - 0.4 * np.exp(
            -0.5 * ((wavelength - wavelength.mean()) / 0.08) ** 2
        )
        error = np.full_like(reference, 0.01)
        for expected in (0.0, 30.0):
            target = rotbroad_help.broaden(wavelength, reference, expected, u1=0.3)
            objective = core.PairwiseRSSFunctionVsiniPolynomial(
                wavelength, target, error, reference, error,
                maxvsini=30.0,
            )
            fit = core.FitTargetRefStarVsiniPolynomial(objective)
            fit.minimize_AMOEBA(verbose=0)
            self.assertAlmostEqual(fit.min_pv[0], expected, delta=0.05)
            self.assertTrue(fit.optimizer_diagnostics["success"])

    def test_high_vsini_gap_mask_and_variance_propagation_are_finite(self):
        wavelength = np.arange(8564.5, 8594.5, 0.01)
        reference = 1.0 - 0.2 * np.exp(
            -0.5 * ((wavelength - 8568.0) / 0.1) ** 2
        )
        reference[1400:1430] = np.nan
        error = np.full_like(reference, 0.02)
        error[1400:1430] = np.nan
        flux, propagated_error = core._prepare_reference_model(
            wavelength, reference, error, [1, 0, 0, 0, 0, 0], 100.0
        )
        valid = np.isfinite(flux) & np.isfinite(propagated_error)
        self.assertGreater(np.count_nonzero(valid), 100)
        self.assertTrue(np.all(propagated_error[valid] > 0))
        self.assertTrue(np.all(~valid[1400:1430]))

    def test_rotbroad_compatibility_module_reexports_neidspec(self):
        from neidspec import rotbroad_help as canonical
        self.assertIs(rotbroad_help.broaden, canonical.broaden)
        self.assertIs(rotbroad_help.broaden_variance, canonical.broaden_variance)


if __name__ == "__main__":
    unittest.main()
