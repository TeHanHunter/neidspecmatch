# Changelog

## 0.2.0b1 (2026-08-31)

This beta validates the corrected software and release pipeline, but it does
not yet provide publication-valid stellar-parameter uncertainties. For
ordinary free-`v sin i` fits, those require the planned v1.5.x
reference-library cross-validation before the stable 0.2.0 release. Fixed and
bounded fits also require independent labeled broad-lined validation.

- Remove the silently applied HPF order-5 metallicity coefficients. Raw
  empirical-library [Fe/H] is now the default; calibration requires a
  validated, order/library/population-specific artifact.
- Audit all 49 archived full-library order results: every cross-fitted
  calibration increased [Fe/H] scatter (order 102: 0.128 to 0.154 dex), so no
  coefficient artifact is bundled or recommended. Near-singular inverse
  calibrations are marked unvalidated and cannot be applied.
- Cross-fit calibration coefficients and emit per-order error-metric,
  calibration, residual, and provenance products. Record the correlation
  limitation of archived spectral leave-one-out predictions.
- Make cross-validation plotting and detailed logging opt-in.
- Add exact or bounded external-vsini constraints without carrying one
  reference star's fit into the next reference-star prior. These modes execute
  but their atmospheric parameters remain exploratory pending independent
  labeled broad-lined validation; synthetic broadening is a stress test only.
- Replace local/stochastic optimizers with profiled continuum least squares
  plus deterministic one-dimensional `v sin i` search and an enumerated
  simplex active-set solve. Freeze one common pairwise pixel mask and reject
  non-identifiable composite solutions.
- Move the large library out of site-packages, add an environment/explicit path
  API and custom reprocessed-library manifests, and verify pinned source size,
  archive checksum, exact ZIP layout, NEID L2 schema, and every payload file
  against a release-owned SHA-256 allowlist.
- Bind result/CV receipts to full matcher and `neidspec` source fingerprints,
  numerical versions, raw/checkpoint hashes, DRP/blaze state, and DQ
  composition.
- Treat calibration cross-fitting as diagnostic only: correlated spectral-LOO
  inputs are not strict nested validation, so 0.2 artifacts cannot be applied.
- Remove source-code archive credentials and require environment variables.
- Remove the McDonald Observatory default from legacy helpers.
- Bind custom-CCF wavelength media to provenance. The exact bundled ESPRESSO
  M3 mask is checksum-registered as air and converted to the NEID vacuum grid;
  unregistered masks require an explicit `air` or `vacuum` declaration. Saved
  results retain the path-free mask policy and per-order CCF diagnostics.
- Require `neidspec>=0.2.1`, the coordinated release that provides the
  checksum-bound CCF mask-medium API.
- Follow the NEID DRP compatibility policy: accept mixed patch releases within
  one major/minor series while recording every exact version, and reject mixed
  minor series by default.
- Use the science/sky blaze arrays embedded in each Level-2 product; the
  historical external v1.3 instrument-response file is no longer required.
- Remove unresolved-rights legacy FITS/catalog assets, output-filled notebooks,
  generated plots, obsolete runners, and HPF-branded documentation from the
  maintained release tree.
- Replace machine-specific runners with argument-driven command-line tools and
  add regression tests for calibration, quiet imports, archive credentials,
  cross-validation products, and vsini constraints.
