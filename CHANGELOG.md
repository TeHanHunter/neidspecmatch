# Changelog

## 0.2.0 (2026-09-03)

This release corrects the software and release pipeline and adds a current
DRP-1.5 empirical reference library with order- and population-specific
leave-one-out validation for the ordinary free-`v sin i` estimator. Orders
55, 101, 102, and 103 each report the full library and both populations split
at 4500 K. Order 55 is recommended for hotter stars and order 102 is the
conservative default for cool/M-dwarf targets. Atmospheric parameters retain raw [Fe/H]. Fixed and bounded
fits still require independent labeled broad-lined validation.

- Remove the silently applied HPF order-5 metallicity coefficients. Raw
  empirical-library [Fe/H] is now the default; calibration requires a
  validated, order/library/population-specific artifact.
- Audit all 49 historical DRP-1.3/20250226-library order results: every
  cross-fitted calibration increased [Fe/H] scatter (order 102: 0.128 to
  0.154 dex), so no coefficient artifact is bundled or recommended.
  Near-singular inverse calibrations are marked unvalidated and cannot be
  applied.
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
- Bind result/CV receipts to the estimator and `neidspec` source fingerprints,
  exact numerical-library versions, raw/checkpoint hashes, DRP/blaze state,
  and DQ composition. Presentation-only figure code is separately hashed and
  does not invalidate numerical validation evidence.
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
- Add the installed `neidspecmatch-crossval-figure` command, which fails closed
  unless all four raw/checkpoint/summary products match the supplied library
  manifest and writes an independently hashed provenance receipt.
