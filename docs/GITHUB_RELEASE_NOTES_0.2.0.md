# NEIDSpecMatch 0.2.0

NEIDSpecMatch 0.2.0 is the first release intended for reproducible science use
with current NEID DRP 1.5 Level-2 spectra. It replaces implicit legacy
assumptions with explicit validation and provenance checks. Users should read
the validation scope below before adopting publication-facing parameters.

## Highlights

- **Current NEID reductions.** The default empirical library contains the
  same 78 labeled stars reprocessed in the DRP 1.5 series (77 v1.5.3 and one
  v1.5.2). Targets and references must share a DRP major/minor series. The
  matcher now uses the blaze arrays embedded in each Level-2 product; the old
  external DRP-1.3 instrument-response file is no longer required.
- **Correct order identity.** An `order` is always the zero-based row of the
  122-row NEID spectral array. Loading a restricted row range no longer shifts
  requested order numbers.
- **Scientifically safer metallicities.** The HPF order-5 detrending
  coefficients formerly applied to NEID spectra have been removed. Raw
  empirical-library `[Fe/H]` is now the default. No metallicity calibration is
  bundled because the tested cross-fitted corrections did not improve held-out
  performance and do not constitute independent validation.
- **Explicit validation scope.** Ordinary free-`v sin i` fits are validated
  for array-row order 55 in the hot population and orders 101, 102, and 103 in
  the cool population, using a 4500 K boundary and the exact DRP-1.5 library,
  runtime, and source fingerprints recorded in the release bundle. Order 102
  is the conservative default for cool/M-dwarf targets. Validation errors are
  empirical population-level predictive metrics, not automatic per-target
  Gaussian uncertainties.
- **Deterministic fitting.** Pairwise continuum/rotation fitting and the
  five-reference composite solve are deterministic, use a common pixel mask,
  and report optimizer, boundary, and identifiability diagnostics. The
  five-reference selection rule remains the established best-five-by-pairwise-
  fit procedure.
- **RV and provenance fixes.** McDonald Observatory defaults are removed.
  NEID barycentric metadata and structurally screened DRP RVs are read from the
  product, and custom CCF masks require a recorded air/vacuum wavelength
  policy. Results bind the library, spectra, DRP/blaze state, numerical stack,
  matcher source, and `neidspec` source with cryptographic receipts.
- **Safer distribution.** Credentials and large data are absent from the
  package. The external library is verified by pinned archive size/MD5, exact
  ZIP layout, release-owned per-file SHA-256 allowlisting, FITS schema, and
  data-quality checks. Imports are quiet and do not download data.

## Important limitations

- The release validation applies only to the ordinary free-`v sin i`
  estimator and the exact order/population/library combinations above.
- Fixed or bounded external-`v sin i` fits execute consistently, but their
  atmospheric parameters remain exploratory pending validation on independent
  labeled broad-lined standards.
- The fitted rotational broadening is not independently calibrated as a
  physical publication-quality `v sin i` measurement.
- A successful optimizer is not by itself evidence of publication validity;
  use the matching validation summary and inspect support/quality warnings.

## Compatibility and reproducibility

Python 3.10, 3.11, and 3.12 are tested. For exact agreement with the release
validation, install with the supplied `constraints/validated-runtime.txt`.
The full validation report, raw cross-validation predictions, fold receipts,
summary tables, and research-note-style figure are linked from the release.

The DRP-1.5 library is archived under DOI
[`10.5281/zenodo.22262405`](https://doi.org/10.5281/zenodo.22262405). Do not
publish or advertise that dataset as available until NEID redistribution terms
have been confirmed and the Zenodo draft has been released.

See [Scientific validation and publication use](https://github.com/TeHanHunter/neidspecmatch/blob/v0.2.0/docs/SCIENCE_VALIDATION.md),
the [validation report](https://github.com/TeHanHunter/neidspecmatch/blob/v0.2.0/docs/VALIDATION_REPORT_0.2.0.md),
and the [changelog](https://github.com/TeHanHunter/neidspecmatch/blob/v0.2.0/CHANGELOG.md)
for details.
