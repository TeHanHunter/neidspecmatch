# NEIDSpecMatch

[![DOI](https://zenodo.org/badge/769285239.svg)](https://doi.org/10.5281/zenodo.14991480)

NEIDSpecMatch estimates stellar parameters by comparing a high-resolution NEID
spectrum with an empirical library. It is based on
[HPFSpecMatch](https://gummiks.github.io/hpfspecmatch/).

Version 0.2.0b1 is a beta science release. Publication use of the ordinary
free-`v sin i` estimator requires matching the target and library DRP series,
choosing orders appropriate for the target, and using the cross-validation
product for exactly that library/order/population. The package does not turn an
unvalidated fit into a publication-quality result automatically.

## Installation

Use Python 3.10, 3.11, or 3.12 in a clean environment:

```bash
python -m pip install 'neidspecmatch==0.2.0b1'
```

Optional features are isolated from the core install:

```bash
python -m pip install 'neidspecmatch[archive]==0.2.0b1'  # archive/catalog helpers
python -m pip install 'neidspecmatch[dev]==0.2.0b1'      # tests and builds
```

For each reference, the pairwise stage solves six continuum coefficients by
linear least squares and profiles a deterministic one-dimensional `v sin i`
grid, refining every sampled local-minimum basin. The five-reference composite
stage enumerates all simplex faces and solves constrained weighted least
squares in a stable null-space basis. Neither stage uses PyTransit/PyDE,
Powell, or stochastic differential evolution.

## Library location and download

The 78-star v1 library is archived under CC BY 4.0 on
[Zenodo (DOI 10.5281/zenodo.14947454)](https://doi.org/10.5281/zenodo.14947454)
by Te Han. Reuse requires attribution. The archive remains external because it
is 7.3 GB and because its exact manifest/reduction provenance must stay
explicit—not because its reuse rights are unresolved. NEIDSpecMatch never
downloads it at import time and never writes into site-packages. Select its
versioned directory in one of three ways:

1. pass `library_path=...`;
2. set `NEIDSPECMATCH_LIBRARY` to that directory; or
3. use the platform-specific user-data default reported by
   `neidspecmatch.resolve_library_path()`.

An explicit download pins the HTTPS Zenodo record URL and byte size, verifies
Zenodo's published MD5, and authenticates every catalog/FITS file against the
release-owned SHA-256 allowlist shipped in the wheel. It rejects unsafe ZIP
members, extra files, and resource-exhaustion archives, validates the catalog
plus NEID L2/HR/FITS/DQ schema for all 78 spectra, and writes a per-file
size/SHA-256 manifest with the allowlist receipt. The source-published MD5 is
an integrity check only, not cryptographic authentication:

```python
import neidspecmatch

neidspecmatch.get_library()  # explicit network and 7.3 GB disk operation
neidspecmatch.validate_library()
```

## Fit a target

The installed command contains no machine-specific paths. Version 0.2 is for
NEID Level-2 high-resolution (HR) products only. Order 102 is the conservative
default for the cool/M-dwarf use case:

```bash
neidspecmatch-fit TARGET_L2.fits 'TIC 123456' results \
  --orders 102 --library-path /data/20250226_specmatch_nir -v
```

Plots are opt-in with `--plot`. `-v` reports per-order progress; omit it for a
quiet batch run. A supplied rotational velocity is fixed exactly by default:

```bash
neidspecmatch-fit TARGET_L2.fits 'TIC 123456' results \
  --vsini 18.2 --vsini-window 0
```

Every reference-star comparison receives the same immutable constraint. A
nonzero window fits within `vsini ± window`; it is not recentered on the
previous reference-star result.

The default RV is the structurally screened NEID DRP value. If a reviewed
external line mask is needed instead, select it explicitly and declare its
input wavelength medium:

```bash
neidspecmatch-fit TARGET_L2.fits 'TIC 123456' results \
  --rv-source custom_ccf --ccf-mask /data/reviewed-mask.txt \
  --ccf-mask-medium vacuum
```

NEID wavelengths and all CCF coordinates are vacuum wavelengths. The exact
bundled ESPRESSO M3 mask has a checksum-bound air-to-vacuum policy, so its
medium option may be omitted; every unregistered mask requires `air` or
`vacuum`. The mask digest and applied conversion policy are retained in the
result provenance. This numerical path does not by itself validate a mask's
astrophysical suitability for a target.

Fixed and bounded constraints are supported execution modes, not validated
atmospheric-parameter modes in this release. The ordinary leave-one-out run
uses unmodified, predominantly narrow-lined library stars as targets and fits
`v sin i` freely. Applying the same nonzero constraint to those unbroadened
held-out targets would test a deliberately mismatched problem, not the
broad-lined science case. Consequently, `Teff`, `[Fe/H]`, and `log(g)` from a
fixed or bounded run remain exploratory pending validation on independent,
labeled broad-lined standards representative of the science targets.

Order 55 may be selected explicitly for the previously defined hot-star
population. Other valid order indices remain available for exploratory work,
but are marked unvalidated unless a current, exactly matching validation
product is supplied for the ordinary free-`v sin i` mode. Publication
validation is population- and order-specific: run one order per invocation.
Do not median order 55 and 102 and describe that combination as validated
without a new combined-order validation product.

Normal input is an untouched NEID L2 product and the matcher applies its blaze
correction. Use `--input-is-deblazed` only when HDU 1 is already corrected;
the flag means “bypass L2 blaze correction,” not “please deblaze.”

## [Fe/H] calibration: safe default and provenance

The safe default is the raw empirical-library weighted [Fe/H]. Version 0.1
silently applied HPF order-5 coefficients to NEID results; those coefficients
have been removed. `calibrate_feh=True` without a validated artifact now raises
an error.

A calibration artifact is:

- specific to a NEID order, named library, and explicit population;
- fitted from `recovered - true` residuals;
- cross-fitted so each row's calibration excludes that row;
- accompanied by raw and calibrated scatter, bias, sample size, source file,
  residual definition, and validation limitations.

No numerical calibration coefficients are embedded in the source release.
They are generated only from supplied cross-validation evidence. Population
cuts use `teff_true` in validation, so bounded populations must be selected
explicitly at science time; recovered Teff is not silently used as an
unvalidated classifier.

For the archived 20250226 library results, an all-star calibration increased
the cross-fitted [Fe/H] scatter in every one of the 49 tested orders. For
example, order 102 increased from 0.128 dex raw to 0.154 dex calibrated. Some
orders were also nearly singular. Consequently, version 0.2 bundles no
calibration artifact and recommends the raw value. Diagnostic artifact
generation remains available for reviewed order/population analyses. Inverse
calibrations that are non-monotonic, have `1 + slope < 0.1` in any cross-fit
fold, or worsen scatter, RMSE, or absolute bias are rejected.

The archived spectral leave-one-out predictions are correlated: another
star's spectral fit can use the calibration-held-out star as a reference.
Calibration cross-fitting therefore avoids in-sample detrending but is not a
strict nested/leave-two-out validation. Because of that residual dependence,
0.2 marks every calibration produced by this workflow unvalidated and refuses
to apply it. A future applicable artifact requires strict nested/leave-two-out
or independent validation and exact manifest, DRP, blaze, matcher, and
`neidspec` source-hash agreement. Legacy pickle output is opt-in and must be
treated as trusted-input-only because loading a pickle can execute code.

## Cross-validation

Cross-validation defaults to no plots and raw [Fe/H]:

```bash
neidspecmatch-crossval --orders 102 \
  --library-path /data/20250226_specmatch_nir \
  --library-id 20250226_specmatch_nir --output validation -v
```

This command performs ordinary spectral leave-one-out validation with free
pairwise `v sin i`; it validates only that estimator on the unmodified library
population. It does not validate a science fit made with `--vsini` or
`--vsini-window`.

This writes a raw result CSV, an atomic fold checkpoint, and a machine-readable
per-order summary. Diagnostic calibration cross-fitting is opt-in and remains
explicitly non-publication-valid:

```bash
neidspecmatch-crossval --orders 102 --fit-feh-calibration \
  --population 'cool::4000' \
  --library-path /data/20250226_specmatch_nir \
  --library-id 20250226_specmatch_nir --output validation
```

Existing archived CV results can be summarized without rerunning thousands of
spectral fits, but the product is archival/exploratory and cannot create an
applicable calibration because it lacks current fold evidence:

```bash
neidspecmatch-crossval \
  --summarize-existing crossvalidation_results_o102.csv \
  --library-id 20250226_specmatch_nir \
  --output validation
```

For multiple orders, `crossvalidation_summary.csv` contains one row per
order/population. The residual RMSE is a predictive-error score that includes
bias; sample scatter (`ddof=1`), bias, finite count, and coverage are retained
alongside it. RMSE is not automatically a Gaussian or 68%-coverage `1 sigma`
uncertainty for an individual target.

Synthetically broadening each held-out spectrum, with consistent variance and
mask propagation, can be a useful stress test of numerical behavior and the
loss of information with line broadening. It is not publication validation:
the injected spectrum and matcher share the same kernel assumptions, and the
experiment does not reproduce the activity, spots, magnetic broadening,
intrinsic rotation, or noise properties of observed broad-lined targets.

Long runs checkpoint atomically after each successful fold and resume only
when the order, library/manifest, DRP set, blaze source, random seed, software
versions, and pipeline fingerprint match. Use `--no-resume` for a fresh run.
Fold pairwise and composite optimizer/identifiability diagnostics are stored;
a summary containing any failed, non-finite, or non-identifiable fit cannot
validate a science result. Calibration
cross-fitting is not strict nested/leave-two-out validation, and oversampling
does not create independent spectral pixels.

The uncertainty equations reported for the 0.1 analysis do not automatically
validate the changed 0.2 L2-blaze/RV/continuum pipeline. An ordinary fit with
free pairwise `v sin i` remains marked exploratory until a current summary
matches the exact result schema, pipeline fingerprint, neidspec version,
deep-verified library manifest, DRP series, blaze source, order, population,
and (when used) calibration artifact hash. Fixed and bounded fits remain
exploratory even when an ordinary leave-one-out summary is available.

## DRP and barycentric metadata

Do not treat successful execution as validation across NEID DRP minor
versions. Patch versions within one minor series are compatible under the NEID
DRP policy, and their exact values remain recorded in provenance. Use targets
and a library processed within one minor series and rerun cross-validation
after changing the library reduction.

NEID L2 spectra already provide per-order `SSBJDnnn` and `SSBRVnnn`. The
`neidspec` dependency should consume those values directly. The legacy local
barycentric helper no longer defaults to McDonald Observatory and refuses to
guess an observatory.

## NEID Archive credentials

Never put credentials in scripts or notebooks. The optional archive helper
uses only:

```bash
export NEID_ARCHIVE_USER='your-user'
export NEID_ARCHIVE_PASSWORD='your-password'
python -m neidspecmatch.neid_archive /private/cookie-directory
```

The cookie is permission-restricted. Keep `.env`, cookies, and debug logs out
of version control. PyPI releases 0.1.1, 0.1.2, 0.1.3, and 0.1.4 each contained
a now-reset hard-coded archive credential. Those releases should be yanked
after 0.2.0 is available. The credential is not reproduced here; cached copies
must be treated as compromised even though the password has been reset.

## Reference

```bibtex
@ARTICLE{2025RNAAS...9...63H,
  author = {{Han}, Te and {Robertson}, Paul and {Ca{\~n}as}, Caleb I. and others},
  title = {NEIDSpecMatch: Stellar Parameter Estimation with NEID Spectra Using an Empirical Library},
  journal = {Research Notes of the American Astronomical Society},
  year = 2025,
  volume = 9,
  number = 3,
  pages = 63,
  doi = {10.3847/2515-5172/adc264}
}
```
