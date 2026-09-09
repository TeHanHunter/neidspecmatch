# NEIDSpecMatch

[![DOI](https://zenodo.org/badge/769285239.svg)](https://doi.org/10.5281/zenodo.14991480)

NEIDSpecMatch estimates stellar parameters by comparing a high-resolution NEID
spectrum with an empirical library. It is based on
[HPFSpecMatch](https://gummiks.github.io/hpfspecmatch/).

Version 0.2.0 is a science release. Publication use of the ordinary
free-`v sin i` estimator requires matching the target and library DRP series,
choosing orders appropriate for the target, and using the cross-validation
product for exactly that library/order/population. The package does not turn an
unvalidated fit into a publication-quality result automatically.

## Installation

Use Python 3.10, 3.11, or 3.12 in a clean environment. For results intended
to match the 0.2.0 validation product, install the exact numerical stack used
for that validation:

```bash
python -m pip install \
  --constraint https://raw.githubusercontent.com/TeHanHunter/neidspecmatch/v0.2.0/constraints/validated-runtime.txt \
  'neidspecmatch==0.2.0'
```

A normal `python -m pip install 'neidspecmatch==0.2.0'` remains supported for
exploration. Every result records its runtime versions, and a different
NumPy/SciPy/Astropy/pandas stack will deliberately not match the published
cross-validation receipt.

Optional features are isolated from the core install:

```bash
python -m pip install 'neidspecmatch[archive]==0.2.0'  # archive/catalog helpers
python -m pip install 'neidspecmatch[dev]==0.2.0'      # tests and builds
```

For each reference, the pairwise stage solves six continuum coefficients by
linear least squares and profiles a deterministic one-dimensional `v sin i`
grid, refining every sampled local-minimum basin. The five-reference composite
stage enumerates all simplex faces and solves constrained weighted least
squares in a stable null-space basis. Neither stage uses PyTransit/PyDE,
Powell, or stochastic differential evolution.

## Library location and download

The historical 78-star v1 library is available from
[Zenodo (DOI 10.5281/zenodo.14947454)](https://doi.org/10.5281/zenodo.14947454).
Although that uploader-supplied record is labeled CC BY 4.0, NEIDSpecMatch does
not treat the label as evidence of authority to relicense the underlying NEID
Level-2 FITS. The replacement DRP-1.5 library is published with the 0.2.0
release at
[Zenodo (DOI 10.5281/zenodo.22262405)](https://doi.org/10.5281/zenodo.22262405).
Its 78 Level-2 FITS files are unmodified public NEID Archive products. The MIT
license applies to the NEIDSpecMatch software only and does not relicense those
archive products; see `DATA_LICENSES.md`.
The archive remains external because it is multi-gigabyte and because its
exact manifest and reduction provenance must stay explicit. External download
does not, by itself, resolve data rights. NEIDSpecMatch never downloads it at
import time and never writes into site-packages. Select its versioned directory
in one of three ways:

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

neidspecmatch.get_library()  # explicit 6.1 GiB download; about 8.6 GiB extracted
neidspecmatch.validate_library()
```

## Fit a target

The installed command contains no machine-specific paths. Version 0.2 is for
NEID Level-2 high-resolution (HR) products only. Order 102 is the conservative
default for the cool/M-dwarf use case:

```bash
neidspecmatch-fit TARGET_L2.fits 'TIC 123456' results \
  --orders 102 --library-path /data/20260831_specmatch_neid_drp15 -v
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

The 0.2.0 release bundle reports the full, cool, and hot populations for each
of orders 55, 101, 102, and 103, with the population boundary at 4500 K.
Order 55 is recommended for hotter stars, while order 102 remains the
conservative default for cool/M-dwarf targets. Other
valid order indices remain available for exploratory work, but are marked
unvalidated unless a current, exactly matching validation product is supplied
for the ordinary free-`v sin i` mode. Publication validation is population-
and order-specific: run one order per publication-facing science-fit
invocation. The cross-validation command may process multiple independent
orders in one invocation. Do not combine results from multiple orders and
describe that estimator as validated without a separate combined-order
validation product.

The complete four-order validation products, including raw predictions, fold
receipts, population summaries, provenance, and the validation figure, are
available from the
[v0.2.0 GitHub release](https://github.com/TeHanHunter/neidspecmatch/releases/tag/v0.2.0)
as `neidspecmatch-0.2.0-four-order-validation.zip`.

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

An independent temperature is not required to run the spectral fit and is not
an input to its atmospheric-parameter optimization. It is used afterward only
to establish membership when the applicable population error metrics are
selected. For this release, use cool for an independent temperature below
4500 K and hot at or above 4500 K; no reference star lies exactly on the
boundary. Generic population bounds are stored inclusively, so a science
target exactly at 4500 K must explicitly select the release convention `hot`.

In the historical DRP-1.3/20250226-library audit, an all-star calibration
increased the cross-fitted [Fe/H] scatter in every one of the 49 tested orders.
For example, order 102 increased from 0.128 dex raw to 0.154 dex calibrated.
Some orders were also nearly singular. Consequently, version 0.2 bundles no
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
neidspecmatch-crossval --orders 55 101 102 103 \
  --population all --population 'cool::4500' --population 'hot:4500' \
  --library-path /data/20260831_specmatch_neid_drp15 \
  --library-id 20260831_specmatch_neid_drp15 --output validation -v
```

This command performs ordinary spectral leave-one-out validation with free
pairwise `v sin i`; it validates only that estimator on the unmodified library
population. It does not validate a science fit made with `--vsini` or
`--vsini-window`.

The 0.2.0 dataset preserves the completed products under their original
0.2.0b1 fingerprints and separately supplies a numerical-continuity record to
the final 0.2.0 source. It does not rewrite historical receipts as exact final-
source products. Consequently, the final package's stricter automatic
runtime-attached validation flag remains conservative for those archived
summaries; use the release report when quoting their population-level errors.

This writes a raw result CSV, an atomic fold checkpoint, and a machine-readable
per-order summary. The installed fail-closed figure command first verifies all
four raw/checkpoint/summary products and the exact library manifest, then
writes the research-note-style figure, a combined summary, and a provenance
receipt:

```bash
neidspecmatch-crossval-figure \
  --validation-dir validation \
  --library-manifest /data/20260831_specmatch_neid_drp15/library_manifest.json
```

For a science fit, supply the matching per-order summary and the independently
established population membership. For example:

```bash
neidspecmatch-fit TARGET_L2.fits 'TIC 123456' results --orders 102 \
  --validation-summary validation/o102_crossval/crossvalidation_summary_o102.csv \
  --validation-population cool --validation-population-teff 3800
```

Diagnostic calibration cross-fitting is opt-in and remains explicitly
non-publication-valid:

```bash
neidspecmatch-crossval --orders 102 --fit-feh-calibration \
  --population 'cool::4500' \
  --library-path /data/20260831_specmatch_neid_drp15 \
  --library-id 20260831_specmatch_neid_drp15 --output validation
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

Never relabel an archived result with a newer library identifier.

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

The uncertainty equations reported for the 0.1 analysis do not validate the
changed 0.2 L2-blaze/RV/continuum pipeline. The 0.2.0 release bundle supplies
78-fold evidence for orders 55, 101, 102, and 103 in both named 4500 K
populations, as well as the full library. An ordinary fit remains marked exploratory unless its summary
matches the exact result schema, pipeline fingerprint, numerical-library
versions, `neidspec` source, deep-verified library manifest, DRP series, blaze
source, order, population, and (when used) calibration artifact hash. Fixed
and bounded fits remain exploratory even when an ordinary leave-one-out
summary is available.

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
