# Scientific validation and publication use

NEIDSpecMatch 0.2 is science software, but installation success and optimizer
success are not evidence that a result is publication quality.  An ordinary
free-`v sin i` result is eligible for atmospheric-parameter publication use
only when its JSON receipt matches a complete cross-validation product from
the same analysis state.

The required match is intentionally strict:

- exact NEIDSpecMatch and `neidspec` source fingerprints and result schema;
- exact deep-hashed empirical-library manifest and reference pool;
- one consistent NEID DRP version and the validated data-quality composition;
- the same Level-2 blaze treatment, NEID array-row order index, and wavelength
  grid;
- the same explicitly defined temperature population; and
- the exact [Fe/H] calibration artifact, if a calibration was applied.

If any item is absent or differs, the result remains useful for diagnosis or
exploration but is marked unvalidated.  Overrides such as mixed DRP versions
do not turn an exploratory result into a validated one.

The standard spectral leave-one-out procedure holds out each unmodified,
predominantly narrow-lined library spectrum and fits pairwise `v sin i` freely.
It therefore validates only the free-`v sin i` estimator on that library
population.  It does not validate atmospheric parameters from a fit in which
pairwise `v sin i` was externally fixed or bounded.

## What changed from 0.1

### Order identifiers

An NEIDSpecMatch `order` is a zero-based row in the 122-row NEID DRP arrays,
not a physical echelle number.  For example, array row 55 is echelle order
118 and array row 102 is echelle order 71.  The mapping is explicit in
`neidspec`; slicing an array no longer changes the meaning of the requested
order.

### Blaze and DRP compatibility

The default uses the science and sky blaze arrays in each Level-2 product.
The historical external instrument-response FITS file is neither required nor
silently selected.  An old response can be requested only as an explicit,
unvalidated comparison path.

Targets and reference spectra must have the same DRP major/minor series.  NEID
declares patch releases within a minor series compatible, so v1.5.2 and v1.5.3
may be validated together while v1.4.x and v1.5.x may not. Exact patch versions
remain recorded in provenance. A target reduced with a newer minor DRP should
not be downgraded merely to satisfy the code; the defensible solution is to
reprocess the reference library and targets consistently, then rerun
cross-validation.  The archived 78-star v1 library contains DRP v1.3.0
spectra, so it cannot validate results from later-minor-DRP targets without a
specific new study.

Level-2 high-resolution NEID science inputs are checked for instrument,
observation type, product level, observing mode, and DRP quality summary.
Fail/reject products are refused.  A warning product may be loaded for an
explicit audit, but a warning target is not publication validated.  Reference
warnings are acceptable only if the exact quality composition was represented
and bound into the matching empirical cross-validation product.

### RV and barycentric metadata

Per-order `SSBJDnnn` and `SSBRVnnn` values are read directly from the NEID
product.  The old McDonald Observatory default is gone.  The default stellar
RV is the product's weighted `CCFRVMOD` after numerical structure checks; this
means "structurally screened," not independently proven astrophysically
correct.  A failed or fast-rotator CCF requires an explicit, independently
justified RV.  The generic bundled M3 mask is not an automatic science
fallback.  If a custom mask is used, its wavelength medium is part of the RV
provenance: unregistered masks require an explicit `air` or `vacuum`
declaration, and air coordinates are converted to the NEID vacuum grid.  The
exact bundled M3 bytes have a SHA-256-bound air-to-vacuum policy; this fixes a
roughly 83 km/s coordinate-medium zero point but does not validate the mask's
astrophysical applicability.

### [Fe/H] calibration

The HPF order-5 detrending coefficients formerly applied by default were not
NEID order-specific validation products and have been removed.  Raw
empirical-library [Fe/H] is now the default. Applying a correction requires a
tamper-evident artifact built for the exact order, library, population, code,
and cross-validation run, plus strict nested or independent validation. The
0.2 calibration-cross-fit workflow is diagnostic and its artifacts are not
applicable.

Calibration cross-fitting prevents direct in-sample fitting of each row, but
the spectral leave-one-out predictions are still correlated: another
calibration row's spectral fit can include the held-out star as a reference.
This is not strict nested leave-two-out validation. Version 0.2 therefore
marks those products unvalidated even when their statistical gates improve.

### Error metrics and uncertainties

No uncertainty from the 0.1 analysis is silently reused. The exact current
leave-one-out run reports residual bias, sample scatter (`ddof=1`), predictive
RMSE, finite count, and coverage for each order and population. RMSE is useful
because it includes nonzero bias, but it is a population-level predictive-error
score. It is not automatically a Gaussian standard deviation, a 68%-coverage
interval, or a per-target `1 sigma` uncertainty. A failed, incomplete, or
non-identifiable fold invalidates the product.

The resampled grid is oversampled and neighboring samples share detector
pixels.  Rotational broadening creates further covariance.  The optimizer's
diagonal weighted loss is therefore a ranking/optimization loss, not a
calibrated chi-square distribution.  Publication uncertainties must remain
empirical.

### Rotation

The rotational kernel now has velocity-complete dynamic support and zero
rotation is an exact identity.  Every pairwise candidate is compared on one
fixed common pixel mask, preventing high trial velocities from appearing
better merely because more edge pixels were discarded.

The fitted broadening is still not independently calibrated as a physical
target `v sin i`, particularly because reference stars can have intrinsic
rotation.  Version 0.2 therefore keeps `v sin i` publication-unvalidated even
when the atmospheric parameters have a matching validation product.

Externally fixed and bounded pairwise-`v sin i` modes execute, but their
atmospheric parameters remain exploratory in this release. Simply forcing the
same nonzero value on an unbroadened held-out library star compares a sharp
target to broadened reference models and does not represent a broad-lined
science target. A synthetic-injection study must instead broaden the held-out
target, propagate its variance, erode masks around gaps and order edges, and
keep the original star excluded from the reference pool. Even that is a stress
test rather than publication validation: it reuses the matcher's kernel model
and does not reproduce the astrophysical and observational diversity of real
rapid rotators. Publication validation requires independent broad-lined stars
with independently determined atmospheric labels and coverage of the relevant
target population.

## Minimum release evidence

A publication-facing validation bundle should preserve:

1. the raw per-star cross-validation CSV and its SHA-256;
2. the per-fold optimizer and identifiability diagnostics;
3. the summary CSV and any [Fe/H] calibration JSON;
4. the deep library manifest, catalog, DRP and data-quality composition;
5. exact software/artifact hashes and the release tag;
6. representative failure-mode and real-spectrum integration results; and
7. for fixed or bounded broadening, independent labeled broad-lined standards.

Until that bundle exists for a selected analysis state, describe outputs as
exploratory rather than quoting the historical equations as current errors.
