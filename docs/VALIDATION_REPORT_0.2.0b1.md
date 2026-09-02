# NEIDSpecMatch 0.2.0b1 beta validation report

Status: **beta release candidate; maintainer approved**  
Validation date: 2026-08-31

This report records software, artifact, and archive-availability evidence. It
does not claim a publication-valid stellar-parameter scale. For the ordinary
free-`v sin i` estimator, that claim requires new cross-validation using the
current v1.5 reference reductions. Fixed and bounded pairwise-`v sin i` fits
remain exploratory pending independent labeled broad-lined validation;
synthetic broadening alone is a stress test, not publication validation.

## Automated regression matrix

The same 58-test suite was run with warnings treated as errors. Each
environment was installed from its Python-specific, fully hashed release lock
with `pip --require-hashes`; `pip check` reported no broken requirements.

| Python | Representative numerical environment | Result |
| --- | --- | --- |
| 3.10.14 | NumPy 1.26.4; SciPy 1.12.0 | 58 passed |
| 3.11.11 | NumPy 1.26.4; SciPy 1.12.0 | 58 passed |
| 3.12.2 | NumPy 1.26.4; SciPy 1.12.0 | 58 passed |

The lock input is identical across the two coordinated repositories. The
locks include the Linux-only keyring dependency chain used by GitHub Actions.
Their SHA-256 values are:

- lock input: `800f4b9ad02401fffd503cc95f1de4200510a8fcda7cea9de269d02b7c4f555c`
- Python 3.10: `25692d843eb13f3bd675bef6566e068e2f0efeb1d48e699a952d26833c5b6436`
- Python 3.11: `d41b8b1a82e9034a03d09519b8e6b4ef23e8d8c9a93ea55ed31ee8e4738b8752`
- Python 3.12: `78ce7b0f88f2637bae4375f37bde58bd89d67c466ef954a54b24e62ca4bdb1c3`

Coverage includes deterministic composite matching, optimizer failure and
non-uniqueness handling, raw and optional calibrated metallicity provenance,
cross-validation checkpoints and sidecar hashes, DRP compatibility, RV-source
screening, library-manifest verification, archive extraction safety, and
target/reference DQ gates.

The DRP guard follows the NEID release policy: patch releases within one
major/minor series are compatible and their exact values remain recorded;
different minor series are rejected by default. Thus v1.5.2 and v1.5.3 may be
validated together, while v1.4.x and v1.5.x may not.

## Live archive availability audit

The public NEID TAP metadata service was queried by exact Level-2 filename on
2026-08-31.

- All 78 reference-library exposures have v1.5.x products: 77 at v1.5.3 and
  one at v1.5.2. The one local filename ending in `_corrected.fits` maps to its
  original archive exposure, which is available at v1.5.3.
- All 78 unique Level-2 filenames currently present under the SURFSUP `data/`
  tree have v1.5.x products: 73 at v1.5.2 and five at v1.5.3.
- All 78 queried SURFSUP products have archive `rejected=0` and `flagged=0`.
- SHA-256 of the sorted 78-row reference metadata CSV returned by the query:
  `285693cac65f16bd9f7ad171782534899bfc913bb8d20b0e0ecdd41c34b69cea`.
- SHA-256 of the 78-row SURFSUP metadata CSV returned by the query:
  `1c2523c7d8143f86fcefdf70f0436d0d97ad0311bdff95861290187ab3ab0c8d`.

The query results establish availability, not possession or scientific
validation. The reference files are public. SURFSUP downloads may require the
maintainer's archive credentials; credentials must be supplied only through
the environment/private cookie workflow.

## One-file current-DRP comparison

The public reference exposure `neidL2_20210116T093238.fits` was downloaded from
the live archive. The returned v1.5.3 file had MD5
`ffe8dd0ff23448c78332b46fe096416d`, matching the archive metadata, and SHA-256
`f5ac87431a4506f325e666b8fe4d0337066386772f70cb198301f4b96dc2f444`.

Compared with the archived v1.3.0 copy:

- extracted flux, variance, and blaze arrays were identical on rows 55 and
  102 for this exposure;
- the wavelength solution changed, reaching an absolute difference of about
  465 m/s on row 55 and 4 m/s on row 102;
- the combined DRP CCF RV changed by about 9 m/s; and
- the v1.5.3 product adds stitched-spectrum HDUs 18--20.

This single comparison suggests that some inputs are stable, but it cannot
replace full empirical cross-validation. Differential wavelength changes are
precisely why an old-v1.3 library must not supply publication uncertainties
for v1.5 targets.

## Artifact rehearsal and validation boundary

A disposable source copy produced one wheel and one sdist. `twine check`,
`check-wheel-contents`, and exact package-data allowlisting passed. The exact
wheel was then installed outside the checkout in all three locked Python
environments: all 58 tests and `pip check` passed, version 0.2.0b1 imported, and
all three command-line entry points displayed help. The rehearsal wheel was
about 86 KiB and the sdist about 241 KiB. Final artifacts must be rebuilt after
this inspection and must contain the release-owned reference-library SHA-256
allowlist but no spectra, legacy catalogs, notebooks, plots, credentials, or
local paths.

The archived v1.3 cross-validation results and paper equations are historical
only. Current free-`v sin i` performance metrics require downloading the v1.5
reference products, producing a new deep library manifest, and rerunning
leave-one-out validation with this exact code. Residual RMSE is not
automatically a Gaussian or per-target `1 sigma` uncertainty. Fixed or bounded
broadening, fitted rotational velocity, and empirical metallicity calibration
remain unvalidated unless separately supported by appropriate independent
evidence.
