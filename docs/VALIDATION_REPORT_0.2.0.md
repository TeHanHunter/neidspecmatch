# NEIDSpecMatch 0.2.0 validation report

Status: **release candidate; four-order validation complete**
Validation date: 2026-09-03

This report records software, artifact, archive, and cross-validation
evidence. The release validation applies only to the ordinary free-`v sin i`
estimator for the exact library/order/population combinations reported below.
Fixed and bounded pairwise-`v sin i` fits remain exploratory pending
independent labeled broad-lined validation; synthetic broadening alone is a
stress test, not publication validation.

## Automated regression matrix

The same 86-test suite was run with warnings treated as errors. Each
environment was installed from its Python-specific, fully hashed release lock
with `pip --require-hashes`; `pip check` reported no broken requirements.

| Python | Representative numerical environment | Result |
| --- | --- | --- |
| 3.10.14 | NumPy 1.26.4; SciPy 1.12.0 | 86 passed |
| 3.11.11 | NumPy 1.26.4; SciPy 1.12.0 | 86 passed |
| 3.12.2 | NumPy 1.26.4; SciPy 1.12.0 | 86 passed |

The lock input is identical across the two coordinated repositories. The
locks include the Linux-only keyring dependency chain used by GitHub Actions.
Their SHA-256 values are:

- lock input: `800f4b9ad02401fffd503cc95f1de4200510a8fcda7cea9de269d02b7c4f555c`
- publication-validation runtime:
  `3336c8922cc74c329e3a11c04dcf5179cb9e84b7a64405fbbc26382becf4dc24`
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
validated together, while v1.4.x and v1.5.x may not. The coordinated
`neidspec 0.2.1` dependency was confirmed available from PyPI on 2026-09-03.

## Frozen DRP-1.5 empirical library

The default release library identifier is
`20260831_specmatch_neid_drp15`. It contains one catalog and 78 NEID Level-2
HR spectra: 77 reduced with DRP v1.5.3 and one with v1.5.2. Seventy-two
reference products pass the package DQ screen and six carry explicitly
accepted warning status; that exact composition is bound into the validation
receipt.

- Reserved dataset DOI: `10.5281/zenodo.22262405`
- ZIP size: `6,536,877,943` bytes
- ZIP MD5, verified locally and returned by Zenodo after upload:
  `b8abc3ade339074cbd137a9fd2749297`
- ZIP SHA-256:
  `4edbd37203236ac4969d3f679d79a131d012417ad3d59a8628670625717979e8`
- Packaged release allowlist SHA-256:
  `2f7a5ef3de977c38e16670769476d29dcee5317cafbda282c5daa7b2f0bd9396`
- Installed deep library-manifest SHA-256:
  `55a1be296e756859598377b62b1dc4a654ba61aabe53368c49a22431a5f29006`

The archive is present in an unpublished Zenodo draft. The draft must remain
unpublished until an authorized NEID representative confirms redistribution
terms in writing. Upload and checksum verification establish integrity, not
permission to redistribute.

## Four-order cross-validation

Leave-one-spectrum-out validation used the ordinary free-`v sin i` estimator
for array-row orders 55, 101, 102, and 103. All 78 folds completed at every
order (312 fits total). Every order reports the full library, 46 cool stars,
and 32 hot stars, divided at 4500 K. Raw `[Fe/H]` is reported; no empirical
detrending is applied.

| Order | Population | N | Teff RMSE (K) | [Fe/H] RMSE (dex) | log(g) RMSE (dex) |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 55 | cool | 46 | 67.25 | 0.1628 | 0.0474 |
| 55 | hot | 32 | 103.98 | 0.0866 | 0.0887 |
| 101 | cool | 46 | 72.09 | 0.1756 | 0.0578 |
| 101 | hot | 32 | 150.14 | 0.0711 | 0.0810 |
| 102 | cool | 46 | 60.92 | 0.1628 | 0.0387 |
| 102 | hot | 32 | 152.11 | 0.0546 | 0.0782 |
| 103 | cool | 46 | 56.59 | 0.1627 | 0.0381 |
| 103 | hot | 32 | 137.42 | 0.0527 | 0.0935 |

The deposited summary CSV also reports all-star metrics, residual bias,
sample scatter, finite coverage, and unrounded values. The completed run was
made with the immediately preceding 0.2.0b1 candidate and records, without
relabeling, these fingerprints:

- NEIDSpecMatch 0.2.0b1 source:
  `fd21efbd772ac88efb5e06b5fb289f49ec26b3b3126087a5ab874b0529cb1186`
- NEIDSpec 0.2.0 source:
  `c494c1af665aa71b29c89ef03ab8a1c7f42f0439c89077c1f0e5526b9bcb13bf`
- library manifest:
  `1b52fe513ffdd7c198bdc95800fea6956a7ad99148d30cabec217b7d5c8e03f6`

After the release-only provenance, packaging, and coordinated NEIDSpec mask-
medium updates, 129 overlapping cool-star folds were rerun with the final
candidate. Median absolute changes were 0.77--1.63 K in Teff,
0.0016--0.0039 dex in `[Fe/H]`, and 0.0006--0.0010 dex in `log(g)`; maxima
were 12.7 K, 0.0248 dex, and 0.0160 dex. These are numerical-continuity checks,
not replacement cross-validation. Their row-level data and provenance are in
the release bundle. The final candidate fingerprints are:

- NEIDSpecMatch 0.2.0 source:
  `4c2dcc6bca85d7ac55f987db5879f7a6a15080202cbfb5d40a5d363e14b88215`
- NEIDSpec 0.2.1 source:
  `f49070f7f18499e85c0c9cb37596a6a9ece3c18eef56bb80ea6e617a599f6ba4`

The original raw predictions and checkpoints are preserved with their true
0.2.0b1 provenance. They are not rewritten to masquerade as an exact final-
source run; the final package's strict runtime-attached validation gate
therefore remains conservative.

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

A frozen release-candidate source copy produced one wheel and one sdist.
`twine check`,
`check-wheel-contents`, and exact package-data allowlisting passed. The exact
wheel was then installed outside the checkout in all three locked Python
environments: all release tests and `pip check` passed, version 0.2.0 imported, and
all four command-line entry points displayed help. The current wheel is about
120 KiB and the sdist about 315 KiB. Final artifacts will be rebuilt once after
the validation bundle is frozen and must contain the release-owned reference-library SHA-256
allowlist but no spectra, legacy catalogs, notebooks, plots, credentials, or
local paths.

The archived v1.3 cross-validation results and paper equations are historical
only. Residual RMSE is not automatically a Gaussian or per-target `1 sigma`
uncertainty. Fixed or bounded broadening, fitted rotational velocity, and
empirical metallicity calibration remain unvalidated unless separately
supported by appropriate independent evidence.
