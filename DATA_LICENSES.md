# Data and model licenses

The MIT license in `LICENSE` applies to the NEIDSpecMatch software source. It
does **not** grant rights to redistribute spectra, catalogs, calibration data,
instrument-response files, masks, figures derived from third-party material,
or other external scientific data.

## Empirical reference library

The approximately 7.3 GB NEIDSpecMatch empirical-library archive is hosted at
Zenodo record [10.5281/zenodo.14947454](https://doi.org/10.5281/zenodo.14947454).
Zenodo identifies Te Han as the creator, reports an archive size of
7,337,841,681 bytes, and licenses the deposit under
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
Reuse and redistribution must preserve the required attribution and license
notice. The archive remains an external, separately downloaded dataset and is
not included in the Python wheel or source distribution.

The historical MD5 value
`e54e203610e948512e641b6e30530570` can detect accidental corruption of that
specific archive. It is integrity metadata only: it is not evidence of
authenticity, scientific provenance, or fitness for analysis. A future archive
release should publish an archive-level SHA-256 digest alongside the license
metadata. Release 0.2.0b1 packages a versioned, release-owned SHA-256 allowlist
for the catalog and all 78 extracted reference FITS files; the downloader uses
that allowlist to authenticate the extracted dataset before it is accepted.

## Other scientific assets

Licensing and redistribution permission have not yet been resolved for the
following bundled or historically distributed inputs:

- catalogs and tables under `lib/`, including literature-derived metadata;
- the NEID master instrument-response FITS file;
- literature tables or spectra under `library/`, including VizieR/ApJ-source
  material and `asu.fit`;
- `combined.csv`, masks, tutorial inputs, and any external spectra copied into
  local working directories; and
- plots or figures that reproduce or are derived from third-party material.

Do not assume that repository visibility or a checksum grants reuse rights.
Before the next public artifact is uploaded, identify the source, citation,
license, redistribution terms, version, and cryptographic digest for every
non-code file. Exclude any file whose status cannot be documented.

Calibration products created from the empirical library require their own
provenance record. That record should identify the input-library version,
selection rules, order, cross-validation procedure, software commit, and
applicable data license. Generating a derivative does not erase restrictions
on its inputs.
