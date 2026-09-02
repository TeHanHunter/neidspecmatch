# Data and model licenses

The MIT license in `LICENSE` applies to the NEIDSpecMatch software source. It
does **not** grant rights to redistribute spectra, catalogs, calibration data,
instrument-response files, masks, figures derived from third-party material,
or other external scientific data.

## Empirical reference libraries

The historical NEIDSpecMatch empirical-library archive is hosted at Zenodo
record [10.5281/zenodo.14947454](https://doi.org/10.5281/zenodo.14947454).
That record is labeled CC BY 4.0, but record-level license metadata supplied by
an uploader is not, by itself, evidence that the uploader was authorized to
relicense third-party NEID Level-2 FITS files. NEIDSpecMatch therefore does not
rely on that label as a grant of redistribution or sublicensing rights.

The proposed 20260831 library contains 78 unmodified public NEID Level-2 FITS
files reduced with DRP 1.5.2 or 1.5.3. Public download availability and an
acknowledgement requirement do not necessarily grant permission to republish
or apply a Creative Commons license to those files. The replacement Zenodo
version must remain an **unpublished draft** until the NEID Archive or another
authorized NEID representative confirms the permitted redistribution terms in
writing. See `docs/NEID_DATA_RIGHTS_DECISION.md` and the ready-to-send request
in `docs/NEID_DATA_PERMISSION_REQUEST.txt`.

Do not apply a deposit-wide CC BY 4.0 statement to the NEID FITS. CC BY 4.0 may
cover wholly depositor-authored documentation, manifests, and checksum lists,
but only to the extent that those files contain original material and with an
explicit exclusion for the NEID FITS, archive-derived metadata, and
third-party stellar labels. The software remains under the MIT license. Any
eventual data deposit remains external and is not included in the Python wheel
or source distribution.

The historical MD5 value
`e54e203610e948512e641b6e30530570` can detect accidental corruption of that
specific archive. It is integrity metadata only: it is not evidence of
authenticity, scientific provenance, or fitness for analysis. A future archive
release should publish an archive-level SHA-256 digest alongside the license
metadata. Release 0.2.0 packages a versioned, release-owned SHA-256 allowlist
for the catalog and all 78 extracted reference FITS files; the downloader uses
that allowlist to authenticate the extracted dataset before it is accepted.

## Excluded historical scientific assets

Licensing and redistribution permission was not sufficiently documented for
the following historically distributed inputs, so they are intentionally
excluded from the 0.2.0 wheel and source distribution:

- catalogs and tables under `lib/`, including literature-derived metadata;
- the NEID master instrument-response FITS file;
- literature tables or spectra under `library/`, including VizieR/ApJ-source
  material and `asu.fit`;
- `combined.csv`, masks, tutorial inputs, and any external spectra copied into
  local working directories; and
- plots or figures that reproduce or are derived from third-party material.

Do not assume that repository visibility or a checksum grants reuse rights.
Before any future release reintroduces one of these assets, identify its
source, citation, license, redistribution terms, version, and cryptographic
digest. Continue to exclude any file whose status cannot be documented.

Calibration products created from the empirical library require their own
provenance record. That record should identify the input-library version,
selection rules, order, cross-validation procedure, software commit, and
applicable data license. Generating a derivative does not erase restrictions
on its inputs.
