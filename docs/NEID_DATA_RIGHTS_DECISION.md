# NEID reference-library rights decision

**Decision: hold publication of the replacement Zenodo library.**

Review date: 2026-09-02

The 20260831 candidate library contains 78 unmodified NEID Level-2 FITS files
that are publicly downloadable from the NEID Archive and were reduced with
DRP 1.5.2 or 1.5.3. Public availability establishes access, but the release
review did not find an explicit grant allowing a third party to republish or
sublicense those FITS files. The existing Zenodo record's CC BY 4.0 label is
uploader-supplied metadata and is not sufficient evidence of such authority.

## Evidence reviewed

The official [NEID Data FAQ](https://neid.ipac.caltech.edu/docs/NEID-DRP/faq.html)
states that data with no proprietary period, or whose proprietary period has
expired, can be downloaded without an account. The official
[NEID Archive help page](https://neid.ipac.caltech.edu/help.php) gives the
required publication acknowledgement and identifies
`neid-help@ipac.caltech.edu` as the contact. Neither page states a data license
or grants third parties permission to redistribute or sublicense Level-2 FITS
files. Access, acknowledgement, redistribution, and licensing are distinct
questions.

The new-version Zenodo draft may be populated for staging, but it must remain
unpublished. Do not assign a deposit-wide CC BY license that purports to cover
the NEID FITS. The following materials may be identified as CC BY 4.0 only to
the extent they are wholly depositor-authored:

- original README and data-rights prose;
- original manifests and checksum lists; and
- original validation summaries, figures, and provenance records.

That grant must explicitly exclude the NEID FITS, archive-derived metadata,
and third-party catalog labels. The NEIDSpecMatch software source is separately
licensed under MIT.

## Publication gate

Send `NEID_DATA_PERMISSION_REQUEST.txt` to `neid-help@ipac.caltech.edu` and
retain the written response in the private release record. Publish the FITS
archive only if an authorized NEID representative confirms redistribution and
provides the applicable rights statement, acknowledgement, and citations.

If permission is denied or remains unclear, do not publish the FITS archive.
Instead, distribute a verified manifest and a reproducible download procedure
that lets each user retrieve the public products directly from the NEID
Archive, subject to its access and acknowledgement terms.

Changing a DOI, checksum, catalog, FITS filename, or software license does not
resolve this gate.
