# Security policy

## Supported versions

The 0.2.x line is the only line intended to receive security fixes. The
published 0.1.x releases are unsupported; users should not treat them as a
secure basis for new analysis.

## Historical credential exposure

An audit of every NEIDSpecMatch source distribution published on PyPI found a
nonempty hard-coded NEID Archive password in versions 0.1.1, 0.1.2, 0.1.3,
and 0.1.4. The exposed password has been revoked. Do not use or copy the
credential from those artifacts or from repository history. After 0.2.0 is
available, the maintainers intend to yank all four affected PyPI releases
with the reason "revoked hard-coded NEID Archive credential; upgrade to
0.2.0". Yanking discourages new installation while preserving immutable
release records; it does not erase existing downloads or git history.

## Reporting a vulnerability

Please report suspected vulnerabilities privately through GitHub's private
vulnerability-reporting feature for this repository. If that feature is not
available, contact the maintainers at `tehanhunter@gmail.com`. Do not open a
public issue for an unpatched vulnerability or include secrets, credentials,
private spectra, or proprietary data in a report.

Include the affected version or commit, operating system and Python version,
the smallest safe reproduction, and the expected impact. The maintainers will
acknowledge the report and coordinate disclosure after a fix is available.

## Security-sensitive surfaces

Treat downloaded library archives and FITS files as untrusted input. Archive
installation must reject path traversal and must verify the documented
checksum and expected layout before replacing an existing library. Python
pickle files are executable serialization and must never be loaded from an
untrusted source.

Archive credentials must be supplied through the documented environment
variables or another local secret store. Never place passwords, API tokens,
cookies, or private download URLs in source, examples, test fixtures, release
artifacts, issue reports, or git history.

Package releases use PyPI Trusted Publishing. No long-lived PyPI token should
be stored in GitHub Actions; the `pypi` environment should require maintainer
approval and should be limited to the release workflow.
