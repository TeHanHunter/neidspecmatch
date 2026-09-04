# Release checklist

Use this checklist for every NEIDSpecMatch release. A release is a scientific
artifact: the source commit, distributions, reference library, calibrations,
and reported validation results must agree.

## Scope and dependencies

- [ ] Freeze the intended source commit and release notes; update
  `CHANGELOG.md` from "unreleased" only after the contents are final.
- [ ] Release and verify the compatible `neidspec>=0.2.1,<0.3` dependency on
  PyPI before creating a NEIDSpecMatch release tag.
- [ ] Confirm supported Python versions (3.10, 3.11, and 3.12) and dependency
  bounds against clean environments.
- [ ] Freeze reviewed, hash-pinned build/test constraints for each supported
  Python minor and archive them with the release record. Keep broad runtime
  compatibility metadata separate from the exact release environment.
- [ ] Confirm optional archive and utility dependencies are imported
  only by the features that need them.
- [ ] Verify package version, git tag (`v<version>`), changelog, citation
  metadata, and documentation all agree.
- [ ] Protect the release-tag pattern in GitHub; the publish job also checks
  `github.ref_protected` and will skip an unprotected tag.

## Security and data rights

- [ ] Run secret scanning over the full git history and every distribution.
  Confirm any historically exposed credential remains revoked or rotated.
- [ ] Confirm the revoked NEID Archive password is absent from the 0.2.0
  source tree, wheel, and sdist. The verified affected PyPI releases are 0.1.1,
  0.1.2, 0.1.3, and 0.1.4; do not print or reuse the historical literal.
- [ ] After 0.2.0 is available, yank PyPI releases 0.1.1 through 0.1.4 with
  the reason "revoked hard-coded NEID Archive credential; upgrade to 0.2.0".
  Record the action and verify the PyPI API flags. Do not confuse yanking with
  deletion or history removal.
- [ ] Exercise the malicious-ZIP, checksum-failure, and partial-install tests.
- [ ] Reject symlinks, absolute paths, `..` traversal, unexpected layout, and
  executable serialization from untrusted library archives.
- [ ] Resolve and record the license, source, citation, version, and digest for
  every data file. Keep the replacement empirical-library Zenodo version as an
  unpublished draft until an authorized NEID representative confirms in
  writing that the Level-2 FITS may be redistributed and states the applicable
  terms. Do not infer that permission from public download access or the old
  Zenodo record's uploader-supplied CC BY label. Do not apply deposit-wide
  CC BY terms to the NEID FITS. See `DATA_LICENSES.md` and
  `docs/NEID_DATA_RIGHTS_DECISION.md`.
- [ ] Do not publish the response file, catalogs, or literature-derived assets
  while their redistribution rights remain unresolved.
- [ ] Confirm `git ls-files lib library tests/img` is empty. Removing unresolved data
  from the release branch does not remove it from existing history; treat
  history rewrite as a separate destructive remediation decision.
- [ ] Keep only maintained, path-agnostic notebooks with empty outputs and null
  execution counts. Tutorials must demonstrate bounded one-order-at-a-time
  runs, not eagerly load the full library across every order.
- [ ] Treat the historical Zenodo MD5 as accidental-corruption detection only.
  Recompute every catalog/FITS SHA-256 from the reviewed source and verify the
  packaged release-owned allowlist before building; update its ID whenever the
  archive payload changes.
- [ ] Inspect the wheel and sdist for credentials, local paths, `.DS_Store`,
  caches, raw spectra, target files, and unintended large data.

## Scientific validation

- [ ] Freeze the empirical-library identifier and manifest, NEID DRP versions,
  embedded Level-2 blaze provenance, order definitions, and software commit.
- [ ] Re-run ordinary free-`vsini` cross-validation for every supported
  order/library/population and archive raw predictions, fold assignments,
  metrics, and logs. Do not apply those results to fixed or bounded fits.
- [ ] Keep raw [Fe/H] as the default unless an order-, library-, and
  population-specific calibration improves held-out performance and ships
  with validation provenance. Never reuse HPF order-5 coefficients for NEID.
- [ ] Report bias, sample scatter, RMSE, finite coverage, and their exact
  pipeline/calibration provenance. Do not label RMSE as an automatic Gaussian,
  68%-coverage, or per-target `1 sigma` uncertainty.
- [ ] Verify externally constrained `vsini` is immutable across all reference
  stars and that result metadata records the supplied center and window. Keep
  atmospheric parameters from fixed and bounded modes exploratory until an
  independent labeled broad-lined validation set covers the science domain.
- [ ] Treat target-broadening injections with propagated variance and masks as
  stress tests only; do not let synthetic recovery unlock publication-valid
  status for fixed or bounded fits.
- [ ] Compare a frozen set of representative spectra against independently
  reviewed expected results, including edge cases and failure modes.
- [ ] Record target/library RV sources, barycentric-correction provenance, and
  instrument/site metadata in machine-readable outputs.

## Build and test

- [ ] Pin every third-party GitHub Action to a reviewed full commit SHA and
  record the review/update date; mutable major-version tags are not a final
  supply-chain control.
- [ ] For each pin update, resolve and dereference the official signed release
  tag, review its changelog, then update the SHA, version comment, and review
  date together. Do not auto-merge an action-pin change.
- [ ] Run `python -m unittest discover -s tests -p 'test_*.py' -v` on Python
  3.10, 3.11, and 3.12 in clean environments.
- [ ] Build once with `python -m build`; never rebuild between testing and
  publication.
- [ ] Run `python -m twine check dist/*` and `check-wheel-contents dist/*.whl`.
- [ ] Confirm every artifact remains below the workflow's 5 MiB software-only
  ceiling. Raise that ceiling only after a reviewed code-size justification,
  never to accommodate the external reference library or other scientific data.
- [ ] Install the exact wheel into clean environments, run the tests, run
  `python -m pip check`, and import the package outside the source checkout.
- [ ] Smoke-test all four installed commands outside the checkout:
  `neidspecmatch-fit`, `neidspecmatch-crossval`,
  `neidspecmatch-crossval-figure`, and `neidspecmatch-library`.
- [ ] Compare local, GitHub-release, and PyPI artifacts by filename, size, and
  SHA-256 digest; verify wheel `METADATA`, `RECORD`, licenses, and source commit.
- [ ] Upload the tested distributions to the GitHub Actions artifact store for
  the protected publish job. Do not rebuild in the publish job.

### Exact artifact-parity check

After publication, download rather than rebuild all three copies. The PyPI
step verifies its API-declared SHA-256 before the final byte comparison:

```bash
release_audit_dir=$(mktemp -d)
run_id=REPLACE_WITH_WORKFLOW_RUN_ID
version=REPLACE_WITH_VERSION
mkdir -p "$release_audit_dir/ci" "$release_audit_dir/github" "$release_audit_dir/pypi"
gh run download "$run_id" --name distributions --dir "$release_audit_dir/ci"
gh release download "v$version" --dir "$release_audit_dir/github"
python - "$version" "$release_audit_dir/pypi" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

version, destination = sys.argv[1], Path(sys.argv[2])
payload = json.load(urlopen("https://pypi.org/pypi/neidspecmatch/json", timeout=30))
for item in payload["releases"][version]:
    if item["packagetype"] not in {"bdist_wheel", "sdist"}:
        continue
    data = urlopen(item["url"], timeout=60).read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != item["digests"]["sha256"]:
        raise SystemExit(f"PyPI digest mismatch: {item['filename']}")
    (destination / item["filename"]).write_bytes(data)
PY
python - "$release_audit_dir/ci" "$release_audit_dir/github" "$release_audit_dir/pypi" <<'PY'
import hashlib
from pathlib import Path
import sys

directories = [Path(value) for value in sys.argv[1:]]
names = [{item.name for item in path.iterdir() if item.is_file()} for path in directories]
if not all(value == names[0] for value in names[1:]):
    raise SystemExit(f"Artifact filename sets differ: {names}")
for name in sorted(names[0]):
    records = []
    for directory in directories:
        data = (directory / name).read_bytes()
        records.append((len(data), hashlib.sha256(data).hexdigest()))
    print(name, records[0])
    if not all(value == records[0] for value in records[1:]):
        raise SystemExit(f"Artifact bytes differ for {name}: {records}")
PY
```

## Trusted publication and post-release checks

- [ ] Configure a PyPI Trusted Publisher for this repository and
  `.github/workflows/release.yml`, scoped to the protected `pypi` environment.
- [ ] Require maintainer approval on the GitHub `pypi` environment and restrict
  who can deploy. Keep `id-token: write` only on the publish job.
- [ ] Create the signed `v<version>` tag only after all gates pass. Do not use a
  branch build as a release candidate for publication.
- [ ] Approve the protected publish job, then verify the PyPI page, provenance,
  hashes, and a clean installation from PyPI on every supported Python version.
- [ ] Create the matching GitHub/Zenodo release and confirm citation metadata.
- [ ] Retain the validation report, artifact digests, workflow URL, commit, tag,
  and approver as the permanent release record.
