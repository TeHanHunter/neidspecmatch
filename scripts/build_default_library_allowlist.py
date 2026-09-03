#!/usr/bin/env python3
"""Build the packaged receipt for the versioned default reference library."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


PROVENANCE_FIELDS = (
    "archive_products",
    "reference_dq_records",
    "source_archive",
    "source_catalog",
)


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--doi", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--allowlist-id", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    missing = [field for field in PROVENANCE_FIELDS if field not in manifest]
    if missing:
        raise SystemExit(f"Manifest is missing provenance fields: {missing}")

    allowlist = {
        "allowlist_id": args.allowlist_id,
        "catalog": manifest["catalog"],
        "files": manifest["files"],
        "fits_count": manifest["fits_count"],
        "header_identity_exceptions": manifest.get(
            "header_identity_exceptions", []
        ),
        "library_id": manifest["library_id"],
        "library_provenance": {
            field: manifest[field] for field in PROVENANCE_FIELDS
        },
        "schema_version": manifest["schema_version"],
        "source_archive_integrity": {
            "algorithm": "md5",
            "scope": "integrity_only_not_authentication",
            "value": file_digest(args.archive, "md5"),
        },
        "source_archive_size_bytes": args.archive.stat().st_size,
        "source_doi": args.doi,
        "source_rights": (
            "NEID Level-2 redistribution guidance pending; "
            "see DATA_LICENSES.md"
        ),
        "source_url": args.url,
    }
    args.output.write_text(
        json.dumps(allowlist, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
