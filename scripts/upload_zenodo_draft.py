#!/usr/bin/env python3
"""Upload one file to an existing unpublished Zenodo deposition draft.

The access token is read from ``ZENODO_TOKEN`` or from a hidden prompt.  This
tool uploads files only; it deliberately has no operation that publishes a
record.
"""

from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen


CHUNK_BYTES = 8 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("deposition_id", help="numeric Zenodo draft identifier")
    parser.add_argument("file", type=Path)
    parser.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org")
    parser.add_argument(
        "--inspect-only",
        action="store_true",
        help="authenticate and inspect the draft without uploading",
    )
    return parser.parse_args()


def token() -> str:
    value = os.environ.get("ZENODO_TOKEN")
    if value is None:
        value = getpass.getpass("Zenodo token (hidden): ")
    value = value.strip()
    if not value:
        raise SystemExit("A Zenodo token is required.")
    return value


def get_deposition(base_url: str, deposition_id: str, access_token: str) -> dict:
    url = f"{base_url}/api/deposit/depositions/{quote(deposition_id)}"
    request = Request(
        url,
        headers={"Authorization": f"Bearer {access_token}"},
    )
    try:
        with urlopen(request, timeout=60) as response:
            return json.load(response)
    except (HTTPError, URLError) as exc:
        raise SystemExit(f"Could not inspect Zenodo draft {deposition_id}: {exc}") from exc


def local_md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def upload(bucket_url: str, path: Path, access_token: str) -> dict:
    target = urlparse(f"{bucket_url.rstrip('/')}/{quote(path.name)}")
    if target.scheme != "https" or not target.hostname:
        raise SystemExit(f"Zenodo returned an invalid bucket URL: {bucket_url!r}")
    curl = shutil.which("curl")
    if curl is None:
        raise SystemExit("curl is required for robust large-file upload.")

    # Keep the bearer token out of the command line and process listings by
    # sending the sensitive header through curl's stdin config. Curl retries
    # transient transport failures automatically; Zenodo does not provide a
    # documented byte-range resume mechanism, so each retry starts the PUT
    # from byte zero.
    curl_config = f'header = "Authorization: Bearer {access_token}"\n'
    command = [
        curl,
        "--config",
        "-",
        "--http1.1",
        "--fail-with-body",
        "--show-error",
        "--progress-bar",
        "--retry",
        "8",
        "--retry-all-errors",
        "--retry-delay",
        "10",
        "--connect-timeout",
        "60",
        "--header",
        "Content-Type: application/octet-stream",
        "--upload-file",
        str(path),
        "--url",
        target.geturl(),
    ]
    completed = subprocess.run(
        command,
        input=curl_config,
        text=True,
        stdout=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(f"Zenodo upload failed after retries (curl exit {completed.returncode}).")
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"Zenodo returned a non-JSON response: {completed.stdout[:1000]!r}"
        ) from exc


def main() -> None:
    args = parse_args()
    path = args.file.expanduser().resolve()
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"Not a regular file: {path}")
    base_url = "https://sandbox.zenodo.org" if args.sandbox else "https://zenodo.org"
    access_token = token()
    deposition = get_deposition(base_url, args.deposition_id, access_token)
    bucket = deposition.get("links", {}).get("bucket")
    if not bucket:
        raise SystemExit("Zenodo draft response did not contain a writable bucket URL.")
    print(
        f"Draft {args.deposition_id}: {deposition.get('metadata', {}).get('title', '(untitled)')}"
    )
    if args.inspect_only:
        print("Draft is accessible; no file was uploaded.")
        return

    expected_md5 = local_md5(path)
    print(f"Uploading {path.name} ({path.stat().st_size:,} bytes; MD5 {expected_md5})")
    uploaded = upload(bucket, path, access_token)
    remote_checksum = str(uploaded.get("checksum", "")).removeprefix("md5:")
    if remote_checksum and remote_checksum.lower() != expected_md5.lower():
        raise SystemExit(
            f"Upload completed but checksum differs: {remote_checksum} != {expected_md5}"
        )
    print(f"Upload complete; Zenodo checksum: {uploaded.get('checksum', '(not returned)')}")


if __name__ == "__main__":
    main()
