"""Filesystem-safe naming helpers.

Scientific object names are metadata, not paths.  This module deliberately
keeps the two concepts separate.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import re
import unicodedata


_UNSAFE_COMPONENT = re.compile(r"[^A-Za-z0-9._-]+")


def safe_filename_component(value, *, max_length=80):
    """Return a deterministic, traversal-safe component for arbitrary text.

    A readable ASCII slug is followed by a short hash of the exact original
    string.  The hash prevents collisions such as ``"A/B"`` and ``"A B"``
    while the slug keeps output directories recognizable to a person.
    """
    original = str(value)
    digest = hashlib.sha256(original.encode("utf-8")).hexdigest()[:10]
    ascii_value = unicodedata.normalize("NFKD", original).encode(
        "ascii", "ignore"
    ).decode("ascii")
    slug = _UNSAFE_COMPONENT.sub("_", ascii_value).strip("._-")
    if not slug:
        slug = "target"
    suffix = "-" + digest
    maximum = int(max_length)
    if maximum < len(suffix) + 1:
        raise ValueError("max_length is too short for a safe filename component")
    slug = slug[:maximum - len(suffix)].rstrip("._-") or "target"
    return slug + suffix


def sha256_file(filename, *, chunk_size=1024 * 1024):
    """Return the SHA-256 digest of a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(filename).open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
