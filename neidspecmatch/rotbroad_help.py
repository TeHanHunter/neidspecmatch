"""Compatibility re-export of the canonical :mod:`neidspec` broadening API.

NEIDSpecMatch and NEIDSpec intentionally use one scientific implementation.
The functions live in the required ``neidspec>=0.2`` dependency; this module
only preserves the historical NEIDSpecMatch import path.
"""

from neidspec.rotbroad_help import (
    broaden,
    broaden_variance,
    rot,
    rotational_kernel,
    rotmacro,
    rotmacro_ft,
)

__all__ = [
    "broaden",
    "broaden_variance",
    "rot",
    "rotational_kernel",
    "rotmacro",
    "rotmacro_ft",
]
