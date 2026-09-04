"""Public API for NEIDSpecMatch.

Heavy spectral and plotting modules are loaded only when their public
functions are requested. Import is quiet and performs no network or writes.
"""

from importlib import import_module

from .calibration import (
    CalibrationError,
    FeHCalibration,
    FeHCalibrationSet,
    PopulationRule,
    UnvalidatedCalibrationError,
    build_calibration_products,
    fit_crossvalidated_feh_calibration,
)
from .config import (
    DEFAULT_LIBRARY_CATALOG,
    DEFAULT_LIBRARY_ID,
    LIBRARY_ENVVAR,
    resolve_library_path,
)
from .paths import safe_filename_component
from .version import __version__


_LAZY_ATTRS = {
    "run_specmatch": (".neidspecmatch", "run_specmatch"),
    "run_specmatch_for_orders": (".neidspecmatch", "run_specmatch_for_orders"),
    "run_crossvalidation_for_orders": (
        ".neidspecmatch", "run_crossvalidation_for_orders"
    ),
    "write_crossvalidation_products": (
        ".neidspecmatch", "write_crossvalidation_products"
    ),
    "load_reference_library": (".neidspecmatch", "load_reference_library"),
    "detrend_feh": (".neidspecmatch", "detrend_feh"),
    "get_library": (".utils", "get_library"),
    "validate_library": (".utils", "validate_library"),
    "build_library_manifest": (".utils", "build_library_manifest"),
}


def __getattr__(name):
    try:
        module_name, attribute = _LAZY_ATTRS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name, __name__), attribute)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_ATTRS))


__all__ = [
    "__version__",
    "CalibrationError",
    "UnvalidatedCalibrationError",
    "PopulationRule",
    "FeHCalibration",
    "FeHCalibrationSet",
    "fit_crossvalidated_feh_calibration",
    "build_calibration_products",
    "DEFAULT_LIBRARY_ID",
    "DEFAULT_LIBRARY_CATALOG",
    "LIBRARY_ENVVAR",
    "resolve_library_path",
    "safe_filename_component",
    *_LAZY_ATTRS,
]
