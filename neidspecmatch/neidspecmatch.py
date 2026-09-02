import neidspec
import numpy as np
import pandas as pd
import astropy.io
import astropy
import pickle
import json
import os
import contextlib
import io
import tempfile
import re
import hashlib
import itertools
from neidspec import stats_help
import glob
from collections import Counter
from pathlib import Path, PurePath

import scipy.optimize
import scipy.linalg
import scipy
import sys
from neidspecmatch import rotbroad_help
from neidspec import utils
import matplotlib.pyplot as plt
import astropy.modeling
from neidspecmatch.priors import PriorSet, UP
from neidspecmatch.likelihood import ll_normal_ev_py
import neidspecmatch.config as config
from matplotlib.gridspec import GridSpec
from neidspecmatch.calibration import (
    CalibrationError,
    CROSSVALIDATION_METRIC_SCHEMA_VERSION,
    CROSSVALIDATION_RESIDUAL_METHOD,
    NO_AUTOMATIC_UNCERTAINTY_SEMANTICS,
    PopulationRule,
    UnvalidatedCalibrationError,
    build_calibration_products,
    coerce_calibration_set,
    raw_crossvalidation_summary,
)
from neidspecmatch.paths import safe_filename_component, sha256_file
from neidspecmatch.version import __version__
import neidspecmatch.utils as package_utils


RESULT_SCHEMA_VERSION = 3
ESTIMATOR_CONFIG_SCHEMA_VERSION = 1
TOP_K = 5
CONTINUUM_DEGREE = 5
LIMB_DARKENING_U1 = 0.3
SIMPLEX_ACTIVE_WEIGHT_TOL = 1e-10
MIN_PAIRWISE_PIXELS = 50
MIN_PAIRWISE_PIXEL_FRACTION = 0.25


def _canonical_estimator_config_json(estimator_config):
    """Return the machine-independent representation used for validation."""
    if not isinstance(estimator_config, dict):
        raise TypeError("estimator_config must be a mapping.")
    return json.dumps(
        estimator_config, sort_keys=True, separators=(',', ':'),
        allow_nan=False,
    )


def _hash_estimator_config(estimator_config):
    return hashlib.sha256(
        _canonical_estimator_config_json(estimator_config).encode('utf-8')
    ).hexdigest()


def _reference_rv_state(references):
    """Return a canonical, path-free identity for the library rest frames."""
    records = []
    for reference in references:
        object_id = str(getattr(reference, 'object', ''))
        rv = pd.to_numeric(
            pd.Series([getattr(reference, 'rv', np.nan)]), errors='coerce'
        ).iloc[0]
        if not object_id or not np.isfinite(rv):
            raise ValueError(
                "Every reference requires a nonempty identity and finite RV."
            )
        records.append({
            'object_id': object_id,
            'rv_kms': float(rv),
            'rv_source': _shareable_rv_source(
                getattr(reference, 'rv_source', 'spectrum object')
            ),
        })
    records.sort(key=lambda item: item['object_id'])
    if len({item['object_id'] for item in records}) != len(records):
        raise ValueError("Reference RV state contains duplicate object identities.")
    canonical = json.dumps(
        records, sort_keys=True, separators=(',', ':'), allow_nan=False
    )
    return {
        'reference_count': len(records),
        'records': records,
        'object_rv_source_sha256': hashlib.sha256(
            canonical.encode('utf-8')
        ).hexdigest(),
    }


def _reference_blaze_state(references):
    """Return the realized, identity-bound blaze state of the library."""
    records = []
    for reference in references:
        object_id = str(getattr(reference, 'object', ''))
        source = getattr(reference, 'blaze_source_used', None)
        if not object_id or source is None or not str(source).strip():
            raise ValueError(
                "Every reference requires a nonempty identity and a realized "
                "blaze source."
            )
        records.append({
            'object_id': object_id,
            'blaze_source': str(source),
        })
    records.sort(key=lambda item: item['object_id'])
    if len({item['object_id'] for item in records}) != len(records):
        raise ValueError(
            "Reference blaze state contains duplicate object identities."
        )
    canonical = json.dumps(
        records, sort_keys=True, separators=(',', ':'), allow_nan=False
    )
    return {
        'reference_count': len(records),
        'sources': sorted({item['blaze_source'] for item in records}),
        'records': records,
        'object_blaze_source_sha256': hashlib.sha256(
            canonical.encode('utf-8')
        ).hexdigest(),
    }


def _support_diagnostics_status(support):
    """Validate the self-consistency of one empirical-support receipt."""
    def invalid(reason):
        return {'valid': False, 'reason': reason}

    if not isinstance(support, dict):
        return invalid('not_a_mapping')
    required = {
        'simplex', 'parameters', 'selected_components',
        'any_full_library_edge', 'vertex_at_full_library_edge',
        'publication_support_validated', 'warnings',
    }
    missing = sorted(required.difference(support))
    if missing:
        return invalid('missing_fields:' + ','.join(missing))
    simplex = support['simplex']
    components = support['selected_components']
    parameters = support['parameters']
    if not isinstance(simplex, dict) or not isinstance(components, list):
        return invalid('invalid_simplex_or_components')
    if not isinstance(parameters, dict):
        return invalid('invalid_parameters')
    if len(components) != TOP_K:
        return invalid('selected_component_count_mismatch')
    try:
        selected_count = int(simplex['selected_count'])
        active_count = int(simplex['active_count'])
        tolerance = float(simplex['active_weight_tolerance'])
        effective_count = float(simplex['effective_component_count'])
        maximum_weight = float(simplex['maximum_weight'])
    except (KeyError, TypeError, ValueError, OverflowError):
        return invalid('invalid_simplex_numeric_fields')
    if (
        selected_count != TOP_K or not 1 <= active_count <= TOP_K
        or not np.isfinite(tolerance) or tolerance <= 0.0
        or not np.isclose(
            tolerance, SIMPLEX_ACTIVE_WEIGHT_TOL, rtol=0.0, atol=0.0
        )
    ):
        return invalid('invalid_simplex_dimensions_or_tolerance')
    bool_fields = (
        'is_vertex', 'is_boundary', 'optimizer_nonunique'
    )
    if any(not isinstance(simplex.get(name), (bool, np.bool_))
           for name in bool_fields):
        return invalid('invalid_simplex_boolean_fields')
    component_ids = [str(item.get('object_id', '')) for item in components]
    if any(not value for value in component_ids) or len(set(component_ids)) != TOP_K:
        return invalid('invalid_or_duplicate_component_identity')
    try:
        ranks = [int(item['pairwise_rank']) for item in components]
        weights = np.asarray(
            [float(item['weight']) for item in components], dtype=float
        )
    except (KeyError, TypeError, ValueError, OverflowError):
        return invalid('invalid_component_rank_or_weight')
    if ranks != list(range(1, TOP_K + 1)):
        return invalid('component_rank_sequence_mismatch')
    if (
        not np.all(np.isfinite(weights)) or np.min(weights) < -1e-12
        or not np.isclose(weights.sum(), 1.0, rtol=0.0, atol=1e-10)
    ):
        return invalid('component_weights_not_on_simplex')
    active = weights > tolerance
    measured_active = int(np.count_nonzero(active))
    if (
        measured_active != active_count
        or bool(simplex['is_vertex']) != (active_count == 1)
        or bool(simplex['is_boundary']) != (active_count < TOP_K)
        or not np.isclose(
            effective_count, 1.0 / np.sum(weights ** 2),
            rtol=1e-12, atol=1e-12,
        )
        or not np.isclose(
            maximum_weight, np.max(weights), rtol=0.0, atol=1e-12
        )
    ):
        return invalid('simplex_summary_inconsistent_with_weights')
    for index, component in enumerate(components):
        if not isinstance(component.get('active'), (bool, np.bool_)):
            return invalid('component_active_flag_invalid')
        if bool(component['active']) != bool(active[index]):
            return invalid('component_active_flag_mismatch')

    component_fields = {
        'teff': 'teff_k', 'feh_raw': 'feh_raw_dex', 'logg': 'logg_dex'
    }
    measured_edges = []
    for name, component_field in component_fields.items():
        details = parameters.get(name)
        if not isinstance(details, dict):
            return invalid(f'{name}_support_missing')
        numeric_fields = (
            'raw_value', 'full_library_min', 'full_library_max',
            'selected_min', 'selected_max',
        )
        try:
            numeric = {
                field: float(details[field]) for field in numeric_fields
            }
            labels = np.asarray([
                float(item[component_field]) for item in components
            ], dtype=float)
        except (KeyError, TypeError, ValueError, OverflowError):
            return invalid(f'{name}_support_numeric_fields_invalid')
        if not all(np.isfinite(value) for value in numeric.values()):
            return invalid(f'{name}_support_nonfinite')
        if not np.all(np.isfinite(labels)):
            return invalid(f'{name}_component_labels_nonfinite')
        if (
            numeric['full_library_min'] > numeric['selected_min']
            or numeric['selected_min'] > numeric['selected_max']
            or numeric['selected_max'] > numeric['full_library_max']
        ):
            return invalid(f'{name}_support_bounds_inconsistent')
        scale = max(
            1.0, abs(numeric['full_library_min']),
            abs(numeric['full_library_max'])
        )
        atol = 32.0 * np.finfo(float).eps * scale
        expected = {
            'at_full_library_lower_edge': bool(np.all(
                np.abs(labels[active] - numeric['full_library_min']) <= atol
            )),
            'at_full_library_upper_edge': bool(np.all(
                np.abs(labels[active] - numeric['full_library_max']) <= atol
            )),
            'at_selected_lower_edge': bool(np.all(
                np.abs(labels[active] - np.min(labels)) <= atol
            )),
            'at_selected_upper_edge': bool(np.all(
                np.abs(labels[active] - np.max(labels)) <= atol
            )),
        }
        if any(
            not isinstance(details.get(field), (bool, np.bool_))
            or bool(details[field]) != value
            for field, value in expected.items()
        ):
            return invalid(f'{name}_edge_flags_inconsistent')
        expected_limited = (
            expected['at_full_library_lower_edge']
            or expected['at_full_library_upper_edge']
        )
        if (
            not isinstance(details.get('support_limited'), (bool, np.bool_))
            or bool(details['support_limited']) != expected_limited
            or not np.isclose(
                numeric['raw_value'], np.dot(weights, labels),
                rtol=1e-12, atol=atol,
            )
            or not np.isclose(
                numeric['selected_min'], np.min(labels), rtol=0.0, atol=atol
            )
            or not np.isclose(
                numeric['selected_max'], np.max(labels), rtol=0.0, atol=atol
            )
        ):
            return invalid(f'{name}_support_summary_inconsistent')
        measured_edges.append(expected_limited)
    any_edge = bool(any(measured_edges))
    if (
        not isinstance(support['any_full_library_edge'], (bool, np.bool_))
        or bool(support['any_full_library_edge']) != any_edge
        or not isinstance(
            support['vertex_at_full_library_edge'], (bool, np.bool_)
        )
        or bool(support['vertex_at_full_library_edge'])
        != (bool(simplex['is_vertex']) and any_edge)
        or not isinstance(
            support['publication_support_validated'], (bool, np.bool_)
        )
        or bool(support['publication_support_validated']) != (not any_edge)
        or not isinstance(support['warnings'], list)
    ):
        return invalid('overall_support_flags_inconsistent')
    return {'valid': True, 'reason': 'self_consistent'}


def _vsini_configuration(*, maxvsini, vsini=None, vsini_window=None):
    maximum = float(maxvsini)
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise ValueError("maxvsini must be finite and positive.")
    if vsini is None:
        if vsini_window not in (None, 0, 0.0):
            raise ValueError("vsini_window requires a supplied vsini value.")
        mode = 'free'
        value = None
        window = None
        bounds = [0.0, maximum]
    else:
        value = float(vsini)
        if not np.isfinite(value) or not 0.0 <= value <= maximum:
            raise ValueError("vsini must be finite and between zero and maxvsini.")
        if vsini_window in (None, 0, 0.0):
            mode = 'fixed'
            window = 0.0
            bounds = [value, value]
        else:
            window = float(vsini_window)
            if not np.isfinite(window) or window < 0.0:
                raise ValueError("vsini_window must be finite and non-negative.")
            lower = max(0.0, value - window)
            upper = min(maximum, value + window)
            if not lower < upper:
                raise ValueError("The constrained vsini interval is empty.")
            mode = 'bounded'
            bounds = [lower, upper]
    return {
        'mode': mode,
        'value_kms': value,
        'max_kms': maximum,
        'window_kms': window,
        'effective_bounds_kms': bounds,
    }


def _build_estimator_config(
        *, order, wavelength, maxvsini, vsini=None, vsini_window=None,
        library_id=None, library_manifest=None, deblazed=False,
        blaze_source=None, reference_rv_state=None,
        reference_blaze_state=None):
    """Describe every runtime choice that changes the atmospheric estimator.

    Star-specific target filenames and target RVs are ordinary data
    provenance, not estimator settings.  The configuration binds both the
    rest-frame policy and the identity-keyed reference-library RV realization:
    every target/reference RV is established exactly once before the common-
    mask pairwise ranking and is never changed for the composite stage.
    """
    wavelength = np.asarray(wavelength, dtype='<f8')
    _scaled_chebyshev_coordinate(wavelength)
    manifest = library_manifest or {}
    wavelength_sha256 = hashlib.sha256(
        np.ascontiguousarray(wavelength).tobytes()
    ).hexdigest()
    return {
        'schema_version': ESTIMATOR_CONFIG_SCHEMA_VERSION,
        'spectral_grid': {
            'order': int(order),
            'wavelength_min_angstrom': float(wavelength[0]),
            'wavelength_max_angstrom': float(wavelength[-1]),
            'n_pixels': int(wavelength.size),
            'wavelength_float64_le_sha256': wavelength_sha256,
        },
        'library': {
            'library_id': None if library_id is None else str(library_id),
            'manifest_sha256': manifest.get('manifest_sha256'),
            'reference_domain': 'complete_authenticated_manifest_library',
        },
        'rv': {
            'input_domain': 'finite_absolute_rest_frame_rv',
            'application_policy': (
                'establish_target_and_reference_rvs_once_before_pairwise_'
                'ranking_and_reuse_without_recentering'
            ),
            'reference_library_state': reference_rv_state,
        },
        'blaze': {
            'input_is_deblazed': bool(deblazed),
            'target_blaze_source': blaze_source,
            'reference_policy': (
                'realized_per_reference_source;publication_validation_'
                'requires_all_l2'
            ),
            'reference_library_state': reference_blaze_state,
            'variance_policy': 'propagate_l2_variance_through_blaze_and_broadening',
        },
        'mask': {
            'policy': (
                'one_target_and_full_reference_library_intersection_frozen_'
                'at_vsini_support_limit_and_reused_by_composite'
            ),
            'minimum_common_pixels': MIN_PAIRWISE_PIXELS,
            'minimum_common_fraction': MIN_PAIRWISE_PIXEL_FRACTION,
            'gap_safe_resampling': True,
        },
        'pairwise': {
            'top_k': TOP_K,
            'unique_references': True,
            'selection_metric': 'unweighted_residual_sum_of_squares',
            'continuum_basis': 'chebyshev_scaled_to_minus1_plus1',
            'continuum_degree': CONTINUUM_DEGREE,
            'limb_darkening_u1': LIMB_DARKENING_U1,
            'solver': 'profiled_linear_coefficients_deterministic_1d_vsini',
        },
        'composite': {
            'top_k': TOP_K,
            'constraints': 'nonnegative_weights_sum_to_one',
            'solver': 'enumerated_active_set_weighted_least_squares',
            'loss': (
                'fixed_diagonal_target_errors_only_no_reference_error_or_'
                'resampling_covariance'
            ),
            'reuses_exact_pairwise_mask': True,
        },
        'vsini': _vsini_configuration(
            maxvsini=maxvsini, vsini=vsini, vsini_window=vsini_window
        ),
        'atmospheric_labels': {
            'method': 'weighted_empirical_reference_labels',
            'feh_scale': 'raw_no_empirical_detrending',
        },
    }

def _verbosity_level(verbose):
    """Normalize bool/int verbosity without treating ``True`` as a flood."""
    if isinstance(verbose, bool):
        return 1 if verbose else 0
    return max(0, int(verbose))


def _log(verbose, message, level=1):
    if _verbosity_level(verbose) >= level:
        print(message)


def _strict_bool_value(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered == 'true':
            return True
        if lowered == 'false':
            return False
    raise ValueError(f"Expected a strict boolean value, received {value!r}.")


@contextlib.contextmanager
def _quiet_external_calls(verbose):
    """Suppress stdout from legacy dependencies when quiet mode is requested."""
    if _verbosity_level(verbose) > 0:
        yield
    else:
        with contextlib.redirect_stdout(io.StringIO()):
            yield


def _set_spectrum_rv(spectrum, rv, source):
    """Apply one explicit RV and keep spectrum/provenance state consistent."""
    value = float(rv)
    if not np.isfinite(value):
        raise ValueError("Explicit RV must be finite.")
    spectrum.redshift(rv=value)
    # neidspec.redshift updates wavelength arrays, not necessarily these
    # bookkeeping attributes.  Set them explicitly for real and test objects.
    spectrum.rv = value
    spectrum.rv_source = source
    provenance = getattr(spectrum, 'provenance', None)
    if isinstance(provenance, dict):
        provenance['rv_source'] = source
        provenance['rv_override_kms'] = value


def _apply_matching_rv_overrides(target, references, *, absrv=None,
                                 reference_rvs=None):
    """Apply authorized RV overrides once, before validation and matching.

    Reference-library state receipts normally detect any mutation after the
    spectra are loaded.  An explicit RV override is an authorized,
    provenance-bearing state transition, so first verify every existing
    receipt, apply all requested shifts, and then refresh only those receipts.
    This prevents an override from blessing an unrelated prior mutation while
    allowing the subsequent manifest-state check to validate the exact state
    used by both matching stages.
    """
    references = list(references)
    requested = []
    if absrv is not None:
        requested.append((target, float(absrv)))
    if reference_rvs is not None:
        values = np.asarray(reference_rvs, dtype=float)
        if values.ndim != 1 or len(values) != len(references):
            raise ValueError(
                "reference_rvs must match the full reference list."
            )
        requested.extend(zip(references, values.tolist()))
    if any(not np.isfinite(value) for _, value in requested):
        raise ValueError("Explicit RV overrides must be finite.")

    # Check every pre-override receipt before mutating any spectrum, avoiding
    # a partial state transition if one reference was already modified.
    receipts = []
    for spectrum, _ in requested:
        receipt = getattr(
            spectrum, '_neidspecmatch_loaded_state_sha256', None
        )
        if (
            receipt is not None
            and receipt != _spectrum_state_fingerprint(spectrum)
        ):
            raise ValueError(
                "Loaded spectrum was mutated before its authorized RV "
                f"override: {getattr(spectrum, 'object', '<unknown>')}"
            )
        receipts.append(receipt)

    for (spectrum, value), receipt in zip(requested, receipts):
        _set_spectrum_rv(
            spectrum, value, "supplied to NEIDSpecMatch"
        )
        if receipt is not None:
            spectrum._neidspecmatch_loaded_state_sha256 = (
                _spectrum_state_fingerprint(spectrum)
            )


def _close_spectrum(value):
    close = getattr(value, 'close', None)
    if callable(close):
        close()


def _close_spectrum_collection(value):
    close = getattr(value, 'close', None)
    if callable(close):
        close()
        return
    for spectrum in getattr(value, 'splist', value):
        _close_spectrum(spectrum)


def _drp_version(spectrum):
    value = getattr(spectrum, 'drp_version', 'unknown')
    if value is None or not str(value).strip():
        return 'unknown'
    return str(value).strip()


def _drp_minor_series(version):
    """Return ``major.minor`` for a semantic NEID DRP version."""
    match = re.fullmatch(
        r"[vV]?(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)"
        r"\.(?P<patch>0|[1-9]\d*)(?:[-+][0-9A-Za-z.-]+)?",
        str(version).strip(),
    )
    if match is None:
        return None
    return "{}.{}".format(match.group("major"), match.group("minor"))


def _spectrum_quality_status(spectrum):
    """Normalize NEID DQ provenance without treating unrelated bit flags as fail."""
    direct = str(getattr(spectrum, 'dq_status', '') or '').strip().lower()
    if direct == 'pass':
        return 'pass'
    assessments = getattr(spectrum, 'dq_assessments', None)
    if not isinstance(assessments, dict):
        provenance = getattr(spectrum, 'provenance', {})
        assessments = provenance.get('dq_assessments') if isinstance(provenance, dict) else None
    statuses = []
    if isinstance(assessments, dict):
        statuses = [
            str(item.get('status', 'unknown')).strip().lower()
            for item in assessments.values() if isinstance(item, dict)
        ]
    for status in ('invalid', 'fail', 'warning', 'unknown'):
        if status in statuses:
            return status
    if direct in {'invalid', 'fail', 'warning', 'unknown', 'unvalidated'}:
        return direct
    return 'unknown'


def _quality_counts(spectra):
    return dict(sorted(Counter(
        _spectrum_quality_status(item) for item in spectra
    ).items()))


def _spectrum_quality_record(spectrum):
    provenance = getattr(spectrum, 'provenance', {})
    assessments = getattr(spectrum, 'dq_assessments', None)
    if assessments is None and isinstance(provenance, dict):
        assessments = provenance.get('dq_assessments')
    return {
        'dq_status': _spectrum_quality_status(spectrum),
        'dq_assessments': assessments if isinstance(assessments, dict) else None,
        'instrument': getattr(spectrum, 'instrument', provenance.get('instrument') if isinstance(provenance, dict) else None),
        'observation_type': getattr(spectrum, 'observation_type', provenance.get('observation_type') if isinstance(provenance, dict) else None),
        'data_level': getattr(spectrum, 'data_level', provenance.get('data_level') if isinstance(provenance, dict) else None),
        'observing_mode': getattr(spectrum, 'observing_mode', provenance.get('observing_mode') if isinstance(provenance, dict) else None),
    }


def _shareable_rv_source(value):
    text = str(value)
    match = re.search(r'mask=(.*?)\s+order_indices=', text)
    if match:
        mask_name = Path(match.group(1)).name
        text = text[:match.start(1)] + mask_name + text[match.end(1):]
    return text


_OMIT_FROM_SHAREABLE_PROVENANCE = object()


def _shareable_provenance_value(value, *, key=None):
    """Return a detached, JSON-safe value without retaining local paths."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, os.PathLike):
        return Path(os.fspath(value)).name
    if value is None or isinstance(value, (str, bool, int)):
        if (
            isinstance(value, str)
            and key is not None
            and str(key).lower().endswith(('_path', '_filename'))
        ):
            return Path(value).name
        return value
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        result = {}
        for nested_key, nested_value in value.items():
            safe_value = _shareable_provenance_value(
                nested_value, key=nested_key
            )
            if safe_value is not _OMIT_FROM_SHAREABLE_PROVENANCE:
                result[str(nested_key)] = safe_value
        return result
    if isinstance(value, (list, tuple)):
        result = []
        for item in value:
            safe_item = _shareable_provenance_value(item)
            if safe_item is not _OMIT_FROM_SHAREABLE_PROVENANCE:
                result.append(safe_item)
        return result
    return _OMIT_FROM_SHAREABLE_PROVENANCE


def _shareable_rv_diagnostics(target_provenance):
    """Copy all JSON-safe RV provenance while stripping machine-local paths."""
    diagnostics = {}
    if not isinstance(target_provenance, dict):
        return diagnostics
    for key, value in target_provenance.items():
        key = str(key)
        if not key.startswith('rv_'):
            continue
        if key == 'rv_source':
            value = _shareable_rv_source(value)
        output_key = key
        if key.endswith('_path'):
            output_key = key[:-5] + '_basename'
        safe_value = _shareable_provenance_value(value, key=key)
        if safe_value is not _OMIT_FROM_SHAREABLE_PROVENANCE:
            diagnostics[output_key] = safe_value
    return diagnostics


def validate_drp_compatibility(target, references, *, allow_mixed_drp=False):
    """Require one NEID DRP minor series unless explicitly waived.

    The NEID DRP team declares patch releases compatible within a minor
    series. Exact versions remain recorded for provenance, but for example
    v1.5.2 and v1.5.3 may be validated together; v1.4.x and v1.5.x may not.
    """
    versions = [_drp_version(target), *[_drp_version(item) for item in references]]
    counts = dict(sorted(Counter(versions).items()))
    series_values = [_drp_minor_series(version) for version in versions]
    unknown = any(value is None for value in series_values)
    series_counts = dict(sorted(Counter(
        value if value is not None else 'unknown' for value in series_values
    ).items()))
    mixed = len(series_counts) != 1
    mixed_patch_versions = len(counts) != 1 and not mixed and not unknown
    invalid = mixed or unknown
    if invalid and not allow_mixed_drp:
        raise ValueError(
            "Target/reference NEID DRP minor series are mixed or unknown: "
            "versions={}, series={}. "
            "Reprocess the target and library consistently or rerun with "
            "allow_mixed_drp=True only for explicitly unvalidated exploratory "
            "work.".format(counts, series_counts)
        )
    return {
        'status': (
            'mixed_or_unknown_unvalidated' if invalid else 'consistent'
        ),
        'versions': counts,
        'minor_series': series_counts,
        'mixed_drp': mixed,
        'mixed_patch_versions': mixed_patch_versions,
        'unknown_drp': unknown,
        'allow_mixed_drp': bool(allow_mixed_drp),
        'publication_validated': not invalid,
    }


def _file_identity(spectrum, dataset_id=None):
    filename = getattr(spectrum, 'filename', None)
    if not filename:
        return {
            'basename': None,
            'safe_identifier': safe_filename_component(
                getattr(spectrum, 'object', 'target')
            ),
            'sha256': None,
            'dataset_id': dataset_id,
        }
    path = Path(filename)
    digest = getattr(spectrum, '_neidspecmatch_sha256', None)
    if digest is None and path.is_file():
        digest = sha256_file(path)
        try:
            spectrum._neidspecmatch_sha256 = digest
        except Exception:
            pass
    return {
        'basename': path.name,
        'safe_identifier': safe_filename_component(path.name),
        'sha256': digest,
        'dataset_id': dataset_id,
    }


def _spectrum_state_fingerprint(spectrum):
    """Fingerprint loaded scientific arrays to detect post-load mutation."""
    digest = hashlib.sha256()
    for name in (
        'w', 'w_shifted', 'f', 'e', 'f_sci', 'f_sky', 'e_sci', 'e_sky',
        'error_scale', 'f_degrade',
    ):
        value = getattr(spectrum, name, None)
        if value is None:
            continue
        array = np.ascontiguousarray(np.asarray(value))
        digest.update(name.encode('utf-8'))
        digest.update(str(array.dtype).encode('ascii'))
        digest.update(json.dumps(array.shape).encode('ascii'))
        digest.update(array.view(np.uint8))
    for name in (
        'object', 'start', 'end', 'rv', 'rv_source', 'blaze_source',
        'legacy_response_path', 'sky_scaling_factor', 'degrade_snr',
    ):
        digest.update(name.encode('utf-8'))
        digest.update(repr(getattr(spectrum, name, None)).encode('utf-8'))
    return digest.hexdigest()


def _refresh_matching_blaze_state(target, references, *, input_is_deblazed):
    """Recreate every derived blaze array from its immutable L2 inputs.

    This is deliberately performed after manifest/state verification and
    before any ranking so a stale or externally modified deblaze cache cannot
    enter either matching stage.
    """
    references = list(references)
    if input_is_deblazed:
        target_source = 'input_predeblazed_hdu1_variance_hdu4'
    else:
        target.deblaze(force=True)
        target_source = getattr(target, 'blaze_source_used', None)
        if target_source is None or not str(target_source).strip():
            raise ValueError("Target deblazing did not record its blaze source.")
        target_source = str(target_source)
    for reference in references:
        reference.deblaze(force=True)
        source = getattr(reference, 'blaze_source_used', None)
        if source is None or not str(source).strip():
            raise ValueError(
                "Reference deblazing did not record its blaze source: "
                f"{getattr(reference, 'object', '<unknown>')}"
            )
    return target_source, _reference_blaze_state(references)


def _library_manifest_provenance(
        catalog_path, *, tolerate_invalid=False, deep=True,
        expected_library_id=None, validate_fits_schema=False):
    manifest_path = Path(catalog_path).resolve().parent / 'library_manifest.json'
    if not manifest_path.is_file():
        return {
            'manifest_present': False,
            'library_id': None,
            'manifest_sha256': None,
            'release_allowlist': None,
            'recognized': False,
            'deep_verified': False,
            'reference_pool_verified': False,
        }
    try:
        validation = package_utils.validate_library(
            manifest_path.parent, require_manifest=True, deep=deep,
            expected_library_id=expected_library_id,
            catalog_name=Path(catalog_path).name,
            validate_fits_schema=validate_fits_schema,
        )
        with manifest_path.open(encoding='utf-8') as stream:
            manifest = json.load(stream)
    except (OSError, ValueError, TypeError):
        if not tolerate_invalid:
            raise
        return {
            'manifest_present': True,
            'library_id': None,
            'manifest_sha256': None,
            'release_allowlist': None,
            'recognized': False,
            'deep_verified': False,
            'reference_pool_verified': False,
        }
    schema_version = manifest.get('schema_version')
    return {
        'manifest_present': True,
        'library_id': manifest.get('library_id'),
        'manifest_schema_version': schema_version,
        'manifest_sha256': sha256_file(manifest_path),
        'source_archive_integrity': manifest.get('source_archive_integrity'),
        'release_allowlist': manifest.get('release_allowlist'),
        'catalog': manifest.get('catalog'),
        'fits_count': manifest.get('fits_count'),
        'recognized': bool(validation.get('manifest') is not None),
        'deep_verified': bool(deep),
        'reference_pool_verified': False,
    }


def _verify_reference_pool_against_manifest(target, references, library_path):
    """Bind loaded reference objects to the exact deep-verified manifest files."""
    root = Path(library_path).resolve()
    manifest_path = root / 'library_manifest.json'
    with manifest_path.open(encoding='utf-8') as stream:
        manifest = json.load(stream)
    expected = {
        Path(record['path']).name: record
        for record in manifest.get('files', [])
        if PurePath(record['path']).parts[:1] == ('FITS',)
    }
    if not expected:
        raise ValueError("Library manifest contains no reference FITS records.")
    spectra = list(references)
    target_name = Path(str(getattr(target, 'filename', ''))).name
    if target_name in expected:
        spectra.append(target)
    provided = [
        Path(str(getattr(spectrum, 'filename', ''))).name
        for spectrum in spectra
    ]
    if any(not name for name in provided):
        raise ValueError(
            "Strict library validation requires a filename for every reference spectrum."
        )
    if len(provided) != len(set(provided)):
        raise ValueError("Loaded reference spectra contain duplicate FITS basenames.")
    if set(provided) != set(expected):
        missing = sorted(set(expected).difference(provided))
        extra = sorted(set(provided).difference(expected))
        raise ValueError(
            "Loaded reference pool does not match the library manifest "
            f"(missing={missing}, extra={extra})."
        )
    for spectrum in spectra:
        filename = Path(str(spectrum.filename)).resolve()
        basename = filename.name
        expected_filename = (root / 'FITS' / basename).resolve()
        if not filename.is_file():
            raise FileNotFoundError(f"Reference spectrum is missing: {basename}")
        # The library tree was already deep-hashed.  An external object with the
        # same basename must independently match the manifest digest.
        if filename != expected_filename:
            digest = sha256_file(filename)
            if digest != expected[basename]['sha256']:
                raise ValueError(
                    f"Reference spectrum SHA-256 does not match manifest: {basename}"
                )
        state_receipt = getattr(
            spectrum, '_neidspecmatch_loaded_state_sha256', None
        )
        if state_receipt is None:
            raise ValueError(
                "Publication validation requires reference spectra loaded by "
                "NEIDSpecMatch with immutable-state receipts."
            )
        if state_receipt != _spectrum_state_fingerprint(spectrum):
            raise ValueError(
                f"Loaded reference spectrum was mutated after loading: {basename}"
            )
    return {
        'reference_pool_verified': True,
        'reference_pool_fits_count': len(provided),
        'reference_pool_verification': (
            'exact_basenames_and_manifest_sha256_deep_verified'
        ),
    }


def _crossvalidation_evidence_status(
        raw_path, checkpoint_path, row, *, expected_order,
        expected_library_id, expected_manifest_sha, expected_blaze_source,
        expected_drp_version_set, expected_quality_counts,
        expected_estimator_json, expected_estimator_sha):
    """Reconstruct one summary row from its raw, fold-level evidence."""
    def invalid(reason):
        return {'valid': False, 'reason': 'crossvalidation_evidence_' + reason}

    def strict_integer(value):
        if isinstance(value, (bool, np.bool_)):
            raise ValueError
        number = float(value)
        if not np.isfinite(number) or not number.is_integer():
            raise ValueError
        return int(number)

    def parse_counts(value):
        if isinstance(value, str):
            value = json.loads(value)
        if not isinstance(value, dict):
            raise ValueError
        output = {}
        for key, count in value.items():
            parsed = strict_integer(count)
            if parsed < 0 or not str(key).strip():
                raise ValueError
            output[str(key)] = parsed
        return dict(sorted(output.items()))

    def optional_finite(value):
        number = pd.to_numeric(
            pd.Series([value]), errors='coerce'
        ).iloc[0]
        return float(number) if np.isfinite(number) else None

    def summary_value_matches(observed, expected):
        if expected is None:
            return bool(pd.isna(observed))
        if isinstance(expected, (str, bool, np.bool_)):
            return observed == expected
        try:
            expected_number = float(expected)
            observed_number = float(observed)
        except (TypeError, ValueError, OverflowError):
            return observed == expected
        if np.isnan(expected_number):
            return bool(np.isnan(observed_number))
        return bool(np.isclose(
            observed_number, expected_number, rtol=2e-12, atol=2e-12
        ))

    try:
        raw = pd.read_csv(raw_path)
    except (OSError, ValueError, pd.errors.ParserError):
        return invalid('raw_csv_unreadable')
    try:
        with Path(checkpoint_path).open(encoding='utf-8') as stream:
            checkpoint = json.load(stream)
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return invalid('checkpoint_unreadable')
    if not isinstance(checkpoint, dict):
        return invalid('checkpoint_not_a_mapping')
    required_checkpoint_fields = {
        'key', 'status', 'folds', 'result_csv', 'result_csv_sha256',
        'result_csv_hash_scope', 'optimizer',
    }
    if not required_checkpoint_fields.issubset(checkpoint):
        return invalid('checkpoint_incomplete')
    if checkpoint['status'] != 'complete':
        return invalid('checkpoint_not_complete')
    if (
        checkpoint['result_csv'] != Path(raw_path).name
        or checkpoint['result_csv_sha256'] != sha256_file(raw_path)
        or checkpoint['result_csv_sha256']
        != str(row['source_cv_results_sha256'])
        or checkpoint['result_csv_hash_scope'] != 'file_bytes'
    ):
        return invalid('checkpoint_raw_csv_binding_mismatch')

    key = checkpoint['key']
    if not isinstance(key, dict):
        return invalid('checkpoint_key_not_a_mapping')
    required_key_fields = {
        'checkpoint_schema_version', 'result_schema_version',
        'neidspecmatch_version', 'neidspec_version',
        'neidspec_source_sha256', 'pipeline_source_sha256',
        'selection_metric', 'order', 'library_id',
        'library_manifest_sha256', 'catalog_sha256', 'drp_versions',
        'blaze_source', 'estimator_config_schema_version',
        'estimator_config_sha256', 'estimator_config_json', 'random_seed',
        'spectrum_identities_in_order', 'reference_dq_status_counts',
    }
    if not required_key_fields.issubset(key):
        return invalid('checkpoint_key_incomplete')
    try:
        integer_key_values = {
            'checkpoint_schema_version': strict_integer(
                key['checkpoint_schema_version']
            ),
            'result_schema_version': strict_integer(
                key['result_schema_version']
            ),
            'order': strict_integer(key['order']),
            'estimator_config_schema_version': strict_integer(
                key['estimator_config_schema_version']
            ),
            'random_seed': strict_integer(key['random_seed']),
        }
    except (TypeError, ValueError, OverflowError):
        return invalid('checkpoint_key_integer_invalid')
    if (
        integer_key_values['checkpoint_schema_version'] != 4
        or integer_key_values['result_schema_version'] != RESULT_SCHEMA_VERSION
        or integer_key_values['order'] != int(expected_order)
        or integer_key_values['estimator_config_schema_version']
        != ESTIMATOR_CONFIG_SCHEMA_VERSION
        or key['neidspecmatch_version'] != __version__
        or key['neidspec_version']
        != getattr(neidspec, '__version__', 'unknown')
        or key['neidspec_source_sha256'] != _neidspec_source_fingerprint()
        or key['pipeline_source_sha256'] != _pipeline_source_fingerprint()
        or key['selection_metric'] != 'unweighted_residual_sum_of_squares'
        or str(key['library_id']) != str(expected_library_id)
        or key['library_manifest_sha256'] != expected_manifest_sha
        or key['blaze_source'] != expected_blaze_source
        or key['estimator_config_sha256'] != expected_estimator_sha
        or key['estimator_config_json'] != expected_estimator_json
    ):
        return invalid('checkpoint_key_runtime_mismatch')
    if re.fullmatch(r'[0-9a-f]{64}', str(key['catalog_sha256'])) is None:
        return invalid('checkpoint_catalog_hash_invalid')
    try:
        parsed_config = json.loads(key['estimator_config_json'])
    except (TypeError, ValueError, json.JSONDecodeError):
        return invalid('checkpoint_estimator_config_invalid')
    try:
        canonical_config = _canonical_estimator_config_json(parsed_config)
    except (TypeError, ValueError):
        return invalid('checkpoint_estimator_config_invalid')
    if (
        canonical_config != key['estimator_config_json']
        or _hash_estimator_config(parsed_config) != expected_estimator_sha
    ):
        return invalid('checkpoint_estimator_config_noncanonical')

    try:
        checkpoint_drp_counts = parse_counts(key['drp_versions'])
        summary_drp_counts = parse_counts(row['source_drp_versions'])
        checkpoint_quality_counts = parse_counts(
            key['reference_dq_status_counts']
        )
        summary_quality_counts = parse_counts(
            row['source_reference_dq_status_counts']
        )
        runtime_quality_counts = parse_counts(expected_quality_counts)
    except (TypeError, ValueError, OverflowError, json.JSONDecodeError):
        return invalid('checkpoint_composition_counts_invalid')
    if (
        checkpoint_drp_counts != summary_drp_counts
        or sorted(checkpoint_drp_counts) != sorted(expected_drp_version_set)
    ):
        return invalid('checkpoint_drp_composition_mismatch')
    if not (
        checkpoint_quality_counts == summary_quality_counts
        == runtime_quality_counts
    ):
        return invalid('checkpoint_dq_composition_mismatch')

    identities = key['spectrum_identities_in_order']
    folds = checkpoint['folds']
    optimizer = checkpoint['optimizer']
    if (
        not isinstance(identities, list)
        or not isinstance(folds, list)
        or not isinstance(optimizer, dict)
        or optimizer.get('all_folds_success') is not True
        or optimizer.get('failed_fold_indices') != []
        or len(identities) != len(raw)
        or len(folds) != len(raw)
    ):
        return invalid('checkpoint_fold_or_optimizer_proof_incomplete')
    try:
        _validate_checkpoint_folds(folds, identities)
    except (TypeError, ValueError, OverflowError):
        return invalid('checkpoint_fold_validation_failed')
    if [item.get('index') for item in folds] != list(range(len(raw))):
        return invalid('checkpoint_fold_order_incomplete')

    required_raw_columns = (
        'teff', 'feh', 'logg', 'vsini', 'd_teff', 'd_feh', 'd_logg',
        'teff_true', 'feh_true', 'logg_true', 'targetname',
    )
    if any(column not in raw.columns for column in required_raw_columns):
        return invalid('raw_csv_columns_incomplete')
    value_columns = required_raw_columns[:7]
    truth_columns = required_raw_columns[7:10]
    numeric = raw.loc[:, value_columns + truth_columns].apply(
        pd.to_numeric, errors='coerce'
    ).to_numpy(dtype=float)
    if not np.all(np.isfinite(numeric)):
        return invalid('raw_csv_contains_nonfinite_values')
    for index, (fold, identity) in enumerate(zip(folds, identities)):
        fold_optimizer = fold.get('optimizer')
        required_optimizer_fields = {
            'success', 'reported_success', 'finite_objective', 'message',
            'nfev', 'nit', 'objective_log_likelihood',
            'composite_nonunique_optimum',
            'pairwise_all_references_success',
            'pairwise_failed_reference_count',
        }
        if (
            not isinstance(identity, dict)
            or not isinstance(fold_optimizer, dict)
            or not required_optimizer_fields.issubset(fold_optimizer)
            or fold_optimizer['success'] is not True
            or fold_optimizer['reported_success'] is not True
            or fold_optimizer['finite_objective'] is not True
            or fold_optimizer['composite_nonunique_optimum'] is not False
            or fold_optimizer['pairwise_all_references_success'] is not True
        ):
            return invalid('fold_optimizer_receipt_incomplete')
        try:
            nfev = strict_integer(fold_optimizer['nfev'])
            nit = strict_integer(fold_optimizer['nit'])
            failed_references = strict_integer(
                fold_optimizer['pairwise_failed_reference_count']
            )
            objective = float(fold_optimizer['objective_log_likelihood'])
        except (TypeError, ValueError, OverflowError):
            return invalid('fold_optimizer_receipt_invalid')
        if (
            nfev < 0 or nit < 0 or failed_references != 0
            or not np.isfinite(objective)
            or not str(fold_optimizer['message']).strip()
        ):
            return invalid('fold_optimizer_receipt_invalid')
        identity_required = {
            'object_id', 'basename', 'file_sha256', 'loaded_state_sha256',
            'drp_version', 'dq_status', 'loaded_order_range',
        }
        if (
            not identity_required.issubset(identity)
            or not str(identity['object_id']).strip()
            or not str(identity['basename']).strip()
            or re.fullmatch(r'[0-9a-f]{64}', str(identity['file_sha256'])) is None
            or re.fullmatch(
                r'[0-9a-f]{64}', str(identity['loaded_state_sha256'])
            ) is None
            or not isinstance(identity['loaded_order_range'], list)
            or len(identity['loaded_order_range']) != 2
        ):
            return invalid('fold_identity_receipt_invalid')
        if (
            str(fold.get('targetname')) != str(raw.iloc[index]['targetname'])
            or str(fold.get('targetname')) != str(identity['object_id'])
        ):
            return invalid('fold_target_identity_mismatch')
        fold_values = np.asarray(fold['values'], dtype=float)
        raw_values = numeric[index, :7]
        if not np.allclose(
            fold_values, raw_values, rtol=2e-12, atol=2e-12
        ):
            return invalid('fold_values_do_not_match_raw_csv')
    identity_drp_counts = dict(sorted(Counter(
        str(item['drp_version']) for item in identities
    ).items()))
    identity_quality_counts = dict(sorted(Counter(
        str(item['dq_status']) for item in identities
    ).items()))
    if identity_drp_counts != checkpoint_drp_counts:
        return invalid('fold_identity_drp_counts_mismatch')
    if identity_quality_counts != checkpoint_quality_counts:
        return invalid('fold_identity_dq_counts_mismatch')

    expected_truth = np.column_stack((
        numeric[:, 0] - numeric[:, 4],
        numeric[:, 1] - numeric[:, 5],
        numeric[:, 2] - numeric[:, 6],
    ))
    truth_scale = np.maximum(1.0, np.abs(expected_truth))
    if not np.all(
        np.abs(numeric[:, 7:10] - expected_truth)
        <= 2e-12 * truth_scale
    ):
        return invalid('raw_truth_columns_inconsistent')

    if str(row['calibration_status']) != 'not_fitted':
        return invalid('calibrated_summary_not_strictly_reconstructible')
    try:
        population = PopulationRule(
            name=str(row['population']),
            teff_min=optional_finite(row['teff_min']),
            teff_max=optional_finite(row['teff_max']),
        )
        reconstructed = raw_crossvalidation_summary(
            raw,
            order=int(expected_order),
            population=population,
            library_id=str(expected_library_id),
            cv_source=Path(raw_path).name,
        )
    except (TypeError, ValueError, CalibrationError):
        return invalid('summary_reconstruction_failed')
    for field, expected in reconstructed.items():
        if field not in row.index:
            return invalid('summary_reconstructed_field_missing')
        if not summary_value_matches(row[field], expected):
            return invalid('summary_not_derived_from_raw_csv')

    population_mask = population.mask(raw).to_numpy(dtype=bool)
    supports = [item['support'] for item in folds]
    recomputed_support_counts = {
        'source_support_limited_folds': int(np.count_nonzero([
            mask and support['publication_support_validated'] is False
            for mask, support in zip(population_mask, supports)
        ])),
        'source_library_edge_folds': int(np.count_nonzero([
            mask and support['any_full_library_edge'] is True
            for mask, support in zip(population_mask, supports)
        ])),
        'source_simplex_vertex_folds': int(np.count_nonzero([
            mask and support['simplex']['is_vertex'] is True
            for mask, support in zip(population_mask, supports)
        ])),
    }
    for field, expected in recomputed_support_counts.items():
        try:
            observed = strict_integer(row[field])
        except (TypeError, ValueError, OverflowError):
            return invalid('summary_support_count_invalid')
        if observed != expected:
            return invalid('summary_support_count_not_derived_from_folds')
    return {
        'valid': True,
        'reason': 'raw_checkpoint_and_summary_reconstructed',
    }


def _validation_product_status(
        validation_summary, *, order, library_id, library_manifest,
        drp_compatibility, blaze_source, validation_population=None,
        validation_population_teff=None, calibration_metadata=None,
        reference_quality_counts=None, reference_drp_versions=None,
        estimator_config=None):
    calibration_metadata = calibration_metadata or {}
    if validation_summary is None:
        return {
            'matched': False,
            'reason': 'no_current_crossvalidation_product',
            'summary_sha256': None,
        }
    if estimator_config is None:
        return {
            'matched': False,
            'reason': 'runtime_estimator_configuration_missing',
            'summary_sha256': None,
        }
    expected_estimator_json = _canonical_estimator_config_json(
        estimator_config
    )
    expected_estimator_sha = _hash_estimator_config(estimator_config)
    summary_hash = None
    summary_path = None
    if isinstance(validation_summary, (str, os.PathLike, Path)):
        summary_path = Path(validation_summary)
        frame = pd.read_csv(summary_path)
        summary_hash = sha256_file(summary_path)
    elif isinstance(validation_summary, pd.DataFrame):
        frame = validation_summary.copy()
    elif isinstance(validation_summary, dict):
        frame = pd.DataFrame([validation_summary])
    else:
        raise TypeError("validation_summary must be a CSV path, DataFrame, or mapping")
    vsini_configuration = estimator_config.get('vsini')
    vsini_mode = (
        vsini_configuration.get('mode')
        if isinstance(vsini_configuration, dict) else None
    )
    if vsini_mode not in {'free', 'fixed', 'bounded'}:
        return {
            'matched': False,
            'reason': 'runtime_estimator_vsini_mode_missing_or_unrecognized',
            'summary_sha256': summary_hash,
        }
    if vsini_mode != 'free':
        return {
            'matched': False,
            'reason': (
                'fixed_or_bounded_vsini_requires_independent_'
                'broadlined_validation'
            ),
            'summary_sha256': summary_hash,
        }
    required = {
        'order', 'library_id', 'source_pipeline_version',
        'source_drp_compatibility_status', 'source_library_manifest_sha256',
        'source_neidspec_version', 'source_result_schema_version',
        'source_blaze_source', 'source_drp_versions',
        'source_drp_version_set', 'population', 'teff_min', 'teff_max',
        'error_metric_schema_version', 'residual_bias_statistic',
        'residual_scatter_statistic', 'residual_scatter_ddof',
        'predictive_error_statistic', 'predictive_error_ddof',
        'error_residual_method', 'uncertainty_semantics', 'n_stars',
        'teff_residual_bias', 'teff_sample_scatter',
        'teff_predictive_rmse', 'teff_n_valid', 'teff_finite_coverage',
        'feh_raw_residual_bias', 'feh_raw_sample_scatter',
        'feh_raw_predictive_rmse', 'feh_raw_n_valid',
        'feh_raw_finite_coverage', 'feh_calibrated_residual_bias',
        'feh_calibrated_sample_scatter', 'feh_calibrated_predictive_rmse',
        'feh_calibrated_n_valid', 'feh_calibrated_finite_coverage',
        'logg_residual_bias', 'logg_sample_scatter',
        'logg_predictive_rmse',
        'logg_n_valid', 'logg_finite_coverage',
        'calibration_status', 'feh_calibration_artifact_id',
        'feh_calibration_artifact_sha256', 'source_pipeline_sha256',
        'source_neidspec_source_sha256', 'source_cv_results_sha256',
        'source_cv_results_hash_scope', 'source_reference_dq_status_counts',
        'source_checkpoint', 'source_checkpoint_sha256',
        'source_checkpoint_hash_scope', 'cv_source',
        'source_all_folds_optimizer_success', 'source_failed_optimizer_folds',
        'source_estimator_config_schema_version',
        'source_estimator_config_sha256', 'source_estimator_config_json',
        'source_support_limited_folds', 'source_library_edge_folds',
        'source_simplex_vertex_folds',
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        return {
            'matched': False,
            'reason': 'validation_product_missing_fields:' + ','.join(missing),
            'summary_sha256': summary_hash,
        }
    expected_drp_versions = sorted(
        str(value) for value in (
            reference_drp_versions
            if reference_drp_versions is not None
            else drp_compatibility.get('versions', {})
        )
    )
    expected_drp_version_set = json.dumps(expected_drp_versions)
    expected_manifest_sha = library_manifest.get('manifest_sha256')
    matches = frame.loc[
        (pd.to_numeric(frame['order'], errors='coerce') == int(order))
        & (frame['library_id'].astype(str) == str(library_id))
        & (frame['source_pipeline_version'].astype(str) == __version__)
        & (frame['source_drp_compatibility_status'].astype(str) == 'consistent')
        & (frame['source_library_manifest_sha256'].astype(str)
           == str(expected_manifest_sha))
        & (frame['source_neidspec_version'].astype(str)
           == str(getattr(neidspec, '__version__', 'unknown')))
        & (frame['source_pipeline_sha256'].astype(str)
           == _pipeline_source_fingerprint())
        & (frame['source_neidspec_source_sha256'].astype(str)
           == _neidspec_source_fingerprint())
        & (pd.to_numeric(frame['source_result_schema_version'], errors='coerce')
           == RESULT_SCHEMA_VERSION)
        & (frame['source_blaze_source'].astype(str) == str(blaze_source))
        & (frame['source_drp_version_set'].astype(str)
           == expected_drp_version_set)
        & (
            pd.to_numeric(
                frame['source_estimator_config_schema_version'],
                errors='coerce',
            ) == ESTIMATOR_CONFIG_SCHEMA_VERSION
        )
        & (frame['source_estimator_config_sha256'].astype(str)
           == expected_estimator_sha)
        & (frame['source_estimator_config_json'].astype(str)
           == expected_estimator_json)
    ]
    if validation_population is None:
        lower = pd.to_numeric(matches['teff_min'], errors='coerce')
        upper = pd.to_numeric(matches['teff_max'], errors='coerce')
        matches = matches.loc[lower.isna() & upper.isna()]
    else:
        matches = matches.loc[
            matches['population'].astype(str) == str(validation_population)
        ]
    if len(matches) > 1:
        return {
            'matched': False,
            'reason': 'multiple_validation_population_rows_match',
            'summary_sha256': summary_hash,
        }
    if matches.empty:
        return {
            'matched': False,
            'reason': 'no_matching_current_order_library_drp_population_product',
            'summary_sha256': summary_hash,
        }
    row = matches.iloc[0]
    cv_results_hash = str(row['source_cv_results_sha256'])
    if (
        re.fullmatch(r'[0-9a-f]{64}', cv_results_hash) is None
        or str(row['source_cv_results_hash_scope']) != 'file_bytes'
    ):
        return {
            'matched': False,
            'reason': 'crossvalidation_source_hash_missing_or_unrecognized',
            'summary_sha256': summary_hash,
        }
    checkpoint_hash = str(row['source_checkpoint_sha256'])
    if (
        re.fullmatch(r'[0-9a-f]{64}', checkpoint_hash) is None
        or str(row['source_checkpoint_hash_scope']) != 'file_bytes'
    ):
        return {
            'matched': False,
            'reason': 'crossvalidation_checkpoint_hash_missing_or_unrecognized',
            'summary_sha256': summary_hash,
        }
    if summary_path is None:
        return {
            'matched': False,
            'reason': 'validation_sidecars_cannot_be_verified_from_in_memory_summary',
            'summary_sha256': summary_hash,
        }
    def local_regular_sidecar(value):
        name = str(value)
        if (
            not name or name in {'.', '..'} or '/' in name or '\\' in name
            or Path(name).name != name
        ):
            return None
        candidate = summary_path.parent / name
        if candidate.is_symlink() or not candidate.is_file():
            return None
        try:
            if candidate.resolve().parent != summary_path.parent.resolve():
                return None
        except OSError:
            return None
        return candidate

    raw_path = local_regular_sidecar(row['cv_source'])
    checkpoint_path = local_regular_sidecar(row['source_checkpoint'])
    if (
        raw_path is None
        or sha256_file(raw_path) != cv_results_hash
        or checkpoint_path is None
        or sha256_file(checkpoint_path) != checkpoint_hash
    ):
        return {
            'matched': False,
            'reason': 'crossvalidation_raw_or_checkpoint_sidecar_hash_mismatch',
            'summary_sha256': summary_hash,
        }
    expected_quality_counts = json.dumps(
        reference_quality_counts or {}, sort_keys=True
    )
    try:
        recorded_quality_counts = json.dumps(
            json.loads(str(row['source_reference_dq_status_counts'])),
            sort_keys=True,
        )
    except (TypeError, json.JSONDecodeError):
        recorded_quality_counts = None
    if recorded_quality_counts != expected_quality_counts:
        return {
            'matched': False,
            'reason': 'reference_dq_composition_does_not_match_validation_product',
            'summary_sha256': summary_hash,
        }
    try:
        optimizer_success = _strict_bool_value(
            row['source_all_folds_optimizer_success']
        )
    except ValueError:
        optimizer_success = False
    failed_optimizer_folds = pd.to_numeric(
        row['source_failed_optimizer_folds'], errors='coerce'
    )
    if not optimizer_success or failed_optimizer_folds != 0:
        return {
            'matched': False,
            'reason': 'crossvalidation_contains_nonconverged_optimizer_folds',
            'summary_sha256': summary_hash,
        }
    n_stars = pd.to_numeric(row['n_stars'], errors='coerce')
    if not np.isfinite(n_stars) or int(n_stars) < 8:
        return {
            'matched': False,
            'reason': 'validation_population_has_fewer_than_8_stars',
            'summary_sha256': summary_hash,
        }
    for field in (
        'source_support_limited_folds', 'source_library_edge_folds',
        'source_simplex_vertex_folds',
    ):
        count = pd.to_numeric(row[field], errors='coerce')
        if (
            not np.isfinite(count)
            or int(count) != count
            or not 0 <= int(count) <= int(n_stars)
        ):
            return {
                'matched': False,
                'reason': 'crossvalidation_support_diagnostics_incomplete',
                'summary_sha256': summary_hash,
            }
    metric_schema = pd.to_numeric(
        row['error_metric_schema_version'], errors='coerce'
    )
    scatter_ddof = pd.to_numeric(
        row['residual_scatter_ddof'], errors='coerce'
    )
    predictive_ddof = pd.to_numeric(
        row['predictive_error_ddof'], errors='coerce'
    )
    if (
        not np.isfinite(metric_schema)
        or int(metric_schema) != CROSSVALIDATION_METRIC_SCHEMA_VERSION
        or str(row['residual_bias_statistic']) != 'arithmetic_mean'
        or str(row['residual_scatter_statistic'])
        != 'sample_standard_deviation'
        or not np.isfinite(scatter_ddof) or int(scatter_ddof) != 1
        or str(row['predictive_error_statistic'])
        != 'root_mean_square_error'
        or not np.isfinite(predictive_ddof) or int(predictive_ddof) != 0
        or str(row['error_residual_method'])
        != CROSSVALIDATION_RESIDUAL_METHOD
        or str(row['uncertainty_semantics'])
        != NO_AUTOMATIC_UNCERTAINTY_SEMANTICS
    ):
        return {
            'matched': False,
            'reason': 'unrecognized_crossvalidation_error_metric_definition',
            'summary_sha256': summary_hash,
        }
    lower = pd.to_numeric(pd.Series([row['teff_min']]), errors='coerce').iloc[0]
    upper = pd.to_numeric(pd.Series([row['teff_max']]), errors='coerce').iloc[0]
    bounded = np.isfinite(lower) or np.isfinite(upper)
    if bounded:
        if validation_population_teff is None:
            return {
                'matched': False,
                'reason': 'bounded_population_requires_external_membership_teff',
                'summary_sha256': summary_hash,
            }
        membership_teff = float(validation_population_teff)
        if not np.isfinite(membership_teff):
            return {
                'matched': False,
                'reason': 'validation_population_membership_teff_must_be_finite',
                'summary_sha256': summary_hash,
            }
        if ((np.isfinite(lower) and membership_teff < lower)
                or (np.isfinite(upper) and membership_teff > upper)):
            return {
                'matched': False,
                'reason': 'target_not_in_selected_validation_population',
                'summary_sha256': summary_hash,
            }

    def finite_or_none(value):
        number = pd.to_numeric(pd.Series([value]), errors='coerce').iloc[0]
        return float(number) if np.isfinite(number) else None

    required_metric_prefixes = ['teff', 'feh_raw', 'logg']
    if calibration_metadata.get('status') == 'applied':
        required_metric_prefixes.append('feh_calibrated')
    for prefix in required_metric_prefixes:
        count = finite_or_none(row[f'{prefix}_n_valid'])
        coverage = finite_or_none(row[f'{prefix}_finite_coverage'])
        if (
            count is None or int(count) < 8 or int(count) != int(n_stars)
            or coverage is None or not np.isclose(coverage, 1.0)
        ):
            return {
                'matched': False,
                'reason': f'{prefix}_residuals_incomplete_or_insufficient',
                'summary_sha256': summary_hash,
            }
    error_metrics = {
        'schema_version': int(row['error_metric_schema_version']),
        'residual_bias': {
            'teff_k': finite_or_none(row['teff_residual_bias']),
            'feh_raw_dex': finite_or_none(row['feh_raw_residual_bias']),
            'feh_calibrated_dex': finite_or_none(
                row.get('feh_calibrated_residual_bias')
            ),
            'logg_dex': finite_or_none(row['logg_residual_bias']),
            'statistic': str(row['residual_bias_statistic']),
        },
        'sample_scatter': {
            'teff_k': finite_or_none(row['teff_sample_scatter']),
            'feh_raw_dex': finite_or_none(row['feh_raw_sample_scatter']),
            'feh_calibrated_dex': finite_or_none(
                row.get('feh_calibrated_sample_scatter')
            ),
            'logg_dex': finite_or_none(row['logg_sample_scatter']),
            'statistic': str(row['residual_scatter_statistic']),
            'ddof': int(row['residual_scatter_ddof']),
        },
        'predictive_rmse': {
            'teff_k': finite_or_none(row['teff_predictive_rmse']),
            'feh_raw_dex': finite_or_none(row['feh_raw_predictive_rmse']),
            'feh_calibrated_dex': finite_or_none(
                row.get('feh_calibrated_predictive_rmse')
            ),
            'logg_dex': finite_or_none(row['logg_predictive_rmse']),
            'statistic': str(row['predictive_error_statistic']),
            'ddof': int(row['predictive_error_ddof']),
        },
        'residual_method': str(row['error_residual_method']),
        'uncertainty_semantics': str(row['uncertainty_semantics']),
    }
    feh_metric = (
        'feh_calibrated_dex'
        if calibration_metadata.get('status') == 'applied'
        else 'feh_raw_dex'
    )
    relevant_metrics = []
    for metric_name in ('sample_scatter', 'predictive_rmse'):
        metric = error_metrics[metric_name]
        relevant_metrics.extend([
            metric['teff_k'], metric['logg_dex'], metric[feh_metric],
        ])
    if any(value is None or value <= 0 for value in relevant_metrics):
        return {
            'matched': False,
            'reason': 'crossvalidation_error_metrics_must_be_finite_and_positive',
            'summary_sha256': summary_hash,
        }
    evidence_status = _crossvalidation_evidence_status(
        raw_path,
        checkpoint_path,
        row,
        expected_order=order,
        expected_library_id=library_id,
        expected_manifest_sha=expected_manifest_sha,
        expected_blaze_source=str(blaze_source),
        expected_drp_version_set=expected_drp_versions,
        expected_quality_counts=reference_quality_counts or {},
        expected_estimator_json=expected_estimator_json,
        expected_estimator_sha=expected_estimator_sha,
    )
    if not evidence_status['valid']:
        return {
            'matched': False,
            'reason': evidence_status['reason'],
            'summary_sha256': summary_hash,
        }
    if calibration_metadata.get('status') == 'applied':
        calibration_rule = calibration_metadata.get('population_rule') or {}
        calibration_membership = calibration_metadata.get(
            'population_membership_teff'
        )
        calibration_bounds = (
            finite_or_none(calibration_rule.get('teff_min')),
            finite_or_none(calibration_rule.get('teff_max')),
        )
        validation_bounds = (
            finite_or_none(row['teff_min']), finite_or_none(row['teff_max'])
        )
        calibration_bounded = any(value is not None for value in calibration_bounds)
        membership_matches = (
            not calibration_bounded
            or (
                calibration_membership is not None
                and validation_population_teff is not None
                and np.isfinite(float(calibration_membership))
                and np.isclose(
                    float(calibration_membership),
                    float(validation_population_teff), rtol=0.0, atol=1e-9,
                )
            )
        )
        if (
            str(row['calibration_status']) != 'validated'
            or str(row['feh_calibration_artifact_id'])
            != str(calibration_metadata.get('artifact_id'))
            or str(row['feh_calibration_artifact_sha256'])
            != str(calibration_metadata.get('artifact_sha256'))
            or str(row['population']) != str(calibration_metadata.get('population'))
            or calibration_bounds != validation_bounds
            or not membership_matches
        ):
            return {
                'matched': False,
                'reason': 'calibration_artifact_not_bound_to_validation_row',
                'summary_sha256': summary_hash,
            }
    return {
        'matched': True,
        'reason': 'matched_current_product',
        'summary_sha256': summary_hash,
        'population': {
            'name': str(row['population']),
            'teff_min': finite_or_none(row['teff_min']),
            'teff_max': finite_or_none(row['teff_max']),
            'membership_teff': (
                float(validation_population_teff)
                if validation_population_teff is not None else None
            ),
            'selection': (
                'explicit_external_membership_teff' if bounded
                else 'unbounded_population'
            ),
        },
        'crossvalidation_error_metrics': error_metrics,
        'estimator_config_sha256': expected_estimator_sha,
        'evidence_validation': evidence_status,
        'fold_support_diagnostics': {
            'support_limited_folds': int(
                row['source_support_limited_folds']
            ),
            'library_edge_folds': int(row['source_library_edge_folds']),
            'simplex_vertex_folds': int(
                row['source_simplex_vertex_folds']
            ),
        },
    }


def build_result_provenance(
        *, target, references, rv_metadata, order, wavelength, deblazed,
        vsini, vsini_window, random_seed, feh_calibration_metadata,
        drp_compatibility, library_id=None, library_manifest=None,
        dataset_id=None, validation_summary=None,
        validation_population=None, validation_population_teff=None,
        composite_optimizer_diagnostics=None,
        pairwise_optimizer_diagnostics=None, estimator_config=None,
        support_diagnostics=None):
    """Build the schema-versioned, shareable provenance stored in result JSON."""
    target_provenance = getattr(target, 'provenance', {})
    rv_diagnostics = _shareable_rv_diagnostics(target_provenance)
    wavelength = np.asarray(wavelength, dtype=float)
    reference_versions = dict(sorted(Counter(
        _drp_version(item) for item in references
    ).items()))
    target_quality = _spectrum_quality_record(target)
    reference_quality_counts = _quality_counts(references)
    target_state_receipt = getattr(
        target, '_neidspecmatch_loaded_state_sha256', None
    )
    target_state_fingerprint = _spectrum_state_fingerprint(target)
    target_state_verified = bool(
        target_state_receipt is not None
        and target_state_receipt == target_state_fingerprint
    )
    calibration_population_details = feh_calibration_metadata.get('population')
    calibration = {
        'status': feh_calibration_metadata.get('status'),
        'artifact_id': feh_calibration_metadata.get('artifact_id'),
        'artifact_sha256': feh_calibration_metadata.get('artifact_sha256'),
        'order': feh_calibration_metadata.get('order'),
        'population': (
            calibration_population_details.get('name')
            if isinstance(calibration_population_details, dict)
            else calibration_population_details
        ),
        'population_rule': (
            calibration_population_details
            if isinstance(calibration_population_details, dict) else None
        ),
        'population_membership_teff': feh_calibration_metadata.get(
            'population_membership_teff'
        ),
        'validated': feh_calibration_metadata.get('validated'),
    }
    manifest = library_manifest or {'manifest_present': False}
    blaze_source = (
        'input_predeblazed_hdu1_variance_hdu4'
        if deblazed else getattr(
            target, 'blaze_source_used', target_provenance.get('blaze_source')
        )
    )
    reference_blaze_sources = dict(sorted(Counter(
        str(getattr(item, 'blaze_source_used', 'unknown'))
        for item in references
    ).items()))
    validation_product = _validation_product_status(
        validation_summary, order=order, library_id=library_id,
        library_manifest=manifest, drp_compatibility=drp_compatibility,
        blaze_source=blaze_source,
        validation_population=validation_population,
        validation_population_teff=validation_population_teff,
        calibration_metadata=calibration,
        reference_quality_counts=reference_quality_counts,
        reference_drp_versions=reference_versions,
        estimator_config=estimator_config,
    )
    validation_reasons = []
    if not drp_compatibility.get('publication_validated', False):
        validation_reasons.append('mixed_or_unknown_drp_versions')
    if (
        str(blaze_source).lower() != 'l2'
        or any(str(source).lower() != 'l2'
               for source in reference_blaze_sources)
    ):
        validation_reasons.append('non_l2_blaze_treatment_unvalidated')
    if not manifest.get('manifest_present') or not manifest.get('recognized'):
        validation_reasons.append('library_manifest_absent_or_unrecognized')
    elif str(manifest.get('library_id')) != str(library_id):
        validation_reasons.append('library_manifest_id_mismatch')
    elif not manifest.get('deep_verified'):
        validation_reasons.append('library_manifest_not_deep_hash_verified')
    elif not manifest.get('reference_pool_verified'):
        validation_reasons.append('reference_pool_not_bound_to_library_manifest')
    if not validation_product['matched']:
        validation_reasons.append(validation_product['reason'])
    if calibration['status'] == 'applied' and calibration['validated'] is not True:
        validation_reasons.append('feh_calibration_not_validated')
    if calibration['status'] == 'applied':
        validation_reasons.append(
            'feh_calibration_crossfit_not_strictly_nested'
        )
    if target_quality['dq_status'] != 'pass':
        validation_reasons.append('target_dq_status_not_pass')
    if str(target_quality.get('instrument')).upper() != 'NEID':
        validation_reasons.append('target_instrument_not_confirmed_neid')
    if str(target_quality.get('observation_type')).upper() != 'SCI':
        validation_reasons.append('target_observation_type_not_confirmed_science')
    if target_quality.get('data_level') != 2:
        validation_reasons.append('target_data_level_not_confirmed_l2')
    if str(target_quality.get('observing_mode')).upper() != 'HR':
        validation_reasons.append('target_mode_not_confirmed_hr')
    if not target_state_verified:
        validation_reasons.append('target_loaded_state_not_verified')
    optimizer_diagnostics = composite_optimizer_diagnostics or {}
    if optimizer_diagnostics.get('success') is not True:
        validation_reasons.append('composite_optimizer_nonconverged')
    if optimizer_diagnostics.get('nonunique_optimum') is True:
        validation_reasons.append('composite_optimizer_nonunique_solution')
    pairwise_diagnostics = pairwise_optimizer_diagnostics or {}
    if (
        not isinstance(pairwise_diagnostics, dict)
        or pairwise_diagnostics.get('all_references_success') is not True
    ):
        validation_reasons.append('pairwise_optimizer_failure')
    support = support_diagnostics if isinstance(support_diagnostics, dict) else {}
    support_status = _support_diagnostics_status(support)
    if support_status['valid']:
        support_parameters = support['parameters']
        support_reasons = [
            reason for reason in support['warnings']
            if reason.endswith('_full_reference_library_lower_boundary')
            or reason.endswith('_full_reference_library_upper_boundary')
            or reason == 'simplex_vertex_at_full_reference_library_boundary'
        ]
    else:
        support_parameters = {}
        support_reasons = []
        validation_reasons.append(
            'reference_support_diagnostics_invalid:'
            + support_status['reason']
        )
    base_atmospheric_validated = not validation_reasons
    error_metrics = validation_product.get(
        'crossvalidation_error_metrics', {}
    )
    residual_bias = error_metrics.get('residual_bias', {})
    sample_scatter = error_metrics.get('sample_scatter', {})
    predictive_rmse = error_metrics.get('predictive_rmse', {})
    feh_metric_name = (
        'feh_calibrated_dex'
        if calibration['status'] == 'applied' else 'feh_raw_dex'
    )

    def parameter_support(name):
        details = support_parameters.get(name, {})
        return details, bool(details.get('support_limited', False))

    teff_support, teff_limited = parameter_support('teff')
    feh_support, feh_limited = parameter_support('feh_raw')
    logg_support, logg_limited = parameter_support('logg')
    parameter_validation = {
        'teff': {
            'publication_validated': (
                base_atmospheric_validated and not teff_limited
            ),
            'reference_support_limited': teff_limited,
            'reference_support': teff_support,
            'crossvalidation_residual_bias_k': residual_bias.get('teff_k'),
            'crossvalidation_sample_scatter_k': sample_scatter.get('teff_k'),
            'crossvalidation_predictive_rmse_k': predictive_rmse.get('teff_k'),
        },
        'feh': {
            'publication_validated': (
                base_atmospheric_validated and not feh_limited
            ),
            'reference_support_limited': feh_limited,
            'reference_support': feh_support,
            'crossvalidation_residual_bias_dex': residual_bias.get(
                feh_metric_name
            ),
            'crossvalidation_sample_scatter_dex': sample_scatter.get(
                feh_metric_name
            ),
            'crossvalidation_predictive_rmse_dex': predictive_rmse.get(
                feh_metric_name
            ),
            'estimator': 'calibrated' if calibration['status'] == 'applied' else 'raw',
        },
        'logg': {
            'publication_validated': (
                base_atmospheric_validated and not logg_limited
            ),
            'reference_support_limited': logg_limited,
            'reference_support': logg_support,
            'crossvalidation_residual_bias_dex': residual_bias.get('logg_dex'),
            'crossvalidation_sample_scatter_dex': sample_scatter.get('logg_dex'),
            'crossvalidation_predictive_rmse_dex': predictive_rmse.get('logg_dex'),
        },
        'vsini': {
            'publication_validated': False,
            'reason': 'no_independent_reference_truth_in_crossvalidation',
            'crossvalidation_residual_bias_kms': None,
            'crossvalidation_sample_scatter_kms': None,
            'crossvalidation_predictive_rmse_kms': None,
        },
    }
    atmospheric_validated = all(
        parameter_validation[name]['publication_validated']
        for name in ('teff', 'feh', 'logg')
    )
    overall_reasons = [
        *validation_reasons, *support_reasons,
        'vsini_no_independent_reference_truth',
    ]
    return {
        'schema_version': RESULT_SCHEMA_VERSION,
        'software': {
            'neidspecmatch': __version__,
            'neidspec': getattr(neidspec, '__version__', 'unknown'),
            'neidspecmatch_source_sha256': _pipeline_source_fingerprint(),
            'neidspec_source_sha256': _neidspec_source_fingerprint(),
            'python': sys.version.split()[0],
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'astropy': astropy.__version__,
            'pandas': pd.__version__,
        },
        'target': {
            'name': str(target.object),
            'file': _file_identity(target, dataset_id=dataset_id),
            'drp_version': _drp_version(target),
            'data_level': getattr(target, 'data_level', None),
            'observatory': getattr(target, 'observatory', None),
            'telescope': getattr(target, 'telescope', None),
            'blaze': {
                'input_is_deblazed': bool(deblazed),
                'source': blaze_source,
                'drift_corrected': getattr(target, 'drift_corrected', None),
            },
            'rv': rv_metadata['target'],
            'rv_diagnostics': rv_diagnostics,
            'quality': target_quality,
            'loaded_state': {
                'receipt_sha256': target_state_receipt,
                'current_sha256': target_state_fingerprint,
                'verified': target_state_verified,
            },
            'resampling': target_provenance.get('last_resampling'),
            'loaded_order_range': {
                'first_order_index': getattr(target, 'start', None),
                'stop_order_index_exclusive': getattr(target, 'end', None),
            },
        },
        'reference_pool': {
            'count': len(references),
            'drp_versions': reference_versions,
            'blaze_sources': reference_blaze_sources,
            'dq_status_counts': reference_quality_counts,
            'loaded_order_ranges': sorted({
                (
                    getattr(item, 'start', None),
                    getattr(item, 'end', None),
                )
                for item in references
            }, key=str),
        },
        'selected_components': {
            'count': len(rv_metadata['references']),
            'rv': rv_metadata['references'],
            'fit': support.get('selected_components', []),
        },
        'drp_compatibility': drp_compatibility,
        'library': {
            'library_id': library_id,
            'manifest': manifest,
        },
        'spectral_grid': {
            'order': int(order),
            'wavelength_min_angstrom': float(np.min(wavelength)),
            'wavelength_max_angstrom': float(np.max(wavelength)),
            'wavelength_step_angstrom': float(np.median(np.diff(wavelength))),
            'n_pixels': int(wavelength.size),
        },
        'estimator_configuration': {
            'configuration': estimator_config,
            'sha256': (
                _hash_estimator_config(estimator_config)
                if estimator_config is not None else None
            ),
        },
        'vsini_constraint': {
            'value_kms': None if vsini is None else float(vsini),
            'window_kms': None if vsini_window is None else float(vsini_window),
            'mode': 'free' if vsini is None else (
                'fixed' if vsini_window in (None, 0, 0.0) else 'bounded'
            ),
            'max_kms': (
                estimator_config.get('vsini', {}).get('max_kms')
                if isinstance(estimator_config, dict) else None
            ),
        },
        'optimization': {
            'pairwise_stage': {
                'name': 'profiled_linear_chebyshev_plus_bounded_1d_vsini',
                'selection_metric': 'unweighted_residual_sum_of_squares',
                'uncertainties_used_for_selection': False,
                'continuum_basis': 'degree_5_Chebyshev_on_order_scaled_minus1_plus1',
                'selected_reference_diagnostics': (
                    pairwise_diagnostics.get('selected_references', [])
                ),
                'all_reference_summary': pairwise_diagnostics,
                'comparison_mask_policy': (
                    'one_target_and_all_reference_intersection_frozen_at_maximum_vsini'
                ),
                'minimum_common_pixels': MIN_PAIRWISE_PIXELS,
                'minimum_common_fraction': MIN_PAIRWISE_PIXEL_FRACTION,
            },
            'composite_stage': {
                'name': 'enumerated_active_set_weighted_least_squares',
                'seed': None,
                'diagnostics': optimizer_diagnostics,
                'loss_model': (
                    'fixed_diagonal_target_errors_only; reference_errors_and_'
                    'resampling_broadening_covariance_not_in_objective'
                ),
            },
        },
        'feh_calibration': calibration,
        'reference_support': support,
        'reference_support_validation': support_status,
        'crossvalidation_product': validation_product,
        'parameter_validation': parameter_validation,
        'overall_validation': {
            'status': (
                'atmospheric_parameters_validated_vsini_unvalidated'
                if atmospheric_validated else 'exploratory_unvalidated'
            ),
            'publication_validated': False,
            'atmospheric_parameters_publication_validated': atmospheric_validated,
            'reasons': overall_reasons,
        },
    }


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return os.fspath(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _pipeline_source_fingerprint():
    """Hash estimator sources and the numerical-library runtime identity.

    Exact runtime versions are recorded explicitly in every result receipt.
    NumPy, SciPy, Astropy, and pandas are also folded into this digest because
    they participate in the numerical result.  Python itself is not: the
    release tests one locked numerical stack on every supported Python minor.
    Publication-figure code is excluded because it cannot change a fit and
    records its own source hash in every figure provenance receipt.
    """
    digest = hashlib.sha256()
    numerical_versions = {
        'numpy': np.__version__,
        'scipy': scipy.__version__,
        'astropy': astropy.__version__,
        'pandas': pd.__version__,
    }
    digest.update(json.dumps(
        numerical_versions, sort_keys=True, separators=(',', ':')
    ).encode('utf-8'))
    digest.update(b'\0')
    package_directory = Path(__file__).resolve().parent
    for path in sorted(package_directory.rglob('*.py')):
        name = path.relative_to(package_directory).as_posix()
        if name == 'crossvalidation_figure.py':
            continue
        digest.update(name.encode('utf-8'))
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


def _neidspec_source_fingerprint():
    """Hash the imported NEIDSpec Python implementation, not only its version."""
    package_directory = Path(neidspec.__file__).resolve().parent
    digest = hashlib.sha256()
    for path in sorted(package_directory.rglob('*.py')):
        name = path.relative_to(package_directory).as_posix()
        digest.update(name.encode('utf-8'))
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


def _validate_calibration_runtime_context(
        calibration, *, library_id, library_manifest, drp_compatibility,
        blaze_source):
    """Fail before applying calibration evidence from a different pipeline."""
    artifact = coerce_calibration_set(calibration)
    if not artifact.source_evidence_validated:
        raise UnvalidatedCalibrationError(
            "Calibration artifact lacks complete successful source evidence."
        )
    if not drp_compatibility.get('publication_validated', False):
        raise UnvalidatedCalibrationError(
            "Calibrations cannot be applied to mixed or unknown DRP data."
        )
    expected_drp = tuple(sorted(str(value) for value in drp_compatibility['versions']))
    checks = {
        'library manifest SHA-256': (
            artifact.source_library_manifest_sha256,
            library_manifest.get('manifest_sha256'),
        ),
        'matcher pipeline SHA-256': (
            artifact.source_pipeline_sha256, _pipeline_source_fingerprint(),
        ),
        'NEIDSpec source SHA-256': (
            artifact.source_neidspec_sha256, _neidspec_source_fingerprint(),
        ),
        'NEIDSpec version': (
            artifact.source_neidspec_version,
            getattr(neidspec, '__version__', 'unknown'),
        ),
        'result schema': (
            artifact.source_result_schema_version, RESULT_SCHEMA_VERSION,
        ),
        'blaze source': (artifact.source_blaze_source, blaze_source),
        'DRP version set': (artifact.source_drp_version_set, expected_drp),
    }
    mismatches = [name for name, values in checks.items() if values[0] != values[1]]
    if mismatches:
        raise UnvalidatedCalibrationError(
            "Calibration evidence does not match this run: " + ', '.join(mismatches)
        )
    if not any(item.library_id == str(library_id) for item in artifact.calibrations):
        raise UnvalidatedCalibrationError(
            "Calibration artifact is not bound to this library_id."
        )
    return artifact


def _scaled_chebyshev_coordinate(wavelength):
    wavelength = np.asarray(wavelength, dtype=float)
    if (
        wavelength.ndim != 1
        or wavelength.size < 2
        or not np.all(np.isfinite(wavelength))
        or np.any(np.diff(wavelength) <= 0.0)
    ):
        raise ValueError(
            "Wavelength grid must be finite, one-dimensional, and strictly increasing."
        )
    midpoint = 0.5 * (wavelength[0] + wavelength[-1])
    half_range = 0.5 * (wavelength[-1] - wavelength[0])
    return (wavelength - midpoint) / half_range


def _prepare_reference_model(
        wavelength, flux, error, coefficients, vsini,
        *, limb_darkening=LIMB_DARKENING_U1):
    """Apply the first-stage continuum/vsini model identically for stage two."""
    coordinate = _scaled_chebyshev_coordinate(wavelength)
    continuum = np.polynomial.chebyshev.chebval(coordinate, coefficients)
    if not np.all(np.isfinite(continuum)):
        raise ValueError("Chebyshev continuum evaluated to non-finite values.")
    transformed_flux = np.asarray(flux, dtype=float) * continuum
    transformed_error = np.asarray(error, dtype=float) * np.abs(continuum)
    transformed_flux = rotbroad_help.broaden(
        np.asarray(wavelength, dtype=float), transformed_flux, float(vsini),
        u1=float(limb_darkening),
    )
    transformed_variance = rotbroad_help.broaden_variance(
        np.asarray(wavelength, dtype=float), transformed_error ** 2,
        float(vsini), u1=float(limb_darkening),
    )
    return transformed_flux, np.sqrt(transformed_variance)


def get_data_ready(H1, Hrefs, w, v, polyvals=None, vsinis=None, order=101,
                   plot=False, deblazed=False, verbose=1, absrv=None,
                   reference_rvs=None, return_rv_metadata=False):
    """
    Get data ready for MCMC
    
    INPUT:
        H1 - target spectrum (HPFSpectrum object)
        Hrefs - reference spectra (HPFSpectraList object)
        w - wavelength grid to interpolate on (array)
        v - velocities in km/s to use for absolute RV consideration (array)
        polyvals - polynomial coefficients (array)
        vsinis - vsini values km/s (array)
        plot - (boolean)

    OUTPUT:
        f1 - spectrum
        e1 - spectrum flux error
        ffrefs - spectra  (array)
        eerefs - spectra flux errors (array)
       
    EXAMPLE:
        "Target={}, rv={:0.3f}km/s, rvmed={:0.3f}km/s".format(H1.target.name,H1.rv,np.median(rabs))
        files = sorted(glob.glob('20200209_ad_leos/AD_Leo/*/*.pkl'))
        summarize_values_from_orders(files,'AD_Leo')
    """
    if absrv is not None:
        with _quiet_external_calls(verbose):
            _set_spectrum_rv(H1, absrv, "supplied to NEIDSpecMatch")
        target_rv_source = "supplied to NEIDSpecMatch"
    else:
        if not hasattr(H1, "rv") or not np.isfinite(H1.rv):
            raise ValueError(
                "Target spectrum has no finite RV. Initialize NEIDSpectrum with "
                "rv_source='drp' or supply absrv explicitly."
            )
        target_rv_source = _shareable_rv_source(
            getattr(H1, "rv_source", "spectrum object")
        )
    target_rv = float(H1.rv)
    with _quiet_external_calls(verbose):
        if not deblazed:
            H1.deblaze()
        f1, e1 = H1.resample_order(w, plot=plot, order=order, deblazed=deblazed)
    _log(verbose, "Target={}, rv={:0.3f} km/s ({})".format(
        H1.object, target_rv, target_rv_source), level=2)

    ffrefs = []
    eerefs = []
    if reference_rvs is not None and len(reference_rvs) != len(Hrefs.splist):
        raise ValueError("reference_rvs must have one value per reference spectrum.")
    reference_metadata = []
    for i, H in enumerate(Hrefs.splist):
        if reference_rvs is not None:
            reference_rv = float(reference_rvs[i])
            reference_source = "supplied to NEIDSpecMatch"
            with _quiet_external_calls(verbose):
                _set_spectrum_rv(H, reference_rv, reference_source)
        else:
            if not hasattr(H, "rv") or not np.isfinite(H.rv):
                raise ValueError(
                    "Reference {!r} has no finite RV; initialize the library "
                    "spectrum with an explicit RV source.".format(H.object)
                )
            reference_rv = float(H.rv)
            reference_source = _shareable_rv_source(
                getattr(H, "rv_source", "spectrum object")
            )
        with _quiet_external_calls(verbose):
            H.deblaze()
            if polyvals is None and vsinis is None:
                _f, _e = H.resample_order(w, order=order)
            else:
                if polyvals is None or vsinis is None:
                    raise ValueError(
                        "polyvals and vsinis must be supplied together for references."
                    )
                _f, _e = H.resample_order(w, order=order)
                _f, _e = _prepare_reference_model(
                    w, _f, _e, polyvals[i], vsinis[i]
                )
        ffrefs.append(_f)
        eerefs.append(_e)
        reference_metadata.append({
            "object": H.object,
            "rv_kms": reference_rv,
            "rv_source": reference_source,
        })
        _log(verbose, "Reference={}, rv={:0.3f} km/s ({})".format(
            H.object, reference_rv, reference_source), level=2)
    # plt.show()
    rv_metadata = {
        "target": {
            "object": H1.object,
            "rv_kms": target_rv,
            "rv_source": target_rv_source,
        },
        "references": reference_metadata,
        "policy": "preserve_neidspec_rv_unless_explicitly_supplied",
    }
    if return_rv_metadata:
        return f1, e1, ffrefs, eerefs, rv_metadata
    return f1, e1, ffrefs, eerefs


class LPFunctionLinComb(object):
    """
    Log Likelihood function for optimizing a combination of 5 top spectra
    """

    def __init__(self, w, f1, e1, f2s, e2s, valid_mask=None):
        self.w = w
        self.data_target = {'f': f1,
                            'e': e1}
        self.data_refs = {'f': f2s,
                          'e': e2s}
        self.num_refs = len(self.data_refs['f'])
        target = np.asarray(f1, dtype=float)
        error = np.asarray(e1, dtype=float)
        references = np.asarray(f2s, dtype=float)
        if references.shape != (self.num_refs, target.size):
            raise ValueError("Composite reference array has an invalid shape.")
        finite_mask = (
            np.isfinite(target) & np.isfinite(error) & (error > 0)
            & np.all(np.isfinite(references), axis=0)
        )
        if valid_mask is None:
            self.valid_mask = finite_mask
        else:
            requested_mask = np.asarray(valid_mask, dtype=bool)
            if requested_mask.shape != target.shape:
                raise ValueError(
                    "Composite valid_mask must match the target spectrum."
                )
            if np.any(requested_mask & ~finite_mask):
                raise ValueError(
                    "Composite valid_mask includes unsupported target/reference "
                    "pixels."
                )
            self.valid_mask = requested_mask.copy()
        if np.count_nonzero(self.valid_mask) < self.num_refs:
            raise ValueError("Too few common valid pixels for composite fitting.")

        self.priors = [UP(0., 1., 'c1', '$c1$', priortype="model"),
                       UP(0., 1., 'c2', '$c2$', priortype="model"),
                       UP(0., 1., 'c3', '$c3$', priortype="model"),
                       UP(0., 1., 'c4', '$c4$', priortype="model")]
        self.ps = PriorSet(self.priors)

    def get_pv_all(self, pv):
        """Return a canonical full simplex vector.

        The public parameterization stores four weights and derives the fifth.
        Floating-point summation can make an exactly active boundary appear as
        ``-2e-16``.  Roundoff at that scale is canonicalized, while material
        simplex violations still fail explicitly.
        """
        partial = np.asarray(pv, dtype=float)
        if partial.shape != (self.num_refs - 1,) or not np.all(np.isfinite(partial)):
            raise ValueError("Composite parameters must be finite and have length four.")
        weights = np.r_[partial, 1.0 - float(np.sum(partial))]
        tolerance = 1e-12
        if np.min(weights) < -tolerance or np.max(weights) > 1.0 + tolerance:
            raise ValueError("Composite parameters lie outside the simplex.")
        weights[np.abs(weights) <= tolerance] = 0.0
        weights[np.abs(weights - 1.0) <= tolerance] = 1.0
        total = float(np.sum(weights))
        if not np.isfinite(total) or total <= 0.0:
            raise ValueError("Composite weights do not have a positive finite sum.")
        weights /= total
        # Put the final floating-point closure error on the largest component;
        # this keeps exact zeros exact and makes the full vector sum to unity.
        closure = 1.0 - float(np.sum(weights))
        weights[int(np.argmax(weights))] += closure
        return weights

    def compute_model(self, pv):
        pv_all = self.get_pv_all(pv)
        # Do not evaluate zero-weight components: IEEE arithmetic makes
        # ``0 * NaN`` a NaN, so unsupported pixels in an inactive reference
        # would otherwise create artificial gaps in the saved composite.
        # The fit itself remains restricted to the frozen all-reference
        # common mask constructed in ``__init__``.
        active = pv_all > 0.0
        references = np.asarray(self.data_refs['f'], dtype=float)
        return pv_all[active] @ references[active]

    def __call__(self, pv):
        try:
            weights = self.get_pv_all(pv)
        except ValueError:
            return -np.inf
        flux_model = np.asarray(self.data_refs['f'], dtype=float).T @ weights
        flux_target = np.asarray(self.data_target['f'])[self.valid_mask]
        error_target = np.asarray(self.data_target['e'])[self.valid_mask]
        flux_model = np.asarray(flux_model)[self.valid_mask]
        log_of_model = ll_normal_ev_py(flux_target, flux_model, error_target)
        # The simplex prior is uniform including its boundary.  Historical
        # UniformPrior used open intervals, spuriously penalizing exact zero
        # weights even though the explicit constraints allow them.
        return log_of_model


class FitLinCombSpec(object):
    """
    A class to fit 5 lin-comb spectra together. Note: look at LPFunctionLinComb
    """

    def __init__(self, LPFunctionLinComb, teffs=None, fehs=None, loggs=None, vsinis=None,
                 teff_known=None, tefferr_known=None,
                 feh_known=None, feherr_known=None,
                 logg_known=None, loggerr_known=None, targetname='',
                 calibrate_feh=False, feh_calibration=None, order=None,
                 library_id=None, calibration_population=None,
                 calibration_population_teff=None, verbose=1):
        self.lpf = LPFunctionLinComb
        self.teffs = [] if teffs is None else teffs
        self.fehs = [] if fehs is None else fehs
        self.loggs = [] if loggs is None else loggs
        self.vsinis = [] if vsinis is None else vsinis
        self.teff_known = teff_known
        self.tefferr_known = tefferr_known
        self.feh_known = feh_known
        self.feherr_known = feherr_known
        self.logg_known = logg_known
        self.loggerr_known = loggerr_known
        self.targetname = targetname
        self.calibrate_feh = calibrate_feh
        self.feh_calibration = feh_calibration
        self.order = order
        self.library_id = library_id
        self.calibration_population = calibration_population
        self.calibration_population_teff = calibration_population_teff
        self.verbose = verbose
        self.feh_calibration_metadata = {
            'status': 'not_applied',
            'reason': 'safe_default_raw_metallicity',
        }

    def calculate_stellar_parameters(self, weights):
        if len(self.teffs) > 0:
            self.teff = weighted_value(self.teffs, weights)
        else:
            self.teff = np.nan
        if len(self.fehs) > 0:
            self.feh_raw = weighted_value(self.fehs, weights)
            self.feh = self.feh_raw
            if self.calibrate_feh or self.feh_calibration is not None:
                if self.feh_calibration is None:
                    raise UnvalidatedCalibrationError(
                        "calibrate_feh=True no longer selects embedded coefficients. "
                        "Supply a validated FeHCalibrationSet artifact."
                    )
                if self.order is None:
                    raise CalibrationError("An order is required to select [Fe/H] calibration.")
                calibration_set = coerce_calibration_set(self.feh_calibration)
                selected = calibration_set.select(
                    self.order,
                    teff=self.calibration_population_teff,
                    population=self.calibration_population,
                    library_id=self.library_id,
                )
                self.feh = selected.apply(self.feh_raw)
                self.feh_calibration_metadata = {
                    'status': 'applied',
                    'artifact_id': calibration_set.artifact_id,
                    'artifact_sha256': calibration_set.canonical_sha256,
                    'population_membership_teff': self.calibration_population_teff,
                    **selected.to_dict(),
                }
                _log(self.verbose, 'Calibrating [Fe/H]: {:0.3f} -> {:0.3f} '
                     '(order {}, population {!r})'.format(
                         self.feh_raw, self.feh, self.order,
                         selected.population.name))
        else:
            self.feh_raw = np.nan
            self.feh = np.nan
        if len(self.loggs) > 0:
            self.logg = weighted_value(self.loggs, weights)
        else:
            self.logg = np.nan
        if len(self.vsinis) > 0:
            self.vsini = weighted_value(self.vsinis, weights)
        else:
            self.vsini = np.nan
        _log(self.verbose, 'Stellar parameters: Teff={:0.1f} K, [Fe/H]={:0.3f} '
             'dex, logg={:0.3f} dex, vsini={:0.3f} km/s'.format(
                 self.teff, self.feh, self.logg, self.vsini))

    def plot_model(self, pv):
        pv_all = self.lpf.get_pv_all(pv)
        fig, (ax, bx) = plt.subplots(nrows=2, dpi=200, sharex=True, gridspec_kw={'height_ratios': [5, 2]})
        ax.plot(self.lpf.w, self.lpf.data_target['f'], color='black', label='Target', lw=1)
        ff = self.lpf.compute_model(pv)

        ax.plot(self.lpf.w, ff, color='crimson', label='Composite', alpha=0.5, lw=1)
        ax.legend(fontsize=10)
        bx.plot(self.lpf.w, self.lpf.data_target['f'] - ff, lw=1)
        bx.set_xlabel('Wavelength [A]', fontsize=12)
        bx.set_ylabel('Residual', fontsize=12)

        self.calculate_stellar_parameters(pv_all)

        title = 'log_ln={}\n'.format(self.lpf(pv))
        title += ' '.join(['c_{}={:0.5f}'.format(i, pv_all[i]) for i in range(5)])
        title += '\nTeff={:0.3f}, Fe/H={:0.3f}, logg={:0.3f}'.format(self.teff, self.feh, self.logg)
        if self.teff_known is not None and self.feh_known is not None and self.logg_known is not None:
            title += '\nTeff_0={:0.3f}+-{:0.1f}K, Fe/H_0={:0.3f}+-{:0.3f}, logg_0={:0.3f}+-{:0.3f}'.format(
                self.teff_known, self.tefferr_known, self.feh_known, self.feherr_known, self.logg_known,
                self.loggerr_known)
        ax.set_title(title, fontsize=10)
        for xx in (ax, bx):
            utils.ax_apply_settings(xx, ticksize=10)
        fig.subplots_adjust(hspace=0.05)

    def plot_model_with_components(self, pv, fig=None, ax=None, names=None, savename='compositeComparison.pdf',
                                   title='', scaleres=1., mode='HR'):
        """
        INPUT:
            scaleres - amount to scale residuals from composite spectrum (default 1)
        Make a plot of the 5 best stars and compare to targets spectrum.
        Compare composite spectrum to target spectrum and plot residuals.
        """
        w = self.lpf.w
        pv_all = self.lpf.get_pv_all(pv)
        self.calculate_stellar_parameters(pv_all)
        if fig is None and ax is None:
            fig, ax = plt.subplots(dpi=200)
        ax.plot(w, self.lpf.data_target['f'], color='black', label='Target', lw=1)
        ff = self.lpf.compute_model(pv)
        ax.plot(w, ff, color='crimson', label='Composite', alpha=0.5, lw=1.5)

        for i in range(5):
            ax.plot(w, self.lpf.data_refs['f'][i] + (5.0 - i), lw=1, color='black')
            if names is not None:
                label = '{}, $c_{}$={:0.5f}'.format(names[i].replace('_', ' '), i + 1, self.lpf.get_pv_all(pv)[i])
                label += ', Teff={:0.0f}K'.format(self.teffs[i])
                label += ', [Fe/H]={:0.2f}'.format(self.fehs[i])
            else:
                label = '$c_{}$={:0.5f}'.format(i + 1, pv_all[i])
                label += ', Teff={:0.0f}K'.format(self.teffs[i])
                label += ', [Fe/H]={:0.2f}'.format(self.fehs[i])
            ax.text(w[0], (6.12 - i), label, color='black', fontsize=8)

        ax.text(w[0], 1.15, 'Target Spectrum (Black), Composite Spectrum (Red)', fontsize=8)
        ax.text(w[0], 0.15, 'Residual: Target - Composite (Scale: {:0.0f}x)'.format(scaleres), fontsize=8)
        if mode == 'HR':
            if self.vsini < 2.8:
                title += 'Target={}, Teff={:0.3f}, Fe/H={:0.3f}, logg={:0.3f}, vsini<{:0.3f}km/s'.format(
                    self.targetname,self.teff, self.feh,self.logg, 2.8)
            else:
                title += 'Target={}, Teff={:0.3f}, Fe/H={:0.3f}, logg={:0.3f}, vsini={:0.3f}km/s'.format(
                    self.targetname,self.teff, self.feh,self.logg, self.vsini)
        elif mode == 'HE':
            if self.vsini < 4.5:
                title += 'Target={}, Teff={:0.3f}, Fe/H={:0.3f}, logg={:0.3f}, vsini<{:0.3f}km/s'.format(
                    self.targetname,self.teff, self.feh,self.logg, 4.5)
            else:
                title += 'Target={}, Teff={:0.3f}, Fe/H={:0.3f}, logg={:0.3f}, vsini={:0.3f}km/s'.format(
                    self.targetname,self.teff, self.feh,self.logg, self.vsini)

        ax.set_title(title, fontsize=10)

        ax.plot(w, (self.lpf.data_target['f'] - ff) * scaleres, color='black', lw=1)
        ax.set_xlabel('Wavelength [A]', fontsize=12, labelpad=2)
        ax.set_ylabel('Flux (+offset)', fontsize=12, labelpad=2)
        ax.set_ylim(-1,7)
        # ax.set_xlim(8645, 8650)
        # ax.set_title('log_ln={}, c={}'.format(self.lpf(pv),str(self.lpf.get_pv_all(pv))),fontsize=10)
        utils.ax_apply_settings(ax, ticksize=10)
        fig.tight_layout()
        fig.savefig(savename, dpi=200)
        _log(self.verbose, 'Saved to {}'.format(savename))

    def plot_model_with_component_save(self, pv, fig=None, ax=None, names=None, savename='compositeComparison.pdf',
                                   title='', scaleres=1., mode='HR', datafile='plot_data.pkl'):
        """
        INPUT:
            scaleres - amount to scale residuals from composite spectrum (default 1)
        Make a plot of the 5 best stars and compare to targets spectrum.
        Compare composite spectrum to target spectrum and plot residuals.
        """
        w = self.lpf.w
        pv_all = self.lpf.get_pv_all(pv)
        self.calculate_stellar_parameters(pv_all)
        target_spectrum = self.lpf.data_target['f']
        composite_spectrum = self.lpf.compute_model(pv)
        ref_spectra = [self.lpf.data_refs['f'][i] for i in range(5)]
        labels = []
        for i in range(5):
            if names is not None:
                label = '{}, $c_{}$={:0.5f}'.format(names[i].replace('_', ' '), i + 1, pv_all[i])
                label += ', Teff={:0.0f}K'.format(self.teffs[i])
                label += ', [Fe/H]={:0.2f}'.format(self.fehs[i])
            else:
                label = '$c_{}$={:0.5f}'.format(i + 1, pv_all[i])
                label += ', Teff={:0.0f}K'.format(self.teffs[i])
                label += ', [Fe/H]={:0.2f}'.format(self.fehs[i])
            labels.append(label)

        data_path = Path(datafile).with_suffix('.npz')
        metadata_path = data_path.with_suffix('.json')
        np.savez_compressed(
            data_path,
            wavelength=np.asarray(w),
            target_spectrum=np.asarray(target_spectrum),
            composite_spectrum=np.asarray(composite_spectrum),
            reference_spectra=np.asarray(ref_spectra),
            weights=np.asarray(pv_all),
        )
        plot_metadata = {
            'schema_version': 1,
            'array_file': data_path.name,
            'labels': labels,
            'title': title,
            'scaleres': scaleres,
            'mode': mode,
            'vsini': self.vsini,
            'targetname': self.targetname,
            'teff': self.teff,
            'feh': self.feh,
            'logg': self.logg
        }
        with metadata_path.open('w', encoding='utf-8') as stream:
            json.dump(plot_metadata, stream, indent=2, sort_keys=True,
                      default=_json_default)
            stream.write('\n')
        _log(self.verbose, f'Saved plot data to {data_path} and {metadata_path}')

        # ... (rest of your plotting code here)

    def minimize_PyDE(self, npop=100, de_iter=200, mc_iter=1000, mcmc=False,
                      threads=1, maximize=True, plot_priors=True,
                      sample_ball=False, random_seed=0, tolerance=1e-7,
                      verbose=None):
        """Solve the five-weight simplex problem exactly by active-set search.

        The historical method name is retained for API compatibility.  The
        composite spectrum is linear in its five non-negative weights, so all
        31 nonempty active faces are solved by equality-constrained weighted
        least squares; stochastic differential evolution is neither needed
        nor used.  The fixed loss uses only the target's diagonal errors.
        Reference errors and resampling/broadening covariance are not included,
        so this is an optimization loss rather than a chi-square distribution.
        """
        verbose = self.verbose if verbose is None else verbose
        _log(verbose, "Solving exact non-negative composite weights")
        target = np.asarray(self.lpf.data_target['f'], dtype=float)
        error = np.asarray(self.lpf.data_target['e'], dtype=float)
        references = np.asarray(self.lpf.data_refs['f'], dtype=float)
        if references.shape != (self.lpf.num_refs, target.size):
            raise ValueError("Composite reference array has an invalid shape.")
        valid = np.asarray(self.lpf.valid_mask, dtype=bool)
        design = references[:, valid].T / error[valid, None]
        observed = target[valid] / error[valid]
        best = None
        feasible_candidates = []
        evaluated = 0
        feasibility_tolerance = 1e-12
        for size in range(1, self.lpf.num_refs + 1):
            for active in itertools.combinations(range(self.lpf.num_refs), size):
                evaluated += 1
                columns = design[:, active]
                particular = np.full(size, 1.0 / size)
                if size == 1:
                    active_weights = particular
                    equality_rank = 0
                    equality_singular = np.array([1.0])
                else:
                    # scipy.linalg.helmert is a deterministic orthonormal
                    # basis for the sum-zero subspace.
                    null_basis = scipy.linalg.helmert(
                        size, full=False
                    ).T
                    reduced_design = columns @ null_basis
                    reduced_target = observed - columns @ particular
                    coordinates, _, equality_rank, equality_singular = np.linalg.lstsq(
                        reduced_design, reduced_target, rcond=None
                    )
                    active_weights = particular + null_basis @ coordinates
                if (
                    not np.all(np.isfinite(active_weights))
                    or np.min(active_weights) < -feasibility_tolerance
                ):
                    continue
                # Only remove roundoff-level negatives; a materially negative
                # equality solution belongs to a lower-dimensional face.
                active_weights[np.abs(active_weights) < feasibility_tolerance] = 0.0
                total = float(active_weights.sum())
                if not np.isfinite(total) or total <= 0 or abs(total - 1.0) > 1e-8:
                    continue
                active_weights /= total
                weights = np.zeros(self.lpf.num_refs, dtype=float)
                weights[list(active)] = active_weights
                closure = 1.0 - float(np.sum(weights))
                weights[int(np.argmax(weights))] += closure
                residual = observed - design @ weights
                weighted_chi2 = float(residual @ residual)
                candidate = (
                    weighted_chi2, weights, active, int(equality_rank),
                    equality_singular,
                )
                feasible_candidates.append(candidate)
                if best is None or candidate[0] < best[0]:
                    best = candidate
        if best is None:
            raise RuntimeError("No feasible composite-weight active set was found.")
        optimum_tolerance = max(1e-12, 1e-10 * max(1.0, abs(best[0])))
        tied_candidates = [
            item for item in feasible_candidates
            if item[0] <= best[0] + optimum_tolerance
        ]
        best = min(
            tied_candidates,
            key=lambda item: (len(item[2]), tuple(item[2]), tuple(item[1])),
        )
        weighted_chi2, weights, active, equality_rank, equality_singular = best
        # Nested active sets can encode the exact same boundary solution.  Only
        # materially distinct weight vectors establish non-identifiability.
        distinct_tied_weights = []
        for candidate in tied_candidates:
            candidate_weights = candidate[1]
            if not any(
                np.max(np.abs(candidate_weights - prior)) <= 1e-9
                for prior in distinct_tied_weights
            ):
                distinct_tied_weights.append(candidate_weights.copy())
        effective_active = tuple(np.flatnonzero(
            weights > SIMPLEX_ACTIVE_WEIGHT_TOL
        ).tolist())
        if not effective_active:
            raise RuntimeError("Composite solution has no positive component.")
        self.optimal_weights = weights.copy()
        self.min_pv = weights[:-1].copy()
        model_valid = design @ weights
        self.min_pv_lnval = float(ll_normal_ev_py(
            target[valid], model_valid * error[valid], error[valid]
        ))
        gradient = 2.0 * design.T @ (design @ weights - observed)
        lagrange = -float(np.mean(gradient[list(effective_active)]))
        stationarity = float(np.max(np.abs(
            gradient[list(effective_active)] + lagrange
        )))
        inactive = sorted(
            set(range(self.lpf.num_refs)).difference(effective_active)
        )
        dual_min = (
            float(np.min(gradient[inactive] + lagrange))
            if inactive else None
        )
        positive_singular = np.asarray(equality_singular)[
            np.asarray(equality_singular) > np.finfo(float).eps
        ]
        condition = (
            float(np.max(positive_singular) / np.min(positive_singular))
            if positive_singular.size else None
        )
        success = bool(
            np.all(np.isfinite(weights))
            and np.isfinite(self.min_pv_lnval)
            and abs(weights.sum() - 1.0) <= 1e-10
            and np.min(weights) >= 0.0
            and stationarity <= 1e-6 * max(1.0, np.max(np.abs(gradient)))
            and (dual_min is None or dual_min >= -1e-6 * max(
                1.0, np.max(np.abs(gradient))
            ))
        )
        if len(distinct_tied_weights) > 1:
            tied_matrix = np.asarray(distinct_tied_weights)
            maximum_weight_spread = float(np.max(
                np.ptp(tied_matrix, axis=0)
            ))
            maximum_l1_weight_distance = float(max(
                np.sum(np.abs(left - right))
                for left in tied_matrix for right in tied_matrix
            ))
        else:
            maximum_weight_spread = 0.0
            maximum_l1_weight_distance = 0.0
        parameter_spreads = {}
        for name, values in (
            ('teff_k', self.teffs), ('feh_dex', self.fehs),
            ('logg_dex', self.loggs), ('vsini_kms', self.vsinis),
        ):
            array = np.asarray(values, dtype=float)
            if array.shape == (self.lpf.num_refs,) and np.all(np.isfinite(array)):
                recovered = np.asarray([
                    float(item @ array) for item in distinct_tied_weights
                ])
                parameter_spreads[name] = float(np.ptp(recovered))
        nonunique = len(distinct_tied_weights) > 1
        self.optimizer_diagnostics = {
            'name': 'enumerated_active_set_weighted_least_squares',
            'success': success,
            'message': (
                'exact simplex optimum found for fixed diagonal target-error loss'
                if success else 'KKT check failed'
            ),
            'nfev': evaluated,
            'nit': 0,
            'active_reference_indices': list(effective_active),
            'active_set_size': len(effective_active),
            'subsets_evaluated': evaluated,
            'feasible_subsets': len(feasible_candidates),
            'objective_tied_subsets': len(tied_candidates),
            'distinct_tied_weight_vectors': len(distinct_tied_weights),
            'nonunique_optimum': nonunique,
            'maximum_tied_weight_spread': maximum_weight_spread,
            'maximum_tied_l1_weight_distance': maximum_l1_weight_distance,
            'tied_stellar_parameter_spreads': parameter_spreads,
            'equality_null_space_rank': equality_rank,
            'equality_null_space_condition': condition,
            'constraint_residual': float(abs(weights.sum() - 1.0)),
            'minimum_weight': float(np.min(weights)),
            'stationarity_residual': stationarity,
            'inactive_dual_minimum': dual_min,
            'weighted_residual_sum_of_squares': weighted_chi2,
            'weighted_chi2_legacy_name': weighted_chi2,
            'loss_interpretation': (
                'fixed_diagonal_target_error_weighted_loss_not_a_chi_square_distribution'
            ),
            'reference_errors_used': False,
            'resampling_covariance_used': False,
            'valid_pixel_count': int(np.count_nonzero(valid)),
            'objective_log_likelihood': float(self.min_pv_lnval),
        }
        if not success:
            raise RuntimeError(
                "Composite active-set weighted least-squares KKT check failed."
            )
        _log(verbose, "Optimization complete; log likelihood={:0.6g}".format(
            self.min_pv_lnval))
        if _verbosity_level(verbose) >= 2:
            self.print_param_diagnostics(self.min_pv)
        if mcmc:
            raise NotImplementedError(
                "Composite-weight MCMC is disabled in 0.2: boundary-safe "
                "simplex sampling and the correlated spectral likelihood have "
                "not been validated. Use the deterministic optimum and fresh "
                "empirical cross-validation uncertainties."
            )

    def print_param_diagnostics(self, pv):
        """
        A function to print nice parameter diagnostics.
        """
        self.df_diagnostics = pd.DataFrame(
            zip(self.lpf.ps.labels, self.lpf.ps.centers, self.lpf.ps.bounds[:, 0], self.lpf.ps.bounds[:, 1], pv,
                self.lpf.ps.centers - pv), columns=["labels", "centers", "lower", "upper", "pv", "center_dist"])
        print(self.df_diagnostics.to_string())
        return self.df_diagnostics

    # def plot_chains(self,labels=None,burn=0,thin=1):
    #    print("Plotting chains")
    #    if labels==None:
    #        labels = self.lpf.ps.descriptions
    #    mcFunc.plot_chains(self.sampler.chain,labels=labels,burn=burn,thin=thin)

    # def plot_corner(self,labels=None,burn=0,thin=1,title_fmt='.5f',**kwargs):
    #    if labels==None:
    #        labels = self.lpf.ps.descriptions
    #    self.fig = mcFunc.plot_corner(self.sampler.chain,labels=labels,burn=burn,thin=thin,title_fmt=title_fmt,**kwargs)


# def sample(df_chain, N=500):
#     spectra = []
#     for i in range(N):
#         c = df_chain.iloc[np.random.randint(len(df_chain))].values[0:4]
#         spectra.append(self.lpf.data_target['f'] - LCS.lpf.compute_model(c))
#     std = np.std(spectra, axis=0)  # np.std(spectra,axis=0)
#     return std


class PairwiseRSSFunctionVsiniPolynomial(object):
    """Unweighted pairwise residual-sum-of-squares objective.

    The historical implementation supplied unit errors to a chi-square helper,
    making the numerical objective RSS rather than an uncertainty-weighted
    chi-square.  Version 0.2 preserves that validated behavior and names it
    explicitly; changing the weights requires fresh cross-validation.
    """
    def __init__(self, w, f1, e1, f2, e2, maxvsini, vsini=None,
                 vsini_window=None, comparison_mask=None):
        self.w = np.asarray(w, dtype=float)
        self.chebyshev_coordinate = _scaled_chebyshev_coordinate(self.w)
        self.wavelength_midpoint = 0.5 * (self.w[0] + self.w[-1])
        self.wavelength_half_range = 0.5 * (self.w[-1] - self.w[0])
        self.data_target = {'f': f1,
                            'e': e1}
        self.data_ref = {'f': f2,
                         'e': e2}
        if vsini is None:
            self.fixed_vsini = None
            vsini_min, vsini_max = 0.0, float(maxvsini)
        elif vsini_window is None or float(vsini_window) == 0.0:
            if not 0.0 <= float(vsini) <= float(maxvsini):
                raise ValueError("Fixed vsini must lie between zero and maxvsini.")
            self.fixed_vsini = float(vsini)
            # PriorSet still describes the full parameter vector.  Exact fixing
            # is implemented by optimizing only the polynomial coefficients.
            delta = max(np.spacing(abs(self.fixed_vsini)), np.finfo(float).eps) * 2.0
            vsini_min = self.fixed_vsini - delta
            vsini_max = self.fixed_vsini + delta
        else:
            if float(vsini_window) < 0:
                raise ValueError("vsini_window cannot be negative.")
            self.fixed_vsini = None
            vsini_min = max(0.0, float(vsini) - float(vsini_window))
            vsini_max = min(float(maxvsini), float(vsini) + float(vsini_window))
            if not vsini_min < vsini_max:
                raise ValueError("The constrained vsini interval is empty.")
        self.vsini_bounds = (vsini_min, vsini_max)
        support_vsini = (
            self.fixed_vsini if self.fixed_vsini is not None else vsini_max
        )
        upper_support = rotbroad_help.broaden(
            self.w, np.asarray(f2, dtype=float), support_vsini,
            u1=LIMB_DARKENING_U1,
        )
        upper_variance = rotbroad_help.broaden_variance(
            self.w, np.asarray(e2, dtype=float) ** 2,
            support_vsini, u1=LIMB_DARKENING_U1,
        )
        support_mask = (
            np.isfinite(np.asarray(f1, dtype=float))
            & np.isfinite(np.asarray(e1, dtype=float))
            & (np.asarray(e1, dtype=float) > 0.0)
            & np.isfinite(upper_support)
            & np.isfinite(upper_variance)
            & (upper_variance > 0.0)
        )
        if comparison_mask is not None:
            comparison_mask = np.asarray(comparison_mask, dtype=bool)
            if comparison_mask.shape != self.w.shape:
                raise ValueError("Pairwise comparison_mask has an invalid shape.")
            if np.any(comparison_mask & ~support_mask):
                raise ValueError(
                    "Pairwise comparison_mask includes unsupported target/reference pixels."
                )
            support_mask = comparison_mask
        self.comparison_mask = support_mask
        self.comparison_pixel_count = int(np.count_nonzero(support_mask))
        self.comparison_pixel_fraction = float(np.mean(support_mask))
        if self.comparison_pixel_count <= 6:
            raise ValueError("Too few frozen common pixels for pairwise fitting.")
        self.priors = [UP(vsini_min, vsini_max, 'vsini', r'$v \sin i$', priortype="model"),
                       UP(-1e10, 1e10, 'c0', 'c_0', priortype="model"),
                       UP(-1e10, 1e10, 'c1', 'c_1', priortype="model"),
                       UP(-1e10, 1e10, 'c2', 'c_2', priortype="model"),
                       UP(-1e10, 1e10, 'c3', 'c_3', priortype="model"),
                       UP(-1e10, 1e10, 'c4', 'c_4', priortype="model"),
                       UP(-1e10, 1e10, 'c5', 'c_5', priortype="model")]
        self.ps = PriorSet(self.priors)

    def compute_model(self, pv, eps=LIMB_DARKENING_U1):
        """
        Multiply reference by polynomial and then rotationally broaden.
        """
        vsini = pv[0]
        coeffs = pv[1:]
        ff_ref, _ = _prepare_reference_model(
            self.w, self.data_ref['f'], self.data_ref['e'], coeffs, vsini,
            limb_darkening=eps,
        )
        return ff_ref

    def solve_coefficients(self, vsini, eps=LIMB_DARKENING_U1):
        """Return the exact linear Chebyshev solution at one ``vsini``."""
        basis = np.polynomial.chebyshev.chebvander(
            self.chebyshev_coordinate, CONTINUUM_DEGREE
        )
        design = np.column_stack([
            rotbroad_help.broaden(
                self.w, np.asarray(self.data_ref['f'], dtype=float) * basis[:, k],
                float(vsini), u1=float(eps),
            )
            for k in range(basis.shape[1])
        ])
        target = np.asarray(self.data_target['f'], dtype=float)
        valid = self.comparison_mask
        if (
            not np.all(np.isfinite(target[valid]))
            or not np.all(np.isfinite(design[valid]))
        ):
            raise ValueError("Frozen pairwise mask is invalid at this vsini.")
        coefficients, _, rank, singular = np.linalg.lstsq(
            design[valid], target[valid], rcond=None
        )
        residual = target[valid] - design[valid] @ coefficients
        rss = float(residual @ residual)
        positive = singular[singular > np.finfo(float).eps]
        condition = (
            float(np.max(positive) / np.min(positive))
            if positive.size else None
        )
        if rank < basis.shape[1]:
            raise RuntimeError("Pairwise Chebyshev design is rank deficient.")
        if condition is None or not np.isfinite(condition) or condition > 1e12:
            raise RuntimeError("Pairwise Chebyshev design is ill-conditioned.")
        return coefficients, rss, {
            'rank': int(rank),
            'condition': condition,
            'valid_pixel_count': int(np.count_nonzero(valid)),
        }

    def __call__(self, pv, verbose=False):
        if any(pv < self.ps.pmins) or any(pv > self.ps.pmaxs):
            return np.inf
        flux_model = self.compute_model(pv)
        flux_target = self.data_target['f']
        valid = self.comparison_mask
        if not np.all(np.isfinite(np.asarray(flux_model)[valid])):
            return np.inf
        residual = np.asarray(flux_target)[valid] - np.asarray(flux_model)[valid]
        return float(np.sum(residual ** 2))


# Compatibility name for downstream imports. The objective is RSS, not a
# statistically error-weighted chi-square.
Chi2FunctionVsiniPolynomial = PairwiseRSSFunctionVsiniPolynomial


class FitTargetRefStarVsiniPolynomial(object):
    """
    A class to fit 5 lin-comb spectra together. Note: look at LPFunctionLinComb
    """

    def __init__(self, Chi2FunctionVsiniPolynomial):
        self.chi2f = Chi2FunctionVsiniPolynomial

    def plot_model(self, pv):
        vsini = pv[0]
        coeffs = pv[1:]
        fig, (ax, bx) = plt.subplots(nrows=2, dpi=200, sharex=True, gridspec_kw={'height_ratios': [5, 2]})
        ax.plot(self.chi2f.w, self.chi2f.data_target['f'], color='black', label='Target', lw=1)
        ff = self.chi2f.compute_model(pv)
        ax.plot(
            self.chi2f.w,
            np.polynomial.chebyshev.chebval(
                self.chi2f.chebyshev_coordinate, coeffs
            ),
        )

        ax.plot(self.chi2f.w, ff, color='crimson', label='Reference (vsini={:0.3f}km/s)'.format(vsini), alpha=0.5, lw=1)
        ax.legend(fontsize=10)
        bx.plot(self.chi2f.w, self.chi2f.data_target['f'] - ff, lw=1)
        bx.set_xlabel('Wavelength [A]', fontsize=12)
        bx.set_ylabel('Residual', fontsize=12)

        title = r'$\chi^2$={}, coeffs={}'.format(self.chi2f(pv), coeffs)
        ax.set_title(title, fontsize=10)
        for xx in (ax, bx):
            utils.ax_apply_settings(xx, ticksize=10)
        fig.subplots_adjust(hspace=0.05)

    def minimize_AMOEBA(self, centers=None, verbose=0):
        """Profile exact linear coefficients over a deterministic 1-D vsini fit."""
        lower, upper = self.chi2f.vsini_bounds
        evaluations = 0

        def solve(value):
            nonlocal evaluations
            evaluations += 1
            coefficients, rss, diagnostics = self.chi2f.solve_coefficients(value)
            return float(rss), np.asarray(coefficients), diagnostics

        if self.chi2f.fixed_vsini is not None:
            best_vsini = self.chi2f.fixed_vsini
            best_rss, best_coefficients, linear_diagnostics = solve(best_vsini)
            refinement_success = True
            coarse_size = 1
        else:
            width = upper - lower
            coarse_size = max(3, int(np.ceil(width / 1.0)) + 1)
            coarse_grid = np.linspace(lower, upper, coarse_size)
            coarse = [solve(value) for value in coarse_grid]
            coarse_rss = np.array([item[0] for item in coarse])
            candidates = [
                (float(value), *result)
                for value, result in zip(coarse_grid, coarse)
            ]
            refinement_success = True
            local_minima = [
                index for index in range(1, coarse_size - 1)
                if (
                    coarse_rss[index] <= coarse_rss[index - 1]
                    and coarse_rss[index] <= coarse_rss[index + 1]
                )
            ]
            for best_index in local_minima:
                bracket = (float(coarse_grid[best_index - 1]),
                           float(coarse_grid[best_index + 1]))

                def profile_rss(value):
                    return solve(value)[0]

                with _quiet_external_calls(verbose):
                    refinement = scipy.optimize.minimize_scalar(
                        profile_rss, bounds=bracket, method='bounded',
                        options={'xatol': 1e-4, 'maxiter': 100},
                    )
                refinement_success &= bool(refinement.success)
                refined_rss, refined_coefficients, refined_diagnostics = solve(
                    float(refinement.x)
                )
                candidates.append((
                    float(refinement.x), refined_rss,
                    refined_coefficients, refined_diagnostics,
                ))
            minimum_rss = min(item[1] for item in candidates)
            tie_tolerance = max(1e-12, 1e-10 * max(1.0, abs(minimum_rss)))
            tied = [
                item for item in candidates
                if item[1] <= minimum_rss + tie_tolerance
            ]
            best_vsini, best_rss, best_coefficients, linear_diagnostics = min(
                tied, key=lambda item: item[0]
            )
            vsini_unresolved = (
                max(item[0] for item in tied) - min(item[0] for item in tied)
                > 1e-3
            )

        self.min_pv = np.r_[best_vsini, best_coefficients]
        objective_value = float(self.chi2f(self.min_pv))
        continuum = np.polynomial.chebyshev.chebval(
            self.chi2f.chebyshev_coordinate, best_coefficients
        )
        continuum_minimum = float(np.min(
            continuum[self.chi2f.comparison_mask]
        ))
        if self.chi2f.fixed_vsini is not None:
            vsini_unresolved = False
        success = bool(
            refinement_success
            and np.all(np.isfinite(self.min_pv))
            and np.isfinite(objective_value)
            and np.isclose(objective_value, best_rss, rtol=1e-9, atol=1e-12)
            and continuum_minimum > 0.0
        )
        self.optimizer_diagnostics = {
            'name': 'profiled_linear_chebyshev_plus_deterministic_1d_vsini_search',
            'success': success,
            'message': (
                'deterministic coarse profile and all sampled-local-minimum '
                'basin refinements complete; not a proof of the continuous '
                'global vsini optimum'
            ),
            'nfev': evaluations,
            'coarse_grid_size': coarse_size,
            'vsini_bounds_kms': [float(lower), float(upper)],
            'vsini_kms': float(best_vsini),
            'linear_rank': linear_diagnostics['rank'],
            'linear_condition': linear_diagnostics['condition'],
            'valid_pixel_count': linear_diagnostics['valid_pixel_count'],
            'valid_pixel_fraction': self.chi2f.comparison_pixel_fraction,
            'continuum_minimum': continuum_minimum,
            'continuum_positive': continuum_minimum > 0.0,
            'vsini_fixed_by_user': self.chi2f.fixed_vsini is not None,
            'vsini_boundary_hit': bool(
                self.chi2f.fixed_vsini is None
                and (
                    np.isclose(best_vsini, lower, atol=1e-4)
                    or np.isclose(best_vsini, upper, atol=1e-4)
                )
            ),
            'vsini_unresolved_on_profile': bool(vsini_unresolved),
            'objective_rss': objective_value,
        }
        if not success:
            raise RuntimeError("Deterministic pairwise profiled optimization failed.")

    def minimize_PyDE(self, npop=100, de_iter=200, mc_iter=1000, mcmc=False,
                      threads=1, maximize=False, plot_priors=True,
                      sample_ball=False, random_seed=0, tolerance=1e-7,
                      verbose=0):
        """Compatibility alias for the deterministic profiled solver."""
        if mcmc:
            raise NotImplementedError("MCMC is not implemented for the chi-square stage.")
        return self.minimize_AMOEBA(verbose=verbose)


def chi2spectraPolyVsini(ww, H1, H2, rv1=None, rv2=None, plot=False,
                         verbose=0, maxvsini=30., order=101,
                         deblazed=False, vsini=None, vsini_window=None,
                         return_diagnostics=False, prepared=None,
                         comparison_mask=None):
    """
    INPUT:
        ww - wavelength grid to interpolate on (array)
        H1 - target spectrum (HPFSpectrum object)
        H2 - reference spectrum (HPFSpectrum object)
        rv1 - radial velocity H1 km/s (float)
        rv2 - radial velocity H2 km/s (float)
        plot - (boolean)
        verbose - print additional info (boolean)

    OUTPUT:
        chi2 - chi2 values for the comparison
        vsini - 
        coeffs - 
        
    EXAMPLE:
        H1 = HPFSpectrum(df[df.name=='G_9-40'].filename.values[0])
        H2 = HPFSpectrum(df[df.name=='AD_Leo'].filename.values[0])

        wmin = 10280.
        wmax = 10380.
        ww = np.arange(wmin,wmax,0.01)
        chi2spectraPolyVsini(ww,H1,H2,rv1=14.51,plot=True)
        
    """
    if prepared is None:
        with _quiet_external_calls(verbose):
            ff1, ee1 = H1.resample_order(ww, order=order, deblazed=deblazed)
            ff2, ee2 = H2.resample_order(ww, order=order)
    else:
        ff1, ee1, ff2, ee2 = prepared
    if plot:
        plt.figure()
        plt.plot(ww, ff1, label='target')
        plt.plot(ww, ff2, label='reference')
    C = PairwiseRSSFunctionVsiniPolynomial(
        ww, ff1, ee1, ff2, ee2, maxvsini,
        vsini=vsini, vsini_window=vsini_window,
        comparison_mask=comparison_mask,
    )
    FTRSVP = FitTargetRefStarVsiniPolynomial(C)
    FTRSVP.minimize_AMOEBA(verbose=verbose)
    fitted_vsini = FTRSVP.min_pv[0]
    coeffs = FTRSVP.min_pv[1:]
    rss = C(FTRSVP.min_pv)

    if plot:
        FTRSVP.plot_model(FTRSVP.min_pv)

    result = (rss, fitted_vsini, coeffs)
    if return_diagnostics:
        return (*result, FTRSVP.optimizer_diagnostics)
    return result


def chi2spectraPolyLoop(ww, H1, Hrefs, plot_all=False, plot_chi=False,
                        verbose=1, maxvsini=30., order=101,
                        deblazed=False, vsini=None, vsini_window=None,
                        return_comparison_mask=False):
    """
    Calculate chi square - target and list of reference spectra
    
    INPUT:
        ww - wavelength grid to interpolate on (array)
        H1 - target spectrum (HPFSpectrum object)
        Hrefs - reference spectra (HPFSpectraList object)
        plot_all - creates many additional plots (boolean)
        plot_chi = plot chi2 H1 vs other stars (boolean)
        verbose - print additional info (boolean)
    
    OUTPUT:
        df - dataframe of all reference stars sorted by chi2 (chi2, poly_params, and vsini values)
        df_best - dataframe with 5 best reference stars (chi2, poly_params, and vsini values)
        Hrefs_best - 5 best fitting reference stars
        
    EXAMPLE:
        
    """
    chis = []
    poly_params = []
    vsinis = []
    optimizer_diagnostics = []
    with _quiet_external_calls(verbose):
        target_flux, target_error = H1.resample_order(
            ww, order=order, deblazed=deblazed
        )
        prepared_references = [
            H2.resample_order(ww, order=order) for H2 in Hrefs
        ]
    support_masks = []
    for (reference_flux, reference_error) in prepared_references:
        provisional = PairwiseRSSFunctionVsiniPolynomial(
            ww, target_flux, target_error, reference_flux, reference_error,
            maxvsini, vsini=vsini, vsini_window=vsini_window,
        )
        support_masks.append(provisional.comparison_mask)
    common_mask = np.logical_and.reduce(support_masks)
    common_count = int(np.count_nonzero(common_mask))
    if (
        common_count < MIN_PAIRWISE_PIXELS
        or common_count / len(common_mask) < MIN_PAIRWISE_PIXEL_FRACTION
    ):
        raise ValueError(
            "Too few pixels remain in the frozen target/all-reference mask: "
            f"{common_count}/{len(common_mask)}; require at least "
            f"{MIN_PAIRWISE_PIXELS} and fraction "
            f"{MIN_PAIRWISE_PIXEL_FRACTION:.2f}."
        )
    for i, (H2, (reference_flux, reference_error)) in enumerate(
            zip(Hrefs, prepared_references)):
        if i == 0:
            _log(verbose, 'Matching target star to all library stars')

        # ``vsini`` is the immutable user constraint.  Never replace it with
        # one reference star's fitted value before fitting the next reference.
        rss, fitted_vsini, p, diagnostics = chi2spectraPolyVsini(
            ww, H1, H2, plot=plot_all, verbose=max(0, _verbosity_level(verbose) - 2),
            maxvsini=maxvsini, order=order, deblazed=deblazed,
            vsini=vsini, vsini_window=vsini_window,
            return_diagnostics=True,
            prepared=(
                target_flux, target_error, reference_flux, reference_error
            ),
            comparison_mask=common_mask,
        )
        if diagnostics.get('success') is not True or not np.isfinite(rss):
            raise RuntimeError(
                f"Pairwise optimization failed for reference {H2.object!r}."
            )
        if _verbosity_level(verbose) >= 2:
            print(
                '{:3d}/{:2d}, Target = {:18s} Library Star = {:18s} RSS = {:6.3f}'.format(i + 1, len(Hrefs), H1.object,
                                                                                           H2.object, rss))
        chis.append(rss)
        poly_params.append(p)
        vsinis.append(fitted_vsini)
        optimizer_diagnostics.append(diagnostics)

    refnames = [H.object for H in Hrefs]
    df = pd.DataFrame(
        list(zip(refnames, chis, poly_params, vsinis, optimizer_diagnostics)),
        columns=[
            'OBJECT_ID', 'rss', 'poly_params', 'vsini',
            'optimizer_diagnostics',
        ],
    )
    # Explicitly labeled legacy alias retained for 0.1 result readers.
    df['chi2_legacy_unit_error'] = df['rss']
    df = df.sort_values(
        ['rss', 'vsini', 'OBJECT_ID'], kind='mergesort'
    )
    df = df.reset_index(drop=False)
    df_best = df.iloc[:TOP_K].copy()
    Hrefs_best = neidspec.NEIDSpecList(np.array(Hrefs)[df_best['index'].values]);  # SEJ

    if plot_chi:
        fig, ax = plt.subplots(dpi=200)
        ax.plot(df.rss, lw=0, marker='h')
        ax.set_xticks(range(len(chis)))
        ax.set_ylabel('Unweighted RSS')
        utils.ax_apply_settings(ax)
        ax.set_xticklabels(df['OBJECT_ID'])
        ax.tick_params('x', labelsize=8)
        ax.set_title('{} vs other stars'.format(H1.object))
        for lab in ax.get_xticklabels():
            lab.set_rotation(90.)
        ax.set_yscale('log')
        # x.set_ylim(0.01,1e3)
    result = (df, df_best, Hrefs_best)
    if return_comparison_mask:
        return (*result, common_mask.copy())
    return result


def weighted_value(values, weights):
    """
    Array of weighted values
    """
    return np.dot(values, weights)


def _build_support_diagnostics(
        df_library, selected_references, weights, *, optimizer_nonunique=False):
    """Describe empirical-library support without changing the fit.

    Full-library boundaries are detected structurally from the active
    reference labels, rather than by comparing a rounded reported result.
    A simplex vertex is recorded but is not by itself considered outside the
    validated library domain; an exact full-library edge is.
    """
    weights = np.asarray(weights, dtype=float)
    selected = selected_references.reset_index(drop=True).copy()
    if len(selected) != TOP_K or weights.shape != (TOP_K,):
        raise ValueError(
            f"Support diagnostics require {TOP_K} selected weights."
        )
    if (
        not np.all(np.isfinite(weights))
        or np.min(weights) < -1e-12
        or not np.isclose(np.sum(weights), 1.0, rtol=0.0, atol=1e-10)
    ):
        raise ValueError("Support weights must be a finite simplex vector.")
    active = weights > SIMPLEX_ACTIVE_WEIGHT_TOL
    if not np.any(active):
        raise ValueError("Support diagnostics require an active reference.")

    parameter_columns = {
        'teff': ('Teff', 'K'),
        'feh_raw': ('[Fe/H]', 'dex'),
        'logg': ('log(g)', 'dex'),
    }
    parameters = {}
    edge_reasons = []
    for output_name, (column, unit) in parameter_columns.items():
        if column not in df_library or column not in selected:
            raise ValueError(f"Library support metadata lack {column!r}.")
        full_values = pd.to_numeric(
            df_library[column], errors='coerce'
        ).to_numpy(dtype=float)
        selected_values = pd.to_numeric(
            selected[column], errors='coerce'
        ).to_numpy(dtype=float)
        if (
            not np.all(np.isfinite(full_values))
            or not np.all(np.isfinite(selected_values))
        ):
            raise ValueError(
                f"Library support metadata contain non-finite {column!r}."
            )
        full_min = float(np.min(full_values))
        full_max = float(np.max(full_values))
        selected_min = float(np.min(selected_values))
        selected_max = float(np.max(selected_values))
        active_values = selected_values[active]
        scale = max(1.0, abs(full_min), abs(full_max))
        tolerance = 32.0 * np.finfo(float).eps * scale
        at_full_lower = bool(np.all(np.abs(active_values - full_min) <= tolerance))
        at_full_upper = bool(np.all(np.abs(active_values - full_max) <= tolerance))
        at_selected_lower = bool(
            np.all(np.abs(active_values - selected_min) <= tolerance)
        )
        at_selected_upper = bool(
            np.all(np.abs(active_values - selected_max) <= tolerance)
        )
        raw_value = float(np.dot(weights, selected_values))
        support_limited = bool(at_full_lower or at_full_upper)
        if at_full_lower:
            edge_reasons.append(
                f'{output_name}_at_full_reference_library_lower_boundary'
            )
        if at_full_upper:
            edge_reasons.append(
                f'{output_name}_at_full_reference_library_upper_boundary'
            )
        parameters[output_name] = {
            'raw_value': raw_value,
            'unit': unit,
            'full_library_min': full_min,
            'full_library_max': full_max,
            'selected_min': selected_min,
            'selected_max': selected_max,
            'at_full_library_lower_edge': at_full_lower,
            'at_full_library_upper_edge': at_full_upper,
            'at_selected_lower_edge': at_selected_lower,
            'at_selected_upper_edge': at_selected_upper,
            'support_limited': support_limited,
        }

    active_count = int(np.count_nonzero(active))
    effective_count = float(1.0 / np.sum(weights ** 2))
    is_vertex = active_count == 1
    any_full_edge = any(
        item['support_limited'] for item in parameters.values()
    )
    components = []
    for rank, row in selected.iterrows():
        record = {
            'pairwise_rank': int(rank + 1),
            'object_id': str(row['OBJECT_ID']),
            'weight': float(weights[rank]),
            'active': bool(active[rank]),
            'teff_k': float(row['Teff']),
            'feh_raw_dex': float(row['[Fe/H]']),
            'logg_dex': float(row['log(g)']),
        }
        for input_name, output_name in (
            ('vsini', 'pairwise_vsini_kms'),
            ('rss', 'pairwise_rss'),
        ):
            value = pd.to_numeric(
                pd.Series([row.get(input_name)]), errors='coerce'
            ).iloc[0]
            record[output_name] = float(value) if np.isfinite(value) else None
        components.append(record)
    warnings = []
    if is_vertex:
        warnings.append('composite_solution_is_single_reference_vertex')
    warnings.extend(edge_reasons)
    if is_vertex and any_full_edge:
        warnings.append('simplex_vertex_at_full_reference_library_boundary')
    return {
        'simplex': {
            'selected_count': TOP_K,
            'active_weight_tolerance': SIMPLEX_ACTIVE_WEIGHT_TOL,
            'active_count': active_count,
            'effective_component_count': effective_count,
            'maximum_weight': float(np.max(weights)),
            'is_vertex': is_vertex,
            'is_boundary': active_count < TOP_K,
            'optimizer_nonunique': bool(optimizer_nonunique),
        },
        'parameters': parameters,
        'selected_components': components,
        'any_full_library_edge': bool(any_full_edge),
        'vertex_at_full_library_edge': bool(is_vertex and any_full_edge),
        'publication_support_validated': not any_full_edge,
        'warnings': warnings,
    }


def run_specmatch(Htarget, Hrefs, ww, v, df_library, df_target=None,
                  plot=False, savefolder='out/', maxvsini=30.,
                  calibrate_feh=False, feh_calibration=None,
                  calibration_population=None,
                  calibration_population_teff=None, library_id=None,
                  scaleres=1., order=101, deblazed=False, mode='HR',
                  save_plot_data=False, save_legacy_pickle=False,
                  vsini=None, vsini_window=None, absrv=None,
                  reference_rvs=None, verbose=1, random_seed=0,
                  allow_mixed_drp=False, drp_compatibility=None,
                  library_manifest=None, dataset_id=None,
                  validation_summary=None, validation_population=None,
                  validation_population_teff=None,
                  allow_unmanifested_library=False, library_path=None):
    """
    Second chi2 loop, creates composite spectrum 
    
    INPUT:
        Htarget - target spectrum (HPFSpectrum object)
        Hrefs - reference spectra (HPFSpectraList object)
        ww - wavelength grid to interpolate on (array)
        v - velocities in km/s to use for absolute RV consideration (array)
        df_library - dataframe with info on Teff/FeH/logg for the library stars
        df_target - dataframe with target parameter info
        plot - save SpecMatch plots (boolean)
        savefolder - output directory name (String)
        scaleres - amount to scale residuals from composite spectrum (default 1)
    
    OUTPUT:
        stellar parameters teff, feh, logg, vsini, and their errors
        df_chi_total, LCS  
    
    EXAMPLE:
        
    """
    _log(verbose, 'Saving results to {}'.format(savefolder))
    utils.make_dir(savefolder, verbose=_verbosity_level(verbose) >= 2)
    targetname = Htarget.object
    vsini_constraint = vsini
    target_component = safe_filename_component(targetname)
    result_prefix = Path(savefolder) / target_component
    target_state_receipt = getattr(
        Htarget, '_neidspecmatch_loaded_state_sha256', None
    )
    if (
        target_state_receipt is not None
        and target_state_receipt != _spectrum_state_fingerprint(Htarget)
    ):
        raise ValueError(
            "Target spectrum was mutated after its NEIDSpecMatch load receipt."
        )
    # Establish the exact target/reference rest frames before manifest-state
    # verification and before the full-library ranking.  Stage two consumes
    # this state without applying another shift.
    _apply_matching_rv_overrides(
        Htarget, Hrefs, absrv=absrv, reference_rvs=reference_rvs
    )
    if not allow_unmanifested_library:
        if library_path is None:
            raise ValueError(
                "library_path is required so the Python API can validate the "
                "catalog, FITS list, sizes, and manifest."
            )
        library_root = Path(library_path).resolve()
        manifest_path = library_root / 'library_manifest.json'
        reusable_receipt = (
            isinstance(library_manifest, dict)
            and library_manifest.get('recognized') is True
            and library_manifest.get('deep_verified') is True
            and manifest_path.is_file()
            and library_manifest.get('manifest_sha256') == sha256_file(manifest_path)
        )
        if not reusable_receipt:
            with manifest_path.open(encoding='utf-8') as stream:
                manifest_catalog = json.load(stream).get('catalog')
            if not manifest_catalog:
                raise ValueError("Library manifest does not name its catalog.")
            library_manifest = _library_manifest_provenance(
                library_root / manifest_catalog, deep=True,
                expected_library_id=library_id,
            )
        library_manifest = dict(library_manifest)
        library_manifest.update(_verify_reference_pool_against_manifest(
            Htarget, Hrefs, library_root
        ))
    target_blaze_source, reference_blaze_state = (
        _refresh_matching_blaze_state(
            Htarget, Hrefs, input_is_deblazed=deblazed
        )
    )
    reference_rv_state = _reference_rv_state(Hrefs)
    if drp_compatibility is None:
        drp_compatibility = validate_drp_compatibility(
            Htarget, Hrefs, allow_mixed_drp=allow_mixed_drp
        )
    ##############################

    _log(verbose, 'Running chi-square library comparison')
    ##############################
    # STEP 1: Chi2 Loop
    df_chi, df_chi_best, Hbest, pairwise_common_mask = chi2spectraPolyLoop(
        ww, Htarget, Hrefs, plot_all=False, plot_chi=plot,
        verbose=verbose, maxvsini=maxvsini, order=order,
        deblazed=deblazed, vsini=vsini, vsini_window=vsini_window,
        return_comparison_mask=True,
    )
    ##############################
    # Combine best data
    # print(df_chi_best['OBJECT_ID'])
    # print(df_library['OBJECT_ID'])
    if (
        len(df_chi_best) != TOP_K
        or df_chi_best['OBJECT_ID'].duplicated().any()
    ):
        raise ValueError(
            f"SpecMatch requires exactly {TOP_K} unique best references."
        )
    if df_library['OBJECT_ID'].duplicated().any():
        raise ValueError("Library metadata OBJECT_ID values must be unique.")
    best_ids = df_chi_best['OBJECT_ID'].astype(str).tolist()
    df_chi_best_total = pd.merge(
        df_chi_best, df_library, on='OBJECT_ID', how='left', sort=False,
        validate='one_to_one',
    )
    if df_chi_best_total['OBJECT_ID'].astype(str).tolist() != best_ids:
        raise RuntimeError("Best-reference metadata merge changed component order.")
    required_parameters = ['Teff', '[Fe/H]', 'log(g)', 'vsini']
    missing_parameters = sorted(
        column for column in required_parameters
        if column not in df_chi_best_total.columns
    )
    if missing_parameters:
        raise ValueError(f"Library metadata are missing columns: {missing_parameters}")
    numeric_parameters = df_chi_best_total[required_parameters].apply(
        pd.to_numeric, errors='coerce'
    )
    if not np.all(np.isfinite(numeric_parameters.to_numpy(dtype=float))):
        bad_rows = []
        finite = np.isfinite(numeric_parameters.to_numpy(dtype=float))
        for row_index in np.flatnonzero(~np.all(finite, axis=1)):
            bad_rows.append({
                'OBJECT_ID': str(
                    df_chi_best_total.iloc[row_index]['OBJECT_ID']
                ),
                'fields': [
                    required_parameters[column_index]
                    for column_index in np.flatnonzero(~finite[row_index])
                ],
            })
        raise ValueError(
            "Best-reference metadata contain missing/non-finite parameters: "
            f"{bad_rows}"
        )
    # print(df_chi_best_total['OBJECT_ID'])
    df_chi_total = pd.merge(
        df_chi, df_library, on='OBJECT_ID', how='left', sort=False,
        validate='one_to_one',
    )
    if plot:
        plot_chi_teff_feh_logg_panel(df_chi_total.rss,
                                     df_chi_total.Teff,
                                     df_chi_total['[Fe/H]'],
                                     df_chi_total['log(g)'],
                                     savename=os.fspath(result_prefix) + '_chi2_3panel.png')
        plot_teff_feh_logg_corr_panel(df_chi_total.Teff.values,
                                      df_chi_total['[Fe/H]'].values,
                                      df_chi_total['log(g)'].values,
                                      df_chi_total.rss.values,
                                      e_teff=df_chi_total['e_Teff'],
                                      e_feh=df_chi_total['e_[Fe/H]'],
                                      e_logg=df_chi_total['e_log(g)'],
                                      savename=os.fspath(result_prefix) + '_chi2_corrpanel.png')
    if df_target is not None:
        teff_known = df_target.Teff.values[0]
        tefferr_known = df_target.e_Teff.values[0]
        feh_known = df_target['[Fe/H]'].values[0]
        feherr_known = df_target['e_[Fe/H]'].values[0]
        logg_known = df_target['log(g)'].values[0]
        loggerr_known = df_target['e_log(g)'].values[0]
    else:
        teff_known = np.nan
        tefferr_known = np.nan
        feh_known = np.nan
        feherr_known = np.nan
        logg_known = np.nan
        loggerr_known = np.nan

    _log(verbose, 'Performing linear-combination optimization')

    ##############################
    # STEP 2 LINEAR COMBINATION
    ##############################
    f1, e1, ffrefs, eerefs, rv_metadata = get_data_ready(
        Htarget, Hbest, ww, v, polyvals=df_chi_best.poly_params.values,
        vsinis=df_chi_best.vsini.values, order=order,
        deblazed=deblazed, plot=plot, verbose=verbose,
        return_rv_metadata=True,
    )
    if feh_calibration is not None:
        feh_calibration = _validate_calibration_runtime_context(
            feh_calibration,
            library_id=library_id,
            library_manifest=(library_manifest or {}),
            drp_compatibility=drp_compatibility,
            blaze_source=getattr(
                Htarget, 'blaze_source_used',
                getattr(Htarget, 'provenance', {}).get('blaze_source'),
            ),
        )
    L = LPFunctionLinComb(
        ww, f1, e1, ffrefs, eerefs, valid_mask=pairwise_common_mask
    )
    LCS = FitLinCombSpec(L, df_chi_best_total.Teff.values,
                         df_chi_best_total['[Fe/H]'].values,
                         df_chi_best_total['log(g)'].values,
                         df_chi_best_total['vsini'].values,
                         teff_known=teff_known,
                         tefferr_known=tefferr_known,
                         feh_known=feh_known,
                         feherr_known=feherr_known,
                         logg_known=logg_known,
                         loggerr_known=loggerr_known,
                         targetname=targetname,
                         calibrate_feh=calibrate_feh,
                         feh_calibration=feh_calibration,
                         order=order,
                         library_id=library_id,
                         calibration_population=calibration_population,
                         calibration_population_teff=calibration_population_teff,
                         verbose=verbose)
    LCS.minimize_PyDE(mcmc=False, random_seed=random_seed, verbose=verbose)
    all_pairwise_diagnostics = df_chi['optimizer_diagnostics'].tolist()
    selected_pairwise_diagnostics = df_chi_best['optimizer_diagnostics'].tolist()
    LCS.pairwise_optimizer_diagnostics = {
        'all_references_success': bool(all(
            item.get('success') is True
            and item.get('linear_rank') == 6
            and np.isfinite(item.get('objective_rss', np.nan))
            for item in all_pairwise_diagnostics
        )),
        'reference_count': len(all_pairwise_diagnostics),
        'failed_reference_count': int(sum(
            not (
                item.get('success') is True
                and item.get('linear_rank') == 6
                and np.isfinite(item.get('objective_rss', np.nan))
            )
            for item in all_pairwise_diagnostics
        )),
        'selected_references': selected_pairwise_diagnostics,
        'common_valid_pixel_count': int(
            selected_pairwise_diagnostics[0].get('valid_pixel_count', 0)
        ) if selected_pairwise_diagnostics else 0,
        'common_valid_pixel_fraction': float(
            selected_pairwise_diagnostics[0].get('valid_pixel_fraction', 0.0)
        ) if selected_pairwise_diagnostics else 0.0,
    }
    LCS.rv_metadata = rv_metadata
    weights = LCS.lpf.get_pv_all(LCS.min_pv)
    LCS.calculate_stellar_parameters(weights)
    LCS.support_diagnostics = _build_support_diagnostics(
        df_chi_total,
        df_chi_best_total,
        weights,
        optimizer_nonunique=LCS.optimizer_diagnostics.get(
            'nonunique_optimum', False
        ),
    )
    estimator_config = _build_estimator_config(
        order=order,
        wavelength=ww,
        maxvsini=maxvsini,
        vsini=vsini_constraint,
        vsini_window=vsini_window,
        library_id=library_id,
        library_manifest=library_manifest,
        deblazed=deblazed,
        blaze_source=target_blaze_source,
        reference_rv_state=reference_rv_state,
        reference_blaze_state=reference_blaze_state,
    )
    # LCS.plot_model(LCS.min_pv)# SEJ
    if plot:
        LCS.plot_model_with_components(LCS.min_pv, names=df_chi_best['OBJECT_ID'].values, title="", mode=mode,
                                       savename=os.fspath(result_prefix) + '_compositecomparison.png', scaleres=scaleres)
    if save_plot_data:
        LCS.plot_model_with_component_save(LCS.min_pv, names=df_chi_best['OBJECT_ID'].values, title="", mode=mode,
                                       datafile=os.fspath(result_prefix) + '_compositecomparison.npz', scaleres=scaleres)

    teff = LCS.teff
    feh = LCS.feh
    logg = LCS.logg
    vsini = LCS.vsini
    teff_delta = teff - teff_known
    feh_delta = feh - feh_known
    logg_delta = logg - logg_known

    provenance = build_result_provenance(
        target=Htarget,
        references=Hrefs,
        rv_metadata=rv_metadata,
        order=order,
        wavelength=ww,
        deblazed=deblazed,
        vsini=vsini_constraint,
        vsini_window=vsini_window,
        random_seed=random_seed,
        feh_calibration_metadata=LCS.feh_calibration_metadata,
        drp_compatibility=drp_compatibility,
        library_id=library_id,
        library_manifest=library_manifest,
        dataset_id=dataset_id,
        validation_summary=validation_summary,
        validation_population=validation_population,
        validation_population_teff=validation_population_teff,
        composite_optimizer_diagnostics=LCS.optimizer_diagnostics,
        pairwise_optimizer_diagnostics=LCS.pairwise_optimizer_diagnostics,
        estimator_config=estimator_config,
        support_diagnostics=LCS.support_diagnostics,
    )
    results = {'schema_version': RESULT_SCHEMA_VERSION,
               'target_name': str(targetname),
               'target_safe_identifier': target_component,
               'teff': teff,
               'logg': logg,
               'feh': feh,
               'feh_raw': LCS.feh_raw,
               'feh_calibration': LCS.feh_calibration_metadata,
               'order': int(order),
               'library_id': library_id,
               'rv_provenance': rv_metadata,
               'provenance': provenance,
               'reference_support': LCS.support_diagnostics,
               'estimator_config_sha256': _hash_estimator_config(
                   estimator_config
               ),
               'vsini': vsini,
               'weights': weights}
    json_path = os.fspath(result_prefix) + '_results.json'
    with open(json_path, 'w', encoding='utf-8') as stream:
        json.dump(results, stream, indent=2, sort_keys=True, default=_json_default)
        stream.write('\n')
    if save_legacy_pickle:
        pickle_path = os.fspath(result_prefix) + '_results.pkl'
        with open(pickle_path, "wb") as savefile:
            pickle.dump(results, savefile, protocol=pickle.HIGHEST_PROTOCOL)
        _log(verbose, 'Saved trusted-files-only legacy pickle to {}'.format(pickle_path))
    _log(verbose, 'Saved results to {}'.format(json_path))

    _log(verbose, 'Saving chi-square dataframe to {}'.format(
        os.fspath(result_prefix) + '_chi2results.csv'), level=2)
    df_chi_total.to_csv(os.fspath(result_prefix) + '_chi2results.csv', index=False)
    df_chi_best_total['weights'] = weights
    df_chi_best_total.to_csv(
        os.fspath(result_prefix) + '_chi2results_best.csv', index=False
    )
    return teff, feh, logg, vsini, teff_delta, feh_delta, logg_delta, df_chi_total, LCS


def plot_chi_teff_feh_logg_panel(chis, teff, feh, logg, savename='chi2panel.pdf', fig=None, title=''):
    """
    3panel plot chi2 vs Teff, FeH, logg for 5 best fit stars among all library stars
    
    INPUT:
        chis - (array)
        teff - (array)
        feh - (array)
        logg - (array)
        savename - directory name to save plot (str)
        title - plot title (str)

    OUTPUT:
        saves 3panel plot
            
    EXAMPLE:
        plot_chi_teff_feh_logg_panel(df_chi.chi2,df_chi.Teff,df_chi['[Fe/H]'],df_chi['log(g)'],savename='chi2.pdf')
    """
    if fig is None:
        fig, (ax, bx, cx) = plt.subplots(dpi=200, ncols=3, sharey=True, figsize=(6, 2))
    df = pd.DataFrame(list(zip(chis, teff, feh, logg)), columns=['chi2', 'teff', 'feh', 'logg'])
    df_best = df.sort_values('chi2').reset_index(drop=True)[0:5]
    ax.plot(teff, chis, marker='o', lw=0)
    ax.plot(df_best.teff, df_best.chi2, color='crimson', marker='h', lw=0)
    # ax.set_xlim(2950.,4050)
    # ax.set_xticks(np.arange(3000,4100,200))

    ax.set_yscale('log')
    ax.set_ylabel(r'$\chi^2$', fontsize=10, labelpad=2)
    ax.set_xlabel(r'$T_{\mathrm{eff}}$ [K]', fontsize=10, labelpad=2)
    bx.set_xlabel('[Fe/H] ', fontsize=10, labelpad=2)
    cx.set_xlabel(r'$\log(g)$', fontsize=10, labelpad=2)

    bx.plot(feh, chis, marker='o', lw=0)
    bx.plot(df_best.feh, df_best.chi2, color='crimson', marker='h', lw=0)
    # bx.set_xlim(-0.55,0.55)
    # bx.set_xticks(np.arange(-0.5,0.55,0.25))

    cx.plot(logg, chis, marker='o', lw=0)
    cx.plot(df_best.logg, df_best.chi2, color='crimson', marker='h', lw=0)
    # cx.set_xlim(4.7,5.1)

    for xx in (ax, bx, cx):
        utils.ax_apply_settings(xx, ticksize=7)
        xx.tick_params(pad=2)
        # xx.set_ylim(3e3,1.5e5)
    bx.set_title(title, fontsize=14)
    fig.subplots_adjust(wspace=0.1, right=0.97, left=0.08, bottom=0.15, top=0.95)
    fig.savefig(savename, dpi=200)
    print('Saved to: {}'.format(savename))


def plot_teff_feh_logg_corr_panel(teff, feh, logg, chis, e_teff=None, e_feh=None, e_logg=None, fig=None, ax=None,
                                  bx=None,
                                  savename='corr_panel.pdf', scale_factor_for_points=1500.):
    """
    Plot of teff vs feh, logg for 5 best fit stars among all library stars
    
    INPUT:
        teff - (array)
        feh - (array)
        logg - (array)
        chis - (array)
        savename - directory name to save plot (str)
        title - plot title (str)

    OUTPUT:
        saves correlation plot
        
    """
    if ax is None and bx is None and fig is None:
        fig, (ax, bx) = plt.subplots(ncols=2, dpi=200, sharey=True, figsize=(10, 3))
    colors = utils.get_cmap_colors(p=chis)
    df_chi = pd.DataFrame(zip(teff, feh, logg, chis), columns=['Teff', '[Fe/H]', 'log(g)', 'chi2'])
    df_chi = df_chi.sort_values('chi2').reset_index(drop=True)
    df_chi_best = df_chi[0:5]

    # First plot:
    ax.scatter(df_chi['[Fe/H]'], df_chi.Teff, s=scale_factor_for_points / df_chi.chi2.values ** 2., color='black')
    ax.scatter(df_chi_best['[Fe/H]'], df_chi_best.Teff, lw=0, color='crimson',
               s=scale_factor_for_points / df_chi_best.chi2.values ** 2.)

    if e_feh is not None and e_logg is not None and e_teff is not None:
        median_xerr = np.median(e_feh)
        median_yerr = np.median(e_teff)
        ax.errorbar(-0.4, 3950, xerr=median_xerr, yerr=median_yerr, capsize=3, elinewidth=0.5, mew=0.5, color='black')
        ax.text(-0.40 + 0.0, 3820, 'Median error', fontsize=10, horizontalalignment='center')

    # Second plot
    bx.scatter(df_chi['log(g)'], df_chi.Teff, s=scale_factor_for_points / df_chi.chi2.values ** 2., color='black',
                label=r'Library spectra:' + '\n' + r'(radius$\propto 1/\chi^2$)')
    bx.scatter(df_chi_best['log(g)'], df_chi_best.Teff, lw=0, color='crimson',
               s=scale_factor_for_points / df_chi_best.chi2.values ** 2., label='5 Best-fit spectra')
    bx.legend(loc='upper right', fontsize=9, frameon=False)

    if e_feh is not None and e_logg is not None and e_teff is not None:
        median_xerr = np.median(e_logg)
        bx.errorbar(4.68, 3230, xerr=median_xerr, yerr=median_yerr, capsize=3, elinewidth=0.5, mew=0.5, color='black')
        bx.text(4.685, 3100, 'Median error', fontsize=10, horizontalalignment='center')

    ax.set_ylabel(r'$T_{\mathrm{eff}}$ [K]', fontsize=12, labelpad=2)
    ax.set_xlabel('Fe/H [dex]', fontsize=12, labelpad=1)
    bx.set_xlabel(r'$\log (g)$ [dex]', fontsize=12, labelpad=1)

    fig.subplots_adjust(wspace=0.05)

    for xx in (ax, bx):
        utils.ax_apply_settings(xx, ticksize=9)
        xx.grid(lw=0.4, alpha=0.3)
        utils.ax_set_linewidth(xx, 1.3)

    fig.savefig(savename, dpi=200)
    print('Saved to: {}'.format(savename))


def summarize_values_from_orders(
        files, targetname, *, trusted_pickle=False, output_directory=None,
        verbose=1):
    """
    Summarize values from different orders from SpecMatch analysis

    INPUT:
        Names of JSON result files. Legacy pickle files require the explicit
        ``trusted_pickle=True`` acknowledgement because pickle can execute
        arbitrary code.

    OUTPUT:
        saves two .csv files:
            1) a .csv file with the Teff, Fe/H, and logg values from all of the orders
            2) a .csv file with the median and the standard deviation from different orders

    Example:
        files = sorted(glob.glob('results/TIC*/*_results.json'))
        summarize_values_from_orders(files,'AD_Leo')
    """
    if not files:
        raise ValueError("At least one result file is required.")
    paths = [Path(filename) for filename in files]
    savefolder = (
        Path(output_directory) if output_directory is not None
        else paths[0].resolve().parent.parent
    )
    savefolder.mkdir(parents=True, exist_ok=True)
    target = safe_filename_component(targetname)
    params = ['teff', 'logg', 'feh', 'vsini']
    results = []
    for filename in paths:
        if filename.suffix.lower() == '.json':
            with filename.open(encoding='utf-8') as stream:
                res = json.load(stream)
        elif filename.suffix.lower() in {'.pkl', '.pickle'}:
            if not trusted_pickle:
                raise ValueError(
                    "Refusing legacy pickle input without trusted_pickle=True."
                )
            with filename.open('rb') as stream:
                res = pickle.load(stream, encoding='latin1')
        else:
            raise ValueError(f"Unsupported result format: {filename}")
        results.append(res)
    df = pd.DataFrame(results)
    df['filenames'] = [os.fspath(path) for path in paths]
    missing = sorted(set(params).difference(df.columns))
    if missing:
        raise ValueError(f"Result files are missing parameters: {missing}")
    medians = [pd.to_numeric(df[name], errors='coerce').median() for name in params]
    stds = [pd.to_numeric(df[name], errors='coerce').std() for name in params]

    df_med = pd.DataFrame(list(zip(params, medians, stds)), columns=['parameters', 'median', 'std'])

    overview_path = savefolder / f'{target}_overview.csv'
    median_path = savefolder / f'{target}_med.csv'
    df.to_csv(overview_path, index=False)
    df_med.to_csv(median_path, index=False)

    _log(verbose, 'Saved to {}'.format(overview_path))
    _log(verbose, 'Saved to {}'.format(median_path))
    return df, df_med


def load_reference_library(
        df_library, path_df_lib_fits, verbose=1, orders=None):
    """Load reference spectra in catalog order with explicit RV provenance.

    DRP RV is preferred.  If and only if that RV fails validation, a finite
    catalog ``rv`` value is used.  Files are joined by basename rather than by
    position or target aliases.
    """
    fits_paths = sorted(glob.glob(os.path.join(path_df_lib_fits, '*.fits')))
    spectrum_order_kwargs = {}
    if orders is not None:
        selected_orders = [int(value) for value in orders]
        if not selected_orders:
            raise ValueError("orders cannot be empty when selecting a library slice")
        spectrum_order_kwargs = {
            'first_order_index': min(selected_orders),
            'stop_order_index': max(selected_orders) + 1,
        }
    by_basename = {}
    for path in fits_paths:
        basename = os.path.basename(path)
        if basename in by_basename:
            raise ValueError(f"Duplicate library FITS basename: {basename}")
        by_basename[basename] = path
    spectra = []
    try:
        for _, row in df_library.reset_index(drop=True).iterrows():
            basename = None
            for column in ('basenames', 'basename', 'filename'):
                value = row.get(column)
                if isinstance(value, str) and value.strip():
                    basename = os.path.basename(value.strip())
                    break
            if basename is None:
                raise ValueError(
                    "Every library row needs 'basenames', 'basename', or 'filename'."
                )
            if basename not in by_basename:
                raise FileNotFoundError(
                    f"Catalog row {row.get('OBJECT_ID')!r} maps to missing FITS {basename!r}."
                )
            filename = by_basename[basename]
            object_id = str(row['OBJECT_ID'])
            try:
                spectrum = neidspec.NEIDSpectrum(
                    filename, targetname=object_id, rv_source='drp',
                    verbose=_verbosity_level(verbose) >= 2,
                    **spectrum_order_kwargs,
                )
            except ValueError as exc:
                message = str(exc)
                catalog_rv = pd.to_numeric(
                    pd.Series([row.get('rv')]), errors='coerce'
                ).iloc[0]
                if (
                    not np.isfinite(catalog_rv)
                    or ('DRP RV' not in message and 'explicit rv' not in message)
                ):
                    raise
                spectrum = neidspec.NEIDSpectrum(
                    filename, targetname=object_id, rv_source='supplied',
                    rv=float(catalog_rv), verbose=_verbosity_level(verbose) >= 2,
                    **spectrum_order_kwargs,
                )
                spectrum.rv_source = 'library_catalog'
                spectrum.rv_uncertainty = np.nan
                provenance = getattr(spectrum, 'provenance', None)
                if isinstance(provenance, dict):
                    provenance['rv_source'] = 'library_catalog'
                    provenance['rv_catalog_column'] = 'rv'
                _log(verbose, 'Using catalog RV for {} ({})'.format(
                    object_id, basename
                ))
            spectrum.library_object_id = object_id
            # The catalog is the reviewed reference-library identity source.
            # Preserve the untouched archive header value separately; this is
            # essential for explicitly declared archive OBJECT typos such as
            # GJ3778 versus catalog object GJ 3378.
            spectrum.fits_object = str(getattr(spectrum, 'object', ''))
            spectrum.object = object_id
            provenance = getattr(spectrum, 'provenance', None)
            if isinstance(provenance, dict):
                provenance['library_object_id'] = object_id
                provenance['library_catalog_basename'] = basename
                provenance['fits_object'] = spectrum.fits_object
                provenance['library_identity_source'] = 'catalog_OBJECT_ID'
            spectrum._neidspecmatch_loaded_state_sha256 = (
                _spectrum_state_fingerprint(spectrum)
            )
            spectra.append(spectrum)
        if len(spectra) != len(df_library):
            raise RuntimeError("Library loader did not preserve the catalog row count.")
        return neidspec.NEIDSpecList(splist=spectra)
    except BaseException:
        for spectrum in spectra:
            _close_spectrum(spectrum)
        raise


def run_specmatch_for_orders(targetfile, targetname, outputdirectory='specmatch_results', HLS=None,
                             path_df_lib=None, path_df_lib_fits=None,
                             orders=None, maxvsini=30., calibrate_feh=False,
                             feh_calibration=None, calibration_population=None,
                             calibration_population_teff=None,
                             library_id=None, scaleres=1., deblazed=False,
                             mode='HR', save_plot_data=False, plot=False,
                             vsini=None, vsini_window=None, verbose=1,
                             absrv=None, rv_source='drp', reference_rvs=None,
                             random_seed=0, allow_mixed_drp=False,
                             dataset_id=None, save_legacy_pickle=False,
                             ccf_mask_path=None, ccf_mask_medium=None,
                             validation_summary=None,
                             allow_unmanifested_library=False,
                             ccf_orders=(55, 56, 91),
                             validation_population=None,
                             validation_population_teff=None):
    """
    run neidspecmatch for a given target file and orders
    
    INPUT:
        targetfile - name of target file
        targetname - target name, queried via simbad or tic ('GJ_251' or TIC 68581262)
        outputdirectory - folder to save overall results and plots
        HLS - refence stars as an HPFSpecList object, defaults to normal library
        path_df_lib - path to .csv file containing info on Teff/FeH/logg for all library stars
                    - defaults to config.PATH_LIBRARY_DB
        orders - hpf orders to run (orders 4, 5, 6, 14, 15, 16, and 17
                    recommended as they are the cleanest orders with minimal tellurics)
        maxvsini - maximum vsini to consider (default = 30 km/s)
    
    OUTPUT:
        result files will be saved to outputdirectory
    
    EXAMPLE:
        filename = '../input/20201020_hpf_gto_targets/Slope-20200114T091114_R01.optimal.fits'
        targetname = 'GJ_251'
        HLS = hpfspec.HPFSpecList(filelist=library_fitsfiles)
        outputdir = '20201020_hpf_gto_targets/GJ_251'
        run_specmatch_for_orders(filename, targetname , outputdir, HLS)
    
    NOTES:
        targetname will be queried via simbad or tic which saves a configuration file to target config directory
    
    """
    orders = ['102'] if orders is None else orders
    analysis_orders = [int(value) for value in orders]
    unknown_orders = sorted(
        value for value in analysis_orders if str(value) not in config.BOUNDS
    )
    if unknown_orders:
        raise ValueError(f"Unsupported NEID order indices: {unknown_orders}")
    if validation_summary is not None and len(analysis_orders) != 1:
        raise ValueError(
            "Publication validation is order/population specific. Run one "
            "order per invocation when supplying validation_summary."
        )
    if mode != 'HR':
        raise ValueError(
            "NEIDSpecMatch 0.2 supports NEID Level-2 high-resolution (HR) "
            "spectra only; HE mode has not been validated."
        )
    path_df_lib = path_df_lib or config.PATH_LIBRARY_DB
    path_df_lib_fits = path_df_lib_fits or config.PATH_LIBRARY_FITS
    if library_id is None:
        library_id = os.path.basename(os.path.dirname(os.path.abspath(path_df_lib)))

    if rv_source == 'supplied' and absrv is None:
        raise ValueError("rv_source='supplied' requires absrv.")
    if absrv is not None and rv_source != 'supplied':
        raise ValueError("absrv requires rv_source='supplied'.")
    if rv_source in {'custom_ccf', 'mask_ccf'} and ccf_mask_path is None:
        raise ValueError("An explicit ccf_mask_path is required for a custom CCF RV.")
    if ccf_mask_path is not None and rv_source not in {'custom_ccf', 'mask_ccf'}:
        raise ValueError(
            "ccf_mask_path is only valid with rv_source='custom_ccf'."
        )
    if ccf_mask_medium is not None and rv_source not in {'custom_ccf', 'mask_ccf'}:
        raise ValueError(
            "ccf_mask_medium is only valid with rv_source='custom_ccf'."
        )
    target_kwargs = {
        'targetname': targetname,
        'rv_source': rv_source,
        'verbose': _verbosity_level(verbose) >= 2,
    }
    target_loaded_orders = list(analysis_orders)
    if rv_source in {'custom_ccf', 'mask_ccf'}:
        target_loaded_orders.extend(int(value) for value in ccf_orders)
        target_kwargs['ccf_orders'] = tuple(int(value) for value in ccf_orders)
    target_kwargs['first_order_index'] = min(target_loaded_orders)
    target_kwargs['stop_order_index'] = max(target_loaded_orders) + 1
    if rv_source == 'supplied':
        target_kwargs['rv'] = float(absrv)
    if ccf_mask_path is not None:
        target_kwargs['ccf_mask_path'] = os.fspath(ccf_mask_path)
    if ccf_mask_medium is not None:
        target_kwargs['ccf_mask_medium'] = ccf_mask_medium
    # NEIDSpectrum establishes a single, provenance-bearing RV here.  The
    # matcher preserves it through every order instead of running another CCF.
    owned_hls = HLS is None
    Htarget = None
    try:
        Htarget = neidspec.NEIDSpectrum(targetfile, **target_kwargs)
        Htarget._neidspecmatch_loaded_state_sha256 = (
            _spectrum_state_fingerprint(Htarget)
        )
        _log(verbose, 'Reading library database from: {}'.format(path_df_lib))
        df_lib = pd.read_csv(path_df_lib)
        library_manifest = _library_manifest_provenance(
            path_df_lib, tolerate_invalid=allow_unmanifested_library,
            expected_library_id=library_id,
            validate_fits_schema=True,
        )
        if not allow_unmanifested_library and not library_manifest.get('recognized'):
            raise ValueError(
                "A recognized, per-file-checksummed library manifest is required."
            )

        vsinis = []
        target_component = safe_filename_component(Htarget.object)
        for o in orders:
            o = str(o)
            order_hls = HLS
            if owned_hls:
                _log(verbose, 'Loading reference library only for order {}'.format(o))
                order_hls = load_reference_library(
                    df_lib, path_df_lib_fits, verbose=verbose, orders=[int(o)]
                )
            try:
                drp_compatibility = validate_drp_compatibility(
                    Htarget, order_hls.splist,
                    allow_mixed_drp=allow_mixed_drp,
                )
                if not drp_compatibility['publication_validated']:
                    _log(
                        verbose,
                        'WARNING: mixed/unknown target/library DRP versions; '
                        'results are marked unvalidated.',
                    )
                _log(verbose, "Order {}".format(o))
                wmin, wmax = config.BOUNDS[o]
                ww = np.arange(wmin, wmax, 0.01)
                v = np.linspace(-175, 175, 2501)
                savefolder = Path(outputdirectory) / f'{target_component}_o{o}'
                _, _, _, fitted_vsini, _, _, _, _, _ = run_specmatch(
                    Htarget, order_hls.splist, ww, v, df_lib,
                    savefolder=os.fspath(savefolder), maxvsini=maxvsini,
                    calibrate_feh=calibrate_feh,
                    feh_calibration=feh_calibration,
                    calibration_population=calibration_population,
                    calibration_population_teff=calibration_population_teff,
                    library_id=library_id, scaleres=scaleres, order=int(o),
                    deblazed=deblazed, mode=mode,
                    save_plot_data=save_plot_data,
                    save_legacy_pickle=save_legacy_pickle, plot=plot,
                    vsini=vsini, vsini_window=vsini_window,
                    absrv=None, reference_rvs=reference_rvs,
                    verbose=verbose, random_seed=random_seed,
                    allow_mixed_drp=allow_mixed_drp,
                    drp_compatibility=drp_compatibility,
                    library_manifest=library_manifest, dataset_id=dataset_id,
                    validation_summary=validation_summary,
                    validation_population=validation_population,
                    validation_population_teff=validation_population_teff,
                    allow_unmanifested_library=allow_unmanifested_library,
                    library_path=Path(path_df_lib).resolve().parent,
                )
                vsinis.append(fitted_vsini)
            finally:
                if owned_hls and order_hls is not None:
                    _close_spectrum_collection(order_hls)
        return vsinis
    finally:
        if Htarget is not None:
            _close_spectrum(Htarget)

def plot_crossvalidation_results_1d(order, df_crossval, savefolder):
    """
    
    """
    fig, (ax, bx, cx) = plt.subplots(nrows=3, sharex=True, dpi=200)
    x = range(len(df_crossval))
    label = 'std={:0.2f}K'.format(np.std(df_crossval.d_teff, ddof=1))
    ax.plot(df_crossval.d_teff, 'k.', markersize=8)
    ax.set_ylim(ax.get_ylim()[0] * 1.6, ax.get_ylim()[1] * 2.2)

    label = 'std={:0.2f}dex'.format(np.std(df_crossval.d_feh, ddof=1))
    bx.plot(df_crossval.d_feh, 'k.', markersize=8)
    bx.set_ylim(bx.get_ylim()[0] * 2.0, bx.get_ylim()[1] * 2.2)

    label = 'std={:0.2f}dex'.format(np.std(df_crossval.d_logg, ddof=1))
    cx.plot(df_crossval.d_logg, 'k.', markersize=8)
    cx.set_ylim(cx.get_ylim()[0] * 2.2, cx.get_ylim()[1] * 2.2)
    for xx in (ax, bx, cx):
        utils.ax_apply_settings(xx)
        xx.tick_params(labelsize=9, pad=2)
        xx.grid(lw=0.3)
        # xx.legend(bbox_to_anchor=(1.,1.),fontsize=12)
    fig.subplots_adjust(hspace=0.05)
    ax.set_ylabel(r'$\Delta$Teff [K]', fontsize=12, labelpad=-4)
    bx.set_ylabel(r'$\Delta$Fe/H [dex]', fontsize=12, labelpad=-2)
    cx.set_ylabel(r'$\Delta$logg [dex]', fontsize=12, labelpad=-2)
    cx.set_xlabel('Library spectrum #', fontsize=12)
    ax.set_title('Library Performance (NEID Order {})'.format(order))
    fig.savefig(Path(savefolder) / f'crossvalidation_o{order}_plot1D.png')
    plt.close(fig)


def plot_crossvalidation_results_2d(order, df_crossval, savefolder):
    """
    """
    dd = df_crossval[["d_teff", "d_feh", "d_logg"]]
    fig, axx = plt.subplots(nrows=3, ncols=3, dpi=200, figsize=(6, 6))

    for xx in [axx[0, 1], axx[0, 2], axx[1, 2]]:
        xx.axes.set_axis_off()

    labels = [r"$\Delta T_{\mathrm{eff}}$", r"$\Delta$[Fe/H]", r"$\Delta \log g$"]
    xlabels = [r"$\Delta T_{\mathrm{eff}}$ [K]", r"$\Delta$[Fe/H] [dex]", r"$\Delta \log g$ [dex]"]

    # xlims = [xx.get_xlim() for xx in [axx[1,0],axx[2,0],axx[2,1]]]
    xlims = []
    diag = [axx[0, 0], axx[1, 1], axx[2, 2]]
    for i, xx in enumerate(diag):
        xx.set_title(labels[i], fontsize=12)
        _y, _x, _ = xx.hist(dd.iloc[:, i].values, color="black", histtype="step", bins=7,
                            range=[dd.iloc[:, i].values.min(), dd.iloc[:, i].values.max()])
        g = astropy.modeling.models.Gaussian1D(stddev=np.std(dd.iloc[:, i].values))
        xmax = np.max(np.abs(dd.iloc[:, i].values))
        x = np.linspace(-xmax, xmax, 1000)
        y = (g(x) / (np.max(g(x)))) * (np.max(_y))
        xx.plot(x, y, color="crimson", lw=1, ls="--")

        if i == 0:
            xx.xaxis.set_visible(False)
            xx.yaxis.set_visible(False)
            # xx.set_xlim(xlims[0][0],xlims[0][1])
        if i == 1:
            xx.xaxis.set_visible(False)
            xx.yaxis.set_visible(False)
            # xx.set_xlim(xlims[2][0],xlims[2][1])
        if i == 2:
            xx.yaxis.set_visible(False)
        xlims.append(xx.get_xlim())

    axx[1, 0].plot(df_crossval.d_teff.values, df_crossval.d_feh.values, marker="o", lw=0, color="k", alpha=0.8)
    axx[1, 0].set_xlim(xlims[0][0], xlims[0][1])
    axx[2, 0].plot(df_crossval.d_teff.values, df_crossval.d_logg.values, marker="o", lw=0, color="k", alpha=0.8)
    axx[2, 0].set_xlim(xlims[0][0], xlims[0][1])
    axx[2, 1].plot(df_crossval.d_feh.values, df_crossval.d_logg.values, marker="o", lw=0, color="k", alpha=0.8)
    axx[2, 1].set_xlim(xlims[1][0], xlims[1][1])

    axx[2, 1].yaxis.set_visible(False)
    axx[1, 0].xaxis.set_visible(False)

    axx[1, 0].set_ylabel(xlabels[1], fontsize=12, labelpad=3)
    axx[2, 0].set_ylabel(xlabels[2], fontsize=12, labelpad=-2)
    axx[2, 0].set_xlabel(xlabels[0], fontsize=12)
    axx[2, 1].set_xlabel(xlabels[1], fontsize=12)
    axx[2, 2].set_xlabel(xlabels[2], fontsize=12)

    for xx in axx.flatten():
        utils.ax_apply_settings(xx, ticksize=10)
        xx.grid(lw=0)
        xx.tick_params(pad=1)

    fig.subplots_adjust(wspace=0.02, hspace=0.02)
    fig.savefig(Path(savefolder) / f'crossvalidation_o{order}_plot2D.png')
    plt.close(fig)


def run_crossvalidation_for_orders(
        order, df_lib=None, HLS=None, outputdir=None,
        path_df_lib_fits=None, plot_results=False, plot_each_target=False,
        calibrate_feh=False, fit_feh_calibration=False,
        population_rules=None, library_id=None, scaleres=1., verbose=1,
        min_calibration_stars=8, random_seed=0, allow_mixed_drp=False,
        resume=True, allow_unmanifested_library=False, maxvsini=30.0):
    """Run spectral leave-one-out validation for one NEID order.

    Spectral fits are always performed on the raw empirical [Fe/H] scale.
    Setting ``fit_feh_calibration=True`` cross-fits a leave-one-out calibration
    layer, writes its validation residuals and metadata artifact, and reports
    the corresponding calibrated scatter.  The artifact records that the
    spectral-LOO predictions are correlated; this is not strict nested CV.

    Plotting is opt-in.  ``verbose=0`` is suitable for long batch jobs.
    """
    if calibrate_feh:
        raise UnvalidatedCalibrationError(
            "Cross-validation must recover raw metallicities first. Use "
            "fit_feh_calibration=True to validate the complete estimator."
        )
    order = str(order)
    df_lib_path = None
    df_lib_source = config.PATH_LIBRARY_DB if df_lib is None else df_lib
    if isinstance(df_lib_source, (str, os.PathLike)):
        df_lib_path = os.path.abspath(os.fspath(df_lib_source))
        df_lib = pd.read_csv(df_lib_path)
        if library_id is None:
            library_id = os.path.basename(os.path.dirname(df_lib_path))
    else:
        df_lib = df_lib_source.copy()
    if not library_id:
        library_id = "unspecified-library"
    path_df_lib_fits = path_df_lib_fits or config.PATH_LIBRARY_FITS
    outputdir = outputdir or config.PATH_LIBRARY_CROSSVAL
    outputdir = os.path.join(os.fspath(outputdir), 'o{}_crossval'.format(order))
    os.makedirs(outputdir, exist_ok=True)
    owned_hls = HLS is None
    try:
        if HLS is None:
            _log(verbose, 'Loading reference spectra from {}'.format(path_df_lib_fits))
            HLS = load_reference_library(
                df_lib, path_df_lib_fits, verbose=verbose,
                orders=[int(order)],
            )
        return _run_crossvalidation_loaded(
            order=order, df_lib=df_lib, df_lib_path=df_lib_path, HLS=HLS,
            outputdir=outputdir, plot_results=plot_results,
            plot_each_target=plot_each_target,
            fit_feh_calibration=fit_feh_calibration,
            population_rules=population_rules, library_id=library_id,
            scaleres=scaleres, verbose=verbose,
            min_calibration_stars=min_calibration_stars,
            random_seed=random_seed, allow_mixed_drp=allow_mixed_drp,
            resume=resume,
            allow_unmanifested_library=allow_unmanifested_library,
            maxvsini=maxvsini,
        )
    finally:
        if owned_hls and HLS is not None:
            _close_spectrum_collection(HLS)


def _atomic_json_write(filename, payload):
    path = Path(filename)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f'.{path.name}-', suffix='.tmp', dir=path.parent
    )
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, default=_json_default)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _atomic_dataframe_write(frame, filename):
    path = Path(filename)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f'.{path.name}-', suffix='.tmp', dir=path.parent
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary_name, index=False)
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _checkpoint_spectrum_identity(spectrum):
    filename = getattr(spectrum, 'filename', None)
    path = Path(filename) if filename else None
    return {
        'object_id': str(getattr(spectrum, 'object', '')),
        'basename': path.name if path is not None else None,
        'file_sha256': sha256_file(path) if path is not None and path.is_file() else None,
        'loaded_state_sha256': _spectrum_state_fingerprint(spectrum),
        'drp_version': _drp_version(spectrum),
        'dq_status': _spectrum_quality_status(spectrum),
        'loaded_order_range': [
            getattr(spectrum, 'start', None), getattr(spectrum, 'end', None)
        ],
    }


def _validate_checkpoint_folds(records, identities):
    if not isinstance(records, list) or not isinstance(identities, list):
        raise ValueError(
            "Cross-validation checkpoint folds and identities must be lists."
        )
    if any(not isinstance(item, dict) for item in records):
        raise ValueError(
            "Every cross-validation checkpoint fold must be a mapping."
        )
    indices = [item.get('index') for item in records]
    if len(indices) != len(set(indices)):
        raise ValueError("Cross-validation checkpoint contains duplicate folds.")
    for record in records:
        index = record.get('index')
        if not isinstance(index, int) or not 0 <= index < len(identities):
            raise ValueError("Cross-validation checkpoint has an invalid fold index.")
        if record.get('target_identity') != identities[index]:
            raise ValueError(
                "Cross-validation checkpoint fold identity/order does not match "
                "the loaded library."
            )
        values = record.get('values')
        if (
            not isinstance(values, list) or len(values) != 7
            or not np.all(np.isfinite(np.asarray(values, dtype=float)))
        ):
            raise ValueError(
                "Cross-validation checkpoint fold must contain seven finite values."
            )
        optimizer = record.get('optimizer')
        if (
            not isinstance(optimizer, dict)
            or optimizer.get('success') is not True
        ):
            raise ValueError(
                "Cross-validation checkpoint contains an unsuccessful optimizer fold."
            )
        support_status = _support_diagnostics_status(record.get('support'))
        if not support_status['valid']:
            raise ValueError(
                "Cross-validation checkpoint fold has incomplete support "
                f"diagnostics ({support_status['reason']})."
            )


def _run_crossvalidation_loaded(
        *, order, df_lib, df_lib_path, HLS, outputdir, plot_results,
        plot_each_target, fit_feh_calibration, population_rules, library_id,
        scaleres, verbose, min_calibration_stars, random_seed,
        allow_mixed_drp, resume, allow_unmanifested_library, maxvsini):
    if len(HLS.splist) != len(df_lib):
        raise ValueError(
            "The spectrum list and library table must contain the same number "
            "of stars in the same order."
        )
    if not HLS.splist:
        raise ValueError("Cross-validation requires a nonempty spectrum list.")
    drp_compatibility = validate_drp_compatibility(
        HLS.splist[0], HLS.splist[1:], allow_mixed_drp=allow_mixed_drp
    )
    if not drp_compatibility['publication_validated']:
        _log(verbose, 'WARNING: mixed/unknown DRP metadata; CV is unvalidated.')
    library_manifest = (
        _library_manifest_provenance(
            df_lib_path, deep=True, expected_library_id=library_id,
            validate_fits_schema=True,
        )
        if df_lib_path is not None
        else {'manifest_present': False, 'library_id': library_id,
              'manifest_sha256': None, 'recognized': False,
              'deep_verified': False, 'reference_pool_verified': False}
    )
    if not allow_unmanifested_library and not library_manifest.get('recognized'):
        raise ValueError(
            "A recognized, per-file-checksummed library manifest is required "
            "for cross-validation."
        )
    if not allow_unmanifested_library:
        library_manifest = dict(library_manifest)
        library_manifest.update(_verify_reference_pool_against_manifest(
            HLS.splist[0], HLS.splist[1:], Path(df_lib_path).parent
        ))
    catalog_sha256 = (
        sha256_file(df_lib_path) if df_lib_path is not None else None
    )
    spectrum_identities = [
        _checkpoint_spectrum_identity(item) for item in HLS.splist
    ]
    wmin, wmax = config.BOUNDS[str(order)]
    ww = np.arange(wmin, wmax, 0.01)
    v = np.linspace(-175, 175, 2501)
    blaze_source = getattr(
        HLS.splist[0], 'blaze_source_used',
        getattr(HLS.splist[0], 'blaze_source', 'unknown'),
    )
    reference_rv_state = _reference_rv_state(HLS.splist)
    reference_blaze_state = _reference_blaze_state(HLS.splist)
    estimator_config = _build_estimator_config(
        order=int(order), wavelength=ww, maxvsini=maxvsini,
        vsini=None, vsini_window=None, library_id=library_id,
        library_manifest=library_manifest, deblazed=False,
        blaze_source=blaze_source,
        reference_rv_state=reference_rv_state,
        reference_blaze_state=reference_blaze_state,
    )
    estimator_config_json = _canonical_estimator_config_json(
        estimator_config
    )
    estimator_config_sha256 = _hash_estimator_config(estimator_config)
    checkpoint_key = {
        'checkpoint_schema_version': 4,
        'result_schema_version': RESULT_SCHEMA_VERSION,
        'neidspecmatch_version': __version__,
        'neidspec_version': getattr(neidspec, '__version__', 'unknown'),
        'neidspec_source_sha256': _neidspec_source_fingerprint(),
        'pipeline_source_sha256': _pipeline_source_fingerprint(),
        'selection_metric': 'unweighted_residual_sum_of_squares',
        'order': int(order),
        'library_id': str(library_id),
        'library_manifest_sha256': library_manifest.get('manifest_sha256'),
        'catalog_sha256': catalog_sha256,
        'drp_versions': drp_compatibility['versions'],
        'blaze_source': blaze_source,
        'estimator_config_schema_version': ESTIMATOR_CONFIG_SCHEMA_VERSION,
        'estimator_config_sha256': estimator_config_sha256,
        'estimator_config_json': estimator_config_json,
        'random_seed': int(random_seed),
        'spectrum_identities_in_order': spectrum_identities,
        'reference_dq_status_counts': _quality_counts(HLS.splist),
    }
    checkpoint_path = Path(outputdir) / f'crossvalidation_checkpoint_o{order}.json'
    fold_records = {}
    if resume and checkpoint_path.is_file():
        with checkpoint_path.open(encoding='utf-8') as stream:
            checkpoint = json.load(stream)
        if checkpoint.get('key') != checkpoint_key:
            raise ValueError(
                "Cross-validation checkpoint is incompatible with this order, "
                "library/manifest, DRP set, seed, or code version. Use "
                "resume=False/--no-resume to start a fresh run."
            )
        checkpoint_folds = checkpoint.get('folds', [])
        _validate_checkpoint_folds(checkpoint_folds, spectrum_identities)
        fold_records = {int(item['index']): item for item in checkpoint_folds}
        _log(verbose, f'Resuming {len(fold_records)}/{len(df_lib)} completed folds')

    columns = ['teff', 'feh', 'logg', 'vsini', 'd_teff', 'd_feh', 'd_logg']
    _log(verbose, 'Running order {} for {} leave-one-out folds'.format(
        order, len(df_lib)))
    for index in range(len(df_lib)):
        if index in fold_records:
            continue
        target = HLS.splist[index]
        _log(verbose, '[{}/{}] {}'.format(index + 1, len(df_lib), target.object))
        target_row = df_lib[df_lib['OBJECT_ID'] == target.object]
        if len(target_row) != 1:
            raise ValueError(
                "Expected exactly one metadata row for {!r}; found {}.".format(
                    target.object, len(target_row))
            )
        references = np.delete(np.array(HLS.splist), index)
        result = run_specmatch(
            target, references, ww, v, df_lib, target_row,
            plot=plot_each_target,
            savefolder=os.fspath(Path(outputdir) / 'plots'),
            calibrate_feh=False, library_id=library_id,
            scaleres=scaleres, order=int(order), maxvsini=maxvsini,
            verbose=max(0, _verbosity_level(verbose) - 1),
            random_seed=random_seed + index,
            allow_mixed_drp=allow_mixed_drp,
            drp_compatibility=drp_compatibility,
            library_manifest=library_manifest,
            allow_unmanifested_library=allow_unmanifested_library,
            library_path=(Path(df_lib_path).parent if df_lib_path else None),
        )
        optimizer = getattr(result[8], 'optimizer_diagnostics', {})
        pairwise_optimizer = getattr(
            result[8], 'pairwise_optimizer_diagnostics', {}
        )
        optimizer_success = (
            optimizer.get('success') is True
            and optimizer.get('nonunique_optimum') is not True
            and pairwise_optimizer.get('all_references_success') is True
        )
        optimizer_objective = optimizer.get('objective_log_likelihood')
        optimizer_finite = (
            optimizer_objective is not None
            and np.isfinite(float(optimizer_objective))
        )
        fold_optimizer = {
            'success': bool(optimizer_success and optimizer_finite),
            'reported_success': bool(optimizer_success),
            'finite_objective': bool(optimizer_finite),
            'message': str(optimizer.get('message', 'missing diagnostics')),
            'nfev': int(optimizer.get('nfev', -1)),
            'nit': int(optimizer.get('nit', -1)),
            'objective_log_likelihood': (
                float(optimizer_objective) if optimizer_finite else None
            ),
            'composite_nonunique_optimum': bool(
                optimizer.get('nonunique_optimum', False)
            ),
            'pairwise_all_references_success': bool(
                pairwise_optimizer.get('all_references_success', False)
            ),
            'pairwise_failed_reference_count': int(
                pairwise_optimizer.get('failed_reference_count', -1)
            ),
        }
        fold_support = getattr(result[8], 'support_diagnostics', {})
        fold_support_status = _support_diagnostics_status(fold_support)
        if not fold_support_status['valid']:
            raise RuntimeError(
                "Cross-validation fold produced invalid support diagnostics "
                f"({fold_support_status['reason']})."
            )
        fold_records[index] = {
            'index': index,
            'targetname': str(target.object),
            'target_identity': spectrum_identities[index],
            'values': [float(value) for value in result[:7]],
            'optimizer': fold_optimizer,
            'support': fold_support,
        }
        _atomic_json_write(checkpoint_path, {
            'key': checkpoint_key,
            'status': 'running',
            'folds': [fold_records[key] for key in sorted(fold_records)],
        })

    ordered = [fold_records[index] for index in range(len(df_lib))]
    failed_optimizer_folds = [
        item['index'] for item in ordered
        if item.get('optimizer', {}).get('success') is not True
    ]
    support_records = [item.get('support', {}) for item in ordered]
    support_limited_folds = int(sum(
        item.get('publication_support_validated') is False
        for item in support_records
    ))
    library_edge_folds = int(sum(
        item.get('any_full_library_edge') is True
        for item in support_records
    ))
    simplex_vertex_folds = int(sum(
        item.get('simplex', {}).get('is_vertex') is True
        for item in support_records
    ))
    frame = pd.DataFrame([item['values'] for item in ordered], columns=columns)
    frame['teff_true'] = frame.teff - frame.d_teff
    frame['feh_true'] = frame.feh - frame.d_feh
    frame['logg_true'] = frame.logg - frame.d_logg
    frame['targetname'] = [item['targetname'] for item in ordered]
    result_path = Path(outputdir) / f'crossvalidation_results_o{order}.csv'
    _atomic_dataframe_write(frame, result_path)
    result_csv_sha256 = sha256_file(result_path)
    _atomic_json_write(checkpoint_path, {
        'key': checkpoint_key,
        'status': 'complete',
        'folds': ordered,
        'result_csv': result_path.name,
        'result_csv_sha256': result_csv_sha256,
        'result_csv_hash_scope': 'file_bytes',
        'optimizer': {
            'all_folds_success': not failed_optimizer_folds,
            'failed_fold_indices': failed_optimizer_folds,
        },
    })
    checkpoint_sha256 = sha256_file(checkpoint_path)
    _log(verbose, 'Saved raw cross-validation results to {}'.format(result_path))

    product_paths = write_crossvalidation_products(
        frame, order=int(order), outputdir=outputdir, library_id=library_id,
        cv_source=os.fspath(result_path),
        fit_feh_calibration=fit_feh_calibration,
        population_rules=population_rules,
        min_calibration_stars=min_calibration_stars,
        source_pipeline_version=__version__,
        source_drp_compatibility_status=drp_compatibility['status'],
        source_library_manifest_sha256=library_manifest.get('manifest_sha256'),
        source_neidspec_version=getattr(neidspec, '__version__', 'unknown'),
        source_neidspec_source_sha256=_neidspec_source_fingerprint(),
        source_blaze_source=blaze_source,
        source_drp_versions=json.dumps(
            drp_compatibility['versions'], sort_keys=True
        ),
        source_pipeline_sha256=_pipeline_source_fingerprint(),
        source_cv_results_sha256=result_csv_sha256,
        source_cv_results_hash_scope='file_bytes',
        source_checkpoint=checkpoint_path.name,
        source_checkpoint_sha256=checkpoint_sha256,
        source_checkpoint_hash_scope='file_bytes',
        source_reference_dq_status_counts=json.dumps(
            _quality_counts(HLS.splist), sort_keys=True
        ),
        source_all_folds_optimizer_success=not failed_optimizer_folds,
        source_failed_optimizer_folds=len(failed_optimizer_folds),
        source_estimator_config=estimator_config,
        source_support_limited_folds=support_limited_folds,
        source_library_edge_folds=library_edge_folds,
        source_simplex_vertex_folds=simplex_vertex_folds,
        source_fold_support=support_records,
    )
    frame.attrs['products'] = product_paths
    frame.attrs['drp_compatibility'] = drp_compatibility
    frame.attrs['checkpoint'] = os.fspath(checkpoint_path)
    if plot_results:
        fig, ax = plt.subplots(dpi=200)
        plot_crossval_feh_delta_feh(
            frame['feh_true'].values, frame['d_feh'].values, ax=ax
        )
        fig.savefig(
            Path(outputdir) / f'crossvalidation_feh_delta_feh_o{order}.png',
            dpi=200,
        )
        plt.close(fig)
        plot_crossvalidation_results_1d(order, frame, outputdir)
        plot_crossvalidation_results_2d(order, frame, outputdir)
        plot_crossvalidation_results_main(order, frame, outputdir)
    return frame


def write_crossvalidation_products(
        df_crossval, *, order, outputdir, library_id, cv_source,
        fit_feh_calibration=False, population_rules=None,
        min_calibration_stars=8, source_pipeline_version=None,
        source_drp_compatibility_status='unknown',
        source_library_manifest_sha256=None,
        source_neidspec_version=None, source_blaze_source='unknown',
        source_drp_versions='{}', source_pipeline_sha256=None,
        source_neidspec_source_sha256=None,
        source_cv_results_sha256=None,
        source_cv_results_hash_scope=None,
        source_checkpoint=None, source_checkpoint_sha256=None,
        source_checkpoint_hash_scope=None,
        source_reference_dq_status_counts='{}',
        source_all_folds_optimizer_success=False,
        source_failed_optimizer_folds=None, source_estimator_config=None,
        source_support_limited_folds=None,
        source_library_edge_folds=None,
        source_simplex_vertex_folds=None, source_fold_support=None):
    """Write machine-readable per-order uncertainty/calibration products."""
    outputdir = os.fspath(outputdir)
    os.makedirs(outputdir, exist_ok=True)
    populations = tuple(
        rule if isinstance(rule, PopulationRule) else PopulationRule.from_mapping(rule)
        for rule in (population_rules or (PopulationRule(),))
    )
    summary_path = os.path.join(
        outputdir, 'crossvalidation_summary_o{}.csv'.format(order))
    paths = {'summary': summary_path}
    calibration_artifact_id = None
    calibration_artifact_sha256 = None
    try:
        source_drp_counts = json.loads(source_drp_versions)
    except (TypeError, json.JSONDecodeError):
        source_drp_counts = {}
    source_drp_version_set = tuple(sorted(
        str(value) for value in source_drp_counts
    ))
    if source_estimator_config is None:
        estimator_config_json = None
        estimator_config_sha256 = None
        estimator_config_schema_version = np.nan
    else:
        estimator_config_json = _canonical_estimator_config_json(
            source_estimator_config
        )
        estimator_config_sha256 = _hash_estimator_config(
            source_estimator_config
        )
        estimator_config_schema_version = source_estimator_config.get(
            'schema_version'
        )
    if fit_feh_calibration:
        if library_id == 'unspecified-library':
            raise CalibrationError(
                "A non-placeholder library_id is required to publish calibration metadata."
            )
        artifact, summary, validation = build_calibration_products(
            df_crossval,
            order=int(order),
            library_id=library_id,
            cv_source=Path(cv_source).name,
            populations=populations,
            min_stars=min_calibration_stars,
            cv_source_sha256=(
                sha256_file(cv_source) if Path(cv_source).is_file() else None
            ),
            cv_source_hash_scope=(
                'file_bytes' if Path(cv_source).is_file()
                else 'canonical_dataframe_csv'
            ),
            source_library_manifest_sha256=source_library_manifest_sha256,
            source_pipeline_sha256=(
                source_pipeline_sha256 or _pipeline_source_fingerprint()
            ),
            source_all_folds_optimizer_success=(
                source_all_folds_optimizer_success
            ),
            source_failed_optimizer_folds=source_failed_optimizer_folds,
            source_neidspec_version=source_neidspec_version,
            source_neidspec_sha256=(
                source_neidspec_source_sha256 or _neidspec_source_fingerprint()
            ),
            source_result_schema_version=RESULT_SCHEMA_VERSION,
            source_blaze_source=source_blaze_source,
            source_drp_version_set=source_drp_version_set,
        )
        artifact_path = os.path.join(
            outputdir, 'feh_calibration_o{}.json'.format(order))
        validation_path = os.path.join(
            outputdir, 'crossvalidation_calibrated_o{}.csv'.format(order))
        artifact.save(artifact_path)
        calibration_artifact_id = artifact.artifact_id
        calibration_artifact_sha256 = artifact.canonical_sha256
        _atomic_dataframe_write(validation, validation_path)
        paths.update(artifact=artifact_path, calibrated_validation=validation_path)
    else:
        summary = pd.DataFrame([
            raw_crossvalidation_summary(
                df_crossval,
                order=int(order),
                population=population,
                library_id=library_id,
                cv_source=Path(cv_source).name,
            )
            for population in populations
        ])
    summary['source_pipeline_version'] = (
        source_pipeline_version if source_pipeline_version is not None
        else 'unknown_archived_pipeline'
    )
    summary['source_result_schema_version'] = RESULT_SCHEMA_VERSION
    summary['source_pipeline_sha256'] = (
        source_pipeline_sha256 or _pipeline_source_fingerprint()
    )
    summary['source_neidspec_source_sha256'] = (
        source_neidspec_source_sha256 or _neidspec_source_fingerprint()
    )
    summary['source_cv_results_sha256'] = (
        source_cv_results_sha256
        or (sha256_file(cv_source) if Path(cv_source).is_file() else None)
    )
    summary['source_cv_results_hash_scope'] = (
        source_cv_results_hash_scope
        or ('file_bytes' if Path(cv_source).is_file() else 'unverified')
    )
    summary['source_checkpoint'] = source_checkpoint
    summary['source_checkpoint_sha256'] = source_checkpoint_sha256
    summary['source_checkpoint_hash_scope'] = source_checkpoint_hash_scope
    summary['source_reference_dq_status_counts'] = (
        source_reference_dq_status_counts
    )
    summary['source_drp_compatibility_status'] = (
        source_drp_compatibility_status
    )
    summary['source_library_manifest_sha256'] = (
        source_library_manifest_sha256
    )
    summary['source_neidspec_version'] = (
        source_neidspec_version or 'unknown'
    )
    summary['source_blaze_source'] = source_blaze_source
    summary['source_drp_versions'] = source_drp_versions
    summary['source_drp_version_set'] = json.dumps(list(source_drp_version_set))
    summary['source_all_folds_optimizer_success'] = bool(
        source_all_folds_optimizer_success
    )
    summary['source_failed_optimizer_folds'] = (
        int(source_failed_optimizer_folds)
        if source_failed_optimizer_folds is not None else np.nan
    )
    summary['source_estimator_config_schema_version'] = (
        estimator_config_schema_version
    )
    summary['source_estimator_config_sha256'] = estimator_config_sha256
    summary['source_estimator_config_json'] = estimator_config_json
    if source_fold_support is not None:
        support_records = list(source_fold_support)
        if len(support_records) != len(df_crossval):
            raise ValueError(
                "source_fold_support must contain one record per CV fold."
            )
        invalid_support = [
            (index, status['reason'])
            for index, item in enumerate(support_records)
            for status in [_support_diagnostics_status(item)]
            if not status['valid']
        ]
        if invalid_support:
            raise ValueError(
                "source_fold_support contains invalid diagnostics: "
                f"{invalid_support[:5]}"
            )
        support_limited = np.asarray([
            item.get('publication_support_validated') is False
            for item in support_records
        ], dtype=bool)
        library_edge = np.asarray([
            item.get('any_full_library_edge') is True
            for item in support_records
        ], dtype=bool)
        simplex_vertex = np.asarray([
            item.get('simplex', {}).get('is_vertex') is True
            for item in support_records
        ], dtype=bool)
        for row_index, population in enumerate(populations):
            population_mask = population.mask(df_crossval).to_numpy(
                dtype=bool
            )
            summary.loc[
                row_index, 'source_support_limited_folds'
            ] = int(np.count_nonzero(support_limited & population_mask))
            summary.loc[
                row_index, 'source_library_edge_folds'
            ] = int(np.count_nonzero(library_edge & population_mask))
            summary.loc[
                row_index, 'source_simplex_vertex_folds'
            ] = int(np.count_nonzero(simplex_vertex & population_mask))
    else:
        summary['source_support_limited_folds'] = (
            int(source_support_limited_folds)
            if source_support_limited_folds is not None else np.nan
        )
        summary['source_library_edge_folds'] = (
            int(source_library_edge_folds)
            if source_library_edge_folds is not None else np.nan
        )
        summary['source_simplex_vertex_folds'] = (
            int(source_simplex_vertex_folds)
            if source_simplex_vertex_folds is not None else np.nan
        )
    summary['error_metric_schema_version'] = (
        CROSSVALIDATION_METRIC_SCHEMA_VERSION
    )
    summary['residual_bias_statistic'] = 'arithmetic_mean'
    summary['residual_scatter_statistic'] = 'sample_standard_deviation'
    summary['residual_scatter_ddof'] = 1
    summary['predictive_error_statistic'] = 'root_mean_square_error'
    summary['predictive_error_ddof'] = 0
    summary['error_residual_method'] = CROSSVALIDATION_RESIDUAL_METHOD
    summary['uncertainty_semantics'] = NO_AUTOMATIC_UNCERTAINTY_SEMANTICS
    # Retain the original columns so older CSV readers fail explicitly rather
    # than silently treating RMSE as one sigma.  The quantitative legacy
    # aliases (``*_bias``, ``*_sigma``, ``*_rmse``) remain readable.
    summary['uncertainty_statistic'] = 'not_automatically_defined'
    summary['uncertainty_ddof'] = np.nan
    summary['uncertainty_residual_method'] = CROSSVALIDATION_RESIDUAL_METHOD
    summary['feh_calibration_artifact_id'] = calibration_artifact_id
    summary['feh_calibration_artifact_sha256'] = calibration_artifact_sha256
    _atomic_dataframe_write(summary, summary_path)
    return paths


def detrend_feh(feh, *, order=None, calibration=None, teff=None,
                population=None, library_id=None):
    """
    Calibrate trend in (Delta Fe/H = Fe/H_recovered - Fe/H_true) vs Fe/H_true
    
    INPUT:
        recovered feh value from specmatch
        
    OUTPUT:
        calibrated feh after removing trend
        
    This compatibility name no longer contains an HPF/order-5 coefficient.
    A validated NEID calibration artifact and order are mandatory.
    """
    if calibration is None or order is None:
        raise UnvalidatedCalibrationError(
            "detrend_feh requires order= and a validated calibration= artifact; "
            "no embedded coefficients are available."
        )
    calibration_set = coerce_calibration_set(calibration)
    selected = calibration_set.select(
        int(order), teff=teff, population=population, library_id=library_id)
    return selected.apply(feh)


def plot_crossvalidation_results_main(order, df_crossval, savefolder):
    """
    Main crossvalidation results plot
    """
    PW = 10.
    PH = 3.
    fig = plt.figure(figsize=(PW, PH), dpi=200)
    gs0 = GridSpec(1, 2)
    gs0.update(top=0.92, bottom=0.11, wspace=0.01, left=0.05, right=0.6)

    gs1 = GridSpec(1, 1)
    gs1.update(top=0.92, bottom=0.11, hspace=0.18, left=0.68, right=0.98)

    ax = plt.subplot(gs0[0, 0])
    bx = plt.subplot(gs0[0, 1])
    cx = plt.subplot(gs1[0, 0])

    ax.plot(df_crossval['feh_true'], df_crossval['teff_true'], marker='o', lw=0, markeredgecolor='black')
    bx.plot(df_crossval['logg_true'], df_crossval['teff_true'], marker='o', lw=0, markeredgecolor='black')
    for i in range(len(df_crossval)):
        _x = [df_crossval['feh_true'].values[i], df_crossval['feh'].values[i]]
        _y = [df_crossval['teff_true'].values[i], df_crossval['teff'].values[i]]
        ax.plot(_x, _y, color='crimson', lw=0.5, zorder=-10)

        _x = [df_crossval['logg_true'].values[i], df_crossval['logg'].values[i]]
        if i == 0:
            bx.plot(_x, _y, color='crimson', lw=0.5, label='Crossvalidation Value', zorder=-10)
        else:
            bx.plot(_x, _y, color='crimson', lw=0.5, zorder=-10)

    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel('[Fe/H]', labelpad=0)
    ax.set_ylabel(r'$T_{\mathrm{eff}}$ [K]', labelpad=0)
    bx.set_xlabel(r'$\log (g)$', labelpad=0)
    ax.minorticks_on()
    bx.minorticks_on()
    bx.yaxis.set_ticks([])
    bx.set_ylim(*ax.get_ylim())
    fig.subplots_adjust(wspace=0.05)
    bx.legend(loc='upper right')
    fig.suptitle('Library Performance (NEID Order {})'.format(order), y=0.98)
    plot_crossval_feh_delta_feh(df_crossval.feh_true.values, df_crossval.d_feh.values, ax=cx)
    fig.savefig('{}/crossvalidation_o{}_main.png'.format(savefolder, order))


def plot_crossval_feh_delta_feh(feh_true, d_feh, ax=None):
    """
    INPUT:
        feh_true - True Fe/H
        d_feh - delta Fe/H = Fe/H_specmatch - Fe/H_True

    EXAMPLE:
        plot_crossval_feh_delta_feh(df_res['[Fe/H]'].values,df_res['d_feh'].values)
    """
    if ax is None:
        fig, ax = plt.subplots(dpi=200)
    ax.plot(feh_true, d_feh, marker='o', lw=0, color='black')

    p = np.polyfit(feh_true, d_feh, deg=1)
    xx = np.linspace(-0.5, 0.5, 200)
    yy = np.polyval(p, xx)
    ax.plot(xx, yy, color='crimson', label='Linear fit\n$p_1$={:0.4f}\n$p_2$={:0.4f}'.format(p[0], p[1]))
    ax.legend(loc='upper right')
    ax.set_xlabel('[Fe/H]', labelpad=0)
    ax.set_ylabel(r'$\Delta$[Fe/H] = $\mathrm{[Fe/H]}_{\mathrm{Recovered}}$ - $\mathrm{[Fe/H]}_{\mathrm{True}}$')
    ax.minorticks_on()
