"""Validated metallicity calibration for NEIDSpecMatch.

The spectral matcher returns an empirical-library weighted metallicity.  A
calibration is an *optional* second stage.  No coefficients live in source
code: calibration artifacts are generated from a named cross-validation data
set and carry enough metadata to reject order, population, or library
mismatches.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


CALIBRATION_SCHEMA_VERSION = 2
RESIDUAL_DEFINITION = "feh_recovered_minus_feh_true"
VALIDATION_METHOD = "spectral_leave_one_out_plus_calibration_leave_one_out"
POPULATION_SELECTION_BASIS = "teff_true_explicit_population_required_at_application"
VALIDATION_LIMITATION = (
    "calibration folds use correlated spectral-LOO predictions; training-fold "
    "spectral fits may include the calibration-held-out star as a reference"
)
STRICT_NESTED_VALIDATION = False
MIN_ABS_CALIBRATION_DENOMINATOR = 0.1
PERFORMANCE_GATE = (
    "calibrated_crossfit_sigma_and_rmse_improve_and_abs_bias_does_not_worsen"
)
CROSSVALIDATION_METRIC_SCHEMA_VERSION = 2
CROSSVALIDATION_RESIDUAL_METHOD = (
    "spectral_leave_one_out_recovered_minus_reference_truth"
)
NO_AUTOMATIC_UNCERTAINTY_SEMANTICS = (
    "crossvalidation_reports_residual_bias_sample_scatter_and_predictive_rmse_"
    "separately;none_is_automatically_a_gaussian_1sigma_uncertainty"
)


class CalibrationError(ValueError):
    """Base exception for invalid or incompatible calibration artifacts."""


class UnvalidatedCalibrationError(CalibrationError):
    """Raised when an unvalidated or unspecified calibration is requested."""


@dataclass(frozen=True)
class PopulationRule:
    """Temperature-domain definition for one calibration population.

    Bounds are inclusive.  ``None`` means that the population is unbounded on
    that side.  The only built-in population is ``all``; science-specific
    temperature cuts must be supplied explicitly and are recorded verbatim.
    """

    name: str = "all"
    teff_min: Optional[float] = None
    teff_max: Optional[float] = None

    def __post_init__(self) -> None:
        if not self.name or not str(self.name).strip():
            raise CalibrationError("A calibration population needs a non-empty name.")
        if (
            self.teff_min is not None
            and self.teff_max is not None
            and self.teff_min > self.teff_max
        ):
            raise CalibrationError("teff_min cannot exceed teff_max.")
        for name, value in (("teff_min", self.teff_min), ("teff_max", self.teff_max)):
            if value is not None and not np.isfinite(value):
                raise CalibrationError(f"{name} must be finite when supplied.")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "PopulationRule":
        return cls(
            name=str(value.get("name", "all")),
            teff_min=_optional_float(value.get("teff_min")),
            teff_max=_optional_float(value.get("teff_max")),
        )

    def contains(self, teff: Optional[float]) -> bool:
        if self.teff_min is None and self.teff_max is None:
            return True
        if teff is None or not np.isfinite(teff):
            return False
        return (
            (self.teff_min is None or teff >= self.teff_min)
            and (self.teff_max is None or teff <= self.teff_max)
        )

    def mask(self, frame: pd.DataFrame) -> pd.Series:
        if "teff_true" not in frame:
            if self.teff_min is None and self.teff_max is None:
                return pd.Series(True, index=frame.index)
            raise CalibrationError(
                "Temperature-bounded populations require a 'teff_true' column."
            )
        values = pd.to_numeric(frame["teff_true"], errors="coerce")
        mask = pd.Series(True, index=frame.index)
        if self.teff_min is not None:
            mask &= values >= self.teff_min
        if self.teff_max is not None:
            mask &= values <= self.teff_max
        return mask

    @property
    def width(self) -> float:
        low = -np.inf if self.teff_min is None else self.teff_min
        high = np.inf if self.teff_max is None else self.teff_max
        return float(high - low)


@dataclass(frozen=True)
class FeHCalibration:
    """One order-specific linear [Fe/H] calibration candidate."""

    order: int
    population: PopulationRule
    slope: float
    intercept: float
    library_id: str
    cv_source: str
    n_stars: int
    raw_sigma: float
    calibrated_sigma: float
    raw_bias: float
    calibrated_bias: float
    raw_rmse: float
    calibrated_rmse: float
    min_abs_denominator: float
    performance_gate: str
    performance_passed: bool
    validation_reason: str
    validated: bool
    validation_method: str = VALIDATION_METHOD
    validation_limitation: str = VALIDATION_LIMITATION
    population_selection_basis: str = POPULATION_SELECTION_BASIS
    residual_definition: str = RESIDUAL_DEFINITION
    created_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    def __post_init__(self) -> None:
        if int(self.order) != self.order:
            raise CalibrationError("Calibration order must be an integer.")
        if not self.library_id.strip():
            raise CalibrationError("Calibration metadata requires a library_id.")
        if not self.cv_source.strip():
            raise CalibrationError("Calibration metadata requires a cv_source.")
        if self.residual_definition != RESIDUAL_DEFINITION:
            raise CalibrationError(
                f"Unsupported residual definition: {self.residual_definition!r}."
            )
        numeric_metrics = {
            "slope": self.slope,
            "intercept": self.intercept,
            "raw_sigma": self.raw_sigma,
            "calibrated_sigma": self.calibrated_sigma,
            "raw_bias": self.raw_bias,
            "calibrated_bias": self.calibrated_bias,
            "raw_rmse": self.raw_rmse,
            "calibrated_rmse": self.calibrated_rmse,
            "min_abs_denominator": self.min_abs_denominator,
        }
        invalid_metrics = [
            name for name, value in numeric_metrics.items()
            if not np.isfinite(value)
        ]
        if invalid_metrics:
            raise CalibrationError(
                "Calibration metrics must be finite: " + ", ".join(invalid_metrics)
            )
        if self.n_stars < 8:
            raise CalibrationError("A calibration artifact requires at least 8 stars.")
        if self.raw_sigma <= 0 or self.calibrated_sigma < 0:
            raise CalibrationError("Calibration scatter values are invalid.")
        if self.raw_rmse < 0 or self.calibrated_rmse < 0:
            raise CalibrationError("Calibration RMSE values cannot be negative.")
        if self.min_abs_denominator < 0:
            raise CalibrationError("Calibration conditioning cannot be negative.")
        denominator = 1.0 + self.slope
        if denominator <= 0.0:
            raise CalibrationError(
                "Calibration inverse must remain positively monotonic."
            )
        if self.min_abs_denominator > denominator + 1e-12:
            raise CalibrationError(
                "Calibration conditioning metadata exceeds the fitted denominator."
            )
        if self.performance_gate != PERFORMANCE_GATE:
            raise CalibrationError(
                f"Unsupported calibration performance gate: {self.performance_gate!r}."
            )
        performance_passed = (
            self.calibrated_sigma < self.raw_sigma
            and self.calibrated_rmse < self.raw_rmse
            and abs(self.calibrated_bias) <= abs(self.raw_bias) + 1e-12
        )
        if bool(self.performance_passed) != performance_passed:
            raise CalibrationError(
                "Calibration performance flag disagrees with cross-fit scatter."
            )
        conditioning_passed = (
            self.min_abs_denominator >= MIN_ABS_CALIBRATION_DENOMINATOR
        )
        # Calibration leave-one-out is performed on correlated spectral-LOO
        # predictions.  It is useful diagnostic cross-fitting, but it is not a
        # strict nested/leave-two-out validation of the calibrated estimator.
        expected_validated = (
            performance_passed and conditioning_passed
            and STRICT_NESTED_VALIDATION
        )
        if bool(self.validated) != expected_validated:
            raise CalibrationError(
                "Calibration validated flag disagrees with the declared gates."
            )
        if self.validated and self.validation_reason != "passed":
            raise CalibrationError("A validated calibration must record reason='passed'.")
        if not self.validated and self.validation_reason == "passed":
            raise CalibrationError("An unvalidated calibration cannot record reason='passed'.")
        if self.validation_method != VALIDATION_METHOD:
            raise CalibrationError(
                f"Unsupported calibration validation method: {self.validation_method!r}."
            )
        if self.population_selection_basis != POPULATION_SELECTION_BASIS:
            raise CalibrationError("Unsupported calibration population-selection basis.")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "FeHCalibration":
        population = value.get("population", {"name": "all"})
        if not isinstance(population, PopulationRule):
            population = PopulationRule.from_mapping(population)
        return cls(
            order=int(value["order"]),
            population=population,
            slope=float(value["slope"]),
            intercept=float(value["intercept"]),
            library_id=str(value["library_id"]),
            cv_source=str(value["cv_source"]),
            n_stars=int(value["n_stars"]),
            raw_sigma=float(value["raw_sigma"]),
            calibrated_sigma=float(value["calibrated_sigma"]),
            raw_bias=float(value["raw_bias"]),
            calibrated_bias=float(value["calibrated_bias"]),
            raw_rmse=float(value.get("raw_rmse", np.nan)),
            calibrated_rmse=float(value["calibrated_rmse"]),
            min_abs_denominator=float(
                value.get("min_abs_denominator", abs(1.0 + float(value["slope"])))
            ),
            performance_gate=str(value.get("performance_gate", PERFORMANCE_GATE)),
            performance_passed=_strict_bool(
                value.get("performance_passed"), "performance_passed"
            ),
            validation_reason=str(value.get("validation_reason", "not recorded")),
            validated=_strict_bool(value.get("validated"), "validated"),
            validation_method=str(value.get("validation_method", VALIDATION_METHOD)),
            validation_limitation=str(
                value.get("validation_limitation", VALIDATION_LIMITATION)
            ),
            population_selection_basis=str(
                value.get("population_selection_basis", POPULATION_SELECTION_BASIS)
            ),
            residual_definition=str(value.get("residual_definition", RESIDUAL_DEFINITION)),
            created_utc=str(value.get("created_utc", "unknown")),
        )

    def apply(self, feh_recovered: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Invert the fitted recovered-minus-true trend."""
        if (
            not self.validated
            or not self.performance_passed
            or self.min_abs_denominator < MIN_ABS_CALIBRATION_DENOMINATOR
        ):
            raise UnvalidatedCalibrationError(
                "Refusing to apply an unvalidated or ill-conditioned calibration."
            )
        calibrated = (np.asarray(feh_recovered) - self.intercept) / (1.0 + self.slope)
        if np.ndim(feh_recovered) == 0:
            return float(calibrated)
        return calibrated

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FeHCalibrationSet:
    """Collection of explicitly versioned order/population calibrations."""

    calibrations: Tuple[FeHCalibration, ...]
    artifact_id: str
    cv_source_sha256: Optional[str] = None
    cv_source_hash_scope: Optional[str] = None
    source_library_manifest_sha256: Optional[str] = None
    source_pipeline_sha256: Optional[str] = None
    source_all_folds_optimizer_success: bool = False
    source_failed_optimizer_folds: Optional[int] = None
    source_neidspec_version: Optional[str] = None
    source_neidspec_sha256: Optional[str] = None
    source_result_schema_version: Optional[int] = None
    source_blaze_source: Optional[str] = None
    source_drp_version_set: Tuple[str, ...] = ()
    schema_version: int = CALIBRATION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != CALIBRATION_SCHEMA_VERSION:
            raise CalibrationError(
                f"Unsupported calibration schema {self.schema_version}; "
                f"expected {CALIBRATION_SCHEMA_VERSION}."
            )
        if not self.artifact_id.strip():
            raise CalibrationError("Calibration artifacts require an artifact_id.")
        keys = [(item.order, item.population.name) for item in self.calibrations]
        if len(keys) != len(set(keys)):
            raise CalibrationError("Duplicate order/population calibration entries.")
        for name in (
            "cv_source_sha256", "source_library_manifest_sha256",
            "source_pipeline_sha256", "source_neidspec_sha256",
        ):
            value = getattr(self, name)
            if value is not None and not re.fullmatch(r"[0-9a-f]{64}", str(value)):
                raise CalibrationError(f"{name} must be a lowercase SHA-256 digest.")
        if self.cv_source_hash_scope not in {
            None, "file_bytes", "canonical_dataframe_csv"
        }:
            raise CalibrationError("Unsupported calibration CV hash scope.")
        if (
            self.source_failed_optimizer_folds is not None
            and self.source_failed_optimizer_folds < 0
        ):
            raise CalibrationError("Failed optimizer fold count cannot be negative.")
        if self.source_result_schema_version is not None:
            if int(self.source_result_schema_version) != self.source_result_schema_version:
                raise CalibrationError("Source result schema must be an integer.")
        if any(not str(value).strip() for value in self.source_drp_version_set):
            raise CalibrationError("Source DRP versions cannot be blank.")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "FeHCalibrationSet":
        artifact = cls(
            calibrations=tuple(
                FeHCalibration.from_mapping(item)
                for item in value.get("calibrations", [])
            ),
            artifact_id=str(value["artifact_id"]),
            cv_source_sha256=value.get("cv_source_sha256"),
            cv_source_hash_scope=value.get("cv_source_hash_scope"),
            source_library_manifest_sha256=value.get(
                "source_library_manifest_sha256"
            ),
            source_pipeline_sha256=value.get("source_pipeline_sha256"),
            source_all_folds_optimizer_success=_strict_bool(
                value.get("source_all_folds_optimizer_success", False),
                "source_all_folds_optimizer_success",
            ),
            source_failed_optimizer_folds=(
                int(value["source_failed_optimizer_folds"])
                if value.get("source_failed_optimizer_folds") is not None else None
            ),
            source_neidspec_version=value.get("source_neidspec_version"),
            source_neidspec_sha256=value.get("source_neidspec_sha256"),
            source_result_schema_version=(
                int(value["source_result_schema_version"])
                if value.get("source_result_schema_version") is not None else None
            ),
            source_blaze_source=value.get("source_blaze_source"),
            source_drp_version_set=tuple(
                str(item) for item in value.get("source_drp_version_set", ())
            ),
            schema_version=int(value.get("schema_version", CALIBRATION_SCHEMA_VERSION)),
        )
        recorded_hash = value.get("artifact_sha256")
        if recorded_hash is None:
            raise CalibrationError(
                "Calibration schema 2 artifacts require artifact_sha256."
            )
        if str(recorded_hash) != artifact.canonical_sha256:
            raise CalibrationError("Calibration artifact SHA-256 does not match its content.")
        return artifact

    @classmethod
    def load(cls, filename: Union[str, Path]) -> "FeHCalibrationSet":
        with Path(filename).open(encoding="utf-8") as stream:
            return cls.from_mapping(json.load(stream))

    def save(self, filename: Union[str, Path]) -> Path:
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}-", suffix=".tmp", dir=path.parent
        )
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                json.dump(self.to_dict(), stream, indent=2, sort_keys=True)
                stream.write("\n")
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
        return path

    def to_dict(self) -> dict[str, Any]:
        payload = self._canonical_payload()
        return {
            **payload,
            "artifact_sha256": self.canonical_sha256,
        }

    def _canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "artifact_id": self.artifact_id,
            "cv_source_sha256": self.cv_source_sha256,
            "cv_source_hash_scope": self.cv_source_hash_scope,
            "source_library_manifest_sha256": self.source_library_manifest_sha256,
            "source_pipeline_sha256": self.source_pipeline_sha256,
            "source_all_folds_optimizer_success": (
                self.source_all_folds_optimizer_success
            ),
            "source_failed_optimizer_folds": self.source_failed_optimizer_folds,
            "source_neidspec_version": self.source_neidspec_version,
            "source_neidspec_sha256": self.source_neidspec_sha256,
            "source_result_schema_version": self.source_result_schema_version,
            "source_blaze_source": self.source_blaze_source,
            "source_drp_version_set": list(self.source_drp_version_set),
            "calibrations": [item.to_dict() for item in self.calibrations],
        }

    @property
    def canonical_sha256(self) -> str:
        encoded = json.dumps(
            self._canonical_payload(), sort_keys=True, separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    @property
    def source_evidence_validated(self) -> bool:
        digests = (
            self.cv_source_sha256,
            self.source_library_manifest_sha256,
            self.source_pipeline_sha256,
            self.source_neidspec_sha256,
        )
        return bool(
            all(re.fullmatch(r"[0-9a-f]{64}", str(value or "")) for value in digests)
            and self.cv_source_hash_scope in {
                "file_bytes", "canonical_dataframe_csv"
            }
            and bool(self.source_all_folds_optimizer_success)
            and self.source_failed_optimizer_folds == 0
            and self.source_neidspec_version not in {None, "", "unknown"}
            and self.source_result_schema_version == 3
            and self.source_blaze_source not in {None, "", "unknown"}
            and len(self.source_drp_version_set) == 1
            and self.source_drp_version_set[0].strip().lower() != "unknown"
        )

    @property
    def publication_validation_eligible(self) -> bool:
        """Current calibration cross-fitting is not strict nested validation."""
        return False

    def select(
        self,
        order: int,
        *,
        teff: Optional[float] = None,
        population: Optional[str] = None,
        library_id: Optional[str] = None,
    ) -> FeHCalibration:
        candidates = [item for item in self.calibrations if item.order == int(order)]
        if population is not None:
            candidates = [
                item for item in candidates if item.population.name == population
            ]
        else:
            # Populations were defined with reference ``teff_true`` during
            # validation.  Selecting one automatically with recovered Teff
            # would add an unvalidated classifier at science time.  Only an
            # unbounded population can therefore be selected implicitly.
            candidates = [
                item for item in candidates
                if item.population.teff_min is None and item.population.teff_max is None
            ]
        if library_id is not None:
            candidates = [item for item in candidates if item.library_id == library_id]
        if not candidates:
            context = f"order={order}, population={population!r}, library={library_id!r}"
            raise CalibrationError(f"No compatible validated [Fe/H] calibration for {context}.")
        if len(candidates) > 1:
            candidates.sort(key=lambda item: item.population.width)
            if np.isclose(candidates[0].population.width, candidates[1].population.width):
                raise CalibrationError(
                    "Multiple equally specific calibrations match; select a population explicitly."
                )
        selected = candidates[0]
        bounded = (
            selected.population.teff_min is not None
            or selected.population.teff_max is not None
        )
        if bounded and teff is None:
            raise CalibrationError(
                "A bounded calibration population requires an independently "
                "supplied membership Teff. Recovered SpecMatch Teff must not "
                "silently select its own calibration population."
            )
        if bounded and not selected.population.contains(teff):
            raise CalibrationError(
                f"Membership Teff {teff!r} is outside population "
                f"{selected.population.name!r}."
            )
        if not selected.validated:
            raise UnvalidatedCalibrationError(
                f"Calibration {self.artifact_id}/{selected.population.name}/o{order} "
                "is marked unvalidated."
            )
        if not self.publication_validation_eligible:
            raise UnvalidatedCalibrationError(
                "NEIDSpecMatch 0.2 does not apply calibrations from the current "
                "correlated spectral-LOO workflow. Strict nested leave-two-out "
                "or independent validation is required."
            )
        if not self.source_evidence_validated:
            raise UnvalidatedCalibrationError(
                f"Calibration {self.artifact_id!r} lacks complete, successful "
                "CV/manifest/pipeline evidence and cannot be applied."
            )
        return selected


CalibrationLike = Union[
    FeHCalibrationSet,
    FeHCalibration,
    Mapping[str, Any],
    str,
    Path,
]


def coerce_calibration_set(value: CalibrationLike) -> FeHCalibrationSet:
    if isinstance(value, FeHCalibrationSet):
        return value
    if isinstance(value, FeHCalibration):
        return FeHCalibrationSet(
            calibrations=(value,),
            artifact_id=f"single-o{value.order}-{value.population.name}",
        )
    if isinstance(value, (str, Path)):
        return FeHCalibrationSet.load(value)
    if isinstance(value, Mapping):
        if "calibrations" in value:
            return FeHCalibrationSet.from_mapping(value)
        return coerce_calibration_set(FeHCalibration.from_mapping(value))
    raise TypeError(f"Unsupported calibration value: {type(value).__name__}.")


def fit_crossvalidated_feh_calibration(
    df_crossval: pd.DataFrame,
    *,
    order: int,
    population: PopulationRule = PopulationRule(),
    library_id: str,
    cv_source: str,
    min_stars: int = 8,
) -> tuple[FeHCalibration, pd.DataFrame]:
    """Fit and cross-fit a linear metallicity calibration.

    The input must be the output of spectral leave-one-out validation.  For
    every validation star, calibration coefficients are fitted using all
    *other* spectral-LOO predictions in that population.  This tests applying
    the calibration without fitting its coefficients on that row, but it is
    not strict nested CV: other spectral-LOO predictions can have used the
    calibration-held-out star as a reference.  The limitation is recorded in
    the artifact.
    """
    required = {"feh", "feh_true"}
    missing = sorted(required.difference(df_crossval.columns))
    if missing:
        raise CalibrationError(f"Cross-validation data are missing columns: {missing}.")
    selected = df_crossval.loc[population.mask(df_crossval)].copy()
    selected["feh"] = pd.to_numeric(selected["feh"], errors="coerce")
    selected["feh_true"] = pd.to_numeric(selected["feh_true"], errors="coerce")
    selected = selected.loc[
        np.isfinite(selected["feh"]) & np.isfinite(selected["feh_true"])
    ].copy()
    n_stars = len(selected)
    if n_stars < min_stars:
        raise CalibrationError(
            f"Population {population.name!r} has {n_stars} valid stars; "
            f"at least {min_stars} are required."
        )
    true = selected["feh_true"].to_numpy(dtype=float)
    recovered = selected["feh"].to_numpy(dtype=float)
    raw_residual = recovered - true
    slope, intercept = np.polyfit(true, raw_residual, deg=1)
    if 1.0 + slope <= 0.0:
        raise CalibrationError(
            "Fitted calibration reverses metallicity rank and is unsupported."
        )

    calibrated = np.empty(n_stars, dtype=float)
    fold_denominators = []
    for index in range(n_stars):
        train = np.arange(n_stars) != index
        fold_slope, fold_intercept = np.polyfit(
            true[train], raw_residual[train], deg=1
        )
        fold_denominator = 1.0 + fold_slope
        fold_denominators.append(fold_denominator)
        if fold_denominator <= 0.0:
            raise CalibrationError(
                f"Calibration fold {index} reverses metallicity rank."
            )
        calibrated[index] = (
            recovered[index] - fold_intercept
        ) / fold_denominator
    calibrated_residual = calibrated - true
    min_abs_denominator = min(1.0 + slope, *fold_denominators)
    raw_sigma = _sample_std(raw_residual)
    calibrated_sigma = _sample_std(calibrated_residual)
    raw_bias = float(np.mean(raw_residual))
    calibrated_bias = float(np.mean(calibrated_residual))
    raw_rmse = float(np.sqrt(np.mean(raw_residual**2)))
    calibrated_rmse = float(np.sqrt(np.mean(calibrated_residual**2)))
    conditioning_passed = bool(
        min_abs_denominator >= MIN_ABS_CALIBRATION_DENOMINATOR)
    performance_passed = bool(
        calibrated_sigma < raw_sigma
        and calibrated_rmse < raw_rmse
        and abs(calibrated_bias) <= abs(raw_bias) + 1e-12
    )
    finite_passed = bool(np.all(np.isfinite(calibrated_residual)))
    validated = bool(
        conditioning_passed and performance_passed and finite_passed
        and STRICT_NESTED_VALIDATION
    )
    failed = []
    if not conditioning_passed:
        failed.append("conditioning")
    if not performance_passed:
        failed.append("crossfit scatter/RMSE/bias performance gate failed")
    if not finite_passed:
        failed.append("non-finite crossfit residuals")
    if not STRICT_NESTED_VALIDATION:
        failed.append("calibration cross-fitting is not strict nested validation")
    validation_reason = "passed" if not failed else "; ".join(failed)

    output = selected.copy()
    output["order"] = int(order)
    output["calibration_population"] = population.name
    output["feh_raw"] = recovered
    output["d_feh_raw"] = raw_residual
    output["feh_calibrated_crossfit"] = calibrated
    output["d_feh_calibrated_crossfit"] = calibrated_residual

    calibration = FeHCalibration(
        order=int(order),
        population=population,
        slope=float(slope),
        intercept=float(intercept),
        library_id=library_id,
        cv_source=cv_source,
        n_stars=n_stars,
        raw_sigma=raw_sigma,
        calibrated_sigma=calibrated_sigma,
        raw_bias=raw_bias,
        calibrated_bias=calibrated_bias,
        raw_rmse=raw_rmse,
        calibrated_rmse=calibrated_rmse,
        min_abs_denominator=float(min_abs_denominator),
        performance_gate=PERFORMANCE_GATE,
        performance_passed=performance_passed,
        validation_reason=validation_reason,
        validated=validated,
    )
    return calibration, output


def build_calibration_products(
    df_crossval: pd.DataFrame,
    *,
    order: int,
    library_id: str,
    cv_source: str,
    populations: Optional[Sequence[PopulationRule]] = None,
    artifact_id: Optional[str] = None,
    min_stars: int = 8,
    cv_source_sha256: Optional[str] = None,
    cv_source_hash_scope: Optional[str] = None,
    source_library_manifest_sha256: Optional[str] = None,
    source_pipeline_sha256: Optional[str] = None,
    source_all_folds_optimizer_success: bool = False,
    source_failed_optimizer_folds: Optional[int] = None,
    source_neidspec_version: Optional[str] = None,
    source_neidspec_sha256: Optional[str] = None,
    source_result_schema_version: Optional[int] = None,
    source_blaze_source: Optional[str] = None,
    source_drp_version_set: Sequence[str] = (),
) -> tuple[FeHCalibrationSet, pd.DataFrame, pd.DataFrame]:
    """Build an artifact, summary table, and long-form validation table."""
    populations = tuple(populations or (PopulationRule(),))
    calibrations = []
    validations = []
    summaries = []
    for population in populations:
        calibration, validation = fit_crossvalidated_feh_calibration(
            df_crossval,
            order=order,
            population=population,
            library_id=library_id,
            cv_source=cv_source,
            min_stars=min_stars,
        )
        calibrations.append(calibration)
        validations.append(validation)
        summaries.append(calibration_summary_row(calibration, validation))
    artifact = FeHCalibrationSet(
        calibrations=tuple(calibrations),
        artifact_id=artifact_id or f"{library_id}-crossval-o{int(order)}",
        cv_source_sha256=(
            cv_source_sha256 or _canonical_dataframe_sha256(df_crossval)
        ),
        cv_source_hash_scope=(
            cv_source_hash_scope or "canonical_dataframe_csv"
        ),
        source_library_manifest_sha256=source_library_manifest_sha256,
        source_pipeline_sha256=source_pipeline_sha256,
        source_all_folds_optimizer_success=source_all_folds_optimizer_success,
        source_failed_optimizer_folds=source_failed_optimizer_folds,
        source_neidspec_version=source_neidspec_version,
        source_neidspec_sha256=source_neidspec_sha256,
        source_result_schema_version=source_result_schema_version,
        source_blaze_source=source_blaze_source,
        source_drp_version_set=tuple(sorted(str(value) for value in source_drp_version_set)),
    )
    return (
        artifact,
        pd.DataFrame(summaries),
        pd.concat(validations, ignore_index=True),
    )


def raw_crossvalidation_summary(
    df_crossval: pd.DataFrame,
    *,
    order: int,
    population: PopulationRule = PopulationRule(),
    library_id: str,
    cv_source: str,
) -> dict[str, Any]:
    selected = df_crossval.loc[population.mask(df_crossval)].copy()
    teff_metrics = _residual_metrics(selected, "d_teff")
    feh_metrics = _residual_metrics(selected, "d_feh")
    calibrated_feh_metrics = _residual_metrics(
        selected, "d_feh_calibrated_crossfit"
    )
    logg_metrics = _residual_metrics(selected, "d_logg")
    return {
        "order": int(order),
        "population": population.name,
        "teff_min": population.teff_min,
        "teff_max": population.teff_max,
        "library_id": library_id,
        "cv_source": cv_source,
        "n_stars": len(selected),
        "calibration_status": "not_fitted",
        "calibration_validation_method": "none",
        **_summary_metric_fields("teff", teff_metrics),
        **_summary_metric_fields("feh_raw", feh_metrics),
        **_summary_metric_fields(
            "feh_calibrated", calibrated_feh_metrics
        ),
        **_summary_metric_fields("logg", logg_metrics),
        "feh_calibration_slope": np.nan,
        "feh_calibration_intercept": np.nan,
    }


def calibration_summary_row(
    calibration: FeHCalibration, validation: pd.DataFrame
) -> dict[str, Any]:
    teff_metrics = _residual_metrics(validation, "d_teff")
    raw_feh_metrics = _residual_metrics(validation, "d_feh_raw")
    calibrated_feh_metrics = _residual_metrics(
        validation, "d_feh_calibrated_crossfit"
    )
    logg_metrics = _residual_metrics(validation, "d_logg")
    return {
        "order": calibration.order,
        "population": calibration.population.name,
        "teff_min": calibration.population.teff_min,
        "teff_max": calibration.population.teff_max,
        "library_id": calibration.library_id,
        "cv_source": calibration.cv_source,
        "n_stars": calibration.n_stars,
        "calibration_status": (
            "statistical_gates_passed_exploratory_non_nested"
            if calibration.validated else "unvalidated"
        ),
        "calibration_validation_reason": calibration.validation_reason,
        "calibration_validation_method": calibration.validation_method,
        "calibration_validation_limitation": calibration.validation_limitation,
        "population_selection_basis": calibration.population_selection_basis,
        **_summary_metric_fields("teff", teff_metrics),
        **_summary_metric_fields("feh_raw", raw_feh_metrics),
        **_summary_metric_fields(
            "feh_calibrated", calibrated_feh_metrics
        ),
        **_summary_metric_fields("logg", logg_metrics),
        "feh_calibration_performance_gate": calibration.performance_gate,
        "feh_calibration_performance_passed": calibration.performance_passed,
        "feh_calibration_min_abs_denominator": calibration.min_abs_denominator,
        "feh_calibration_conditioning_threshold": MIN_ABS_CALIBRATION_DENOMINATOR,
        "feh_calibration_slope": calibration.slope,
        "feh_calibration_intercept": calibration.intercept,
        "residual_definition": calibration.residual_definition,
    }


def _sample_std(values: Iterable[float]) -> float:
    values = np.asarray(tuple(values), dtype=float)
    values = values[np.isfinite(values)]
    return float(np.std(values, ddof=1)) if len(values) > 1 else np.nan


def _residual_metrics(frame: pd.DataFrame, column: str) -> dict[str, Any]:
    if column not in frame:
        values = np.asarray([], dtype=float)
    else:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
        values = values[np.isfinite(values)]
    total = len(frame)
    return {
        "n_valid": int(values.size),
        "coverage": float(values.size / total) if total else 0.0,
        "bias": float(np.mean(values)) if values.size else np.nan,
        "sample_scatter": _sample_std(values),
        "predictive_rmse": (
            float(np.sqrt(np.mean(values ** 2))) if values.size else np.nan
        ),
    }


def _summary_metric_fields(
    prefix: str, metrics: Mapping[str, Any]
) -> dict[str, Any]:
    """Return explicit residual metrics plus backward-readable aliases.

    ``*_sigma``, ``*_bias``, and ``*_rmse`` existed in the first 0.2
    validation products.  They remain byte-level CSV fields for readers of
    those products, but the canonical names state the estimand and never imply
    that predictive RMSE is a Gaussian one-sigma uncertainty.
    """
    bias = metrics["bias"]
    sample_scatter = metrics["sample_scatter"]
    predictive_rmse = metrics["predictive_rmse"]
    return {
        f"{prefix}_residual_bias": bias,
        f"{prefix}_sample_scatter": sample_scatter,
        f"{prefix}_predictive_rmse": predictive_rmse,
        f"{prefix}_n_valid": metrics["n_valid"],
        f"{prefix}_finite_coverage": metrics["coverage"],
        # Backward-readable aliases.  Semantics are declared by the summary's
        # metric metadata, not by an automatic one-sigma interpretation.
        f"{prefix}_bias": bias,
        f"{prefix}_sigma": sample_scatter,
        f"{prefix}_rmse": predictive_rmse,
    }


def _optional_float(value: Any) -> Optional[float]:
    return None if value is None else float(value)


def _strict_bool(value: Any, name: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    raise CalibrationError(f"Calibration field {name!r} must be a JSON boolean.")


def _canonical_dataframe_sha256(frame: pd.DataFrame) -> str:
    canonical = frame.to_csv(index=False, lineterminator="\n", float_format="%.17g")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
