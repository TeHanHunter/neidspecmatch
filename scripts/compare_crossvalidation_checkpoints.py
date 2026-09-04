#!/usr/bin/env python3
"""Compare completed validation folds with a later continuity checkpoint.

This utility does not create or replace cross-validation.  It records the
numerical change for targets that were rerun after a release-candidate source
or dependency update, while preserving both checkpoints' original provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PARAMETERS = ("teff", "feh_raw", "logg", "vsini")
VALUE_INDEX = {"teff": 0, "feh_raw": 1, "logg": 2, "vsini": 3}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_checkpoint(path: Path, *, require_complete: bool) -> dict:
    path = path.resolve()
    value = json.loads(path.read_text(encoding="utf-8"))
    required = {"key", "status", "folds"}
    if not isinstance(value, dict) or not required.issubset(value):
        raise ValueError(f"Malformed checkpoint: {path}")
    if require_complete and value["status"] != "complete":
        raise ValueError(f"Baseline checkpoint is not complete: {path}")
    if not isinstance(value["folds"], list) or not value["folds"]:
        raise ValueError(f"Checkpoint has no folds: {path}")
    value["_path"] = path
    value["_sha256"] = sha256(path)
    return value


def fold_map(checkpoint: dict) -> dict[str, dict]:
    output = {}
    for fold in checkpoint["folds"]:
        name = str(fold.get("targetname", "")).strip()
        values = fold.get("values")
        if not name or name in output or not isinstance(values, list) or len(values) < 7:
            raise ValueError("Checkpoint has malformed or duplicate fold records")
        output[name] = fold
    return output


def parse_pair(value: str) -> tuple[int, Path, Path]:
    try:
        order_text, baseline, candidate = value.split("=", 2)
        order = int(order_text)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            "comparison must be ORDER=BASELINE_CHECKPOINT=CANDIDATE_CHECKPOINT"
        ) from exc
    return order, Path(baseline), Path(candidate)


def compare(order: int, baseline_path: Path, candidate_path: Path):
    baseline = load_checkpoint(baseline_path, require_complete=True)
    candidate = load_checkpoint(candidate_path, require_complete=False)
    for label, checkpoint in (("baseline", baseline), ("candidate", candidate)):
        observed = int(checkpoint["key"].get("order", -1))
        if observed != order:
            raise ValueError(f"{label} checkpoint order {observed} != {order}")

    baseline_folds = fold_map(baseline)
    candidate_folds = fold_map(candidate)
    shared = [name for name in baseline_folds if name in candidate_folds]
    if not shared:
        raise ValueError(f"No shared targets for order {order}")

    rows = []
    for name in shared:
        old = np.asarray(baseline_folds[name]["values"], dtype=float)
        new = np.asarray(candidate_folds[name]["values"], dtype=float)
        teff_true = float(old[0] - old[4])
        row = {
            "order": order,
            "targetname": name,
            "population": "cool" if teff_true < 4500.0 else "hot",
            "teff_true": teff_true,
        }
        for parameter in PARAMETERS:
            index = VALUE_INDEX[parameter]
            row[f"baseline_{parameter}"] = float(old[index])
            row[f"candidate_{parameter}"] = float(new[index])
            row[f"delta_{parameter}"] = float(new[index] - old[index])
            row[f"abs_delta_{parameter}"] = abs(row[f"delta_{parameter}"])
        rows.append(row)

    frame = pd.DataFrame(rows)
    summaries = []
    for population in ("all", "cool", "hot"):
        selected = frame if population == "all" else frame.loc[
            frame["population"] == population
        ]
        if selected.empty:
            continue
        row = {"order": order, "population": population, "n_compared": len(selected)}
        for parameter in PARAMETERS:
            values = selected[f"abs_delta_{parameter}"].to_numpy(float)
            row[f"{parameter}_median_abs_delta"] = float(np.median(values))
            row[f"{parameter}_max_abs_delta"] = float(np.max(values))
        summaries.append(row)

    provenance = {
        "order": order,
        "baseline": {
            "path": baseline_path.name,
            "sha256": baseline["_sha256"],
            "status": baseline["status"],
            "fold_count": len(baseline["folds"]),
            "key": baseline["key"],
        },
        "candidate": {
            "path": candidate_path.name,
            "sha256": candidate["_sha256"],
            "status": candidate["status"],
            "fold_count": len(candidate["folds"]),
            "key": candidate["key"],
        },
        "shared_fold_count": len(shared),
    }
    return frame, pd.DataFrame(summaries), provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison", action="append", type=parse_pair, required=True,
        help="ORDER=BASELINE_CHECKPOINT=CANDIDATE_CHECKPOINT; repeat per order",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    detail_frames = []
    summary_frames = []
    provenance = {
        "schema_version": 1,
        "purpose": "post-validation numerical continuity check; not replacement CV",
        "population_boundary_k": 4500.0,
        "comparisons": [],
    }
    for order, baseline, candidate in args.comparison:
        detail, summary, receipt = compare(order, baseline, candidate)
        detail_frames.append(detail)
        summary_frames.append(summary)
        provenance["comparisons"].append(receipt)

    detail = pd.concat(detail_frames, ignore_index=True)
    summary = pd.concat(summary_frames, ignore_index=True)
    detail_path = args.output_dir / "crossvalidation_release_continuity_rows.csv"
    summary_path = args.output_dir / "crossvalidation_release_continuity_summary.csv"
    receipt_path = args.output_dir / "crossvalidation_release_continuity_provenance.json"
    detail.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    receipt_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
