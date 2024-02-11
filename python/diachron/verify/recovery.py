"""Recompute known-truth coverage from the archived per-dataset draws."""

from __future__ import annotations

from collections import defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .verifier import _read_json


def verify_recovery(summary_path: Path, rows_path: Path, expected_commit: str) -> dict:
    summary = _read_json(summary_path, 1024 * 1024)
    if rows_path.stat().st_size > 32 * 1024 * 1024:
        raise ValueError("oversized compressed recovery rows")
    with rows_path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != summary["raw_rows_sha256"] or summary["source_commit"] != expected_commit:
        raise ValueError("recovery rows or source identity differ")
    matched = summary.get("protocol", 2) == 3
    data_base, chain_base = (52000, 53000) if matched else (20000, 30000)
    borrow_field = "borrow_regime" if matched else "borrow_prob"
    group_key = "by_borrow_regime" if matched else "by_borrow_probability"
    date_field = "n_dated_tips" if matched else "n_attestations"
    probes = set(np.linspace(0, summary["n_datasets"] - 1, 12, dtype=int)) if matched else set()
    histogram = [0] * 10
    completed = []
    groups = defaultdict(list)
    regimes = defaultdict(set)
    failed = count = 0
    sample_count = summary["chains"] * (summary["iterations"] // summary["thin"])
    with gzip.open(rows_path, "rb") as stream:
        while line := stream.readline(256 * 1024 + 1):
            if len(line) > 256 * 1024 or count >= 5000:
                raise ValueError("recovery row or dataset budget exceeded")
            row = json.loads(line)
            if row["i"] != count or row["data_seed"] != data_base + count or row["chain_seed"] != chain_base + count:
                raise ValueError("recovery seeds or dataset order differ")
            count += 1
            if row["status"] == "failed":
                failed += 1
                continue
            if row["status"] != "complete":
                raise ValueError("unrecognized recovery outcome")
            regime = row["regime"]
            draws, masses = row["root_samples"], row["true_clade_masses"]
            if (len(draws) != sample_count or not draws or len(masses) != regime["n_languages"] - 2
                    or not masses or any(not math.isfinite(value) or value < 0 for value in draws)
                    or any(not 0 <= mass <= 1 for mass in masses) or row["root_true"] != regime["root_age"]):
                raise ValueError("recovery samples or truth dimensions differ")
            low, high = np.quantile(draws, [.05, .95])
            covered = bool(low <= row["root_true"] <= high)
            recovered = sum(mass > .5 for mass in masses) / len(masses)
            if (not np.allclose([low, high], row["root_ci"], atol=1e-9, rtol=0)
                    or row["root_in_90ci"] != covered or abs(row["clade_recovery"] - recovered) > 1e-12):
                raise ValueError("reported recovery row differs from its draws")
            if matched:
                if row["source_commit"] != expected_commit or row["chain_seeds"] != [
                    row["chain_seed"] + 7919 * (chain + 1) for chain in range(4)
                ]:
                    raise ValueError("matched recovery source or chain seeds differ")
                rank = sum(value < row["root_true"] for value in draws)
                if row["root_rank"] != rank:
                    raise ValueError("matched recovery rank differs from raw samples")
                histogram[min(9, 10 * rank // (len(draws) + 1))] += 1
                if row["i"] in probes:
                    from ..sensitivity.matched import generate_matched_case

                    _, _, _, expected_regime, truth = generate_matched_case(row["i"])
                    if expected_regime != regime or truth != row["truth"]:
                        raise ValueError("matched generator did not reproduce the saved truth")
            completed.append((covered, recovered))
            groups[str(regime[borrow_field])].append((covered, recovered))
            for key in ("n_languages", borrow_field, "missing_frac"):
                regimes[key].add(regime[key])
            regimes[date_field].add(regime[date_field] if matched else row[date_field])
    if not completed:
        raise ValueError("recovery has no completed datasets")
    covered = sum(value[0] for value in completed)
    coverage = covered / len(completed)
    recovery = sum(value[1] for value in completed) / len(completed)
    valid = (count >= 200 and count == summary["n_datasets"]
             and len(completed) == summary["completed"] and failed == summary["failed"]
             and covered == summary["covered"] and abs(coverage - summary["root_age_90ci_coverage"]) < 1e-12
             and abs(recovery - summary["mean_clade_recovery"]) < 1e-12
             and summary["coverage_target"] == .9 and all(len(values) > 1 for values in regimes.values())
             and set(summary[group_key]) == set(groups))
    if matched:
        valid &= (summary["root_rank_histogram"] == histogram and summary["chains"] == 4
                  and summary["warmup"] == 800 and summary["iterations"] == 3200 and summary["thin"] == 10)
    for key, outcomes in groups.items():
        measured = summary[group_key][key]
        valid &= (measured["datasets"] == len(outcomes)
                  and abs(measured["coverage"] - sum(row[0] for row in outcomes) / len(outcomes)) < 1e-12
                  and abs(measured["mean_clade_recovery"] - sum(row[1] for row in outcomes) / len(outcomes)) < 1e-12)
    return {"artifact_valid": bool(valid), "source_commit": expected_commit, "datasets": count,
            "completed": len(completed), "failed": failed, "covered": covered, "coverage": coverage,
            "coverage_target": .9, "clade_recovery": recovery, "calibration_target_met": coverage >= .9,
            "matched_design": matched, "generator_probes": len(probes),
            "regimes": {key: sorted(values) for key, values in regimes.items()}}
