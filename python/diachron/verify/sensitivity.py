"""Recompute sensitivity summaries and identify incomplete observed designs."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import numpy as np

from ..sensitivity.run import _declared_configs
from .campaign import _checkpoint_matches, _diagnose
from .verifier import _read_json


def verify_sensitivity(root: Path) -> list[dict]:
    checks = []

    def record(name, passed, detail=None):
        checks.append({"check": f"sensitivity:{name}", "ok": bool(passed), "detail": detail})

    report = _read_json(root / "sensitivity.json", 8 * 1024 * 1024)
    plan = _read_json(root / "SENSITIVITY_PLAN.json", 2 * 1024 * 1024)
    settings = _declared_configs()[:15] if report["is_synthetic"] else _declared_configs()
    shard = plan.get("shard")
    if shard is not None:
        settings = [setting for setting in settings if setting.name == shard]
    expected = [{"name": setting.name, "lever": setting.lever, "model": asdict(setting.build())}
                for setting in settings]
    names = [setting.name for setting in settings]
    record("predeclared_settings", plan["configurations"] == expected
           and [row["config"] for row in report["configs"]] == names
           and report["n_configs"] == len(names)
           and (len(names) == 1 if shard is not None else len(names) >= 12))
    record("source", report["source_commit"] == plan["source_commit"])
    if plan["thin"] != 10 or not 1 <= plan["chains"] <= 8 or plan["warmup"] < 0 or plan["iterations"] < 10:
        record("budget_shape", False)
        return checks
    declared = {row["name"]: row["model"] for row in expected}
    baseline = None
    computed = {}
    converged = []
    for row in report["configs"]:
        name = row["config"]
        if name not in declared or row["model"] != declared[name]:
            record(f"{name}:declared_model", False)
            return checks
        model_name = declared[name]["name"]
        raw = _read_json(root / name / model_name / "bundle.json", 128 * 1024 * 1024)
        roots = np.array([[draw["root_age"] for draw in chain] for chain in raw["chains"]]).flatten()
        mean, sd = float(roots.mean()), float(roots.std())
        quantiles = np.quantile(roots, [.05, .5, .95])
        computed[name] = mean
        if name == "pop_size_3000":
            baseline = mean
        record(f"{name}:summary", abs(row["root_age_mean"] - mean) < 1e-9
               and abs(row["root_age_sd"] - sd) < 1e-9
               and np.allclose(row["root_quantiles_05_50_95"], quantiles, atol=1e-9, rtol=0))
        record(f"{name}:model", row["model"] == raw["context"]["model"]["config"]
               and row["intervention"] == raw["context"]["sensitivity"]
               and row["n_chars"] == len(raw["context"]["model"]["columns"])
               and row["n_attestations"] == len(raw["context"]["attestations"])
               and raw["source_commit"] == report["source_commit"]
               and raw["is_synthetic"] == report["is_synthetic"])
        configs = raw["context"]["chain_configurations"]
        expected_draws = plan["iterations"] // plan["thin"]
        record(f"{name}:chain_budgets", len(configs) == raw["n_chains"] == plan["chains"]
               and len({config["seed"] for config in configs}) == plan["chains"]
               and all(all(config[key] == plan[key] for key in (
                   "warmup", "iterations", "thin", "initialization",
               )) for config in configs)
               and all(len(chain) == expected_draws and all(
                   draw["iteration"] == plan["warmup"] + plan["thin"] * (index + 1)
                   for index, draw in enumerate(chain)
               ) for chain in raw["chains"]))
        if plan["iterations"] >= 12000 and plan["iterations"] % plan["thin"] == 0:
            record(f"{name}:final_checkpoints", all(
                _checkpoint_matches(_read_json(
                    root / name / "checkpoints" / f"{model_name}_chain{index}.json", 128 * 1024 * 1024,
                ), raw, index, plan) for index in range(plan["chains"])
            ))
        diagnostic = _diagnose(raw)
        converged.append(diagnostic["converged"])
        record(f"{name}:convergence_claim", diagnostic["converged"] == row["converged"])
    if baseline is None:
        record("baseline", shard is not None and shard != "pop_size_3000" and report["baseline"] is None
               and not report["levers_that_shift_root_age"]
               and all(row["root_shift_vs_baseline"] is None for row in report["configs"]))
    else:
        record("baseline", report["baseline"] == "pop_size_3000")
    if baseline is not None:
        record("shifts", all(abs(row["root_shift_vs_baseline"] - abs(computed[row["config"]] - baseline)) < 1e-9
                             for row in report["configs"]))
        record("material_shifts", report["levers_that_shift_root_age"] == [
            name for name in names if abs(computed[name] - baseline) > .15 * max(baseline, 1)
        ])
    record("aggregate_convergence_claim", report["all_converged"] == all(converged))
    return checks
