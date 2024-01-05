"""Recompute observed campaign completion from draws and checkpoints."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from ..diagnostics import convergence as cv, topology as tp
from ..models.posterior import ModelConfig
from .verifier import _model, _read_json

_MODELS = ("tree_only", "dated", "geography", "borrowing")
_ENDPOINTS = ("root_age", "gain", "loss", "clock_rate")


def _diagnose(raw: dict) -> dict:
    """Use full precision diagnostics, independent of the CLI's rounded report."""
    scalars = {}
    for name in _ENDPOINTS:
        values = np.array([[draw[name] for draw in chain] for chain in raw["chains"]])
        scalars[name] = {
            "split_rhat": float(cv.split_rhat(values)),
            "rank_rhat": float(cv.rank_normalized_rhat(values)),
            "bulk_ess": float(cv.bulk_ess(values)),
            "tail_ess": float(cv.tail_ess(values)),
        }
    newicks = [[draw["newick"] for draw in chain] for chain in raw["chains"]]
    agreement = tp.between_chain_clades(newicks)
    mixing = [tp.topology_mixing_summary(chain) for chain in newicks]
    scalar_pass = all(
        all(math.isfinite(value) for value in endpoint.values())
        and endpoint["split_rhat"] < 1.01 and endpoint["rank_rhat"] < 1.01
        and endpoint["bulk_ess"] > 400 and endpoint["tail_ess"] > 400
        for endpoint in scalars.values()
    )
    topology_pass = (
        all(chain["distinct"] > 1 for chain in mixing)
        and agreement["mean_frequency_sd"] is not None
        and agreement["mean_frequency_sd"] < 0.01
        and agreement["maximum_chain_delta"] < 0.1
    )
    return {
        "scalars": scalars, "between_chain_clades": agreement,
        "topology_per_chain": mixing,
        "converged": bool(len(raw["chains"]) >= 4 and scalar_pass and topology_pass),
    }


def _source_matches(root: Path, commit: str) -> bool:
    """Bind the source hash inventory to trusted local Git objects."""
    record = _read_json(root / "SOURCE_INPUTS.json", 2 * 1024 * 1024)
    if record["source_commit"] != commit:
        return False
    tracked = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "-z", commit],
    ).decode().rstrip("\0").split("\0")
    if not tracked or len(tracked) > 4096 or set(record["files"]) != set(tracked):
        return False
    for name in tracked:
        content = subprocess.check_output(["git", "show", f"{commit}:{name}"])
        if hashlib.sha256(content).hexdigest() != record["files"][name]:
            return False
    return True


def _checkpoint_matches(checkpoint: dict, raw: dict, index: int, plan: dict) -> bool:
    config = raw["context"]["chain_configurations"][index]
    total = config["warmup"] + config["iterations"]
    final = raw["chains"][index][-1]
    model = _model(raw["context"]["model"])
    canonical = json.dumps(
        {key: value for key, value in asdict(model).items() if key in raw["context"]["model"]}, sort_keys=True, separators=(",", ":"), allow_nan=False)
    for name in ("parent", "heights", "taxon"):
        if checkpoint[name] != final[name]:
            return False
    state = checkpoint["state"]
    return (
        checkpoint["schema"] == 3 and checkpoint["source_commit"] == plan["source_commit"]
        and checkpoint["model_digest"] == hashlib.sha256(canonical.encode()).hexdigest()
        and checkpoint["config"] == config and checkpoint["iteration"] == total
        and checkpoint["draws"] == raw["chains"][index]
        and checkpoint["n_draws"] == len(raw["chains"][index])
        and checkpoint["proposals"] == total
        and sum(checkpoint["move_proposals"].values()) == total
        and 0 <= checkpoint["accepts"] <= total
        and sum(checkpoint["move_accepts"].values()) == checkpoint["accepts"]
        and all(0 <= count <= checkpoint["move_proposals"].get(name, -1)
                for name, count in checkpoint["move_accepts"].items())
        and checkpoint["log_post"] == final["log_posterior"]
        and all(state[name] == final[name] for name in (
            "gain", "loss", "clock_rate", "geo_sigma2", "borrow_rates", "clock_multipliers",
        ))
        and len(checkpoint["rng_state"]) == 4 and any(checkpoint["rng_state"])
    )


def verify_inference_campaign(root: Path, expected_commit: str | None) -> list[dict]:
    checks = []

    def record(name, ok, detail=None):
        checks.append({"check": f"inference:{name}", "ok": bool(ok), "detail": detail})

    plan = _read_json(root / "RUN_PLAN.json", 1024 * 1024)
    record("trusted_source", expected_commit is not None and plan["source_commit"] == expected_commit)
    record("frozen_plan", (
        plan["schema"] == 1 and plan["protocol"] == 3 and plan["component"] == "observed_inference"
        and plan["models"] == list(_MODELS) and plan["study"] == "iecor"
        and plan["chains"] >= 4 and plan["warmup"] >= 3000 and plan["iterations"] >= 12000
        and plan["thin"] == 10 and plan["initialization"] == "distance"
    ))
    if expected_commit is not None:
        record("source_inventory", _source_matches(root, expected_commit))
    status = _read_json(root / "INFERENCE.json", 1024 * 1024)
    commands = [json.loads(line) for line in (root / "COMMANDS.jsonl").read_text().splitlines()]
    record("status_identity", status["source_commit"] == plan["source_commit"]
           and status["component"] == plan["component"] and set(status["models"]) == set(_MODELS)
           and status["status"] == "completed" and status["full_project_acceptance"] is False)
    converged = []
    for name in _MODELS:
        directory = root / name
        raw = _read_json(directory / name / "bundle.json", 128 * 1024 * 1024)
        configs = raw["context"]["chain_configurations"]
        record(f"{name}:model", raw["model_name"] == name
               and raw["context"]["model"]["config"] == asdict(getattr(ModelConfig, name)())
               and raw["source_commit"] == plan["source_commit"]
               and raw["dataset_name"] == "iecor" and raw["is_synthetic"] is False)
        expected_iterations = list(range(plan["warmup"] + plan["thin"],
                                         plan["warmup"] + plan["iterations"] + 1, plan["thin"]))
        budgets_match = (
            len(configs) == raw["n_chains"] == plan["chains"]
            and len({config["seed"] for config in configs}) == plan["chains"]
            and all(all(config[key] == plan[key] for key in (
                "warmup", "iterations", "thin", "initialization",
            )) for config in configs)
            and all([draw["iteration"] for draw in chain] == expected_iterations for chain in raw["chains"])
        )
        record(f"{name}:complete_budgets", budgets_match)
        checkpoint_matches = budgets_match and all(
            _checkpoint_matches(_read_json(
                directory / "checkpoints" / f"{name}_chain{index}.json", 128 * 1024 * 1024,
            ), raw, index, plan) for index in range(plan["chains"])
        )
        record(f"{name}:final_checkpoints", checkpoint_matches)
        for action in ("sample", "diagnostics", "verify"):
            label = f"{name}-{action}"
            attempts = [command for command in commands if command["label"] == label]
            record(f"{label}:command", bool(attempts) and attempts[-1]["exit"] == 0
                   and attempts[-1]["source_commit"] == plan["source_commit"]
                   and attempts[-1]["seconds"] > 0)
        diagnostic = _diagnose(raw)
        # Keep artifact validity separate from the scientific acceptance condition.
        converged.append(diagnostic["converged"])
        record(f"{name}:diagnostic_claim", status["models"][name]["converged"] == diagnostic["converged"])
        record(f"{name}:convergence", diagnostic["converged"], diagnostic)
    record("aggregate_claim", status["all_converged"] == all(converged))
    return checks
