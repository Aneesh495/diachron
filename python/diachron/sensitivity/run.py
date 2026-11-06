"""Predeclared sensitivity configurations (gate G08).

The brief requires at least twelve predeclared sensitivity configurations that
test whether conclusions are stable or depend on prior and calibration choices.
Each configuration changes one lever and reruns a short inference, then records a
scalar summary (the posterior mean root age and the mean within-clade sharing) so
stable findings can be separated from prior- or calibration-dependent ones.

The configurations are declared here, before any run, and cover: coalescent
population size (tree prior strength), clock-rate prior mean, ascertainment on
versus off, singleton inclusion, loan inclusion versus removal, calibration
subsets, spatial resolution for the geographic model, and the four observation
models. The point is not that every lever moves the answer, but that the analysis
reports which ones do.
"""

from __future__ import annotations

import json
from pathlib import Path
import hashlib
import os
from dataclasses import asdict, dataclass, replace

import numpy as np

from ..io.corpus import build_primary_corpus
from ..infer.campaign import run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..util.gitmeta import source_commit, working_tree_clean
from .interventions import prepare_sensitivity_study


@dataclass
class SensitivityConfig:
    name: str
    lever: str
    build: callable


def _declared_configs() -> list[SensitivityConfig]:
    protocol = json.loads(Path("data/manifests/sensitivity-protocol.json").read_text())
    configurations = protocol["configurations"]
    if protocol["schema"] != 1 or len({row["name"] for row in configurations}) != len(configurations):
        raise ValueError("invalid sensitivity protocol")
    allowed = {"tree_only", "dated", "geography", "borrowing", "relaxed"}
    result = []
    for row in configurations:
        if row["base_model"] not in allowed:
            raise ValueError("unsupported sensitivity base model")
        factory = getattr(ModelConfig, row["base_model"])
        result.append(SensitivityConfig(
            row["name"], row["lever"],
            lambda factory=factory, overrides=row["overrides"]: replace(factory(), **overrides),
        ))
    return result


def run_sensitivity(
    out_dir: str = "runs/sensitivity", warmup: int = 500, iterations: int = 1500,
    *, dataset=None, attestations=None, frozen_study=None, n_chains: int = 2,
    jobs: int = 1, resume: bool = False, only_configuration: str | None = None,
) -> dict:
    if dataset is None:
        corpus = build_primary_corpus()
        dataset, attestations = corpus.dataset, corpus.attestations
    attestations = list(attestations or [])
    if not dataset.is_synthetic and not working_tree_clean():
        raise ValueError("observed sensitivity requires a clean committed checkout")
    commit = source_commit()
    configurations = _declared_configs()
    # Keep the earlier synthetic diagnostic protocol reproducible.
    if dataset.is_synthetic:
        configurations = configurations[:15]
    if only_configuration is not None:
        configurations = [setting for setting in configurations if setting.name == only_configuration]
        if not configurations:
            raise ValueError("unknown sensitivity configuration")
    fingerprint = json.dumps(
        {"dataset": asdict(dataset), "attestations": [asdict(row) for row in attestations]},
        sort_keys=True, separators=(",", ":"), allow_nan=False,
    )
    plan = {
        "schema": 1, "source_commit": commit, "study": dataset.name, "shard": only_configuration,
        "input_sha256": hashlib.sha256(fingerprint.encode()).hexdigest(),
        "frozen_study": frozen_study, "chains": n_chains, "warmup": warmup,
        "iterations": iterations, "thin": 10, "jobs": jobs, "initialization": "distance",
        "configurations": [{"name": setting.name, "lever": setting.lever,
                            "model": asdict(setting.build())} for setting in configurations],
    }
    out = Path(out_dir)
    plan_path = out / "SENSITIVITY_PLAN.json"
    if out.exists():
        if (out / "SHARDS.json").exists():
            raise ValueError("merged evidence is immutable; start a new campaign")
        if not resume or not plan_path.is_file():
            raise FileExistsError("sensitivity output exists; resume with the same plan")
        if json.loads(plan_path.read_text()) != plan:
            raise ValueError("sensitivity source, input or budget differs")
    else:
        out.mkdir(parents=True)
        _write(plan_path, plan)
    rows = []
    for setting in configurations:
        study, intervention = prepare_sensitivity_study(dataset, attestations, setting.name)
        config = setting.build()
        if config.use_calibration and study.calibrations:
            calibration = study.calibrations[0]
            start_age = (calibration.age_low + calibration.age_high) / 2
        else:
            start_age = 8000.0 if dataset.name == "iecor" else 6000.0
        directory = out / setting.name
        chains = run_campaign(
            study, config, n_chains=n_chains, warmup=warmup, iterations=iterations, thin=10,
            root_age_init=start_age, base_seed=40000, jobs=jobs, resume=resume,
            checkpoint_dir=str(directory / "checkpoints"),
        )
        bundle = campaign_bundle(chains, study, config, source_commit=commit)
        bundle.context["sensitivity"] = intervention
        if frozen_study is not None:
            bundle.context["frozen_study"] = frozen_study
        bundle.save(str(directory / config.name))
        from ..cli.sampling import run_diagnose

        diagnostic = run_diagnose(str(directory / config.name))
        _write(directory / "diagnostics.json", diagnostic)
        root = np.array(bundle.scalar("root_age")).flatten()
        rows.append({
            "config": setting.name, "lever": setting.lever,
            "root_age_mean": float(root.mean()), "root_age_sd": float(root.std()),
            "n_chars": len(study.columns), "n_attestations": len(study.attestations),
            "model": asdict(config), "intervention": intervention,
            "root_quantiles_05_50_95": np.quantile(root, [0.05, 0.5, 0.95]).tolist(),
            "converged": diagnostic["converged"],
        })
        _write(out / "PROGRESS.json", {
            "source_commit": commit, "completed_configurations": [row["config"] for row in rows],
            "status": "running", "required_configurations": len(configurations),
        })
    result = summarize_sensitivity(rows, dataset.name, dataset.is_synthetic, commit)
    _write(out / "sensitivity.json", result)
    _write(out / "PROGRESS.json", {
        "source_commit": commit, "status": "completed",
        "completed_configurations": [row["config"] for row in rows],
    })
    return result


def summarize_sensitivity(rows, study_name, is_synthetic, commit):
    baseline = next((row for row in rows if row["config"] == "pop_size_3000"), None)
    for row in rows:
        row["root_shift_vs_baseline"] = (
            abs(row["root_age_mean"] - baseline["root_age_mean"]) if baseline else None
        )
    return {
        "gate": "G08", "source_commit": commit, "is_synthetic": is_synthetic,
        "study": study_name, "profile": "saved-chain sensitivity campaign",
        "n_configs": len(rows), "configs": rows, "baseline": "pop_size_3000" if baseline else None,
        "all_converged": bool(rows) and all(row["converged"] for row in rows),
        "levers_that_shift_root_age": [row["config"] for row in rows if baseline and
            row["root_shift_vs_baseline"] > .15 * max(baseline["root_age_mean"], 1)],
        "note": "Shifts are descriptive. Unconverged fits do not establish robustness.",
    }


def merge_sensitivity_shards(shards: list[Path], out: Path):
    """Assemble complete shards without duplicating large raw draw files."""
    import shutil

    if out.exists():
        raise FileExistsError("merged sensitivity output already exists")
    common = None
    rows = {}
    sources = {}
    for path in shards:
        plan = json.loads((path / "SENSITIVITY_PLAN.json").read_text())
        report = json.loads((path / "sensitivity.json").read_text())
        name = plan["shard"]
        if not name or name in rows or report["n_configs"] != 1:
            raise ValueError("duplicate or invalid sensitivity shard")
        identity = {key: value for key, value in plan.items() if key not in {"shard", "configurations"}}
        if common is not None and common != identity:
            raise ValueError("shards have different source, input, or budgets")
        common = identity
        row = report["configs"][0]
        if row["config"] != name or report["source_commit"] != plan["source_commit"]:
            raise ValueError("shard report identity differs")
        rows[name] = row
        sources[name] = path
    settings = _declared_configs()
    if set(rows) != {setting.name for setting in settings}:
        raise ValueError("observed sensitivity merge requires every predeclared configuration")
    out.mkdir(parents=True)
    common["configurations"] = [{"name": setting.name, "lever": setting.lever,
                                "model": asdict(setting.build())} for setting in settings]
    common["shard"] = None
    _write(out / "SENSITIVITY_PLAN.json", common)
    for name, source in sources.items():
        if not (source / name).is_dir() or (source / name).is_symlink():
            raise ValueError("missing shard model directory")
        # Hardlinks preserve identical immutable inputs without a second allocation.
        # The assembled evidence is read-only until a new campaign is generated.
        shutil.copytree(source / name, out / name, copy_function=os.link)
    result = summarize_sensitivity(
        [rows[setting.name] for setting in settings], common["study"], False, common["source_commit"],
    )
    _write(out / "sensitivity.json", result)
    _write(out / "SHARDS.json", {name: str(path) for name, path in sources.items()})
    return result


def _write(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)
