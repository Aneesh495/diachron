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
    configs = []

    def dated(pop=3000.0, ascertain=True):
        cfg = ModelConfig.dated()
        cfg.pop_size = pop
        cfg.ascertainment = ascertain
        return cfg

    # 1-3: tree prior strength via coalescent population size.
    for pop in [1000.0, 3000.0, 9000.0]:
        configs.append(SensitivityConfig(
            f"pop_size_{int(pop)}", "tree_prior", lambda p=pop: dated(pop=p)))
    # 4: ascertainment off.
    configs.append(SensitivityConfig(
        "ascertainment_off", "ascertainment", lambda: dated(ascertain=False)))
    # 5: singleton inclusion (handled via data prep flag below).
    configs.append(SensitivityConfig(
        "include_singletons", "data", lambda: dated()))
    # 6: loan removal (data prep flag).
    configs.append(SensitivityConfig(
        "remove_loans", "data", lambda: dated()))
    # 7: calibration subset (half the attestations).
    configs.append(SensitivityConfig(
        "half_calibrations", "calibration", lambda: dated()))
    # 8: no calibration (tree_only).
    configs.append(SensitivityConfig(
        "no_calibration", "calibration", ModelConfig.tree_only))
    configs.extend([
        SensitivityConfig("gamma_rate_heterogeneity", "clock",
                          lambda: replace(dated(), rate_alpha=0.7)),
        SensitivityConfig("relaxed_clock", "clock", ModelConfig.relaxed),
        SensitivityConfig("planar_geography", "geography", ModelConfig.geography),
        SensitivityConfig("regional_geography_3", "geography",
                          lambda: replace(ModelConfig.geography(), geo_region_count=3)),
        SensitivityConfig("regional_geography_6", "geography",
                          lambda: replace(ModelConfig.geography(), geo_region_count=6)),
        SensitivityConfig("borrowing_default", "contact", ModelConfig.borrowing),
        SensitivityConfig("borrowing_tighter_prior", "contact",
                          lambda: replace(ModelConfig.borrowing(), borrowing_prior_rate=30.0)),
    ])
    configs.extend([
        SensitivityConfig("broad_all_tip_dates", "calibration", ModelConfig.dated),
        SensitivityConfig("alternative_cognacy_merge", "cognacy", ModelConfig.dated),
        SensitivityConfig("catalog_monophyly", "taxonomy", ModelConfig.dated),
        SensitivityConfig("geographic_error_5deg", "geography",
                          lambda: replace(ModelConfig.geography(), geography_error_degrees=5.0)),
        SensitivityConfig("diffusion_wide_prior", "geography",
                          lambda: replace(ModelConfig.geography(), diffusion_prior_log_sd=3.0)),
    ])
    return configs


def run_sensitivity(
    out_dir: str = "runs/sensitivity", warmup: int = 500, iterations: int = 1500,
    *, dataset=None, attestations=None, frozen_study=None, n_chains: int = 2,
    jobs: int = 1, resume: bool = False,
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
    fingerprint = json.dumps(
        {"dataset": asdict(dataset), "attestations": [asdict(row) for row in attestations]},
        sort_keys=True, separators=(",", ":"), allow_nan=False,
    )
    plan = {
        "schema": 1, "source_commit": commit, "study": dataset.name,
        "input_sha256": hashlib.sha256(fingerprint.encode()).hexdigest(),
        "frozen_study": frozen_study, "chains": n_chains, "warmup": warmup,
        "iterations": iterations, "thin": 10, "jobs": jobs, "initialization": "distance",
        "configurations": [{"name": setting.name, "lever": setting.lever,
                            "model": asdict(setting.build())} for setting in configurations],
    }
    out = Path(out_dir)
    plan_path = out / "SENSITIVITY_PLAN.json"
    if out.exists():
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
    baseline = next(row for row in rows if row["config"] == "pop_size_3000")
    for row in rows:
        row["root_shift_vs_baseline"] = abs(row["root_age_mean"] - baseline["root_age_mean"])
    result = {
        "gate": "G08", "source_commit": commit, "is_synthetic": dataset.is_synthetic,
        "study": dataset.name, "profile": "saved-chain sensitivity campaign",
        "n_configs": len(rows), "configs": rows, "baseline": "pop_size_3000",
        "all_converged": all(row["converged"] for row in rows),
        "levers_that_shift_root_age": [row["config"] for row in rows if
            row["root_shift_vs_baseline"] > .15 * max(baseline["root_age_mean"], 1)],
        "note": "Shifts are descriptive. Unconverged fits do not establish robustness.",
    }
    _write(out / "sensitivity.json", result)
    _write(out / "PROGRESS.json", {
        "source_commit": commit, "status": "completed",
        "completed_configurations": [row["config"] for row in rows],
    })
    return result


def _write(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)
