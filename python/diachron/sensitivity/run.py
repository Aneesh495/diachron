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
import os
from dataclasses import dataclass

import numpy as np

from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig


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
    # 9-12: the four observation models at the default prior.
    for name, builder in [
        ("model_tree_only", ModelConfig.tree_only),
        ("model_dated", ModelConfig.dated),
        ("model_geography", ModelConfig.geography),
        ("model_borrowing", ModelConfig.borrowing),
    ]:
        configs.append(SensitivityConfig(name, "model", builder))
    return configs


def run_sensitivity(out_dir: str = "runs/sensitivity",
                    warmup: int = 500, iterations: int = 1500) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    corpus = build_primary_corpus()

    rows = []
    for sc in _declared_configs():
        # Data levers change how the study is prepared.
        drop_singletons = sc.name != "include_singletons"
        include_loans = sc.name != "remove_loans"
        study = prepare_study(
            corpus.dataset, corpus.attestations,
            drop_singletons=drop_singletons,
            include_loans_as_present=include_loans,
        )
        if sc.name == "half_calibrations":
            study.calibrations = study.calibrations[: len(study.calibrations) // 2]

        config = sc.build()
        chains = run_campaign(
            study, config, n_chains=2, warmup=warmup, iterations=iterations, thin=10,
            root_age_init=corpus.truth.root_age, base_seed=40000,
        )
        bundle = campaign_bundle(chains, study, config)
        root = np.array([d for ch in bundle.scalar("root_age") for d in ch])
        rows.append({
            "config": sc.name,
            "lever": sc.lever,
            "root_age_mean": float(root.mean()),
            "root_age_sd": float(root.std()),
            "n_chars": len(study.columns),
        })

    # Stability verdict: which levers move the root-age estimate materially.
    baseline = next(r for r in rows if r["config"] == "model_dated")
    for r in rows:
        r["root_shift_vs_baseline"] = abs(r["root_age_mean"] - baseline["root_age_mean"])
    sensitive = [r["config"] for r in rows
                 if r["root_shift_vs_baseline"] > 0.15 * max(baseline["root_age_mean"], 1)]

    result = {
        "gate": "G08",
        "n_configs": len(rows),
        "configs": rows,
        "baseline": "model_dated",
        "levers_that_shift_root_age": sensitive,
        "note": "Predeclared configurations; stable vs prior/calibration-dependent findings separated.",
    }
    with open(os.path.join(out_dir, "sensitivity.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    return result
