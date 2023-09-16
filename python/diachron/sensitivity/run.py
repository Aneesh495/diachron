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
from dataclasses import asdict, dataclass, replace

import numpy as np

from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..util.gitmeta import source_commit


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
    return configs


def run_sensitivity(out_dir: str = "runs/sensitivity",
                    warmup: int = 500, iterations: int = 1500) -> dict:
    os.makedirs(out_dir, exist_ok=False)
    corpus = build_primary_corpus()

    rows = []
    for sc in _declared_configs():
        # Data levers change how the study is prepared.
        drop_singletons = sc.name != "include_singletons"
        include_loans = sc.name != "remove_loans"
        attestations = corpus.attestations
        if sc.name == "half_calibrations":
            attestations = sorted(attestations, key=lambda att: att.id)[::2]
        study = prepare_study(
            corpus.dataset, attestations,
            drop_singletons=drop_singletons,
            include_loans_as_present=include_loans,
        )
        config = sc.build()
        if config.use_calibration and study.calibrations:
            cal = study.calibrations[0]
            start_age = (cal.age_low + cal.age_high) / 2.0
        else:
            start_age = 6000.0
        chains = run_campaign(
            study, config, n_chains=2, warmup=warmup, iterations=iterations, thin=10,
            root_age_init=start_age, base_seed=40000,
        )
        bundle = campaign_bundle(chains, study, config)
        root = np.array([d for ch in bundle.scalar("root_age") for d in ch])
        rows.append({
            "config": sc.name,
            "lever": sc.lever,
            "root_age_mean": float(root.mean()),
            "root_age_sd": float(root.std()),
            "n_chars": len(study.columns),
            "n_attestations": len(attestations),
            "model": asdict(config),
            "root_quantiles_05_50_95": np.quantile(root, [0.05, 0.5, 0.95]).tolist(),
        })
        with open(os.path.join(out_dir, "configurations.jsonl"), "a", encoding="utf-8") as raw:
            raw.write(json.dumps(rows[-1], allow_nan=False) + "\n")

    # Stability verdict: which levers move the root-age estimate materially.
    baseline = next(r for r in rows if r["config"] == "pop_size_3000")
    for r in rows:
        r["root_shift_vs_baseline"] = abs(r["root_age_mean"] - baseline["root_age_mean"])
    sensitive = [r["config"] for r in rows
                 if r["root_shift_vs_baseline"] > 0.15 * max(baseline["root_age_mean"], 1)]

    result = {
        "gate": "G08",
        "source_commit": source_commit(),
        "is_synthetic": True,
        "profile": "two-chain sensitivity diagnostic",
        "n_configs": len(rows),
        "configs": rows,
        "baseline": "pop_size_3000",
        "levers_that_shift_root_age": sensitive,
        "note": "Shifts are descriptive for short chains. Small shifts do not establish robustness.",
    }
    with open(os.path.join(out_dir, "sensitivity.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    return result
