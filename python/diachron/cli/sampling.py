"""CLI glue for running and diagnosing sampling campaigns."""

from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path
import time

import numpy as np

from ..io.corpus import build_primary_corpus
from ..models.posterior import ModelConfig
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..infer.artifacts import PosteriorBundle
from ..diagnostics import convergence as cv
from ..diagnostics import topology as tp
from ..util.gitmeta import source_commit


_MODEL_CONFIGS = {
    "tree_only": ModelConfig.tree_only,
    "dated": ModelConfig.dated,
    "relaxed": ModelConfig.relaxed,
    "geography": ModelConfig.geography,
    "borrowing": ModelConfig.borrowing,
}


def _profile_budget(profile: str, warmup: int, iterations: int) -> tuple[int, int]:
    """CPU diagnostic profile uses shorter runs than the full acceptance profile.

    The profiles are named distinctly so a quick diagnostic run can never be
    mistaken for the full campaign.
    """
    if profile == "cpu":
        return warmup, iterations
    if profile == "full":
        return max(warmup, 3000), max(iterations, 12000)
    return warmup, iterations


def run_sample(
    *,
    model_name: str,
    n_chains: int,
    warmup: int,
    iterations: int,
    thin: int,
    out_dir: str,
    profile: str = "cpu",
    resume: bool = False,
    jobs: int = 1,
    initialization: str = "distance",
    height_resampled_spr: bool = False,
    temperature_ladder: tuple[float, ...] = (1.0,),
    attachment_resampling: bool = False,
    geography_error: float | None = None,
    observed_study: str | None = None,
    observed_source: Path | None = None,
    catalog: Path | None = None,
    edh_source: Path | None = None,
    edh_geography: Path | None = None,
) -> dict:
    if model_name not in _MODEL_CONFIGS:
        raise ValueError(f"unknown model {model_name}; choose from {list(_MODEL_CONFIGS)}")
    if n_chains < 1 or n_chains > 8 or thin < 1 or iterations < 1 or warmup < 0:
        raise ValueError("invalid bounded chain budget")
    config = _MODEL_CONFIGS[model_name]()
    if geography_error is not None:
        if not config.use_geography:
            raise ValueError("--geography-error applies to the geographic model")
        config.geography_error_degrees = geography_error
    if observed_study is None:
        if any((observed_source, catalog, edh_source, edh_geography)):
            raise ValueError("source and catalog arguments require --study")
        corpus = build_primary_corpus()
        study = prepare_study(corpus.dataset, corpus.attestations)
        root_age_init = corpus.truth.root_age
    else:
        from ..io.freeze import load_study
        from ..io.types import Attestation

        if observed_source is None or catalog is None:
            raise ValueError("--study requires --source and --catalog")
        if observed_study == "abvd" and config.use_calibration:
            raise ValueError("the ABVD subset supports tree_only without independent dates")
        if bool(edh_source) != bool(edh_geography) or (edh_source and observed_study != "iecor"):
            raise ValueError("EDH requires both source inputs and the IE-CoR study")
        dataset, _, report = load_study(observed_study, observed_source, catalog)
        attestations = []
        if edh_source:
            from ..io.epidoc import load_edh

            evidence = load_edh(edh_source, edh_geography)
            attestations = [Attestation(**record) for record in evidence["attestations"]]
        study = prepare_study(dataset, attestations)
        root_age_init = 8000.0 if observed_study == "iecor" else 3500.0
    warmup, iterations = _profile_budget(profile, warmup, iterations)

    t0 = time.time()
    chains = run_campaign(
        study,
        config,
        n_chains=n_chains,
        warmup=warmup,
        iterations=iterations,
        thin=thin,
        root_age_init=root_age_init,
        checkpoint_dir=os.path.join(out_dir, "checkpoints"),
        resume=resume,
        jobs=jobs,
        initialization=initialization,
        height_resampled_spr=height_resampled_spr,
        temperature_ladder=temperature_ladder,
        attachment_resampling=attachment_resampling,
    )
    bundle = campaign_bundle(chains, study, config, source_commit=source_commit())
    if observed_study:
        bundle.context["frozen_study"] = report
    bundle_dir = os.path.join(out_dir, model_name)
    bundle.save(bundle_dir)
    elapsed = time.time() - t0

    accept = [round(c.acceptance_rate, 3) for c in chains]
    return {
        "model": model_name,
        "profile": profile,
        "study": study.dataset.name,
        "is_synthetic": study.dataset.is_synthetic,
        "dated_tips": len(chains[0].model.tip_dates),
        "chains": n_chains,
        "workers": min(n_chains, jobs),
        "initialization": initialization,
        "height_resampled_spr": height_resampled_spr,
        "attachment_resampling": attachment_resampling,
        "inverse_temperatures": list(temperature_ladder),
        "warmup": warmup,
        "iterations": iterations,
        "draws_per_chain": bundle.draws_per_chain,
        "acceptance_rates": accept,
        "bundle_dir": bundle_dir,
        "seconds": round(elapsed, 1),
    }


def run_diagnose(bundle_dir: str) -> dict:
    bundle = PosteriorBundle.load(bundle_dir)
    out: dict = {"model": bundle.model_name, "n_chains": bundle.n_chains}

    def finite(value: float) -> float | None:
        return float(value) if np.isfinite(value) else None

    scalars = {}
    for name in ["root_age", "gain", "loss", "clock_rate", "log_posterior"]:
        chains = np.array(bundle.scalar(name))
        scalars[name] = {
            "mean": float(np.mean(chains)),
            "split_rhat": finite(cv.split_rhat(chains)),
            "rank_rhat": finite(cv.rank_normalized_rhat(chains)),
            "bulk_ess": float(cv.bulk_ess(chains)),
            "tail_ess": float(cv.tail_ess(chains)),
            "mcse_mean": finite(float(cv.mcse_mean(chains))),
        }
    out["scalars"] = scalars

    topo = []
    for c in range(bundle.n_chains):
        topo.append(tp.topology_mixing_summary(bundle.chain_newicks(c)))
    out["topology_per_chain"] = topo
    agreement = tp.between_chain_clades([
        bundle.chain_newicks(chain) for chain in range(bundle.n_chains)
    ])
    out["between_chain_clades"] = agreement
    out["move_statistics"] = bundle.context.get("move_statistics")
    if "replica_exchange" in bundle.context:
        out["replica_exchange"] = bundle.context["replica_exchange"]

    out["clade_agreement"] = [
        {"chains": [a, b], "score": round(
            tp.clade_agreement(bundle.chain_newicks(a), bundle.chain_newicks(b)), 4,
        )}
        for a, b in combinations(range(bundle.n_chains), 2)
    ]

    endpoints = ["root_age", "gain", "loss", "clock_rate"]
    rhat_ok = all(
        scalars[e][metric] is not None and scalars[e][metric] < 1.01
        for e in endpoints for metric in ("split_rhat", "rank_rhat")
    )
    ess_ok = all(
        scalars[e]["bulk_ess"] > 400 and scalars[e]["tail_ess"] > 400 for e in endpoints
    )
    topology_moves = all(chain["distinct"] > 1 for chain in topo)
    clades_agree = (
        agreement["mean_frequency_sd"] is not None
        and agreement["mean_frequency_sd"] < 0.01
        and agreement["maximum_chain_delta"] < 0.1
    )
    out["converged"] = bool(
        bundle.n_chains >= 4 and rhat_ok and ess_ok and topology_moves and clades_agree
    )
    out["thresholds"] = {
        "split_rhat_below": 1.01, "ess_above": 400,
        "mean_clade_frequency_sd_below": 0.01, "maximum_chain_clade_delta_below": 0.1,
    }
    return out
