"""CLI glue for running and diagnosing sampling campaigns."""

from __future__ import annotations

import os
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

    scalars = {}
    for name in ["root_age", "gain", "loss", "clock_rate", "log_posterior"]:
        chains = np.array(bundle.scalar(name))
        scalars[name] = {
            "mean": float(np.mean(chains)),
            "split_rhat": round(cv.split_rhat(chains), 4),
            "rank_rhat": round(cv.rank_normalized_rhat(chains), 4),
            "bulk_ess": round(cv.bulk_ess(chains), 1),
            "tail_ess": round(cv.tail_ess(chains), 1),
            "mcse_mean": float(cv.mcse_mean(chains)),
        }
    out["scalars"] = scalars

    topo = []
    for c in range(bundle.n_chains):
        topo.append(tp.topology_mixing_summary(bundle.chain_newicks(c)))
    out["topology_per_chain"] = topo
    if bundle.n_chains >= 2:
        out["clade_agreement_01"] = round(
            tp.clade_agreement(bundle.chain_newicks(0), bundle.chain_newicks(1)), 3
        )

    # Convergence verdict against the declared thresholds.
    endpoints = ["root_age", "gain", "loss", "clock_rate"]
    rhat_ok = all(scalars[e]["split_rhat"] < 1.01 for e in endpoints)
    ess_ok = all(
        scalars[e]["bulk_ess"] > 400 and scalars[e]["tail_ess"] > 400 for e in endpoints
    )
    out["converged"] = bool(rhat_ok and ess_ok)
    out["thresholds"] = {"split_rhat_below": 1.01, "ess_above": 400}
    return out
