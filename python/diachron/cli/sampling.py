"""CLI glue for running and diagnosing sampling campaigns."""

from __future__ import annotations

import os
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
) -> dict:
    if model_name not in _MODEL_CONFIGS:
        raise ValueError(f"unknown model {model_name}; choose from {list(_MODEL_CONFIGS)}")
    corpus = build_primary_corpus()
    study = prepare_study(corpus.dataset, corpus.attestations)
    config = _MODEL_CONFIGS[model_name]()
    warmup, iterations = _profile_budget(profile, warmup, iterations)

    t0 = time.time()
    chains = run_campaign(
        study,
        config,
        n_chains=n_chains,
        warmup=warmup,
        iterations=iterations,
        thin=thin,
        root_age_init=corpus.truth.root_age,
        checkpoint_dir=os.path.join(out_dir, "checkpoints"),
    )
    bundle = campaign_bundle(chains, study, config, source_commit=source_commit())
    bundle_dir = os.path.join(out_dir, model_name)
    bundle.save(bundle_dir)
    elapsed = time.time() - t0

    accept = [round(c.acceptance_rate, 3) for c in chains]
    return {
        "model": model_name,
        "profile": profile,
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
