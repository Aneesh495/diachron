"""Known-truth recovery study (gate G05).

We generate many synthetic datasets spanning tree shape, missingness, borrowing,
and dating regimes, fit a short chain to each, and measure two things: whether the
inferred tree recovers the true clades (clade recovery rate) and whether the true
root age falls inside the posterior credible interval (credible-region coverage).

Coverage is the honesty check: if the 90% credible interval contains the truth
about 90% of the time across datasets, the uncertainty is calibrated. Systematic
over- or under-coverage is reported rather than hidden, because a posterior that
is narrower than it should be is a failure even when its point estimate is close.

Each dataset is labeled synthetic. The study runs a bounded number of short
chains; the count and per-dataset budget are arguments so a pilot and the full
200-dataset run use the same code path.
"""

from __future__ import annotations

import json
import os

import numpy as np

from ..io.simulate import simulate_dataset, simulate_attestations
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..diagnostics import topology as tp
from ..models.tree import tree_from_truth


def _true_clades(truth) -> set:
    tree = tree_from_truth(truth)
    nw = tree.newick()
    return tp._clade_signature(nw)


def _regime(i: int) -> dict:
    """Vary the generative regime across datasets to span the required space."""
    rng = np.random.default_rng(10_000 + i)
    return {
        "n_languages": int(rng.integers(8, 20)),
        "n_concepts": int(rng.integers(40, 90)),
        "root_age": float(rng.uniform(3000, 8000)),
        "gain": float(rng.uniform(0.0002, 0.0006)),
        "loss": float(rng.uniform(0.0006, 0.0015)),
        "borrow_prob": float(rng.choice([0.0, 0.15, 0.3])),
        "missing_frac": float(rng.choice([0.0, 0.08, 0.2])),
    }


def run_recovery_study(n_datasets: int = 200, out_dir: str = "runs/recovery",
                       warmup: int = 400, iterations: int = 1200) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    clade_recoveries = []
    coverage_hits = 0
    coverage_total = 0
    per_dataset = []

    for i in range(n_datasets):
        reg = _regime(i)
        dataset, truth = simulate_dataset(seed=20000 + i, name=f"rec{i}", **reg)
        attestations = simulate_attestations(dataset, truth, n_records=20, seed=20000 + i)
        study = prepare_study(dataset, attestations)

        chains = run_campaign(
            study, ModelConfig.dated(), n_chains=2,
            warmup=warmup, iterations=iterations, thin=10,
            root_age_init=truth.root_age, base_seed=30000 + i,
        )
        bundle = campaign_bundle(chains, study, ModelConfig.dated())

        # Clade recovery: MAP clades vs true clades.
        newicks = [nw for c in range(bundle.n_chains) for nw in bundle.chain_newicks(c)]
        post_clades = tp.clade_posteriors(newicks)
        true_clades = _true_clades(truth)
        if true_clades:
            recovered = sum(1 for cl in true_clades if post_clades.get(cl, 0) > 0.5)
            clade_recoveries.append(recovered / len(true_clades))

        # Root-age credible interval coverage (dated model).
        root_samples = np.array([d for ch in bundle.scalar("root_age") for d in ch])
        lo, hi = np.quantile(root_samples, [0.05, 0.95])
        coverage_total += 1
        if lo <= truth.root_age <= hi:
            coverage_hits += 1
        per_dataset.append({
            "i": i, "regime": reg,
            "clade_recovery": clade_recoveries[-1] if clade_recoveries else None,
            "root_in_90ci": bool(lo <= truth.root_age <= hi),
            "root_true": truth.root_age,
            "root_ci": [float(lo), float(hi)],
        })

    result = {
        "gate": "G05",
        "n_datasets": n_datasets,
        "mean_clade_recovery": float(np.mean(clade_recoveries)) if clade_recoveries else None,
        "root_age_90ci_coverage": coverage_hits / coverage_total if coverage_total else None,
        "coverage_target": 0.90,
        "note": "All datasets labeled synthetic known-truth.",
    }
    with open(os.path.join(out_dir, "recovery.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": result, "per_dataset": per_dataset}, f, indent=2)
    return result
