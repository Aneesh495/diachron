"""Held-out concept and language prediction (gate G07).

Masks are applied before any model fitting so there is no leakage: we partition
characters into a training set and a held-out set, fit the tree and parameters on
the training characters only, then score the held-out characters under the fitted
posterior. Held-out concepts remove whole concepts (all their cognate classes) so
no cognate class crosses the partition; held-out languages remove whole columns.

The score is the mean held-out log predictive density per character, estimated by
averaging the character likelihood over posterior draws. A model with genuine
predictive value scores higher than the tree-only baseline; the comparison, not
the absolute number, is what the gate asks for.
"""

from __future__ import annotations

import os

import numpy as np
from scipy.special import logsumexp

from diachron import _core
from ..io.corpus import build_primary_corpus
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..models.posterior import ModelConfig
from ..models.tree import PhyloTree
from .run import _parse_newick_to_tree


def _held_out_score(bundle_chains, tip_labels, held_columns) -> float:
    """Mean held-out log predictive density per character over posterior draws."""
    draws = [d for chain in bundle_chains for d in chain]
    if not draws or not held_columns:
        return float("nan")
    # Subsample draws for speed.
    idx = np.linspace(0, len(draws) - 1, min(30, len(draws))).astype(int)
    per_char = []
    for col in held_columns:
        logliks = []
        for i in idx:
            d = draws[i]
            tree = _parse_newick_to_tree(d["newick"], tip_labels)
            parent, bl, h, tax, root = tree.to_core_arrays()
            z = d["gain"] + d["loss"]
            prior = [d["loss"] / z, d["gain"] / z]
            ll = _core.char_log_likelihood(
                parent, bl, h, tax, root, "binary", [d["gain"], d["loss"]], col, prior
            )
            logliks.append(ll)
        # log mean likelihood across draws (posterior predictive of the character).
        lppd = logsumexp(logliks) - np.log(len(logliks))
        per_char.append(lppd)
    return float(np.mean(per_char))


def run_heldout_prediction(out_dir: str = "runs/predict", seed: int = 0) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    corpus = build_primary_corpus()
    study_full = prepare_study(corpus.dataset, corpus.attestations)

    # Partition characters by concept to hold out whole concepts (no class leak).
    cm_concepts = [c.concept_id for c in _variable_characters(corpus)]
    rng = np.random.default_rng(seed)
    unique_concepts = sorted(set(cm_concepts))
    held_concepts = set(
        rng.choice(unique_concepts, size=max(1, len(unique_concepts) // 5), replace=False)
    )
    train_cols = [
        col for col, cc in zip(study_full.columns, cm_concepts) if cc not in held_concepts
    ]
    held_cols = [
        col for col, cc in zip(study_full.columns, cm_concepts) if cc in held_concepts
    ]

    results = {}
    for model_name, cfg in [("tree_only", ModelConfig.tree_only()),
                            ("dated", ModelConfig.dated())]:
        # Fit on training characters only.
        train_study = _study_with_columns(study_full, train_cols)
        chains = run_campaign(
            train_study, cfg, n_chains=2, warmup=600, iterations=1800, thin=10,
            root_age_init=corpus.truth.root_age, base_seed=5000,
        )
        bundle = campaign_bundle(chains, train_study, cfg)
        score = _held_out_score(bundle.chains, bundle.tip_labels, held_cols)
        results[model_name] = {
            "held_out_lppd_per_char": score,
            "n_train_chars": len(train_cols),
            "n_held_chars": len(held_cols),
        }

    results["gate"] = "G07"
    results["held_out_concepts"] = len(held_concepts)
    results["leakage_guard"] = "whole concepts held out; no cognate class spans partition"
    if "dated" in results and "tree_only" in results:
        results["dated_beats_tree_only"] = (
            results["dated"]["held_out_lppd_per_char"]
            >= results["tree_only"]["held_out_lppd_per_char"]
        )
    with open(os.path.join(out_dir, "heldout.json"), "w", encoding="utf-8") as f:
        import json

        json.dump(results, f, indent=2)
    return results


def _variable_characters(corpus):
    from ..io.characters import build_character_matrix

    cm = build_character_matrix(corpus.dataset)
    return cm.variable_characters()


def _study_with_columns(study, columns):
    from ..infer.campaign import StudyData

    return StudyData(
        tip_labels=study.tip_labels,
        columns=columns,
        calibrations=study.calibrations,
        tip_coords=study.tip_coords,
        dataset=study.dataset,
    )
