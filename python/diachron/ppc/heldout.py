"""Pre-fit concept and language partitions with conditional predictive scoring."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from ..infer.campaign import prepare_study, run_campaign
from ..io.characters import build_character_matrix
from ..io.corpus import build_primary_corpus
from ..models.posterior import ModelConfig
from ..models.tree import PhyloTree

_MODELS = (
    ModelConfig.tree_only, ModelConfig.dated, ModelConfig.geography, ModelConfig.borrowing,
)


def _filter_lexical(dataset, predicate):
    lexemes = [form for form in dataset.lexemes if predicate(form)]
    identifiers = {form.id for form in lexemes}
    return replace(dataset, lexemes=lexemes, judgments=[
        judgment for judgment in dataset.judgments if judgment.lexeme_id in identifiers
    ])


def _select(values, rng):
    unique = sorted(values)
    n = max(1, len(unique) // 5)
    return set(rng.choice(unique, size=n, replace=False).tolist())


def _concept_split(dataset, attestations, held):
    train = _filter_lexical(dataset, lambda form: form.concept_id not in held)
    test = _filter_lexical(dataset, lambda form: form.concept_id in held)
    study = prepare_study(train, attestations, include_unobserved_languages=True)
    test_study = prepare_study(test, include_unobserved_languages=True)
    return study, test_study.character_ids, test_study.columns, None


def _language_split(dataset, attestations, held):
    train = _filter_lexical(dataset, lambda form: form.language_id not in held)
    study = prepare_study(train, attestations, include_unobserved_languages=True)

    # The evaluation vocabulary is the set selected on training taxa alone.
    # Newly encountered classes have no prediction target in this design.
    truth = build_character_matrix(
        dataset, drop_singletons=False, include_unobserved_languages=True,
    )
    by_id = {(c.concept_id, c.cognate_class): c.states for c in truth.characters}
    held_indices = [i for i, label in enumerate(study.tip_labels) if label in held]
    ids, patterns, conditioners = [], [], []
    for char_id, train_pattern in zip(study.character_ids, study.columns):
        test_pattern = by_id[char_id]
        if not any(test_pattern[i] >= 0 for i in held_indices):
            continue
        combined = list(train_pattern)
        for i in held_indices:
            combined[i] = test_pattern[i]
        ids.append(char_id)
        patterns.append(combined)
        conditioners.append(train_pattern)
    return study, ids, patterns, conditioners


def _predictive_log_scores(chains, patterns, conditioning):
    results = []
    model = chains[0].model
    if conditioning is not None:
        model = replace(model, config=replace(model.config, ascertainment=False))
    for chain in chains:
        draws = chain.draws
        for index in np.linspace(0, len(draws) - 1, min(20, len(draws)), dtype=int):
            draw = draws[index]
            tree = PhyloTree(draw.parent, draw.heights, draw.taxon)
            params = (draw.gain, draw.loss, draw.clock_rate)
            scores = np.asarray(model.character_scores(tree, *params, patterns))
            if conditioning is not None:
                scores -= model.character_scores(tree, *params, conditioning)
            results.append(scores)
    if not results or not patterns:
        raise ValueError("predictive scoring requires saved draws and held-out patterns")
    scores = logsumexp(results, axis=0) - np.log(len(results))
    if not np.isfinite(scores).all():
        raise ValueError("nonfinite held-out prediction score")
    return scores.tolist()


def _summary(character_ids, scores, baseline, rng):
    by_concept = defaultdict(list)
    for (concept, _), model_score, base in zip(character_ids, scores, baseline):
        by_concept[concept].append(model_score - base)
    shifts = np.asarray([np.mean(values) for values in by_concept.values()])
    bootstrap = rng.choice(shifts, size=(1000, shifts.size), replace=True).mean(axis=1)
    return {
        "mean_log_score_per_character": float(np.mean(scores)),
        "mean_concept_delta_vs_tree": float(shifts.mean()),
        "concept_bootstrap_interval_95": np.quantile(bootstrap, [0.025, 0.975]).tolist(),
    }


def run_heldout_prediction(out_dir: str = "runs/predict", seed: int = 0) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    corpus = build_primary_corpus()
    dataset = corpus.dataset
    rng = np.random.default_rng(seed)
    experiments = {
        "concept": (_concept_split, _select({c.id for c in dataset.concepts}, rng)),
        "language": (_language_split, _select({l.id for l in dataset.languages}, rng)),
    }
    report = {
        "study": dataset.name, "is_synthetic": dataset.is_synthetic, "split_seed": seed,
        "profile": "two-chain predictive diagnostic",
        "experiments": {},
    }
    for split, (builder, held) in experiments.items():
        study, char_ids, outcomes, conditioners = builder(dataset, corpus.attestations, held)
        results = {}
        scores = {}
        split_out = out / split
        split_out.mkdir()
        (split_out / "partition.json").write_text(json.dumps({
            "held_out_ids": sorted(held), "train_character_ids": study.character_ids,
            "test_character_ids": char_ids,
        }, indent=2))
        for build in _MODELS:
            config = build()
            chains = run_campaign(
                study, config, n_chains=2, warmup=600, iterations=1800, thin=10,
                root_age_init=corpus.truth.root_age, base_seed=5000,
                checkpoint_dir=str(split_out / "checkpoints"),
            )
            values = _predictive_log_scores(chains, outcomes, conditioners)
            scores[config.name] = values
            results[config.name] = _summary(char_ids, values, scores["tree_only"], rng)
            (split_out / f"{config.name}.json").write_text(
                json.dumps({"model": config.name, "log_scores": values}, indent=2)
            )
        report["experiments"][split] = {
            "n_train_characters": len(study.columns), "n_scored_characters": len(char_ids),
            "held_out_ids": sorted(held), "models": results,
            "scope": (
                "Conditioned on held-out characters meeting the same selection rule."
                if split == "concept" else
                "Predicts held-out taxa on training-defined classes, with external geography and dates available."
            ),
        }
    (out / "heldout.json").write_text(json.dumps(report, indent=2))
    return report
