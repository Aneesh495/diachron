"""Independent checks of held-out score arrays and fixed-probe predictive replay."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, replace
import json
import math
from pathlib import Path

import numpy as np

from ..io.characters import build_character_matrix
from ..models.posterior import ModelConfig
from ..models.tree import PhyloTree

_MODELS = ("tree_only", "dated", "geography", "borrowing")


def _read(path: Path):
    if path.stat().st_size > 128 * 1024 * 1024:
        raise ValueError(f"oversized predictive artifact: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def _masked_dataset(dataset, kind, held, keep_held):
    forms = [form for form in dataset.lexemes
             if (((form.concept_id if kind == "concept" else form.language_id) in held) == keep_held)]
    form_ids = {form.id for form in forms}
    return replace(dataset, lexemes=forms,
                   judgments=[j for j in dataset.judgments if j.lexeme_id in form_ids])


def _test_patterns(dataset, partition):
    """Build evaluation outcomes independently from the saved score producer."""
    if partition["kind"] == "concept":
        evidence = _masked_dataset(dataset, "concept", set(partition["identifiers"]), True)
        matrix = build_character_matrix(evidence, include_unobserved_languages=True)
        chars = matrix.variable_characters()
    else:
        matrix = build_character_matrix(dataset, drop_singletons=False, include_unobserved_languages=True)
        chars = matrix.characters
    return {(char.concept_id, char.cognate_class): char.states for char in chars}


def _replay_scores(bundle, indices, character_ids, test_lookup, partition):
    # The engine and model contract are shared with saved-posterior replay.
    # Score aggregation, conditioning and probe selection are implemented here.
    from .verifier import _model

    model = _model(bundle["context"]["model"])
    original_columns = model.columns
    training_ids = {
        tuple(identifier): index
        for index, identifier in enumerate(partition["character_ids"])
    }
    outcomes = [test_lookup[tuple(character_ids[index])] for index in indices]
    conditioning = None
    if partition["kind"] == "language":
        training = [original_columns[training_ids[tuple(character_ids[index])]] for index in indices]
        held_positions = {
            index for index, label in enumerate(bundle["tip_labels"])
            if label in partition["identifiers"]
        }
        outcomes = [
            [test[taxon] if taxon in held_positions else train[taxon]
             for taxon in range(len(train))]
            for train, test in zip(training, outcomes)
        ]
        conditioning = training
        model = replace(model, config=replace(model.config, ascertainment=False))

    samples = []
    for chain in bundle["chains"]:
        for index in np.linspace(0, len(chain) - 1, min(20, len(chain)), dtype=int):
            draw = chain[index]
            tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"])
            parameters = (
                tree, draw["gain"], draw["loss"], draw["clock_rate"],
            )
            loglikelihoods = model.character_scores(
                *parameters, outcomes, draw["borrow_rates"], draw.get("clock_multipliers"),
            )
            if conditioning is not None:
                normalizers = model.character_scores(
                    *parameters, conditioning, draw["borrow_rates"], draw.get("clock_multipliers"),
                )
                loglikelihoods = [a - b for a, b in zip(loglikelihoods, normalizers)]
            samples.append(loglikelihoods)
    per_character = np.array(samples).T
    scores = []
    for values in per_character:
        largest = max(values)
        scores.append(largest + math.log(sum(math.exp(v - largest) for v in values) / len(values)))
    return scores


def verify_prediction(directory: Path, dataset, probe_characters: int = 16) -> list[dict]:
    checks = []

    def record(name, ok, detail=None):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    report = _read(directory / "heldout.json")
    if report["study"] != dataset.name or report["is_synthetic"] != dataset.is_synthetic:
        record("prediction:study", False)
        return checks
    plan_path = directory / "PREDICTION_PLAN.json"
    plan = _read(plan_path) if plan_path.is_file() else None
    budget = report["budget"]
    budget_ok = (all(type(budget[key]) is int for key in ("chains", "warmup", "iterations", "thin"))
                 and 1 <= budget["chains"] <= 8 and budget["warmup"] >= 0
                 and budget["iterations"] >= budget["thin"] == 10)
    record("prediction:budget_shape", budget_ok)
    if not budget_ok:
        return checks
    if plan is not None:
        record("prediction:frozen_plan", plan["schema"] == 1 and plan["source_commit"] == report["source_commit"]
               and plan["seed"] == report["split_seed"] and plan["study"] == dataset.name
               and plan["models"] == list(_MODELS) and all(plan[key] == value for key, value in budget.items()))
    rng = np.random.default_rng(report["split_seed"])
    expected = {}
    for kind, choices in (
        ("concept", sorted(c.id for c in dataset.concepts)),
        ("language", sorted(language.id for language in dataset.languages)),
    ):
        expected[kind] = sorted(rng.choice(
            choices, size=max(1, len(choices) // 5), replace=False,
        ).tolist())

    for kind in ("concept", "language"):
        split = directory / kind
        partition_manifest = _read(split / "partition.json")
        ids = partition_manifest["test_character_ids"]
        record(f"prediction:{kind}:partition",
               partition_manifest["held_out_ids"] == expected[kind]
               and report["experiments"][kind]["held_out_ids"] == expected[kind])
        record(f"prediction:{kind}:character_count",
               len(ids) == report["experiments"][kind]["n_scored_characters"]
               and len({tuple(key) for key in ids}) == len(ids))
        baseline = _read(split / "tree_only.json")["log_scores"]
        if len(baseline) != len(ids):
            record(f"prediction:{kind}:baseline_count", False)
            continue
        partition = {
            "kind": kind, "identifiers": expected[kind],
            "character_ids": partition_manifest["train_character_ids"],
        }
        training = build_character_matrix(
            _masked_dataset(dataset, kind, set(expected[kind]), False), include_unobserved_languages=True,
        )
        training_ids = [(character.concept_id, character.cognate_class)
                        for character in training.variable_characters()]
        training_columns = [character.states for character in training.variable_characters()]
        record(f"prediction:{kind}:training_vocabulary",
               partition_manifest["train_character_ids"] == [list(key) for key in training_ids]
               and report["experiments"][kind]["n_train_characters"] == len(training_ids))
        lookup = _test_patterns(dataset, partition)
        if kind == "concept":
            scored_ids = list(lookup)
        else:
            positions = [index for index, label in enumerate(training.language_ids) if label in expected[kind]]
            scored_ids = [key for key in training_ids if any(lookup[key][index] >= 0 for index in positions)]
        exact_vocabulary = ids == [list(key) for key in scored_ids]
        record(f"prediction:{kind}:evaluation_vocabulary", exact_vocabulary)
        if not exact_vocabulary:
            continue
        probes = np.linspace(0, len(ids) - 1, min(probe_characters, len(ids)), dtype=int)
        for model_name in _MODELS:
            raw = _read(split / f"{model_name}.json")
            values = raw["log_scores"]
            valid = len(values) == len(ids) and all(math.isfinite(value) for value in values)
            record(f"prediction:{kind}:{model_name}:finite_scores", valid)
            if not valid:
                continue

            summary = report["experiments"][kind]["models"][model_name]
            by_concept = defaultdict(list)
            for (concept, _), value, base in zip(ids, values, baseline):
                by_concept[concept].append(value - base)
            deltas = np.array([sum(group) / len(group) for group in by_concept.values()])
            bootstrap_indices = rng.integers(0, len(deltas), size=(1000, len(deltas)))
            quantiles = np.quantile(deltas[bootstrap_indices].mean(axis=1), [0.025, 0.975])
            record(
                f"prediction:{kind}:{model_name}:summary",
                abs(sum(values) / len(values) - summary["mean_log_score_per_character"]) < 1e-10
                and abs(float(deltas.mean()) - summary["mean_concept_delta_vs_tree"]) < 1e-10
                and np.allclose(quantiles, summary["concept_bootstrap_interval_95"], atol=1e-10, rtol=0),
            )
            bundle = _read(split / model_name / "bundle.json")
            record(f"prediction:{kind}:{model_name}:source",
                   bundle["source_commit"] == report["source_commit"] == raw["source_commit"])
            stored = bundle["context"]
            record(f"prediction:{kind}:{model_name}:training_inputs",
                   bundle["tip_labels"] == training.language_ids
                   and stored["model"]["columns"] == training_columns
                   and stored["model"]["config"] == asdict(getattr(ModelConfig, model_name)())
                   and stored["held_out_partition"] == {**partition, "seed": report["split_seed"]})
            configurations = stored["chain_configurations"]
            expected_draws = budget["iterations"] // budget["thin"]
            record(f"prediction:{kind}:{model_name}:chain_budgets",
                   len(configurations) == bundle["n_chains"] == budget["chains"]
                   and len({configuration["seed"] for configuration in configurations}) == budget["chains"]
                   and all(all(configuration[key] == budget[key] for key in ("warmup", "iterations", "thin"))
                           for configuration in configurations)
                   and all(len(chain) == expected_draws and all(
                       draw["iteration"] == budget["warmup"] + budget["thin"] * (index + 1)
                       for index, draw in enumerate(chain)
                   ) for chain in bundle["chains"]))
            if budget["iterations"] >= 12000 and budget["iterations"] % budget["thin"] == 0:
                from .campaign import _checkpoint_matches

                record(f"prediction:{kind}:{model_name}:final_checkpoints", all(
                    _checkpoint_matches(_read(split / "checkpoints" / f"{model_name}_chain{index}.json"),
                                        bundle, index, report)
                    for index in range(budget["chains"])
                ))
            reproduced = _replay_scores(bundle, probes, ids, lookup, partition)
            record(f"prediction:{kind}:{model_name}:predictive_replay",
                   np.allclose(reproduced, np.array(values)[probes], atol=1e-9, rtol=0),
                   {"characters": len(probes), "draws_per_chain_at_most": 20})
    return checks
