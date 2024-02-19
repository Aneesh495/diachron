"""Held-out region prediction on known trees under two spatial truth processes."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

from ..models.tree import PhyloTree, random_starting_tree
from ..util.gitmeta import source_commit
from .marginals import SpatialKernel, state_marginals

_AXES = np.concatenate([np.eye(3), -np.eye(3)])
_SCALES = (.25, .5, 1., 2., 4.)
_RATE = 0.00016
_TIP_COUNT = 12


def fibonacci_sphere(size: int) -> np.ndarray:
    indices = np.arange(size)
    height = 1 - 2 * (indices + .5) / size
    phi = np.pi * (3 - np.sqrt(5)) * indices
    radius = np.sqrt(1 - height ** 2)
    return np.column_stack([radius * np.cos(phi), radius * np.sin(phi), height])


def _generator(points: np.ndarray, local: bool) -> np.ndarray:
    dots = np.clip(points @ points.T, -1, 1)
    angles = np.arccos(dots)
    if local:
        weights = np.exp(-angles ** 2 / (2 * .34 ** 2)) * (angles <= .65)
    else:
        weights = 1 / np.maximum(angles, .01)
    np.fill_diagonal(weights, 0.)
    # Match average instantaneous loss of directional correlation per year.
    weights *= _RATE / np.mean(np.sum(weights * (1 - dots), axis=1))
    return weights - np.diag(weights.sum(axis=1))


def _geodesic_step(point: np.ndarray, variance: float, rng: np.random.Generator) -> np.ndarray:
    tangent = rng.normal(size=3)
    tangent -= np.dot(tangent, point) * point
    tangent *= np.sqrt(variance)
    length = np.linalg.norm(tangent)
    moved = np.cos(length) * point + np.sinc(length / np.pi) * tangent
    return moved / np.linalg.norm(moved)


def _simulate_sphere(point: np.ndarray, years: float, rng) -> np.ndarray:
    steps = max(1, math.ceil(_RATE * years / .0025))
    for _ in range(steps):
        point = _geodesic_step(point, _RATE * years / steps, rng)
    return point


def _simulate_regions(state: int, years: float, generator: np.ndarray, rng) -> int:
    elapsed = 0.
    while True:
        rate = -generator[state, state]
        elapsed += rng.exponential(1 / rate)
        if elapsed >= years:
            return state
        weights = generator[state].copy()
        weights[state] = 0.
        state = int(rng.choice(len(generator), p=weights / rate))


@dataclass(frozen=True)
class _SpatialCase:
    tree: PhyloTree
    root_region: int
    tip_regions: list[int]
    withheld: list[int]


def _case(seed: int, truth: str, region_generator: np.ndarray) -> _SpatialCase:
    rng = np.random.default_rng(seed)
    tree = random_starting_tree([str(taxon) for taxon in range(_TIP_COUNT)], 4500., rng)
    # Freeze which tips will be hidden before generating geography.
    withheld = sorted(rng.choice(_TIP_COUNT, size=3, replace=False).tolist())
    children = tree.children()
    order = [tree.root]
    regions = [0] * tree.n_nodes
    points = np.zeros((tree.n_nodes, 3))
    if truth == "sphere":
        root = rng.normal(size=3)
        points[tree.root] = root / np.linalg.norm(root)
    else:
        regions[tree.root] = int(rng.integers(6))
    durations = tree.branch_lengths()
    for node in order:
        order.extend(children[node])
        parent = tree.parent[node]
        if truth == "sphere":
            if parent >= 0:
                points[node] = _simulate_sphere(points[parent], durations[node], rng)
            regions[node] = int(np.argmax(_AXES @ points[node]))
        elif parent >= 0:
            regions[node] = _simulate_regions(regions[parent], durations[node], region_generator, rng)
    tip_regions = [regions[tree.taxon.index(taxon)] for taxon in range(_TIP_COUNT)]
    return _SpatialCase(tree, regions[tree.root], tip_regions, withheld)


def _root_interval(probabilities, truth):
    ordered = np.argsort(-probabilities, kind="stable")
    size = int(np.searchsorted(probabilities[ordered].cumsum(), .9)) + 1
    return {"size": size, "covered": bool(truth in ordered[:size])}


def _predict(case, kernel, mapping, scale):
    masks = np.array([mapping == region for region in range(6)], dtype=float)
    evidence = {taxon: masks[region] for taxon, region in enumerate(case.tip_regions) if taxon not in case.withheld}
    marginals = state_marginals(case.tree, kernel, evidence, scale).probabilities
    regional = marginals @ masks.T
    root = regional[case.tree.root]
    predictions = [regional[case.tree.taxon.index(taxon)].tolist() for taxon in case.withheld]
    return {"root": root.tolist(), "root_90set": _root_interval(root, case.root_region), "heldout": predictions}


def _score_rows(rows):
    statistics = {}
    for model in ("regions", "sphere_mesh"):
        units = []
        covered = sizes = 0
        for row in rows:
            prediction = row["predictions"][model]
            probabilities = [probabilities[row["tip_regions"][tip]]
                             for probabilities, tip in zip(prediction["heldout"], row["withheld"])]
            units.append(-float(np.mean(np.log(probabilities))))
            covered += prediction["root_90set"]["covered"]
            sizes += prediction["root_90set"]["size"]
        statistics[model] = {"heldout_nll": float(np.mean(units)), "unit_nll": units,
                             "root_90set_coverage": covered / len(rows), "root_90set_size": sizes / len(rows)}
    paired = np.subtract(statistics["sphere_mesh"].pop("unit_nll"), statistics["regions"].pop("unit_nll"))
    rng = np.random.default_rng(65327)
    boot = paired[rng.integers(0, len(paired), (1000, len(paired)))].mean(axis=1)
    statistics["paired_nll_mesh_minus_regions"] = {"mean": float(paired.mean()),
                                                   "bootstrap_95ci": np.quantile(boot, [.025, .975]).tolist()}
    return statistics


def run_spatial_study(out: Path, train_count: int = 40, test_count: int = 120) -> dict:
    if not 1 <= train_count <= 40 or not 1 <= test_count <= 120:
        raise ValueError("invalid spatial study budget")
    out.mkdir(parents=True)
    region_generator = _generator(_AXES, local=False)
    mesh = fibonacci_sphere(96)
    models = {"regions": (SpatialKernel(region_generator), np.arange(6)),
              "sphere_mesh": (SpatialKernel(_generator(mesh, local=True)), np.argmax(mesh @ _AXES.T, axis=1))}
    results = {}
    with (out / "cases.jsonl").open("x") as stream:
        for truth_index, truth in enumerate(("sphere", "regions")):
            base_seed = 64000 + 1000 * truth_index
            training = [_case(base_seed + i, truth, region_generator) for i in range(train_count)]
            scales = {}
            training_scores = {}
            for name, (kernel, mapping) in models.items():
                masks = np.array([mapping == region for region in range(6)], dtype=float)
                scores = [sum(state_marginals(case.tree, kernel,
                                               {taxon: masks[region] for taxon, region in enumerate(case.tip_regions)},
                                               scale).log_likelihood for case in training) for scale in _SCALES]
                scales[name] = _SCALES[int(np.argmax(scores))]
                training_scores[name] = dict(zip(_SCALES, scores))
            rows = []
            for i in range(test_count):
                seed = base_seed + 200 + i
                case = _case(seed, truth, region_generator)
                row = {"truth": truth, "seed": seed, "root_region": case.root_region, "tip_regions": case.tip_regions,
                       "withheld": case.withheld, "parent": case.tree.parent, "heights": case.tree.heights,
                       "predictions": {name: _predict(case, kernel, mapping, scales[name])
                                       for name, (kernel, mapping) in models.items()}}
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                rows.append(row)
            results[truth] = {"train_seed_start": base_seed, "test_seed_start": base_seed + 200,
                              "selected_rate_scales": scales, "training_log_likelihood": training_scores,
                              "scores": _score_rows(rows)}
    summary = {"schema": 1, "is_synthetic": True, "source_commit": source_commit(), "train_per_truth": train_count,
               "test_per_truth": test_count, "regions": 6, "mesh_states": 96, "test_tips_per_tree": 3,
               "results": results}
    (out / "spatial-study.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary
