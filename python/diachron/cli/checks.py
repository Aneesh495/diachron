"""Runnable correctness and performance checks behind the CLI.

`run_likelihood_oracle` is the command-line face of gate G02: it compares the
Rust pruning likelihood against the independent enumeration on many random tiny
trees and reports the worst disagreement. `run_bench` and `run_stress` measure
likelihood throughput and exercise large trees to confirm the scaling stays
numerically stable.
"""

from __future__ import annotations

import math
import time

import numpy as np

from diachron import _core
from ..models.tree import random_starting_tree


def _random_tiny_tree(rng: np.random.Generator):
    n_tips = int(rng.integers(2, 6))
    parent = [-1] * (2 * n_tips - 1)
    heights = [0.0] * (2 * n_tips - 1)
    taxon = [-1] * (2 * n_tips - 1)
    for t in range(n_tips):
        taxon[t] = t
    active = list(range(n_tips))
    nxt = n_tips
    age = 0.0
    while len(active) > 1:
        age += float(rng.uniform(0.1, 1.0))
        a = active.pop(int(rng.integers(0, len(active))))
        b = active.pop(int(rng.integers(0, len(active))))
        node = nxt
        nxt += 1
        parent[a] = node
        parent[b] = node
        heights[node] = age
        active.append(node)
    root = active[0]
    bl = [0.0] * len(parent)
    for i, p in enumerate(parent):
        if p >= 0:
            bl[i] = heights[p] - heights[i]
    return parent, bl, heights, taxon, root, n_tips


def _enumerate_ll(parent, bl, taxon, root, gain, loss, column, root_prior):
    import itertools

    n_nodes = len(parent)

    def trans(t):
        z = gain + loss
        e = math.exp(-z * t)
        pi0, pi1 = loss / z, gain / z
        return [[pi0 + pi1 * e, pi1 - pi1 * e], [pi0 - pi0 * e, pi1 + pi0 * e]]

    tmats = [trans(bl[i]) for i in range(n_nodes)]
    internal = [i for i in range(n_nodes) if taxon[i] < 0]
    missing = [i for i in range(n_nodes) if taxon[i] >= 0 and column[taxon[i]] < 0]
    free = internal + missing
    total = 0.0
    for assignment in itertools.product(range(2), repeat=len(free)):
        state = [-1] * n_nodes
        for idx, node in enumerate(free):
            state[node] = assignment[idx]
        for i in range(n_nodes):
            if taxon[i] >= 0 and column[taxon[i]] >= 0:
                state[i] = column[taxon[i]]
        term = root_prior[state[root]]
        for i in range(n_nodes):
            p = parent[i]
            if p >= 0:
                term *= tmats[i][state[p]][state[i]]
        total += term
    return math.log(total) if total > 0 else float("-inf")


def run_likelihood_oracle(n_cases: int = 1000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    worst = 0.0
    checked = 0
    t0 = time.time()
    for _ in range(n_cases):
        parent, bl, heights, taxon, root, n_tips = _random_tiny_tree(rng)
        gain = float(rng.uniform(0.1, 2.0))
        loss = float(rng.uniform(0.1, 2.0))
        z = gain + loss
        prior = [loss / z, gain / z]
        column = [
            -1 if rng.random() < 0.15 else int(rng.integers(0, 2))
            for _ in range(n_tips)
        ]
        fast = _core.char_log_likelihood(
            parent, bl, heights, taxon, root, "binary", [gain, loss], column, prior
        )
        ref = _enumerate_ll(parent, bl, taxon, root, gain, loss, column, prior)
        if math.isinf(ref) and math.isinf(fast):
            continue
        worst = max(worst, abs(fast - ref))
        checked += 1
    return {
        "gate": "G02",
        "cases": n_cases,
        "checked": checked,
        "worst_abs_diff": worst,
        "tolerance": 1e-9,
        "passed": worst < 1e-9,
        "seconds": round(time.time() - t0, 2),
    }


def run_bench() -> dict:
    """Likelihood throughput on a moderately sized tree."""
    rng = np.random.default_rng(42)
    labels = [f"t{i}" for i in range(60)]
    tree = random_starting_tree(labels, 6000.0, rng)
    parent, bl, h, tax, root = tree.to_core_arrays()
    gain, loss = 0.0005, 0.001
    z = gain + loss
    prior = [loss / z, gain / z]
    columns = [
        [int(rng.integers(0, 2)) for _ in range(60)] for _ in range(500)
    ]
    t0 = time.time()
    reps = 20
    for _ in range(reps):
        _core.block_log_likelihood(
            parent, bl, h, tax, root, "binary", [gain, loss], columns, prior
        )
    dt = time.time() - t0
    evals = reps * len(columns)
    return {
        "tree_tips": 60,
        "characters": len(columns),
        "block_evals": reps,
        "seconds": round(dt, 3),
        "char_likelihoods_per_sec": round(evals / dt),
    }


def run_stress() -> dict:
    """Large-tree numerical stability: likelihood stays finite as depth grows."""
    rng = np.random.default_rng(7)
    results = []
    for n in [50, 200, 500]:
        labels = [f"t{i}" for i in range(n)]
        tree = random_starting_tree(labels, 6000.0, rng)
        parent, bl, h, tax, root = tree.to_core_arrays()
        gain, loss = 0.0005, 0.001
        z = gain + loss
        prior = [loss / z, gain / z]
        column = [int(rng.integers(0, 2)) for _ in range(n)]
        ll = _core.char_log_likelihood(
            parent, bl, h, tax, root, "binary", [gain, loss], column, prior
        )
        results.append({"tips": n, "log_likelihood": ll, "finite": math.isfinite(ll)})
    return {"stress": results, "all_finite": all(r["finite"] for r in results)}
