"""Gate G02: the Rust pruning likelihood matches an independent enumerator.

We generate many tiny random trees and character patterns, compute the log
likelihood both with the production Rust binding and with the pure-Python
enumeration oracle, and require agreement to a tight tolerance. The count exceeds
the 1000-case floor the acceptance profile sets.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from diachron import _core
from oracle import enumerate_log_likelihood


def _random_tiny_tree(rng: np.random.Generator):
    """Build a random rooted tree with 2..5 tips, returning core arrays."""
    n_tips = int(rng.integers(2, 6))
    # Coalescent-style merge to make a valid rooted tree with ordered times.
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
    branch_len = [0.0] * len(parent)
    for i, p in enumerate(parent):
        if p >= 0:
            branch_len[i] = heights[p] - heights[i]
    return parent, branch_len, heights, taxon, root, n_tips


@pytest.mark.parametrize("seed", range(1200))
def test_pruning_matches_enumeration(seed):
    rng = np.random.default_rng(seed)
    parent, branch_len, heights, taxon, root, n_tips = _random_tiny_tree(rng)
    gain = float(rng.uniform(0.1, 2.0))
    loss = float(rng.uniform(0.1, 2.0))
    z = gain + loss
    root_prior = [loss / z, gain / z]

    # Random pattern with occasional missing states.
    column = []
    for _ in range(n_tips):
        if rng.random() < 0.15:
            column.append(-1)  # missing
        else:
            column.append(int(rng.integers(0, 2)))

    fast = _core.char_log_likelihood(
        parent, branch_len, heights, taxon, root, "binary", [gain, loss], column, root_prior
    )
    ref = enumerate_log_likelihood(
        parent, branch_len, taxon, root, gain, loss, column, root_prior
    )

    if math.isinf(ref) and math.isinf(fast):
        return
    assert abs(fast - ref) < 1e-9, f"seed {seed}: fast={fast} ref={ref}"
