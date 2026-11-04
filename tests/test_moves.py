"""Tree-move integrity tests.

Topology moves must never corrupt the tree: no cycles, no orphaned nodes, a
stable node count, and valid node-time ordering after every accepted move. A
regression here once caused SPR to create a cycle that made likelihood traversal
blow up, so these tests run many moves and check structural invariants each time.
"""

from __future__ import annotations

import numpy as np

from diachron import _core
from diachron.models.tree import random_starting_tree
from diachron.infer.moves import MoveSet


def _has_no_cycle(tree) -> bool:
    for i in range(tree.n_nodes):
        node = i
        steps = 0
        seen = set()
        while tree.parent[node] >= 0:
            node = tree.parent[node]
            steps += 1
            if node in seen or steps > tree.n_nodes:
                return False
            seen.add(node)
    return True


def _exactly_one_root(tree) -> bool:
    return sum(1 for p in tree.parent if p < 0) == 1


def test_spr_preserves_tree_integrity():
    rng = np.random.default_rng(3)
    labels = [f"t{i}" for i in range(16)]
    tree = random_starting_tree(labels, 6000.0, rng)
    ms = MoveSet()
    pr = _core.PyRng(42)
    n0 = tree.n_nodes
    for _ in range(3000):
        snap_p, snap_h = list(tree.parent), list(tree.heights)
        res = ms.spr(tree, pr)
        assert _has_no_cycle(tree), "SPR created a cycle"
        assert _exactly_one_root(tree), "SPR produced not-exactly-one root"
        assert tree.n_nodes == n0, "SPR changed node count"
        if not res.valid:
            tree.parent, tree.heights = snap_p, snap_h
        else:
            assert tree.times_valid()


def test_nni_preserves_tree_integrity():
    rng = np.random.default_rng(5)
    labels = [f"t{i}" for i in range(12)]
    tree = random_starting_tree(labels, 6000.0, rng)
    ms = MoveSet()
    pr = _core.PyRng(7)
    n0 = tree.n_nodes
    for _ in range(3000):
        snap_p, snap_h = list(tree.parent), list(tree.heights)
        res = ms.nni(tree, pr)
        assert _has_no_cycle(tree)
        assert _exactly_one_root(tree)
        assert tree.n_nodes == n0
        if not res.valid:
            tree.parent, tree.heights = snap_p, snap_h


def test_node_slide_stays_in_window():
    rng = np.random.default_rng(9)
    labels = [f"t{i}" for i in range(10)]
    tree = random_starting_tree(labels, 6000.0, rng)
    ms = MoveSet()
    pr = _core.PyRng(11)
    for _ in range(2000):
        internals = tree.internal_nodes()
        node = internals[pr.below(len(internals))]
        snap = list(tree.heights)
        res = ms.slide_node_time(tree, node, pr)
        if res.valid:
            assert tree.times_valid()
        else:
            tree.heights = snap
