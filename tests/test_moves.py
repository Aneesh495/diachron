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


def test_spr_returns_to_exact_original_state_at_same_parent_age():
    class Choice:
        def __init__(self, *choices):
            self.choices = iter(choices)

        def below(self, n):
            choice = next(self.choices)
            assert 0 <= choice < n
            return choice

    from diachron.models.tree import PhyloTree

    tree = PhyloTree([4, 4, 5, 5, 6, 6, -1], [0, 0, 0, 0, 2, 4, 8], [0, 1, 2, 3, -1, -1, -1])
    original = tree.copy()
    moves = MoveSet()
    # Prune child 0 at age 2 and insert its parent above tip 2. The
    # detached backbone has target edges above tips 1, 2, and 3.
    assert moves.spr(tree, Choice(0, 1)).valid
    assert tree.parent != original.parent
    assert tree.heights == original.heights
    assert moves.spr(tree, Choice(0, 0)).valid
    assert tree == original


def test_root_slide_uses_same_window_in_both_directions():
    class Step:
        def __init__(self, value):
            self.value = value

        def uniform(self):
            return self.value

    from diachron.models.tree import PhyloTree

    tree = PhyloTree([2, 2, 4, 4, -1], [0, 0, 100, 0, 200], [0, 1, -1, 2, -1])
    moves = MoveSet()
    assert moves.slide_node_time(tree, 4, Step(0.75)).valid
    assert moves.slide_node_time(tree, 4, Step(0.25)).valid
    assert tree.root_age() == 200


def test_local_time_window_remains_useful_between_close_events():
    class HalfStep:
        def uniform(self):
            return 0.75

    from diachron.models.tree import PhyloTree

    tree = PhyloTree([3, 3, 4, 4, -1], [999.0, 998.0, 0.0, 999.5, 1000.0], [0, 1, 2, -1, -1])
    result = MoveSet().slide_node_time(tree, 3, HalfStep())
    assert result.valid
    assert tree.heights[3] == 999.55
