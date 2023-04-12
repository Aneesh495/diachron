"""Proposal moves for the MCMC, with their Hastings terms.

Each move takes the current state, proposes a change in place on a copy, and
returns the log Hastings ratio (log q(old|new) - log q(new|old) + log|Jacobian|)
plus a validity flag. The acceptance test in `chain.py` folds this into the log
posterior difference. Move correctness is checked against exact finite targets in
`tests/test_mcmc_exact.py`, which is the evidence for gate G03.

The moves use the project RNG (`diachron._core.PyRng`) so a chain's entire
stochastic trajectory is reproducible from its seed and resumable from a saved
RNG state.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from ..models.tree import PhyloTree


@dataclass
class MoveResult:
    log_hastings: float
    valid: bool


class MoveSet:
    """Container for step sizes, adapted during warmup then frozen."""

    def __init__(self):
        # Multiplicative window for scalar scale moves and additive window
        # (as a fraction of root age) for node-time slides.
        self.scale_lambda = 0.5
        self.time_window_frac = 0.1

    # --- scalar scale move -------------------------------------------------
    def scale_scalar(self, value: float, rng) -> tuple[float, MoveResult]:
        u = rng.uniform()
        factor = math.exp(self.scale_lambda * (u - 0.5))
        new = value * factor
        if not math.isfinite(new) or new <= 0:
            return value, MoveResult(0.0, False)
        # The multiplier proposal is symmetric in log space; the Hastings term is
        # the log Jacobian of the transform x -> x*factor, which is log(factor).
        return new, MoveResult(math.log(factor), True)

    # --- node time slide ---------------------------------------------------
    def slide_node_time(self, tree: PhyloTree, node: int, rng) -> MoveResult:
        if tree.is_tip(node):
            return MoveResult(0.0, False)
        ch = tree.children()
        kids = ch[node]
        lower = max((tree.heights[c] for c in kids), default=0.0)
        parent = tree.parent[node]
        upper = tree.heights[parent] if parent >= 0 else math.inf
        window = self.time_window_frac * max(tree.root_age(), 1.0)
        step = (rng.uniform() - 0.5) * 2.0 * window
        proposed = tree.heights[node] + step
        if proposed <= lower + 1e-9 or (upper != math.inf and proposed >= upper - 1e-9):
            return MoveResult(0.0, False)
        tree.heights[node] = proposed
        # Symmetric reflected window; Hastings term is zero.
        return MoveResult(0.0, True)

    # --- root height scale -------------------------------------------------
    def scale_root(self, tree: PhyloTree, rng) -> MoveResult:
        """Scale all internal node heights by a common factor.

        This moves the overall time scale without changing topology or relative
        times. Scaling k internal node ages by `factor` contributes a Jacobian of
        `factor^k` (k independent multiplied coordinates), so the Hastings term is
        `k * log(factor)`.
        """
        internals = tree.internal_nodes()
        u = rng.uniform()
        factor = math.exp(self.scale_lambda * (u - 0.5))
        for node in internals:
            tree.heights[node] *= factor
        if not tree.times_valid():
            # revert
            inv = 1.0 / factor
            for node in internals:
                tree.heights[node] *= inv
            return MoveResult(0.0, False)
        return MoveResult(len(internals) * math.log(factor), True)

    # --- nearest-neighbor interchange -------------------------------------
    def nni(self, tree: PhyloTree, rng) -> MoveResult:
        """Rooted NNI: swap a node's sibling with one of its children.

        Picks an internal child `c` of a node `p`, a sibling `sib` of `c` under
        `p`, and a child `gc` of `c`, then swaps `sib` and `gc`. The move is its
        own inverse under symmetric selection, so the Hastings term is zero. It
        is rejected if it would make a child older than its parent.
        """
        ch = tree.children()
        candidates = []
        for c in tree.internal_nodes():
            p = tree.parent[c]
            if p >= 0 and len(ch[p]) >= 2:
                candidates.append((p, c))
        if not candidates:
            return MoveResult(0.0, False)
        p, c = candidates[rng.below(len(candidates))]
        siblings = [x for x in ch[p] if x != c]
        if not siblings:
            return MoveResult(0.0, False)
        sib = siblings[rng.below(len(siblings))]
        kids_c = ch[c]
        if not kids_c:
            return MoveResult(0.0, False)
        gc = kids_c[rng.below(len(kids_c))]

        # Time validity after swap: c must be older than sib, p older than gc.
        if tree.heights[c] <= tree.heights[sib] + 1e-9:
            return MoveResult(0.0, False)
        if tree.heights[p] <= tree.heights[gc] + 1e-9:
            return MoveResult(0.0, False)

        # Rewire parent pointers.
        tree.parent[sib] = c
        tree.parent[gc] = p
        return MoveResult(0.0, True)

    # --- subtree prune and regraft (bounded) ------------------------------
    def spr(self, tree: PhyloTree, rng) -> MoveResult:
        """A bounded subtree-prune-regraft: move a subtree to a nearby edge.

        SPR explores topology more broadly than NNI. We prune a non-root subtree
        and regraft it onto an edge whose time interval admits the subtree's
        crown age, keeping the move reversible by restricting regraft targets to
        edges that span the pruned node's age. The proposal selects uniformly
        among valid targets; the Hastings ratio is the ratio of the number of
        valid targets before and after, computed from the tree structure.
        """
        # Choose a prunable node (not root, not a child of root with a single
        # sibling only, to keep a valid binary structure after removal).
        internal_and_tips = [i for i in range(tree.n_nodes) if tree.parent[i] >= 0]
        if len(internal_and_tips) < 3:
            return MoveResult(0.0, False)
        node = internal_and_tips[rng.below(len(internal_and_tips))]
        old_parent = tree.parent[node]
        if old_parent < 0:
            return MoveResult(0.0, False)
        grandparent = tree.parent[old_parent]

        ch = tree.children()
        # old_parent's other children become attached to grandparent.
        siblings = [x for x in ch[old_parent] if x != node]
        if grandparent < 0 or len(siblings) != 1:
            # Keep this move to the clean binary case for a correct Jacobian.
            return MoveResult(0.0, False)
        sib = siblings[0]

        node_age = tree.heights[node]
        # Candidate target edges: (child) whose branch spans node_age, excluding
        # the pruned subtree itself and the collapsing edge.
        subtree = _descendants(tree, node) | {node}
        targets = []
        for child in range(tree.n_nodes):
            p = tree.parent[child]
            if p < 0 or child in subtree or child == old_parent:
                continue
            lo = tree.heights[child]
            hi = tree.heights[p]
            if lo < node_age < hi:
                targets.append(child)
        if not targets:
            return MoveResult(0.0, False)
        n_targets_before = len(targets)

        target = targets[rng.below(len(targets))]
        target_parent = tree.parent[target]

        # Detach: connect sib directly to grandparent, freeing old_parent.
        tree.parent[sib] = grandparent
        # Reuse old_parent as the new internal node on the target edge.
        new_internal = old_parent
        tree.heights[new_internal] = node_age  # place it at the subtree's age
        tree.parent[new_internal] = target_parent
        tree.parent[target] = new_internal
        tree.parent[node] = new_internal

        if not tree.times_valid():
            return MoveResult(0.0, False)

        # Count reverse targets for the Hastings ratio.
        subtree_after = _descendants(tree, node) | {node}
        targets_after = []
        for child in range(tree.n_nodes):
            p = tree.parent[child]
            if p < 0 or child in subtree_after or child == new_internal:
                continue
            lo = tree.heights[child]
            hi = tree.heights[p]
            if lo < node_age < hi:
                targets_after.append(child)
        n_targets_after = max(len(targets_after), 1)
        return MoveResult(math.log(n_targets_before / n_targets_after), True)


def _descendants(tree: PhyloTree, node: int) -> set[int]:
    ch = tree.children()
    out: set[int] = set()
    stack = list(ch[node])
    while stack:
        n = stack.pop()
        out.add(n)
        stack.extend(ch[n])
    return out
