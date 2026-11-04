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
        window = self.time_window_frac * max(lower if parent < 0 else tree.root_age(), 1.0)
        step = (rng.uniform() - 0.5) * 2.0 * window
        proposed = tree.heights[node] + step
        if proposed <= lower + 1e-9 or (upper != math.inf and proposed >= upper - 1e-9):
            return MoveResult(0.0, False)
        tree.heights[node] = proposed
        # The window is invariant under this move; out-of-bounds steps reject.
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
        scaled = [tree.heights[i] * factor if not tree.is_tip(i) else tree.heights[i]
                  for i in range(tree.n_nodes)]
        if any(not math.isfinite(age) for age in scaled) or any(
            p >= 0 and scaled[p] < scaled[i] for i, p in enumerate(tree.parent)
        ):
            return MoveResult(0.0, False)
        tree.heights = scaled
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

    def spr(self, tree: PhyloTree, rng) -> MoveResult:
        """Time-preserving SPR with symmetric selection on the detached backbone.

        Choose an ordered prunable child and preserve its parent's age. Suppress
        the parent and draw uniformly from backbone edges crossing that age.
        The original sibling edge is in the same set and is a self-transition.
        Given the same child on reversal, both the backbone and age are
        identical, so the proposal probabilities cancel. Binary node degrees
        and the number of children beneath the root are invariant.
        """
        children = tree.children()
        prunable = [
            child for child, parent in enumerate(tree.parent)
            if parent >= 0 and tree.parent[parent] >= 0 and len(children[parent]) == 2
        ]
        if not prunable:
            return MoveResult(0.0, False)

        child = prunable[rng.below(len(prunable))]
        parent = tree.parent[child]
        grandparent = tree.parent[parent]
        sibling = next(node for node in children[parent] if node != child)
        age = tree.heights[parent]
        if tree.heights[child] >= age or not tree.heights[sibling] < age < tree.heights[grandparent]:
            return MoveResult(0.0, False)

        detached = _descendants(tree, child) | {child, parent}
        targets = []
        for node, ancestor in enumerate(tree.parent):
            if node in detached or ancestor < 0:
                continue
            backbone_parent = grandparent if node == sibling else ancestor
            if tree.heights[node] < age < tree.heights[backbone_parent]:
                targets.append((node, backbone_parent))
        target, target_parent = targets[rng.below(len(targets))]
        if target == sibling:
            return MoveResult(0.0, False)

        tree.parent[sibling] = grandparent
        tree.parent[parent] = target_parent
        tree.parent[target] = parent
        return MoveResult(0.0, True)


def _descendants(tree: PhyloTree, node: int) -> set[int]:
    ch = tree.children()
    out: set[int] = set()
    stack = list(ch[node])
    while stack:
        n = stack.pop()
        out.add(n)
        stack.extend(ch[n])
    return out
