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

    def rate_clock_ridge(
        self, gain: float, loss: float, clock: float, rng,
    ) -> tuple[tuple[float, float, float], MoveResult]:
        """Move along constant effective gain and loss with a log-symmetric step."""
        factor = math.exp(self.scale_lambda * (rng.uniform() - 0.5))
        values = (gain / factor, loss / factor, clock * factor)
        if any(not math.isfinite(value) or value <= 0 for value in values):
            return (gain, loss, clock), MoveResult(0.0, False)
        # In log coordinates the proposal translates by (-d,-d,d). The
        # product of the three original coordinates changes by exp(-d).
        return values, MoveResult(-math.log(factor), True)

    def root_age_scale(self, tree: PhyloTree, rng) -> MoveResult:
        """Scale only the root's positive gap above its oldest child."""
        node = tree.root
        lower = max(tree.heights[child] for child in tree.children()[node])
        gap = tree.heights[node] - lower
        if gap <= 0:
            return MoveResult(0.0, False)
        factor = math.exp(self.scale_lambda * (rng.uniform() - 0.5))
        proposed = lower + gap * factor
        if not math.isfinite(proposed) or proposed <= lower:
            return MoveResult(0.0, False)
        tree.heights[node] = proposed
        return MoveResult(math.log(factor), True)

    def time_clock_scale(self, tree: PhyloTree, clock: float, rng) -> tuple[float, MoveResult]:
        """Scale all internal ages and inversely scale the global clock."""
        old_root = tree.root_age()
        result = self.scale_root(tree, rng)
        if not result.valid:
            return clock, result
        factor = tree.root_age() / old_root
        proposed = clock / factor
        return proposed, MoveResult(result.log_hastings - math.log(factor), True)

    # --- node time slide ---------------------------------------------------
    def slide_node_time(self, tree: PhyloTree, node: int, rng) -> MoveResult:
        if tree.is_tip(node):
            return MoveResult(0.0, False)
        ch = tree.children()
        kids = ch[node]
        lower = max((tree.heights[c] for c in kids), default=0.0)
        parent = tree.parent[node]
        upper = tree.heights[parent] if parent >= 0 else math.inf
        span = upper - lower if parent >= 0 else max(lower, 1.0)
        window = self.time_window_frac * span
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

    def spr(self, tree: PhyloTree, rng, *, resample_height: bool = False) -> MoveResult:
        """Use the native attachment-age move when requested.

        The original time-preserving proposal draws uniformly from backbone
        edges crossing the fixed attachment age. Its reverse set is identical.
        The height-resampled kernel also moves the attachment age, using the
        ratio of the destination and original edge spans.
        """
        if resample_height:
            proposed = rng.spr_with_height(tree.parent, tree.heights, tree.taxon, tree.root)
            if proposed is None:
                return MoveResult(0.0, False)
            tree.parent, tree.heights, log_hastings = proposed
            return MoveResult(log_hastings, True)

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
