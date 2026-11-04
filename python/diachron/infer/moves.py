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

        SPR explores topology more broadly than NNI. We prune a node `v` together
        with its parent edge. Its parent `u` has exactly one other child `sib` in
        the binary case; detaching `v` collapses `u`, so `sib` reconnects directly
        to `u`'s parent `gp`, and `u` is freed for reuse. We then regraft by
        reinserting `u` on a target edge `(w, parent_w)` whose time interval spans
        `v`'s crown age, making `u` the new parent of both `v` and `w`.

        The move is kept reversible and safe by three guards: the target edge is
        chosen from edges that span `v`'s age and that are not inside the pruned
        subtree (so no cycle can form), `sib` is never chosen as the target (that
        would be the identity after collapse), and the resulting times are checked.
        The Hastings ratio is the ratio of the number of valid targets before and
        after the move, computed from the tree structure.

        Correctness note: an earlier version reattached `u` without excluding the
        collapsed edge and the pruned subtree, which could create a cycle. The
        exclusions below are what prevent that; `tests/test_moves.py` checks that
        many SPR moves never corrupt the tree.
        """
        ch = tree.children()

        # Prunable nodes: any non-root node whose parent is binary (exactly two
        # children) so collapsing the parent yields a valid binary tree.
        prunable = []
        for v in range(tree.n_nodes):
            u = tree.parent[v]
            if u < 0:
                continue
            gp = tree.parent[u]
            if gp < 0:
                continue  # parent is the root; skip to keep the root binary
            if len(ch[u]) != 2:
                continue
            prunable.append(v)
        if not prunable:
            return MoveResult(0.0, False)

        v = prunable[rng.below(len(prunable))]
        u = tree.parent[v]
        gp = tree.parent[u]
        sib = next(x for x in ch[u] if x != v)

        node_age = tree.heights[v]
        subtree = _descendants(tree, v) | {v, u}

        # Valid target edges: child-node `w` whose edge (w -> parent_w) strictly
        # spans node_age, excluding anything in the pruned subtree, the collapsed
        # node `u`, and the sibling edge (reattaching above sib is the identity).
        def valid_targets() -> list[int]:
            out = []
            for w in range(tree.n_nodes):
                pw = tree.parent[w]
                if pw < 0 or w in subtree or w == u or w == sib:
                    continue
                if tree.heights[w] < node_age < tree.heights[pw]:
                    out.append(w)
            return out

        targets_before = valid_targets()
        if not targets_before:
            return MoveResult(0.0, False)
        n_before = len(targets_before)

        w = targets_before[rng.below(n_before)]
        parent_w = tree.parent[w]

        # Prune: reconnect sib to gp, freeing u.
        tree.parent[sib] = gp
        # Regraft: u becomes the parent of w and v, inserted on edge (parent_w, w).
        tree.parent[u] = parent_w
        tree.parent[w] = u
        tree.parent[v] = u
        # u sits at the subtree's crown age on the target edge.
        tree.heights[u] = node_age

        if not tree.times_valid():
            return MoveResult(0.0, False)

        # Reverse-move target count for the Hastings ratio. After the move, the
        # pruned subtree is {v, u} again; recompute from the new structure.
        ch2 = tree.children()
        sib_after = sib  # sib now hangs under gp
        subtree_after = _descendants(tree, v) | {v, u}

        def valid_targets_after() -> int:
            count = 0
            for x in range(tree.n_nodes):
                px = tree.parent[x]
                if px < 0 or x in subtree_after or x == u or x == sib_after:
                    continue
                if tree.heights[x] < node_age < tree.heights[px]:
                    count += 1
            return count

        n_after = max(valid_targets_after(), 1)
        _ = ch2
        return MoveResult(math.log(n_before / n_after), True)


def _descendants(tree: PhyloTree, node: int) -> set[int]:
    ch = tree.children()
    out: set[int] = set()
    stack = list(ch[node])
    while stack:
        n = stack.pop()
        out.add(n)
        stack.extend(ch[n])
    return out
