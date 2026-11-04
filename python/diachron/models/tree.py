"""Python-side timed phylogenetic tree.

This mirrors the Rust `Tree` but lives in Python because the MCMC driver, the
priors, and the serialization all need to read and mutate tree structure. The
flat-array form (parent/heights/taxon) is what crosses to the Rust likelihood;
this class maintains it and offers the queries the sampler needs.

Heights are node ages in years before present. Tips are usually at age 0, but a
dated or extinct tip can sit at a positive age. Branch length above a node is
its parent's height minus its own; the tree enforces that children are younger
than parents so branch lengths are nonnegative.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import copy

import numpy as np


@dataclass
class PhyloTree:
    # parent[i] is the parent node index, or -1 for the root.
    parent: list[int]
    # heights[i] is node age in years before present.
    heights: list[float]
    # taxon[i] is the taxon index for a tip, or -1 for an internal node.
    taxon: list[int]
    # Optional labels for tips, keyed by taxon index.
    tip_labels: dict[int, str] = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return len(self.parent)

    @property
    def root(self) -> int:
        for i, p in enumerate(self.parent):
            if p < 0:
                return i
        raise ValueError("tree has no root")

    def children(self) -> dict[int, list[int]]:
        ch: dict[int, list[int]] = {i: [] for i in range(self.n_nodes)}
        for i, p in enumerate(self.parent):
            if p >= 0:
                ch[p].append(i)
        return ch

    def is_tip(self, i: int) -> bool:
        return self.taxon[i] >= 0

    def internal_nodes(self) -> list[int]:
        return [i for i in range(self.n_nodes) if not self.is_tip(i)]

    def tip_nodes(self) -> list[int]:
        return [i for i in range(self.n_nodes) if self.is_tip(i)]

    def branch_lengths(self) -> list[float]:
        """Edge length above each node; the root's entry is 0."""
        bl = [0.0] * self.n_nodes
        for i, p in enumerate(self.parent):
            if p >= 0:
                bl[i] = self.heights[p] - self.heights[i]
        return bl

    def times_valid(self) -> bool:
        for i, p in enumerate(self.parent):
            if p >= 0 and self.heights[p] < self.heights[i] - 1e-9:
                return False
        return True

    def root_age(self) -> float:
        return self.heights[self.root]

    def copy(self) -> "PhyloTree":
        return PhyloTree(
            parent=list(self.parent),
            heights=list(self.heights),
            taxon=list(self.taxon),
            tip_labels=dict(self.tip_labels),
        )

    def to_core_arrays(self) -> tuple[list[int], list[float], list[float], list[int], int]:
        """Return (parent, branch_len, height, taxon, root) for the Rust binding."""
        return (
            list(self.parent),
            self.branch_lengths(),
            list(self.heights),
            list(self.taxon),
            self.root,
        )

    def newick(self) -> str:
        """Serialize to a Newick string with branch lengths in years."""
        ch = self.children()
        bl = self.branch_lengths()

        def render(node: int) -> str:
            if self.is_tip(node):
                label = self.tip_labels.get(self.taxon[node], f"t{self.taxon[node]}")
                return f"{label}:{bl[node]:.4f}"
            inner = ",".join(render(c) for c in ch[node])
            length = bl[node]
            return f"({inner}):{length:.4f}"

        return render(self.root) + ";"


def tree_from_truth(truth) -> PhyloTree:
    """Build a PhyloTree from a simulator TrueHistory."""
    labels = {t: lid for t, lid in enumerate(truth.tip_language_ids)}
    return PhyloTree(
        parent=list(truth.parent),
        heights=list(truth.heights),
        taxon=list(truth.taxon),
        tip_labels=labels,
    )


def random_starting_tree(
    tip_labels: list[str], root_age: float, rng: np.random.Generator,
    tip_ages: dict[int, float] | None = None,
) -> PhyloTree:
    """A random ultrametric starting tree for a chain, independent of the truth.

    Built by the same coalescent merge used in simulation but with the inference
    chain's own RNG, so chains start from genuinely different topologies.
    """
    n = len(tip_labels)
    if n < 2:
        raise ValueError("a study requires at least two tips")
    tip_ages = tip_ages or {}
    oldest = max(tip_ages.values(), default=0.0)
    if oldest >= root_age:
        raise ValueError("starting root age must exceed all dated tips")
    parent = [-1] * (2 * n - 1)
    heights = [0.0] * (2 * n - 1)
    taxon = [-1] * (2 * n - 1)
    for t in range(n):
        taxon[t] = t
    active = list(range(n))
    next_node = n
    age = 0.0
    k = n
    while len(active) > 1:
        rate = k * (k - 1) / 2.0
        age += rng.exponential(1.0 / rate) if rate > 0 else root_age
        i = int(rng.integers(0, len(active)))
        a = active.pop(i)
        j = int(rng.integers(0, len(active)))
        b = active.pop(j)
        node = next_node
        next_node += 1
        parent[a] = node
        parent[b] = node
        heights[node] = age
        active.append(node)
        k -= 1
    root = active[0]
    scale = (root_age - oldest) / heights[root]
    heights = [oldest + h * scale for h in heights]
    for tip in range(n):
        heights[tip] = tip_ages.get(tip, 0.0)
    labels = {t: tip_labels[t] for t in range(n)}
    return PhyloTree(parent=parent, heights=heights, taxon=taxon, tip_labels=labels)
