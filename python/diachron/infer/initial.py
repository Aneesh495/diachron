"""Construct independent starting trees from the permitted lexical observations."""

from __future__ import annotations

import numpy as np

from ..models.tree import PhyloTree


def lexical_distances(columns: list[list[int]]) -> np.ndarray:
    """Observed-mask Jaccard distances, used only to initialize the chain."""
    data = np.asarray(columns, dtype=np.int32)
    if data.ndim != 2 or not data.shape[0] or data.shape[1] < 2:
        raise ValueError("distance initialization requires a nonempty character matrix")
    observed = (data >= 0).astype(np.int32)
    present = (data == 1).astype(np.int32)
    shared = present.T @ present
    positive_observed = present.T @ observed
    union = positive_observed + positive_observed.T - shared
    distance = np.ones(union.shape, dtype=float)
    np.divide(shared, union, out=distance, where=union > 0)
    distance[union > 0] = 1.0 - distance[union > 0]
    np.fill_diagonal(distance, 0.0)
    return distance


def distance_starting_tree(
    labels: list[str], columns: list[list[int]], root_age: float,
    rng: np.random.Generator, tip_ages: dict[int, float],
) -> PhyloTree:
    """Average-linkage seed from a bootstrap of the supplied training columns.

    Distances set only the initial arrangement. The posterior still evaluates
    every state with the production likelihood and the declared priors.
    """
    if not columns or len(labels) != len(columns[0]):
        raise ValueError("initialization labels disagree with character width")
    n = len(labels)
    if not 2 <= n <= 1000:
        raise ValueError("distance initialization supports two to 1000 taxa")
    if not np.isfinite(root_age) or root_age <= max(tip_ages.values(), default=0.0):
        raise ValueError("initial root must predate every sampled tip")
    bootstrap = rng.integers(0, len(columns), size=len(columns))
    distance = lexical_distances([columns[index] for index in bootstrap])
    # Tiny independent tie-breaking noise preserves the order of distinct distances.
    jitter = rng.uniform(0.0, 1e-10, (n, n))
    distance += (jitter + jitter.T) / 2
    np.fill_diagonal(distance, np.inf)

    parents = [-1] * (2 * n - 1)
    taxa = list(range(n)) + [-1] * (n - 1)
    cluster_node = list(range(n))
    sizes = np.ones(n, dtype=np.int64)
    raw = [0.0] * (2 * n - 1)
    oldest = [tip_ages.get(taxon, 0.0) for taxon in range(n)] + [0.0] * (n - 1)
    active = np.ones(n, dtype=bool)
    for node in range(n, 2 * n - 1):
        a, b = np.unravel_index(np.argmin(distance), distance.shape)
        if not np.isfinite(distance[a, b]):
            raise ValueError("initial clustering exhausted finite distances")
        left, right = cluster_node[a], cluster_node[b]
        parents[left] = parents[right] = node
        raw[node] = max(distance[a, b] / 2, raw[left] + 1e-8, raw[right] + 1e-8)
        oldest[node] = max(oldest[left], oldest[right])
        merged = (sizes[a] * distance[a] + sizes[b] * distance[b]) / (sizes[a] + sizes[b])
        sizes[a] += sizes[b]
        cluster_node[a] = node
        active[b] = False
        distance[a, active] = merged[active]
        distance[active, a] = merged[active]
        distance[b, :] = distance[:, b] = np.inf
        distance[a, a] = np.inf

    root = 2 * n - 2
    heights = [tip_ages.get(taxon, 0.0) for taxon in range(n)]
    heights.extend(
        oldest[node] + raw[node] / raw[root] * (root_age - oldest[node])
        for node in range(n, 2 * n - 1)
    )
    result = PhyloTree(parents, heights, taxa, dict(enumerate(labels)))
    if not result.times_valid():
        raise ValueError("initial tree violates a terminal date")
    return result
