"""Conditional importance resampling of a subtree's attachment edge and age."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math

import numpy as np

from ..models.tree import PhyloTree


@dataclass(frozen=True)
class AttachmentEdge:
    child: int
    low: float
    high: float
    mass: float
    tail_scale: float | None = None

    def log_density(self, age: float) -> float:
        if self.tail_scale is not None:
            return math.log(self.mass / self.tail_scale) - (age - self.low) / self.tail_scale
        width = self.high - self.low
        fraction = (age - self.low) / width
        density = 0.25 + 2.25 * (1 - fraction) ** 2
        return math.log(self.mass) + math.log(density) - math.log(width)

    def draw_age(self, rng) -> float:
        if self.tail_scale is not None:
            return self.low - self.tail_scale * math.log1p(-rng.uniform())
        fraction = rng.uniform()
        if rng.uniform() >= 0.25:
            fraction = min(fraction, rng.uniform(), rng.uniform())
        return self.low + (self.high - self.low) * fraction


def _choose(weights: list[float], rng) -> int:
    threshold = rng.uniform() * math.fsum(weights)
    for index, weight in enumerate(weights):
        threshold -= weight
        if threshold < 0:
            return index
    return max(index for index, weight in enumerate(weights) if weight > 0)


def _attachment_edges(tree: PhyloTree, child: int, junction: int, sibling: int,
                      children: dict[int, list[int]], distances: np.ndarray,
                      stem_scale: float | None) -> tuple[list[int], list[AttachmentEdge]]:
    backbone = list(tree.parent)
    backbone[sibling] = backbone[junction]

    detached = []
    removed = {junction}
    stack = [child]
    while stack:
        node = stack.pop()
        removed.add(node)
        if tree.is_tip(node):
            detached.append(tree.taxon[node])
        stack.extend(children[node])
    nearest = distances[sorted(detached)].mean(axis=0)
    # The same labeled backbone is seen from every position of this junction.
    grandparent = backbone[junction]
    if grandparent >= 0:
        anchor = children[grandparent]
        anchor.remove(junction)
        anchor.append(sibling)
        anchor.sort()
    root = sibling if grandparent < 0 else tree.root
    order = [root]
    for node in order:
        order.extend(children[node])
    totals = {}
    counts = {}
    for node in reversed(order):
        if tree.is_tip(node):
            totals[node], counts[node] = float(nearest[tree.taxon[node]]), 1
        else:
            totals[node] = math.fsum(totals[sub] for sub in children[node])
            counts[node] = sum(counts[sub] for sub in children[node])
    edges = []
    for node, parent in enumerate(backbone):
        if parent >= 0 and node not in removed:
            low = max(tree.heights[node], tree.heights[child])
            high = tree.heights[parent]
            if low < high:
                distance = totals[node] / counts[node]
                edges.append((node, low, high, math.exp(-12 * distance)))
    total = math.fsum(edge[3] for edge in edges)
    stem_mass = (0.18 if edges else 1.0) if stem_scale is not None else 0.0
    proposals = [AttachmentEdge(node, low, high, (1 - stem_mass) * (0.1 / len(edges) + 0.9 * weight / total))
                 for node, low, high, weight in edges]
    if stem_scale is not None:
        proposals.append(AttachmentEdge(root, max(tree.heights[root], tree.heights[child]),
                                        math.inf, stem_mass, stem_scale))
    return backbone, proposals


def resample_attachment(tree: PhyloTree, log_target: float, evaluate: Callable[[PhyloTree], float],
                        distances: np.ndarray, rng, *, stem_scale: float | None = None) -> tuple[PhyloTree, float]:
    """Keep the current attachment plus four independent proposals, then resample."""
    child = rng.below(tree.n_nodes)
    junction = tree.parent[child]
    if junction < 0 or (tree.parent[junction] < 0 and stem_scale is None):
        return tree, log_target
    children = tree.children()
    if len(children[junction]) != 2:
        return tree, log_target
    sibling = next(node for node in children[junction] if node != child)
    backbone, edges = _attachment_edges(tree, child, junction, sibling, children, distances, stem_scale)
    current = next((edge for edge in edges if edge.child == sibling), None)
    if current is None or not current.low < tree.heights[junction] < current.high:
        return tree, log_target
    candidates = [(tree, log_target)]
    log_weights = [log_target - current.log_density(tree.heights[junction])]
    masses = [edge.mass for edge in edges]
    for _ in range(4):
        edge = edges[_choose(masses, rng)]
        age = edge.draw_age(rng)
        if not edge.low < age < edge.high:
            continue
        proposed = tree.copy()
        proposed.parent = list(backbone)
        proposed.parent[junction] = backbone[edge.child]
        proposed.parent[edge.child] = junction
        proposed.heights[junction] = age
        score = evaluate(proposed)
        if math.isfinite(score):
            candidates.append((proposed, score))
            log_weights.append(score - edge.log_density(age))
    maximum = max(log_weights)
    return candidates[_choose([math.exp(weight - maximum) for weight in log_weights], rng)]
