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

    def log_density(self, age: float) -> float:
        width = self.high - self.low
        fraction = (age - self.low) / width
        density = 0.25 + 2.25 * (1 - fraction) ** 2
        return math.log(self.mass) + math.log(density) - math.log(width)

    def draw_age(self, rng) -> float:
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
                      children: dict[int, list[int]], distances: np.ndarray) -> tuple[list[int], list[AttachmentEdge]]:
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
    anchor = children[backbone[junction]]
    anchor.remove(junction)
    anchor.append(sibling)
    anchor.sort()
    order = [tree.root]
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
    return backbone, [AttachmentEdge(node, low, high, 0.1 / len(edges) + 0.9 * weight / total)
                      for node, low, high, weight in edges]


def resample_attachment(tree: PhyloTree, log_target: float, evaluate: Callable[[PhyloTree], float],
                        distances: np.ndarray, rng) -> tuple[PhyloTree, float]:
    """Keep the current attachment plus four independent proposals, then resample."""
    child = rng.below(tree.n_nodes)
    junction = tree.parent[child]
    if junction < 0 or tree.parent[junction] < 0:
        return tree, log_target
    children = tree.children()
    if len(children[junction]) != 2:
        return tree, log_target
    sibling = next(node for node in children[junction] if node != child)
    backbone, edges = _attachment_edges(tree, child, junction, sibling, children, distances)
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
