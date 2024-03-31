"""Gradient-guided moves of dated tips, internal ages, and the clock."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..models.posterior import PosteriorModel
from ..models.tip_dates import TipDate
from ..models.tree import PhyloTree

if TYPE_CHECKING:
    from .chain import ChainState


class TimeChart:
    """Log positive gaps above the oldest child or a terminal date offset."""

    def __init__(self, tree: PhyloTree, dates: dict[int, TipDate]):
        self.tree = tree
        self.children = tree.children()
        traversal = [tree.root]
        for node in traversal:
            traversal.extend(self.children[node])
        self.nodes = [i for i in reversed(traversal) if self.children[i] or tree.taxon[i] in dates]
        self.dates = {i: dates[tree.taxon[i]] for i in self.nodes if tree.taxon[i] in dates}
        # Tip ages can have much tighter independent measurement uncertainty.
        widths = []
        for node in self.nodes:
            date = self.dates.get(node)
            if date is None:
                widths.append(1.0)
            else:
                relative = date.sd if date.kind == "offset_lognormal" else date.sd / date.mean
                widths.append(min(1.0, max(.005, 3 * relative)))
        self.widths = np.array([*widths, .12])

    def _origin(self, node: int, heights: list[float]) -> tuple[float, int | None]:
        children = self.children[node]
        if children:
            anchor = max(children, key=heights.__getitem__)
            return heights[anchor], anchor
        return self.dates[node].offset, None

    def encode(self, clock: float) -> tuple[np.ndarray, list[int | None]]:
        gaps, anchors = [], []
        for node in self.nodes:
            origin, anchor = self._origin(node, self.tree.heights)
            gaps.append(math.log(self.tree.heights[node] - origin))
            anchors.append(anchor)
        return np.array([*gaps, math.log(clock)]), anchors

    def decode(self, coordinates: np.ndarray) -> tuple[PhyloTree, float, list[int | None]]:
        heights = list(self.tree.heights)
        anchors = []
        for node, coordinate in zip(self.nodes, coordinates):
            origin, anchor = self._origin(node, heights)
            value = origin + math.exp(coordinate)
            if not origin < value < math.inf:
                raise ValueError("nonpositive or nonfinite proposed time gap")
            heights[node] = value
            anchors.append(anchor)
        clock = math.exp(coordinates[-1])
        if not 0 < clock < math.inf:
            raise ValueError("nonpositive or nonfinite proposed clock")
        return PhyloTree(self.tree.parent, heights, self.tree.taxon, self.tree.tip_labels), clock, anchors

    def pullback(self, age_gradient: np.ndarray, coordinates: np.ndarray, anchors: list[int | None]) -> np.ndarray:
        # Older nodes can move with their anchor. Pull their sensitivities back
        # before reading the child's gap derivative.
        result = np.ones(len(coordinates))
        for j in reversed(range(len(self.nodes))):
            node, anchor = self.nodes[j], anchors[j]
            force = age_gradient[node]
            result[j] += force * math.exp(coordinates[j])
            if anchor is not None:
                age_gradient[anchor] += force
        return result


def age_prior_gradient(tree: PhyloTree, pop_size: float, dates: dict[int, TipDate]) -> np.ndarray:
    """Derivative of the serial coalescent and measured terminal date densities."""
    gradient = np.zeros(tree.n_nodes)
    lineages = 0
    for node in sorted(range(tree.n_nodes), key=lambda i: (tree.heights[i], tree.taxon[i] < 0)):
        taxon = tree.taxon[node]
        # Moving the event adds time with the previous lineage count and
        # removes the same time with the next count.
        if taxon < 0:
            lineages -= 1
            gradient[node] = -lineages / pop_size
        else:
            gradient[node] = lineages / pop_size
            lineages += 1
            if taxon in dates:
                date = dates[taxon]
                age = tree.heights[node] - date.offset
                if date.kind == "normal":
                    gradient[node] -= (age - date.mean) / date.sd ** 2
                else:
                    log_mean = math.log(date.mean) - date.sd ** 2 / 2
                    gradient[node] -= (1 + (math.log(age) - log_mean) / date.sd ** 2) / age
    return gradient


def propose_time_block(
    model: PosteriorModel, tree: PhyloTree, state: ChainState, inverse_temperature: float,
    step: float, rng,
) -> tuple[PhyloTree, float, float] | None:
    """Return the Gaussian proposal with its reverse density and chart Jacobian."""
    dates = model.tip_dates if model.config.use_calibration else {}
    chart = TimeChart(tree, dates)
    widths = step * chart.widths

    def drift(candidate: PhyloTree, clock: float, coordinates: np.ndarray, anchors) -> np.ndarray:
        branches = np.asarray(model.character_branch_gradient(
            candidate, state.gain, state.loss, clock, state.borrow_rates, state.clock_multipliers,
        )) * inverse_temperature
        ages = age_prior_gradient(candidate, model.config.pop_size, dates)
        for child, parent in enumerate(candidate.parent):
            if parent >= 0:
                ages[parent] += branches[child]
                ages[child] -= branches[child]
        gradient = chart.pullback(ages, coordinates, anchors)
        gradient[-1] = (branches @ candidate.branch_lengths()
                        - coordinates[-1] / model.config.clock_prior_log_sd ** 2)
        if not np.isfinite(gradient).all():
            raise ValueError("nonfinite time proposal force")
        return np.clip(.5 * widths ** 2 * gradient, -widths, widths)

    try:
        current, anchors = chart.encode(state.clock_rate)
        forward = drift(tree, state.clock_rate, current, anchors)
        noise = np.array([rng.normal() for _ in widths])
        proposed = current + forward + widths * noise
        candidate, clock, new_anchors = chart.decode(proposed)
        reverse = (current - proposed - drift(candidate, clock, proposed, new_anchors)) / widths
    except (ValueError, OverflowError):
        return None
    log_ratio = .5 * (noise @ noise - reverse @ reverse)
    log_jacobian = float((proposed - current).sum())
    return candidate, clock, float(log_ratio + log_jacobian)
