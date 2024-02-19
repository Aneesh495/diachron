"""Scaled conditional location messages for a symmetric spatial CTMC."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..models.tree import PhyloTree


class SpatialKernel:
    """Diagonalize a fixed per-year symmetric generator once per spatial model."""

    def __init__(self, generator):
        rates = np.asarray(generator, dtype=float)
        size = len(rates)
        if (not 2 <= size <= 128 or rates.shape != (size, size) or not np.isfinite(rates).all()
                or not np.allclose(rates, rates.T, rtol=0, atol=1e-14)
                or not np.allclose(rates.sum(axis=1), 0, rtol=0, atol=1e-14)
                or np.any(rates[~np.eye(size, dtype=bool)] < 0)):
            raise ValueError("invalid symmetric spatial generator")
        eigenvalues, self._vectors = np.linalg.eigh(rates)
        self._rates = np.minimum(eigenvalues, 0.0)
        self.size = size

    def transitions(self, durations) -> np.ndarray:
        lengths = np.asarray(durations, dtype=float)
        if (lengths.ndim != 1 or not np.isfinite(lengths).all() or np.any(lengths < 0)
                or len(lengths) * self.size * self.size > 8_000_000):
            raise ValueError("invalid spatial transition durations or matrix budget")
        decays = np.exp(lengths[:, None] * self._rates)
        matrices = (self._vectors * decays[:, None, :]) @ self._vectors.T
        if matrices.min(initial=0.0) < -1e-10:
            raise ValueError("spatial transition matrix has negative probabilities")
        np.maximum(matrices, 0.0, out=matrices)
        matrices /= matrices.sum(axis=2, keepdims=True)
        matrices[lengths == 0] = np.eye(self.size)
        return matrices


@dataclass(frozen=True)
class SpatialMarginals:
    log_likelihood: float
    probabilities: np.ndarray


def _scaled(values: np.ndarray) -> tuple[np.ndarray, float]:
    total = float(values.sum())
    if not math.isfinite(total) or total <= 0:
        raise ValueError("spatial observations have zero support")
    return values / total, math.log(total)


def state_marginals(tree: PhyloTree, kernel: SpatialKernel, tip_emissions: dict[int, np.ndarray],
                    rate_scale: float = 1.0) -> SpatialMarginals:
    """Compute every node's state probabilities given all observed tips.

    The emission vector contains likelihoods in state order. Missing tips
    supply an all-ones vector. Each branch uses the same per-year generator.
    """
    roots = [node for node, parent in enumerate(tree.parent) if parent == -1]
    if (len(roots) != 1 or any(parent < -1 or parent >= tree.n_nodes for parent in tree.parent)
            or len(tree.heights) != tree.n_nodes or len(tree.taxon) != tree.n_nodes
            or not math.isfinite(rate_scale) or rate_scale <= 0):
        raise ValueError("invalid tree or spatial rate scale")
    root = roots[0]
    children = tree.children()
    order = [root]
    for node in order:
        order.extend(children[node])
    if len(order) != tree.n_nodes:
        raise ValueError("spatial tree is disconnected")
    tips = [tree.taxon[node] for node in range(tree.n_nodes) if not children[node]]
    taxa = set(tips)
    if (any(children[node] for node in tree.tip_nodes()) or any(tip < 0 for tip in taxa)
            or len(taxa) != len(tips) or not set(tip_emissions) <= taxa):
        raise ValueError("spatial emissions require known terminal taxa")
    transitions = kernel.transitions(rate_scale * np.asarray(tree.branch_lengths()))
    emissions = np.ones((tree.n_nodes, kernel.size))
    for node in tree.tip_nodes():
        taxon = tree.taxon[node]
        if taxon in tip_emissions:
            observation = np.asarray(tip_emissions[taxon], dtype=float)
            if observation.shape != (kernel.size,) or not np.isfinite(observation).all() or np.any(observation < 0):
                raise ValueError("invalid spatial emission")
            emissions[node] = observation

    up = np.empty_like(emissions)
    transported = np.empty_like(emissions)
    log_likelihood = -math.log(kernel.size)
    for node in reversed(order):
        message, scale = _scaled(emissions[node])
        log_likelihood += scale
        for child in children[node]:
            message, scale = _scaled(message * transported[child])
            log_likelihood += scale
        up[node] = message
        if node != root:
            transported[node] = transitions[node] @ message

    down = np.empty_like(emissions)
    down[root] = 1.0 / kernel.size
    probabilities = np.empty_like(emissions)
    for node in order:
        probabilities[node], _ = _scaled(up[node] * down[node])
        prefix = [_scaled(down[node] * emissions[node])[0]]
        for child in children[node]:
            prefix.append(_scaled(prefix[-1] * transported[child])[0])
        suffix = np.ones(kernel.size)
        for index in range(len(children[node]) - 1, -1, -1):
            child = children[node][index]
            complement, _ = _scaled(prefix[index] * suffix)
            down[child], _ = _scaled(transitions[child].T @ complement)
            suffix, _ = _scaled(suffix * transported[child])
    return SpatialMarginals(log_likelihood, probabilities)
