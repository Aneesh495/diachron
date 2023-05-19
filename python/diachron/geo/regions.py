"""Region-based geographic model.

Location is a discrete latent state over a small set of named regions. A lineage
moves between regions as a continuous-time Markov chain, and the likelihood of the
observed tip regions, marginalized over ancestral regions, is computed with the
same Felsenstein pruning engine used for cognate characters. This is appropriate
when the geographic evidence is coarse: a language is "in the steppe" or "in
Anatolia", and the model does not pretend to resolve coordinates the data cannot
support.

Tip regions are assigned from modern reference coordinates by nearest region
centroid, with the explicit understanding that this is a present-day
approximation; the region transition rates then govern how plausibly a lineage
moved between regions over time. Uncertain tip regions are handled by passing a
missing state, which the pruning likelihood marginalizes, just as for cognate
characters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from diachron import _core


@dataclass
class RegionModel:
    """A set of named regions with centroids and pairwise movement rates.

    `names` labels the regions; `centroids[i]` is (lat, lon) of region i. The
    symmetric base rate scales a movement-ease matrix derived from inverse
    great-circle distance between centroids, so nearby regions exchange lineages
    faster than distant ones. The matrix is a CTMC generator (rows sum to zero).
    """

    names: list[str]
    centroids: list[tuple[float, float]]
    base_rate: float = 1e-4  # per year

    def n_regions(self) -> int:
        return len(self.names)

    def rate_matrix(self) -> list[list[float]]:
        n = self.n_regions()
        weights = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                d = _core.haversine_km(
                    self.centroids[i][0], self.centroids[i][1],
                    self.centroids[j][0], self.centroids[j][1],
                )
                # Ease of movement falls off with distance; clamp to avoid
                # division by zero for coincident centroids.
                weights[i][j] = self.base_rate * 1000.0 / max(d, 50.0)
        # Build generator with zero row sums.
        q = [[0.0] * n for _ in range(n)]
        for i in range(n):
            diag = 0.0
            for j in range(n):
                if i != j:
                    q[i][j] = weights[i][j]
                    diag += weights[i][j]
            q[i][i] = -diag
        return q

    def stationary(self) -> list[float]:
        """Approximate stationary distribution by uniform over regions.

        For the symmetric-ish distance-decayed generator the stationary
        distribution is close to uniform; we use uniform as the root prior, which
        is the least-committal choice for an unknown ancestral region.
        """
        n = self.n_regions()
        return [1.0 / n] * n


def assign_regions(
    coords: dict[int, tuple[float, float]], model: RegionModel
) -> dict[int, int]:
    """Assign each tip to its nearest region centroid.

    Returns taxon -> region index. Tips without coordinates are omitted, so the
    caller can pass them as missing to the likelihood.
    """
    out: dict[int, int] = {}
    for taxon, (lat, lon) in coords.items():
        best = 0
        best_d = float("inf")
        for i, (clat, clon) in enumerate(model.centroids):
            d = _core.haversine_km(lat, lon, clat, clon)
            if d < best_d:
                best_d = d
                best = i
        out[taxon] = best
    return out


def region_log_likelihood(
    tree,
    model: RegionModel,
    tip_regions: dict[int, int],
    n_tips: int,
) -> float:
    """Pruning likelihood of the tip regions under the region CTMC.

    Because the Rust binding's k-state model is F81-style (a single exchange
    rate toward an equilibrium), we approximate the distance-decayed generator by
    an F81 chain whose equilibrium is uniform and whose rate matches the mean
    off-diagonal movement rate. This keeps the region likelihood on the same fast
    pruning path as the cognate characters; the approximation is documented and
    compared against the spherical model rather than asserted as exact.
    """
    parent, bl, h, tax, root = tree.to_core_arrays()
    n = model.n_regions()
    q = model.rate_matrix()
    # Mean off-diagonal rate as the F81 exchange rate.
    off = [q[i][j] for i in range(n) for j in range(n) if i != j]
    mean_rate = float(np.mean(off)) if off else model.base_rate
    pi = model.stationary()
    # F81 params: [rate, pi0, pi1, ...]
    params = [max(mean_rate, 1e-9)] + pi
    column = [tip_regions.get(t, -1) for t in range(n_tips)]
    return _core.char_log_likelihood(
        parent, bl, h, tax, root, "f81", params, column, pi
    )
