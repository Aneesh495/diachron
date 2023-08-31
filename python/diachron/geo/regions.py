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
        """The symmetric distance-decayed generator has uniform equilibrium."""
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
    """Prune under the exact distance-decayed generator, in years."""
    parent, bl, h, tax, root = tree.to_core_arrays()
    n = model.n_regions()
    q = model.rate_matrix()
    pi = model.stationary()
    mean_rate = sum(-p * q[i][i] for i, p in enumerate(pi))
    exchange = [q[i][j] / pi[j] if i != j else 0.0 for i in range(n) for j in range(n)]
    params = [float(n), mean_rate, *pi, *exchange]
    column = [tip_regions.get(t, -1) for t in range(n_tips)]
    return _core.char_log_likelihood(
        parent, bl, h, tax, root, "gtr", params, column, pi
    )
