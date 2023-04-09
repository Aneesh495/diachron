"""Prior densities for tree times, rates, and calibrations.

All densities return natural logs. They are written out rather than pulled from
scipy so the exact normalization each term uses is visible and so the posterior
is a single auditable expression; scipy is still used in tests to cross-check
these against an independent implementation.
"""

from __future__ import annotations

import math


def exponential_logpdf(x: float, rate: float) -> float:
    """log density of Exponential(rate) at x >= 0."""
    if x < 0:
        return float("-inf")
    return math.log(rate) - rate * x


def lognormal_logpdf(x: float, mu: float, sigma: float) -> float:
    """log density of LogNormal(mu, sigma) at x > 0."""
    if x <= 0:
        return float("-inf")
    return (
        -math.log(x * sigma * math.sqrt(2 * math.pi))
        - (math.log(x) - mu) ** 2 / (2 * sigma * sigma)
    )


def normal_logpdf(x: float, mu: float, sigma: float) -> float:
    return -0.5 * ((x - mu) / sigma) ** 2 - math.log(sigma * math.sqrt(2 * math.pi))


def uniform_interval_logpdf(x: float, low: float, high: float) -> float:
    """log density of a uniform on [low, high], used for a flat date calibration."""
    if x < low or x > high or high <= low:
        return float("-inf")
    return -math.log(high - low)


def coalescent_log_prior(heights: list[float], parent: list[int], pop_size: float) -> float:
    """Log prior on node times under a constant-size coalescent.

    The coalescent gives the probability of the ordered internal node times given
    the number of lineages between successive coalescences. With `k` lineages the
    pairwise coalescent rate is `k(k-1)/(2 N)`; the waiting time to the next
    coalescence is exponential with that rate. This is a standard tree prior; it
    is implemented here so the density entering the posterior is explicit.

    Returns the log density of the internal node ages. Assumes tips at known
    ages (usually 0) and internal nodes strictly ordered above them.
    """
    n_nodes = len(heights)
    # A node is internal (a coalescence) if it is some other node's parent;
    # otherwise it is a tip (a sampling event).
    is_parent = [False] * n_nodes
    for p in parent:
        if p >= 0:
            is_parent[p] = True
    events = [
        (heights[i], "coal" if is_parent[i] else "sample") for i in range(n_nodes)
    ]
    events.sort(key=lambda e: e[0])

    lp = 0.0
    lineages = 0
    prev_time = 0.0
    for time, kind in events:
        dt = time - prev_time
        if lineages >= 2 and dt > 0:
            rate = lineages * (lineages - 1) / (2.0 * pop_size)
            # Contribution of the waiting interval with `lineages` lineages.
            lp += -rate * dt
        if kind == "sample":
            lineages += 1
        else:  # coalescence reduces lineage count by one and contributes a rate
            if lineages >= 2:
                rate = lineages * (lineages - 1) / (2.0 * pop_size)
                lp += math.log(rate) if rate > 0 else float("-inf")
            lineages -= 1
        prev_time = time
    return lp


def calibration_log_density(node_age: float, low: float, high: float, kind: str = "uniform") -> float:
    """Density of a node-age calibration.

    `uniform` is a flat interval [low, high]; `normal` treats [low, high] as a
    95% interval and uses the matching Gaussian. Calibrations are the external
    temporal evidence (a dated attestation, a defended archaeological bound) and
    are the only thing that makes the strict clock's time scale identifiable.
    """
    if kind == "uniform":
        return uniform_interval_logpdf(node_age, low, high)
    if kind == "normal":
        mu = 0.5 * (low + high)
        # Treat [low, high] as a 95% interval; 1.959964 is the standard normal
        # two-sided 95% quantile.
        sigma = (high - low) / (2 * 1.959963984540054)
        return normal_logpdf(node_age, mu, sigma)
    raise ValueError(f"unknown calibration kind {kind}")
