"""Scalar convergence diagnostics: ESS, split-Rhat, autocorrelation, MCSE.

The definitions follow the rank-normalized split-Rhat and ESS of Vehtari,
Gelman, Simpson, Carpenter, and Burkner (2021), "Rank-normalization, folding, and
localization: An improved R-hat for assessing convergence of MCMC". The algebra
is implemented directly here rather than taken from a library so the estimator
Diachron reports is explicit and testable.

Inputs are arrays of shape (n_chains, n_draws). Rhat near 1 and large ESS are the
signals that chains have converged and mix well; the acceptance profile targets
split-Rhat below 1.01 and bulk/tail ESS above 400 for declared scalar endpoints.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def _as_2d(samples: np.ndarray) -> np.ndarray:
    arr = np.asarray(samples, dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr


def split_chains(samples: np.ndarray) -> np.ndarray:
    """Split each chain in half, doubling the chain count.

    Splitting guards against within-chain trends that a between-chain comparison
    would otherwise miss: a slowly drifting chain looks fine across chains but its
    two halves disagree.
    """
    arr = _as_2d(samples)
    n_chains, n = arr.shape
    half = n // 2
    if half == 0:
        return arr
    first = arr[:, :half]
    second = arr[:, -half:]
    return np.concatenate([first, second], axis=0)


def _rhat_core(arr: np.ndarray) -> float:
    """Potential scale reduction factor on an (n_chains, n_draws) array."""
    m, n = arr.shape
    if n < 2 or m < 2 or not np.isfinite(arr).all():
        return float("nan")
    chain_means = arr.mean(axis=1)
    chain_vars = arr.var(axis=1, ddof=1)
    between = n * chain_means.var(ddof=1)
    within = chain_vars.mean()
    if within == 0:
        return float("inf") if between > 0 else float("nan")
    var_hat = (n - 1) / n * within + between / n
    return float(np.sqrt(var_hat / within))


def split_rhat(samples: np.ndarray) -> float:
    """Split-Rhat: the PSRF computed on split chains."""
    return _rhat_core(split_chains(samples))


def _rank_normalize(arr: np.ndarray) -> np.ndarray:
    """Rank-normalize values to approximate normal scores (Blom transform)."""
    flat = arr.reshape(-1)
    ranks = stats.rankdata(flat).reshape(arr.shape)
    n = flat.size
    # Blom offset gives approximately normal scores.
    z = stats.norm.ppf((ranks - 0.375) / (n + 0.25))
    return z


def rank_normalized_rhat(samples: np.ndarray) -> float:
    """Rank-normalized split-Rhat, robust to heavy tails and non-normality."""
    split = split_chains(samples)
    if not split.size:
        return float("nan")
    location = _rhat_core(_rank_normalize(split))
    spread = _rhat_core(_rank_normalize(np.abs(split - np.median(split))))
    return max(location, spread) if np.isfinite(spread) else spread


def autocorrelation(x: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Autocorrelation function of a 1D series via FFT."""
    x = np.asarray(x, dtype=float)
    n = x.size
    x = x - x.mean()
    if max_lag is None:
        max_lag = n - 1
    # FFT-based autocovariance.
    size = 1
    while size < 2 * n:
        size *= 2
    f = np.fft.rfft(x, size)
    acov = np.fft.irfft(f * np.conjugate(f), size)[:n].real
    acov /= n
    if acov[0] == 0:
        return np.zeros(min(max_lag + 1, n))
    acf = acov / acov[0]
    return acf[: max_lag + 1]


def effective_sample_size(samples: np.ndarray) -> float:
    """Effective sample size using the Geyer initial-monotone-sequence estimator.

    Combines autocorrelations across chains and truncates the sum where the
    paired autocorrelations turn negative, the standard way to get a stable ESS
    without summing noise from the autocorrelation tail.
    """
    arr = _as_2d(samples)
    m, n = arr.shape
    if n < 4 or not np.isfinite(arr).all():
        return 0.0
    variances = np.where(np.ptp(arr, axis=1) == 0.0, 0.0, arr.var(axis=1, ddof=1))
    within = float(variances.mean())
    if within <= 0:
        return 0.0
    between = arr.mean(axis=1).var(ddof=1) if m > 1 else 0.0
    pooled = within * (n - 1) / n + between
    covariances = np.mean([
        autocorrelation(chain) * variance * (n - 1) / n
        for chain, variance in zip(arr, variances)
    ], axis=0)
    rho = 1.0 - (within - covariances) / pooled
    rho[0] = 1.0

    total = 0.0
    previous = np.inf
    for even in range(0, n - 1, 2):
        pair = rho[even] + rho[even + 1]
        if pair <= 0:
            break
        previous = min(pair, previous)
        total += previous
    tau = max(2.0 * total - 1.0, 1.0)
    return float(m * n / tau)


def bulk_ess(samples: np.ndarray) -> float:
    """Bulk-ESS: ESS of rank-normalized draws (sensitive to the bulk)."""
    return effective_sample_size(_rank_normalize(split_chains(samples)))


def tail_ess(samples: np.ndarray) -> float:
    """Tail-ESS: the smaller of the ESS of the 5% and 95% tail indicators."""
    arr = split_chains(samples)
    if not arr.size:
        return 0.0
    flat = arr.reshape(-1)
    q05, q95 = np.quantile(flat, [0.05, 0.95])
    low = (arr <= q05).astype(float)
    high = (arr >= q95).astype(float)
    return float(min(effective_sample_size(low), effective_sample_size(high)))


def mcse_mean(samples: np.ndarray) -> float:
    """Monte Carlo standard error of the mean: sd / sqrt(ESS)."""
    arr = _as_2d(samples)
    sd = arr.std(ddof=1)
    ess = effective_sample_size(arr)
    if ess <= 0:
        return float("nan")
    return float(sd / np.sqrt(ess))
