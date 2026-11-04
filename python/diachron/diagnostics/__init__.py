"""Convergence and mixing diagnostics for MCMC output.

These implement the standard quantities the inference gate requires: effective
sample size (bulk and tail), split-Rhat, autocorrelation, Monte Carlo standard
error, and topology mixing summaries. They follow the rank-normalized,
split-chain definitions from Vehtari et al. (2021), implemented here so the exact
estimator Diachron reports is auditable; the implementations are cross-checked in
tests against hand-computed small cases and against independent formulas.
"""

from .convergence import (
    effective_sample_size,
    split_rhat,
    autocorrelation,
    mcse_mean,
    rank_normalized_rhat,
    bulk_ess,
    tail_ess,
)
from .topology import topology_frequencies, topology_mixing_summary

__all__ = [
    "effective_sample_size",
    "split_rhat",
    "autocorrelation",
    "mcse_mean",
    "rank_normalized_rhat",
    "bulk_ess",
    "tail_ess",
    "topology_frequencies",
    "topology_mixing_summary",
]
