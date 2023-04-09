"""Generative models and the composed log posterior.

The posterior factorizes as

    log p(tree, times, rates, geo, borrow | data)
        = log p(characters | tree, times, rates, ascertainment)   [likelihood]
        + log p(dates | times)                                    [calibration]
        + log p(geo observations | tree, times, geo params)       [geography]
        + log p(tree)                                             [tree prior]
        + log p(times, rates)                                     [clock prior]
        + log p(geo params) + log p(borrow params)                [hyperpriors]
        + const.

Each term is implemented as its own callable so models can be swapped (tree-only,
dated, geography-conditioned, borrowing-aware) while sharing the pieces they have
in common. Limiting cases are documented where a richer model reduces to a
simpler one (e.g. the relaxed clock collapses to the strict clock as its rate
variance goes to zero).
"""

from .tree import PhyloTree, tree_from_truth
from .clock import StrictClock, RelaxedClock
from .priors import coalescent_log_prior, lognormal_logpdf, exponential_logpdf
from .posterior import PosteriorModel, ModelConfig

__all__ = [
    "PhyloTree",
    "tree_from_truth",
    "StrictClock",
    "RelaxedClock",
    "coalescent_log_prior",
    "lognormal_logpdf",
    "exponential_logpdf",
    "PosteriorModel",
    "ModelConfig",
]
