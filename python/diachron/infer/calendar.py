"""Move along a common calendar scale with reciprocal rate changes."""

from __future__ import annotations

from collections.abc import Callable
import math

from ..models.posterior import PosteriorModel
from ..models.priors import coalescent_log_prior, lognormal_logpdf
from ..models.tree import PhyloTree


def slice_offset(score: Callable[[float], float], start_score: float, rng) -> float:
    """Shrink a randomly placed unit interval around the current offset zero."""
    level = start_score + math.log(max(rng.uniform(), 1e-300))
    left = -rng.uniform()
    right = left + 1.0
    for _ in range(64):
        candidate = left + (right - left) * rng.uniform()
        if score(candidate) >= level:
            return candidate
        if candidate < 0:
            left = candidate
        else:
            right = candidate
    return 0.0


def propose_calendar_scale(
    model: PosteriorModel, tree: PhyloTree, clock: float, diffusion: float, rng,
) -> tuple[PhyloTree, float, float, float] | None:
    """Return a prior slice proposal and the reverse/forward correction.

    Absolute time scaling leaves lexical exposure unchanged when the clock
    changes reciprocally. The geography diffusion coefficient has the same
    units relationship. The final MH step uses the full posterior regardless.
    """
    dated = model.tip_dates.keys() if model.config.use_calibration else ()
    if any(age != 0 and taxon >= 0 and taxon not in dated
           for age, taxon in zip(tree.heights, tree.taxon)):
        return None
    degrees = sum(taxon < 0 or taxon in dated for taxon in tree.taxon)
    geographic = model.config.use_geography
    exponent = degrees - 1 - geographic
    proposal = tree.copy()
    cfg = model.config

    def rates(factor: float) -> tuple[float, float]:
        return clock / factor, diffusion / factor if geographic else diffusion

    def scale_prior(offset: float) -> float:
        factor = math.exp(offset)
        proposal.heights = [height * factor for height in tree.heights]
        proposed_clock, proposed_diffusion = rates(factor)
        if (not all(math.isfinite(height) for height in proposal.heights)
                or not 0 < proposed_clock < math.inf or not 0 < proposed_diffusion < math.inf
                or not proposal.times_valid()):
            return -math.inf
        score = (model.calibration_log_density(proposal)
                 + coalescent_log_prior(proposal.heights, proposal.parent, cfg.pop_size)
                 + lognormal_logpdf(proposed_clock, 0.0, cfg.clock_prior_log_sd))
        if geographic:
            score += lognormal_logpdf(proposed_diffusion, math.log(cfg.diffusion_prior_median),
                                      cfg.diffusion_prior_log_sd)
        return score + exponent * offset

    initial = scale_prior(0.0)
    if not math.isfinite(initial):
        return None
    offset = slice_offset(scale_prior, initial, rng)
    if offset == 0:
        return None
    correction = initial - scale_prior(offset) + exponent * offset
    return proposal, *rates(math.exp(offset)), correction
