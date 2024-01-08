"""Independent finite-grid posterior audit for the four observation models.

The reference sums complete binary node states, integrates a dense Gaussian
precision matrix, and uses SciPy probability densities. It calls no production
pruning, geography, calibration or prior helper. The finite chain uses the
production acceptance method. This checks a conditional grid target, not the
mixing of the unrestricted continuous-time tree sampler.
"""

from __future__ import annotations

from dataclasses import replace
from itertools import combinations, product
import math

import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp
from scipy.stats import expon, lognorm, truncnorm


PATTERNS = ((1, 1, 0, 0), (0, 1, 1, 0), (-1, 1, 0, 1), (1, 0, 1, 0))
COORDINATES = {0: (.1, .7), 1: (.6, .2), 3: (.8, .5)}


def _character_probability(tree, gain, loss, clock, observed, contacts, minimum):
    """Jointly sum inheritance and copy histories before conditioning selection."""
    generator = np.array([[-gain, gain], [loss, -loss]]) * clock
    transitions = {
        node: expm(generator * (tree.heights[parent] - tree.heights[node]))
        for node, parent in enumerate(tree.parent) if parent >= 0
    }
    root_prior = (loss / (gain + loss), gain / (gain + loss))
    numerator = denominator = 0.0
    for states in product((0, 1), repeat=len(tree.parent)):
        inherited = root_prior[states[tree.root]]
        for node, parent in enumerate(tree.parent):
            if parent >= 0:
                inherited *= transitions[node][states[parent], states[node]]
        tips = [states[tree.taxon.index(taxon)] for taxon in range(4)]
        for donor, recipient, weight in [(0, 0, 1.0), *contacts]:
            visible = tips.copy()
            visible[recipient] = tips[donor]
            masked = [value for value, original in zip(visible, observed) if original >= 0]
            retained = masked.count(1) >= minimum and 0 in masked
            if not retained:
                continue
            denominator += inherited * weight
            if all(original < 0 or original == value for original, value in zip(observed, visible)):
                numerator += inherited * weight
    return numerator / denominator


def _geographic_integral(tree, sigma2):
    """Integrate free node locations under independent Gaussian edge densities."""
    count = len(tree.parent)
    precision = np.zeros((count, count))
    log_variances = 0.0
    for node, parent in enumerate(tree.parent):
        if parent < 0:
            continue
        variance = sigma2 * (tree.heights[parent] - tree.heights[node])
        log_variances += math.log(variance)
        precision[node, node] += 1 / variance
        precision[parent, parent] += 1 / variance
        precision[node, parent] -= 1 / variance
        precision[parent, node] -= 1 / variance
    observed = [tree.taxon.index(taxon) for taxon in COORDINATES]
    free = sorted(set(range(count)) - set(observed))
    ff = precision[np.ix_(free, free)]
    fo = precision[np.ix_(free, observed)]
    oo = precision[np.ix_(observed, observed)]
    sign, logdet = np.linalg.slogdet(ff)
    if sign <= 0:
        raise ValueError("singular geographic precision in finite audit")
    marginal = oo - fo.T @ np.linalg.solve(ff, fo)
    score = 0.0
    for dimension in range(2):
        coordinates = np.array([COORDINATES[taxon][dimension] for taxon in COORDINATES])
        score += -.5 * (
            (count - 1 - len(free)) * math.log(2 * math.pi)
            + log_variances + logdet + coordinates @ marginal @ coordinates
        )
    return float(score)


def _independent_posterior(tree, config, gain, loss, clock, sigma2, borrowing, dates):
    contacts = [(0, 2, borrowing)] if config.use_borrowing else []
    score = sum(math.log(_character_probability(
        tree, gain, loss, clock, observed, contacts, 2,
    )) for observed in PATTERNS)
    # Integrate the total lineage-pair exposure on each interval. This is
    # separate from the production event-by-event implementation.
    times = sorted(set(tree.heights))
    exposure = 0.0
    for low, high in zip(times, times[1:]):
        midpoint = (low + high) / 2
        lineages = sum(age < midpoint for node, age in enumerate(tree.heights) if tree.taxon[node] >= 0)
        lineages -= sum(age < midpoint for node, age in enumerate(tree.heights) if tree.taxon[node] < 0)
        exposure += lineages * (lineages - 1) / 2 * (high - low)
    score -= 3 * math.log(config.pop_size) + exposure / config.pop_size
    for value in (gain, loss):
        score += lognorm.logpdf(value, s=config.rate_prior_log_sd, scale=config.rate_prior_median)
    score += lognorm.logpdf(clock, s=config.clock_prior_log_sd, scale=1)
    if config.use_calibration:
        for taxon, (mean, sd) in dates.items():
            age = tree.heights[tree.taxon.index(taxon)]
            score += truncnorm.logpdf(age, -mean / sd, np.inf, loc=mean, scale=sd)
    if config.use_geography:
        score += _geographic_integral(tree, sigma2)
        score += lognorm.logpdf(sigma2, s=config.diffusion_prior_log_sd, scale=config.diffusion_prior_median)
    if config.use_borrowing:
        score += expon.logpdf(borrowing, scale=1 / config.borrowing_prior_rate)
    return float(score)


def _ranked_trees():
    from ..models.tree import PhyloTree

    for first in combinations(range(4), 2):
        remaining = sorted(set(range(4)) - set(first)) + [4]
        for second in combinations(remaining, 2):
            parents = [6] * 6 + [-1]
            for node in first:
                parents[node] = 4
            for node in second:
                parents[node] = 5
            for age in (.4, .65):
                tree = PhyloTree(parents.copy(), [0, 0, .25, age, .9, 1.7, 3.2], [0, 1, 2, 3, -1, -1, -1])
                for gain, loss, clock, sigma2, borrowing in ((.4, .8, .8, .5, .1), (.7, 1.2, 1.1, 1.1, .5)):
                    yield tree, gain, loss, clock, sigma2, borrowing


def run_tiny_posterior_audit(iterations: int = 200000, seed: int = 7159) -> dict:
    from diachron import _core
    from ..infer.chain import Chain
    from ..models.posterior import ModelConfig, PosteriorModel
    from ..models.tip_dates import TipDate

    if not 10000 <= iterations <= 2000000:
        raise ValueError("finite audit supports 10000 to 2000000 iterations")
    results = {}
    dates = {2: (.3, .12), 3: (.5, .15)}
    for model_index, builder in enumerate((ModelConfig.tree_only, ModelConfig.dated,
                                          ModelConfig.geography, ModelConfig.borrowing)):
        config = replace(builder(), pop_size=2.5, rate_prior_median=.5,
                         rate_prior_log_sd=.8, diffusion_prior_median=.8, borrowing_prior_rate=3)
        model = PosteriorModel(
            config, [list(pattern) for pattern in PATTERNS], n_tips=4, minimum_present=2,
            tip_coords=COORDINATES, borrow_candidates=[(0, 2)] if config.use_borrowing else [],
            tip_dates={taxon: TipDate("normal", mean, sd) for taxon, (mean, sd) in dates.items()}
            if config.use_calibration else {},
        )
        reference, production = [], []
        for tree, gain, loss, clock, sigma2, borrowing in _ranked_trees():
            reference.append(_independent_posterior(tree, config, gain, loss, clock, sigma2, borrowing, dates))
            production.append(model.log_posterior(
                tree, gain, loss, clock, sigma2, [borrowing] if config.use_borrowing else [],
            ))
        difference = float(np.max(np.abs(np.array(reference) - production)))
        probabilities = np.exp(np.array(reference) - logsumexp(reference))
        # Call the production acceptance implementation with symmetric proposals
        # over this explicitly finite conditional state space.
        chain = object.__new__(Chain)
        chain.rng = _core.PyRng(seed + model_index)
        state = 0
        chain._log_post = production[state]
        counts = np.zeros(len(reference), dtype=int)
        warmup = 1000
        for iteration in range(iterations + warmup):
            candidate = chain.rng.below(len(reference))
            if chain._accept(production[candidate], 0):
                state = candidate
                chain._log_post = production[state]
            if iteration >= warmup:
                counts[state] += 1
        deviation = float(np.max(np.abs(counts / iterations - probabilities)))
        results[config.name] = {
            "grid_states": len(reference), "max_log_posterior_difference": difference,
            "maximum_frequency_deviation": deviation,
            "passed": difference < 1e-9 and deviation < .01,
        }
    return {
        "scope": "finite conditional joint targets; unrestricted mixing remains separate",
        "seed": seed, "sampling_iterations_per_model": iterations,
        "log_posterior_tolerance": 1e-9, "frequency_tolerance": .01,
        "models": results, "passed": all(result["passed"] for result in results.values()),
    }
