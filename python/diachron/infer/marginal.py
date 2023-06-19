"""Marginal likelihood estimation for model comparison (D22).

Comparing models by Bayes factors requires the marginal likelihood (the evidence),
the normalizing constant of the posterior. Raw likelihoods from different
observation models are not comparable and must never be treated as Bayes factors;
the estimator here computes a proper (if approximate) marginal likelihood so a
comparison, when reported, is between comparable quantities.

The method is path sampling (thermodynamic integration) along a temperature
ladder. For an inverse temperature beta in [0, 1], the power posterior is
proportional to prior times likelihood^beta. The log marginal likelihood is the
integral over beta of the expected log likelihood under the power posterior:

    log Z = integral_0^1 E_beta[ log p(data | theta) ] d beta.

We estimate the inner expectation by short MCMC at each beta and integrate with
the trapezoid rule. The estimator is validated on a conjugate problem where the
marginal likelihood is known in closed form, so a reported Bayes factor rests on
a checked estimator rather than an unvalidated one. This is why the brief makes
Bayes factors optional only if the estimator is validated: the validation is the
precondition, and it lives in `tests/test_marginal.py`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PathSamplingResult:
    log_marginal: float
    betas: list[float]
    mean_log_lik: list[float]


def power_posterior_path_sampling(
    log_likelihood_fn,
    log_prior_sampler,
    *,
    n_betas: int = 16,
    n_samples: int = 400,
    rng: np.random.Generator,
) -> PathSamplingResult:
    """Estimate log marginal likelihood by path sampling over a beta ladder.

    `log_likelihood_fn(theta) -> float` is the log likelihood at a parameter.
    `log_prior_sampler(beta, n, rng) -> (thetas, log_liks)` returns `n` samples
    from the power posterior at inverse temperature `beta` and their log
    likelihoods. For the validated conjugate test the power posterior is sampled
    exactly; for the phylogenetic models it is approximated by short MCMC.

    A power ladder (beta_k = (k/K)^5) concentrates points near beta = 0 where the
    integrand changes fastest, which reduces the discretization bias of the
    trapezoid rule.
    """
    betas = [(k / n_betas) ** 5 for k in range(n_betas + 1)]
    mean_ll = []
    for beta in betas:
        _thetas, log_liks = log_prior_sampler(beta, n_samples, rng)
        mean_ll.append(float(np.mean(log_liks)))
        _ = log_likelihood_fn  # kept for callers that pass an explicit evaluator

    # Trapezoid integration of E_beta[log lik] over beta.
    log_z = 0.0
    for i in range(1, len(betas)):
        db = betas[i] - betas[i - 1]
        log_z += 0.5 * (mean_ll[i] + mean_ll[i - 1]) * db
    return PathSamplingResult(log_marginal=log_z, betas=betas, mean_log_lik=mean_ll)


def bayes_factor(log_marginal_a: float, log_marginal_b: float) -> float:
    """log Bayes factor of model A over model B (difference of log marginals).

    Positive favors A. This is only meaningful when both marginals come from the
    validated estimator on comparable data; the caller is responsible for not
    passing raw, unnormalized likelihoods here.
    """
    return log_marginal_a - log_marginal_b
