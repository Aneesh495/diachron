"""Validate the path-sampling marginal likelihood on a conjugate problem.

For a conjugate Normal model the marginal likelihood is known in closed form, so
we can check the path-sampling estimator against the truth. This validation is the
precondition for reporting any Bayes factor: an unvalidated estimator could make
an incomparable comparison look rigorous.

Model: data x_1..x_n ~ Normal(mu, sigma^2) with known sigma^2, and a conjugate
prior mu ~ Normal(mu0, tau^2). The marginal likelihood p(x) is a multivariate
Normal with mean mu0 (per element) and covariance sigma^2 I + tau^2 11^T, which
has a closed form we compute directly.
"""

from __future__ import annotations

import math
import numpy as np

from diachron.infer.marginal import power_posterior_path_sampling


def _closed_form_log_marginal(x, mu0, tau2, sigma2):
    """log p(x) under the conjugate Normal model, in closed form."""
    n = len(x)
    # Covariance Sigma = sigma2 I + tau2 * ones ones^T.
    # Use the matrix determinant lemma and Sherman-Morrison for log-density.
    x = np.asarray(x) - mu0
    # Sigma^{-1} = (1/sigma2) (I - (tau2/(sigma2 + n tau2)) ones ones^T)
    s = np.sum(x)
    quad = (np.dot(x, x) / sigma2) - (tau2 * s * s) / (sigma2 * (sigma2 + n * tau2))
    # log det Sigma = (n-1) log sigma2 + log(sigma2 + n tau2)
    logdet = (n - 1) * math.log(sigma2) + math.log(sigma2 + n * tau2)
    return -0.5 * (n * math.log(2 * math.pi) + logdet + quad)


def test_path_sampling_matches_conjugate_marginal():
    rng = np.random.default_rng(0)
    mu_true = 1.5
    sigma2 = 1.0
    n = 8
    x = rng.normal(mu_true, math.sqrt(sigma2), size=n)
    mu0, tau2 = 0.0, 4.0

    closed = _closed_form_log_marginal(x, mu0, tau2, sigma2)

    # Power-posterior sampler: at inverse temperature beta the posterior for mu is
    # Normal with precision 1/tau2 + beta n / sigma2 and mean weighted
    # accordingly. We sample mu exactly and return its data log likelihood.
    xbar = float(np.mean(x))

    def sampler(beta, n_samp, rng_):
        prec = 1.0 / tau2 + beta * n / sigma2
        var = 1.0 / prec
        mean = var * (mu0 / tau2 + beta * n * xbar / sigma2)
        mus = rng_.normal(mean, math.sqrt(var), size=n_samp)
        log_liks = []
        for mu in mus:
            # full-data log likelihood at mu (beta = 1 likelihood)
            ll = -0.5 * (
                n * math.log(2 * math.pi * sigma2)
                + np.sum((x - mu) ** 2) / sigma2
            )
            log_liks.append(ll)
        return mus, np.array(log_liks)

    result = power_posterior_path_sampling(
        log_likelihood_fn=lambda t: 0.0,
        log_prior_sampler=sampler,
        n_betas=24,
        n_samples=2000,
        rng=rng,
    )
    # Path sampling estimates log Z = log marginal; compare to closed form.
    assert abs(result.log_marginal - closed) < 0.15, (
        f"path sampling {result.log_marginal} vs closed form {closed}"
    )


def test_bayes_factor_sign():
    from diachron.infer.marginal import bayes_factor

    assert bayes_factor(-10.0, -12.0) > 0
    assert bayes_factor(-12.0, -10.0) < 0
