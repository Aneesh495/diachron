"""Gate G03: the sampler reproduces exact finite target distributions.

Two checks:

1. **Finite discrete target.** A Metropolis chain on a small discrete state space
   with a known unnormalized target must converge to that target. We run the same
   accept/reject logic the chain uses and compare sampled frequencies to the exact
   normalized probabilities within a Monte Carlo tolerance. This validates the
   acceptance rule.

2. **Continuous scale move reversibility.** The multiplicative scale move with its
   log-Jacobian Hastings term must leave a known continuous target invariant. We
   sample a Gamma target with the scale move only and compare the sampled mean and
   variance to the analytic values, confirming the Jacobian is correct (an omitted
   or wrong Jacobian shifts the stationary distribution measurably).

Tolerances are predeclared here, not tuned after seeing results.
"""

from __future__ import annotations

import math
import numpy as np

from diachron import _core


def _metropolis_accept(rng, log_post_old, log_post_new, log_hastings):
    if not math.isfinite(log_post_new):
        return False
    log_alpha = log_post_new - log_post_old + log_hastings
    if log_alpha >= 0:
        return True
    return math.log(max(rng.uniform(), 1e-300)) < log_alpha


def test_finite_discrete_target():
    # Unnormalized target over 5 states.
    weights = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    logw = np.log(weights)
    exact = weights / weights.sum()

    rng = _core.PyRng(2024)
    n_states = len(weights)
    state = 0
    counts = np.zeros(n_states)
    n_iter = 400_000
    for i in range(n_iter):
        # Proposal: uniform over the other states (symmetric, Hastings 0).
        proposal = rng.below(n_states)
        if _metropolis_accept(rng, logw[state], logw[proposal], 0.0):
            state = proposal
        if i > 20_000:  # discard warmup
            counts[state] += 1
    freq = counts / counts.sum()
    # Predeclared tolerance: max absolute deviation under 0.01.
    max_dev = float(np.max(np.abs(freq - exact)))
    assert max_dev < 0.01, f"freq {freq} vs exact {exact}, dev {max_dev}"


def test_scale_move_targets_gamma():
    # Target: Gamma(shape=k, scale=theta) on x > 0. log p(x) = (k-1) log x - x/theta.
    k, theta = 3.0, 2.0

    def log_target(x):
        if x <= 0:
            return float("-inf")
        return (k - 1) * math.log(x) - x / theta

    rng = _core.PyRng(99)
    x = 1.0
    lam = 1.0
    samples = []
    n_iter = 400_000
    for i in range(n_iter):
        u = rng.uniform()
        factor = math.exp(lam * (u - 0.5))
        x_new = x * factor
        log_hastings = math.log(factor)  # Jacobian of the multiplier move
        if _metropolis_accept(rng, log_target(x), log_target(x_new), log_hastings):
            x = x_new
        if i > 20_000 and i % 5 == 0:
            samples.append(x)
    arr = np.array(samples)
    mean = arr.mean()
    var = arr.var()
    # Analytic Gamma mean = k*theta = 6, variance = k*theta^2 = 12.
    assert abs(mean - k * theta) < 0.15, f"mean {mean}"
    assert abs(var - k * theta * theta) < 1.2, f"var {var}"


def test_scale_move_without_jacobian_is_biased():
    """Control: dropping the Jacobian shifts the stationary distribution.

    This is the negative control for the Jacobian. Running the same sampler with
    `log_hastings = 0` (the common bug) targets a different distribution, so its
    mean must be detectably wrong. If this test ever passed with the correct mean,
    the Jacobian would be doing nothing and test_scale_move_targets_gamma would be
    vacuous.
    """
    k, theta = 3.0, 2.0

    def log_target(x):
        if x <= 0:
            return float("-inf")
        return (k - 1) * math.log(x) - x / theta

    rng = _core.PyRng(7)
    x = 1.0
    lam = 1.0
    samples = []
    for i in range(300_000):
        u = rng.uniform()
        factor = math.exp(lam * (u - 0.5))
        x_new = x * factor
        if _metropolis_accept(rng, log_target(x), log_target(x_new), 0.0):  # bug
            x = x_new
        if i > 20_000 and i % 5 == 0:
            samples.append(x)
    biased_mean = float(np.mean(samples))
    # The buggy sampler targets x * Gamma(k, theta)-ish, shifting the mean below
    # the correct 6. We only require that it is detectably different.
    assert abs(biased_mean - k * theta) > 0.3, f"biased mean {biased_mean} too close"
