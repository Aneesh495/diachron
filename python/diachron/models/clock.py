"""Molecular clock models relating calendar time to substitution rate.

A clock turns a branch's duration in years into the number of expected
substitutions on that branch, which is what the character likelihood needs. Two
clocks are offered:

* **Strict clock**: one rate for the whole tree. Expected substitutions on a
  branch equal `rate * duration`. Identifiable when at least one calibration
  pins the time scale; without a calibration, rate and time trade off.

* **Relaxed clock (uncorrelated lognormal)**: each branch draws its own rate
  from a lognormal with shared mean and standard deviation. As the standard
  deviation goes to zero the relaxed clock reduces exactly to the strict clock,
  which is the documented limiting case and a test target.

Units: `rate` is in substitutions per year (the character model's time is then
in substitution units). Durations are in years.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class StrictClock:
    rate: float  # substitutions per year

    def branch_substitutions(self, durations: np.ndarray) -> np.ndarray:
        return self.rate * durations

    def log_prior(self) -> float:
        # Exponential prior on the rate with a weakly informative mean; keeps the
        # rate positive and finite without dominating the likelihood.
        from .priors import exponential_logpdf

        mean_rate = 1e-3
        return exponential_logpdf(self.rate, 1.0 / mean_rate)


@dataclass
class RelaxedClock:
    """Uncorrelated lognormal relaxed clock.

    `mean_rate` is the mean of the branch-rate distribution; `sigma` is the
    lognormal standard deviation on the log scale. `branch_rates` holds the
    per-branch multipliers (one per node, root entry ignored); the effective rate
    on a branch is `mean_rate * branch_rates[i]`.
    """

    mean_rate: float
    sigma: float
    branch_rates: np.ndarray

    def branch_substitutions(self, durations: np.ndarray) -> np.ndarray:
        return self.mean_rate * self.branch_rates * durations

    def log_prior(self) -> float:
        from .priors import exponential_logpdf, lognormal_logpdf

        # Hyperpriors: exponential on the mean rate, half-normal on sigma, and a
        # lognormal(0, sigma) on each branch multiplier (mean 1 in log space).
        lp = exponential_logpdf(self.mean_rate, 1.0 / 1e-3)
        # half-normal density on sigma (sigma >= 0)
        lp += -0.5 * (self.sigma / 0.5) ** 2 - np.log(0.5) - 0.5 * np.log(2 * np.pi)
        for r in self.branch_rates:
            if r <= 0:
                return float("-inf")
            lp += lognormal_logpdf(r, mu=0.0, sigma=max(self.sigma, 1e-6))
        return float(lp)

    def reduces_to_strict(self, tol: float = 1e-6) -> bool:
        """True when sigma is small enough that this behaves like a strict clock."""
        return self.sigma < tol
