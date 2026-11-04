"""Tests for the MCMC diagnostics (gate G06/G14 support).

Each diagnostic is checked against a case with a known answer: iid draws give ESS
near the sample size and Rhat near one; a strongly autocorrelated AR(1) series
gives an ESS close to its analytic value; well-separated chains give a large
Rhat; and the Newick topology signature is invariant to branch lengths and node
ordering but distinguishes genuinely different shapes.
"""

from __future__ import annotations

import numpy as np

from diachron.diagnostics import convergence as cv
from diachron.diagnostics import topology as tp


def test_iid_ess_near_n_and_rhat_near_one():
    rng = np.random.default_rng(0)
    x = rng.standard_normal((4, 2000))
    ess = cv.effective_sample_size(x)
    assert ess > 7000, ess
    assert abs(cv.split_rhat(x) - 1.0) < 0.02
    assert abs(cv.rank_normalized_rhat(x) - 1.0) < 0.02


def test_ar1_ess_matches_theory():
    # AR(1) with phi has asymptotic ESS factor (1 - phi) / (1 + phi).
    rng = np.random.default_rng(1)
    phi = 0.9
    n = 4000
    x = np.zeros((4, n))
    for c in range(4):
        for t in range(1, n):
            x[c, t] = phi * x[c, t - 1] + rng.standard_normal()
    ess = cv.effective_sample_size(x)
    theory = 4 * n * (1 - phi) / (1 + phi)
    # Within a factor of 1.6 of theory (Geyer estimator is conservative).
    assert theory / 1.6 < ess < theory * 1.6, f"ess {ess} vs theory {theory}"


def test_separated_chains_flag_nonconvergence():
    rng = np.random.default_rng(2)
    x = np.array([rng.standard_normal(1000) + 5 * k for k in range(4)])
    assert cv.split_rhat(x) > 1.5


def test_mcse_shrinks_with_ess():
    rng = np.random.default_rng(3)
    x = rng.standard_normal((4, 4000))
    mcse = cv.mcse_mean(x)
    # sd ~ 1, ess ~ 16000, so mcse ~ 1/sqrt(16000) ~ 0.008.
    assert 0.003 < mcse < 0.02, mcse


def test_topology_signature_invariances():
    nw1 = "((t0:1,t1:1):1,(t2:1,t3:1):1);"
    nw1_relabeled_bl = "((t1:9,t0:2):5,(t3:1,t2:7):3);"
    nw2 = "((t0:1,t2:1):1,(t1:1,t3:1):1);"
    assert tp._clade_signature(nw1) == tp._clade_signature(nw1_relabeled_bl)
    assert tp._clade_signature(nw1) != tp._clade_signature(nw2)


def test_topology_mixing_summary():
    nw1 = "((t0:1,t1:1):1,(t2:1,t3:1):1);"
    nw2 = "((t0:1,t2:1):1,(t1:1,t3:1):1);"
    summ = tp.topology_mixing_summary([nw1, nw1, nw2, nw1, nw2])
    assert summ["distinct"] == 2
    assert 0.0 < summ["switch_rate"] <= 1.0
    assert 0.0 < summ["map_prob"] <= 1.0


def test_clade_agreement():
    nw1 = "((t0:1,t1:1):1,(t2:1,t3:1):1);"
    nw2 = "((t0:1,t2:1):1,(t1:1,t3:1):1);"
    # Identical chains agree perfectly.
    assert tp.clade_agreement([nw1, nw1], [nw1, nw1]) == 1.0
    # Disjoint topology distributions agree less.
    assert tp.clade_agreement([nw1, nw1], [nw2, nw2]) < 1.0


def test_constant_separated_chains_have_no_effective_samples():
    x = np.broadcast_to(np.arange(4)[:, None], (4, 200))
    assert np.isinf(cv.split_rhat(x))
    assert cv.bulk_ess(x) == 0.0


def test_rank_folding_detects_different_scale():
    rng = np.random.default_rng(2901)
    x = rng.standard_normal((4, 1000)) * np.array([1, 1, 5, 5])[:, None]
    assert cv.split_rhat(x) < 1.01
    assert cv.rank_normalized_rhat(x) > 1.1
