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
import pytest

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


def test_production_topology_proposals_recover_exact_ranked_tree_posterior():
    """Condition on three event ages and enumerate every ranked four-tip tree."""
    from collections import Counter
    from itertools import combinations, product

    from diachron import _core
    from diachron.infer.moves import MoveSet
    from diachron.models.tree import PhyloTree

    taxa = [0, 1, 2, 3, -1, -1, -1]
    heights = [0.0, 0.0, 0.0, 0.0, 0.7, 1.6, 3.2]
    patterns = ([1, 1, 0, 0], [0, 1, 0, 1], [0, 0, 1, 1])
    gain, loss = 0.4, 0.8
    prior = [loss / (gain + loss), gain / (gain + loss)]

    def exact_log_target(parent):
        score = 0.0
        for pattern in patterns:
            mass = 0.0
            for ancestral in product((0, 1), repeat=3):
                states = list(pattern) + list(ancestral)
                term = prior[states[6]]
                for node in range(6):
                    p = parent[node]
                    retain = np.exp(-(gain + loss) * (heights[p] - heights[node]))
                    term *= (1 - retain) * prior[states[node]] + retain * (states[node] == states[p])
                mass += term
            score += np.log(mass)
        return float(score)

    exact = {}
    for first in combinations(range(4), 2):
        available = sorted(set(range(4)) - set(first)) + [4]
        for second in combinations(available, 2):
            parents = [6] * 6 + [-1]
            for node in first:
                parents[node] = 4
            for node in second:
                parents[node] = 5
            exact[tuple(parents)] = exact_log_target(parents)
    assert len(exact) == 18
    scale = sum(np.exp(score) for score in exact.values())
    probabilities = {key: np.exp(value) / scale for key, value in exact.items()}

    move_set = MoveSet()
    rng = _core.PyRng(30922)
    state = next(iter(exact))
    counts = Counter()
    for step in range(61_000):
        candidate = PhyloTree(list(state), list(heights), list(taxa))
        proposal = move_set.nni if rng.uniform() < 0.5 else move_set.spr
        move = proposal(candidate, rng)
        if move.valid:
            dest = tuple(candidate.parent)
            assert dest in exact
            log_ratio = exact[dest] - exact[state] + move.log_hastings
            if np.log(max(rng.uniform(), 1e-300)) < log_ratio:
                state = dest
        if step >= 1000:
            counts[state] += 1
    assert len(counts) == len(exact)
    assert max(abs(counts[key] / 60_000 - p) for key, p in probabilities.items()) < 0.015


def test_rate_clock_move_targets_analytic_conditional_lognormal():
    from diachron import _core
    from diachron.infer.moves import MoveSet
    from diachron.models.priors import lognormal_logpdf

    rng = _core.PyRng(3279)
    moves = MoveSet()
    moves.scale_lambda = 2.5
    means = np.array([-1.0, -0.2, 0.3])
    widths = np.array([0.7, 0.9, 0.5])
    effective_gain, effective_loss = 0.6, 1.2
    precision = 1 / widths ** 2
    expected_mean = np.dot(
        [np.log(effective_gain) - means[0], np.log(effective_loss) - means[1], means[2]],
        precision,
    ) / precision.sum()
    expected_variance = 1 / precision.sum()

    def target(state):
        return sum(lognormal_logpdf(value, mean, width)
                   for value, mean, width in zip(state, means, widths))

    state = (effective_gain, effective_loss, 1.0)
    logp = target(state)
    logs = []
    for iteration in range(31_000):
        proposed, result = moves.rate_clock_ridge(*state, rng)
        score = target(proposed)
        if np.log(max(rng.uniform(), 1e-300)) < score - logp + result.log_hastings:
            state, logp = proposed, score
        if iteration >= 1000:
            logs.append(np.log(state[2]))
    assert abs(np.mean(logs) - expected_mean) < 0.02
    assert abs(np.var(logs) - expected_variance) < 0.02
    assert state[0] * state[2] == pytest.approx(effective_gain, rel=1e-12)
    assert state[1] * state[2] == pytest.approx(effective_loss, rel=1e-12)


def test_height_resampled_spr_targets_uniform_four_tip_time_volume():
    """Balanced trees have twice the valid age volume of pectinate trees."""
    from collections import Counter

    from diachron.infer.moves import MoveSet
    from diachron.models.constraints import clade_masks
    from diachron.models.tree import PhyloTree

    tree = PhyloTree([4, 4, 5, 5, 6, 6, -1], [0, 0, 0, 0, .3, .7, 1], [0, 1, 2, 3, -1, -1, -1])
    rng = _core.PyRng(8541)
    moves = MoveSet()
    draws = 60000
    counts = Counter()
    for index in range(draws + 2000):
        proposed = tree.copy()
        move = moves.spr(proposed, rng, resample_height=True)
        if move.valid and math.log(max(rng.uniform(), 1e-300)) < move.log_hastings:
            assert all(proposed.heights[p] > proposed.heights[n] for n, p in enumerate(proposed.parent) if p >= 0)
            tree = proposed
        if index >= 2000:
            counts[frozenset(clade_masks(tree) - {15})] += 1
    # For each of three balanced topologies the two internal heights fill a
    # unit square. Each of twelve pectinate topologies occupies half a square.
    # The internal node indexing multiplicity is the same for both shapes.
    assert len(counts) == 15
    for clades, count in counts.items():
        volume = 1.0 if all(mask.bit_count() == 2 for mask in clades) else .5
        assert abs(count / draws - volume / 9) < .012


def test_attachment_resampling_corrects_nonuniform_proposal_to_continuous_target():
    from collections import Counter

    from diachron.infer.attachment import resample_attachment
    from diachron.models.constraints import clade_masks
    from diachron.models.tree import PhyloTree

    tree = PhyloTree([4, 4, 5, 5, 6, 6, -1], [0, 0, 0, 0, .3, .7, 1], [0, 1, 2, 3, -1, -1, -1])
    distances = np.array([[0, .1, .6, .7], [.1, 0, .4, .5], [.6, .4, 0, .2], [.7, .5, .2, 0]])
    rng = _core.PyRng(8552)
    draws = 40000
    counts = Counter()
    for iteration in range(draws + 2000):
        tree, log_target = resample_attachment(tree, 0.0, lambda _: 0.0, distances, rng)
        assert log_target == 0.0
        assert all(tree.heights[p] > tree.heights[n] for n, p in enumerate(tree.parent) if p >= 0)
        if iteration >= 2000:
            counts[frozenset(clade_masks(tree) - {15})] += 1
    assert len(counts) == 15
    for clades, count in counts.items():
        volume = 1.0 if all(mask.bit_count() == 2 for mask in clades) else .5
        assert abs(count / draws - volume / 9) < .012


@pytest.mark.parametrize("parsimony", [False, True])
def test_root_stem_resampling_targets_unbounded_four_tip_age_distribution(parsimony):
    from collections import Counter

    from diachron.infer.attachment import resample_attachment
    from diachron.models.constraints import clade_masks
    from diachron.models.tree import PhyloTree

    tree = PhyloTree([4, 4, 5, 5, 6, 6, -1], [0, 0, 0, 0, .3, .7, 1], [0, 1, 2, 3, -1, -1, -1])
    distances = np.array([[0, .1, .6, .7], [.1, 0, .4, .5], [.6, .4, 0, .2], [.7, .5, .2, 0]])
    columns = [[1, 1, 0, 0], [0, -1, 0, 1], [0, 1, 1, 0]] * 5 if parsimony else None
    rng = _core.PyRng(8599)
    counts, roots, ages = Counter(), set(), []

    def target(candidate):
        return -candidate.root_age()

    score = target(tree)
    for iteration in range(52000):
        tree, score = resample_attachment(tree, score, target, distances, rng, stem_scale=1.7, columns=columns)
        assert score == target(tree)
        assert all(tree.heights[p] > tree.heights[n] for n, p in enumerate(tree.parent) if p >= 0)
        if iteration >= 2000:
            counts[frozenset(clade_masks(tree) - {15})] += 1
            roots.add(tree.root)
            ages.append(tree.root_age())
    # Integrate over the other two internal ages at a fixed root age r.
    # Their volume scales as r², so the marginal root law is Gamma(3, 1).
    assert roots == {4, 5, 6} and len(counts) == 15
    assert abs(np.mean(ages) - 3.0) < .14
    assert abs(np.mean(np.array(ages) <= 3) - (1 - 8.5 * np.exp(-3))) < .025
    for clades, count in counts.items():
        volume = 1.0 if all(mask.bit_count() == 2 for mask in clades) else .5
        assert abs(count / len(ages) - volume / 9) < .015
