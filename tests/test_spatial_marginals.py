"""Independent spatial state enumeration and held-out observation checks."""

from dataclasses import replace
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from diachron.geo.marginals import SpatialKernel, state_marginals
from diachron.geo.regions import RegionModel, region_log_likelihood
from diachron.geo.study import _AXES, _case, _generator, _predict, _simulate_sphere, _RATE
from diachron.models.tree import PhyloTree


def test_spatial_messages_match_all_states_and_native_region_likelihood():
    tree = PhyloTree([3, 3, 4, 4, -1], [0, 0, 0, 400, 1200], [0, 1, 2, -1, -1])
    model = RegionModel(["West", "East", "North"], [(45, -15), (42, 30), (67, 12)])
    rates = np.asarray(model.rate_matrix())
    messages = SpatialKernel(rates)
    evidence = {0: np.array([0, 1, 0]), 2: np.array([.1, .2, .8])}
    result = state_marginals(tree, messages, evidence)
    transitions = [expm(rates * years) for years in tree.branch_lengths()]
    weights = np.zeros((tree.n_nodes, 3))
    normalizer = 0.
    for states in product(range(3), repeat=tree.n_nodes):
        mass = evidence[0][states[0]] * evidence[2][states[2]] / 3
        for node, parent in enumerate(tree.parent):
            if parent >= 0:
                mass *= transitions[node][states[parent], states[node]]
        normalizer += mass
        weights[np.arange(tree.n_nodes), states] += mass
    assert result.log_likelihood == pytest.approx(np.log(normalizer), abs=1e-12)
    np.testing.assert_allclose(result.probabilities, weights / normalizer, atol=1e-12)
    exact = state_marginals(tree, messages, {0: np.eye(3)[1], 2: np.eye(3)[2]})
    assert exact.log_likelihood == pytest.approx(region_log_likelihood(tree, model, {0: 1, 2: 2}, 3), abs=1e-10)


def test_predictions_ignore_hidden_tip_states_and_normalize():
    rates = _generator(_AXES, local=False)
    case = _case(65130, "regions", rates)
    heldout_changed = [int((region + 1) % 6) if taxon in case.withheld else region
                       for taxon, region in enumerate(case.tip_regions)]
    kernel = SpatialKernel(rates)
    original = _predict(case, kernel, np.arange(6), 1.)
    changed = _predict(replace(case, tip_regions=heldout_changed), kernel, np.arange(6), 1.)
    assert original == changed
    np.testing.assert_allclose(np.sum(original["heldout"], axis=1), 1.)


def test_spherical_simulator_tracks_first_directional_moment():
    rng = np.random.default_rng(6288)
    point = np.array([1., 0., 0.])
    duration = .10 / _RATE
    results = np.array([_simulate_sphere(point, duration, rng) for _ in range(1500)])
    np.testing.assert_allclose(np.linalg.norm(results, axis=1), 1., atol=1e-12)
    assert results[:, 0].mean() == pytest.approx(np.exp(-.10), abs=.012)
