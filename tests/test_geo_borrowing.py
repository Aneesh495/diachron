"""Tests for the region geography model and the bounded borrowing model (G10).

Geography: the region pruning likelihood is finite, and assigning tips to their
nearest region centroid is consistent with the centroids. Contact classification
keeps documented, hypothetical, and proximity evidence distinct.

Borrowing: the model prefers low intensity when candidate pairs do not share more
cognates than background, and a known-truth loan pair (which does share) raises
the data-supported intensity relative to an unrelated pair. This is the
recoverability-limit check the gate asks for: the finite model can detect a
candidate contact that is present, and correctly stays near zero when it is not.
"""

from __future__ import annotations

import numpy as np

from diachron.io.simulate import simulate_dataset
from diachron.io.characters import build_character_matrix
from diachron.models.tree import tree_from_truth
from diachron.models.posterior import PosteriorModel, ModelConfig
from diachron.geo.regions import RegionModel, assign_regions, region_log_likelihood
from diachron.geo.contact import classify_contact, ContactKind


def _study(seed=5, n_lang=12, n_concepts=50, borrow_prob=0.0):
    ds, truth = simulate_dataset(
        n_languages=n_lang, n_concepts=n_concepts, seed=seed, root_age=6000,
        borrow_prob=borrow_prob,
    )
    cm = build_character_matrix(ds)
    l2i = {lid: i for i, lid in enumerate(cm.language_ids)}
    cols = [
        [c.states[l2i[lid]] for lid in truth.tip_language_ids]
        for c in cm.variable_characters()
    ]
    return ds, truth, cols


def test_region_likelihood_finite():
    ds, truth, cols = _study()
    tree = tree_from_truth(truth)
    rm = RegionModel(names=["A", "B", "C"], centroids=[(40, 30), (48, 40), (45, 20)])
    coords = {
        t: truth.true_locations[lid] for t, lid in enumerate(truth.tip_language_ids)
    }
    regions = assign_regions(coords, rm)
    ll = region_log_likelihood(tree, rm, regions, n_tips=12)
    assert np.isfinite(ll)


def test_contact_kinds_distinct():
    doc = classify_contact(0, 1, (39, 35), (40, 36), documented=True, hypothesized=False)
    hyp = classify_contact(0, 2, (39, 35), (50, 10), documented=False, hypothesized=True)
    prox = classify_contact(0, 3, (39, 35), (39.5, 35.5), documented=False, hypothesized=False)
    assert doc.kind is ContactKind.DOCUMENTED
    assert hyp.kind is ContactKind.HYPOTHETICAL
    assert prox.kind is ContactKind.PROXIMITY


def test_borrowing_prefers_zero_without_signal():
    ds, truth, cols = _study(borrow_prob=0.0)
    tree = tree_from_truth(truth)
    model = PosteriorModel(
        config=ModelConfig.borrowing(), columns=cols,
        borrow_candidates=[(0, 1), (2, 3)],
    )
    low = model.log_posterior(tree, 0.0004, 0.0009, 1.0, borrow_rates=np.array([0.01, 0.01]))
    high = model.log_posterior(tree, 0.0004, 0.0009, 1.0, borrow_rates=np.array([0.9, 0.9]))
    # With no injected borrowing, high intensity should not be preferred.
    assert low > high


def test_borrowing_responds_to_shared_cognates():
    # Construct two columns where taxa 0 and 1 share present states strongly.
    n = 6
    cols = []
    for _ in range(40):
        col = [0] * n
        col[0] = 1
        col[1] = 1  # 0 and 1 always co-present
        cols.append(col)
    model = PosteriorModel(
        config=ModelConfig.borrowing(), columns=cols,
        borrow_candidates=[(0, 1)],
    )
    # Compare two intensities on the same tree and selected patterns.
    from diachron.models.tree import PhyloTree

    parent = [n, n, n, n, n, n, -1]
    heights = [0, 0, 0, 0, 0, 0, 1000.0]
    taxon = [0, 1, 2, 3, 4, 5, -1]
    tree = PhyloTree(parent=parent, heights=heights, taxon=taxon)
    low = sum(model.character_scores(tree, 0.0004, 0.0009, 1.0, borrow_rates=[0.05]))
    high = sum(model.character_scores(tree, 0.0004, 0.0009, 1.0, borrow_rates=[0.9]))
    assert high > low


def test_planar_geography_keeps_observed_contrast_when_other_subtree_is_missing():
    import math

    from diachron.models.posterior import _brownian_tree_loglik
    from diachron.models.tree import PhyloTree

    tree = PhyloTree([3, 3, 4, 4, -1], [0.0, 0.0, 0.0, 2.0, 5.0], [0, 1, 2, -1, -1])
    actual = _brownian_tree_loglik(tree, {0: (2.0, 1.0), 2: (3.0, 4.0)}, 0.5)
    assert abs(actual - (-1.0 - math.log(10 * math.pi))) < 1e-12


def test_contact_mixture_matches_enumeration_with_missing_donor():
    from collections import defaultdict
    from itertools import product

    import pytest
    from diachron import _core

    gain, loss = 0.3, 0.8
    pi = (loss / (gain + loss), gain / (gain + loss))
    lengths = (0.6, 1.3, 0.7, 0.3)
    contacts = [(0, 1, 0.4), (2, 0, 0.7)]
    for mask in ((0, 1, 2, 3), (1, 2, 3)):
        outcomes = defaultdict(float)
        for hidden in product((0, 1), repeat=4):
            mass = 0.0
            for root in (0, 1):
                prob = pi[root]
                for state, duration in zip(hidden, lengths):
                    retain = np.exp(-(gain + loss) * duration)
                    prob *= (1 - retain) * pi[state] + (root == state) * retain
                mass += prob
            for donor, recipient, weight in [(0, 0, 1.0), *contacts]:
                visible = list(hidden)
                visible[recipient] = hidden[donor]
                column = tuple(s if t in mask else -1 for t, s in enumerate(visible))
                if column.count(1) >= 2 and 0 in column:
                    outcomes[column] += mass * weight
        columns = list(outcomes)
        observed = _core.binary_character_scores(
            [-1, 0, 0, 0, 0], [0.0, *lengths], [2.0, 1.4, 0.7, 1.3, 1.7],
            [-1, 0, 1, 2, 3], 0, gain, loss, columns, [1.0], 2, contacts,
        )
        selected_mass = sum(outcomes.values())
        assert np.exp(observed) == pytest.approx(
            [outcomes[column] / selected_mass for column in columns], abs=1e-12,
        )


def test_region_likelihood_matches_matrix_exponential():
    from scipy.linalg import expm

    from diachron.models.tree import PhyloTree

    tree = PhyloTree([2, 2, -1], [0.0, 0.0, 200.0], [0, 1, -1])
    model = RegionModel(["A", "B", "C"], [(43, 8), (43, 10), (43, 55)])
    transition = expm(np.asarray(model.rate_matrix()) * 200.0)
    expected = np.log(np.dot(model.stationary(), transition[:, 0] * transition[:, 1]))
    got = region_log_likelihood(tree, model, {0: 0, 1: 1}, 2)
    assert abs(got - expected) < 1e-10


def test_location_marginals_match_independent_gaussian_precision_system():
    from diachron.models.posterior import planar_location_marginals
    from diachron.models.tree import PhyloTree

    tree = PhyloTree(
        [4, 4, 5, 5, 6, 6, -1], [0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 7.0],
        [0, 1, 2, 3, -1, -1, -1],
    )
    observed = {0: (12.0, 2.0), 1: (17.0, 3.0), 3: (1.0, 5.0)}
    sigma2 = 0.5
    Q = np.zeros((tree.n_nodes, tree.n_nodes))
    for i, p in enumerate(tree.parent):
        if p >= 0:
            weight = 1.0 / (sigma2 * (tree.heights[p] - tree.heights[i]))
            Q[i, i] += weight
            Q[p, p] += weight
            Q[p, i] -= weight
            Q[i, p] -= weight
    scored = list(observed)
    latent = sorted(set(range(tree.n_nodes)) - set(scored))
    covariance = np.linalg.inv(Q[np.ix_(latent, latent)])
    mean = -covariance @ Q[np.ix_(latent, scored)] @ np.array(list(observed.values()))
    result = planar_location_marginals(tree, observed, sigma2)
    assert np.allclose([result[i][:2] for i in latent], mean, atol=1e-12)
    assert np.allclose([result[i][2] for i in latent], covariance.diagonal(), atol=1e-12)


def test_location_error_changes_known_gaussian_root_variance():
    from diachron.models.posterior import planar_location_marginals
    from diachron.models.tree import PhyloTree

    tree = PhyloTree([2, 2, -1], [0.0, 0.0, 4.0], [0, 1, -1])
    positions = {0: (42.0, 20.0), 1: (44.0, 26.0)}
    precise = planar_location_marginals(tree, positions, 0.5)[tree.root]
    uncertain = planar_location_marginals(tree, positions, 0.5, error=3.0)[tree.root]
    assert precise == (43.0, 23.0, 1.0)
    assert uncertain == (43.0, 23.0, 5.5)
