"""Catalog constraints change the target explicitly and preserve valid starts."""

import math

import numpy as np
import pytest

from diachron.infer.campaign import build_model, prepare_study, run_campaign
from diachron.infer.initial import distance_starting_tree
from diachron.io.simulate import simulate_dataset
from diachron.models.constraints import satisfies_clades, validate_clades
from diachron.models.posterior import ModelConfig
from diachron.models.tree import PhyloTree


def test_nested_constraints_and_dated_initialization():
    clades = [[0, 1], [0, 1, 2], [3, 4]]
    masks = validate_clades(clades, 6)
    for seed in range(8):
        tree = distance_starting_tree(
            list('abcdef'), [[1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1]],
            8000, np.random.default_rng(seed), {1: 2000, 4: 4000}, clades,
        )
        assert tree.times_valid() and satisfies_clades(tree, masks)
        assert tree.heights[1] == 2000 and tree.heights[4] == 4000
    with pytest.raises(ValueError, match="incompatible"):
        validate_clades([[0, 1], [1, 2]], 4)


def test_composed_sampler_respects_conditioned_topology_target():
    dataset, _ = simulate_dataset(n_languages=6, n_concepts=20, seed=483, name="clade-test")
    study = prepare_study(dataset)
    study.required_clades = [[0, 1], [2, 3]]
    model = build_model(study, ModelConfig.tree_only())
    chains = run_campaign(study, ModelConfig.tree_only(), n_chains=1, warmup=50, iterations=200)
    for draw in chains[0].draws:
        tree = PhyloTree(draw.parent, draw.heights, draw.taxon)
        assert satisfies_clades(tree, validate_clades(study.required_clades, 6))
        assert math.isfinite(model.log_posterior(tree, draw.gain, draw.loss, draw.clock_rate))
    invalid = chains[0].tree.copy()
    # Exchange leaf identities across the two required groups without changing ages.
    invalid.taxon[0], invalid.taxon[2] = invalid.taxon[2], invalid.taxon[0]
    assert model.log_posterior(invalid, .0005, .001, 1) == -math.inf
