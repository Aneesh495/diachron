"""Tempering preserves the cold target and exact continuation."""

import math

import numpy as np
from scipy.special import softmax

from diachron import _core
from diachron.infer.campaign import build_model, prepare_study, run_campaign, campaign_bundle
from diachron.infer.chain import Chain, ChainConfig, ChainState
from diachron.infer.replica import ReplicaExchange, exchange_log_ratio
from diachron.io.simulate import simulate_dataset
from diachron.models.posterior import ModelConfig
from diachron.models.tree import random_starting_tree


def test_exchange_preserves_independently_enumerated_joint_target():
    rng = _core.PyRng(4721)
    betas = (1., .3)
    lexical = np.array([-1.5, -.9, -4.2, -2.3])
    prior = np.log([.1, .2, .5, .2])
    scores = [prior + beta * lexical for beta in betas]
    counts = np.zeros((4, 4))
    state = [0, 0]
    for iteration in range(102000):
        for index, score in enumerate(scores):
            candidate = rng.below(4)
            if math.log(max(rng.uniform(), 1e-300)) < score[candidate] - score[state[index]]:
                state[index] = candidate
        ratio = exchange_log_ratio(*betas, lexical[state[0]], lexical[state[1]])
        if math.log(max(rng.uniform(), 1e-300)) < ratio:
            state.reverse()
        if iteration >= 2000:
            counts[tuple(state)] += 1
    exact = np.outer(softmax(scores[0]), softmax(scores[1]))
    assert np.max(np.abs(counts / counts.sum() - exact)) < .01


def test_production_exchange_resume_matches_every_replica(tmp_path):
    dataset, _ = simulate_dataset(n_languages=6, n_concepts=16, seed=701, name="replica-fixture")
    study = prepare_study(dataset)
    model = build_model(study, ModelConfig.tree_only())

    def start():
        replicas = []
        for index, beta in enumerate((1., .85, .65)):
            tree = random_starting_tree(study.tip_labels, 6000, np.random.default_rng(237 + index))
            config = ChainConfig(seed=723 + index, warmup=30, iterations=120,
                                 thin=10 if index == 0 else 121, inverse_temperature=beta)
            replicas.append(Chain(model, tree, ChainState(.0005, .001, 1.), config))
        return ReplicaExchange(replicas, seed=1580, swap_every=5)

    continuous = start()
    continuous.run()
    interrupted = start()
    path = tmp_path / "replicas.json"
    interrupted.run(path, checkpoint_every=11, until_iteration=73)
    resumed = ReplicaExchange.resume(model, path)
    resumed.run()
    assert resumed.statistics() == continuous.statistics()
    assert sum(resumed.accepts) > 0
    assert resumed.rng.snapshot() == continuous.rng.snapshot()
    for actual, expected in zip(resumed.replicas, continuous.replicas):
        assert actual.checkpoint_payload() == expected.checkpoint_payload()
    assert len(resumed.cold.draws) == 12
    assert all(not chain.draws for chain in resumed.replicas[1:])


def test_parallel_exchange_exports_only_cold_draws(tmp_path):
    from diachron.verify.verifier import verify_evidence
    from diachron.util.gitmeta import source_commit

    dataset, _ = simulate_dataset(n_languages=6, n_concepts=16, seed=811, name="parallel-replica")
    study = prepare_study(dataset)
    config = ModelConfig.tree_only()
    kwargs = dict(n_chains=2, warmup=20, iterations=60, temperature_ladder=(1., .8, .5))
    serial = run_campaign(study, config, **kwargs)
    parallel = run_campaign(study, config, **kwargs, jobs=2, checkpoint_dir=str(tmp_path / "states"))
    assert [chain.draws for chain in parallel] == [chain.draws for chain in serial]
    assert [chain.exchange_statistics for chain in parallel] == [chain.exchange_statistics for chain in serial]
    bundle = campaign_bundle(parallel, study, config, source_commit=source_commit())
    bundle.save(str(tmp_path / "bundle"))
    assert verify_evidence(str(tmp_path / "bundle"), expected_commit=source_commit())["ok"]
