"""Gate G12: interrupted chains resume with exact CPU state.

We run a chain to completion, then run an equivalent chain that is checkpointed
partway, resumed from the checkpoint, and continued. The two must agree
bit-for-bit on the final log posterior, topology, parameters, and the retained
draws, with no draw duplicated across the interruption. We also check that
resuming onto a mismatched model or a corrupted checkpoint is rejected rather
than silently producing wrong state.
"""

from __future__ import annotations

import json
import numpy as np
import pytest

from diachron.io.simulate import simulate_dataset
from diachron.io.characters import build_character_matrix
from diachron.models.tree import random_starting_tree
from diachron.models.posterior import PosteriorModel, ModelConfig
from diachron.infer.chain import Chain, ChainConfig, ChainState


def _build_model():
    ds, truth = simulate_dataset(n_languages=8, n_concepts=40, seed=3, root_age=6000)
    cm = build_character_matrix(ds)
    lang_to_idx = {lid: i for i, lid in enumerate(cm.language_ids)}
    cols = [
        [c.states[lang_to_idx[lid]] for lid in truth.tip_language_ids]
        for c in cm.variable_characters()
    ]
    model = PosteriorModel(config=ModelConfig.tree_only(), columns=cols)
    return model, truth


@pytest.mark.parametrize("attachment_resampling,root_stem,calendar_scale,time_steps", [
    (False, False, False, 0), (True, False, False, 0), (True, True, False, 0),
    (False, False, True, 0), (False, False, False, 1), (False, False, False, 4),
])
def test_exact_resume(tmp_path, attachment_resampling, root_stem, calendar_scale, time_steps):
    model, truth = _build_model()
    rng = np.random.default_rng(0)
    start = random_starting_tree(truth.tip_language_ids, 6000.0, rng)

    cfg = ChainConfig(seed=4321, warmup=200, iterations=400, thin=10,
                      attachment_resampling=attachment_resampling, attachment_root_stem=root_stem,
                      calendar_scale=calendar_scale, time_block=bool(time_steps), time_block_steps=time_steps or 1)

    full = Chain(model, start.copy(),
                 ChainState(gain=0.0005, loss=0.001, clock_rate=1.0), cfg)
    full.run()

    partial = Chain(model, start.copy(),
                    ChainState(gain=0.0005, loss=0.001, clock_rate=1.0),
                    cfg)
    ckpt = str(tmp_path / "ckpt.json")
    partial.run(checkpoint_path=ckpt, checkpoint_every=50, until_iteration=353)

    resumed = Chain.resume(model, ckpt)
    resumed.run()

    assert full._log_post == resumed._log_post
    assert full.draws == resumed.draws
    assert full.rng.snapshot() == resumed.rng.snapshot()
    assert full.tree.parent == resumed.tree.parent
    assert full.tree.heights == resumed.tree.heights
    if attachment_resampling:
        assert resumed.move_accepts["attachment_resampling"] > 0
    if calendar_scale:
        assert resumed.move_accepts["calendar_scale"] > 0
    if time_steps:
        assert resumed.move_accepts["time_block"] > 0
        assert full.time_block_step == resumed.time_block_step
    assert len(full.draws) == len(resumed.draws)
    # No duplicate retained iterations across the interruption.
    iters = [d.iteration for d in resumed.draws]
    assert len(iters) == len(set(iters))


def test_corrupted_checkpoint_rejected(tmp_path):
    model, truth = _build_model()
    rng = np.random.default_rng(1)
    start = random_starting_tree(truth.tip_language_ids, 6000.0, rng)
    cfg = ChainConfig(seed=1, warmup=50, iterations=50, thin=10)
    chain = Chain(model, start.copy(),
                  ChainState(gain=0.0005, loss=0.001, clock_rate=1.0), cfg)
    ckpt = str(tmp_path / "bad.json")
    chain.checkpoint(ckpt)

    with open(ckpt) as f:
        payload = json.load(f)
    del payload["rng_state"]  # corrupt the state
    with open(ckpt, "w") as f:
        json.dump(payload, f)

    with pytest.raises((KeyError, ValueError)):
        Chain.resume(model, ckpt)


def test_checkpoint_rejects_missingness_change(tmp_path):
    from dataclasses import replace

    model, truth = _build_model()
    tree = random_starting_tree(truth.tip_language_ids, 6000.0, np.random.default_rng(0))
    chain = Chain(model, tree, ChainState(0.0005, 0.001, 1.0), ChainConfig(1))
    path = tmp_path / "ckpt.json"
    chain.checkpoint(path)
    altered = [list(c) for c in model.columns]
    altered[0][0] = -1 if altered[0][0] >= 0 else 0
    with pytest.raises(ValueError, match="model inputs differ"):
        Chain.resume(replace(model, columns=altered), path)


def test_parallel_chain_order_and_draws_match_serial(tmp_path):
    from diachron.infer.campaign import prepare_study, run_campaign

    dataset, truth = simulate_dataset(n_languages=6, n_concepts=12, seed=15, root_age=6000)
    study = prepare_study(dataset)
    args = dict(n_chains=2, warmup=10, iterations=40, thin=5, base_seed=32)
    config = ModelConfig.tree_only()
    serial = run_campaign(study, config, **args)
    parallel = run_campaign(study, config, **args, jobs=2, checkpoint_dir=str(tmp_path))
    assert [chain.draws for chain in serial] == [chain.draws for chain in parallel]
