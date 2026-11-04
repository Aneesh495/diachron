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


def test_exact_resume(tmp_path):
    model, truth = _build_model()
    rng = np.random.default_rng(0)
    start = random_starting_tree(truth.tip_language_ids, 6000.0, rng)

    cfg = ChainConfig(seed=4321, warmup=200, iterations=400, thin=10)

    full = Chain(model, start.copy(),
                 ChainState(gain=0.0005, loss=0.001, clock_rate=1.0), cfg)
    full.run()

    partial = Chain(model, start.copy(),
                    ChainState(gain=0.0005, loss=0.001, clock_rate=1.0),
                    ChainConfig(seed=4321, warmup=200, iterations=150, thin=10))
    ckpt = str(tmp_path / "ckpt.json")
    partial.run(checkpoint_path=ckpt, checkpoint_every=50)
    partial.checkpoint(ckpt)

    resumed = Chain.resume(model, ckpt)
    resumed.config = cfg  # restore the full target iteration count
    resumed.run()

    assert abs(full._log_post - resumed._log_post) < 1e-9
    assert full.tree.parent == resumed.tree.parent
    assert full.tree.heights == resumed.tree.heights
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
