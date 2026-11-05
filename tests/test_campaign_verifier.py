"""Campaign diagnostics reject plausible scalar and topology failures."""

from copy import deepcopy
from dataclasses import asdict
import json

import numpy as np

from diachron.infer.campaign import prepare_study, run_campaign, campaign_bundle
from diachron.io.simulate import simulate_dataset
from diachron.verify.campaign import _checkpoint_matches, _diagnose


def test_full_precision_campaign_diagnostics_detect_isolated_clade_modes():
    rng = np.random.default_rng(2971)
    topologies = [
        "((a:1,b:1):1,(c:1,d:1):1);",
        "((a:1,c:1):1,(b:1,d:1):1);",
    ]
    chains = []
    for _ in range(4):
        order = rng.permutation(np.repeat([0, 1], 1000))
        chains.append([
            {"root_age": rng.normal(6000, 100), "gain": rng.lognormal(-8, .1),
             "loss": rng.lognormal(-7, .1), "clock_rate": rng.lognormal(0, .1),
             "newick": topologies[index]}
            for index in order
        ])
    raw = {"chains": chains}
    assert _diagnose(raw)["converged"]
    bad = deepcopy(raw)
    for index, chain in enumerate(bad["chains"]):
        for draw in chain:
            draw["newick"] = topologies[index % 2]
    assert not _diagnose(bad)["converged"]
    bad = deepcopy(raw)
    for draw in bad["chains"][0]:
        draw["root_age"] += 1000
    assert not _diagnose(bad)["converged"]


def test_final_checkpoint_binds_model_draws_and_move_accounting(tmp_path):
    from diachron.models.posterior import ModelConfig
    from diachron.util.gitmeta import source_commit

    dataset, _ = simulate_dataset(n_languages=6, n_concepts=12, seed=61, name="checkpoint-audit")
    study = prepare_study(dataset)
    config = ModelConfig.tree_only()
    chains = run_campaign(
        study, config, n_chains=1, warmup=10, iterations=30,
        checkpoint_dir=str(tmp_path),
    )
    bundle = campaign_bundle(chains, study, config, source_commit=source_commit())
    raw = json.loads(json.dumps(asdict(bundle)))
    checkpoint = json.loads((tmp_path / "tree_only_chain0.json").read_text())
    plan = {"source_commit": source_commit()}
    assert _checkpoint_matches(checkpoint, raw, 0, plan)
    for mutation in ("model", "draw", "move"):
        altered = deepcopy(checkpoint)
        if mutation == "model":
            altered["model_digest"] = "0" * 64
        elif mutation == "draw":
            altered["draws"][0]["gain"] *= 1.01
        else:
            altered["move_accepts"]["gain"] = altered["move_proposals"].get("gain", 0) + 1
        assert not _checkpoint_matches(altered, raw, 0, plan)
