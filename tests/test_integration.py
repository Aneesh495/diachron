"""End-to-end integration: data -> characters -> campaign -> bundle -> diagnose.

This is the deterministic integration path the demo exercises, run on a small
slice so it finishes quickly in CI. It checks that a short campaign produces a
loadable posterior bundle, that the bundle's draw digest is stable, and that the
verifier accepts an honest bundle and rejects a tampered one (a negative control
for gate G13).
"""

from __future__ import annotations

import json
import numpy as np

from diachron.io.simulate import simulate_dataset, simulate_attestations
from diachron.infer.campaign import prepare_study, run_campaign, campaign_bundle
from diachron.models.posterior import ModelConfig
from diachron.infer.artifacts import PosteriorBundle
from diachron.verify.verifier import verify_evidence


def _small_campaign(tmp_path, seed=1):
    ds, truth = simulate_dataset(
        n_languages=10, n_concepts=50, root_age=6000.0,
        borrow_prob=0.2, missing_frac=0.08, seed=seed, name="diachron_primary",
    )
    atts = simulate_attestations(ds, truth, n_records=30, seed=seed)
    study = prepare_study(ds, atts)
    chains = run_campaign(
        study, ModelConfig.dated(), n_chains=4, warmup=300, iterations=900,
        thin=10, root_age_init=truth.root_age, base_seed=1000,
    )
    bundle = campaign_bundle(chains, study, ModelConfig.dated(), source_commit="testcommit")
    bundle.save(str(tmp_path / "dated"))
    return bundle, study


def test_campaign_produces_loadable_bundle(tmp_path):
    bundle, study = _small_campaign(tmp_path)
    loaded = PosteriorBundle.load(str(tmp_path / "dated"))
    assert loaded.n_chains == 4
    assert loaded.draws_per_chain > 0
    assert loaded.model_name == "dated"
    # Digest is stable across save/load.
    assert loaded.draw_digest() == bundle.draw_digest()
    # Root-age posterior is finite and positive.
    root = np.array(loaded.scalar("root_age"))
    assert np.all(root > 0)


def test_verifier_accepts_honest_bundle(tmp_path):
    _small_campaign(tmp_path)
    result = verify_evidence(str(tmp_path))
    assert result["ok"], result


def test_verifier_rejects_tampered_results(tmp_path):
    _small_campaign(tmp_path)
    # Tamper with a stored draw's root age without updating the manifest digest.
    bundle_path = tmp_path / "dated" / "bundle.json"
    with open(bundle_path) as f:
        payload = json.load(f)
    payload["chains"][0][0]["root_age"] *= 2.0  # altered result
    with open(bundle_path, "w") as f:
        json.dump(payload, f)
    result = verify_evidence(str(tmp_path))
    assert not result["ok"], "verifier failed to catch tampered root age"
    # The failing check must be the draw digest, for a semantic reason.
    digest_checks = [c for c in result["checks"] if c["check"].endswith("draw_digest")]
    assert any(not c["ok"] for c in digest_checks)
