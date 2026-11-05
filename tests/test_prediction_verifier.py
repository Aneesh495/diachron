"""Independent predictive score replay catches changes to stored observations."""

import json

from diachron.io.simulate import simulate_dataset, simulate_attestations
from diachron.ppc.heldout import run_heldout_prediction
from diachron.verify.prediction import verify_prediction


def test_prediction_replay_and_tampered_score(tmp_path):
    dataset, truth = simulate_dataset(
        n_languages=8, n_concepts=30, root_age=6000.0, seed=327, name="prediction-test",
    )
    attestations = simulate_attestations(dataset, truth, n_records=12, seed=27)
    out = tmp_path / "prediction"
    run_heldout_prediction(
        str(out), dataset=dataset, attestations=attestations,
        n_chains=2, warmup=20, iterations=60,
    )
    checks = verify_prediction(out, dataset, probe_characters=3)
    assert all(check["ok"] for check in checks), checks

    path = out / "concept/dated.json"
    payload = json.loads(path.read_text())
    payload["log_scores"][0] -= 0.1
    path.write_text(json.dumps(payload))
    tampered = verify_prediction(out, dataset, probe_characters=3)
    assert any(c["check"].endswith("dated:predictive_replay") and not c["ok"] for c in tampered)


def test_interrupted_prediction_resumes_identically(tmp_path, monkeypatch):
    import pytest
    import diachron.ppc.heldout as producer

    dataset, truth = simulate_dataset(
        n_languages=6, n_concepts=15, root_age=6000, seed=227, name="predict-resume",
    )
    attestations = simulate_attestations(dataset, truth, n_records=10, seed=22)
    kwargs = dict(dataset=dataset, attestations=attestations, n_chains=1, warmup=5, iterations=20)
    original = producer.run_campaign
    calls = 0

    def interrupted(*args, **kw):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise InterruptedError("controlled interruption after a saved model")
        return original(*args, **kw)

    monkeypatch.setattr(producer, "run_campaign", interrupted)
    out = tmp_path / "interrupted"
    with pytest.raises(InterruptedError):
        producer.run_heldout_prediction(str(out), **kwargs)
    assert json.loads((out / "PROGRESS.json").read_text())["completed_models"] == ["tree_only"]
    monkeypatch.setattr(producer, "run_campaign", original)
    result = producer.run_heldout_prediction(str(out), resume=True, **kwargs)
    fresh = tmp_path / "fresh"
    assert result == producer.run_heldout_prediction(str(fresh), **kwargs)
    for kind in ("concept", "language"):
        for model in ("tree_only", "dated", "geography", "borrowing"):
            assert json.loads((out / kind / model / "bundle.json").read_text()) == json.loads(
                (fresh / kind / model / "bundle.json").read_text()
            )
    with pytest.raises(ValueError, match="budget changed"):
        producer.run_heldout_prediction(str(out), resume=True, **{**kwargs, "iterations": 30})
