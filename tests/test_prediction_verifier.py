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
