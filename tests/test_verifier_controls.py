"""Gate G13: project-specific semantic negative controls.

Each control makes a structurally valid copy of an honest evidence bundle, alters
a meaningful result or relationship, rebuilds the ordinary hashes, and then
invokes the real verifier (as a subprocess) which must reject for a semantic
reason. We record the actual exit status and the failing check. A control that
the verifier accepted would be a failure of the verifier.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from diachron.io.simulate import simulate_dataset, simulate_attestations
from diachron.infer.campaign import prepare_study, run_campaign, campaign_bundle
from diachron.models.posterior import ModelConfig
from diachron.verify import controls


def _honest_bundle(tmp_path):
    ds, truth = simulate_dataset(
        n_languages=8, n_concepts=40, root_age=6000.0,
        borrow_prob=0.2, missing_frac=0.08, seed=3, name="diachron_primary",
    )
    atts = simulate_attestations(ds, truth, n_records=20, seed=3)
    study = prepare_study(ds, atts)
    chains = run_campaign(
        study, ModelConfig.dated(), n_chains=4, warmup=200, iterations=600,
        thin=10, root_age_init=truth.root_age, base_seed=900,
    )
    bundle = campaign_bundle(chains, study, ModelConfig.dated(), source_commit="t")
    src = str(tmp_path / "honest")
    bundle.save(src)
    return src


def _run_verifier(evidence_dir: str) -> dict:
    """Invoke the real verifier as a subprocess and return its JSON result."""
    proc = subprocess.run(
        [sys.executable, "-m", "diachron.cli.main", "verify", "--evidence", evidence_dir],
        capture_output=True, text=True, check=False, timeout=30,
    )
    result = json.loads(proc.stdout)
    assert proc.returncode == (0 if result["ok"] else 1), proc.stderr
    return result


def test_honest_bundle_accepted(tmp_path):
    src = _honest_bundle(tmp_path)
    # Verifier expects an evidence dir containing bundle dirs; pass the parent.
    result = _run_verifier(str(tmp_path))
    assert result["ok"], result


@pytest.mark.parametrize(
    "control_fn,expected_substring",
    [
        (controls.control_fabricated_likelihood, "posterior_replay"),
        (controls.control_rehashed_altered_result, "attestation_intervals"),
        (controls.control_missing_to_absent, "posterior_replay"),
        (controls.control_modern_as_ancient, "geographic_honesty"),
        (controls.control_corrupted_digest, "draw_digest"),
    ],
)
def test_control_rejected_for_semantic_reason(tmp_path, control_fn, expected_substring):
    src = _honest_bundle(tmp_path)
    work = str(tmp_path / "work")
    os.makedirs(work, exist_ok=True)
    tampered_dir = control_fn(src, work)
    # The verifier walks for bundle.json, so point it at the tampered bundle's parent.
    parent = os.path.dirname(tampered_dir)
    result = _run_verifier(parent)
    assert not result["ok"], f"verifier accepted a tampered bundle: {result}"
    failing = [c for c in result["checks"] if not c["ok"]]
    assert any(expected_substring in c["check"] for c in failing), (
        f"expected a {expected_substring} rejection, got {failing}"
    )
