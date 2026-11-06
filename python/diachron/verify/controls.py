"""Make independent rehashed edits to the saved evidence for negative controls."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil


def _rehash(value) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _edit(src_dir, work_dir, name, change, *, rehash=True):
    out = Path(work_dir) / name
    shutil.copytree(src_dir, out)
    artifact = out / "bundle.json"
    manifest_path = out / "MANIFEST.json"
    bundle = json.loads(artifact.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    change(bundle)
    if rehash:
        manifest["draw_digest"] = _rehash(bundle["chains"])
        manifest["context_sha256"] = _rehash(bundle["context"])
    artifact.write_text(json.dumps(bundle), encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return str(out)


def control_fabricated_likelihood(src_dir, work_dir):
    def change(bundle):
        bundle["chains"][0][0]["log_posterior"] -= 0.25

    return _edit(src_dir, work_dir, "wrong_posterior", change)


def control_rehashed_altered_result(src_dir, work_dir):
    def change(bundle):
        records = bundle["context"]["attestations"]
        oldest = max(records, key=lambda att: att["age_high"])
        oldest["age_high"] += 15

    return _edit(src_dir, work_dir, "altered_attestation_interval", change)


def control_missing_to_absent(src_dir, work_dir):
    def change(bundle):
        for column in bundle["context"]["model"]["columns"]:
            if -1 in column:
                column[column.index(-1)] = 0
                return
        raise ValueError("control requires a missing cell")

    return _edit(src_dir, work_dir, "missing_to_absent", change)


def control_modern_as_ancient(src_dir, work_dir):
    def change(bundle):
        attestation = bundle["context"]["attestations"][0]
        attestation["location_kind"] = "modern-reference"
        latitude, longitude = next(iter(bundle["context"]["model"]["tip_coords"].values()))
        attestation["latitude"] = latitude
        attestation["longitude"] = longitude

    return _edit(src_dir, work_dir, "modern_as_ancient", change)


def control_corrupted_digest(src_dir, work_dir):
    def change(bundle):
        bundle["chains"][0][0]["gain"] *= 3

    return _edit(src_dir, work_dir, "corrupted_digest", change, rehash=False)


def run_semantic_controls(bundle, out, expected_commit, *, source=None, catalog=None,
                          edh_source=None, edh_geography=None):
    """Execute the real verifier CLI on an intact bundle and four rehashed edits."""
    import subprocess
    import sys
    import tempfile

    from ..util.gitmeta import source_commit

    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    arguments = ["--expected-commit", expected_commit]
    for flag, path in (("--source", source), ("--catalog", catalog),
                       ("--edh-source", edh_source), ("--edh-geography", edh_geography)):
        if path is not None:
            arguments.extend([flag, str(path)])

    def execute(label, directory):
        command = [sys.executable, "-m", "diachron.cli.main", "verify", "--evidence", str(directory), *arguments]
        result = subprocess.run(command, capture_output=True, text=True, timeout=180)
        (out / f"{label}.stdout.json").write_text(result.stdout)
        (out / f"{label}.stderr.log").write_text(result.stderr)
        try:
            payload = json.loads(result.stdout)
        except ValueError:
            payload = {"ok": False, "checks": []}
        return result.returncode, payload

    code, original = execute("original", bundle)
    positive = code == 0 and original["ok"]
    results = []
    for name, edit, expected_check in (
        ("wrong_posterior", control_fabricated_likelihood, "posterior_replay"),
        ("altered_attestation_interval", control_rehashed_altered_result,
         "frozen_observed_inputs" if source else "attestation_intervals"),
        ("missing_to_absent", control_missing_to_absent,
         "frozen_observed_inputs" if source else "posterior_replay"),
        ("modern_as_ancient", control_modern_as_ancient, "geographic_honesty"),
    ):
        with tempfile.TemporaryDirectory(dir=out, prefix=".control-") as work:
            changed = edit(bundle, work)
            code, verification = execute(name, changed)
            checks = {check["check"].rsplit(":", 1)[-1]: check["ok"] for check in verification["checks"]}
            rejected = (code == 1 and not verification["ok"] and checks.get(expected_check) is False
                        and checks.get("draw_digest") is True and checks.get("context_digest") is True
                        and checks.get("source_commit") is True)
            results.append({"control": name, "verifier_exit": code, "expected_check": expected_check,
                            "rejected_semantically": rejected})
    report = {"source_commit": expected_commit, "verifier_source_commit": source_commit(),
              "positive_control_passed": positive, "controls": results,
              "passed": positive and all(row["rejected_semantically"] for row in results)}
    (out / "SEMANTIC_CONTROLS.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
