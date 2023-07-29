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
