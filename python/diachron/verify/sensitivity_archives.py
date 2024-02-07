"""Verify all archived sensitivity shards with one restored fit at a time."""

from __future__ import annotations

from pathlib import Path
import tempfile

from ..io.evidence_archive import restore_evidence
from ..io.types import DataError
from ..sensitivity.run import _declared_configs, summarize_sensitivity
from ..util.gitmeta import source_commit
from .verifier import _read_json, verify_evidence


def verify_sensitivity_archives(directory: Path, asset_manifest: Path, expected_commit: str, **inputs) -> dict:
    manifest = _read_json(asset_manifest, 1024 * 1024)
    if manifest["source_commit"] != expected_commit:
        raise ValueError("archive inventory identifies a different source")
    assets = manifest["assets"]
    rows = []
    records = {}
    identity = None
    for setting in _declared_configs():
        name = setting.name
        filename = f"{name}.tar.gz"
        if filename not in assets or not (directory / filename).is_file():
            records[name] = {"status": "unavailable", "reason": "missing archive"}
            continue
        try:
            with tempfile.TemporaryDirectory(dir=directory, prefix=".sensitivity-replay-") as temporary:
                restored = Path(temporary) / "restored"
                result = restore_evidence(directory / filename, restored, assets[filename], max_bytes=256 * 1024 * 1024)
                job = _read_json(restored / "JOB.json", 1024 * 1024)
                fit = restored / "fit"
                plan = _read_json(fit / "SENSITIVITY_PLAN.json", 2 * 1024 * 1024)
                common = {key: value for key, value in plan.items() if key not in {"shard", "configurations"}}
                if identity is not None and common != identity:
                    raise ValueError("archived fits do not share the same frozen inputs and budgets")
                if (job["source_commit"] != expected_commit or job["configuration"] != name
                        or job["status"] != "executed_and_verified" or plan["shard"] != name
                        or [record["label"] for record in job["commands"]] != ["acquire", "fit", "verify"]
                        or any(record["exit"] != 0 for record in job["commands"])):
                    raise ValueError("archive job execution did not complete")
                verification = verify_evidence(str(fit), expected_commit=expected_commit, **inputs)
                if not verification["ok"] or verification["n_bundles"] != 1:
                    raise ValueError(f"shard replay failed: {[c['check'] for c in verification['checks'] if not c['ok']]}")
                report = _read_json(fit / "sensitivity.json", 8 * 1024 * 1024)
                if report["is_synthetic"] or report["study"] != "iecor" or report["n_configs"] != 1:
                    raise ValueError("archive is outside the observed sensitivity profile")
                if identity is None:
                    identity = common
                row = report["configs"][0]
                rows.append(row)
                records[name] = {"status": "verified", "archive_sha256": result["archive_sha256"],
                                 "restored_files": result["restored_files"], "converged": row["converged"]}
        except (DataError, OSError, ValueError, KeyError, TypeError, IndexError) as exc:
            records[name] = {"status": "failed", "reason": str(exc)}
    complete = len(rows) == len(records)
    full_budget = bool(identity and identity["chains"] >= 4 and identity["warmup"] >= 3000
                       and identity["iterations"] >= 12000 and identity["thin"] == 10)
    summary = summarize_sensitivity(rows, "iecor", False, expected_commit)
    return {"source_commit": expected_commit, "verifier_source_commit": source_commit(),
            "artifact_verification": complete,
            "full_budget": full_budget, "completed": len(rows), "required": len(records),
            "passed": complete and full_budget and summary["all_converged"],
            "archives": records, "summary": summary}
