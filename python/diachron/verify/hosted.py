"""Validate hosted job execution records before replaying their result files."""

from __future__ import annotations

import math
from pathlib import Path
import re

from .campaign import _source_matches
from .verifier import _read_json


def hosted_execution(root: Path, commit: str, kind: str, configuration: str | None = None) -> dict:
    job = _read_json(root / "JOB.json", 1024 * 1024)
    environment = _read_json(root / "ENVIRONMENT.json", 1024 * 1024)
    source = _read_json(root / "SOURCE_INPUTS.json", 2 * 1024 * 1024)
    steps = job["commands"]
    if (job["source_commit"] != commit or job["kind"] != kind or job["configuration"] != configuration
            or job["status"] != "executed_and_verified" or environment["source_commit"] != commit
            or not re.fullmatch(r"[0-9]+", environment.get("hosted_run") or "")
            or not re.fullmatch(r"[0-9a-f]{64}", source["native_binary_sha256"])
            or [step["label"] for step in steps] != ["acquire", "fit", "verify"]):
        raise ValueError("hosted job identity or execution status differs")
    subcommands = ("acquire", "predict" if kind == "prediction" else "sensitivity", "verify")
    for step, command in zip(steps, subcommands):
        if (step["source_commit"] != commit or step["exit"] != 0 or not math.isfinite(step["seconds"])
                or not 0 < step["seconds"] <= 18100 or step["arguments"][1:4] != ["-m", "diachron.cli.main", command]):
            raise ValueError(f"hosted {command} did not record a successful bounded execution")
    if not _source_matches(root, commit):
        raise ValueError("hosted source inventory differs from trusted Git bytes")
    return {"run_id": environment["hosted_run"], "attempt": environment["attempt"],
            "native_binary_sha256": source["native_binary_sha256"], "fit_seconds": steps[1]["seconds"]}
