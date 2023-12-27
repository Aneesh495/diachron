"""Execute and resume a source-bound observed inference campaign."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

from ..util.gitmeta import source_commit, working_tree_clean

_MODELS = ("tree_only", "dated", "geography", "borrowing")


def _write(path: Path, payload) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def run_inference_campaign(
    out: Path, source: Path, catalog: Path,
    edh_source: Path | None, edh_geography: Path | None, jobs: int, resume: bool,
) -> int:
    if not 1 <= jobs <= 4:
        raise ValueError("inference campaign supports one to four workers")
    if bool(edh_source) != bool(edh_geography):
        raise ValueError("provide both EDH inputs")
    if not working_tree_clean():
        raise ValueError("campaign requires committed source in a clean checkout")
    commit = source_commit()
    plan = {
        "schema": 1, "component": "observed_inference", "protocol": 3,
        "source_commit": commit, "models": list(_MODELS),
        "study": "iecor", "source": str(source.resolve()), "catalog": str(catalog.resolve()),
        "edh_source": str(edh_source.resolve()) if edh_source else None,
        "edh_geography": str(edh_geography.resolve()) if edh_geography else None,
        "chains": 4, "warmup": 3000, "iterations": 12000, "thin": 10,
        "initialization": "distance", "jobs": jobs,
    }
    plan_path = out / "RUN_PLAN.json"
    if out.exists():
        if not resume or not plan_path.is_file():
            raise FileExistsError("campaign output exists; use --resume with the same source and plan")
        if json.loads(plan_path.read_text()) != plan:
            raise ValueError("resumed campaign plan or executable source differs")
    else:
        out.mkdir(parents=True)
        _write(plan_path, plan)
        files = subprocess.check_output(["git", "ls-files", "-z"]).decode().split("\0")
        _write(out / "SOURCE_INPUTS.json", {
            "source_commit": commit,
            "files": {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                      for name in files if name and Path(name).is_file()},
        })
        _write(out / "ENVIRONMENT.json", {
            "python": sys.version, "platform": platform.platform(),
            "executable": sys.executable, "processors_used_at_most": jobs,
        })

    records = {}
    status_path = out / "INFERENCE.json"
    if status_path.is_file():
        records = json.loads(status_path.read_text()).get("models", {})
    commands = out / "COMMANDS.jsonl"

    def execute(arguments: list[str], label: str) -> int:
        started = time.time()
        command = [sys.executable, "-m", "diachron.cli.main", *arguments]
        with (out / f"{label}.stdout.json").open("w") as stdout, (
            out / f"{label}.stderr.log"
        ).open("w") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr)
        record = {
            "arguments": command, "label": label, "exit": result.returncode,
            "seconds": time.time() - started, "source_commit": commit,
        }
        with commands.open("a", encoding="utf-8") as log:
            log.write(json.dumps(record, allow_nan=False) + "\n")
            log.flush()
        return result.returncode

    for name in _MODELS:
        sample_args = [
            "sample", "--study", "iecor", "--model", name, "--profile", "full",
            "--chains", "4", "--jobs", str(jobs), "--warmup", "3000",
            "--iterations", "12000", "--thin", "10", "--initialization", "distance",
            "--out", str(out / name), "--source", str(source), "--catalog", str(catalog),
        ]
        if edh_source:
            sample_args.extend(["--edh-source", str(edh_source), "--edh-geography", str(edh_geography)])
        if resume:
            sample_args.append("--resume")
        code = execute(sample_args, f"{name}-sample")
        record = {"sample_exit": code, "status": "failed", "converged": False}
        if code == 0:
            bundle = out / name / name
            diagnose_code = execute(["diagnose", "--bundle", str(bundle)], f"{name}-diagnostics")
            verify_args = [
                "verify", "--evidence", str(bundle), "--expected-commit", commit,
                "--source", str(source), "--catalog", str(catalog),
            ]
            if edh_source:
                verify_args.extend(["--edh-source", str(edh_source), "--edh-geography", str(edh_geography)])
            verify_code = execute(verify_args, f"{name}-verify")
            record.update({"diagnostic_exit": diagnose_code, "verify_exit": verify_code})
            if diagnose_code == verify_code == 0:
                diagnostic = json.loads((out / f"{name}-diagnostics.stdout.json").read_text())
                record.update({"status": "executed", "converged": diagnostic["converged"]})
        records[name] = record
        _write(status_path, {
            "source_commit": commit, "component": "observed_inference", "models": records,
            "status": "completed" if len(records) == len(_MODELS) else "running",
            "all_converged": len(records) == len(_MODELS) and all(r["converged"] for r in records.values()),
            "full_project_acceptance": False,
        })
    return 0 if all(record["converged"] for record in records.values()) else 1
