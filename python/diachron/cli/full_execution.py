"""Run each full-profile workload with durable command evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time

from ..util.gitmeta import source_commit, working_tree_clean


def _write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def execute_full_profile(root: Path, inputs: Path, jobs: int, release_policy: Path | None):
    from diachron import _core
    from ..io.epidoc import load_edh
    from ..io.freeze import load_study

    if not 1 <= jobs <= 4 or not working_tree_clean():
        raise ValueError("full execution requires clean committed source and one to four workers")
    commit = source_commit()
    inputs = inputs.resolve()
    root = root.resolve()
    paths = {"iecor": str(inputs / "iecor"), "abvd": str(inputs / "abvd"),
             "catalog": str(inputs / "glottolog53.zip"), "edh_source": str(inputs / "edh"),
             "edh_geography": str(inputs / "edh-geography.json")}
    _, _, first = load_study("iecor", Path(paths["iecor"]), Path(paths["catalog"]))
    _, _, second = load_study("abvd", Path(paths["abvd"]), Path(paths["catalog"]))
    attestations = load_edh(Path(paths["edh_source"]), Path(paths["edh_geography"]))
    _write(root / "DATASETS.json", {"primary": first, "secondary": second,
                                  "edh": {key: value for key, value in attestations.items() if key != "attestations"}})
    _write(root / "ENVIRONMENT.json", {"platform": platform.platform(), "python": sys.version,
                                       "workers": jobs, "source_commit": commit})
    files = subprocess.check_output(["git", "ls-files", "-z"], text=True).split("\0")
    _write(root / "SOURCE_INPUTS.json", {"source_commit": commit,
        "native_binary_sha256": hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
        "files": {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                  for name in files if name and Path(name).is_file()}})
    models = ("tree_only", "dated", "geography", "borrowing")
    components = root / "components"
    primary_args = ["--source", paths["iecor"], "--catalog", paths["catalog"],
                    "--edh-source", paths["edh_source"], "--edh-geography", paths["edh_geography"]]
    budgets = ["--chains", "4", "--warmup", "3000", "--iterations", "12000", "--jobs", str(jobs)]
    plan = {"schema": 1, "source_commit": commit, "inputs": paths, "execution_mode": "full",
            "execution_record": str(root / "EXECUTION.json"),
            "recovery": {"source_commit": commit, "summary": str(components / "recovery/recovery-summary.json"),
                         "rows": str(components / "recovery/rows.jsonl.gz")},
            "workbench": {"html": str(components / "workbench/Diachron.html"),
                          "qa": str(components / "workbench/qa/workbench-qa.json")},
            "release": {"demo_directory": str(components / "demo"),
                        "private_policy": str(release_policy.resolve()) if release_policy else None,
                        "private_policy_sha256": hashlib.sha256(release_policy.read_bytes()).hexdigest()
                        if release_policy else None}}
    for name in ("inference", "prediction", "sensitivity"):
        plan[name] = {"directory": str(components / name), "source_commit": commit}
    _write(root / "FULL_PLAN.json", plan)
    reports = {}

    def run(name, arguments, timeout=18000):
        start = time.monotonic()
        with (root / f"{name}.stdout.log").open("w") as stdout, (root / f"{name}.stderr.log").open("w") as stderr:
            try:
                process = subprocess.Popen(arguments, stdout=stdout, stderr=stderr, start_new_session=True)
                try:
                    code = process.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    # A timed-out sampler can own multiple worker processes.
                    try:
                        import os
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                    code = 124
            except OSError as exc:
                stderr.write(str(exc) + "\n")
                code = 127
        record = {"label": name, "arguments": arguments, "source_commit": commit, "exit": code,
                  "seconds": time.monotonic() - start,
                  "max_child_rss_platform_units": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
        with (root / "COMMANDS.jsonl").open("a") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
        reports[name] = record
        _write(root / "EXECUTION.json", {"source_commit": commit, "status": "running", "commands": reports})
        return code

    python = [sys.executable, "-m", "diachron.cli.main"]
    run("inference", [*python, "inference-campaign", "--out", plan["inference"]["directory"],
                      "--jobs", str(jobs), *primary_args])
    run("prediction", [*python, "predict", "--study", "iecor", "--out", plan["prediction"]["directory"],
                       *budgets, *primary_args])
    run("sensitivity", [*python, "sensitivity", "--study", "iecor", "--out", plan["sensitivity"]["directory"],
                        *budgets, *primary_args])
    run("secondary", [*python, "sample", "--study", "abvd", "--model", "tree_only", "--profile", "full",
                      "--out", str(components / "secondary"), "--source", paths["abvd"],
                      "--catalog", paths["catalog"], *budgets])
    run("recovery", [*python, "recovery", "--datasets", "200", "--out", str(components / "recovery")])
    run("demo", [*python, "demo", "--out", plan["release"]["demo_directory"]], timeout=600)
    workbench = components / "workbench"
    workbench_data = workbench / "workbench.json"
    flags = [flag for model in models for flag in
             ("--bundle", f"{model}={components / 'inference' / model / model}")]
    run("export-workbench", [*python, "export-workbench", "--source", paths["iecor"], "--catalog", paths["catalog"],
                            "--out", str(workbench_data), *flags], timeout=600)
    run("build-workbench", ["npm", "--prefix", "ui", "run", "build"], timeout=300)
    run("package-workbench", [sys.executable, "scripts/bundle_workbench_html.py", "--data", str(workbench_data),
                             "--out", plan["workbench"]["html"]], timeout=120)
    run("browser", ["node", "ui/qa/workbench.mjs", plan["workbench"]["html"], str(workbench / "qa")], timeout=180)
    run("wheel", [str(Path(sys.executable).parent / "maturin"), "build", "--locked", "--release",
                  "--out", str(components / "wheels")], timeout=600)
    _write(root / "EXECUTION.json", {"source_commit": commit, "status": "completed", "commands": reports})
    return plan
