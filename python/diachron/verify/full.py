"""Evaluate the full profile from pinned inputs and completed component evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from ..io.types import DataError
from ..util.gitmeta import source_commit, working_tree_clean
from .verifier import _read_json, verify_evidence

_GATE_IDS = tuple(f"G{index:02d}" for index in range(1, 16))


def normalize_plan(path: Path) -> dict:
    plan = _read_json(path, 1024 * 1024)
    if plan.get("schema") != 1:
        raise ValueError("unsupported full evidence plan")
    paths = {
        "inputs": ("iecor", "abvd", "catalog", "edh_source", "edh_geography"),
        "inference": ("directory",), "prediction": ("directory",),
        "sensitivity": ("directory", "archive_directory", "asset_manifest"),
        "recovery": ("summary", "rows"), "workbench": ("html", "qa"),
        "release": ("demo_directory", "private_policy"),
    }
    for name, fields in paths.items():
        section = plan.get(name, {})
        if not isinstance(section, dict):
            raise ValueError(f"invalid plan section: {name}")
        for field in fields:
            if section.get(field) is not None:
                section[field] = str((path.parent / section[field]).resolve())
    plan["execution_mode"] = "audit_existing"
    plan.pop("execution_record", None)
    return plan


def _source(section):
    commit = section["source_commit"]
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("component requires a full trusted source commit")
    if subprocess.check_output(["git", "cat-file", "-t", commit], timeout=10).strip() != b"commit":
        raise ValueError("component source must identify a commit")
    return commit


def _tests(*names):
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", *names],
                            capture_output=True, text=True, timeout=180)
    return {"passed": result.returncode == 0, "exit": result.returncode,
            "files": list(names), "result": result.stdout.strip().splitlines()[-1:]}


def _inference_options(plan):
    inputs = plan["inputs"]
    return {"observed_source": Path(inputs["iecor"]), "catalog": Path(inputs["catalog"]),
            "edh_source": Path(inputs["edh_source"]), "edh_geography": Path(inputs["edh_geography"])}


def _check_fits(plan, component):
    section = plan[component]
    if component == "sensitivity" and "archive_directory" in section:
        from .sensitivity_archives import verify_sensitivity_archives

        result = verify_sensitivity_archives(Path(section["archive_directory"]), Path(section["asset_manifest"]),
                                             _source(section), **_inference_options(plan))
        summary = result.pop("summary")
        result["all_converged"] = summary["all_converged"]
        result["material_root_shift_settings"] = summary["levers_that_shift_root_age"]
        return result
    directory = Path(section["directory"])
    expected = _source(section)
    hosted = None
    if component == "prediction" and (directory / "JOB.json").is_file():
        from .hosted import hosted_execution

        hosted = hosted_execution(directory, expected, component)
        directory /= "fit"
    marker = {"inference": "RUN_PLAN.json", "prediction": "heldout.json", "sensitivity": "sensitivity.json"}[component]
    if (directory / "ACCEPTANCE.json").is_file() or not (directory / marker).is_file():
        raise ValueError(f"missing completed {component} component")
    report = verify_evidence(str(directory), expected_commit=expected, **_inference_options(plan))
    failed = [check["check"] for check in report["checks"] if not check["ok"]]
    detail = {"artifact_verification": report["ok"], "source_commit": expected,
              "bundles": report.get("n_bundles", 0), "failed_checks": failed}
    if component == "inference":
        frozen = any(check["check"] == "inference:frozen_plan" and check["ok"] for check in report["checks"])
        detail["passed"] = report["ok"] and frozen and detail["bundles"] == 4
        return detail
    if component == "prediction":
        from .campaign import _diagnose

        predictive = _read_json(directory / "heldout.json", 2 * 1024 * 1024)
        budget = predictive["budget"]
        budget_ok = (budget["chains"] >= 4 and budget["warmup"] >= 3000
                     and budget["iterations"] >= 12000 and budget["thin"] == 10)
        fits = {}
        if report["ok"]:
            for split in ("concept", "language"):
                for model in ("tree_only", "dated", "geography", "borrowing"):
                    raw = _read_json(directory / split / model / "bundle.json", 128 * 1024 * 1024)
                    fits[f"{split}/{model}"] = _diagnose(raw)["converged"]
        detail.update({"full_budget": budget_ok, "fit_convergence": fits,
                       "passed": report["ok"] and budget_ok and len(fits) == 8 and all(fits.values())})
        if hosted is not None:
            detail["hosted_execution"] = hosted
        return detail
    sensitivity = _read_json(directory / "sensitivity.json", 8 * 1024 * 1024)
    budget = _read_json(directory / "SENSITIVITY_PLAN.json", 2 * 1024 * 1024)
    budget_ok = budget["chains"] >= 4 and budget["warmup"] >= 3000 and budget["iterations"] >= 12000
    detail.update({"full_budget": budget_ok, "settings": sensitivity["n_configs"],
                   "all_converged": sensitivity["all_converged"],
                   "passed": report["ok"] and budget_ok and sensitivity["n_configs"] >= 12
                   and sensitivity["all_converged"]})
    return detail


def _observations(plan):
    from ..io.epidoc import load_edh
    from ..io.freeze import load_study

    inputs = _inference_options(plan)
    primary, _, primary_report = load_study("iecor", inputs["observed_source"], inputs["catalog"])
    _, _, second = load_study("abvd", Path(plan["inputs"]["abvd"]), inputs["catalog"])
    edh = load_edh(inputs["edh_source"], inputs["edh_geography"])
    return primary, primary_report, second, edh


def _evidence_profile(primary, primary_report, second, edh):
    records = edh["attestations"]
    qualified = [row for row in records if row["source_ref"] and row["bibliography"] and row["region"]
                 and row["age_low"] <= row["age_high"] and row["location_uncertainty_note"]]
    alternatives = sum(record.proposed_cognate_id is not None for record in primary.cognate_classes)
    passed = (not primary.is_synthetic and primary_report["eligible_languages"] >= 100
              and primary_report["eligible_concepts"] >= 100 and second["eligible_languages"] >= 30
              and len({row["source_ref"] for row in qualified}) >= 100 and alternatives > 0)
    return {"passed": passed, "primary_varieties": primary_report["eligible_languages"],
            "primary_concepts": primary_report["eligible_concepts"], "secondary_varieties": second["eligible_languages"],
            "attestations": len(qualified), "proposed_cognacy_links": alternatives}


def _workbench(plan):
    section = plan["workbench"]
    html = Path(section["html"])
    if html.stat().st_size > 20_000_000:
        raise ValueError("oversized standalone workbench")
    content = html.read_bytes()
    qa = _read_json(Path(section["qa"]), 1024 * 1024)
    raw = re.search(rb'<script id="workbench-data" type="application/json">(.*?)</script>', content, re.DOTALL)
    if raw is None or qa["html_sha256"] != hashlib.sha256(content).hexdigest():
        raise ValueError("workbench snapshot or browser evidence identity differs")
    return json.loads(raw[1]), qa


def _geography(plan, primary, edh):
    from ..cli.workbench import _clade_locations, _retained_subset
    from ..infer.artifacts import PosteriorBundle

    snapshot, _ = _workbench(plan)
    reference = {(language.id, language.latitude, language.longitude) for language in primary.languages
                 if not language.is_historical and language.latitude is not None and language.longitude is not None}
    points = snapshot["modernPoints"]
    modern = {(point["label"], point["lat"], point["lon"]) for point in points}
    historical = {row["id"]: row for row in edh["attestations"]}
    shown = snapshot["attestations"]
    records_match = len(shown) == len(historical) and len({point["id"] for point in shown}) == len(shown)
    for point in shown:
        original = historical[point["id"]]
        records_match &= (point["kind"] == "attestation" and point["sourceRef"] == original["source_ref"]
                          and point["lat"] == original["latitude"] and point["lon"] == original["longitude"]
                          and point["sourceEpoch"] == 1950 and point["ageLow"] == original["age_low"] + 50
                          and point["ageHigh"] == original["age_high"] + 50)
    bundle = PosteriorBundle.load(str(Path(plan["inference"]["directory"]) / "geography/geography"))
    roots, clades = _clade_locations(bundle, _retained_subset(bundle))
    model = snapshot["bundles"]["geography"]
    locations_match = (bool(roots) and model["source_commit"] == bundle.source_commit
                       and model["rootLocations"] == roots and model["cladeLocations"] == clades)
    return {"passed": modern == reference and all(point["kind"] == "modern-reference" for point in points)
            and records_match and locations_match, "modern_points": len(points),
            "attested_records": len(shown), "recomputed_root_locations": len(roots),
            "recomputed_clade_clouds": len(clades)}


def audit_full_plan(plan: dict, expected_commit: str, out: Path) -> dict:
    """Re-run each available gate while retaining incomplete and failed checks."""
    trusted = (plan.get("schema") == 1 and plan.get("source_commit") == expected_commit
               and expected_commit == source_commit() and working_tree_clean())
    if not trusted:
        return {"status": "source_mismatch", "full_profile_executed": False, "source_commit": source_commit(),
                "profile": "full", "ready": False, "gates": {}, "blocking_gates": list(_GATE_IDS)}
    from ..cli.checks import run_likelihood_oracle
    from ..cli.thesis import build_thesis
    from .tiny_posterior import run_tiny_posterior_audit
    from .recovery import verify_recovery

    execution = None
    if plan.get("execution_mode") == "full":
        execution = _read_json(Path(plan["execution_record"]), 2 * 1024 * 1024)
        if execution["source_commit"] != expected_commit or execution["status"] != "completed":
            execution = None
    gates = {}

    def gate(name, function):
        try:
            detail = function()
            gates[name] = {"status": "passed" if detail["passed"] else "failed", "evidence": detail}
        except (OSError, KeyError, ValueError, TypeError, IndexError, DataError, subprocess.SubprocessError) as exc:
            gates[name] = {"status": "unavailable", "reason": str(exc)}

    observations = None

    def evidence():
        nonlocal observations
        observations = _observations(plan)
        return _evidence_profile(*observations)

    gate("G01", evidence)
    gate("G02", lambda: run_likelihood_oracle(n_cases=1200))

    def sampler():
        oracle = run_tiny_posterior_audit()
        moves = _tests("tests/test_mcmc_exact.py", "tests/test_moves.py")
        return {"passed": oracle["passed"] and moves["passed"], "finite_targets": oracle, "moves": moves}

    gate("G03", sampler)
    gate("G04", lambda: _tests("tests/test_ascertainment.py"))

    def recovery():
        section = plan["recovery"]
        detail = verify_recovery(Path(section["summary"]), Path(section["rows"]), _source(section))
        detail["passed"] = detail["artifact_valid"] and detail["failed"] == 0 and detail["calibration_target_met"]
        return detail

    gate("G05", recovery)
    for name, component in (("G06", "inference"), ("G07", "prediction"), ("G08", "sensitivity")):
        gate(name, lambda component=component: _check_fits(plan, component))
    gate("G09", lambda: _geography(plan, observations[0], observations[3]))

    def borrowing():
        from ..infer.campaign import prepare_study

        dataset = observations[0]
        study = prepare_study(dataset)
        example = {key: states for key, states in zip(study.character_ids, study.columns)
                   if key in {("61", "4364"), ("61", "227")}}
        annotations = sum(judgment.cognate_class == "4364" and judgment.is_loan for judgment in dataset.judgments)
        tests = _tests("tests/test_geo_borrowing.py")
        return {"passed": tests["passed"] and len(example) == 2 and annotations > 0,
                "contact_model_checks": tests, "annotated_fruit_forms": annotations,
                "case_class_presence_counts": {key[1]: values.count(1) for key, values in example.items()}}

    gate("G10", borrowing)

    def workbench():
        snapshot, qa = _workbench(plan)
        steps = {step["id"] for step in snapshot["tour"]}
        source = plan["inference"]["source_commit"]
        passed = (qa["passed"] and qa["network_requests"] == 0 and qa["mobile"] == [390, 844]
                  and steps == {"missing", "cognacy", "borrowing", "dating", "geography"}
                  and bool(snapshot["taxonomy"]) and set(snapshot["bundles"]) == {"tree_only", "dated", "geography", "borrowing"}
                  and all(bundle["source_commit"] == source for bundle in snapshot["bundles"].values()))
        return {"passed": bool(passed), "html_sha256": qa["html_sha256"], "tour_steps": sorted(steps),
                "catalog_groups": len(snapshot["taxonomy"]), "browser": qa["browser"]}

    gate("G11", workbench)
    gate("G12", lambda: _tests("tests/test_recovery.py", "tests/test_replica.py"))

    def controls():
        from .controls import run_semantic_controls

        section = plan["inference"]
        options = _inference_options(plan)
        options["source"] = options.pop("observed_source")
        return run_semantic_controls(Path(section["directory"]) / "dated/dated", out / "controls",
                                     _source(section), **options)

    gate("G13", controls)

    def thesis():
        result = build_thesis()
        return {**result, "passed": bool(result["ledger_ok"] and result["rendered_html"])}

    gate("G14", thesis)

    def publication():
        from ..cli.release import run_release
        from .release_policy import private_policy_passes

        section = plan["release"]
        commands = execution["commands"]
        expected_stages = {"inference", "prediction", "sensitivity", "secondary", "recovery", "demo",
                           "export-workbench", "build-workbench", "package-workbench", "browser", "wheel"}
        commands_ok = (set(commands) == expected_stages and all(
            command["exit"] == 0 and command["seconds"] > 0 and command["source_commit"] == expected_commit
            for command in commands.values()
        ))
        policy_ok = False
        if section.get("private_policy"):
            policy_path = Path(section["private_policy"])
            if policy_path.stat().st_size <= 65536:
                policy_ok = (hashlib.sha256(policy_path.read_bytes()).hexdigest() == section["private_policy_sha256"]
                             and private_policy_passes(policy_path, expected_commit))
        demo_path = Path(section["demo_directory"])
        if (demo_path / "ACCEPTANCE.json").is_file():
            raise ValueError("the demonstration must be a standalone component")
        from .campaign import _source_matches

        identities_ok = _source_matches(Path(plan["execution_record"]).parent, expected_commit)
        source_consistent = all(plan[name]["source_commit"] == expected_commit
                                for name in ("inference", "prediction", "sensitivity", "recovery"))
        demo = verify_evidence(str(demo_path), expected_commit=expected_commit)
        remote = subprocess.check_output(["git", "ls-remote", "origin", "refs/heads/main"], text=True, timeout=30)
        hosted = json.loads(subprocess.check_output([
            "gh", "run", "list", "--workflow", "ci.yml", "--commit", expected_commit, "--limit", "5",
            "--json", "headSha,conclusion,url",
        ], text=True, timeout=30))
        exact_ci = [run["url"] for run in hosted if run["headSha"] == expected_commit and run["conclusion"] == "success"]
        package = run_release()
        return {"passed": bool(commands_ok and policy_ok and identities_ok and source_consistent and demo["ok"] and remote.split()[0] == expected_commit and exact_ci
                               and package["source_package_checks_passed"]), "demo_replay": demo["ok"],
                "hosted_ci": exact_ci, "source_package": package["source_package_checks_passed"],
                "execution_records_passed": commands_ok, "release_policy_passed": policy_ok,
                "source_inventory_passed": identities_ok and source_consistent}

    if execution is not None and all(result["status"] == "passed" for result in gates.values()):
        gate("G15", publication)
    else:
        reason = "required evidence gates remain open" if execution else "full execution record required"
        gates["G15"] = {"status": "blocked", "reason": reason}
    pending = [name for name in _GATE_IDS if gates[name]["status"] != "passed"]
    return {"profile": "full", "source_commit": expected_commit, "full_profile_executed": execution is not None,
            "execution_mode": plan.get("execution_mode", "audit_existing"),
            "status": "passed" if not pending else "incomplete", "ready": not pending,
            "blocking_gates": pending, "gates": gates}


def verify_full_acceptance(root: Path, expected_commit: str | None) -> dict:
    if expected_commit is None:
        return {"ok": False, "checks": [], "reason": "provide a trusted full profile source commit"}
    recorded = _read_json(root / "ACCEPTANCE.json", 4 * 1024 * 1024)
    plan = _read_json(root / "FULL_PLAN.json", 1024 * 1024)
    with tempfile.TemporaryDirectory(prefix="diachron-reaudit-") as scratch:
        reproduced = audit_full_plan(plan, expected_commit, Path(scratch))
    def stable(value):
        if isinstance(value, dict):
            return {key: stable(item) for key, item in value.items()
                    if key not in {"seconds", "result", "hosted_ci"}}
        if isinstance(value, list):
            return [stable(item) for item in value]
        return value

    agreement = (recorded["source_commit"] == expected_commit and recorded["status"] == reproduced["status"]
                 and recorded["blocking_gates"] == reproduced["blocking_gates"]
                 and recorded["full_profile_executed"] == reproduced["full_profile_executed"]
                 and stable(recorded["gates"]) == stable(reproduced["gates"]))
    return {"ok": agreement and reproduced["ready"], "checks": [
        {"check": "full_gate_recomputation", "ok": agreement, "detail": reproduced["blocking_gates"]},
        {"check": "full_acceptance", "ok": reproduced["ready"], "detail": reproduced["status"]},
    ], "gates": reproduced["gates"]}
