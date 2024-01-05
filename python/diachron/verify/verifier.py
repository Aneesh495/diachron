"""Check saved posterior structure and replay a bounded set of exact draws."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

from ..models.posterior import Calibration, ModelConfig, PosteriorModel
from ..models.tree import PhyloTree
from ..models.tip_dates import TipDate


def _sha256(value) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def _model(raw: dict) -> PosteriorModel:
    kwargs = dict(raw)
    kwargs["config"] = ModelConfig(**raw["config"])
    kwargs["calibrations"] = [
        Calibration(**{**c, "clade_taxa": tuple(c["clade_taxa"])}) for c in raw["calibrations"]
    ]
    kwargs["tip_coords"] = {int(k): tuple(v) for k, v in raw["tip_coords"].items()}
    kwargs["tip_dates"] = {int(k): TipDate(**v) for k, v in raw.get("tip_dates", {}).items()}
    kwargs["tip_regions"] = {int(k): v for k, v in raw["tip_regions"].items()}
    if kwargs["region_model"] is not None:
        from ..geo.regions import RegionModel

        kwargs["region_model"] = RegionModel(**raw["region_model"])
    return PosteriorModel(**kwargs)


def _validate_geography(bundle: dict) -> bool:
    if "ancient_coordinates" in bundle:
        return False
    context = bundle["context"]
    ids = {a["id"] for a in context["attestations"]}
    if len(ids) != len(context["attestations"]):
        return False
    for attestation in context["attestations"]:
        source = attestation["source_ref"]
        if not source or attestation.get("location_kind") in {"modern-reference", "illustrative-spread"}:
            return False
        if source.startswith("synthetic:") and not bundle["is_synthetic"]:
            return False
        if not bundle["is_synthetic"] and attestation.get("latitude") is not None:
            if attestation.get("location_kind") not in {"reported_origin", "reported_findspot"}:
                return False
    return True


def _calibrations_match(bundle: dict) -> bool:
    raw = bundle["context"]["model"]
    attestations = bundle["context"]["attestations"]
    for att in attestations:
        lo, hi = att["age_low"], att["age_high"]
        if not (math.isfinite(lo) and math.isfinite(hi) and 0 <= lo <= hi):
            return False
    derived = [
        cal for cal in raw["calibrations"] if cal["source_ref"] == "aggregate of dated attestations"
    ]
    if not derived:
        return True
    known = set(bundle["tip_labels"])
    candidates = [
        att["age_high"] for att in attestations if att["language_id"] in known
        and att["language_assignment_confidence"] is not None
        and att["language_assignment_confidence"] >= 0.5
    ]
    if not candidates or len(derived) != 1:
        return False
    oldest = max(candidates)
    cal = derived[0]
    return (
        cal["age_low"] == oldest and cal["age_high"] == oldest * 1.8 + 500.0
        and cal["clade_taxa"] == list(range(len(known))) and cal["kind"] == "uniform"
    )


def _draws_well_formed(bundle: dict) -> bool:
    chains = bundle["chains"]
    tips = bundle["tip_labels"]
    if not tips or len(set(tips)) != len(tips) or not chains or len(chains) != bundle["n_chains"]:
        return False
    if bundle["draws_per_chain"] < 1:
        return False
    for chain in chains:
        if len(chain) != bundle["draws_per_chain"]:
            return False
        last = -1
        for draw in chain:
            fields = ("root_age", "log_posterior", "gain", "loss", "clock_rate", "geo_sigma2")
            if any(not math.isfinite(draw[name]) for name in fields):
                return False
            if draw["iteration"] <= last or draw["root_age"] < 0:
                return False
            last = draw["iteration"]
    return True


def _replay(bundle: dict, probe_draws: int) -> tuple[bool, int]:
    model = _model(bundle["context"]["model"])
    labels = dict(enumerate(bundle["tip_labels"]))
    probes_per_chain = max(1, probe_draws // len(bundle["chains"]))
    count = 0
    for chain in bundle["chains"]:
        stride = max(1, len(chain) // probes_per_chain)
        for draw in chain[::stride][:probes_per_chain]:
            count += 1
            tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"], labels)
            if tree.root_age() != draw["root_age"] or tree.newick() != draw["newick"]:
                return False, count
            score = model.log_posterior(
                tree, draw["gain"], draw["loss"], draw["clock_rate"],
                draw["geo_sigma2"], draw["borrow_rates"], draw.get("clock_multipliers"),
            )
            if not math.isfinite(score) or abs(score - draw["log_posterior"]) > 1e-9:
                return False, count
    return count > 0, count


def _observed_context_matches(
    bundle: dict, source: Path, catalog: Path, edh_source: Path | None, edh_geography: Path | None,
) -> bool:
    from dataclasses import asdict, replace

    import numpy as np

    from ..infer.campaign import prepare_study, build_model
    from ..io.epidoc import load_edh
    from ..io.freeze import load_study
    from ..io.types import Attestation

    name = bundle["dataset_name"]
    dataset, _, manifest = load_study(name, source, catalog)
    observations = []
    if edh_source is not None:
        cohort = load_edh(edh_source, edh_geography)
        observations = cohort["attestations"]
    model = bundle["context"]["model"]
    partition = bundle["context"].get("held_out_partition")
    if partition is not None:
        seed = partition["seed"]
        if not isinstance(seed, int) or not 0 <= seed < 2 ** 32:
            return False
        rng = np.random.default_rng(seed)
        expected = {}
        for kind, identifiers in (
            ("concept", sorted(c.id for c in dataset.concepts)),
            ("language", sorted(lang.id for lang in dataset.languages)),
        ):
            expected[kind] = set(rng.choice(
                identifiers, size=max(1, len(identifiers) // 5), replace=False,
            ).tolist())
        kind = partition["kind"]
        held = expected[kind]
        if partition["identifiers"] != sorted(held):
            return False
        forms = [
            form for form in dataset.lexemes
            if (form.concept_id if kind == "concept" else form.language_id) not in held
        ]
        form_ids = {form.id for form in forms}
        dataset = replace(
            dataset, lexemes=forms,
            judgments=[judgment for judgment in dataset.judgments if judgment.lexeme_id in form_ids],
        )
    sensitivity = bundle["context"].get("sensitivity")
    if sensitivity is not None:
        if partition is not None:
            return False
        from ..sensitivity.interventions import prepare_sensitivity_study
        from ..sensitivity.run import _declared_configs

        settings = {setting.name: setting for setting in _declared_configs()}
        setting = settings[sensitivity["name"]]
        if model["config"] != asdict(setting.build()):
            return False
        rebuilt, changes = prepare_sensitivity_study(
            dataset, [Attestation(**att) for att in observations], setting.name,
        )
        if sensitivity != json.loads(json.dumps(changes)):
            return False
        observations = [asdict(att) for att in rebuilt.attestations]
    else:
        rebuilt = prepare_study(
            dataset, [Attestation(**att) for att in observations],
            drop_singletons=model["minimum_present"] == 2,
            include_unobserved_languages=partition is not None,
        )
    if partition is not None and partition["character_ids"] != [
        list(identifier) for identifier in rebuilt.character_ids
    ]:
        return False
    expected_model = json.loads(json.dumps(asdict(build_model(rebuilt, ModelConfig(**model["config"])))))
    # Older bundles omit subsequently added optional fields. Compare all stored
    # model inputs and let dataclass defaults apply to omitted optional inputs.
    return (
        not bundle["is_synthetic"] and bundle["dataset_version"] == dataset.version
        and bundle["context"]["frozen_study"] == manifest
        and bundle["tip_labels"] == rebuilt.tip_labels
        and model == {key: expected_model[key] for key in model}
        and bundle["context"]["attestations"] == json.loads(json.dumps(observations))
    )


def _read_json(path: Path, max_bytes: int) -> dict:
    if path.stat().st_size > max_bytes:
        raise ValueError(f"artifact exceeds size limit: {path.name}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path.name}")
    return value


def verify_evidence(
    evidence_dir: str, expected_commit: str | None = None, probe_draws: int = 16,
    *, observed_source: Path | None = None, catalog: Path | None = None,
    edh_source: Path | None = None, edh_geography: Path | None = None,
) -> dict:
    from ..io.types import DataError

    root = Path(evidence_dir)
    checks = []

    def record(name: str, passed: bool, detail=None) -> None:
        checks.append({"check": name, "ok": bool(passed), "detail": detail})

    if bool(observed_source) != bool(catalog) or bool(edh_source) != bool(edh_geography):
        return {"ok": False, "reason": "source verification requires both files of each pair", "checks": []}
    if edh_source and not observed_source:
        return {"ok": False, "reason": "EDH source verification requires the lexical inputs", "checks": []}
    if probe_draws < 1:
        return {"ok": False, "reason": "probe budget must be positive", "checks": []}
    bundle_dirs = sorted(path.parent for path in root.rglob("bundle.json"))
    if not bundle_dirs:
        return {"ok": False, "reason": "no saved posterior bundles found", "checks": []}

    for directory in bundle_dirs:
        name = str(directory.relative_to(root)) if directory != root else directory.name
        try:
            bundle = _read_json(directory / "bundle.json", 128 * 1024 * 1024)
            manifest = _read_json(directory / "MANIFEST.json", 1024 * 1024)
            context = bundle.get("context", {})
            if context.get("schema_version") != 2:
                record(f"{name}:schema", False, "exact model inputs and draw state required")
                continue
            record(f"{name}:draw_digest", _sha256(bundle["chains"]) == manifest["draw_digest"])
            record(f"{name}:context_digest", _sha256(context) == manifest["context_sha256"])
            bound_source = bundle["source_commit"]
            record(
                f"{name}:source_commit",
                bound_source == manifest["source_commit"]
                and (expected_commit is None or bound_source == expected_commit),
                {"stored": bound_source, "expected": expected_commit},
            )
            well_formed = _draws_well_formed(bundle)
            record(f"{name}:draw_structure", well_formed)
            record(f"{name}:attestation_intervals", _calibrations_match(bundle))
            record(f"{name}:geographic_honesty", _validate_geography(bundle))
            if observed_source is not None:
                record(
                    f"{name}:frozen_observed_inputs",
                    _observed_context_matches(bundle, observed_source, catalog, edh_source, edh_geography),
                )
            if well_formed:
                replay, count = _replay(bundle, probe_draws)
                record(f"{name}:posterior_replay", replay, {"probes": count})
        except (DataError, OSError, KeyError, IndexError, TypeError, ValueError, RecursionError) as exc:
            record(f"{name}:artifact_error", False, str(exc))

    if (root / "heldout.json").is_file():
        try:
            from .prediction import verify_prediction
            if observed_source is None:
                record("prediction:observed_inputs", False, "provide pinned lexical inputs for score replay")
            else:
                from ..io.freeze import load_study

                name = _read_json(root / "heldout.json", 1024 * 1024)["study"]
                dataset, _, _ = load_study(name, observed_source, catalog)
                checks.extend(verify_prediction(root, dataset))
        except (DataError, OSError, KeyError, IndexError, TypeError, ValueError) as exc:
            record("prediction:artifact_error", False, str(exc))

    if (root / "sensitivity.json").is_file():
        try:
            from .sensitivity import verify_sensitivity

            checks.extend(verify_sensitivity(root))
        except (OSError, KeyError, IndexError, TypeError, ValueError) as exc:
            record("sensitivity:artifact_error", False, str(exc))

    if (root / "RUN_PLAN.json").is_file():
        try:
            from .campaign import verify_inference_campaign

            plan = _read_json(root / "RUN_PLAN.json", 1024 * 1024)
            if plan.get("component") == "observed_inference":
                record("inference:pinned_inputs", observed_source is not None)
                checks.extend(verify_inference_campaign(root, expected_commit))
        except (OSError, KeyError, IndexError, TypeError, ValueError, subprocess.SubprocessError) as exc:
            record("inference:artifact_error", False, str(exc))

    acceptance_path = root / "ACCEPTANCE.json"
    if acceptance_path.is_file():
        try:
            acceptance = _read_json(acceptance_path, 1024 * 1024)
            if acceptance["profile"] == "full":
                record("full_acceptance", False, "the complete campaign verifier is unavailable")
        except (OSError, KeyError, TypeError, ValueError) as exc:
            record("acceptance_manifest", False, str(exc))
    return {"ok": all(check["ok"] for check in checks), "n_bundles": len(bundle_dirs), "checks": checks}
