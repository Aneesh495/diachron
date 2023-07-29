"""Check saved posterior structure and replay a bounded set of exact draws."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from ..models.posterior import Calibration, ModelConfig, PosteriorModel
from ..models.tree import PhyloTree


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
                draw["geo_sigma2"], draw["borrow_rates"],
            )
            if not math.isfinite(score) or abs(score - draw["log_posterior"]) > 1e-9:
                return False, count
    return count > 0, count


def verify_evidence(evidence_dir: str, expected_commit: str | None = None, probe_draws: int = 16) -> dict:
    root = Path(evidence_dir)
    checks = []

    def record(name: str, passed: bool, detail=None) -> None:
        checks.append({"check": name, "ok": bool(passed), "detail": detail})

    if probe_draws < 1:
        return {"ok": False, "reason": "probe budget must be positive", "checks": []}
    bundle_dirs = sorted(path.parent for path in root.rglob("bundle.json"))
    if not bundle_dirs:
        return {"ok": False, "reason": "no saved posterior bundles found", "checks": []}

    for directory in bundle_dirs:
        name = str(directory.relative_to(root)) if directory != root else directory.name
        try:
            bundle = json.loads((directory / "bundle.json").read_text(encoding="utf-8"))
            manifest = json.loads((directory / "MANIFEST.json").read_text(encoding="utf-8"))
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
            if well_formed:
                replay, count = _replay(bundle, probe_draws)
                record(f"{name}:posterior_replay", replay, {"probes": count})
        except (OSError, KeyError, IndexError, TypeError, ValueError, RecursionError) as exc:
            record(f"{name}:artifact_error", False, str(exc))

    acceptance_path = root / "ACCEPTANCE.json"
    if acceptance_path.is_file():
        try:
            acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
            if acceptance["profile"] == "full":
                record("full_acceptance", False, "the complete campaign verifier is unavailable")
        except (OSError, KeyError, TypeError, ValueError) as exc:
            record("acceptance_manifest", False, str(exc))
    return {"ok": all(check["ok"] for check in checks), "n_bundles": len(bundle_dirs), "checks": checks}
