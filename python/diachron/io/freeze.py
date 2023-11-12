"""Freeze verified CLDF inputs into versioned Arrow tables."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import tempfile
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, replace
from pathlib import Path

from .cldf import load_cldf_wordlist
from .types import DataError, Dataset

SCHEMA_VERSION = 2


def verify_inputs(root: Path, files: dict) -> None:
    """Check the actual bytes before parsing any source tables."""
    root = root.resolve()
    for name, expected in files.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise DataError(f"missing pinned input: {name}")
        if path.stat().st_size != expected["bytes"]:
            raise DataError(f"byte count differs for pinned input: {name}")
        with path.open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != expected["sha256"]:
                raise DataError(f"digest differs for pinned input: {name}")


def load_study(study: str, source: Path, catalog: Path) -> tuple[Dataset, list[dict], dict]:
    """Return an observed lexical study and a crosswalk preserving both taxonomies."""
    if study not in ("iecor", "abvd"):
        raise DataError(f"unknown observed study: {study}")
    manifest = json.loads(Path(f"data/manifests/{study}.json").read_text(encoding="utf-8"))
    glottolog = json.loads(Path("data/manifests/glottolog.json").read_text(encoding="utf-8"))
    verify_inputs(source, manifest["files"])
    verify_inputs(catalog.parent, {catalog.name: glottolog["archive"]})
    dataset = load_cldf_wordlist(str(source / manifest["metadata"]), profile=study)
    if dataset.license != manifest["license"]:
        raise DataError("CLDF metadata disagrees with the pinned source license")
    dataset = replace(
        dataset, name=study, version=manifest["version"], doi=manifest.get("doi"),
        source_url=manifest["url"],
    )
    raw_totals = {
        "languages": len(dataset.languages),
        "concepts": len(dataset.concepts),
        "lexemes": len(dataset.lexemes),
        "judgments": len(dataset.judgments),
    }
    dataset, exclusions = select_languages(dataset, manifest["selection"])

    with zipfile.ZipFile(catalog) as archive:
        with archive.open("languoid.csv") as raw:
            catalog_by_id = {
                row["id"]: row for row in csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"))
            }
    crosswalk = []
    for language in dataset.languages:
        entry = catalog_by_id.get(language.glottocode, {})
        crosswalk.append({
            "language_id": language.id,
            "glottocode": language.glottocode,
            "catalog_version": glottolog["version"],
            "catalog_name": entry.get("name"),
            "catalog_level": entry.get("level"),
            "catalog_parent": entry.get("parent_id") or None,
            "catalog_family": entry.get("family_id") or None,
            "source_classification": language.source_classification,
        })
    observed_cells = dataset.observed_cell_count()
    lexemes = {form.id: form for form in dataset.lexemes}
    coded_cells = {
        (lexemes[j.lexeme_id].language_id, j.concept_id) for j in dataset.judgments
    }
    n_languages = dataset.n_eligible_languages
    n_concepts = dataset.n_eligible_concepts
    stats = {
        "schema_version": SCHEMA_VERSION,
        "study": study,
        "is_synthetic": False,
        "source": manifest,
        "catalog": glottolog,
        "source_totals": raw_totals,
        "eligible_languages": n_languages,
        "eligible_concepts": n_concepts,
        "observed_cells": observed_cells,
        "lexically_missing_cells": n_languages * n_concepts - observed_cells,
        "cells_with_cognate_judgments": len(coded_cells),
        "observed_cells_without_judgments": observed_cells - len(coded_cells),
        "lexemes": len(dataset.lexemes),
        "judgments": len(dataset.judgments),
        "loan_judgments": sum(j.is_loan for j in dataset.judgments),
        "doubtful_judgments": sum(j.doubtful is True for j in dataset.judgments),
        "partial_judgments": sum(j.is_partial for j in dataset.judgments),
        "forms_with_multiple_judgments": sum(
            count > 1 for count in Counter(j.lexeme_id for j in dataset.judgments).values()
        ),
        "historical_varieties": sum(lang.is_historical for lang in dataset.languages),
        "resolved_glottocodes": sum(r["catalog_name"] is not None for r in crosswalk),
        "language_exclusions": dict(exclusions),
        "language_ids": [language.id for language in dataset.languages],
    }
    if n_languages < manifest["selection"]["required_varieties"] or n_concepts < 100:
        raise DataError(f"observed study fails its frozen coverage floor: {study}")
    return dataset, crosswalk, stats


def select_languages(dataset: Dataset, rules: dict) -> tuple[Dataset, Counter]:
    """Select by judged coverage and source identity before character construction."""
    lexemes = {form.id: form for form in dataset.lexemes}
    judged_concepts: dict[str, set[str]] = defaultdict(set)
    for j in dataset.judgments:
        judged_concepts[lexemes[j.lexeme_id].language_id].add(j.concept_id)

    exclusions: Counter = Counter()
    candidates = []
    for language in dataset.languages:
        if language.family != rules["family"]:
            exclusions["source_family"] += 1
        elif len(judged_concepts[language.id]) < rules["min_judged_concepts"]:
            exclusions["judged_coverage"] += 1
        elif rules["distinct_glottocodes"] and not language.glottocode:
            exclusions["missing_glottocode"] += 1
        else:
            candidates.append(language)
    if rules["distinct_glottocodes"]:
        distinct = {}
        for language in sorted(candidates, key=lambda l: (-len(judged_concepts[l.id]), l.id)):
            distinct.setdefault(language.glottocode, language)
        exclusions["duplicate_glottocode"] = len(candidates) - len(distinct)
        candidates = list(distinct.values())
    if rules["max_varieties"] is not None:
        seed = rules["selection_seed"]
        candidates.sort(key=lambda l: hashlib.sha256(f"{seed}:{l.id}".encode()).digest())
        exclusions["bounded_sample"] = max(0, len(candidates) - rules["max_varieties"])
        candidates = candidates[:rules["max_varieties"]]
    candidates.sort(key=lambda l: l.id)
    selected_ids = {language.id for language in candidates}
    forms = [form for form in dataset.lexemes if form.language_id in selected_ids]
    form_ids = {form.id for form in forms}
    concept_ids = {form.concept_id for form in forms}
    judgments = [j for j in dataset.judgments if j.lexeme_id in form_ids]
    used_classes = {j.cognate_class for j in judgments}
    selected = replace(
        dataset, languages=candidates, lexemes=forms,
        concepts=[c for c in dataset.concepts if c.id in concept_ids],
        judgments=judgments,
        cognate_classes=[c for c in dataset.cognate_classes if c.id in used_classes],
    )
    selected.validate()
    return selected, exclusions


# Explicit schemas retain types even when a dataset has no annotations of a kind.
def _schema(fields):
    import pyarrow as pa

    kinds = {"s": pa.string(), "b": pa.bool_(), "f": pa.float64(), "i": pa.int32()}
    return pa.schema([(name, kinds[kind]) for name, kind in fields.items()])


def write_study(out: Path, dataset: Dataset, crosswalk: list[dict], report: dict) -> None:
    import pyarrow as pa

    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise DataError(f"frozen output already exists: {out}")
    identity = {
        b"diachron.schema": str(SCHEMA_VERSION).encode(),
        b"diachron.source": dataset.source_url.encode(),
        b"diachron.license": dataset.license.encode(),
        b"diachron.version": dataset.version.encode(),
    }
    tables = {
        "languages": (dataset.languages, {
            "id": "s", "name": "s", "glottocode": "s", "family": "s", "latitude": "f",
            "longitude": "f", "is_extinct": "b", "tip_age": "f", "is_historical": "b",
            "source_classification": "s", "date_distribution": "s", "date_mean": "f",
            "date_sd": "f", "date_offset": "f",
        }),
        "concepts": (dataset.concepts, {
            "id": "s", "gloss": "s", "concepticon_id": "s", "semantic_field": "s",
        }),
        "lexemes": (dataset.lexemes, {
            "id": "s", "language_id": "s", "concept_id": "s", "form": "s", "original": "s",
            "transliteration_source": "s", "source_ref": "s", "source_value": "s",
            "loan_annotation": "b",
        }),
        "judgments": (dataset.judgments, {
            "lexeme_id": "s", "concept_id": "s", "cognate_class": "s", "confidence": "f",
            "is_loan": "b", "is_partial": "b", "annotator": "s", "source_ref": "s",
            "id": "s", "doubtful": "b", "segment_slice": "s", "detection_method": "s",
            "loan_event": "s", "loan_source": "s",
        }),
        "cognate_classes": (dataset.cognate_classes, {
            "id": "s", "root_form": "s", "root_language": "s", "source_ref": "s",
            "proposed_cognate_id": "s", "proposed_scale": "i", "superset_id": "s",
            "parallel_derivation": "b",
        }),
        "taxonomy": (crosswalk, {
            "language_id": "s", "glottocode": "s", "catalog_version": "s", "catalog_name": "s",
            "catalog_level": "s", "catalog_parent": "s", "catalog_family": "s",
            "source_classification": "s",
        }),
    }
    with tempfile.TemporaryDirectory(dir=out.parent, prefix=".freeze-") as tmp:
        stage = Path(tmp) / "study"
        stage.mkdir()
        digests = {}
        for name, (rows, fields) in tables.items():
            schema = _schema(fields).with_metadata(identity)
            table = pa.Table.from_pylist(
                [row if isinstance(row, dict) else asdict(row) for row in rows], schema=schema,
            )
            path = stage / f"{name}.arrow"
            with pa.ipc.new_file(str(path), schema) as writer:
                writer.write_table(table)
            with path.open("rb") as stream:
                digests[path.name] = hashlib.file_digest(stream, "sha256").hexdigest()
        (stage / "DATASETS.json").write_text(
            json.dumps({**report, "arrow_sha256": digests}, indent=2) + "\n", encoding="utf-8",
        )
        stage.rename(out)
