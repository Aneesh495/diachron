"""Evidence examples and catalog comparisons for saved observed draws."""

from __future__ import annotations

from collections import defaultdict
import csv
import io
import json
from pathlib import Path
import zipfile

from ..infer.campaign import prepare_study
from ..io.freeze import load_study
from ..models.constraints import clade_masks
from ..models.tree import PhyloTree


def _taxonomy_groups(dataset, labels, catalog, version):
    with zipfile.ZipFile(catalog) as archive, archive.open("languoid.csv") as raw:
        entries = {row["id"]: row for row in csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"))}
    groups = defaultdict(set)
    display = {}
    positions = {identifier: index for index, identifier in enumerate(labels)}
    for language in dataset.languages:
        if language.id not in positions:
            continue
        index = positions[language.id]
        chain = []
        if language.source_classification:
            for name in language.source_classification.split(";"):
                chain.append(name)
                key = (dataset.name, ";".join(chain))
                groups[key].add(index)
                display[key] = {"label": " / ".join(chain), "url": dataset.source_url,
                                "taxonomy": f"IE-CoR {dataset.version}"}
        identifier = language.glottocode
        visited = set()
        while identifier in entries and identifier not in visited and len(visited) < 100:
            visited.add(identifier)
            entry = entries[identifier]
            key = ("glottolog", identifier)
            groups[key].add(index)
            display[key] = {"label": entry["name"], "url": f"https://glottolog.org/resource/languoid/id/{identifier}",
                            "taxonomy": f"Glottolog {version}"}
            identifier = entry["parent_id"]
    selected = sorted(
        (key for key, taxa in groups.items() if 1 < len(taxa) < len(labels)),
        key=lambda key: (-len(groups[key]), key),
    )[:96]
    return [{"id": ":".join(key), **display[key], "taxa": sorted(groups[key]), "models": {}}
            for key in selected]


def record_taxonomy_frequencies(groups, name, bundle):
    if not groups:
        return
    masks = [sum(1 << taxon for taxon in group["taxa"]) for group in groups]
    counts = dict.fromkeys(masks, 0)
    total = 0
    for chain in bundle.chains:
        for draw in chain:
            tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"])
            for mask in clade_masks(tree).intersection(counts):
                counts[mask] += 1
            total += 1
    for group, mask in zip(groups, masks):
        group["models"][name] = {"clade_draws": counts[mask], "total_draws": total}


def _matrix_examples(study, dataset):
    by_id = {key: index for index, key in enumerate(study.character_ids)}
    concepts = {concept.id: concept.gloss for concept in dataset.concepts}
    names = {language.id: language.name for language in dataset.languages}
    steps = []

    def cell(column, taxon):
        return {"column": column, "taxon": taxon}

    for column, states in enumerate(study.columns):
        if -1 in states and 0 in states:
            steps.append({
                "id": "missing", "title": "Missing evidence and coded absence",
                "body": "Compare an unknown state with a coded absence in the same class. An absence requires classified lexical evidence for that meaning.",
                "cells": [cell(column, states.index(-1)), cell(column, states.index(0))],
                "links": [{"label": "Frozen lexical source", "url": dataset.source_url}],
            })
            break
    for record in dataset.cognate_classes:
        if record.proposed_cognate_id is None:
            continue
        pair = next(((index, by_id[(concept, record.proposed_cognate_id)])
                     for (concept, cognate), index in by_id.items()
                     if cognate == record.id and (concept, record.proposed_cognate_id) in by_id
                     and cognate != record.proposed_cognate_id), None)
        if pair is None:
            continue
        first, second = pair
        meaning = concepts[study.character_ids[first][0]]
        ordinal = (f" The source records ordinal scale {record.proposed_scale}."
                   if record.proposed_scale is not None else "")
        steps.append({
            "id": "cognacy", "title": "An alternate cognacy judgment",
            "body": f"For {meaning}, the source proposes a link between classes {record.id} and "
                    f"{record.proposed_cognate_id}. The baseline scores separate classes.{ordinal}",
            "cells": [cell(index, study.columns[index].index(1)) for index in (first, second)],
            "links": [{"label": "Source cognate annotations", "url": dataset.source_url}],
        })
        break
    if dataset.name == "iecor" and all(("61", identifier) in by_id for identifier in ("4364", "227")):
        loan, source = (by_id[("61", identifier)] for identifier in ("4364", "227"))
        loan_members = {i for i, value in enumerate(study.columns[loan]) if value == 1}
        source_members = {i for i, value in enumerate(study.columns[source]) if value == 1}
        crosses = any(a in loan_members and b in source_members or b in loan_members and a in source_members
                      for a, b in study.borrow_candidates)
        candidate = ("The proximity screen includes a pair across the two groups." if crosses else
                     "The proximity screen misses all pairs across these two groups.")
        steps.append({
            "id": "borrowing", "title": "A documented loan annotation",
            "body": f"The source links fruit loan class 4364 to class 227. {candidate} "
                    "The copy likelihood operates within one binary class.",
            "cells": [cell(loan, min(loan_members)), cell(source, min(source_members))],
            "links": [{"label": "IE-CoR loan source records", "url": dataset.source_url}],
        })
    selected = sorted({cell["column"] for step in steps for cell in step["cells"]})
    displayed = selected + [index for index in range(len(study.columns)) if index not in selected]
    displayed = displayed[:400]
    positions = {column: index for index, column in enumerate(displayed)}
    for step in steps:
        for item in step["cells"]:
            item["column"] = positions[item["column"]]
    matrix = {
        "languageLabels": [names[label] for label in study.tip_labels],
        "columns": [study.columns[index] for index in displayed],
        "characterLabels": [f"{concepts[concept]} · class {identifier}"
                            for concept, identifier in (study.character_ids[index] for index in displayed)],
    }
    return matrix, steps


def build_observed_context(bundle, source, catalog):
    dataset, _, report = load_study(bundle.dataset_name, source, catalog)
    study = prepare_study(dataset, drop_singletons=bundle.context["model"]["minimum_present"] == 2)
    if (report != bundle.context["frozen_study"] or study.tip_labels != bundle.tip_labels
            or study.columns != bundle.context["model"]["columns"]):
        raise ValueError("workbench lexical examples require the same frozen character data")
    matrix, tour = _matrix_examples(study, dataset)
    taxonomy = _taxonomy_groups(dataset, study.tip_labels, catalog, report["catalog"]["version"])
    if dataset.name == "iecor":
        case = json.loads(Path("data/cases/knossos.json").read_text())
        body = case["scope"] + " " + case["context"]
        taxon = {}
        if case["language_id"] in study.tip_labels:
            taxon["taxon"] = study.tip_labels.index(case["language_id"])
        else:
            body += " The frozen lexical cohort excludes Mycenaean Greek. These markers provide calendar context."
        tour.append({
            "id": "dating", "title": "Competing dates at Knossos", "body": body, **taxon,
            "dates": [{**item, "age": 2000 + item["calendar_marker_bce"] - 1}
                      for item in case["alternatives"]],
            "links": [{"label": "Later reassessment of the deposits", "url": case["context_url"]}],
        })
    return {"matrix": matrix, "tour": tour, "taxonomy": taxonomy}
