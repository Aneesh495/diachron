"""Deterministic evidence interventions for observed sensitivity studies."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, replace

from ..infer.campaign import prepare_study
from ..models.tip_dates import TipDate


def merge_proposed_cognates(dataset):
    """An extreme alternative, merging source-linked classes within a concept.

    The source scale is ordinal. This intervention assigns no probability to it.
    Links with missing endpoints or ambiguous concept identity are left unused.
    """
    concepts = defaultdict(set)
    for judgment in dataset.judgments:
        concepts[judgment.cognate_class].add(judgment.concept_id)
    parent = {identifier: identifier for identifier in concepts}

    def root(identifier):
        while parent[identifier] != identifier:
            parent[identifier] = parent[parent[identifier]]
            identifier = parent[identifier]
        return identifier

    used, excluded = [], []
    for record in sorted(dataset.cognate_classes, key=lambda row: row.id):
        target = record.proposed_cognate_id
        if target is None:
            continue
        if record.id not in concepts or target not in concepts:
            excluded.append({"source": record.id, "target": target, "reason": "unused endpoint"})
            continue
        if len(concepts[record.id]) != 1 or concepts[record.id] != concepts[target]:
            excluded.append({"source": record.id, "target": target, "reason": "concept disagreement"})
            continue
        a, b = root(record.id), root(target)
        parent[max(a, b)] = min(a, b)
        used.append({"source": record.id, "target": target, "source_ordinal_scale": record.proposed_scale})
    mapping = {identifier: root(identifier) for identifier in parent if root(identifier) != identifier}
    altered = replace(dataset, judgments=[
        replace(judgment, cognate_class=mapping.get(judgment.cognate_class, judgment.cognate_class))
        for judgment in dataset.judgments
    ])
    return altered, {"mapping": mapping, "used_links": used, "excluded_links": excluded}


def catalog_clades(dataset, labels):
    """Condition on first-level IE-CoR source groups, explicitly as assumptions."""
    positions = {label: index for index, label in enumerate(labels)}
    groups = defaultdict(list)
    for language in dataset.languages:
        if language.source_classification and language.id in positions:
            groups[language.source_classification.split(";")[0]].append(positions[language.id])
    return {name: sorted(taxa) for name, taxa in sorted(groups.items()) if 1 < len(taxa) < len(labels)}


def prepare_sensitivity_study(dataset, attestations, name):
    """Apply one named intervention with a complete, reproducible change record."""
    report = {
        "name": name, "drop_singletons": name != "include_singletons",
        "include_loans_as_present": name != "remove_loans",
        "tip_date_overrides": {}, "catalog_groups": {}, "cognacy_alternative": None,
    }
    if name == "alternative_cognacy_merge":
        dataset, report["cognacy_alternative"] = merge_proposed_cognates(dataset)
    observations = list(attestations)
    if name == "half_calibrations" and dataset.is_synthetic:
        observations = sorted(observations, key=lambda record: record.id)[::2]
    study = prepare_study(
        dataset, observations, drop_singletons=report["drop_singletons"],
        include_loans_as_present=report["include_loans_as_present"],
    )
    if not dataset.is_synthetic and name in {"half_calibrations", "broad_all_tip_dates"}:
        dates = sorted(study.tip_dates)
        selected = dates[1::2] if name == "half_calibrations" else dates
        for taxon in selected:
            # This common proper prior does not reuse the withheld source mean.
            date = TipDate("normal", 4000.0, 4000.0)
            study.tip_dates[taxon] = date
            report["tip_date_overrides"][study.tip_labels[taxon]] = asdict(date)
    if name == "catalog_monophyly":
        report["catalog_groups"] = catalog_clades(dataset, study.tip_labels)
        study.required_clades = list(report["catalog_groups"].values())
    report["dated_tips"] = len(study.tip_dates)
    report["attestations"] = len(study.attestations)
    report["character_ids"] = study.character_ids
    return study, report
