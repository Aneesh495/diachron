"""Sensitivity changes source evidence explicitly and saves auditable chains."""

from dataclasses import replace
import json

from diachron.io.simulate import simulate_dataset
from diachron.io.types import CognateClassRecord
from diachron.sensitivity.interventions import merge_proposed_cognates, prepare_sensitivity_study
from diachron.sensitivity.run import run_sensitivity
from diachron.verify.sensitivity import verify_sensitivity


def test_observed_date_subset_and_catalog_constraints_are_effective():
    dataset, _ = simulate_dataset(n_languages=8, n_concepts=12, seed=193, name="interventions-test")
    languages = [replace(
        language, is_historical=True, date_distribution="normal", date_mean=1000 + 100 * index,
        date_sd=100.0, date_offset=0.0, source_classification="A;sub" if index < 4 else "B;sub",
    ) for index, language in enumerate(dataset.languages)]
    observed = replace(dataset, languages=languages, is_synthetic=False)
    baseline, _ = prepare_sensitivity_study(observed, [], "pop_size_3000")
    subset, report = prepare_sensitivity_study(observed, [], "half_calibrations")
    assert len(report["tip_date_overrides"]) == 4
    assert subset.tip_dates[0] == baseline.tip_dates[0]
    assert subset.tip_dates[1].mean == subset.tip_dates[1].sd == 4000
    assert baseline.tip_dates[1].mean == 1100
    broad, _ = prepare_sensitivity_study(observed, [], "broad_all_tip_dates")
    assert all(date.mean == date.sd == 4000 for date in broad.tip_dates.values())
    constrained, report = prepare_sensitivity_study(observed, [], "catalog_monophyly")
    assert report["catalog_groups"] == {"A": [0, 1, 2, 3], "B": [4, 5, 6, 7]}
    assert constrained.required_clades == list(report["catalog_groups"].values())


def test_alternative_merges_preserve_concept_boundaries_and_original_records():
    dataset, _ = simulate_dataset(n_languages=8, n_concepts=12, seed=193, name="alternative-test")
    first = dataset.judgments[0]
    same = next(row for row in dataset.judgments if row.concept_id == first.concept_id
                and row.cognate_class != first.cognate_class)
    other = next(row for row in dataset.judgments if row.concept_id != first.concept_id)
    dataset.cognate_classes = [
        CognateClassRecord(first.cognate_class, None, None, 'source', same.cognate_class, 2, None, None),
        CognateClassRecord(other.cognate_class, None, None, 'source', first.cognate_class, 1, None, None),
    ]
    altered, report = merge_proposed_cognates(dataset)
    assert len(report["used_links"]) == 1
    assert report["excluded_links"][0]["reason"] == "concept disagreement"
    assert altered.cognate_classes == dataset.cognate_classes
    assert len(altered.judgments) == len(dataset.judgments)
    assert altered.judgments != dataset.judgments


def test_saved_sensitivity_is_recomputable_and_resumable(tmp_path):
    dataset, _ = simulate_dataset(n_languages=6, n_concepts=15, seed=499, name="sensitivity-test")
    out = tmp_path / "sensitivity"
    result = run_sensitivity(str(out), 5, 20, dataset=dataset, n_chains=1)
    assert result["n_configs"] == 15
    assert result == run_sensitivity(str(out), 5, 20, dataset=dataset, n_chains=1, resume=True)
    checks = verify_sensitivity(out)
    assert all(check["ok"] for check in checks), checks
    report = json.loads((out / "sensitivity.json").read_text())
    report["configs"][0]["root_age_mean"] += 10
    (out / "sensitivity.json").write_text(json.dumps(report))
    assert any(not check["ok"] for check in verify_sensitivity(out))
