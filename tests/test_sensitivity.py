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


def test_observed_shards_require_complete_identical_protocol_before_merge(tmp_path, monkeypatch):
    from diachron.sensitivity.run import _declared_configs, merge_sensitivity_shards
    import diachron.sensitivity.run as runner
    import pytest

    dataset, _ = simulate_dataset(n_languages=6, n_concepts=12, seed=490, name="shard-fixture")
    dataset = replace(dataset, is_synthetic=False)
    monkeypatch.setattr(runner, 'working_tree_clean', lambda: True)
    shards = []
    for setting in _declared_configs():
        out = tmp_path / setting.name
        result = run_sensitivity(str(out), 2, 20, dataset=dataset, n_chains=1,
                                 only_configuration=setting.name)
        assert result['n_configs'] == 1
        shards.append(out)
    with pytest.raises(ValueError, match='every predeclared'):
        merge_sensitivity_shards(shards[:-1], tmp_path / 'partial')
    merged = tmp_path / 'complete'
    result = merge_sensitivity_shards(shards, merged)
    assert result['n_configs'] == 20
    assert all(check['ok'] for check in verify_sensitivity(merged))
    assert (shards[0] / 'pop_size_1000/dated/bundle.json').stat().st_ino == (
        merged / 'pop_size_1000/dated/bundle.json'
    ).stat().st_ino
    with pytest.raises(ValueError, match='immutable'):
        run_sensitivity(str(merged), 2, 20, dataset=dataset, n_chains=1, resume=True)
    path = shards[-1] / 'SENSITIVITY_PLAN.json'
    changed = json.loads(path.read_text())
    changed['iterations'] = 30
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match='different source, input, or budgets'):
        merge_sensitivity_shards(shards, tmp_path / 'different')
