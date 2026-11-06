"""Matched recovery generates proper serial trees and persists complete units."""

from pathlib import Path
import json

import numpy as np
import pytest

from diachron.sensitivity.matched import generate_matched_case, run_matched_recovery, serial_coalescent
from diachron.verify.recovery import _validate_time_block_row


def test_serial_coalescent_two_tip_wait_is_exponential_after_oldest_tip():
    rng = np.random.default_rng(1620)
    roots = [serial_coalescent(['a', 'b'], [0., 2500.], 1200., rng).root_age() for _ in range(15000)]
    waits = np.asarray(roots) - 2500
    assert waits.min() > 0
    assert abs(waits.mean() - 1200) < 30
    assert abs(np.mean(waits < 1200) - (1 - np.exp(-1))) < .015


def test_matched_case_masks_and_dated_tips_follow_the_target():
    study, config, tree, regime, truth = generate_matched_case(13)
    assert config.use_borrowing and len(truth['borrow_rates']) == 2
    assert len(study.tip_dates) == regime['n_dated_tips'] == 2
    assert len(study.columns) == regime['n_concepts']
    assert all(column.count(1) >= 2 and 0 in column for column in study.columns)
    assert tree.times_valid() and tree.root_age() > max(tree.heights[:len(study.tip_labels)])


def test_matched_recovery_resumes_without_regenerating_saved_cases(tmp_path):
    out = Path(tmp_path) / 'matched'
    result = run_matched_recovery(2, out, 1, warmup=5, iterations=20)
    assert result['completed'] == 2 and result['failed'] == 0
    saved = (out / 'units/case-000.json').stat().st_mtime_ns
    assert result == run_matched_recovery(2, out, 1, True, warmup=5, iterations=20)
    assert (out / 'units/case-000.json').stat().st_mtime_ns == saved


def test_joint_recovery_row_replays_clade_bits_and_rejects_altered_counts(tmp_path):
    out = tmp_path / 'matched-time'
    result = run_matched_recovery(1, out, 1, protocol=4, warmup=5, iterations=100)
    assert result['completed'] == 1
    row = json.loads((out / 'units/case-000.json').read_text())
    _validate_time_block_row(row, result, row['root_samples'], row['true_clade_masses'])
    row['true_clade_membership'][0] ^= 1
    with pytest.raises(ValueError, match='frequencies differ'):
        _validate_time_block_row(row, result, row['root_samples'], row['true_clade_masses'])
