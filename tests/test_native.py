"""Tests for the native Rust bindings: gamma rates and the native MCMC chain."""

from __future__ import annotations

from diachron import _core


def test_discrete_gamma_mean_one():
    for alpha in (0.1, 0.5, 1.0, 5.0):
        rates = _core.discrete_gamma_rates(alpha, 4)
        assert abs(sum(rates) / len(rates) - 1.0) < 1e-9
        assert all(r > 0 for r in rates)


def test_gamma_heterogeneity_spread():
    # Smaller alpha gives a wider spread of rate categories.
    narrow = _core.discrete_gamma_rates(100.0, 4)
    wide = _core.discrete_gamma_rates(0.2, 4)
    spread_narrow = max(narrow) - min(narrow)
    spread_wide = max(wide) - min(wide)
    assert spread_wide > spread_narrow


def _tiny_tree_inputs():
    parent = [4, 4, 5, 5, 6, 6, -1]
    height = [0.0, 0.0, 0.0, 0.0, 400.0, 600.0, 1000.0]
    taxon = [0, 1, 2, 3, -1, -1, -1]
    cols = [[1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1]]
    return parent, height, taxon, cols


def test_native_chain_runs():
    parent, height, taxon, cols = _tiny_tree_inputs()
    draws = _core.run_native_chain(
        parent, height, taxon, 6, 5e-4, 1e-3, cols, 3000.0, 42, 300, 900, 10
    )
    assert len(draws) == 90
    # Each draw is (iteration, log_post, root_age, gain, loss).
    assert all(len(d) == 5 for d in draws)
    assert all(d[1] == d[1] for d in draws)  # finite log posteriors (not NaN)
    assert all(d[2] > 0 for d in draws)  # positive root ages


def test_native_chain_reproducible():
    parent, height, taxon, cols = _tiny_tree_inputs()
    a = _core.run_native_chain(parent, height, taxon, 6, 5e-4, 1e-3, cols, 3000.0, 7, 200, 400, 5)
    b = _core.run_native_chain(parent, height, taxon, 6, 5e-4, 1e-3, cols, 3000.0, 7, 200, 400, 5)
    assert a == b


def test_native_chain_different_seeds_differ():
    parent, height, taxon, cols = _tiny_tree_inputs()
    a = _core.run_native_chain(parent, height, taxon, 6, 5e-4, 1e-3, cols, 3000.0, 1, 200, 400, 5)
    b = _core.run_native_chain(parent, height, taxon, 6, 5e-4, 1e-3, cols, 3000.0, 2, 200, 400, 5)
    assert a != b
