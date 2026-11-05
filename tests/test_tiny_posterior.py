"""Independent composed-model targets include dates, geography and contact."""

from diachron.verify.tiny_posterior import run_tiny_posterior_audit


def test_composed_posterior_and_production_acceptance_match_exact_grid():
    result = run_tiny_posterior_audit()
    assert result["passed"], result
    assert all(model["grid_states"] == 72 for model in result["models"].values())
