import numpy as np
import pytest

from diachron.infer.time_block import TimeChart, age_prior_gradient
from diachron.models.priors import coalescent_log_prior
from diachron.models.tip_dates import TipDate
from diachron.models.tree import PhyloTree


def test_time_chart_jacobian_and_serial_age_force_match_finite_differences():
    tree = PhyloTree([3, 3, 4, 4, -1], [100., 270., 0., 1000., 2400.], [0, 1, 2, -1, -1])
    dates = {0: TipDate("normal", 120, 40), 1: TipDate("offset_lognormal", 200, .4, 50)}
    chart = TimeChart(tree, dates)
    point, anchors = chart.encode(.9)
    restored, clock, _ = chart.decode(point)
    np.testing.assert_allclose(restored.heights, tree.heights)
    assert clock == pytest.approx(.9)

    def evaluate(coordinates):
        candidate, clock, _ = chart.decode(coordinates)
        log_prior = coalescent_log_prior(candidate.heights, candidate.parent, 3000.)
        log_prior += sum(dates[t].log_density(candidate.heights[i])
                         for i, t in enumerate(candidate.taxon) if t in dates)
        return np.array([*[candidate.heights[i] for i in chart.nodes], clock]), log_prior + coordinates.sum()

    differences = np.eye(len(point)) * 1e-6
    columns = []
    reference = []
    for delta in differences:
        before, lo = evaluate(point - delta)
        after, hi = evaluate(point + delta)
        columns.append((after - before) / 2e-6)
        reference.append((hi - lo) / 2e-6)
    sign, log_jacobian = np.linalg.slogdet(np.column_stack(columns))
    assert sign == 1
    assert log_jacobian == pytest.approx(point.sum(), abs=1e-8)
    force = chart.pullback(age_prior_gradient(tree, 3000., dates), point, anchors)
    np.testing.assert_allclose(force, reference, atol=1e-8)
