import math

import numpy as np
import pytest
from scipy.stats import norm

from diachron import _core
from diachron.infer.calendar import propose_calendar_scale, slice_offset
from diachron.models.posterior import ModelConfig, PosteriorModel
from diachron.models.tip_dates import TipDate
from diachron.models.tree import PhyloTree


def test_slice_recovers_normal_mass_with_finite_bracket():
    rng = _core.PyRng(72090)
    x = 0.0
    draws = []
    for index in range(21000):
        x += slice_offset(lambda delta: -0.5 * ((x + delta - .4) / .19) ** 2,
                          -0.5 * ((x - .4) / .19) ** 2, rng)
        if index >= 1000:
            draws.append(x)
    np.testing.assert_allclose(np.quantile(draws, [.1, .5, .9]),
                               norm.ppf([.1, .5, .9], .4, .19), atol=.012)


@pytest.mark.parametrize("config", [ModelConfig.dated(), ModelConfig.geography(), ModelConfig.borrowing()])
def test_calendar_scale_preserves_exposure_and_corrects_full_score(config):
    model = PosteriorModel(config, [[1, 1, 0], [0, 0, 1], [1, -1, 0]],
                           tip_dates={1: TipDate("normal", 400, 100)},
                           tip_coords={0: (52., 10.), 1: (50., 13.), 2: (44., 25.)},
                           borrow_candidates=[(0, 1)] if config.use_borrowing else [])
    tree = PhyloTree([3, 3, 4, 4, -1], [0., 400., 0., 1100., 2500.], [0, 1, 2, -1, -1])
    gain, loss, clock, diffusion = .0003, .001, .9, .012
    rates = [.08] if config.use_borrowing else []
    proposal = propose_calendar_scale(model, tree, clock, diffusion, _core.PyRng(675))
    assert proposal is not None
    candidate, new_clock, new_diffusion, correction = proposal
    factor = candidate.root_age() / tree.root_age()
    assert not math.isclose(factor, 1.0)
    assert candidate.heights[1] == pytest.approx(400 * factor)
    assert tree.heights[1] == 400
    assert candidate.parent == tree.parent
    np.testing.assert_allclose(np.array(candidate.branch_lengths()) * new_clock,
                               np.array(tree.branch_lengths()) * clock)
    if config.use_geography:
        assert new_diffusion * factor == pytest.approx(diffusion)
    else:
        assert new_diffusion == diffusion
    assert model.log_posterior(candidate, gain, loss, new_clock, new_diffusion, rates) + correction == pytest.approx(
        model.log_posterior(tree, gain, loss, clock, diffusion, rates), abs=1e-9)
