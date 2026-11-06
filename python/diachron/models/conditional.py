"""Sample the same retained binary patterns used by the likelihood."""

import numpy as np

from diachron import _core


def sample_selected_patterns(tree, gain, loss, masks, minimum_present, rng, *,
                             rate_categories=(1.0,), contacts=(), clock_multipliers=()):
    parent, lengths, heights, taxa, root = tree.to_core_arrays()
    if clock_multipliers:
        lengths = [length * multiplier for length, multiplier in zip(lengths, clock_multipliers)]
    return _core.binary_conditioned_draws(
        parent, lengths, heights, taxa, root, gain, loss, np.asarray(masks, dtype=bool).tolist(),
        list(rate_categories), minimum_present, int(rng.integers(2 ** 63)), list(contacts),
    )
