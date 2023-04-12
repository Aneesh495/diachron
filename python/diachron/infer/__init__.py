"""MCMC inference: chain driver, proposals, checkpointing, and resumption.

The sampler is a Metropolis-within-Gibbs cycle. Each iteration applies a randomly
chosen move (scale a scalar, slide a node time, or rewire topology), evaluates
the composed log posterior, and accepts or rejects. Adaptation of proposal step
sizes happens only during a warmup phase and is frozen before any retained
sample is drawn, so the retained chain targets the stationary posterior without
the moving-target bias that online adaptation would introduce.

State is checkpointed atomically (RNG words, tree, parameters, iteration count,
and frozen step sizes) so an interrupted chain resumes bit-for-bit on CPU.
"""

from .chain import Chain, ChainConfig, ChainState
from .moves import MoveSet

__all__ = ["Chain", "ChainConfig", "ChainState", "MoveSet"]
