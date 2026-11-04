"""The MCMC chain: warmup, adaptation freeze, sampling, and checkpointing.

A chain holds a `PosteriorModel`, a current `PhyloTree`, the scalar parameters,
and a project RNG. Running proceeds in two phases:

1. **Warmup** adapts the move step sizes toward a target acceptance rate using a
   simple Robbins-Monro update. No samples are retained in warmup.
2. **Sampling** freezes the step sizes and records a draw every `thin`
   iterations. Because adaptation is off, the retained draws target the
   stationary posterior.

Checkpointing writes the full state (RNG words, tree arrays, parameters, step
sizes, iteration counters, and the number of retained draws) to JSON atomically
via a temp-file rename. `resume` reconstructs an identical chain, so continuing
an interrupted run reproduces the exact trajectory it would have taken, and no
retained draw is duplicated (gate G12).
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
import json
import math
import os
import tempfile

import numpy as np

from diachron import _core
from ..models.tree import PhyloTree
from ..models.posterior import PosteriorModel
from .moves import MoveSet


@dataclass
class ChainConfig:
    seed: int
    warmup: int = 2000
    iterations: int = 8000
    thin: int = 10
    target_accept: float = 0.3


@dataclass
class ChainState:
    gain: float
    loss: float
    clock_rate: float
    geo_sigma2: float = 0.004
    # Per-candidate-edge borrowing intensities; empty unless the borrowing model
    # is active. Stored as a plain list so it serializes in checkpoints.
    borrow_rates: list = field(default_factory=list)


@dataclass
class Draw:
    iteration: int
    log_posterior: float
    root_age: float
    gain: float
    loss: float
    clock_rate: float
    newick: str


class Chain:
    def __init__(
        self,
        model: PosteriorModel,
        tree: PhyloTree,
        state: ChainState,
        config: ChainConfig,
    ):
        self.model = model
        self.tree = tree
        self.state = state
        self.config = config
        self.rng = _core.PyRng(config.seed)
        self.moves = MoveSet()
        self.iteration = 0
        self.draws: list[Draw] = []
        self._accepts = 0
        self._proposals = 0
        self._log_post = self._eval()

    def _eval(self) -> float:
        import numpy as _np

        br = _np.array(self.state.borrow_rates) if self.state.borrow_rates else None
        return self.model.log_posterior(
            self.tree,
            self.state.gain,
            self.state.loss,
            self.state.clock_rate,
            self.state.geo_sigma2,
            borrow_rates=br,
        )

    def _accept(self, log_post_new: float, log_hastings: float) -> bool:
        if not math.isfinite(log_post_new):
            return False
        log_alpha = log_post_new - self._log_post + log_hastings
        if log_alpha >= 0:
            return True
        return math.log(max(self.rng.uniform(), 1e-300)) < log_alpha

    def _step(self, adapt: bool) -> None:
        """One Metropolis-within-Gibbs step: pick a move and accept/reject."""
        move_kind = self.rng.below(6)
        self._proposals += 1
        accepted = False

        if move_kind == 0:  # scale gain
            new, res = self.moves.scale_scalar(self.state.gain, self.rng)
            if res.valid:
                old = self.state.gain
                self.state.gain = new
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.state.gain = old
        elif move_kind == 1:  # scale loss
            new, res = self.moves.scale_scalar(self.state.loss, self.rng)
            if res.valid:
                old = self.state.loss
                self.state.loss = new
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.state.loss = old
        elif move_kind == 2:  # scale clock rate
            new, res = self.moves.scale_scalar(self.state.clock_rate, self.rng)
            if res.valid:
                old = self.state.clock_rate
                self.state.clock_rate = new
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.state.clock_rate = old
        elif move_kind == 3:  # slide an internal node time
            internals = self.tree.internal_nodes()
            node = internals[self.rng.below(len(internals))]
            snapshot = list(self.tree.heights)
            res = self.moves.slide_node_time(self.tree, node, self.rng)
            if res.valid:
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.tree.heights = snapshot
        elif move_kind == 4:  # scale the whole time axis (root move)
            snapshot = list(self.tree.heights)
            res = self.moves.scale_root(self.tree, self.rng)
            if res.valid:
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.tree.heights = snapshot
        else:  # topology: NNI or SPR
            parent_snap = list(self.tree.parent)
            height_snap = list(self.tree.heights)
            if self.rng.uniform() < 0.5:
                res = self.moves.nni(self.tree, self.rng)
            else:
                res = self.moves.spr(self.tree, self.rng)
            if res.valid:
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                    accepted = True
                else:
                    self.tree.parent = parent_snap
                    self.tree.heights = height_snap
            else:
                self.tree.parent = parent_snap
                self.tree.heights = height_snap

        if accepted:
            self._accepts += 1

        # Hyperparameter moves for the active model's extra variables. These run
        # occasionally alongside the main move so the geographic diffusion rate
        # and the borrowing intensities are sampled rather than held fixed.
        self._maybe_sample_hyperparameters()

        if adapt and self._proposals % 100 == 0:
            self._adapt()

    def _maybe_sample_hyperparameters(self) -> None:
        """Sample geo_sigma2 (geography) and borrow_rates (borrowing) in place.

        Each is a scale move with the log-factor Jacobian, accepted against the
        full posterior. They are gated on the model config so a tree-only or dated
        run does nothing here.
        """
        cfg = self.model.config
        if cfg.use_geography and self.rng.uniform() < 0.25:
            new, res = self.moves.scale_scalar(self.state.geo_sigma2, self.rng)
            if res.valid:
                old = self.state.geo_sigma2
                self.state.geo_sigma2 = new
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                else:
                    self.state.geo_sigma2 = old
        if cfg.use_borrowing and self.state.borrow_rates and self.rng.uniform() < 0.25:
            idx = self.rng.below(len(self.state.borrow_rates))
            new, res = self.moves.scale_scalar(self.state.borrow_rates[idx], self.rng)
            if res.valid:
                old = self.state.borrow_rates[idx]
                self.state.borrow_rates[idx] = new
                lp = self._eval()
                if self._accept(lp, res.log_hastings):
                    self._log_post = lp
                else:
                    self.state.borrow_rates[idx] = old

    def _adapt(self) -> None:
        """Robbins-Monro step-size adaptation toward the target acceptance rate.

        Only called during warmup. Grows the proposal windows when acceptance is
        too high and shrinks them when too low, with a shrinking learning rate so
        adaptation settles.
        """
        rate = self._accepts / max(self._proposals, 1)
        lr = 0.5 / (1.0 + self._proposals / 1000.0)
        adjust = math.exp(lr * (rate - self.config.target_accept))
        self.moves.scale_lambda = min(max(self.moves.scale_lambda * adjust, 1e-3), 5.0)
        self.moves.time_window_frac = min(
            max(self.moves.time_window_frac * adjust, 1e-4), 0.9
        )

    def run(self, checkpoint_path: str | None = None, checkpoint_every: int = 1000) -> None:
        """Run warmup then sampling, checkpointing periodically."""
        total = self.config.warmup + self.config.iterations
        while self.iteration < total:
            adapt = self.iteration < self.config.warmup
            self._step(adapt)
            self.iteration += 1
            in_sampling = self.iteration > self.config.warmup
            if in_sampling and (self.iteration - self.config.warmup) % self.config.thin == 0:
                self.draws.append(
                    Draw(
                        iteration=self.iteration,
                        log_posterior=self._log_post,
                        root_age=self.tree.root_age(),
                        gain=self.state.gain,
                        loss=self.state.loss,
                        clock_rate=self.state.clock_rate,
                        newick=self.tree.newick(),
                    )
                )
            if checkpoint_path and self.iteration % checkpoint_every == 0:
                self.checkpoint(checkpoint_path)

    @property
    def acceptance_rate(self) -> float:
        return self._accepts / max(self._proposals, 1)

    # --- checkpointing -----------------------------------------------------
    def checkpoint(self, path: str) -> None:
        payload = {
            "rng_state": self.rng.snapshot(),
            "iteration": self.iteration,
            "parent": self.tree.parent,
            "heights": self.tree.heights,
            "taxon": self.tree.taxon,
            "tip_labels": {str(k): v for k, v in self.tree.tip_labels.items()},
            "state": asdict(self.state),
            "scale_lambda": self.moves.scale_lambda,
            "time_window_frac": self.moves.time_window_frac,
            "log_post": self._log_post,
            "accepts": self._accepts,
            "proposals": self._proposals,
            "n_draws": len(self.draws),
            "config": asdict(self.config),
            "draws": [asdict(d) for d in self.draws],
        }
        # Atomic write: temp file then rename.
        d = os.path.dirname(path) or "."
        os.makedirs(d, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f)
        os.replace(tmp, path)

    @classmethod
    def resume(cls, model: PosteriorModel, path: str) -> "Chain":
        with open(path) as f:
            payload = json.load(f)
        tree = PhyloTree(
            parent=list(payload["parent"]),
            heights=list(payload["heights"]),
            taxon=list(payload["taxon"]),
            tip_labels={int(k): v for k, v in payload["tip_labels"].items()},
        )
        state = ChainState(**payload["state"])
        config = ChainConfig(**payload["config"])
        chain = cls(model, tree, state, config)
        chain.rng = _core.PyRng.from_state(payload["rng_state"])
        chain.iteration = payload["iteration"]
        chain.moves.scale_lambda = payload["scale_lambda"]
        chain.moves.time_window_frac = payload["time_window_frac"]
        chain._log_post = payload["log_post"]
        chain._accepts = payload["accepts"]
        chain._proposals = payload["proposals"]
        chain.draws = [Draw(**d) for d in payload["draws"]]
        return chain

    def scalar_trace(self, name: str) -> np.ndarray:
        return np.array([getattr(d, name) for d in self.draws], dtype=float)
