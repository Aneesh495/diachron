"""Replica exchange that tempers lexical likelihoods and retains cold draws."""

from __future__ import annotations

import json
import math

from diachron import _core

from .chain import Chain, write_checkpoint


def validate_ladder(values) -> tuple[float, ...]:
    ladder = tuple(values)
    if (not 1 <= len(ladder) <= 12 or ladder[0] != 1.0
            or any(not math.isfinite(value) or not .05 <= value <= 1 for value in ladder)
            or any(a <= b for a, b in zip(ladder, ladder[1:]))):
        raise ValueError("temperature ladder must start at 1 and strictly descend to at least 0.05")
    return ladder


def exchange_log_ratio(beta_a: float, beta_b: float, lexical_a: float, lexical_b: float) -> float:
    """The two unheated prior factors cancel under an exchange."""
    return (beta_a - beta_b) * (lexical_b - lexical_a)


class ReplicaExchange:
    def __init__(self, replicas: list[Chain], seed: int, swap_every: int = 20):
        self.ladder = validate_ladder(chain.config.inverse_temperature for chain in replicas)
        if len(replicas) < 2 or not 1 <= swap_every <= 1000:
            raise ValueError("replica exchange requires at least two temperatures and a bounded cadence")
        cold = replicas[0]
        if any(chain.model is not cold.model or chain.config.warmup != cold.config.warmup
               or chain.config.iterations != cold.config.iterations or chain.iteration != cold.iteration
               for chain in replicas):
            raise ValueError("replicas must share a target and iteration budget")
        if any(chain.config.thin <= chain.config.iterations for chain in replicas[1:]):
            raise ValueError("hot replicas must not retain posterior draws")
        if len({chain.config.seed for chain in replicas}) != len(replicas):
            raise ValueError("replicas require independent local RNG seeds")
        self.replicas = replicas
        self.swap_every = swap_every
        self.seed = seed
        self.rng = _core.PyRng(seed)
        self.cycles = 0
        self.attempts = [0] * (len(replicas) - 1)
        self.accepts = [0] * (len(replicas) - 1)
        self.walkers = list(range(len(replicas)))
        self.started_cold = [True] + [False] * (len(replicas) - 1)
        self.visited_hot = [False] * len(replicas)
        self.round_trips = 0

    @property
    def cold(self) -> Chain:
        chain = self.replicas[0]
        chain.exchange_statistics = self.statistics()
        return chain

    def statistics(self) -> dict:
        return {"inverse_temperatures": list(self.ladder), "swap_every": self.swap_every,
                "seed": self.seed, "swap_cycles": self.cycles,
                "adjacent_attempts": list(self.attempts), "adjacent_accepts": list(self.accepts),
                "round_trips": self.round_trips}

    @staticmethod
    def _lexical(chain: Chain) -> float:
        state = chain.state
        return sum(chain.model.character_scores(
            chain.tree, state.gain, state.loss, state.clock_rate,
            borrow_rates=state.borrow_rates, clock_multipliers=state.clock_multipliers,
        ))

    def _exchange(self) -> None:
        for left in range(self.cycles % 2, len(self.replicas) - 1, 2):
            a, b = self.replicas[left:left + 2]
            log_ratio = exchange_log_ratio(self.ladder[left], self.ladder[left + 1],
                                           self._lexical(a), self._lexical(b))
            self.attempts[left] += 1
            if log_ratio >= 0 or math.log(max(self.rng.uniform(), 1e-300)) < log_ratio:
                a.tree, b.tree = b.tree, a.tree
                a.state, b.state = b.state, a.state
                # Re-evaluate in each temperature's normal operation order so
                # restored state scores agree exactly, including roundoff.
                a._log_post = a._eval()
                b._log_post = b._eval()
                self.accepts[left] += 1
                self.walkers[left], self.walkers[left + 1] = self.walkers[left + 1], self.walkers[left]
        self.cycles += 1
        cold, hot = self.walkers[0], self.walkers[-1]
        if self.visited_hot[cold]:
            self.round_trips += 1
            self.visited_hot[cold] = False
        self.started_cold[cold] = True
        if self.started_cold[hot]:
            self.visited_hot[hot] = True

    def run(self, checkpoint_path=None, checkpoint_every: int = 1000, *, until_iteration=None) -> None:
        if checkpoint_every < 1:
            raise ValueError("checkpoint cadence must be positive")
        first = self.replicas[0]
        maximum = first.config.warmup + first.config.iterations
        limit = maximum if until_iteration is None else until_iteration
        if not first.iteration <= limit <= maximum:
            raise ValueError("stop iteration outside the configured replica run")
        while first.iteration < limit:
            next_swap = (first.iteration // self.swap_every + 1) * self.swap_every
            next_checkpoint = (first.iteration // checkpoint_every + 1) * checkpoint_every
            boundary = min(next_swap, next_checkpoint, limit)
            for replica in self.replicas:
                replica.run(until_iteration=boundary)
            # Record the cold draw before the exchange at a shared boundary.
            if boundary == next_swap:
                self._exchange()
            if checkpoint_path is not None and boundary in (next_checkpoint, limit):
                self.checkpoint(checkpoint_path)

    def checkpoint(self, path) -> None:
        write_checkpoint(path, {
            "schema": "diachron-replica-exchange-1", "replicas": [chain.checkpoint_payload() for chain in self.replicas],
            "statistics": self.statistics(), "rng_state": self.rng.snapshot(),
            "walkers": self.walkers, "started_cold": self.started_cold, "visited_hot": self.visited_hot,
        })

    @classmethod
    def resume(cls, model, path):
        with open(path, encoding="utf-8") as stream:
            payload = json.load(stream)
        if payload["schema"] != "diachron-replica-exchange-1" or not 2 <= len(payload["replicas"]) <= 12:
            raise ValueError("unsupported replica exchange checkpoint")
        statistics = payload["statistics"]
        replicas = [Chain.from_checkpoint(model, raw) for raw in payload["replicas"]]
        restored = cls(replicas, statistics["seed"], statistics["swap_every"])
        cycles = replicas[0].iteration // restored.swap_every
        attempts = [(cycles + (index % 2 == 0)) // 2 for index in range(len(replicas) - 1)]
        accepts = statistics["adjacent_accepts"]
        if (list(restored.ladder) != statistics["inverse_temperatures"] or cycles != statistics["swap_cycles"]
                or attempts != statistics["adjacent_attempts"] or len(accepts) != len(attempts)
                or any(not 0 <= accepted <= tried for accepted, tried in zip(accepts, attempts))
                or sorted(payload["walkers"]) != list(range(len(replicas)))
                or any(len(payload[key]) != len(replicas) for key in ("started_cold", "visited_hot"))
                or not 0 <= statistics["round_trips"] <= sum(accepts)
                or not any(payload["rng_state"])):
            raise ValueError("replica checkpoint counters or ladder differ")
        restored.cycles, restored.attempts, restored.accepts = cycles, attempts, accepts
        restored.round_trips = statistics["round_trips"]
        restored.walkers = payload["walkers"]
        restored.started_cold = payload["started_cold"]
        restored.visited_hot = payload["visited_hot"]
        restored.rng = _core.PyRng.from_state(payload["rng_state"])
        return restored
