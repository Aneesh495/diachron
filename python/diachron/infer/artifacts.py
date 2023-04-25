"""Posterior artifact storage.

A completed inference run is saved as a self-describing bundle: per-chain draws,
the model configuration, the dataset identity, and a summary. The workbench and
the bounded serving path read these bundles; they never trigger new inference.
The format is plain JSON plus, for the larger per-draw arrays, a columnar store,
so a saved posterior can be inspected offline without the sampler.

A bundle records the executable source commit separately from the results, so a
later report-only write can publish results while still pointing at the exact
code that produced them.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
import json
import os
import hashlib

from .chain import Chain


@dataclass
class PosteriorBundle:
    model_name: str
    dataset_name: str
    dataset_version: str
    is_synthetic: bool
    n_chains: int
    draws_per_chain: int
    tip_labels: list[str]
    # chains[i] is a list of Draw dicts.
    chains: list[list[dict]]
    source_commit: str = "unknown"
    notes: str = ""
    summary: dict = field(default_factory=dict)

    def save(self, directory: str) -> None:
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "bundle.json"), "w", encoding="utf-8") as f:
            json.dump(asdict(self), f)
        # A separate checksum over the draws, so tampering with results is
        # detectable independently of the bundle's own fields.
        digest = self.draw_digest()
        with open(os.path.join(directory, "MANIFEST.json"), "w", encoding="utf-8") as f:
            json.dump(
                {
                    "model_name": self.model_name,
                    "dataset_name": self.dataset_name,
                    "n_chains": self.n_chains,
                    "draws_per_chain": self.draws_per_chain,
                    "draw_digest": digest,
                    "source_commit": self.source_commit,
                },
                f,
                indent=2,
            )

    def draw_digest(self) -> str:
        """A stable hash of the draw payload, used to detect result tampering."""
        h = hashlib.sha256()
        for chain in self.chains:
            for d in chain:
                h.update(
                    f"{d['iteration']}|{d['log_posterior']:.9g}|{d['root_age']:.9g}|"
                    f"{d['gain']:.9g}|{d['loss']:.9g}|{d['clock_rate']:.9g}|{d['newick']}".encode()
                )
        return h.hexdigest()

    @classmethod
    def load(cls, directory: str) -> "PosteriorBundle":
        with open(os.path.join(directory, "bundle.json"), encoding="utf-8") as f:
            payload = json.load(f)
        return cls(**payload)

    def chain_newicks(self, chain_index: int) -> list[str]:
        return [d["newick"] for d in self.chains[chain_index]]

    def scalar(self, name: str) -> list[list[float]]:
        return [[d[name] for d in chain] for chain in self.chains]


def bundle_from_chains(
    chains: list[Chain],
    model_name: str,
    dataset_name: str,
    dataset_version: str,
    is_synthetic: bool,
    tip_labels: list[str],
    source_commit: str = "unknown",
) -> PosteriorBundle:
    draw_lists = [[asdict(d) for d in c.draws] for c in chains]
    return PosteriorBundle(
        model_name=model_name,
        dataset_name=dataset_name,
        dataset_version=dataset_version,
        is_synthetic=is_synthetic,
        n_chains=len(chains),
        draws_per_chain=min((len(c.draws) for c in chains), default=0),
        tip_labels=tip_labels,
        chains=draw_lists,
        source_commit=source_commit,
    )
