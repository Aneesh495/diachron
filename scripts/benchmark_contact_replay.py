"""Measure lexical scoring on frozen draws from an observed contact bundle."""

import argparse
import hashlib
import json
from pathlib import Path
from statistics import median
import struct
import time

from diachron import _core
from diachron.models.tree import PhyloTree
from diachron.verify.verifier import _model, _read_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bundle", type=Path)
    args = parser.parse_args()
    bundle = _read_json(args.bundle / "bundle.json", 128 * 1024 * 1024)
    model = _model(bundle["context"]["model"])
    if not model.config.use_borrowing or not model.borrow_candidates:
        raise ValueError("benchmark requires a fitted contact model")
    scores, times = [], []
    for chain in bundle["chains"]:
        for draw in (chain[len(chain) // 4], chain[3 * len(chain) // 4]):
            tree = PhyloTree(draw["parent"], draw["heights"], draw["taxon"])
            started = time.perf_counter()
            for _ in range(10):
                result = model.character_scores(tree, draw["gain"], draw["loss"], draw["clock_rate"],
                                                borrow_rates=draw["borrow_rates"],
                                                clock_multipliers=draw["clock_multipliers"])
            times.append((time.perf_counter() - started) / 10)
            scores.extend(result)
    print(json.dumps({"source_commit": bundle["source_commit"], "characters": len(model.columns),
                      "contact_edges": len(model.borrow_candidates), "states": len(times), "repetitions": 10,
                      "native_sha256": hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
                      "score_sha256": hashlib.sha256(struct.pack(f"<{len(scores)}d", *scores)).hexdigest(),
                      "seconds_per_call": times, "median_seconds": median(times)}, indent=2))


if __name__ == "__main__":
    main()
