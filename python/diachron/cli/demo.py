"""The bounded end-to-end demonstration.

The demo exercises the whole pipeline on a small slice of the primary corpus so a
fresh reader can see Diachron work in a couple of minutes without network access:
build data, transform to characters, run a short four-chain MCMC on the dated
model, compute diagnostics, save a posterior bundle, and emit a small summary the
workbench can load. It deliberately uses a reduced language/concept slice and
short chains so it stays fast; it is a demonstration, not the acceptance
campaign, and it labels its output accordingly.
"""

from __future__ import annotations

import json
import os
import time

from ..io.simulate import simulate_dataset, simulate_attestations
from ..models.posterior import ModelConfig
from ..infer.campaign import prepare_study, run_campaign, campaign_bundle
from ..diagnostics import convergence as cv
from ..diagnostics import topology as tp
from ..util.gitmeta import source_commit


def run_demo(out_dir: str = "runs/demo") -> int:
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    # A small, fast slice with known truth.
    dataset, truth = simulate_dataset(
        n_languages=16,
        n_concepts=80,
        root_age=6000.0,
        gain=0.00035,
        loss=0.00085,
        borrow_prob=0.3,
        missing_frac=0.08,
        seed=20240601,
        name="demo",
    )
    attestations = simulate_attestations(dataset, truth, n_records=40, seed=20240601)
    study = prepare_study(dataset, attestations)

    print(f"[demo] dataset: {len(dataset.languages)} varieties, "
          f"{len(dataset.concepts)} concepts, {len(study.columns)} variable characters")
    print(f"[demo] {len(attestations)} dated attestations; synthetic known-truth corpus")

    config = ModelConfig.dated()
    chains = run_campaign(
        study,
        config,
        n_chains=4,
        warmup=800,
        iterations=2400,
        thin=10,
        root_age_init=truth.root_age,
    )
    accept = [round(c.acceptance_rate, 3) for c in chains]
    print(f"[demo] 4 chains done; acceptance rates {accept}")

    bundle = campaign_bundle(chains, study, config, source_commit=source_commit())
    bundle.save(os.path.join(out_dir, "dated"))

    # Diagnostics on the root age endpoint.
    import numpy as np

    root_chains = np.array(bundle.scalar("root_age"))
    rhat = cv.split_rhat(root_chains)
    bulk = cv.bulk_ess(root_chains)
    topo = tp.topology_mixing_summary(bundle.chain_newicks(0))

    summary = {
        "note": "Bounded demonstration on a labeled-synthetic slice. Not the acceptance campaign.",
        "varieties": len(dataset.languages),
        "concepts": len(dataset.concepts),
        "variable_characters": len(study.columns),
        "attestations": len(attestations),
        "chains": 4,
        "acceptance_rates": accept,
        "root_age": {
            "true": truth.root_age,
            "posterior_mean": float(root_chains.mean()),
            "posterior_sd": float(root_chains.std()),
            "split_rhat": round(float(rhat), 4),
            "bulk_ess": round(float(bulk), 1),
        },
        "topology_mixing_chain0": topo,
        "bundle_dir": os.path.join(out_dir, "dated"),
        "seconds": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "demo_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[demo] root age: true {truth.root_age:.0f}, "
          f"posterior mean {root_chains.mean():.0f} +/- {root_chains.std():.0f} years")
    print(f"[demo] split-Rhat {rhat:.3f}, bulk-ESS {bulk:.0f}")
    print(f"[demo] wrote posterior bundle to {out_dir}/dated and summary to "
          f"{out_dir}/demo_summary.json")
    print(f"[demo] finished in {time.time() - t0:.1f}s")
    return 0
