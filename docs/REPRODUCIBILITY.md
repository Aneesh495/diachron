# Reproducibility

## Bounded demonstration

```bash
make bootstrap
make build
make test
make likelihood-check
make demo
make verify EVIDENCE=runs/demo VERIFY_ARGS="--expected-commit $(git rev-parse HEAD)"
make thesis
```

Python tooling is pinned in `requirements-dev.txt`. The build uses
CPython 3.12 to 3.14 and Rust 1.88.0. The workbench uses Node 24.21.0
and the committed npm lockfile. The thesis builder writes Markdown and
standalone HTML under `runs/thesis`.

## Observed diagnostic

Acquire the pinned inputs from [Observed data](OBSERVED_DATA.md), then run

```bash
.venv/bin/diachron sample --study iecor --model dated --chains 4 --jobs 4 \
  --warmup 400 --iterations 1200 --thin 10 --out runs/observed-iecor \
  --source data/cache/iecor --catalog data/cache/glottolog53.zip \
  --edh-source data/cache/edh --edh-geography data/cache/edh-geography.json

.venv/bin/diachron diagnose --bundle runs/observed-iecor/dated
.venv/bin/diachron verify --evidence runs/observed-iecor/dated \
  --expected-commit "$(git rev-parse HEAD)" \
  --source data/cache/iecor --catalog data/cache/glottolog53.zip \
  --edh-source data/cache/edh --edh-geography data/cache/edh-geography.json
```

The bundle stores exact tree arrays, rates, input columns and attestation
context. Verification checks recorded source identity, replays stored scores
and, with the pinned inputs supplied, rebuilds the observation context.
An honest bundle may pass replay and miss convergence.

## Profiles and recovery

`make acceptance RUN_ID=<id> PROFILE=cpu` runs the synthetic diagnostic
into a new evidence directory. The full profile records a blocked status
until the complete observed campaign and verifier are available.

To resume a sampling command, keep its source checkout and its original
arguments, including output directory, and add `--resume`. The chain
rejects stale source, inputs or budget. Per-chain checkpoints preserve
the RNG words, tree state, frozen step sizes and retained draw order.

Run `diachron recovery --datasets 200 --out <new-directory>` for the revised
known-truth protocol. The public results include raw per-dataset draws and
clade masses under [Recovery](results/recovery_200.md).
