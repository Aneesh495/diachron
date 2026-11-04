# Data card

## Corpora shipped with Diachron

### Primary reference corpus (`diachron_primary`)

- Kind: labeled synthetic known-truth. Not observed history.
- Dimensions: 112 language varieties, 110 concepts, 11333 observed cells, 1576
  variable binary characters after singleton exclusion, 140 dated attestations.
- Generation: `io.simulate.simulate_dataset` with seed 20240601, a coalescent
  tree at root age 6000 years, gain 0.00035 and loss 0.00085 per year, borrowing
  probability 0.25 per concept, and 8 percent missing cells.
- License: synthetic, no redistribution restrictions.
- Intended use: recovery and coverage studies, pipeline demonstration, and
  workbench artifacts. Known truth is recorded so inference can be scored against
  it.
- Not intended for: any claim about a real language family's history.

### Second reference corpus (`diachron_second`)

- Kind: labeled synthetic known-truth.
- Dimensions: 36 varieties, 90 concepts, 60 attestations; a younger root (3500
  years) with its own calibration regime.
- Generation: seed 20240602; see `io.corpus.build_second_corpus`.
- Purpose: test generalization to a second family without forcing the primary
  family's calibration assumptions.

## Real-data path

Diachron ingests CLDF collections (for example IE-CoR / Lexibank) that the user
acquires under the collection's own license. The reader
(`io.cldf.load_cldf_wordlist`) is schema-driven and records the dataset license
and DOI. Restricted lexical records are never redistributed by this repository;
the public release ships the acquisition recipe and derived artifacts the license
permits.

## Known biases and limitations

- The synthetic corpora are generated under the inference model's own
  assumptions, which is ideal for recovery but means a clean fit on them does not
  prove fit on real data. Deliberate misspecification is injected in the recovery
  regimes (borrowing, missingness) to probe this.
- Attestations are deliberately sparse and region-biased so that inference must
  reflect the bias rather than assume even coverage.
