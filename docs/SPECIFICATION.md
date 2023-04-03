# Diachron specification

Diachron is a research platform for studying how language families diversify,
borrow words, and shift their geographic footprint under uncertain evidence. It
infers posterior distributions over trees, divergence times, borrowing, and
latent locations from cognate data and dated attestations, and it keeps the
distinction between observation, inference, and illustration visible everywhere.

## Scope and claims

Diachron makes two kinds of statement and never blurs them.

1. **Mathematical statements** are about the algorithms and their properties: the
   pruning recursion computes the phylogenetic likelihood of a character under a
   continuous-time Markov model; a Metropolis-Hastings kernel with a stated
   proposal and acceptance ratio leaves a stated target invariant. These are
   proved against independent tiny-state enumeration and exact finite targets.

2. **Historical statements** are empirical claims with uncertainty. A posterior
   over the Indo-European root age or a latent homeland region is a distribution,
   not a recovered fact. Diachron reports credible regions, sensitivity to
   priors and calibrations, and explicit counterevidence.

A language is not a population. A word lineage is not a language tree. Neither is
a genetic ancestry tree. The product and the thesis state these distinctions.

## Units of inheritance

Observations attach to a (language, concept) cell. A cell may hold multiple
lexemes (synonyms), cognate-class judgments that disagree, partial cognates,
explicit loans, inferred absences, and genuine missingness. The transformation
from lexeme judgments to the character matrix used by the likelihood is explicit,
versioned, and auditable (`diachron.io.characters`). Missing and absent are
distinct states and produce distinct likelihoods (gate G04).

## Models

- **Character evolution (D07, D08):** original continuous-time Markov chain over
  cognate presence/absence classes, composed into a Felsenstein pruning
  likelihood with log-space scaling, multifurcations, missing states, and branch
  length edge cases. Matrix exponential is a permitted numerical primitive; the
  model composition and likelihood are original.
- **Ascertainment (D09):** a conditioning correction matching the character
  selection process (variable-only / no-constant-site selection), derived and
  validated against simulated selection (gate G04).
- **Clock and dates (D10, D11):** strict and relaxed (uncorrelated lognormal)
  clocks with identifiable parameterization; calibration densities on nodes and
  dated tips with calendar-consistent units (years before present).
- **Tree prior:** a birth or coalescent-style prior with a documented density.
- **Borrowing (D16):** a finite contact model on a bounded candidate edge set,
  with identifiability and prior-sensitivity analysis. No unrestricted network
  reconstruction is claimed.
- **Geography (D17, D18):** latent location states with uncertain tip locations
  and a transition model (region-based CTMC and a Brownian-on-sphere alternative),
  distinguishing proximity, documented contact, and hypothetical contact.

## Inference

Original Metropolis-Hastings MCMC (`diachron-core`) with branch-length, node-time,
clock-rate, and topology (NNI/SPR-style) moves, correct proposal ratios and
Jacobians, and constraint validation (D12). Chains checkpoint RNG and all state
for exact CPU continuation (D13). Diagnostics compute ESS, split-Rhat,
autocorrelation, Monte Carlo standard error, and topology mixing (D14).

## Evidence profile (G01)

- Primary study: at least 100 eligible language varieties and 100 concepts.
- Second study: at least 30 eligible varieties, with its own calibration.
- At least 100 individually sourced dated/located attestation records with
  bibliographic references and uncertainty.

Real restricted lexical records are not redistributed. The ingester reads CLDF
exactly and the public release ships a reproducible acquisition recipe plus
derived artifacts that the dataset license permits. Synthetic known-truth corpora
are always labeled synthetic and are used for recovery and calibration, never
presented as observed history.

## Gates

See `docs/BUILD_STATUS.md` for the live gate table (G01-G15).
