# Bibliography

This bibliography attributes the ideas, algorithms, datasets, and standards
Diachron builds on. It does not imply that any published result is reproduced
here; Diachron reimplements algorithms and makes its own model comparisons with
their limitations stated.

## Algorithms and methods

- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
  likelihood approach. Journal of Molecular Evolution. The pruning recursion for
  the phylogenetic likelihood. Diachron implements the discrete-character case
  with log-space scaling and a missing-versus-absent distinction, validated
  against independent enumeration.
  Primary method reference: https://gen-files.princeton.edu/sysbio/Felsenstein81.pdf
- Felsenstein, J. (1985). Phylogenies and the comparative method. The American
  Naturalist. Independent contrasts, which Diachron uses for the Gaussian
  Brownian-on-tree geographic likelihood.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., Burkner, P. (2021).
  Rank-normalization, folding, and localization: an improved R-hat for assessing
  convergence of MCMC. Bayesian Analysis. The rank-normalized split-Rhat and
  bulk/tail ESS definitions Diachron implements.
  https://arxiv.org/abs/1903.08008
- Metropolis, N., Rosenbluth, A., Rosenbluth, M., Teller, A., Teller, E. (1953);
  Hastings, W. K. (1970). The Metropolis-Hastings algorithm underlying the
  sampler.
- Higham, N. J. (2005). The scaling and squaring method for the matrix
  exponential. The approach Diachron's `expm` uses, with a degree-6 Pade
  approximant.
- Moler, C., Van Loan, C. (2003). Nineteen dubious ways to compute the exponential
  of a matrix. Background on matrix-exponential pitfalls.

## Computational historical linguistics

- Gray, R. D., Atkinson, Q. D. (2003). Language-tree divergence times support the
  Anatolian theory of Indo-European origin. Nature. A source of the dated-tree
  framing and a competing hypothesis.
- Bouckaert, R., et al. (2012). Mapping the origins and expansion of the
  Indo-European language family. Science. Phylogeographic inference of
  Indo-European, treated as a competing analysis, not ground truth.
- Heggarty, P., et al. (2023). Language trees with sampled ancestors support a
  hybrid model for the origin of Indo-European languages. Science.
  https://doi.org/10.1126/science.abg0818 . The IE-CoR data framing and a recent
  competing chronology.

## Datasets and standards

- Glottolog (pinned release and DOI recorded per study). Versioned language
  identifiers and reference coordinates. https://glottolog.org/meta/downloads
- CLDF: Cross-Linguistic Data Formats. The lexical data contract Diachron
  ingests. https://cldf.clld.org/
- IE-CoR. Indo-European cognate relationships dataset. https://iecor.clld.org/
  and its source https://github.com/lexibank/iecor . Ingested under its own
  license; not redistributed by Diachron.
- Lexibank. A collection of standardized lexical datasets, ingested under
  per-dataset licenses.

## Attribution note

Where Diachron reports a chronology, a tree, or a homeland, it reports a posterior
with uncertainty and sensitivity, and it cites the calibration and assumptions
that produced it. It does not claim agreement with or reproduction of any of the
published analyses above; those are context and competing interpretations against
which Diachron's own, separately computed results are read.
