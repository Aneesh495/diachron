//! Diachron core: original phylogenetic likelihood and MCMC kernels.
//!
//! This crate is the computational heart of Diachron. It implements, from
//! scratch, the continuous-time Markov substitution models, the Felsenstein
//! pruning likelihood, an independent enumeration oracle used to validate that
//! likelihood, and the Metropolis-Hastings kernels used for posterior sampling.
//! Matrix exponentiation is the only numerical primitive that leans on a
//! standard algorithm (scaling and squaring with a Pade approximant); it too is
//! implemented here rather than pulled from a library, so the whole inference
//! path is auditable.
//!
//! The crate is consumed two ways: as a pure-Rust library (unit and property
//! tests run against it directly) and as a Python extension module via PyO3,
//! which is how the orchestration layer in `python/diachron` drives it.

pub mod linalg;
pub mod subst;
pub mod tree;
pub mod pruning;
pub mod enumerate;
pub mod rng;
pub mod mcmc;
pub mod geo;
pub mod rates;
pub mod chain;
pub mod eigen;

#[cfg(feature = "python")]
mod bindings;
