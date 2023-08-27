//! PyO3 bindings exposing the core to Python.
//!
//! The Python orchestration layer builds trees and character matrices, then
//! calls into this module for the heavy inner loops: the pruning likelihood and
//! the MCMC moves. The binding keeps the interface explicit and flat (indices
//! and arrays, no opaque handles leaking across the boundary) so the contract is
//! easy to audit and serialize.
//!
//! Trees cross the boundary as parallel arrays: `parent[i]` (negative for the
//! root), `branch_len[i]`, `height[i]`, and `taxon[i]` (negative for internal
//! nodes). Character columns are lists of per-taxon states where a negative
//! value means missing.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use crate::pruning::{log_likelihood, log_likelihood_block, TipStates};
use crate::subst::{BinaryGainLoss, F81, SubstModel};
use crate::tree::Tree;
use crate::rng::Rng;
use crate::geo;

fn build_tree(
    parent: &[i64],
    branch_len: &[f64],
    height: &[f64],
    taxon: &[i64],
    root: usize,
) -> PyResult<Tree> {
    let n = parent.len();
    if branch_len.len() != n || height.len() != n || taxon.len() != n {
        return Err(PyValueError::new_err("tree arrays length mismatch"));
    }
    if n == 0 || n > 10_000 || root >= n || parent[root] != -1 {
        return Err(PyValueError::new_err("invalid tree root or node count"));
    }
    if branch_len.iter().chain(height).any(|x| !x.is_finite() || *x < 0.0) {
        return Err(PyValueError::new_err("invalid tree branch length or height"));
    }
    // Construct directly from the arrays, tolerant of arbitrary node ordering
    // (a parent may be listed after its children). We fill nodes first, then
    // wire up the children lists in a second pass.
    let mut tree = Tree::new();
    for i in 0..n {
        tree.nodes.push(crate::tree::Node {
            parent: if parent[i] < 0 { None } else { Some(parent[i] as usize) },
            children: Vec::new(),
            branch_len: branch_len[i],
            taxon: if taxon[i] < 0 { None } else { Some(taxon[i] as usize) },
            height: height[i],
        });
    }
    for i in 0..n {
        if let Some(pp) = tree.nodes[i].parent {
            if pp >= n {
                return Err(PyValueError::new_err("parent index out of range"));
            }
            tree.nodes[pp].children.push(i);
        } else if i != root {
            return Err(PyValueError::new_err("multiple tree roots"));
        }
    }
    tree.root = root;
    if tree.postorder().len() != n {
        return Err(PyValueError::new_err("tree has unreachable nodes or a cycle"));
    }
    let mut tips = vec![false; tree.n_tips()];
    for node in &tree.nodes {
        match (node.children.is_empty(), node.taxon) {
            (true, Some(t)) if t < tips.len() && !tips[t] => tips[t] = true,
            (false, None) => (),
            _ => return Err(PyValueError::new_err("invalid or duplicate tree taxon index")),
        }
    }
    Ok(tree)
}

fn make_model(kind: &str, params: &[f64]) -> PyResult<Box<dyn SubstModel>> {
    match kind {
        "binary" => {
            if params.len() != 2 {
                return Err(PyValueError::new_err("binary model needs [gain, loss]"));
            }
            Ok(Box::new(BinaryGainLoss::new(params[0], params[1])))
        }
        "f81" => {
            if params.len() < 2 {
                return Err(PyValueError::new_err("f81 needs [rate, pi0, pi1, ...]"));
            }
            let rate = params[0];
            let pi = params[1..].to_vec();
            Ok(Box::new(F81::new(pi, rate)))
        }
        _ => Err(PyValueError::new_err(format!("unknown model kind {kind}"))),
    }
}

/// Log-likelihood of a single character column.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn char_log_likelihood(
    parent: Vec<i64>,
    branch_len: Vec<f64>,
    height: Vec<f64>,
    taxon: Vec<i64>,
    root: usize,
    model_kind: String,
    model_params: Vec<f64>,
    column: Vec<i64>,
    root_prior: Vec<f64>,
) -> PyResult<f64> {
    let tree = build_tree(&parent, &branch_len, &height, &taxon, root)?;
    let model = make_model(&model_kind, &model_params)?;
    let states: Vec<Option<usize>> = column
        .iter()
        .map(|&s| if s < 0 { None } else { Some(s as usize) })
        .collect();
    let tips = TipStates { states: &states };
    Ok(log_likelihood(&tree, model.as_ref(), &tips, &root_prior))
}

/// Log-likelihood of a block of character columns sharing one tree and model.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn block_log_likelihood(
    parent: Vec<i64>,
    branch_len: Vec<f64>,
    height: Vec<f64>,
    taxon: Vec<i64>,
    root: usize,
    model_kind: String,
    model_params: Vec<f64>,
    columns: Vec<Vec<i64>>,
    root_prior: Vec<f64>,
) -> PyResult<f64> {
    let tree = build_tree(&parent, &branch_len, &height, &taxon, root)?;
    let model = make_model(&model_kind, &model_params)?;
    let cols: Vec<Vec<Option<usize>>> = columns
        .iter()
        .map(|col| {
            col.iter()
                .map(|&s| if s < 0 { None } else { Some(s as usize) })
                .collect()
        })
        .collect();
    Ok(log_likelihood_block(&tree, model.as_ref(), &cols, &root_prior))
}

/// Binary character scores with observed-mask conditioning and a site rate mixture.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (parent, branch_len, height, taxon, root, gain, loss, columns,
                    rate_categories, min_present, contacts=None))]
fn binary_character_scores(
    parent: Vec<i64>,
    branch_len: Vec<f64>,
    height: Vec<f64>,
    taxon: Vec<i64>,
    root: usize,
    gain: f64,
    loss: f64,
    columns: Vec<Vec<i64>>,
    rate_categories: Vec<f64>,
    min_present: usize,
    contacts: Option<Vec<(usize, usize, f64)>>,
) -> PyResult<Vec<f64>> {
    let tree = build_tree(&parent, &branch_len, &height, &taxon, root)?;
    let n_tips = tree.n_tips();
    let contacts = contacts.unwrap_or_default();
    if contacts.len() > 16
        || contacts.iter().any(|&(donor, recipient, intensity)| {
            donor >= n_tips || recipient >= n_tips || donor == recipient
                || !intensity.is_finite() || intensity < 0.0
        })
        || !contacts.iter().map(|(_, _, w)| w).sum::<f64>().is_finite()
        || min_present > 2
        || rate_categories.is_empty()
        || rate_categories.len() > 32
        || rate_categories.iter().any(|r| !r.is_finite() || *r <= 0.0)
        || !gain.is_finite() || gain <= 0.0
        || !loss.is_finite() || loss <= 0.0
        || !(gain + loss).is_finite()
        || columns.iter().any(|col| {
            col.len() != n_tips || col.iter().any(|s| !(-1..=1).contains(s))
        })
    {
        return Err(PyValueError::new_err("invalid binary mixture rates or character states"));
    }
    let model = BinaryGainLoss::new(gain, loss);
    let cols: Vec<Vec<Option<usize>>> = columns.iter().map(|col| {
        col.iter().map(|&s| usize::try_from(s).ok()).collect()
    }).collect();
    Ok(crate::pruning::binary_mixture_scores(
        &tree, &model, &cols, &model.stationary(), &rate_categories, min_present, &contacts,
    ))
}

/// A checkpointable RNG exposed to Python for reproducible sampling.
#[pyclass]
struct PyRng {
    inner: Rng,
}

#[pymethods]
impl PyRng {
    #[new]
    fn new(seed: u64) -> Self {
        PyRng { inner: Rng::seed_from_u64(seed) }
    }

    #[staticmethod]
    fn from_state(state: Vec<u64>) -> PyResult<Self> {
        if state.len() != 4 {
            return Err(PyValueError::new_err("rng state must have 4 words"));
        }
        Ok(PyRng {
            inner: Rng::from_state([state[0], state[1], state[2], state[3]]),
        })
    }

    fn snapshot(&self) -> Vec<u64> {
        self.inner.snapshot().to_vec()
    }

    fn uniform(&mut self) -> f64 {
        self.inner.next_f64()
    }

    fn normal(&mut self) -> f64 {
        self.inner.next_normal()
    }

    fn below(&mut self, n: usize) -> usize {
        self.inner.next_below(n)
    }
}

/// Brownian-on-sphere log density, exposed for the geographic likelihood.
#[pyfunction]
fn brownian_sphere_logpdf(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    sigma2_per_time: f64,
    t: f64,
) -> f64 {
    geo::brownian_sphere_logpdf(lat1, lon1, lat2, lon2, sigma2_per_time, t)
}

/// Great-circle distance in km.
#[pyfunction]
fn haversine_km(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    geo::haversine_km(lat1, lon1, lat2, lon2)
}

/// Discrete-gamma rate categories (median representatives, mean normalized to 1).
#[pyfunction]
fn discrete_gamma_rates(alpha: f64, k: usize) -> Vec<f64> {
    crate::rates::discrete_gamma_rates(alpha, k)
}

/// Run the native Rust MCMC chain and return thinned scalar draws.
///
/// Returns a list of (iteration, log_post, root_age, gain, loss) tuples. The tree
/// crosses as the usual parallel arrays; node heights carry the times.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn run_native_chain(
    parent: Vec<i64>,
    height: Vec<f64>,
    taxon: Vec<i64>,
    root: usize,
    gain: f64,
    loss: f64,
    columns: Vec<Vec<i64>>,
    pop_size: f64,
    seed: u64,
    warmup: u64,
    iterations: u64,
    thin: u64,
) -> PyResult<Vec<(u64, f64, f64, f64, f64)>> {
    // Build the tree (branch lengths derived from heights).
    let n = parent.len();
    if height.len() != n || taxon.len() != n {
        return Err(PyValueError::new_err("tree arrays length mismatch"));
    }
    let mut tree = crate::tree::Tree::new();
    for i in 0..n {
        tree.nodes.push(crate::tree::Node {
            parent: if parent[i] < 0 { None } else { Some(parent[i] as usize) },
            children: Vec::new(),
            branch_len: 0.0,
            taxon: if taxon[i] < 0 { None } else { Some(taxon[i] as usize) },
            height: height[i],
        });
    }
    for i in 0..n {
        if let Some(p) = tree.nodes[i].parent {
            if p >= n {
                return Err(PyValueError::new_err("parent index out of range"));
            }
            tree.nodes[p].children.push(i);
        }
    }
    tree.root = root;
    tree.set_branch_lengths_from_heights();

    let cols: Vec<Vec<Option<usize>>> = columns
        .iter()
        .map(|c| c.iter().map(|&s| if s < 0 { None } else { Some(s as usize) }).collect())
        .collect();
    let target = crate::chain::Target { columns: cols, pop_size };
    let mut chain = crate::chain::NativeChain::new(tree, gain, loss, seed, target);
    chain.run(warmup, iterations, thin);
    Ok(chain
        .draws
        .iter()
        .map(|d| (d.iteration, d.log_post, d.root_age, d.gain, d.loss))
        .collect())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(char_log_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(block_log_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(binary_character_scores, m)?)?;
    m.add_function(wrap_pyfunction!(brownian_sphere_logpdf, m)?)?;
    m.add_function(wrap_pyfunction!(haversine_km, m)?)?;
    m.add_function(wrap_pyfunction!(discrete_gamma_rates, m)?)?;
    m.add_function(wrap_pyfunction!(run_native_chain, m)?)?;
    m.add_class::<PyRng>()?;
    Ok(())
}
