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
        }
    }
    if root >= n {
        return Err(PyValueError::new_err("root index out of range"));
    }
    tree.root = root;
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

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(char_log_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(block_log_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(brownian_sphere_logpdf, m)?)?;
    m.add_function(wrap_pyfunction!(haversine_km, m)?)?;
    m.add_class::<PyRng>()?;
    Ok(())
}
