//! Substitution models for cognate-class evolution.
//!
//! The observable at a (language, concept) cell, after the character transform,
//! is a discrete state. Diachron uses two model families:
//!
//! * **Binary presence/absence** of a cognate class, with asymmetric gain and
//!   loss rates. This is the workhorse for "does language L belong to cognate
//!   class K for this concept". Gain is the birth of a new form in a lineage;
//!   loss is its replacement. Gain is typically rarer than loss, so the rates
//!   are kept distinct and both are inferred.
//!
//! * **General k-state** with an equilibrium distribution and a single exchange
//!   rate (a Felsenstein-81-style parameterization), used when a concept's forms
//!   are modeled as competing classes rather than independent binary characters.
//!
//! All rate matrices are built so that rows sum to zero and the process has the
//! stated stationary distribution. Time is in the clock's rate units; branch
//! lengths carry the calendar scaling elsewhere.

use crate::eigen::ReversibleModel;
use crate::linalg::{expm, Mat};

/// A substitution model exposes its state count, stationary distribution, and a
/// transition probability matrix over a branch of a given length.
pub trait SubstModel {
    fn n_states(&self) -> usize;
    fn stationary(&self) -> Vec<f64>;
    /// P(t): row i, column j is P(end = j | start = i) over branch length t.
    fn transition(&self, t: f64) -> Mat;
    /// The instantaneous rate matrix Q (rows sum to zero).
    fn rate_matrix(&self) -> Mat;
}

/// Reversible exchangeability model with an explicit expected event rate.
pub struct Gtr {
    model: ReversibleModel,
    generator: Mat,
    rate: f64,
}

impl Gtr {
    pub fn new(pi: Vec<f64>, exchange: &[f64], rate: f64) -> Option<Self> {
        let n = pi.len();
        if !(2..=32).contains(&n)
            || exchange.len() != n * n
            || (pi.iter().sum::<f64>() - 1.0).abs() > 1e-9
            || pi.iter().any(|p| !p.is_finite() || *p <= 0.0)
            || !rate.is_finite()
            || rate <= 0.0
        {
            return None;
        }
        let mut q = Mat::zeros(n);
        for i in 0..n {
            let mut outgoing = 0.0;
            for j in 0..n {
                if i != j {
                    let value = exchange[i * n + j];
                    if !value.is_finite() || value < 0.0 || value != exchange[j * n + i] {
                        return None;
                    }
                    let qij = value * pi[j];
                    q.set(i, j, qij);
                    outgoing += qij;
                }
            }
            q.set(i, i, -outgoing);
        }
        let mean: f64 = pi.iter().enumerate().map(|(i, p)| -p * q.at(i, i)).sum();
        if !mean.is_finite() || mean <= 0.0 {
            return None;
        }
        q.scale(rate / mean);
        Some(Self {
            model: ReversibleModel::new(pi, exchange),
            generator: q,
            rate,
        })
    }
}

impl SubstModel for Gtr {
    fn n_states(&self) -> usize {
        self.model.n_states()
    }

    fn stationary(&self) -> Vec<f64> {
        self.model.stationary().to_vec()
    }

    fn transition(&self, duration: f64) -> Mat {
        let n = self.n_states();
        if duration == 0.0 {
            return Mat::identity(n);
        }
        let data = self.model.transition(duration * self.rate);
        let rows = data
            .chunks_exact(n)
            .map(<[f64]>::to_vec)
            .collect::<Vec<_>>();
        Mat::from_rows(&rows)
    }

    fn rate_matrix(&self) -> Mat {
        self.generator.clone()
    }
}

/// Binary gain/loss model. State 0 = absent, state 1 = present.
///
/// Q = [[-gain, gain], [loss, -loss]]. Stationary dist is proportional to
/// (loss, gain): a class is present at equilibrium in proportion to how fast it
/// is gained relative to lost.
#[derive(Clone, Debug)]
pub struct BinaryGainLoss {
    pub gain: f64,
    pub loss: f64,
}

impl BinaryGainLoss {
    pub fn new(gain: f64, loss: f64) -> Self {
        assert!(gain > 0.0 && loss > 0.0, "rates must be positive");
        BinaryGainLoss { gain, loss }
    }
}

impl SubstModel for BinaryGainLoss {
    fn n_states(&self) -> usize {
        2
    }

    fn stationary(&self) -> Vec<f64> {
        let z = self.gain + self.loss;
        vec![self.loss / z, self.gain / z]
    }

    fn rate_matrix(&self) -> Mat {
        Mat::from_rows(&[vec![-self.gain, self.gain], vec![self.loss, -self.loss]])
    }

    fn transition(&self, t: f64) -> Mat {
        assert!(t >= 0.0, "branch length must be nonnegative");
        if t == 0.0 {
            return Mat::identity(2);
        }
        // Closed form for the 2-state chain; avoids expm and is exact.
        let z = self.gain + self.loss;
        let change = -(-z * t).exp_m1();
        let pi0 = self.loss / z; // stationary absent
        let pi1 = self.gain / z; // stationary present
                                 // P(t) = Pi + e * (I - Pi) structure for a reversible 2-state chain.
        Mat::from_rows(&[
            vec![1.0 - pi1 * change, pi1 * change],
            vec![pi0 * change, 1.0 - pi0 * change],
        ])
    }
}

/// Felsenstein-81-style k-state model: Q_ij = rate * pi_j for i != j, so the
/// chain jumps to state j at a rate proportional to its equilibrium frequency.
#[derive(Clone, Debug)]
pub struct F81 {
    pub pi: Vec<f64>,
    pub rate: f64,
}

impl F81 {
    pub fn new(pi: Vec<f64>, rate: f64) -> Self {
        let s: f64 = pi.iter().sum();
        assert!((s - 1.0).abs() < 1e-9, "stationary dist must sum to 1");
        assert!(pi.iter().all(|&p| p > 0.0), "frequencies must be positive");
        assert!(rate > 0.0);
        F81 { pi, rate }
    }
}

impl SubstModel for F81 {
    fn n_states(&self) -> usize {
        self.pi.len()
    }

    fn stationary(&self) -> Vec<f64> {
        self.pi.clone()
    }

    fn rate_matrix(&self) -> Mat {
        let n = self.pi.len();
        let mut q = Mat::zeros(n);
        for i in 0..n {
            let mut diag = 0.0;
            for j in 0..n {
                if i != j {
                    let r = self.rate * self.pi[j];
                    q.set(i, j, r);
                    diag -= r;
                }
            }
            q.set(i, i, diag);
        }
        q
    }

    fn transition(&self, t: f64) -> Mat {
        assert!(t >= 0.0);
        if t == 0.0 {
            return Mat::identity(self.pi.len());
        }
        // F81 has a closed form: P_ij(t) = e^{-beta t} delta_ij + (1 - e^{-beta t}) pi_j
        // with beta chosen so the expected substitution rate is `rate`.
        let n = self.pi.len();
        let beta = self.rate;
        let e = (-beta * t).exp();
        let mut p = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let kron = if i == j { 1.0 } else { 0.0 };
                p.set(i, j, e * kron + (1.0 - e) * self.pi[j]);
            }
        }
        p
    }
}

/// Fallback transition via matrix exponential, used to cross-check the closed
/// forms above and for models without a closed form.
pub fn transition_expm(model: &dyn SubstModel, t: f64) -> Mat {
    let mut q = model.rate_matrix();
    q.scale(t);
    expm(&q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_closed_form_matches_expm() {
        let m = BinaryGainLoss::new(0.3, 0.9);
        for &t in &[0.1, 0.5, 1.0, 3.0] {
            let closed = m.transition(t);
            let viaexp = transition_expm(&m, t);
            for i in 0..2 {
                for j in 0..2 {
                    assert!(
                        (closed.at(i, j) - viaexp.at(i, j)).abs() < 1e-9,
                        "mismatch t={t} i={i} j={j}"
                    );
                }
            }
        }
    }

    #[test]
    fn f81_closed_form_matches_expm() {
        let m = F81::new(vec![0.2, 0.3, 0.5], 1.4);
        for &t in &[0.2, 0.7, 2.0] {
            let closed = m.transition(t);
            let viaexp = transition_expm(&m, t);
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (closed.at(i, j) - viaexp.at(i, j)).abs() < 1e-8,
                        "mismatch t={t} i={i} j={j}"
                    );
                }
            }
        }
    }

    #[test]
    fn stationary_is_invariant() {
        let m = BinaryGainLoss::new(0.5, 1.5);
        let pi = m.stationary();
        let p = m.transition(2.3);
        // pi P = pi.
        for j in 0..2 {
            let mut acc = 0.0;
            for i in 0..2 {
                acc += pi[i] * p.at(i, j);
            }
            assert!((acc - pi[j]).abs() < 1e-10);
        }
    }
}
