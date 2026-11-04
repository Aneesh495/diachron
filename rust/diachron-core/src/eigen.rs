//! Symmetric eigendecomposition for reversible-model transition matrices.
//!
//! A general time-reversible (GTR) substitution model has a rate matrix Q that is
//! not symmetric, but reversibility means Q is similar to a symmetric matrix via
//! the stationary distribution: with Pi = diag(pi), the matrix
//! S = Pi^{1/2} Q Pi^{-1/2} is symmetric. We diagonalize S with the cyclic
//! Jacobi eigenvalue algorithm (robust and simple for the small dense symmetric
//! matrices here), then recover Q = U Lambda U^{-1} with
//! U = Pi^{-1/2} V, U^{-1} = V^T Pi^{1/2}, where V holds S's eigenvectors.
//!
//! With the eigendecomposition, the transition matrix is
//! P(t) = U exp(Lambda t) U^{-1}, which is cheaper and more stable to evaluate at
//! many branch lengths than repeating scaling-and-squaring, and it makes the
//! stationary distribution and reversibility exact by construction.
//!
//! The Jacobi rotation, the sweep schedule, and the reconstruction are
//! implemented here and validated against the matrix-exponential route.

/// Eigendecomposition of a real symmetric n-by-n matrix (row-major).
/// Returns (eigenvalues, eigenvectors-as-columns row-major).
pub fn jacobi_symmetric(a_in: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = a_in.to_vec();
    // V starts as the identity and accumulates the rotations.
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }

    for _sweep in 0..100 {
        // Off-diagonal norm; stop when negligible.
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // Compute the rotation angle.
                let theta = (aqq - app) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // Apply rotation to rows/cols p and q.
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                // Accumulate eigenvectors.
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }

    let eigvals: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    (eigvals, v)
}

/// A GTR-style reversible model specified by a stationary distribution and a
/// symmetric exchangeability matrix. Precomputes the eigendecomposition so
/// transitions at any branch length are a fast U exp(Lambda t) U^{-1}.
pub struct ReversibleModel {
    n: usize,
    pi: Vec<f64>,
    eigvals: Vec<f64>,
    // U (columns are right eigenvectors of Q) and U_inv, row-major.
    u: Vec<f64>,
    u_inv: Vec<f64>,
}

impl ReversibleModel {
    /// Build from stationary `pi` (sums to 1) and symmetric exchangeabilities
    /// `exch` (n*n, exch[i][j] = exch[j][i], diagonal ignored). The rate matrix is
    /// Q_ij = exch_ij * pi_j for i != j, scaled so the expected substitution rate
    /// is one, which is the standard GTR normalization.
    pub fn new(pi: Vec<f64>, exch: &[f64]) -> Self {
        let n = pi.len();
        assert!((pi.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        let mut q = vec![0.0; n * n];
        for i in 0..n {
            let mut diag = 0.0;
            for j in 0..n {
                if i != j {
                    let r = exch[i * n + j] * pi[j];
                    q[i * n + j] = r;
                    diag += r;
                }
            }
            q[i * n + i] = -diag;
        }
        // Normalize so sum_i pi_i * (-Q_ii) = 1 (one expected substitution / unit).
        let mu: f64 = (0..n).map(|i| -pi[i] * q[i * n + i]).sum();
        if mu > 0.0 {
            for x in q.iter_mut() {
                *x /= mu;
            }
        }

        // Symmetrize: S = Pi^{1/2} Q Pi^{-1/2}.
        let sqrt_pi: Vec<f64> = pi.iter().map(|&p| p.sqrt()).collect();
        let mut s = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                s[i * n + j] = sqrt_pi[i] * q[i * n + j] / sqrt_pi[j];
            }
        }
        // Symmetrize numerically (S should already be symmetric).
        for i in 0..n {
            for j in (i + 1)..n {
                let avg = 0.5 * (s[i * n + j] + s[j * n + i]);
                s[i * n + j] = avg;
                s[j * n + i] = avg;
            }
        }

        let (eigvals, vmat) = jacobi_symmetric(&s, n);
        // U = Pi^{-1/2} V, U_inv = V^T Pi^{1/2}.
        let mut u = vec![0.0; n * n];
        let mut u_inv = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                u[i * n + j] = vmat[i * n + j] / sqrt_pi[i];
                u_inv[i * n + j] = vmat[j * n + i] * sqrt_pi[j];
            }
        }
        ReversibleModel { n, pi, eigvals, u, u_inv }
    }

    pub fn n_states(&self) -> usize {
        self.n
    }

    pub fn stationary(&self) -> &[f64] {
        &self.pi
    }

    /// Transition matrix P(t) = U exp(Lambda t) U^{-1}, row-major.
    pub fn transition(&self, t: f64) -> Vec<f64> {
        let n = self.n;
        let exp_lambda: Vec<f64> = self.eigvals.iter().map(|&l| {
            if l >= 0.0 { 1.0 } else { (l * t).exp() }
        }).collect();
        // P = U * diag(exp_lambda) * U_inv.
        let mut p = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for k in 0..n {
                    acc += self.u[i * n + k] * exp_lambda[k] * self.u_inv[k * n + j];
                }
                // Clamp tiny negatives from roundoff.
                p[i * n + j] = acc.max(0.0);
            }
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{expm, Mat};

    #[test]
    fn jacobi_diagonalizes_known_matrix() {
        // 2x2 symmetric [[2,1],[1,2]] has eigenvalues 1 and 3.
        let (vals, _v) = jacobi_symmetric(&[2.0, 1.0, 1.0, 2.0], 2);
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 1.0).abs() < 1e-10);
        assert!((sorted[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn reversible_transition_matches_expm() {
        let pi = vec![0.2, 0.3, 0.5];
        // Symmetric exchangeabilities.
        let exch = vec![
            0.0, 1.0, 0.5,
            1.0, 0.0, 2.0,
            0.5, 2.0, 0.0,
        ];
        let model = ReversibleModel::new(pi.clone(), &exch);
        // Rebuild Q the same way for the expm cross-check.
        let n = 3;
        let mut q = vec![0.0; n * n];
        for i in 0..n {
            let mut diag = 0.0;
            for j in 0..n {
                if i != j {
                    let r = exch[i * n + j] * pi[j];
                    q[i * n + j] = r;
                    diag += r;
                }
            }
            q[i * n + i] = -diag;
        }
        let mu: f64 = (0..n).map(|i| -pi[i] * q[i * n + i]).sum();
        for x in q.iter_mut() {
            *x /= mu;
        }
        for &t in &[0.1, 0.5, 1.5] {
            let p_eig = model.transition(t);
            let mut qm = Mat { n, data: q.clone() };
            qm.scale(t);
            let p_exp = expm(&qm);
            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (p_eig[i * n + j] - p_exp.at(i, j)).abs() < 1e-8,
                        "t={t} i={i} j={j}: {} vs {}",
                        p_eig[i * n + j],
                        p_exp.at(i, j)
                    );
                }
            }
        }
    }

    #[test]
    fn reversible_rows_sum_to_one() {
        let pi = vec![0.25, 0.25, 0.25, 0.25];
        let exch = vec![
            0.0, 1.0, 1.0, 1.0,
            1.0, 0.0, 1.0, 1.0,
            1.0, 1.0, 0.0, 1.0,
            1.0, 1.0, 1.0, 0.0,
        ];
        let model = ReversibleModel::new(pi, &exch);
        let p = model.transition(0.8);
        for i in 0..4 {
            let s: f64 = (0..4).map(|j| p[i * 4 + j]).sum();
            assert!((s - 1.0).abs() < 1e-9, "row {i} sums to {s}");
        }
    }

    #[test]
    fn reversible_detailed_balance() {
        // pi_i P_ij(t) = pi_j P_ji(t) for a reversible model.
        let pi = vec![0.2, 0.3, 0.5];
        let exch = vec![0.0, 1.2, 0.7, 1.2, 0.0, 1.5, 0.7, 1.5, 0.0];
        let model = ReversibleModel::new(pi.clone(), &exch);
        let p = model.transition(0.6);
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                let lhs = pi[i] * p[i * n + j];
                let rhs = pi[j] * p[j * n + i];
                assert!((lhs - rhs).abs() < 1e-9, "detailed balance i={i} j={j}");
            }
        }
    }
}
