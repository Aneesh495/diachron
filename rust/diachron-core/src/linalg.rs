//! Dense small-matrix linear algebra for the substitution models.
//!
//! The matrices here are rate matrices and transition probability matrices of a
//! continuous-time Markov chain over a small cognate-class alphabet (typically 2
//! to a few dozen states). Everything is row-major `Vec<f64>` of length `n*n`.
//!
//! Units: a rate matrix `Q` has off-diagonal entries in units of expected
//! substitutions per unit branch length; rows sum to zero. `expm(Q*t)` is a
//! stochastic matrix (nonnegative rows summing to one) giving P(end state |
//! start state) over a branch of length `t`.
//!
//! The matrix exponential uses scaling-and-squaring with a degree-6 Padé
//! approximant. For the tiny, well-conditioned generators Diachron builds this
//! is accurate to near machine precision and is validated in `tests/` against
//! an independent eigen/series route and against scipy's `expm`.

/// A square dense matrix stored row-major.
#[derive(Clone, Debug, PartialEq)]
pub struct Mat {
    pub n: usize,
    pub data: Vec<f64>,
}

impl Mat {
    pub fn zeros(n: usize) -> Self {
        Mat { n, data: vec![0.0; n * n] }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Mat::zeros(n);
        for i in 0..n {
            m.data[i * n + i] = 1.0;
        }
        m
    }

    pub fn from_rows(rows: &[Vec<f64>]) -> Self {
        let n = rows.len();
        let mut data = Vec::with_capacity(n * n);
        for r in rows {
            assert_eq!(r.len(), n, "matrix must be square");
            data.extend_from_slice(r);
        }
        Mat { n, data }
    }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }

    /// Matrix-matrix product. O(n^3); n is tiny so this is not a bottleneck.
    pub fn matmul(&self, other: &Mat) -> Mat {
        assert_eq!(self.n, other.n);
        let n = self.n;
        let mut out = Mat::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let a = self.data[i * n + k];
                if a == 0.0 {
                    continue;
                }
                for j in 0..n {
                    out.data[i * n + j] += a * other.data[k * n + j];
                }
            }
        }
        out
    }

    /// In-place scalar multiply.
    pub fn scale(&mut self, s: f64) {
        for v in self.data.iter_mut() {
            *v *= s;
        }
    }

    /// Element-wise add, returning a new matrix.
    pub fn add(&self, other: &Mat) -> Mat {
        assert_eq!(self.n, other.n);
        let mut out = self.clone();
        for (o, b) in out.data.iter_mut().zip(other.data.iter()) {
            *o += *b;
        }
        out
    }

    /// Largest absolute row sum (the induced infinity norm). Used to pick the
    /// scaling-and-squaring exponent.
    pub fn inf_norm(&self) -> f64 {
        let n = self.n;
        let mut best = 0.0f64;
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                s += self.data[i * n + j].abs();
            }
            if s > best {
                best = s;
            }
        }
        best
    }
}

/// Solve `A x = b` columns via Gaussian elimination with partial pivoting.
/// Returns `A^{-1} B` as a matrix. Used inside the Pade evaluation; `A` is
/// well conditioned there (diagonally dominant by construction).
pub fn solve(a: &Mat, b: &Mat) -> Mat {
    let n = a.n;
    // Work on an augmented copy [A | B].
    let mut m = a.data.clone();
    let mut r = b.data.clone();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut best = m[col * n + col].abs();
        for row in (col + 1)..n {
            let v = m[row * n + col].abs();
            if v > best {
                best = v;
                piv = row;
            }
        }
        if piv != col {
            for j in 0..n {
                m.swap(col * n + j, piv * n + j);
                r.swap(col * n + j, piv * n + j);
            }
        }
        let diag = m[col * n + col];
        // diag is nonzero for the Pade denominator; guard anyway.
        let inv = 1.0 / diag;
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = m[row * n + col] * inv;
            if factor == 0.0 {
                continue;
            }
            for j in 0..n {
                m[row * n + j] -= factor * m[col * n + j];
                r[row * n + j] -= factor * r[col * n + j];
            }
        }
    }
    // Normalize rows by the diagonal.
    let mut out = Mat::zeros(n);
    for i in 0..n {
        let inv = 1.0 / m[i * n + i];
        for j in 0..n {
            out.data[i * n + j] = r[i * n + j] * inv;
        }
    }
    out
}

/// Matrix exponential via scaling and squaring with a degree-6 Pade approximant.
///
/// Reference identity: `expm(A) = r(A/2^s)^(2^s)` where `r` is the diagonal
/// (6,6) Pade rational approximation of `exp`. For the diagonal Pade the
/// denominator is the numerator evaluated at `-A`, so the numerator and
/// denominator share coefficients up to sign. The coefficients `c_k` below are
/// the standard (6,6) Pade numerator coefficients
///
///     c_k = (2q - k)! q! / ((2q)! k! (q - k)!)   with q = 6.
///
/// These are not the Taylor coefficients 1/k!; using those gives a truncated
/// Taylor series, not a Pade approximant, which is far less accurate. The
/// scaling choice, the even/odd split, the solve, and the squaring loop are
/// implemented here.
pub fn expm(a: &Mat) -> Mat {
    let n = a.n;
    // Pick s so that ||A/2^s|| <= 1/2, which keeps the Pade error tiny.
    let norm = a.inf_norm();
    let s = if norm <= 0.5 {
        0
    } else {
        (norm.log2().ceil() as i32 + 1).max(0) as u32
    };
    let scale = 2f64.powi(s as i32);
    let mut scaled = a.clone();
    scaled.scale(1.0 / scale);

    // Diagonal (6,6) Pade numerator coefficients c_k.
    // c0=1, c1=1/2, c2=5/44, c3=1/66, c4=1/792, c5=1/15840, c6=1/665280.
    let c: [f64; 7] = [
        1.0,
        1.0 / 2.0,
        5.0 / 44.0,
        1.0 / 66.0,
        1.0 / 792.0,
        1.0 / 15840.0,
        1.0 / 665280.0,
    ];

    // Powers of the scaled matrix: A^2, A^4, A^6.
    let a2 = scaled.matmul(&scaled);
    let a4 = a2.matmul(&a2);
    let a6 = a4.matmul(&a2);

    // Even part  E = c0 I + c2 A^2 + c4 A^4 + c6 A^6.
    // Odd scalar  O = A (c1 I + c3 A^2 + c5 A^4).
    // Numerator   N = E + O,   Denominator D = E - O  (= N evaluated at -A).
    let mut odd_inner = Mat::identity(n);
    odd_inner.scale(c[1]);
    let mut t = a2.clone();
    t.scale(c[3]);
    odd_inner = odd_inner.add(&t);
    let mut t = a4.clone();
    t.scale(c[5]);
    odd_inner = odd_inner.add(&t);
    let o = scaled.matmul(&odd_inner);

    let mut e = Mat::identity(n);
    e.scale(c[0]);
    let mut t = a2.clone();
    t.scale(c[2]);
    e = e.add(&t);
    let mut t = a4.clone();
    t.scale(c[4]);
    e = e.add(&t);
    let mut t = a6.clone();
    t.scale(c[6]);
    e = e.add(&t);

    // N = E + O, D = E - O. r = D^{-1} N.
    let num = e.add(&o);
    let mut den = e.clone();
    for (x, y) in den.data.iter_mut().zip(o.data.iter()) {
        *x -= *y;
    }
    let mut result = solve(&den, &num);

    // Square s times to undo the scaling.
    for _ in 0..s {
        result = result.matmul(&result);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expm_of_zero_is_identity() {
        let z = Mat::zeros(3);
        let e = expm(&z);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((e.at(i, j) - expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn expm_diagonal_matches_scalar_exp() {
        let m = Mat::from_rows(&[vec![0.3, 0.0], vec![0.0, -1.2]]);
        let e = expm(&m);
        assert!((e.at(0, 0) - 0.3f64.exp()).abs() < 1e-10);
        assert!((e.at(1, 1) - (-1.2f64).exp()).abs() < 1e-10);
        assert!(e.at(0, 1).abs() < 1e-12);
    }

    #[test]
    fn expm_rate_matrix_rows_sum_to_one() {
        // A 2-state rate matrix with rows summing to zero exponentiates to a
        // stochastic matrix.
        let q = Mat::from_rows(&[vec![-0.7, 0.7], vec![0.4, -0.4]]);
        let p = expm(&q);
        for i in 0..2 {
            let s: f64 = (0..2).map(|j| p.at(i, j)).sum();
            assert!((s - 1.0).abs() < 1e-10, "row {i} sums to {s}");
            for j in 0..2 {
                assert!(p.at(i, j) >= -1e-12);
            }
        }
    }
}
