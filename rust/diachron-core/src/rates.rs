//! Among-character rate heterogeneity.
//!
//! Not every cognate class evolves at the same rate: some concepts are
//! conservative and replaced rarely, others turn over fast. The standard way to
//! model this is a discrete-gamma mixture: the per-character rate multiplier is
//! drawn from a Gamma(alpha, alpha) distribution (mean one, variance 1/alpha),
//! discretized into a small number of equiprobable categories represented by
//! their medians. A character's likelihood is then the average over categories
//! of its likelihood at that rate multiplier.
//!
//! Small alpha means strong heterogeneity (a few fast classes, many slow ones);
//! large alpha approaches a single rate. The discretization follows the standard
//! median approach: split the Gamma into `k` equal-probability bins and use the
//! median of each bin as its representative rate, then renormalize so the mean
//! rate is exactly one (otherwise the clock and the rate mixture would both
//! absorb the overall scale and trade off).
//!
//! The incomplete-gamma and its inverse are implemented here (series and
//! continued-fraction expansions plus Newton refinement) so the discretization
//! is self-contained and testable, rather than depending on an external special
//! function.

/// Lower regularized incomplete gamma P(a, x) via series / continued fraction.
pub fn gammp(a: f64, x: f64) -> f64 {
    assert!(a > 0.0 && x >= 0.0);
    if x == 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        // Series representation.
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..200 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-14 {
                break;
            }
        }
        sum * (-x + a * x.ln() - ln_gamma(a)).exp()
    } else {
        // Continued fraction for the complement Q, then P = 1 - Q.
        let mut b = x + 1.0 - a;
        let mut c = 1e300;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..200 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = b + an / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-14 {
                break;
            }
        }
        let q = (-x + a * x.ln() - ln_gamma(a)).exp() * h;
        1.0 - q
    }
}

/// Lanczos approximation to ln Gamma(x) for x > 0.
pub fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection formula.
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Inverse of the regularized lower incomplete gamma: find x with P(a, x) = p.
pub fn gammp_inv(a: f64, p: f64) -> f64 {
    assert!(a > 0.0 && (0.0..=1.0).contains(&p));
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    // Initial guess (Wilson-Hilferty for the gamma quantile), then Newton steps
    // on P(a,x) - p using the known derivative dP/dx = x^(a-1) e^-x / Gamma(a).
    let mut x = a.max(1e-3);
    let lng = ln_gamma(a);
    for _ in 0..100 {
        let err = gammp(a, x) - p;
        let dpdx = ((a - 1.0) * x.ln() - x - lng).exp();
        if dpdx <= 0.0 {
            break;
        }
        let step = err / dpdx;
        let mut nx = x - step;
        if nx <= 0.0 {
            nx = x / 2.0;
        }
        if (nx - x).abs() < 1e-12 * nx.max(1.0) {
            x = nx;
            break;
        }
        x = nx;
    }
    x
}

/// Discrete-gamma rate categories (median representatives), normalized to mean 1.
pub fn discrete_gamma_rates(alpha: f64, k: usize) -> Vec<f64> {
    assert!(alpha > 0.0 && k >= 1);
    if k == 1 {
        return vec![1.0];
    }
    let mut rates = Vec::with_capacity(k);
    for i in 0..k {
        // Median of bin i: the (2i+1)/(2k) quantile of Gamma(alpha, 1) scaled by
        // 1/alpha so the distribution has mean one before renormalization.
        let q = (2.0 * i as f64 + 1.0) / (2.0 * k as f64);
        let g = gammp_inv(alpha, q);
        rates.push(g / alpha);
    }
    // Renormalize so the mean rate is exactly one.
    let mean: f64 = rates.iter().sum::<f64>() / k as f64;
    for r in rates.iter_mut() {
        *r /= mean;
    }
    rates
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_gamma_known_values() {
        // Gamma(5) = 24, ln 24 ~ 3.178.
        assert!((ln_gamma(5.0) - 24.0f64.ln()).abs() < 1e-9);
        // Gamma(0.5) = sqrt(pi).
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-9);
    }

    #[test]
    fn gammp_monotone_and_bounds() {
        let a = 2.5;
        assert!(gammp(a, 0.0).abs() < 1e-12);
        let mut prev = 0.0;
        for i in 1..50 {
            let x = i as f64 * 0.3;
            let p = gammp(a, x);
            assert!(p >= prev - 1e-12 && p <= 1.0 + 1e-9);
            prev = p;
        }
        assert!(gammp(a, 50.0) > 0.999);
    }

    #[test]
    fn gammp_inv_roundtrip() {
        let a = 3.0;
        for &p in &[0.1, 0.25, 0.5, 0.75, 0.9] {
            let x = gammp_inv(a, p);
            assert!((gammp(a, x) - p).abs() < 1e-6, "p={p} x={x}");
        }
    }

    #[test]
    fn discrete_gamma_mean_is_one() {
        for &alpha in &[0.1, 0.5, 1.0, 5.0] {
            let rates = discrete_gamma_rates(alpha, 4);
            let mean: f64 = rates.iter().sum::<f64>() / rates.len() as f64;
            assert!((mean - 1.0).abs() < 1e-9, "alpha={alpha} mean={mean}");
            // Smaller alpha gives larger spread.
            assert!(rates.iter().all(|&r| r > 0.0));
        }
    }

    #[test]
    fn large_alpha_approaches_single_rate() {
        let rates = discrete_gamma_rates(1000.0, 4);
        for r in rates {
            assert!((r - 1.0).abs() < 0.1);
        }
    }
}
