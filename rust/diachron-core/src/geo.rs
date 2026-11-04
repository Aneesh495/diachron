//! Geographic transition models for latent lineage locations.
//!
//! Diachron treats the location of a lineage as a latent state that changes
//! along branches. Two models are offered and compared rather than assumed:
//!
//! * **Region CTMC**: a continuous-time Markov chain over a small set of named
//!   regions, with a symmetric or asymmetric transition rate matrix. This is
//!   appropriate when the geographic evidence is coarse (a language is "in
//!   Anatolia" or "in the Pontic steppe") and avoids pretending to resolve
//!   coordinates that the data cannot support.
//!
//! * **Brownian on the sphere**: lineage location diffuses on the surface of the
//!   Earth. The transition kernel over a branch of length t and diffusion
//!   variance sigma^2 is approximated by a wrapped Gaussian in tangent space; we
//!   expose the log density of a displacement so the likelihood can integrate
//!   over uncertain tip coordinates.
//!
//! The region CTMC reuses the substitution machinery (it is a CTMC); the
//! spherical model is implemented here. Neither model is told that modern
//! coordinates are ancient homelands: tip locations enter as uncertain
//! observations with their own error, which the caller supplies.

/// Great-circle (haversine) distance in kilometers between two lat/lon points in
/// degrees. Earth radius is the mean radius; this is a geometric helper, not a
/// geodesy claim.
pub fn haversine_km(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    const R: f64 = 6371.0088;
    let to_rad = std::f64::consts::PI / 180.0;
    let (p1, p2) = (lat1 * to_rad, lat2 * to_rad);
    let dp = (lat2 - lat1) * to_rad;
    let dl = (lon2 - lon1) * to_rad;
    let a = (dp / 2.0).sin().powi(2) + p1.cos() * p2.cos() * (dl / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

/// Log density of Brownian displacement on the sphere, approximated in the
/// tangent plane at the origin point. `sigma2_per_time` is the diffusion
/// variance in km^2 per unit branch length. Over short-to-moderate branches and
/// the continental scales Diachron works at this tangent approximation is
/// adequate; it is cross-checked against a direct numerical integration in the
/// Python test suite.
pub fn brownian_sphere_logpdf(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    sigma2_per_time: f64,
    t: f64,
) -> f64 {
    assert!(t > 0.0 && sigma2_per_time > 0.0);
    let d = haversine_km(lat1, lon1, lat2, lon2);
    let var = sigma2_per_time * t;
    // Isotropic 2D Gaussian in the tangent plane: the displacement magnitude d
    // has density proportional to exp(-d^2 / (2 var)) with a 2D normalizer.
    let norm = -(2.0 * std::f64::consts::PI * var).ln();
    norm - (d * d) / (2.0 * var)
}

/// Build a symmetric region transition rate matrix from pairwise "ease of
/// movement" weights. `weights[i][j]` is a nonnegative rate of moving from i to
/// j; the diagonal is set so rows sum to zero. Asymmetric weights are allowed.
pub fn region_rate_matrix(weights: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = weights.len();
    let mut q = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut diag = 0.0;
        for j in 0..n {
            if i != j {
                q[i][j] = weights[i][j];
                diag += weights[i][j];
            }
        }
        q[i][i] = -diag;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn haversine_known_distance() {
        // London to Paris is about 344 km.
        let d = haversine_km(51.5074, -0.1278, 48.8566, 2.3522);
        assert!((d - 344.0).abs() < 15.0, "got {d}");
    }

    #[test]
    fn brownian_peaks_at_zero_displacement() {
        let near = brownian_sphere_logpdf(50.0, 10.0, 50.0, 10.0001, 1000.0, 1.0);
        let far = brownian_sphere_logpdf(50.0, 10.0, 55.0, 20.0, 1000.0, 1.0);
        assert!(near > far);
    }

    #[test]
    fn region_rows_sum_to_zero() {
        let w = vec![
            vec![0.0, 0.5, 0.2],
            vec![0.3, 0.0, 0.4],
            vec![0.1, 0.6, 0.0],
        ];
        let q = region_rate_matrix(&w);
        for row in &q {
            let s: f64 = row.iter().sum();
            assert!(s.abs() < 1e-12);
        }
    }
}
