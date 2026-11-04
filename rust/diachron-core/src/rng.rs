//! Deterministic, serializable pseudo-random number generator.
//!
//! MCMC checkpoint/resume requires that the RNG state be captured and restored
//! exactly, so chains continue bit-for-bit after an interruption (gate G12).
//! Rust's standard library has no stable, serializable RNG, so Diachron carries
//! its own: a 64-bit `splitmix64` seeder feeding a `xoshiro256**` generator.
//! Both are small, well-studied integer-arithmetic algorithms; the
//! implementation and the serialization are original here.
//!
//! The state is four `u64` words. `snapshot` returns them; `from_state`
//! reconstructs an identical generator. Two generators with the same state
//! produce the same stream.

/// xoshiro256** generator with explicit, serializable state.
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

#[inline]
fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

#[inline]
fn rotl(x: u64, k: u32) -> u64 {
    (x << k) | (x >> (64 - k))
}

impl Rng {
    /// Seed from a single 64-bit value, expanded through splitmix64.
    pub fn seed_from_u64(seed: u64) -> Self {
        let mut x = seed;
        let s = [
            splitmix64(&mut x),
            splitmix64(&mut x),
            splitmix64(&mut x),
            splitmix64(&mut x),
        ];
        Rng { s }
    }

    /// Reconstruct from a saved state.
    pub fn from_state(state: [u64; 4]) -> Self {
        Rng { s: state }
    }

    /// Capture the current state for checkpointing.
    pub fn snapshot(&self) -> [u64; 4] {
        self.s
    }

    /// Next raw 64-bit value (xoshiro256**).
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = rotl(self.s[1].wrapping_mul(5), 7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = rotl(self.s[3], 45);
        result
    }

    /// Uniform in [0, 1). Uses the top 53 bits for a double.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform integer in [0, n).
    #[inline]
    pub fn next_below(&mut self, n: usize) -> usize {
        // Lemire-style rejection to avoid modulo bias; n is tiny in practice.
        assert!(n > 0);
        loop {
            let x = self.next_u64();
            let m = (x as u128) * (n as u128);
            let low = m as u64;
            if low >= (n as u64) {
                return (m >> 64) as usize;
            }
            // Rejection zone; retry.
            let threshold = (0u64.wrapping_sub(n as u64)) % (n as u64);
            if low >= threshold {
                return (m >> 64) as usize;
            }
        }
    }

    /// Standard normal via Box-Muller. Caches the second variate.
    pub fn next_normal(&mut self) -> f64 {
        // No caching across calls keeps the state simple to serialize; the cost
        // of one extra transcendental per draw is negligible for MCMC here.
        let u1 = (self.next_f64()).max(1e-300);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// Exponential with the given rate (mean 1/rate).
    pub fn next_exp(&mut self, rate: f64) -> f64 {
        let u = (self.next_f64()).max(1e-300);
        -u.ln() / rate
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_stream() {
        let mut a = Rng::seed_from_u64(42);
        let mut b = Rng::seed_from_u64(42);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn snapshot_restores_stream() {
        let mut a = Rng::seed_from_u64(7);
        for _ in 0..50 {
            a.next_u64();
        }
        let snap = a.snapshot();
        let future: Vec<u64> = (0..20).map(|_| a.next_u64()).collect();
        let mut b = Rng::from_state(snap);
        let replay: Vec<u64> = (0..20).map(|_| b.next_u64()).collect();
        assert_eq!(future, replay);
    }

    #[test]
    fn uniform_in_range() {
        let mut r = Rng::seed_from_u64(99);
        for _ in 0..10000 {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x));
            let k = r.next_below(13);
            assert!(k < 13);
        }
    }

    #[test]
    fn normal_mean_and_variance() {
        let mut r = Rng::seed_from_u64(123);
        let n = 200_000;
        let mut sum = 0.0;
        let mut sumsq = 0.0;
        for _ in 0..n {
            let z = r.next_normal();
            sum += z;
            sumsq += z * z;
        }
        let mean = sum / n as f64;
        let var = sumsq / n as f64 - mean * mean;
        assert!(mean.abs() < 0.02, "mean {mean}");
        assert!((var - 1.0).abs() < 0.03, "var {var}");
    }
}
