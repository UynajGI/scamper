//! Numerically stable stochastic acceptance decisions.

use rand::{Rng, RngExt};

/// Draw a Metropolis-Hastings decision from a log acceptance ratio.
///
/// This performs the comparison entirely in log space:
/// `ln(u) < min(0, log_acceptance)`. Positive infinity is accepted,
/// negative infinity is rejected, and NaN is rejected.
#[inline]
pub fn accept_log_probability<R: Rng + ?Sized>(log_acceptance: f64, rng: &mut R) -> bool {
    if log_acceptance.is_nan() {
        return false;
    }
    if log_acceptance >= 0.0 {
        return true;
    }
    if log_acceptance == f64::NEG_INFINITY {
        return false;
    }
    rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn handles_non_finite_ratios_without_exponentiation() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        assert!(accept_log_probability(f64::INFINITY, &mut rng));
        assert!(!accept_log_probability(f64::NEG_INFINITY, &mut rng));
        assert!(!accept_log_probability(f64::NAN, &mut rng));
    }

    #[test]
    fn zero_log_ratio_is_always_accepted() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        for _ in 0..100 {
            assert!(accept_log_probability(0.0, &mut rng));
        }
    }
}
