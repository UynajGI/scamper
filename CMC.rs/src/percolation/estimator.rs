//! Shared scalar sample-moment accumulator for report-time estimators.
//!
//! The first-passage, passage-shape, and chemical-distance estimators all
//! accumulate independent realizations the same way: a count, a sum, and a
//! sum of squares, reporting the unbiased sample variance
//! `(sum_sq - sum^2/n) / (n - 1)` clamped at zero and the standard error
//! `sqrt(variance / n)`. This type owns that arithmetic once; the
//! estimators keep their own validation, sample spaces, and report structs
//! because those legitimately differ (a profile entry may be conditioned on
//! reachability, a diameter is always finite, and so on).
//!
//! All methods are allocation-free and branch-free after `add`'s identity
//! check; the order of floating-point operations is fixed to match the
//! expressions the estimators used before this type existed.

/// Running count, sum, and sum of squares of independent samples.
///
/// The estimator order is preserved exactly: `add` accumulates
/// `value` then `value * value`, and the variance expression is evaluated
/// as written in [`Self::sample_variance`], so existing bitwise results do
/// not change.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Moments {
    samples: usize,
    sum: f64,
    sum_sq: f64,
}

impl Moments {
    /// An empty accumulator.
    pub const fn new() -> Self {
        Self {
            samples: 0,
            sum: 0.0,
            sum_sq: 0.0,
        }
    }

    /// Accumulate one finite sample.
    #[inline]
    pub fn add(&mut self, value: f64) {
        self.samples += 1;
        self.sum += value;
        self.sum_sq += value * value;
    }

    /// Samples accumulated.
    #[inline]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// Sample mean, or `None` with no samples.
    #[inline]
    pub fn mean(&self) -> Option<f64> {
        (self.samples != 0).then(|| self.sum / self.samples as f64)
    }

    /// Unbiased sample variance, or `None` with fewer than two samples.
    ///
    /// The subtraction form is retained from the v2 estimators (and clamped
    /// at zero) so numeric results stay bit-identical; it is not a
    /// numerically stable Welford update.
    #[inline]
    pub fn sample_variance(&self) -> Option<f64> {
        if self.samples < 2 {
            return None;
        }
        let samples = self.samples as f64;
        Some(((self.sum_sq - self.sum * self.sum / samples) / (samples - 1.0)).max(0.0))
    }

    /// Standard error of [`Self::mean`], or `None` with fewer than two
    /// samples.
    #[inline]
    pub fn standard_error(&self) -> Option<f64> {
        self.sample_variance()
            .map(|variance| (variance / self.samples as f64).sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_and_single_sample_have_no_variance_or_error() {
        let empty = Moments::new();
        assert_eq!(empty.samples(), 0);
        assert_eq!(empty.mean(), None);
        assert_eq!(empty.sample_variance(), None);
        assert_eq!(empty.standard_error(), None);

        let mut single = Moments::new();
        single.add(2.5);
        assert_eq!(single.samples(), 1);
        assert_eq!(single.mean(), Some(2.5));
        assert_eq!(single.sample_variance(), None);
        assert_eq!(single.standard_error(), None);
    }

    #[test]
    fn known_samples_report_the_unbiased_variance_and_error() {
        let mut moments = Moments::new();
        for value in [1.0, 2.0, 3.0, 4.0, 5.0] {
            moments.add(value);
        }
        assert_eq!(moments.samples(), 5);
        assert_eq!(moments.mean(), Some(3.0));
        assert_eq!(moments.sample_variance(), Some(2.5));
        assert_eq!(moments.standard_error(), Some((2.5_f64 / 5.0).sqrt()));
    }

    #[test]
    fn identical_samples_clamp_to_zero_variance() {
        let mut moments = Moments::new();
        moments.add(7.0);
        moments.add(7.0);
        assert_eq!(moments.sample_variance(), Some(0.0));
        assert_eq!(moments.standard_error(), Some(0.0));
    }
}
