//! Statistical-efficiency estimates for benchmark time series.

/// Integrated autocorrelation and effective-sample summary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StatisticalEfficiency {
    /// Integrated autocorrelation time in units of recorded samples.
    pub integrated_autocorrelation_time: f64,
    /// Effective number of independent samples.
    pub effective_samples: f64,
    /// Effective samples produced per wall-clock second.
    pub effective_samples_per_second: f64,
}

/// Estimate statistical efficiency using Geyer's initial-positive-pair rule.
///
/// The reported convention is `ESS = N / (2 τ_int)`, so independent samples
/// have `τ_int = 0.5` and `ESS = N`.
pub fn statistical_efficiency(samples: &[f64], elapsed_seconds: f64) -> StatisticalEfficiency {
    assert!(
        elapsed_seconds.is_finite() && elapsed_seconds >= 0.0,
        "elapsed time must be finite and non-negative"
    );
    assert!(
        samples.iter().all(|value| value.is_finite()),
        "statistical-efficiency samples must be finite"
    );

    let n = samples.len();
    if n < 2 {
        let effective_samples = n as f64;
        return StatisticalEfficiency {
            integrated_autocorrelation_time: 0.5,
            effective_samples,
            effective_samples_per_second: if elapsed_seconds > 0.0 {
                effective_samples / elapsed_seconds
            } else {
                0.0
            },
        };
    }

    let mean = samples.iter().sum::<f64>() / n as f64;
    let centered = samples.iter().map(|value| value - mean).collect::<Vec<_>>();
    let variance = centered.iter().map(|value| value * value).sum::<f64>() / n as f64;
    if variance <= f64::EPSILON {
        let effective_samples = n as f64;
        return StatisticalEfficiency {
            integrated_autocorrelation_time: 0.5,
            effective_samples,
            effective_samples_per_second: if elapsed_seconds > 0.0 {
                effective_samples / elapsed_seconds
            } else {
                0.0
            },
        };
    }

    let max_lag = (n / 2).max(1);
    let correlation = |lag: usize| {
        centered[..n - lag]
            .iter()
            .zip(&centered[lag..])
            .map(|(left, right)| left * right)
            .sum::<f64>()
            / (n - lag) as f64
            / variance
    };

    let mut tau = 0.5;
    let mut lag = 1usize;
    while lag <= max_lag {
        let first = correlation(lag);
        let second = if lag < max_lag {
            correlation(lag + 1)
        } else {
            0.0
        };
        let pair = first + second;
        if !pair.is_finite() || pair <= 0.0 {
            break;
        }
        tau += pair;
        lag += 2;
    }
    tau = tau.max(0.5);
    let effective_samples = (n as f64 / (2.0 * tau)).clamp(1.0, n as f64);
    StatisticalEfficiency {
        integrated_autocorrelation_time: tau,
        effective_samples,
        effective_samples_per_second: if elapsed_seconds > 0.0 {
            effective_samples / elapsed_seconds
        } else {
            0.0
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_series_is_well_defined() {
        let result = statistical_efficiency(&[3.0; 64], 2.0);
        assert_eq!(result.integrated_autocorrelation_time, 0.5);
        assert_eq!(result.effective_samples, 64.0);
        assert_eq!(result.effective_samples_per_second, 32.0);
    }

    #[test]
    fn positively_correlated_series_has_lower_ess() {
        let samples = (0..256)
            .scan(0.0, |state, index| {
                *state = 0.95 * *state + ((index * 17 % 23) as f64 - 11.0);
                Some(*state)
            })
            .collect::<Vec<_>>();
        let result = statistical_efficiency(&samples, 1.0);
        assert!(result.integrated_autocorrelation_time > 0.5);
        assert!(result.effective_samples < samples.len() as f64);
    }
}
