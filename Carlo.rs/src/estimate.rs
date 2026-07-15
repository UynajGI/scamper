//! Statistical estimates for simulation observables.

use std::fmt;

use serde::{Deserialize, Serialize};

/// Statistical estimate with error analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Estimate {
    /// Mean value.
    pub mean: f64,

    /// Standard error of the mean.
    pub stderr: f64,

    /// Integrated autocorrelation time.
    pub autocorr_time: f64,

    /// Number of bins used.
    pub n_bins: usize,
}

impl Estimate {
    /// Compute estimate from bin means.
    pub fn from_bins(bins: &[f64]) -> Self {
        if bins.is_empty() {
            return Self {
                mean: 0.0,
                stderr: 0.0,
                autocorr_time: 1.0,
                n_bins: 0,
            };
        }

        let n = bins.len() as f64;
        let mean = bins.iter().sum::<f64>() / n;

        // Standard deviation of bin means
        let variance = if n > 1.0 {
            bins.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };

        // Standard error = std / sqrt(n_bins)
        let stderr = variance.sqrt() / n.sqrt();

        Self {
            mean,
            stderr,
            autocorr_time: 1.0, // Placeholder; proper estimation is complex
            n_bins: bins.len(),
        }
    }

    /// Convenience: format as "mean ± stderr".
    pub fn format(&self) -> String {
        format!("{:.6} ± {:.6}", self.mean, self.stderr)
    }
}

/// Complex statistical estimate with error analysis.
///
/// Stores real and imaginary parts as separate estimates.
/// JSON serialization matches Carlo.jl format: `{re, im}`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexEstimate {
    /// Real part estimate.
    pub re: Estimate,
    /// Imaginary part estimate.
    pub im: Estimate,
}

impl ComplexEstimate {
    /// Create from real and imaginary parts.
    pub fn new(re: Estimate, im: Estimate) -> Self {
        Self { re, im }
    }

    /// Create from complex bin values.
    pub fn from_bins(bins_re: &[f64], bins_im: &[f64]) -> Self {
        Self {
            re: Estimate::from_bins(bins_re),
            im: Estimate::from_bins(bins_im),
        }
    }

    /// Format as "(re ± err) + (im ± err)i".
    pub fn format(&self) -> String {
        format!(
            "({:.6} ± {:.6}) + ({:.6} ± {:.6})i",
            self.re.mean, self.re.stderr, self.im.mean, self.im.stderr
        )
    }

    /// Magnitude of the complex mean.
    pub fn magnitude(&self) -> f64 {
        self.re.mean.hypot(self.im.mean)
    }
}

impl fmt::Display for ComplexEstimate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({:.6} ± {:.6}) + ({:.6} ± {:.6})i",
            self.re.mean, self.re.stderr, self.im.mean, self.im.stderr
        )
    }
}
