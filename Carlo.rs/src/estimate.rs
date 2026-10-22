//! Statistical estimates for simulation observables.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::CarloError;

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
        Self::from_bins_with_autocorr(bins, 1.0)
    }

    /// Compute estimate from bin means with a caller-provided autocorrelation time.
    pub fn from_bins_with_autocorr(bins: &[f64], autocorr_time: f64) -> Self {
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

        let variance = if n > 1.0 {
            bins.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };

        let stderr = variance.sqrt() / n.sqrt();

        Self {
            mean,
            stderr,
            autocorr_time: autocorr_time.max(1.0),
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

/// Element-preserving estimate for an array observable.
///
/// Each element carries an independent [`Estimate`] (mean, stderr,
/// autocorrelation time, bin count) computed from that element's own bin
/// column. The `shape` records the logical geometry so a merged or
/// reweighted curve keeps its coordinate structure instead of collapsing
/// to one scalar.
///
/// A scalar observable is represented by [`Estimate`] alone; `ArrayEstimate`
/// with an empty `shape` is rejected.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArrayEstimate {
    /// Logical shape of the observable (non-empty; row-major element order).
    pub shape: Vec<usize>,

    /// One estimate per element, `shape` product entries in row-major order.
    pub elements: Vec<Estimate>,
}

impl ArrayEstimate {
    /// Create from a shape and per-element estimates.
    ///
    /// Returns [`CarloError::ShapeMismatch`] unless `elements.len()` equals
    /// the product of `shape` dimensions. Scalar observables (empty `shape`)
    /// are rejected with [`CarloError::InvalidConfig`] — they belong in an
    /// [`Estimate`].
    pub fn new(shape: Vec<usize>, elements: Vec<Estimate>) -> Result<Self, CarloError> {
        if shape.is_empty() || shape.contains(&0) {
            return Err(CarloError::InvalidConfig {
                field: "ArrayEstimate::new".into(),
                reason: "array estimates require a non-empty shape with positive dimensions; \
                         scalar observables are plain Estimate values"
                    .into(),
            });
        }
        let n = shape.iter().product::<usize>();
        if elements.len() != n {
            return Err(CarloError::ShapeMismatch {
                observable: "ArrayEstimate".into(),
                expected_len: n,
                got_len: elements.len(),
            });
        }
        Ok(Self { shape, elements })
    }

    /// Create from one bin-mean column per element.
    ///
    /// Column `i` holds the bin means of element `i`; every column must have
    /// the same length. Autocorrelation times are estimated per element via
    /// the same variance-ratio method the scalar path uses.
    pub fn from_element_bins(shape: &[usize], columns: &[Vec<f64>]) -> Result<Self, CarloError> {
        let n = shape.iter().product::<usize>();
        if shape.is_empty() || n == 0 {
            return Err(CarloError::InvalidConfig {
                field: "ArrayEstimate::from_element_bins".into(),
                reason: "array estimates require a non-empty shape with positive dimensions".into(),
            });
        }
        if columns.len() != n {
            return Err(CarloError::ShapeMismatch {
                observable: "ArrayEstimate".into(),
                expected_len: n,
                got_len: columns.len(),
            });
        }
        let expected_bins = columns.first().map(Vec::len).unwrap_or(0);
        if columns.iter().any(|c| c.len() != expected_bins) {
            return Err(CarloError::InvalidConfig {
                field: "ArrayEstimate::from_element_bins".into(),
                reason: format!(
                    "element columns must share one bin count, expected {expected_bins}"
                ),
            });
        }

        let elements = columns
            .iter()
            .map(|column| {
                let tau = variance_ratio_autocorr(column);
                Estimate::from_bins_with_autocorr(column, tau)
            })
            .collect();
        Ok(Self {
            shape: shape.to_vec(),
            elements,
        })
    }

    /// Logical shape of the observable.
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// Per-element estimates in row-major order.
    pub fn elements(&self) -> &[Estimate] {
        &self.elements
    }

    /// Number of elements (`shape` product).
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Estimate for one element by flat row-major index.
    pub fn element(&self, index: usize) -> Option<&Estimate> {
        self.elements.get(index)
    }

    /// Estimate for one element by multi-dimensional index.
    pub fn at(&self, index: &[usize]) -> Option<&Estimate> {
        if index.len() != self.shape.len() {
            return None;
        }
        let mut flat = 0usize;
        for (&idx, &dim) in index.iter().zip(self.shape.iter()) {
            if idx >= dim {
                return None;
            }
            flat = flat * dim + idx;
        }
        self.elements.get(flat)
    }
}

/// Variance-ratio integrated autocorrelation time over bin means.
///
/// τ = 0.5 · (σ²_rebin / (σ²_bins / 2) − 1), floored at zero and reported
/// as 1.0 when there are too few bins or no variance to compare. This is
/// the same estimator [`crate::Accumulator::autocorr_time`] applies to
/// scalar bin means, factored out so array elements reuse it directly.
pub(crate) fn variance_ratio_autocorr(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return 1.0;
    }

    let n_f64 = n as f64;
    let mu = values.iter().sum::<f64>() / n_f64;
    let var_no_rebin = values.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (n_f64 - 1.0);

    let rebin_n = n / 2;
    if rebin_n < 2 {
        return 1.0;
    }

    let rebin_means: Vec<f64> = (0..rebin_n)
        .map(|i| (values[2 * i] + values[2 * i + 1]) / 2.0)
        .collect();
    let rebin_mu = rebin_means.iter().sum::<f64>() / rebin_n as f64;
    let var_rebin = rebin_means
        .iter()
        .map(|x| (x - rebin_mu).powi(2))
        .sum::<f64>()
        / (rebin_n as f64 - 1.0);

    let expected_var = var_no_rebin / 2.0;
    if expected_var <= 0.0 {
        return 1.0;
    }

    0.5 * (var_rebin / expected_var - 1.0).max(0.0)
}
