use serde::{Deserialize, Serialize};

use crate::McmcError;

/// Online dense covariance estimate using a matrix-valued Welford recurrence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DenseCovarianceAdaptation {
    dimension: usize,
    mean: Vec<f64>,
    m2: Vec<f64>,
    count: u64,
    regularization: f64,
    frozen: bool,
    #[serde(default)]
    delta: Vec<f64>,
}

impl DenseCovarianceAdaptation {
    pub fn new(dimension: usize, regularization: f64) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "dense covariance adaptation dimension must be positive".to_string(),
            ));
        }
        if !regularization.is_finite() || regularization <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "dense covariance regularization must be finite and positive".to_string(),
            ));
        }
        Ok(Self {
            dimension,
            mean: vec![0.0; dimension],
            m2: vec![0.0; dimension.saturating_mul(dimension)],
            count: 0,
            regularization,
            frozen: false,
            delta: vec![0.0; dimension],
        })
    }

    pub fn observe(&mut self, position: &[f64]) -> Result<(), McmcError> {
        if self.frozen {
            return Err(McmcError::AdaptationFrozen);
        }
        if position.len() != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: position.len(),
            });
        }
        if position.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "dense covariance observation must be finite".to_string(),
            ));
        }
        if self.mean.len() != self.dimension
            || self.m2.len() != self.dimension.saturating_mul(self.dimension)
        {
            return Err(McmcError::InvalidConfig(
                "dense covariance checkpoint has inconsistent workspace dimensions".to_string(),
            ));
        }
        if self.delta.len() != self.dimension {
            self.delta.resize(self.dimension, 0.0);
        }

        self.count = self.count.saturating_add(1);
        let count = self.count as f64;
        for (index, value) in position.iter().copied().enumerate() {
            self.delta[index] = value - self.mean[index];
            self.mean[index] += self.delta[index] / count;
        }
        for (row, &delta_row) in self.delta.iter().enumerate() {
            for (column, (&pos_col, &mean_col)) in position.iter().zip(self.mean.iter()).enumerate()
            {
                self.m2[row * self.dimension + column] += delta_row * (pos_col - mean_col);
            }
        }
        Ok(())
    }

    /// Return the regularized covariance without freezing adaptation.
    pub fn covariance(&self) -> Option<Vec<f64>> {
        if self.count < 2
            || self.mean.len() != self.dimension
            || self.m2.len() != self.dimension.saturating_mul(self.dimension)
        {
            return None;
        }
        let denominator = (self.count - 1) as f64;
        let mut covariance = vec![0.0; self.dimension * self.dimension];
        for row in 0..self.dimension {
            for column in 0..self.dimension {
                let left = self.m2[row * self.dimension + column];
                let right = self.m2[column * self.dimension + row];
                covariance[row * self.dimension + column] = 0.5 * (left + right) / denominator;
            }
            covariance[row * self.dimension + row] += self.regularization;
        }
        Some(covariance)
    }

    /// Freeze adaptation and return a row-major lower Cholesky factor.
    pub fn finalize_cholesky(&mut self) -> Option<Vec<f64>> {
        self.frozen = true;
        let covariance = self.covariance()?;
        regularized_cholesky(&covariance, self.dimension, self.regularization)
    }

    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }

    pub const fn count(&self) -> u64 {
        self.count
    }
}

pub(crate) fn regularized_cholesky(
    matrix: &[f64],
    dimension: usize,
    base_jitter: f64,
) -> Option<Vec<f64>> {
    if dimension == 0 || matrix.len() != dimension.saturating_mul(dimension) {
        return None;
    }
    let mut jitter = 0.0;
    for attempt in 0..8 {
        let mut lower = vec![0.0; matrix.len()];
        let mut success = true;
        for row in 0..dimension {
            for column in 0..=row {
                let mut value =
                    0.5 * (matrix[row * dimension + column] + matrix[column * dimension + row]);
                if row == column {
                    value += jitter;
                }
                for inner in 0..column {
                    value -= lower[row * dimension + inner] * lower[column * dimension + inner];
                }
                if row == column {
                    if !value.is_finite() || value <= 0.0 {
                        success = false;
                        break;
                    }
                    lower[row * dimension + column] = value.sqrt();
                } else {
                    let diagonal = lower[column * dimension + column];
                    if !diagonal.is_finite() || diagonal <= 0.0 {
                        success = false;
                        break;
                    }
                    lower[row * dimension + column] = value / diagonal;
                }
            }
            if !success {
                break;
            }
        }
        if success {
            return Some(lower);
        }
        jitter = if attempt == 0 {
            base_jitter
        } else {
            jitter * 10.0
        };
    }
    None
}
