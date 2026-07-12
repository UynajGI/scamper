use serde::{Deserialize, Serialize};

use crate::McmcError;

/// Online diagonal covariance estimate using Welford's stable recurrence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagonalCovarianceAdaptation {
    mean: Vec<f64>,
    m2: Vec<f64>,
    count: u64,
    minimum_scale: f64,
    frozen: bool,
}

impl DiagonalCovarianceAdaptation {
    pub fn new(dimension: usize, minimum_scale: f64) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "covariance adaptation dimension must be positive".to_string(),
            ));
        }
        if !minimum_scale.is_finite() || minimum_scale <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "minimum covariance scale must be finite and positive".to_string(),
            ));
        }
        Ok(Self {
            mean: vec![0.0; dimension],
            m2: vec![0.0; dimension],
            count: 0,
            minimum_scale,
            frozen: false,
        })
    }

    pub fn observe(&mut self, position: &[f64]) -> Result<(), McmcError> {
        if self.frozen {
            return Err(McmcError::AdaptationFrozen);
        }
        if position.len() != self.mean.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.mean.len(),
                actual: position.len(),
            });
        }
        self.count = self.count.saturating_add(1);
        let count = self.count as f64;
        for ((mean, m2), value) in self
            .mean
            .iter_mut()
            .zip(self.m2.iter_mut())
            .zip(position.iter().copied())
        {
            let delta = value - *mean;
            *mean += delta / count;
            let delta2 = value - *mean;
            *m2 += delta * delta2;
        }
        Ok(())
    }

    pub fn finalize_scales(&mut self) -> Option<Vec<f64>> {
        self.frozen = true;
        if self.count < 2 {
            return None;
        }
        let denominator = (self.count - 1) as f64;
        Some(
            self.m2
                .iter()
                .map(|m2| (m2 / denominator).sqrt().max(self.minimum_scale))
                .collect(),
        )
    }

    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }
}
