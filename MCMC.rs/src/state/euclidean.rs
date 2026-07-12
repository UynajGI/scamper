use serde::{Deserialize, Serialize};

use crate::target::{validate_log_density, LogDensity};
use crate::{ChainState, McmcError};

/// Cache reserved for future gradient-based Euclidean kernels.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EuclideanCache {
    gradient: Vec<f64>,
    gradient_valid: bool,
}

impl EuclideanCache {
    pub fn with_dimension(dimension: usize) -> Self {
        Self {
            gradient: vec![0.0; dimension],
            gradient_valid: false,
        }
    }

    pub fn gradient(&self) -> Option<&[f64]> {
        self.gradient_valid.then_some(self.gradient.as_slice())
    }

    pub fn gradient_buffer_mut(&mut self) -> &mut [f64] {
        self.gradient_valid = false;
        &mut self.gradient
    }

    pub fn mark_gradient_valid(&mut self) {
        self.gradient_valid = true;
    }

    pub fn invalidate_gradient(&mut self) {
        self.gradient_valid = false;
    }
}

pub type EuclideanState = ChainState<Vec<f64>, EuclideanCache>;

impl EuclideanState {
    pub fn initialize<T>(target: &mut T, position: Vec<f64>) -> Result<Self, McmcError>
    where
        T: LogDensity<[f64]>,
    {
        if position.is_empty() {
            return Err(McmcError::InvalidConfig(
                "state dimension must be positive".to_string(),
            ));
        }
        if position.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "initial position must contain only finite values".to_string(),
            ));
        }
        let log_density = validate_log_density(target.log_density(&position))?;
        if log_density == f64::NEG_INFINITY {
            return Err(McmcError::InvalidConfig(
                "initial position lies outside target support".to_string(),
            ));
        }
        let dimension = position.len();
        Ok(Self::new(
            position,
            log_density,
            EuclideanCache::with_dimension(dimension),
        ))
    }

    pub fn dimension(&self) -> usize {
        self.position().len()
    }
}
