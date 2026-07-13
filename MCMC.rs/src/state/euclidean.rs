use serde::{Deserialize, Serialize};

use crate::target::{validate_log_density, LogDensity};
use crate::{ChainState, McmcError};

/// Cached gradient for gradient-based Euclidean kernels.
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

    pub(crate) fn set_gradient(&mut self, gradient: &[f64]) -> Result<(), McmcError> {
        if gradient.len() != self.gradient.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.gradient.len(),
                actual: gradient.len(),
            });
        }
        if gradient.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "gradient cache must contain only finite values".to_string(),
            ));
        }
        self.gradient.copy_from_slice(gradient);
        self.gradient_valid = true;
        Ok(())
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

    pub(crate) fn synchronize_gradient(
        &mut self,
        log_density: f64,
        gradient: &[f64],
    ) -> Result<(), McmcError> {
        let log_density = validate_log_density(log_density)?;
        if log_density == f64::NEG_INFINITY {
            return Err(McmcError::InvalidConfig(
                "accepted state cannot lie outside target support".to_string(),
            ));
        }
        self.cache_mut().set_gradient(gradient)?;
        self.synchronize_log_density(log_density);
        Ok(())
    }

    pub(crate) fn commit_hamiltonian_proposal(
        &mut self,
        proposed_position: &mut Vec<f64>,
        proposed_log_density: f64,
        proposed_gradient: &[f64],
    ) -> Result<(), McmcError> {
        let proposed_log_density = validate_log_density(proposed_log_density)?;
        if proposed_log_density == f64::NEG_INFINITY {
            return Err(McmcError::InvalidConfig(
                "Hamiltonian proposal lies outside target support".to_string(),
            ));
        }
        if proposed_position.len() != self.dimension() {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension(),
                actual: proposed_position.len(),
            });
        }
        if proposed_position.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "Hamiltonian proposal position must be finite".to_string(),
            ));
        }
        self.cache_mut().set_gradient(proposed_gradient)?;
        self.swap_position(proposed_position, proposed_log_density);
        Ok(())
    }

    /// Validate invariants required by every Euclidean transition kernel.
    pub fn validate(&self) -> Result<(), McmcError> {
        let dimension = self.dimension();
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "state dimension must be positive".to_string(),
            ));
        }
        if self.position().iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "state position must contain only finite values".to_string(),
            ));
        }
        let log_density = validate_log_density(self.log_density())?;
        if log_density == f64::NEG_INFINITY {
            return Err(McmcError::InvalidConfig(
                "accepted state cannot lie outside target support".to_string(),
            ));
        }
        if self.cache().gradient.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: self.cache().gradient.len(),
            });
        }
        if self.cache().gradient_valid
            && self.cache().gradient.iter().any(|value| !value.is_finite())
        {
            return Err(McmcError::InvalidConfig(
                "valid gradient cache must contain only finite values".to_string(),
            ));
        }
        Ok(())
    }
}
