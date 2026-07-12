use serde::{Deserialize, Serialize};

use super::{check_lengths, Bijector, TransformError};
use crate::McmcError;

/// Ordered-vector transform with an unconstrained first element and positive gaps.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ordered {
    dimension: usize,
}

impl Ordered {
    pub fn new(dimension: usize) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "ordered transform dimension must be positive".to_string(),
            ));
        }
        Ok(Self { dimension })
    }
}

impl Bijector for Ordered {
    fn unconstrained_dimension(&self) -> usize {
        self.dimension
    }

    fn constrained_dimension(&self) -> usize {
        self.dimension
    }

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(unconstrained, self.dimension, constrained, self.dimension)?;
        constrained[0] = unconstrained[0];
        let mut log_jacobian = 0.0;
        for index in 1..self.dimension {
            let unconstrained_value = unconstrained[index];
            let gap = unconstrained_value.exp();
            let previous = constrained[index - 1];
            let value = previous + gap;
            if !gap.is_finite() || !value.is_finite() || value <= previous {
                return Err(TransformError::NonFinite);
            }
            constrained[index] = value;
            log_jacobian += unconstrained_value;
        }
        Ok(log_jacobian)
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(constrained, self.dimension, unconstrained, self.dimension)?;
        if constrained.iter().any(|value| !value.is_finite()) {
            return Err(TransformError::NonFinite);
        }
        unconstrained[0] = constrained[0];
        let mut inverse_log_jacobian = 0.0;
        for (output, pair) in unconstrained.iter_mut().skip(1).zip(constrained.windows(2)) {
            let gap = pair[1] - pair[0];
            if gap <= 0.0 || !gap.is_finite() {
                return Err(TransformError::OutsideDomain(
                    "ordered vector must be strictly increasing".to_string(),
                ));
            }
            *output = gap.ln();
            inverse_log_jacobian -= *output;
        }
        Ok(inverse_log_jacobian)
    }
}
