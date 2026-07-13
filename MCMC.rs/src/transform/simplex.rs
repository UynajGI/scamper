use serde::{Deserialize, Serialize};

use super::{
    check_lengths, log_one_minus_sigmoid, log_sigmoid, sigmoid, Bijector, DifferentiableBijector,
    TransformError,
};
use crate::McmcError;

/// Stick-breaking transform from `R^(K-1)` to the interior of a K-simplex.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct Simplex {
    dimension: usize,
}

impl Simplex {
    pub fn new(dimension: usize) -> Result<Self, McmcError> {
        if dimension < 2 {
            return Err(McmcError::InvalidConfig(
                "simplex transform dimension must be at least two".to_string(),
            ));
        }
        Ok(Self { dimension })
    }
}

impl Bijector for Simplex {
    fn unconstrained_dimension(&self) -> usize {
        self.dimension - 1
    }

    fn constrained_dimension(&self) -> usize {
        self.dimension
    }

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(
            unconstrained,
            self.dimension - 1,
            constrained,
            self.dimension,
        )?;
        let mut remaining = 1.0;
        let mut log_remaining = 0.0;
        let mut log_jacobian = 0.0;
        for (output, &unconstrained_value) in constrained
            .iter_mut()
            .take(self.dimension - 1)
            .zip(unconstrained)
        {
            let fraction = sigmoid(unconstrained_value);
            let value = remaining * fraction;
            if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 || value <= 0.0 {
                return Err(TransformError::NonFinite);
            }
            *output = value;
            log_jacobian += log_remaining
                + log_sigmoid(unconstrained_value)
                + log_one_minus_sigmoid(unconstrained_value);
            remaining *= 1.0 - fraction;
            log_remaining += log_one_minus_sigmoid(unconstrained_value);
        }
        if !remaining.is_finite() || remaining <= 0.0 {
            return Err(TransformError::NonFinite);
        }
        constrained[self.dimension - 1] = remaining;
        Ok(log_jacobian)
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(
            constrained,
            self.dimension,
            unconstrained,
            self.dimension - 1,
        )?;
        if constrained
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(TransformError::OutsideDomain(
                "simplex coordinates must be finite and strictly positive".to_string(),
            ));
        }
        let total = constrained.iter().sum::<f64>();
        if (total - 1.0).abs() > 1.0e-10 {
            return Err(TransformError::OutsideDomain(
                "simplex coordinates must sum to one".to_string(),
            ));
        }

        let mut remaining = 1.0;
        let mut forward_log_jacobian = 0.0;
        let mut log_remaining = 0.0;
        for (&coordinate, output) in constrained
            .iter()
            .take(self.dimension - 1)
            .zip(unconstrained.iter_mut())
        {
            let fraction = coordinate / remaining;
            if !fraction.is_finite() || fraction <= 0.0 || fraction >= 1.0 {
                return Err(TransformError::OutsideDomain(
                    "simplex stick fraction lies outside (0, 1)".to_string(),
                ));
            }
            *output = fraction.ln() - (-fraction).ln_1p();
            forward_log_jacobian += log_remaining + fraction.ln() + (-fraction).ln_1p();
            remaining -= coordinate;
            log_remaining = remaining.ln();
        }
        Ok(-forward_log_jacobian)
    }
}

impl DifferentiableBijector for Simplex {
    fn pullback(
        &mut self,
        unconstrained: &[f64],
        constrained: &[f64],
        constrained_gradient: &[f64],
        unconstrained_gradient: &mut [f64],
    ) -> Result<(), TransformError> {
        check_lengths(
            unconstrained,
            self.dimension - 1,
            constrained,
            self.dimension,
        )?;
        check_lengths(
            constrained_gradient,
            self.dimension,
            unconstrained_gradient,
            self.dimension - 1,
        )?;
        let mut weighted_tail =
            constrained_gradient[self.dimension - 1] * constrained[self.dimension - 1];
        for index in (0..self.dimension - 1).rev() {
            let fraction = sigmoid(unconstrained[index]);
            let own = constrained_gradient[index] * constrained[index] * (1.0 - fraction);
            unconstrained_gradient[index] = own - fraction * weighted_tail;
            weighted_tail += constrained_gradient[index] * constrained[index];
        }
        if unconstrained_gradient.iter().all(|value| value.is_finite()) {
            Ok(())
        } else {
            Err(TransformError::NonFinite)
        }
    }

    fn log_jacobian_gradient(
        &mut self,
        unconstrained: &[f64],
        output: &mut [f64],
    ) -> Result<(), TransformError> {
        check_lengths(
            unconstrained,
            self.dimension - 1,
            output,
            self.dimension - 1,
        )?;
        for (index, (value, output)) in unconstrained
            .iter()
            .copied()
            .zip(output.iter_mut())
            .enumerate()
        {
            *output = 1.0 - (self.dimension - index) as f64 * sigmoid(value);
        }
        Ok(())
    }
}
