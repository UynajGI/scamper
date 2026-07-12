use serde::{Deserialize, Serialize};

use super::{check_lengths, log_one_minus_sigmoid, log_sigmoid, sigmoid, Bijector, TransformError};
use crate::McmcError;

/// Identity transform over a fixed positive dimension.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct Identity {
    dimension: usize,
}

impl Identity {
    pub fn new(dimension: usize) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "identity transform dimension must be positive".to_string(),
            ));
        }
        Ok(Self { dimension })
    }
}

impl Bijector for Identity {
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
        constrained.copy_from_slice(unconstrained);
        Ok(0.0)
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(constrained, self.dimension, unconstrained, self.dimension)?;
        unconstrained.copy_from_slice(constrained);
        Ok(0.0)
    }
}

/// Scalar positive transform `x = exp(z)`.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Positive;

impl Bijector for Positive {
    fn unconstrained_dimension(&self) -> usize {
        1
    }

    fn constrained_dimension(&self) -> usize {
        1
    }

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(unconstrained, 1, constrained, 1)?;
        let value = unconstrained[0].exp();
        if !value.is_finite() || value <= 0.0 {
            return Err(TransformError::NonFinite);
        }
        constrained[0] = value;
        Ok(unconstrained[0])
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(constrained, 1, unconstrained, 1)?;
        if !constrained[0].is_finite() || constrained[0] <= 0.0 {
            return Err(TransformError::OutsideDomain(
                "positive transform requires x > 0".to_string(),
            ));
        }
        unconstrained[0] = constrained[0].ln();
        Ok(-unconstrained[0])
    }
}

/// Scalar bounded transform using a numerically stable logistic map.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Interval {
    lower: f64,
    upper: f64,
}

impl Interval {
    pub fn new(lower: f64, upper: f64) -> Result<Self, McmcError> {
        if !lower.is_finite() || !upper.is_finite() || lower >= upper {
            return Err(McmcError::InvalidConfig(
                "interval transform requires finite lower < upper".to_string(),
            ));
        }
        Ok(Self { lower, upper })
    }

    pub const fn lower(&self) -> f64 {
        self.lower
    }

    pub const fn upper(&self) -> f64 {
        self.upper
    }
}

impl Bijector for Interval {
    fn unconstrained_dimension(&self) -> usize {
        1
    }

    fn constrained_dimension(&self) -> usize {
        1
    }

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(unconstrained, 1, constrained, 1)?;
        let probability = sigmoid(unconstrained[0]);
        let width = self.upper - self.lower;
        let value = self.lower + width * probability;
        if !value.is_finite() || value <= self.lower || value >= self.upper {
            return Err(TransformError::NonFinite);
        }
        constrained[0] = value;
        Ok(width.ln() + log_sigmoid(unconstrained[0]) + log_one_minus_sigmoid(unconstrained[0]))
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        check_lengths(constrained, 1, unconstrained, 1)?;
        let value = constrained[0];
        if !value.is_finite() || value <= self.lower || value >= self.upper {
            return Err(TransformError::OutsideDomain(
                "interval value must lie strictly inside its bounds".to_string(),
            ));
        }
        let probability = (value - self.lower) / (self.upper - self.lower);
        let unconstrained_value = probability.ln() - (-probability).ln_1p();
        unconstrained[0] = unconstrained_value;
        Ok(-((self.upper - self.lower).ln()
            + log_sigmoid(unconstrained_value)
            + log_one_minus_sigmoid(unconstrained_value)))
    }
}
