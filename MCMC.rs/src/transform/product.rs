use serde::{Deserialize, Serialize};

use super::{Bijector, TransformError};

/// Static product of two independent transforms.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Product<A, B> {
    first: A,
    second: B,
}

impl<A, B> Product<A, B> {
    pub const fn new(first: A, second: B) -> Self {
        Self { first, second }
    }

    pub const fn first(&self) -> &A {
        &self.first
    }

    pub const fn second(&self) -> &B {
        &self.second
    }

    pub fn into_inner(self) -> (A, B) {
        (self.first, self.second)
    }
}

impl<A, B> Bijector for Product<A, B>
where
    A: Bijector,
    B: Bijector,
{
    fn unconstrained_dimension(&self) -> usize {
        self.first.unconstrained_dimension() + self.second.unconstrained_dimension()
    }

    fn constrained_dimension(&self) -> usize {
        self.first.constrained_dimension() + self.second.constrained_dimension()
    }

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        let first_unconstrained = self.first.unconstrained_dimension();
        let first_constrained = self.first.constrained_dimension();
        if unconstrained.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: unconstrained.len(),
            });
        }
        if constrained.len() != self.constrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.constrained_dimension(),
                actual: constrained.len(),
            });
        }
        let (left_input, right_input) = unconstrained.split_at(first_unconstrained);
        let (left_output, right_output) = constrained.split_at_mut(first_constrained);
        let left_log_jacobian = self.first.forward(left_input, left_output)?;
        let right_log_jacobian = self.second.forward(right_input, right_output)?;
        Ok(left_log_jacobian + right_log_jacobian)
    }

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError> {
        let first_constrained = self.first.constrained_dimension();
        let first_unconstrained = self.first.unconstrained_dimension();
        if constrained.len() != self.constrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.constrained_dimension(),
                actual: constrained.len(),
            });
        }
        if unconstrained.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: unconstrained.len(),
            });
        }
        let (left_input, right_input) = constrained.split_at(first_constrained);
        let (left_output, right_output) = unconstrained.split_at_mut(first_unconstrained);
        let left_log_jacobian = self.first.inverse(left_input, left_output)?;
        let right_log_jacobian = self.second.inverse(right_input, right_output)?;
        Ok(left_log_jacobian + right_log_jacobian)
    }
}
