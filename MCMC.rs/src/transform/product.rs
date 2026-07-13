use serde::{Deserialize, Serialize};

use super::{Bijector, DifferentiableBijector, TransformError};

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

impl<A, B> DifferentiableBijector for Product<A, B>
where
    A: DifferentiableBijector,
    B: DifferentiableBijector,
{
    fn pullback(
        &mut self,
        unconstrained: &[f64],
        constrained: &[f64],
        constrained_gradient: &[f64],
        unconstrained_gradient: &mut [f64],
    ) -> Result<(), TransformError> {
        let first_unconstrained = self.first.unconstrained_dimension();
        let first_constrained = self.first.constrained_dimension();
        if unconstrained.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: unconstrained.len(),
            });
        }
        if constrained.len() != self.constrained_dimension()
            || constrained_gradient.len() != self.constrained_dimension()
        {
            return Err(TransformError::DimensionMismatch {
                expected: self.constrained_dimension(),
                actual: constrained.len().min(constrained_gradient.len()),
            });
        }
        if unconstrained_gradient.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: unconstrained_gradient.len(),
            });
        }
        let (left_z, right_z) = unconstrained.split_at(first_unconstrained);
        let (left_x, right_x) = constrained.split_at(first_constrained);
        let (left_gx, right_gx) = constrained_gradient.split_at(first_constrained);
        let (left_gz, right_gz) = unconstrained_gradient.split_at_mut(first_unconstrained);
        self.first.pullback(left_z, left_x, left_gx, left_gz)?;
        self.second.pullback(right_z, right_x, right_gx, right_gz)
    }

    fn log_jacobian_gradient(
        &mut self,
        unconstrained: &[f64],
        output: &mut [f64],
    ) -> Result<(), TransformError> {
        let first_unconstrained = self.first.unconstrained_dimension();
        if unconstrained.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: unconstrained.len(),
            });
        }
        if output.len() != self.unconstrained_dimension() {
            return Err(TransformError::DimensionMismatch {
                expected: self.unconstrained_dimension(),
                actual: output.len(),
            });
        }
        let (left_input, right_input) = unconstrained.split_at(first_unconstrained);
        let (left_output, right_output) = output.split_at_mut(first_unconstrained);
        self.first.log_jacobian_gradient(left_input, left_output)?;
        self.second.log_jacobian_gradient(right_input, right_output)
    }
}
