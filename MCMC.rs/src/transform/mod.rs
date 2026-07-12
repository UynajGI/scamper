mod product;
mod scalar;
mod simplex;
mod target;
mod vector;

pub use product::Product;
pub use scalar::{Identity, Interval, Positive};
pub use simplex::Simplex;
pub use target::TransformedTarget;
pub use vector::Ordered;

use thiserror::Error;

/// Errors produced while mapping between constrained and unconstrained spaces.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum TransformError {
    #[error("transform dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("constrained value is outside the transform domain: {0}")]
    OutsideDomain(String),

    #[error("transform produced a non-finite numerical result")]
    NonFinite,
}

/// A differentiable bijection used to run MCMC in an unconstrained Euclidean space.
///
/// `forward` returns `log |det(dx / dz)|`; `inverse` returns
/// `log |det(dz / dx)|`.
pub trait Bijector: Send {
    fn unconstrained_dimension(&self) -> usize;
    fn constrained_dimension(&self) -> usize;

    fn forward(
        &mut self,
        unconstrained: &[f64],
        constrained: &mut [f64],
    ) -> Result<f64, TransformError>;

    fn inverse(
        &mut self,
        constrained: &[f64],
        unconstrained: &mut [f64],
    ) -> Result<f64, TransformError>;
}

pub(crate) fn check_lengths(
    input: &[f64],
    expected_input: usize,
    output: &[f64],
    expected_output: usize,
) -> Result<(), TransformError> {
    if input.len() != expected_input {
        return Err(TransformError::DimensionMismatch {
            expected: expected_input,
            actual: input.len(),
        });
    }
    if output.len() != expected_output {
        return Err(TransformError::DimensionMismatch {
            expected: expected_output,
            actual: output.len(),
        });
    }
    Ok(())
}

pub(crate) fn sigmoid(value: f64) -> f64 {
    if value >= 0.0 {
        1.0 / (1.0 + (-value).exp())
    } else {
        let exponential = value.exp();
        exponential / (1.0 + exponential)
    }
}

pub(crate) fn log_sigmoid(value: f64) -> f64 {
    -softplus(-value)
}

pub(crate) fn log_one_minus_sigmoid(value: f64) -> f64 {
    -softplus(value)
}

fn softplus(value: f64) -> f64 {
    if value > 0.0 {
        value + (-value).exp().ln_1p()
    } else {
        value.exp().ln_1p()
    }
}
