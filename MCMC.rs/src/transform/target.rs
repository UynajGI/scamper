use super::{Bijector, TransformError};
use crate::target::LogDensity;
use crate::McmcError;

/// Target wrapper that includes a transform Jacobian in unconstrained space.
pub struct TransformedTarget<T, B> {
    target: T,
    bijector: B,
    constrained: Vec<f64>,
}

impl<T, B> TransformedTarget<T, B>
where
    B: Bijector,
{
    pub fn new(target: T, bijector: B) -> Result<Self, McmcError> {
        let unconstrained_dimension = bijector.unconstrained_dimension();
        let constrained_dimension = bijector.constrained_dimension();
        if unconstrained_dimension == 0 || constrained_dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "transform dimensions must be positive".to_string(),
            ));
        }
        Ok(Self {
            target,
            bijector,
            constrained: vec![0.0; constrained_dimension],
        })
    }

    pub const fn target(&self) -> &T {
        &self.target
    }

    pub fn target_mut(&mut self) -> &mut T {
        &mut self.target
    }

    pub const fn bijector(&self) -> &B {
        &self.bijector
    }

    pub fn constrained_position(&self) -> &[f64] {
        &self.constrained
    }

    pub fn into_inner(self) -> (T, B) {
        (self.target, self.bijector)
    }
}

impl<T, B> LogDensity<[f64]> for TransformedTarget<T, B>
where
    T: LogDensity<[f64]>,
    B: Bijector,
{
    fn log_density(&mut self, state: &[f64]) -> f64 {
        if state.len() != self.bijector.unconstrained_dimension() {
            return f64::NAN;
        }
        let log_jacobian = match self.bijector.forward(state, &mut self.constrained) {
            Ok(value) => value,
            Err(TransformError::DimensionMismatch { .. }) => return f64::NAN,
            Err(TransformError::OutsideDomain(_) | TransformError::NonFinite) => {
                return f64::NEG_INFINITY;
            }
        };
        let constrained_log_density = self.target.log_density(&self.constrained);
        if constrained_log_density == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else {
            constrained_log_density + log_jacobian
        }
    }
}
