use super::{Bijector, DifferentiableBijector, TransformError};
use crate::target::{DifferentiableLogDensity, LogDensity};
use crate::McmcError;

/// Target wrapper that includes a transform Jacobian in unconstrained space.
pub struct TransformedTarget<T, B> {
    target: T,
    bijector: B,
    constrained: Vec<f64>,
    constrained_gradient: Vec<f64>,
    jacobian_gradient: Vec<f64>,
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
            constrained_gradient: vec![0.0; constrained_dimension],
            jacobian_gradient: vec![0.0; unconstrained_dimension],
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

impl<T, B> DifferentiableLogDensity for TransformedTarget<T, B>
where
    T: DifferentiableLogDensity,
    B: DifferentiableBijector,
{
    fn log_density_and_gradient(&mut self, state: &[f64], gradient: &mut [f64]) -> f64 {
        let unconstrained_dimension = self.bijector.unconstrained_dimension();
        if state.len() != unconstrained_dimension || gradient.len() != unconstrained_dimension {
            gradient.fill(f64::NAN);
            return f64::NAN;
        }
        let log_jacobian = match self.bijector.forward(state, &mut self.constrained) {
            Ok(value) => value,
            Err(TransformError::DimensionMismatch { .. }) => {
                gradient.fill(f64::NAN);
                return f64::NAN;
            }
            Err(TransformError::OutsideDomain(_) | TransformError::NonFinite) => {
                gradient.fill(f64::NAN);
                return f64::NEG_INFINITY;
            }
        };
        let constrained_log_density = self
            .target
            .log_density_and_gradient(&self.constrained, &mut self.constrained_gradient);
        if !constrained_log_density.is_finite()
            || self
                .constrained_gradient
                .iter()
                .any(|value| !value.is_finite())
        {
            gradient.fill(f64::NAN);
            return constrained_log_density;
        }
        if self
            .bijector
            .pullback(
                state,
                &self.constrained,
                &self.constrained_gradient,
                gradient,
            )
            .is_err()
            || self
                .bijector
                .log_jacobian_gradient(state, &mut self.jacobian_gradient)
                .is_err()
        {
            gradient.fill(f64::NAN);
            return f64::NAN;
        }
        for (gradient, jacobian) in gradient
            .iter_mut()
            .zip(self.jacobian_gradient.iter().copied())
        {
            *gradient += jacobian;
        }
        constrained_log_density + log_jacobian
    }
}
