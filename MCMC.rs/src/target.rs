/// Log target density, up to an additive normalization constant.
///
/// `NEG_INFINITY` denotes a point outside the support. `NaN` and positive
/// infinity are invalid. Implementations may mutate internal workspaces, but
/// repeated evaluation at the same state must describe the same density.
pub trait LogDensity<S: ?Sized>: Send {
    fn log_density(&mut self, state: &S) -> f64;
}

/// Log target density with a reusable, allocation-free gradient interface.
pub trait DifferentiableLogDensity: LogDensity<[f64]> {
    /// Write `gradient = ∇ log π(position)` and return `log π(position)`.
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64;
}

/// Adapter from a mutable closure to [`LogDensity`].
pub struct FnLogDensity<F> {
    function: F,
}

impl<F> FnLogDensity<F> {
    pub const fn new(function: F) -> Self {
        Self { function }
    }

    pub fn into_inner(self) -> F {
        self.function
    }
}

impl<F> LogDensity<[f64]> for FnLogDensity<F>
where
    F: FnMut(&[f64]) -> f64 + Send,
{
    fn log_density(&mut self, state: &[f64]) -> f64 {
        (self.function)(state)
    }
}

/// Validate the target-density convention used throughout the crate.
pub(crate) fn validate_log_density(value: f64) -> Result<f64, crate::McmcError> {
    if value.is_nan() || value == f64::INFINITY {
        Err(crate::McmcError::InvalidLogDensity { value })
    } else {
        Ok(value)
    }
}
