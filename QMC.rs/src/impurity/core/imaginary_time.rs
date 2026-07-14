//! Periodic imaginary-time helpers shared by impurity solvers.

use crate::impurity::ImpurityError;

/// Direction in which a loop head propagates along the periodic worldline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropagationDirection {
    /// Increasing imaginary time.
    Forward,
    /// Decreasing imaginary time.
    Backward,
}

/// Validate an inverse temperature.
pub fn validate_beta(beta: f64) -> Result<(), ImpurityError> {
    if !beta.is_finite() || beta <= 0.0 {
        return Err(ImpurityError::parameter(
            "beta",
            format!("must be finite and positive, got {beta}"),
        ));
    }
    Ok(())
}

/// Wrap a time into `[0, beta)`.
pub fn wrap_tau(tau: f64, beta: f64) -> f64 {
    tau.rem_euclid(beta)
}

/// Directed periodic distance from `from` to `to`.
pub fn directed_distance(from: f64, to: f64, beta: f64, direction: PropagationDirection) -> f64 {
    match direction {
        PropagationDirection::Forward => (to - from).rem_euclid(beta),
        PropagationDirection::Backward => (from - to).rem_euclid(beta),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directed_distances_wrap_periodically() {
        assert!(
            (directed_distance(3.0, 1.0, 4.0, PropagationDirection::Forward) - 2.0).abs() < 1e-14
        );
        assert!(
            (directed_distance(1.0, 3.0, 4.0, PropagationDirection::Backward) - 2.0).abs() < 1e-14
        );
    }
}
