//! Local Hilbert-space primitives for impurity solvers.

use crate::impurity::ImpurityError;

/// Spin state stored on a spin-1/2 worldline leg (`-1` or `+1`, representing
/// the eigenvalue of `sigma_z`).
pub type Spin = i8;

/// Validate the compact spin-1/2 encoding.
pub fn validate_spin(spin: Spin, field: &str) -> Result<(), ImpurityError> {
    if !matches!(spin, -1 | 1) {
        return Err(ImpurityError::parameter(field, "must be -1 or +1"));
    }
    Ok(())
}
