//! Error type for the variational QMC family.

use thiserror::Error;

/// Errors produced by variational constructors, kernels and checkpoints.
///
/// Invalid user input is always reported as [`VariationalError::InvalidConfig`]
/// (input-validation criterion G: never panic, never silently accept garbage),
/// and unreadable snapshots as [`VariationalError::CheckpointCorrupted`].
#[derive(Debug, Error)]
pub enum VariationalError {
    /// A configuration field has an invalid or missing value.
    #[error("invalid variational configuration: field '{field}' - {reason}")]
    InvalidConfig {
        /// Field that failed validation.
        field: String,
        /// Why the value was rejected.
        reason: String,
    },

    /// A checkpoint snapshot is corrupt or uses an unknown format.
    #[error("variational checkpoint corrupted: {detail}")]
    CheckpointCorrupted {
        /// What exactly failed to parse or validate.
        detail: String,
    },
}

impl VariationalError {
    /// Shorthand for [`VariationalError::InvalidConfig`].
    pub fn invalid(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidConfig {
            field: field.into(),
            reason: reason.into(),
        }
    }

    /// Shorthand for [`VariationalError::CheckpointCorrupted`].
    pub fn checkpoint(detail: impl Into<String>) -> Self {
        Self::CheckpointCorrupted {
            detail: detail.into(),
        }
    }

    /// Require a finite, strictly positive parameter (criterion G).
    pub fn require_positive(field: &str, value: f64) -> Result<(), Self> {
        if value.is_finite() && value > 0.0 {
            Ok(())
        } else {
            Err(Self::invalid(
                field,
                format!("must be finite and strictly positive, got {value}"),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers_build_the_documented_variants() {
        assert!(matches!(
            VariationalError::invalid("alpha", "non-finite"),
            VariationalError::InvalidConfig { .. }
        ));
        assert!(matches!(
            VariationalError::checkpoint("bad tag"),
            VariationalError::CheckpointCorrupted { .. }
        ));
        assert!(VariationalError::require_positive("b", 0.5).is_ok());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(VariationalError::require_positive("b", bad).is_err());
        }
    }
}
