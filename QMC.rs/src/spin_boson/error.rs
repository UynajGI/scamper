//! Errors produced by the continuous-time spin-boson engine.

use thiserror::Error;

/// Spin-boson QMC construction and runtime errors.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum SpinBosonError {
    /// A physical or algorithmic parameter is outside its valid domain.
    #[error("invalid parameter `{field}`: {reason}")]
    InvalidParameter {
        /// Parameter name.
        field: String,
        /// Human-readable reason.
        reason: String,
    },

    /// A sampled operator configuration violates worldline invariants.
    #[error("invalid wormhole configuration: {0}")]
    InvalidConfiguration(String),

    /// The directed loop exceeded its safety limit without closing.
    #[error("directed loop did not close after {steps} steps (limit {limit})")]
    LoopDidNotClose {
        /// Steps executed.
        steps: usize,
        /// Configured safety limit.
        limit: usize,
    },

    /// A tabulated bath was malformed.
    #[error("invalid tabulated bath: {0}")]
    InvalidBathTable(String),
}

impl SpinBosonError {
    /// Convenience constructor for invalid parameters.
    pub fn parameter(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidParameter {
            field: field.into(),
            reason: reason.into(),
        }
    }
}
