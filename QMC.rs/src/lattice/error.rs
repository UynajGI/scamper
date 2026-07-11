//! Error types for continuous-time lattice QMC.

use thiserror::Error;

use crate::graph::GraphError;
use crate::local_space::LocalSpaceError;

/// Runtime, model-construction, or configuration error.
#[derive(Debug, Error)]
pub enum LatticeQmcError {
    /// Invalid user parameter.
    #[error("invalid parameter `{field}`: {reason}")]
    InvalidParameter {
        /// Parameter name.
        field: String,
        /// Validation message.
        reason: String,
    },
    /// A model cannot be represented with positive local weights.
    #[error("invalid or non-stoquastic model: {0}")]
    InvalidModel(String),
    /// A sampled configuration violates a representation invariant.
    #[error("invalid continuous-time configuration: {0}")]
    InvalidConfiguration(String),
    /// No diagonal matrix element exists for a proposed insertion.
    #[error("term {term} has no diagonal vertex for local states {states:?}")]
    MissingDiagonalVertex {
        /// Term index.
        term: usize,
        /// Local basis states.
        states: Vec<u16>,
    },
    /// A directed loop exceeded its safety limit.
    #[error("directed loop did not close after {steps} steps (limit {limit})")]
    LoopDidNotClose {
        /// Performed steps.
        steps: usize,
        /// Configured limit.
        limit: usize,
    },
    /// No valid worm discontinuity could be inserted.
    #[error("no valid directed-loop start state exists in the current configuration")]
    NoLoopStart,
    /// Graph-construction error.
    #[error(transparent)]
    Graph(#[from] GraphError),
    /// Local-space error.
    #[error(transparent)]
    LocalSpace(#[from] LocalSpaceError),
}

impl LatticeQmcError {
    /// Construct an invalid-parameter error.
    pub fn parameter(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidParameter {
            field: field.into(),
            reason: reason.into(),
        }
    }
}
