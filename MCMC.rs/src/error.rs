use thiserror::Error;

/// Errors produced by target evaluation, kernel configuration, diagnostics or I/O.
#[derive(Debug, Error)]
pub enum McmcError {
    #[error("invalid MCMC configuration: {0}")]
    InvalidConfig(String),

    #[error("dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("target returned invalid log density {value}")]
    InvalidLogDensity { value: f64 },

    #[error("warmup adaptation was used after it had been frozen")]
    AdaptationFrozen,

    #[error("diagnostics require at least two chains with at least four draws each")]
    InsufficientDraws,

    #[error("all chains must have the same parameter dimension")]
    InconsistentTraceDimension,

    #[error("checkpoint format mismatch: expected {expected}, found {found}")]
    CheckpointFormat { expected: String, found: String },

    #[error("checkpoint target fingerprint does not match the reconstructed target")]
    TargetMismatch,

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("serialization error: {0}")]
    Serialization(#[from] serde_json::Error),

    #[cfg(feature = "hdf5")]
    #[error("HDF5 error: {0}")]
    Hdf5(#[from] hdf5::Error),
}
