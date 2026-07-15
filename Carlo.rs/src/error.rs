//! Error types for the Carlo.rs framework.

use thiserror::Error;

/// Errors returned by Carlo.rs operations.
#[derive(Debug, Error)]
pub enum CarloError {
    /// A configuration field has an invalid or missing value.
    #[error("Invalid configuration: field '{field}' - {reason}")]
    InvalidConfig { field: String, reason: String },

    /// HDF5 read/write failure (only available with the `hdf5` feature).
    #[cfg(feature = "hdf5")]
    #[error("HDF5 I/O error: {0}")]
    Hdf5Error(#[from] hdf5::Error),

    /// File-system I/O error with the offending path.
    #[error("I/O error for path '{path}': {source}")]
    IoError {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// A requested observable was not found in the measurement set.
    #[error("Measurement '{name}' not found")]
    MeasurementNotFound { name: String },

    /// A checkpoint file is corrupt or incompatible with the expected schema.
    #[error("Checkpoint corrupted: {detail}")]
    CheckpointCorrupted { detail: String },

    /// An adaptive run did not converge within the allotted sweep budget.
    #[error("Convergence not reached after {sweeps} sweeps")]
    ConvergenceTimeout { sweeps: u64 },

    /// (De)serialization failure (typically JSON or HDF5 attribute parsing).
    #[error("Serialization error: {0}")]
    SerializationError(#[from] serde_json::Error),

    /// MPI collective or topology failure (only with the `mpi` feature).
    #[cfg(feature = "mpi")]
    #[error("MPI error: {0}")]
    MpiError(#[from] crate::backend::MpiError),
}
