use thiserror::Error;

#[derive(Debug, Error)]
pub enum CarloError {
    #[error("Invalid configuration: field '{field}' - {reason}")]
    InvalidConfig { field: String, reason: String },

    #[cfg(feature = "hdf5")]
    #[error("HDF5 I/O error: {0}")]
    Hdf5Error(#[from] hdf5::Error),

    #[error("I/O error for path '{path}': {source}")]
    IoError {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("Measurement '{name}' not found")]
    MeasurementNotFound { name: String },

    #[error("Checkpoint corrupted: {detail}")]
    CheckpointCorrupted { detail: String },

    #[error("Convergence not reached after {sweeps} sweeps")]
    ConvergenceTimeout { sweeps: u64 },

    #[error("Serialization error: {0}")]
    SerializationError(#[from] serde_json::Error),

    #[cfg(feature = "mpi")]
    #[error("MPI error: {0}")]
    MpiError(#[from] crate::backend::MpiError),
}
