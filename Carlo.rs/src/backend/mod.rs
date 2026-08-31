//! Parallel execution backends for Monte Carlo simulations.
//!
//! This module provides the [`Backend`] trait for abstracting parallel execution,
//! allowing simulations to run on different parallel architectures.
//!
//! # Available Backends
//!
//! - [`RayonBackend`]: Thread-parallel execution using Rayon (default)
//!
//! ## MPI Backend (requires `mpi` feature)
//!
//! When compiled with `--features mpi`, an MPI backend is available for
//! distributed execution across multiple nodes.
//!
//! # Example
//!
//! ```rust
//! use carlo_rs::{Backend, RayonBackend, Scheduler, RunConfig};
//!
//! // Create a Rayon backend for thread-parallel execution
//! let backend = RayonBackend::new(4); // 4 threads
//!
//! // Use with scheduler
//! let config = RunConfig {
//!     thermalization_sweeps: 1000,
//!     measurement_sweeps: 10000,
//!     binsize: 100,
//!     base_seed: 42,
//!     progress_interval: 1000,
//!     checkpoint_interval: 0,
//! };
//! let scheduler = Scheduler::new(backend, config);
//! ```
//!
//! # MPI Backend
//!
//! The MPI backend requires the `mpi` feature and MPI installation:
//!
//! ```bash
//! # Install MPI (Ubuntu/Debian)
//! sudo apt-get install libopenmpi-dev openmpi-bin
//!
//! # Build with MPI
//! cargo build --features mpi
//!
//! # Run with mpirun
//! mpirun -np 16 ./carlo-rs run
//! ```

use rand_core::Rng;
use rand_core::SeedableRng;

/// Parallel execution backend abstraction.
pub trait Backend: Clone + Send + Sync {
    /// RNG type for this backend.
    type Rng: Rng + SeedableRng + Send;

    /// Spawn n tasks in parallel, each with isolated RNG.
    fn spawn_tasks<F>(&self, n_tasks: usize, base_seed: u64, f: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync;

    /// Wait for all tasks to complete.
    fn barrier(&self);
}

mod rayon;
pub use rayon::RayonBackend;

#[cfg(feature = "mpi")]
mod mpi;
#[cfg(feature = "mpi")]
pub use mpi::{
    run_distributed, run_distributed_compat, DistributedConfig, Done, Idle, MpiBackend, MpiError,
    MpiRng, MpiRunConfig, ResultsAggregator, Running, SchedulerTask, TaskSpec, TaskStream,
    TimeLimits, Worker, WorkerState,
};
