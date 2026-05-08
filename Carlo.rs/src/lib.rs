//! # Carlo.rs - Monte Carlo Simulation Framework
//!
//! Carlo.rs is a Rust implementation of the [Carlo.jl](https://github.com/lukas-weber/Carlo.jl)
//! framework for developing high-performance, distributed Monte Carlo simulations.
//!
//! ## Overview
//!
//! The framework handles model-independent tasks:
//! - Autocorrelation and error analysis
//! - Monte-Carlo-aware MPI scheduling
//! - Checkpointing and result merging
//!
//! while leaving all flexibility of implementing Monte Carlo updates and estimators to you.
//!
//! ## Quick Start
//!
//! Implement the [`MonteCarlo`] trait for your model:
//!
//! ```rust,ignore
//! use carlo_rs::{MonteCarlo, Context, CarloError, FromParams, Params};
//! use rand_xoshiro::Xoshiro256PlusPlus;
//! use rand_core::Rng;
//!
//! // Your Monte Carlo model
//! struct IsingModel {
//!     lattice: Vec<i8>,
//!     beta: f64,  // inverse temperature
//! }
//!
//! impl MonteCarlo for IsingModel {
//!     type Rng = Xoshiro256PlusPlus;
//!
//!     fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
//!         // Perform one Monte Carlo sweep (update configuration)
//!         for i in 0..self.lattice.len() {
//!             // Metropolis update
//!             let neighbor_sum = self.lattice[(i + 1) % self.lattice.len()]
//!                 + self.lattice[(i - 1 + self.lattice.len()) % self.lattice.len()];
//!             let delta_e = -2.0 * self.lattice[i] as f64 * neighbor_sum as f64 * self.beta;
//!             if delta_e < 0.0 || ctx.rng.gen::<f64>() < (-delta_e).exp() {
//!                 self.lattice[i] *= -1;
//!             }
//!         }
//!     }
//!
//!     fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
//!         // Measure observables after each sweep
//!         let magnetization = self.lattice.iter().sum::<i8>() as f64 / self.lattice.len() as f64;
//!         ctx.measure("Magnetization", magnetization);
//!     }
//! }
//!
//! impl FromParams for IsingModel {
//!     fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
//!         let size = params.get::<usize>("L").unwrap_or(100);
//!         let beta = params.get::<f64>("beta").unwrap_or(1.0);
//!         Ok(Self {
//!             lattice: vec![1; size],
//!             beta,
//!         })
//!     }
//! }
//! ```
//!
//! ## Features
//!
//! - `hdf5`: HDF5 checkpoint and measurement file support
//! - `mpi`: MPI distributed backend for large-scale parallel simulations
//! - `strict-repro`: Strict reproducibility mode using jump sequences for RNG
//!
//! ## Architecture
//!
//! - [`MonteCarlo`]: Core trait for implementing simulation algorithms
//! - [`Backend`]: Parallel execution abstraction (Rayon or MPI)
//! - [`Scheduler`]: Orchestrates thermalization and measurement phases
//! - [`Context`]: Runtime state including RNG and measurement collection
//! - [`Results`]: Final simulation results with metadata
//!
//! For result analysis, see [`merge`] and [`evaluable`] modules.

pub mod backend;
pub mod cli;
mod context;
mod error;
mod estimate;
pub mod evaluable;
pub mod job;
pub mod lattice;
mod measurements;
pub mod merge;
mod monte_carlo;
pub mod output;
pub mod parallel_tempering;
mod params;
pub mod progress;
mod results;
mod rng_checkpoint;
pub mod run;
mod scheduler;
mod version;

#[cfg(feature = "mpi")]
pub use backend::{run_distributed, MpiBackend, MpiRunConfig, SchedulerTask};
pub use backend::{Backend, RayonBackend};
pub use context::{Context, ContextCheckpoint};
pub use error::CarloError;
pub use estimate::{ComplexEstimate, Estimate};
pub use evaluable::{jackknife, Evaluator, MultiplexEvaluator};
pub use job::{parse_duration, JobInfo, TaskInfo, TaskMaker};
pub use measurements::{
    Accumulator, ComplexAccumulator, ComplexValue, Measurements,
};
pub use merge::{
    AddSamplesState, calc_rebin_count, calc_rebin_length, compute_decorrelated_autocorr_time,
    compute_regular_autocorr_time, cov_of_mean, list_meas_files,
    MergeOptions, ObservableType, ResultObservable,
};
#[cfg(feature = "hdf5")]
pub use merge::{iterate_measfile_observables, merge_results, merge_results_from_files};
#[cfg(feature = "hdf5")]
pub use monte_carlo::MonteCarloCheckpoint;
pub use monte_carlo::{FromParams, MonteCarlo, MonteCarloExt};
pub use progress::{print_status_table, spinner, MultiTaskProgress, SimProgress};
pub use rng_checkpoint::{RngCheckpointHdf5, RNG_TYPE, RNG_VERSION};
pub use run::{Run, RunId, TaskId, timing};
pub use version::Version;

// Re-export CLI types for library users
pub use cli::run as cli_run;
pub use output::{save_hdf5, save_json};
pub use parallel_tempering::{
    ParallelTemperingCompatible, ParallelTemperingConfig, ParallelTemperingMC,
};
pub use params::Params;
pub use results::{ComplexResult, Metadata, Results};
pub use scheduler::{RunConfig, Scheduler};
