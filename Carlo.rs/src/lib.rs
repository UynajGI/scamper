//! # Carlo.rs — Monte Carlo Simulation Framework
//!
//! Carlo.rs handles the **model-independent** concerns of Monte Carlo simulations:
//! scheduling, measurement accumulation, error analysis, checkpointing, and
//! parallel execution. You implement the physics — Carlo.rs runs it.
//!
//! ## Trait hierarchy
//!
//! ```text
//! MonteCarlo          — you MUST implement: sweep() + type Rng
//!  ├─ measure()       — you SHOULD implement: record observables via ctx
//!  ├─ name()          — optional: label for logging
//!  └─ MonteCarloExt   — blanket impl: init(), register_evaluables()
//!
//! FromParams: MonteCarlo   — you MUST implement: construct model from Params
//! ```
//!
//! ## Execution flow
//!
//! ```text
//! Params (HashMap<String, String>)
//!   │
//!   ▼ FromParams::from_params(params, rng)
//! YourModel : MonteCarlo
//!   │
//!   ▼ Scheduler::new(backend, RunConfig)::run_one::<YourModel>(params)
//!   │
//!   ├─ thermalization (N sweeps):  sweep() + advance_sweep()
//!   │
//!   ├─ measurement (M sweeps):     sweep() + measure() + advance_sweep()
//!   │
//!   └─ finalize_measurements() → Results
//!        │
//!        └─ Estimates { mean, stderr, autocorr_time, n_bins }
//! ```
//!
//! ## Quick start
//!
//! A complete 1D Ising chain simulation from definition to results:
//!
//! ```rust,ignore
//! use carlo_rs::{
//!     accept_log_probability, MonteCarlo, FromParams, Context, Params, CarloError,
//!     Scheduler, RunConfig, RayonBackend, Backend,
//! };
//! use rand_xoshiro::Xoshiro256PlusPlus;
//!
//! // ── Step 1: define your model ──
//! struct Ising1D {
//!     spins: Vec<i8>,
//!     beta: f64,
//! }
//!
//! // ── Step 2: implement MonteCarlo ──
//! impl MonteCarlo for Ising1D {
//!     type Rng = Xoshiro256PlusPlus;
//!
//!     fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
//!         let n = self.spins.len();
//!         for i in 0..n {
//!             let left  = self.spins[(i + n - 1) % n] as f64;
//!             let right = self.spins[(i + 1) % n] as f64;
//!             let dE = 2.0 * self.beta * self.spins[i] as f64 * (left + right);
//!             if accept_log_probability(-dE, &mut ctx.rng) {
//!                 self.spins[i] *= -1;
//!             }
//!         }
//!     }
//!
//!     fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
//!         let m = self.spins.iter().map(|&s| s as f64).sum::<f64>() / self.spins.len() as f64;
//!         ctx.measure("Magnetization", m);
//!     }
//! }
//!
//! // ── Step 3: implement FromParams ──
//! impl FromParams for Ising1D {
//!     fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
//!         let n = params.get::<usize>("L").unwrap_or(100);
//!         // Random initial configuration
//!         let spins: Vec<i8> = (0..n)
//!             .map(|_| if rng.gen::<bool>() { 1 } else { -1 })
//!             .collect();
//!         let beta = params.get::<f64>("beta").unwrap_or(1.0);
//!         Ok(Self { spins, beta })
//!     }
//! }
//!
//! // ── Step 4: run ──
//! let mut params = Params::new();
//! params.set("L", 128);
//! params.set("beta", 0.5);
//!
//! let config = RunConfig {
//!     thermalization_sweeps: 1000,
//!     measurement_sweeps: 10_000,
//!     binsize: 100,
//!     base_seed: 42,
//!     ..Default::default()
//! };
//!
//! let backend = RayonBackend::new(1);
//! let scheduler = Scheduler::new(backend, config);
//! let results: Results = scheduler.run_one::<Ising1D>(&params);
//!
//! // ── Step 5: read results ──
//! if let Some(est) = results.get("Magnetization") {
//!     println!("Magnetization = {}", est.format());
//! }
//! ```
//!
//! ## Key types
//!
//! | Type | Role |
//! |------|------|
//! | [`MonteCarlo`] | You implement this. One sweep = one MC update pass. |
//! | [`FromParams`] | You implement this. Construct model from `Params` dict. |
//! | [`Context`] | Passed to `sweep()` / `measure()`. Holds RNG + measurement buffer. |
//! | [`Params`] | `HashMap<String, String>` parameter bag. `set()` / `get::<T>()`. |
//! | [`Scheduler`] | Runs thermalization → measurement loops. |
//! | [`RunConfig`] | Sweep counts, binsize, seed, progress/checkpoint intervals. |
//! | [`Backend`] | Parallel execution. [`RayonBackend`] for threads, `MpiBackend` for MPI. |
//! | [`Results`] | `HashMap<String, Estimate>`. Serialize to JSON, merge across tasks. |
//! | [`Estimate`] | `{ mean, stderr, autocorr_time, n_bins }` |
//!
//! ## Features
//!
//! | Feature | Effect |
//! |---------|--------|
//! | `hdf5` | HDF5 checkpoint + measurement files |
//! | `mpi` | MPI distributed backend (`mpirun -np N ./carlo run`) |
//! | `strict-repro` | Jump-sequence RNG for exact reproducibility across task counts |
//!
//! ## Module map
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`monte_carlo`] | `MonteCarlo`, `FromParams`, `MonteCarloExt` traits |
//! | [`context`] | `Context<R>` — RNG + measurements + sweep counter |
//! | [`scheduler`] | `Scheduler<B>`, `RunConfig` — run orchestration |
//! | [`backend`] | `Backend` trait, `RayonBackend`, `MpiBackend` |
//! | [`measurements`] | `Measurements`, `Accumulator` — binned sample collection |
//! | [`merge`] | Rebinning, autocorr time, `merge_results` |
//! | [`evaluable`] | Jackknife resampling, `Evaluator` |
//! | [`results`] | `Results`, `Metadata`, `Estimate`, `ComplexEstimate` |
//! | [`params`] | `Params` — typed key-value parameter store |
//! | [`lattice`] | `LatticeParams` — basic 2D lattice helpers |
//! | [`run`] | `Run`, `RunId`, `TaskId` — single-run lifecycle |
//! | [`output`] | JSON/HDF5 save/load, `dataframe()` |
//! | [`parallel_tempering`] | PT MC with chain scheduling |
//! | [`job`] | `JobInfo`, `TaskInfo`, `TaskMaker` — multi-run job management |
//! | [`cli`] | `carlo run/status/merge/delete` CLI |
//! | [`progress`] | Progress bars and status tables |

mod acceptance;
pub mod backend;
pub mod cli;
mod clock;
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
mod phase;
pub mod progress;
mod results;
mod rng_checkpoint;
pub mod rng_stream;
pub mod run;
mod run_control;
mod scheduler;
mod version;

pub use acceptance::accept_log_probability;
#[cfg(feature = "mpi")]
pub use backend::{
    run_distributed, run_distributed_compat, DistributedConfig, MpiBackend, MpiError, MpiRng,
    MpiRunConfig, SchedulerTask, TaskSpec,
};
pub use backend::{Backend, RayonBackend};
pub use clock::SimulationClock;
pub use context::{Context, ContextCheckpoint};
pub use error::CarloError;
pub use estimate::{ComplexEstimate, Estimate};
pub use evaluable::{jackknife, Evaluator, MultiplexEvaluator};
pub use job::{parse_duration, JobInfo, TaskInfo, TaskMaker};
pub use measurements::{Accumulator, ComplexAccumulator, ComplexValue, Measurements};
pub use merge::{
    calc_rebin_count, calc_rebin_length, compute_decorrelated_autocorr_time,
    compute_regular_autocorr_time, cov_of_mean, list_meas_files, AddSamplesState, MergeOptions,
    ObservableType, ResultObservable,
};
#[cfg(feature = "hdf5")]
pub use merge::{
    iterate_measfile_observables, merge_results, merge_results_from_files, merge_task_results,
};
#[cfg(feature = "hdf5")]
pub use monte_carlo::MonteCarloCheckpoint;
pub use monte_carlo::{FromParams, MonteCarlo, MonteCarloExt};
pub use progress::{print_status_table, spinner, MultiTaskProgress, SimProgress};
pub use rng_checkpoint::{RngCheckpointHdf5, RNG_TYPE, RNG_VERSION};
pub use run::{timing, Run, RunId, TaskId};
pub use version::Version;

// Re-export CLI types for library users
pub use cli::run as cli_run;
pub use output::{
    dataframe, make_scalar, make_scalar_owned, measurement_from_obs, recursive_stack, save_hdf5,
    save_json, ResultRow,
};
pub use parallel_tempering::{
    ParallelTemperingCompatible, ParallelTemperingConfig, ParallelTemperingMC,
};
pub use params::Params;
pub use phase::RunPhase;
pub use results::{ComplexResult, Metadata, Results};
pub use rng_stream::{RngPhase, RngStreamKey};
pub use run_control::{AdaptiveRunControl, RunDecision};
pub use scheduler::{RunConfig, Scheduler};
