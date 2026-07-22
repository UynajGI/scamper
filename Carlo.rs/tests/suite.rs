//! Carlo.rs integration test suite.
//!
//! Each module lives in a subdirectory and is wired via explicit `#[path]`
//! attributes so the directory tree mirrors the test taxonomy without
//! requiring `mod.rs` files.

// ── Unit ──────────────────────────────────────────────────────────────────

#[path = "unit/accumulator.rs"]
mod accumulator;
#[path = "unit/autocorr_reference.rs"]
mod autocorr_reference;
#[path = "unit/clock_phase.rs"]
mod clock_phase;
#[path = "unit/complex.rs"]
mod complex;
#[path = "unit/context.rs"]
mod context;
#[path = "unit/error.rs"]
mod error;
#[path = "unit/estimate.rs"]
mod estimate;
#[path = "unit/evaluable.rs"]
mod evaluable;
#[path = "unit/measurements.rs"]
mod measurements;
#[path = "unit/merge.rs"]
mod merge;
#[path = "unit/params.rs"]
mod params;
#[path = "unit/results.rs"]
mod results;
#[path = "unit/version.rs"]
mod version;

// ── Integration ───────────────────────────────────────────────────────────

#[path = "integration/backend.rs"]
mod backend;
#[path = "integration/checkpoint.rs"]
mod checkpoint;
#[path = "integration/cli.rs"]
mod cli;
#[path = "integration/lifecycle.rs"]
mod lifecycle;
#[path = "integration/monte_carlo.rs"]
mod monte_carlo;
#[path = "integration/parallel_tempering.rs"]
mod parallel_tempering;
#[path = "integration/reproducibility.rs"]
mod reproducibility;
#[path = "integration/run.rs"]
mod run;
#[path = "integration/run_control.rs"]
mod run_control;
#[path = "integration/scheduler.rs"]
mod scheduler;
#[path = "integration/workflow.rs"]
mod workflow;

// ── I/O ───────────────────────────────────────────────────────────────────

#[cfg(feature = "hdf5")]
#[path = "io/checkpoint_hdf5.rs"]
mod checkpoint_hdf5;
#[path = "io/dataframe.rs"]
mod dataframe;
#[path = "io/job.rs"]
mod job;
#[path = "io/merge_io.rs"]
mod merge_io;
#[path = "io/output_io.rs"]
mod output_io;

// ── MPI (feature-gated) ───────────────────────────────────────────────────

#[path = "mpi/distributed.rs"]
mod mpi_distributed;
#[path = "mpi/mpi_test.rs"]
mod mpi_test;

// ── Performance ───────────────────────────────────────────────────────────

#[path = "perf/perf_test.rs"]
mod perf_test;
