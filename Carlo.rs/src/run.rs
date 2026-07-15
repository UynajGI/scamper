//! Run struct for simulation lifecycle management.
//!
//! [`Run`] encapsulates a single Monte Carlo simulation run,
//! providing checkpoint/restart capability and lifecycle tracking.
//!
//! # Usage
//!
//! ```rust,ignore
//! let mut run = Run::new(params, task_id, run_id, &config, seed)?;
//! while !run.is_complete() {
//!     run.step();
//!     if should_checkpoint {
//!         run.write_checkpoint(&path)?;
//!     }
//! }
//! let results = run.finalize(base_seed);
//! ```

use rand_core::{Rng, SeedableRng};
use std::time::Instant;

#[cfg(feature = "mpi")]
use mpi::topology::SimpleCommunicator;
#[cfg(feature = "hdf5")]
use std::path::Path;

use crate::{
    CarloError, Context, FromParams, Metadata, MonteCarlo, Params, Results, RunConfig, RunPhase,
};

#[cfg(feature = "hdf5")]
use hdf5::File as Hdf5File;

/// Internal timing observables (prefixed with `_ll_` matching Carlo.jl).
pub mod timing {
    /// Wall-clock time of one sweep, in seconds.
    pub const SWEEP_TIME: &str = "_ll_sweep_time";
    /// Wall-clock time of one measurement pass, in seconds.
    pub const MEASURE_TIME: &str = "_ll_measure_time";
    /// Wall-clock time to read a checkpoint, in seconds.
    pub const CHECKPOINT_READ_TIME: &str = "_ll_checkpoint_read_time";
    /// Wall-clock time to write a checkpoint, in seconds.
    pub const CHECKPOINT_WRITE_TIME: &str = "_ll_checkpoint_write_time";
}

/// Run ID for tracking multiple runs of the same task.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunId(pub u64);

impl RunId {
    /// Create a new run ID.
    pub fn new(id: u64) -> Self {
        Self(id)
    }

    /// Underlying `u64` value.
    pub fn as_u64(&self) -> u64 {
        self.0
    }
}

/// Task ID for identifying parameter sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TaskId(pub usize);

impl TaskId {
    /// Create a new task ID.
    pub fn new(id: usize) -> Self {
        Self(id)
    }

    /// Underlying `usize` value.
    pub fn as_usize(&self) -> usize {
        self.0
    }
}

/// Single Monte Carlo simulation run.
pub struct Run<MC: MonteCarlo, R: Rng + SeedableRng> {
    /// Runtime context (RNG, measurements, sweep counter).
    context: Context<R>,

    /// Monte Carlo implementation.
    mc: MC,

    /// Task ID (parameter set index).
    task_id: TaskId,

    /// Run ID (unique within task).
    run_id: RunId,

    /// Target sweeps for this run.
    target_sweeps: u64,

    /// Completed sweeps.
    sweeps_done: u64,

    /// Run configuration.
    config: RunConfig,
}

impl<MC: MonteCarlo<Rng = R>, R: Rng + SeedableRng + Send> Run<MC, R> {
    /// Create a new run.
    pub fn new(
        params: &Params,
        task_id: TaskId,
        run_id: RunId,
        config: &RunConfig,
        seed: u64,
    ) -> Result<Self, CarloError>
    where
        MC: FromParams<Rng = R>,
    {
        let rng = R::seed_from_u64(seed);
        let mut context =
            Context::new_with_binsize(rng, config.thermalization_sweeps, config.binsize);
        let mut mc = MC::from_params(params, &mut context.rng)?;
        let initial_phase = if config.thermalization_sweeps == 0 {
            RunPhase::Measurement
        } else {
            RunPhase::Thermalization
        };
        context.enter_phase(initial_phase);
        mc.on_phase_start(initial_phase, &mut context);

        Ok(Self {
            context,
            mc,
            task_id,
            run_id,
            target_sweeps: config.measurement_sweeps,
            sweeps_done: 0,
            config: config.clone(),
        })
    }

    /// Create a configured run from an already constructed Monte Carlo value.
    ///
    /// This avoids forcing models with closures, shared datasets or other
    /// non-`Params` construction paths to implement [`crate::FromParams`].
    pub fn from_parts(
        mut context: Context<R>,
        mut mc: MC,
        task_id: TaskId,
        run_id: RunId,
        config: RunConfig,
    ) -> Self {
        let initial_phase = if config.thermalization_sweeps == 0 {
            RunPhase::Measurement
        } else {
            RunPhase::Thermalization
        };
        context.enter_phase(initial_phase);
        mc.on_phase_start(initial_phase, &mut context);
        Self {
            context,
            mc,
            task_id,
            run_id,
            target_sweeps: config.measurement_sweeps,
            sweeps_done: 0,
            config,
        }
    }

    /// Create run from existing context and MC (backward compatibility).
    pub fn from_context(context: Context<R>, mc: MC) -> Self {
        Self {
            context,
            mc,
            task_id: TaskId::new(0),
            run_id: RunId::new(0),
            target_sweeps: 0,
            sweeps_done: 0,
            config: RunConfig::default(),
        }
    }

    /// Execute one Monte Carlo step (sweep + optional measurement).
    /// Returns number of measurement sweeps (1 if thermalized, 0 otherwise).
    ///
    /// Automatically records sweep_time and measure_time as observables
    /// (prefixed with `_ll_` matching Carlo.jl convention).
    pub fn step(&mut self) -> u64 {
        let desired_phase = if self.context.sweep_count() < self.config.thermalization_sweeps {
            RunPhase::Thermalization
        } else {
            RunPhase::Measurement
        };
        if self.context.phase() != desired_phase {
            let previous = self.context.phase();
            self.mc.on_phase_end(previous, &mut self.context);
            self.context.enter_phase(desired_phase);
            self.mc.on_phase_start(desired_phase, &mut self.context);
        }

        let sweep_start = Instant::now();
        self.mc.sweep(&mut self.context);
        let sweep_time = sweep_start.elapsed().as_secs_f64();
        let collect = self.context.phase().collects_measurements();
        self.context.advance_sweep();

        if collect {
            let measure_start = Instant::now();
            self.mc.measure(&mut self.context);
            let measure_time = measure_start.elapsed().as_secs_f64();

            // Record timing observables
            self.context.measure(timing::SWEEP_TIME, sweep_time);
            self.context.measure(timing::MEASURE_TIME, measure_time);

            self.sweeps_done += 1;
            1
        } else {
            0
        }
    }

    /// Execute one Monte Carlo step with MPI communicator for multi-rank coordination.
    /// Default behavior (via MonteCarlo::sweep_with_comm) delegates to regular sweep.
    #[cfg(feature = "mpi")]
    pub fn step_with_comm(&mut self, comm: &SimpleCommunicator) -> u64 {
        let desired_phase = if self.context.sweep_count() < self.config.thermalization_sweeps {
            RunPhase::Thermalization
        } else {
            RunPhase::Measurement
        };
        if self.context.phase() != desired_phase {
            let previous = self.context.phase();
            self.mc.on_phase_end(previous, &mut self.context);
            self.context.enter_phase(desired_phase);
            self.mc.on_phase_start(desired_phase, &mut self.context);
        }

        let sweep_start = Instant::now();
        self.mc.sweep_with_comm(&mut self.context, comm);
        let sweep_time = sweep_start.elapsed().as_secs_f64();
        let collect = self.context.phase().collects_measurements();
        self.context.advance_sweep();

        if collect {
            let measure_start = Instant::now();
            self.mc.measure_with_comm(&mut self.context, comm);
            let measure_time = measure_start.elapsed().as_secs_f64();

            self.context.measure(timing::SWEEP_TIME, sweep_time);
            self.context.measure(timing::MEASURE_TIME, measure_time);

            self.sweeps_done += 1;
            1
        } else {
            0
        }
    }

    /// Run multiple sweeps.
    pub fn run(&mut self, sweeps: u64) {
        for _ in 0..sweeps {
            self.step();
        }
    }

    /// Check if thermalized.
    pub fn is_thermalized(&self) -> bool {
        self.context.is_thermalized()
    }

    /// Get sweeps completed.
    pub fn sweeps_done(&self) -> u64 {
        self.sweeps_done
    }

    /// Get total sweep count (including thermalization).
    pub fn sweep_count(&self) -> u64 {
        self.context.sweep_count()
    }

    /// Check if run is complete.
    pub fn is_complete(&self) -> bool {
        self.sweeps_done >= self.target_sweeps
    }

    /// Get configured measurement-sweep target.
    pub fn target_sweeps(&self) -> u64 {
        self.target_sweeps
    }

    /// Get remaining measurement sweeps.
    pub fn remaining_sweeps(&self) -> u64 {
        self.target_sweeps.saturating_sub(self.sweeps_done)
    }

    /// Get task ID.
    pub fn task_id(&self) -> TaskId {
        self.task_id
    }

    /// Get run ID.
    pub fn run_id(&self) -> RunId {
        self.run_id
    }

    /// Get context reference.
    pub fn context(&self) -> &Context<R> {
        &self.context
    }

    /// Get mutable context reference.
    pub fn context_mut(&mut self) -> &mut Context<R> {
        &mut self.context
    }

    /// Get MC implementation reference.
    pub fn mc(&self) -> &MC {
        &self.mc
    }

    /// Get mutable MC implementation reference.
    pub fn mc_mut(&mut self) -> &mut MC {
        &mut self.mc
    }

    /// Finalize run and return results.
    pub fn finalize(self, base_seed: u64) -> Results {
        self.finalize_with_mc(base_seed).0
    }

    /// Finalize a run while returning ownership of the Monte Carlo value.
    ///
    /// Statistical samplers use this to recover posterior traces that are
    /// intentionally not flattened into Carlo.rs scalar measurements.
    pub fn finalize_with_mc(mut self, base_seed: u64) -> (Results, MC) {
        let previous = self.context.phase();
        self.mc.on_phase_end(previous, &mut self.context);
        self.context.enter_phase(RunPhase::Finished);
        self.mc
            .on_phase_start(RunPhase::Finished, &mut self.context);
        let estimates = self.context.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: chrono::Utc::now(),
            base_seed,
            thermalization_sweeps: self.config.thermalization_sweeps,
            measurement_sweeps: self.sweeps_done,
            n_tasks: 1,
        });
        (results, self.mc)
    }
}

#[cfg(feature = "hdf5")]
impl<MC: MonteCarlo<Rng = R>, R: Rng + SeedableRng + Send + crate::RngCheckpointHdf5> Run<MC, R>
where
    MC: crate::MonteCarloCheckpoint,
{
    /// Write checkpoint to HDF5 file.
    ///
    /// Automatically records checkpoint_write_time as observable.
    pub fn write_checkpoint(&mut self, path: &Path) -> Result<(), CarloError> {
        let checkpoint_start = Instant::now();

        let file = Hdf5File::create(path).map_err(|e| CarloError::InvalidConfig {
            field: "checkpoint".into(),
            reason: format!("Cannot create checkpoint file {}: {}", path.display(), e),
        })?;

        // Write version info
        let version_group =
            file.create_group("version")
                .map_err(|e| CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot create version group: {}", e),
                })?;

        version_group
            .create_dataset_simple("carlo_version", &[1], &env!("CARGO_PKG_VERSION"))
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write carlo_version: {}", e),
            })?;

        // Write context
        let context_group =
            file.create_group("context")
                .map_err(|e| CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot create context group: {}", e),
                })?;

        // Create rank-specific subgroup (for multi-rank runs)
        let rank_group = context_group
            .create_group(&format!("rank{:04}", 0))
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot create rank group: {}", e),
            })?;

        self.context.write_checkpoint_hdf5(&mut rank_group)?;

        // Write MC-specific state
        let mc_group = file
            .create_group("simulation")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot create simulation group: {}", e),
            })?;

        self.mc.write_checkpoint(&mut mc_group)?;

        // Write run metadata
        let meta_group = file
            .create_group("metadata")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot create metadata group: {}", e),
            })?;

        meta_group
            .create_dataset_simple("task_id", &[1], &(self.task_id.as_usize() as u64))
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write task_id: {}", e),
            })?;

        meta_group
            .create_dataset_simple("run_id", &[1], &self.run_id.as_u64())
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write run_id: {}", e),
            })?;

        meta_group
            .create_dataset_simple("sweeps_done", &[1], &self.sweeps_done)
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write sweeps_done: {}", e),
            })?;

        meta_group
            .create_dataset_simple("target_sweeps", &[1], &self.target_sweeps)
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write target_sweeps: {}", e),
            })?;

        // Record checkpoint write time
        let checkpoint_time = checkpoint_start.elapsed().as_secs_f64();
        self.context
            .measure(timing::CHECKPOINT_WRITE_TIME, checkpoint_time);

        Ok(())
    }

    /// Read checkpoint from HDF5 file and restore run state.
    ///
    /// Automatically records checkpoint_read_time as observable.
    pub fn read_checkpoint(
        path: &Path,
        params: &Params,
        config: &RunConfig,
        seed: u64,
    ) -> Result<Option<Self>, CarloError>
    where
        MC: FromParams<Rng = R>,
    {
        if !path.exists() {
            return Ok(None);
        }

        let checkpoint_start = Instant::now();

        let file = Hdf5File::open(path).map_err(|e| CarloError::InvalidConfig {
            field: "checkpoint".into(),
            reason: format!("Cannot open checkpoint file {}: {}", path.display(), e),
        })?;

        // Read metadata
        let meta_group = file
            .group("metadata")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot open metadata group: {}", e),
            })?;

        let task_id: u64 = meta_group
            .dataset("task_id")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read task_id: {}", e),
            })?
            .read_1d()
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse task_id: {}", e),
            })?[0];

        let run_id: u64 = meta_group
            .dataset("run_id")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read run_id: {}", e),
            })?
            .read_1d()
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse run_id: {}", e),
            })?[0];

        let sweeps_done: u64 = meta_group
            .dataset("sweeps_done")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read sweeps_done: {}", e),
            })?
            .read_1d()
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse sweeps_done: {}", e),
            })?[0];

        let target_sweeps: u64 = meta_group
            .dataset("target_sweeps")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read target_sweeps: {}", e),
            })?
            .read_1d()
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse target_sweeps: {}", e),
            })?[0];

        // Read context
        let context_group = file
            .group("context")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot open context group: {}", e),
            })?;

        let rank_group = context_group.group(&format!("rank{:04}", 0)).map_err(|e| {
            CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot open rank group: {}", e),
            }
        })?;

        let mut context = Context::read_checkpoint_hdf5_full(&rank_group, config.binsize)?;

        // Read MC state (if available)
        let mc_group = file
            .group("simulation")
            .map_err(|e| CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot open simulation group: {}", e),
            })?;

        // Create MC from params then restore state. Re-entering the restored
        // phase lets lifecycle-aware kernels rebuild/freeze phase-local state.
        let mut mc = MC::from_params(params, &mut context.rng)?;
        mc.read_checkpoint(&mc_group)?;
        mc.on_phase_start(context.phase(), &mut context);

        // Record checkpoint read time
        let checkpoint_time = checkpoint_start.elapsed().as_secs_f64();
        context.measure(timing::CHECKPOINT_READ_TIME, checkpoint_time);

        Ok(Some(Self {
            context,
            mc,
            task_id: TaskId::new(task_id as usize),
            run_id: RunId::new(run_id),
            target_sweeps,
            sweeps_done,
            config: config.clone(),
        }))
    }
}
