//! Monte Carlo scheduler for orchestrating simulation runs.
//!
//! The [`Scheduler`] coordinates the execution of Monte Carlo simulations,
//! managing thermalization and measurement phases, and collecting results.
//!
//! # Run Configuration
//!
//! [`RunConfig`] specifies simulation parameters:
//! - `thermalization_sweeps`: Number of sweeps before measurement
//! - `measurement_sweeps`: Number of sweeps for data collection
//! - `binsize`: Bin size for accumulating measurements
//! - `base_seed`: Base RNG seed (each task gets a derived seed)
//! - `progress_interval`: How often to report progress
//! - `checkpoint_interval`: How often to checkpoint (0 = disabled)
//!
//! # Example
//!
//! ```rust
//! use carlo_rs::{Scheduler, RayonBackend, RunConfig};
//!
//! let backend = RayonBackend::new(4);
//! let config = RunConfig {
//!     thermalization_sweeps: 1000,
//!     measurement_sweeps: 10000,
//!     binsize: 100,
//!     base_seed: 42,
//!     progress_interval: 1000,
//!     checkpoint_interval: 0,
//! };
//! let scheduler = Scheduler::new(backend, config);
//!
//! // Run a single simulation
//! // let results = scheduler.run_one::<MyModel>(&params);
//!
//! // Run multiple parallel simulations
//! // let results = scheduler.run_parallel::<MyModel>(8, &params);
//! ```

use std::time::Instant;

use rand_core::{Rng, SeedableRng};

use crate::{
    AdaptiveRunControl, Backend, CarloError, Context, FromParams, Metadata, Params, Results,
    RunDecision, RunPhase,
};

/// Configuration for a Monte Carlo run.
#[derive(Debug, Clone)]
pub struct RunConfig {
    /// Number of thermalization sweeps.
    pub thermalization_sweeps: u64,

    /// Number of measurement sweeps.
    pub measurement_sweeps: u64,

    /// Binsize for measurements.
    pub binsize: usize,

    /// Base RNG seed.
    pub base_seed: u64,

    /// Progress report interval (sweeps).
    pub progress_interval: u64,

    /// Checkpoint interval (sweeps). 0 = disabled.
    pub checkpoint_interval: u64,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            thermalization_sweeps: 1000,
            measurement_sweeps: 10000,
            binsize: 100,
            base_seed: 0,
            progress_interval: 1000,
            checkpoint_interval: 0, // disabled by default
        }
    }
}

/// Monte Carlo scheduler. Orchestrates simulation runs.
pub struct Scheduler<B: Backend> {
    backend: B,
    config: RunConfig,
}

impl<B: Backend> Scheduler<B> {
    /// Create a new scheduler with given backend and config.
    pub fn new(backend: B, config: RunConfig) -> Self {
        Self { backend, config }
    }

    /// Run a single simulation task.
    pub fn run_one<MC: FromParams>(&self, params: &Params) -> Results {
        let start = Instant::now();

        // Initialize RNG and context
        let rng = MC::Rng::seed_from_u64(self.config.base_seed);
        let mut ctx =
            Context::new_with_binsize(rng, self.config.thermalization_sweeps, self.config.binsize);

        // Create model from params
        let mut mc =
            MC::from_params(params, &mut ctx.rng).expect("Failed to create model from params");

        // Thermalization phase. Zero-warmup runs enter production directly.
        if self.config.thermalization_sweeps > 0 {
            ctx.enter_phase(RunPhase::Thermalization);
            mc.on_phase_start(RunPhase::Thermalization, &mut ctx);
            for _ in 0..self.config.thermalization_sweeps {
                mc.sweep(&mut ctx);
                ctx.advance_sweep();
            }
            mc.on_phase_end(RunPhase::Thermalization, &mut ctx);
        }

        // Measurement phase. Entering it explicitly freezes adaptive kernels
        // before the first production sweep, including zero-warmup runs.
        ctx.enter_phase(RunPhase::Measurement);
        mc.on_phase_start(RunPhase::Measurement, &mut ctx);
        for sweep in 0..self.config.measurement_sweeps {
            mc.sweep(&mut ctx);
            mc.measure(&mut ctx);
            ctx.advance_sweep();

            // Progress reporting (placeholder)
            if self.config.progress_interval > 0 && sweep % self.config.progress_interval == 0 {
                // Could emit tracing event here
            }
        }
        mc.on_phase_end(RunPhase::Measurement, &mut ctx);
        ctx.enter_phase(RunPhase::Finished);
        mc.on_phase_start(RunPhase::Finished, &mut ctx);

        // Finalize measurements
        let estimates = ctx.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);

        // Set metadata
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: chrono::Utc::now(),
            base_seed: self.config.base_seed,
            thermalization_sweeps: self.config.thermalization_sweeps,
            measurement_sweeps: self.config.measurement_sweeps,
            n_tasks: 1,
        });

        let _duration = start.elapsed();
        // Could log duration here

        results
    }

    /// Run a single simulation whose phase transition and stopping point are
    /// decided by an algorithm-specific controller.
    ///
    /// The fixed-count [`Self::run_one`] path remains the default. This method
    /// is intended for adaptive workflows whose thermalization length is known
    /// only at runtime.
    pub fn run_controlled<MC, C>(
        &self,
        params: &Params,
        mut control: C,
    ) -> Result<Results, CarloError>
    where
        MC: FromParams,
        C: AdaptiveRunControl<MC>,
    {
        let rng = MC::Rng::seed_from_u64(self.config.base_seed);
        let mut ctx =
            Context::new_with_binsize(rng, self.config.thermalization_sweeps, self.config.binsize);
        let mut mc = MC::from_params(params, &mut ctx.rng)?;
        let initial_phase = control.initial_phase();
        if !matches!(
            initial_phase,
            RunPhase::Thermalization | RunPhase::Measurement
        ) {
            return Err(CarloError::InvalidConfig {
                field: "run_control.initial_phase".into(),
                reason: "expected Thermalization or Measurement".into(),
            });
        }

        ctx.enter_phase(initial_phase);
        mc.on_phase_start(initial_phase, &mut ctx);
        let mut thermalization_sweeps = 0_u64;
        let mut measurement_sweeps = 0_u64;

        loop {
            mc.sweep(&mut ctx);
            match ctx.phase() {
                RunPhase::Thermalization => thermalization_sweeps += 1,
                RunPhase::Measurement => {
                    mc.measure(&mut ctx);
                    measurement_sweeps += 1;
                }
                RunPhase::Initialization | RunPhase::Finished => {
                    return Err(CarloError::InvalidConfig {
                        field: "run_control.phase".into(),
                        reason: "controller entered a non-runnable phase".into(),
                    });
                }
            }
            ctx.advance_sweep();

            match (ctx.phase(), control.after_sweep(&mc, &ctx)) {
                (RunPhase::Thermalization, RunDecision::ContinueAdaptation)
                | (RunPhase::Measurement, RunDecision::ContinueProduction) => {}
                (RunPhase::Thermalization, RunDecision::BeginProduction) => {
                    mc.on_phase_end(RunPhase::Thermalization, &mut ctx);
                    ctx.enter_phase(RunPhase::Measurement);
                    mc.on_phase_start(RunPhase::Measurement, &mut ctx);
                }
                (_, RunDecision::Stop) => break,
                (RunPhase::Thermalization, RunDecision::ContinueProduction)
                | (RunPhase::Measurement, RunDecision::ContinueAdaptation)
                | (RunPhase::Measurement, RunDecision::BeginProduction) => {
                    return Err(CarloError::InvalidConfig {
                        field: "run_control.decision".into(),
                        reason: "decision does not match the active phase".into(),
                    });
                }
                (RunPhase::Initialization | RunPhase::Finished, _) => unreachable!(),
            }
        }

        let previous = ctx.phase();
        mc.on_phase_end(previous, &mut ctx);
        ctx.enter_phase(RunPhase::Finished);
        mc.on_phase_start(RunPhase::Finished, &mut ctx);

        let estimates = ctx.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: chrono::Utc::now(),
            base_seed: self.config.base_seed,
            thermalization_sweeps,
            measurement_sweeps,
            n_tasks: 1,
        });
        Ok(results)
    }

    /// Run an adaptive simulation and return both the finalized model state and
    /// ordinary measurement results.
    ///
    /// Density-of-states estimators and other adaptive algorithms often keep
    /// their primary result in the Monte Carlo state rather than in scalar
    /// measurements. This variant preserves that state after scheduler
    /// finalization, including the terminal lifecycle hook.
    pub fn run_controlled_with_state<MC, C>(
        &self,
        params: &Params,
        mut control: C,
    ) -> Result<(MC, Results), CarloError>
    where
        MC: FromParams,
        C: AdaptiveRunControl<MC>,
    {
        let rng = MC::Rng::seed_from_u64(self.config.base_seed);
        let mut ctx =
            Context::new_with_binsize(rng, self.config.thermalization_sweeps, self.config.binsize);
        let mut mc = MC::from_params(params, &mut ctx.rng)?;
        let initial_phase = control.initial_phase();
        if !matches!(
            initial_phase,
            RunPhase::Thermalization | RunPhase::Measurement
        ) {
            return Err(CarloError::InvalidConfig {
                field: "run_control.initial_phase".into(),
                reason: "expected Thermalization or Measurement".into(),
            });
        }

        ctx.enter_phase(initial_phase);
        mc.on_phase_start(initial_phase, &mut ctx);
        let mut thermalization_sweeps = 0_u64;
        let mut measurement_sweeps = 0_u64;

        loop {
            mc.sweep(&mut ctx);
            match ctx.phase() {
                RunPhase::Thermalization => thermalization_sweeps += 1,
                RunPhase::Measurement => {
                    mc.measure(&mut ctx);
                    measurement_sweeps += 1;
                }
                RunPhase::Initialization | RunPhase::Finished => {
                    return Err(CarloError::InvalidConfig {
                        field: "run_control.phase".into(),
                        reason: "controller entered a non-runnable phase".into(),
                    });
                }
            }
            ctx.advance_sweep();

            match (ctx.phase(), control.after_sweep(&mc, &ctx)) {
                (RunPhase::Thermalization, RunDecision::ContinueAdaptation)
                | (RunPhase::Measurement, RunDecision::ContinueProduction) => {}
                (RunPhase::Thermalization, RunDecision::BeginProduction) => {
                    mc.on_phase_end(RunPhase::Thermalization, &mut ctx);
                    ctx.enter_phase(RunPhase::Measurement);
                    mc.on_phase_start(RunPhase::Measurement, &mut ctx);
                }
                (_, RunDecision::Stop) => break,
                (RunPhase::Thermalization, RunDecision::ContinueProduction)
                | (RunPhase::Measurement, RunDecision::ContinueAdaptation)
                | (RunPhase::Measurement, RunDecision::BeginProduction) => {
                    return Err(CarloError::InvalidConfig {
                        field: "run_control.decision".into(),
                        reason: "decision does not match the active phase".into(),
                    });
                }
                (RunPhase::Initialization | RunPhase::Finished, _) => unreachable!(),
            }
        }

        let previous = ctx.phase();
        mc.on_phase_end(previous, &mut ctx);
        ctx.enter_phase(RunPhase::Finished);
        mc.on_phase_start(RunPhase::Finished, &mut ctx);

        let estimates = ctx.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: chrono::Utc::now(),
            base_seed: self.config.base_seed,
            thermalization_sweeps,
            measurement_sweeps,
            n_tasks: 1,
        });
        Ok((mc, results))
    }

    /// Run multiple parallel simulation tasks.
    pub fn run_parallel<MC: FromParams>(&self, n_tasks: usize, params: &Params) -> Vec<Results> {
        use std::sync::Mutex;

        let results: Mutex<Vec<Results>> = Mutex::new(Vec::with_capacity(n_tasks));
        let config = self.config.clone();

        self.backend
            .spawn_tasks(n_tasks, config.base_seed, |task_id, rng| {
                let mut ctx = Context::new_with_binsize(
                    MC::Rng::seed_from_u64(rng.next_u64()),
                    config.thermalization_sweeps,
                    config.binsize,
                );

                let mut mc = MC::from_params(params, &mut ctx.rng)
                    .expect("Failed to create model from params");

                // Thermalization. Zero-warmup tasks enter production directly.
                if config.thermalization_sweeps > 0 {
                    ctx.enter_phase(RunPhase::Thermalization);
                    mc.on_phase_start(RunPhase::Thermalization, &mut ctx);
                    for _ in 0..config.thermalization_sweeps {
                        mc.sweep(&mut ctx);
                        ctx.advance_sweep();
                    }
                    mc.on_phase_end(RunPhase::Thermalization, &mut ctx);
                }

                // Measurement
                ctx.enter_phase(RunPhase::Measurement);
                mc.on_phase_start(RunPhase::Measurement, &mut ctx);
                for _ in 0..config.measurement_sweeps {
                    mc.sweep(&mut ctx);
                    mc.measure(&mut ctx);
                    ctx.advance_sweep();
                }
                mc.on_phase_end(RunPhase::Measurement, &mut ctx);
                ctx.enter_phase(RunPhase::Finished);
                mc.on_phase_start(RunPhase::Finished, &mut ctx);

                let estimates = ctx.finalize_measurements();
                let mut result = Results::from_measurements(&estimates);
                result.set_metadata(Metadata {
                    version: env!("CARGO_PKG_VERSION").to_string(),
                    timestamp: chrono::Utc::now(),
                    base_seed: config.base_seed.wrapping_add(task_id as u64),
                    thermalization_sweeps: config.thermalization_sweeps,
                    measurement_sweeps: config.measurement_sweeps,
                    n_tasks,
                });

                results.lock().expect("results mutex poisoned").push(result);
            });

        self.backend.barrier();
        results.into_inner().expect("results mutex poisoned")
    }
}
