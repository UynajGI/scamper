//! Runtime context for Monte Carlo simulations.
//!
//! [`Context`] holds the state during simulation execution:
//! - Random number generator
//! - Measurement accumulator
//! - Sweep counter and thermalization tracking
//!
//! # Usage
//!
//! The context is passed to [`MonteCarlo::sweep()`] and [`MonteCarlo::measure()`]
//! methods, providing access to RNG and measurement recording.
//!
//! # Checkpointing
//!
//! With the `hdf5` feature, context can be saved/loaded from HDF5 files
//! for simulation restart capability.

use rand_core::Rng;
use rand_core::SeedableRng;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::{ComplexEstimate, Estimate, Measurements, RunPhase, SimulationClock};

/// Checkpoint state for serialization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContextCheckpoint {
    /// Sweeps completed at checkpoint time.
    pub sweep_count: u64,
    /// Configured thermalization sweep count.
    pub thermalization_sweeps: u64,
    /// Whether the warmup boundary has been crossed.
    pub thermalized: bool,
    /// Explicit lifecycle phase (defaults to `Initialization` for legacy checkpoints).
    #[serde(default)]
    pub phase: RunPhase,
    /// Algorithm-reported attempted elementary transitions.
    #[serde(default)]
    pub attempted_updates: u64,
    /// Algorithm-reported accepted or executed transitions.
    #[serde(default)]
    pub accepted_moves: u64,
    /// Physical/event time advanced by kinetic Monte Carlo kernels.
    #[serde(default)]
    pub event_time: f64,
}

/// Runtime context for Monte Carlo simulation.
pub struct Context<R: Rng + SeedableRng> {
    /// Random number generator (public access).
    pub rng: R,

    /// Measurement collector.
    measurements: Measurements,

    /// Optional namespace prepended to observable names. This is transient
    /// scheduler state used by algorithms such as parallel tempering.
    measurement_namespace: Option<String>,

    /// Current sweep count.
    sweep_count: u64,

    /// Thermalization sweeps threshold.
    thermalization_sweeps: u64,

    /// Whether the warmup boundary has been crossed.
    thermalized: bool,

    /// Explicit scheduler-owned lifecycle phase.
    phase: RunPhase,

    /// Algorithm-reported attempted elementary transitions.
    attempted_updates: u64,

    /// Algorithm-reported accepted or executed transitions.
    accepted_moves: u64,

    /// Physical/event time advanced by kinetic Monte Carlo kernels.
    event_time: f64,
}

impl<R: Rng + SeedableRng> Context<R> {
    /// Create new context.
    pub fn new(rng: R, thermalization_sweeps: u64) -> Self {
        Self {
            rng,
            measurements: Measurements::new(100), // default binsize
            measurement_namespace: None,
            sweep_count: 0,
            thermalization_sweeps,
            thermalized: false,
            phase: RunPhase::Initialization,
            attempted_updates: 0,
            accepted_moves: 0,
            event_time: 0.0,
        }
    }

    /// Create new context with custom binsize.
    pub fn new_with_binsize(rng: R, thermalization_sweeps: u64, binsize: usize) -> Self {
        Self {
            rng,
            measurements: Measurements::new(binsize),
            measurement_namespace: None,
            sweep_count: 0,
            thermalization_sweeps,
            thermalized: false,
            phase: RunPhase::Initialization,
            attempted_updates: 0,
            accepted_moves: 0,
            event_time: 0.0,
        }
    }

    /// Finalize measurements and return estimates.
    pub fn finalize_measurements(&self) -> HashMap<String, Estimate> {
        self.measurements.finalize()
    }

    /// Finalize complex measurements and return estimates.
    pub fn finalize_complex_measurements(&self) -> HashMap<String, ComplexEstimate> {
        self.measurements.finalize_complex()
    }

    fn qualified_measurement_name(&self, name: &str) -> String {
        match &self.measurement_namespace {
            Some(namespace) if !namespace.is_empty() => format!("{namespace}/{name}"),
            _ => name.to_string(),
        }
    }

    /// Set a transient namespace for subsequently recorded observables.
    ///
    /// Schedulers should clear the namespace after the scoped operation. It is
    /// deliberately not persisted in checkpoints.
    pub fn set_measurement_namespace(&mut self, namespace: Option<String>) {
        self.measurement_namespace = namespace;
    }

    /// Return the current transient measurement namespace.
    pub fn measurement_namespace(&self) -> Option<&str> {
        self.measurement_namespace.as_deref()
    }

    /// Record an observable sample (scalar).
    pub fn measure(&mut self, name: &str, value: f64) {
        let name = self.qualified_measurement_name(name);
        self.measurements.add_sample(&name, value);
    }

    /// Record an array observable sample.
    /// The shape is determined by the first call for each observable name.
    pub fn measure_array(&mut self, name: &str, values: &[f64]) {
        let name = self.qualified_measurement_name(name);
        self.measurements.add_sample_array(&name, values);
    }

    /// Record a complex observable sample.
    /// Real and imaginary parts are accumulated separately.
    pub fn measure_complex(&mut self, name: &str, re: f64, im: f64) {
        let name = self.qualified_measurement_name(name);
        self.measurements.add_sample_complex(&name, re, im);
    }

    /// Register a scalar observable with custom binsize.
    pub fn register_observable(&mut self, name: &str, binsize: usize) {
        let name = self.qualified_measurement_name(name);
        self.measurements.register(&name, binsize);
    }

    /// Register an array observable with custom binsize and shape.
    pub fn register_observable_with_shape(&mut self, name: &str, binsize: usize, shape: &[usize]) {
        let name = self.qualified_measurement_name(name);
        self.measurements.register_array(&name, binsize, shape);
    }

    /// Check if thermalized.
    pub fn is_thermalized(&self) -> bool {
        match self.phase {
            RunPhase::Thermalization => false,
            RunPhase::Measurement | RunPhase::Finished => true,
            RunPhase::Initialization => self.thermalized,
        }
    }

    /// Current explicit run phase.
    #[inline]
    pub const fn phase(&self) -> RunPhase {
        self.phase
    }

    /// Enter a scheduler phase.
    ///
    /// This is public so custom schedulers and [`crate::Run`] can preserve the
    /// same lifecycle semantics as [`crate::Scheduler`].
    pub fn enter_phase(&mut self, phase: RunPhase) {
        self.phase = phase;
        if matches!(phase, RunPhase::Measurement | RunPhase::Finished) {
            self.thermalized = true;
        }
    }

    /// Get sweep count.
    pub fn sweep_count(&self) -> u64 {
        self.sweep_count
    }

    /// Advance sweep counter.
    pub fn advance_sweep(&mut self) {
        self.sweep_count += 1;
        if self.sweep_count >= self.thermalization_sweeps {
            self.thermalized = true;
        }
    }

    /// Add attempted elementary transitions to the explicit attempt clock.
    #[inline]
    pub fn record_attempts(&mut self, count: u64) {
        self.attempted_updates = self.attempted_updates.saturating_add(count);
    }

    /// Add accepted/executed transitions to the explicit accepted-move clock.
    #[inline]
    pub fn record_accepted_moves(&mut self, count: u64) {
        self.accepted_moves = self.accepted_moves.saturating_add(count);
    }

    /// Advance physical/event time.  The increment must be finite and non-negative.
    #[inline]
    pub fn advance_event_time(&mut self, delta: f64) {
        assert!(
            delta.is_finite() && delta >= 0.0,
            "event-time increment must be finite and non-negative"
        );
        self.event_time += delta;
        assert!(self.event_time.is_finite(), "event-time clock overflowed");
    }

    /// Total number of attempted elementary transitions.
    #[inline]
    pub const fn attempted_updates(&self) -> u64 {
        self.attempted_updates
    }

    /// Total number of accepted or executed transitions.
    #[inline]
    pub const fn accepted_moves(&self) -> u64 {
        self.accepted_moves
    }

    /// Accumulated physical/event time.
    #[inline]
    pub const fn event_time(&self) -> f64 {
        self.event_time
    }

    /// Return all clock readings without conflating their units.
    pub fn simulation_clocks(&self) -> [SimulationClock; 4] {
        [
            SimulationClock::Sweeps(self.sweep_count),
            SimulationClock::Attempts(self.attempted_updates),
            SimulationClock::AcceptedMoves(self.accepted_moves),
            SimulationClock::EventTime(self.event_time),
        ]
    }

    /// Snapshot the current context state for serialization.
    pub fn checkpoint_state(&self) -> ContextCheckpoint {
        ContextCheckpoint {
            sweep_count: self.sweep_count,
            thermalization_sweeps: self.thermalization_sweeps,
            thermalized: self.thermalized,
            phase: self.phase,
            attempted_updates: self.attempted_updates,
            accepted_moves: self.accepted_moves,
            event_time: self.event_time,
        }
    }

    /// Reconstruct a context from checkpoint state, inferring the lifecycle
    /// phase from counters when the checkpoint predates explicit phases.
    pub fn restore_from_checkpoint(checkpoint: ContextCheckpoint, rng: R, binsize: usize) -> Self {
        // Old serialized checkpoints do not contain `phase` and deserialize it
        // as `Initialization`. Infer the active phase from their counters while
        // preserving a genuine pre-run initialization checkpoint.
        let phase = if checkpoint.phase == RunPhase::Initialization {
            if checkpoint.thermalized || checkpoint.sweep_count >= checkpoint.thermalization_sweeps
            {
                RunPhase::Measurement
            } else if checkpoint.sweep_count > 0 {
                RunPhase::Thermalization
            } else {
                RunPhase::Initialization
            }
        } else {
            checkpoint.phase
        };
        Self {
            rng,
            measurements: Measurements::new(binsize),
            measurement_namespace: None,
            sweep_count: checkpoint.sweep_count,
            thermalization_sweeps: checkpoint.thermalization_sweeps,
            thermalized: checkpoint.thermalized
                || matches!(phase, RunPhase::Measurement | RunPhase::Finished),
            phase,
            attempted_updates: checkpoint.attempted_updates,
            accepted_moves: checkpoint.accepted_moves,
            event_time: checkpoint.event_time,
        }
    }

    /// Get thermalization sweeps setting.
    pub fn thermalization_sweeps(&self) -> u64 {
        self.thermalization_sweeps
    }
}

#[cfg(feature = "hdf5")]
use hdf5::Group;

#[cfg(feature = "hdf5")]
use crate::RngCheckpointHdf5;

#[cfg(feature = "hdf5")]
impl<R: Rng + SeedableRng + RngCheckpointHdf5> Context<R> {
    /// Write context state to HDF5 group (includes RNG state).
    pub fn write_checkpoint_hdf5(&self, group: &mut Group) -> Result<(), crate::CarloError> {
        // Write sweep state
        group
            .create_dataset_simple("sweep_count", &[1], &self.sweep_count.to_ne_bytes())
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write sweep_count: {}", e),
            })?;
        group
            .create_dataset_simple(
                "thermalization_sweeps",
                &[1],
                &self.thermalization_sweeps.to_ne_bytes(),
            )
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write thermalization_sweeps: {}", e),
            })?;

        // Write RNG state
        let rng_group =
            group
                .create_group("rng")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot create rng group: {}", e),
                })?;
        self.rng.write_checkpoint(&mut rng_group)?;

        // Write measurements
        let meas_group =
            group
                .create_group("measurements")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot create measurements group: {}", e),
                })?;
        self.measurements.write_checkpoint_hdf5(&mut meas_group)?;

        Ok(())
    }

    /// Read context from HDF5 group (includes RNG state).
    pub fn read_checkpoint_hdf5_full(
        group: &Group,
        binsize: usize,
    ) -> Result<Self, crate::CarloError> {
        let sweep_count = group
            .dataset("sweep_count")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read sweep_count: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse sweep_count: {}", e),
            })?[0];

        let thermalization_sweeps = group
            .dataset("thermalization_sweeps")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read thermalization_sweeps: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse thermalization_sweeps: {}", e),
            })?[0];

        // Read RNG
        let rng_group = group
            .group("rng")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot open rng group: {}", e),
            })?;
        let rng = R::read_checkpoint(&rng_group)?;

        // Read measurements
        let meas_group =
            group
                .group("measurements")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot open measurements group: {}", e),
                })?;
        let measurements = crate::Measurements::read_checkpoint_hdf5(&meas_group)?;

        Ok(Self {
            rng,
            measurements,
            measurement_namespace: None,
            sweep_count,
            thermalization_sweeps,
            thermalized: sweep_count >= thermalization_sweeps,
            phase: if sweep_count >= thermalization_sweeps {
                RunPhase::Measurement
            } else {
                RunPhase::Thermalization
            },
            attempted_updates: 0,
            accepted_moves: 0,
            event_time: 0.0,
        })
    }
}
