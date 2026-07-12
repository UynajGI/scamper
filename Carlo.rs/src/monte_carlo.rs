//! Core trait for Monte Carlo algorithms.
//!
//! The [`MonteCarlo`] trait is the central abstraction for implementing
//! Monte Carlo simulations in this framework. Users implement the `sweep()`
//! method to update configurations and optionally `measure()` for observables.
//!
//! # Required Methods
//!
//! - [`MonteCarlo::sweep()`]: Execute one Monte Carlo sweep (update configuration)
//! - [`MonteCarlo::Rng`]: Specify the RNG type (must implement `Rng + SeedableRng + Send`)
//!
//! # Optional Methods
//!
//! - [`MonteCarlo::measure()`]: Measure observables (called after thermalization)
//! - [`MonteCarlo::name()`]: Return algorithm name (for logging/metadata)
//! - [`MonteCarloCheckpoint`]: Save/load state to HDF5 (requires `hdf5` feature)

use rand_core::Rng;
use rand_core::SeedableRng;

#[cfg(feature = "mpi")]
use mpi::topology::SimpleCommunicator;

use crate::CarloError;
use crate::Context;
use crate::Params;
use crate::RunPhase;

/// Core trait for Monte Carlo algorithms.
/// Users implement `sweep()` and optionally override other methods.
pub trait MonteCarlo: Sized {
    /// Execute one sweep (update configuration).
    /// Users may call `ctx.measure()` inside sweep.
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);

    /// RNG type used by this model.
    type Rng: Rng + SeedableRng + Send;

    /// Optional: measure observables (default: empty).
    fn measure(&mut self, _ctx: &mut Context<Self::Rng>) {}

    /// Lifecycle hook called after the context enters a new run phase.
    ///
    /// Adaptive algorithms can use this to reset warmup statistics or freeze
    /// transition-kernel parameters at the start of production.
    fn on_phase_start(&mut self, _phase: RunPhase, _ctx: &mut Context<Self::Rng>) {}

    /// Lifecycle hook called immediately before leaving a run phase.
    fn on_phase_end(&mut self, _phase: RunPhase, _ctx: &mut Context<Self::Rng>) {}

    /// Optional: sweep with MPI communicator for multi-rank coordination.
    /// Default: delegate to [`sweep`](MonteCarlo::sweep).
    #[cfg(feature = "mpi")]
    fn sweep_with_comm(
        &mut self,
        ctx: &mut Context<Self::Rng>,
        _comm: &mpi::topology::SimpleCommunicator,
    ) {
        self.sweep(ctx);
    }

    /// Optional: measure with MPI communicator for multi-rank coordination.
    /// Default: delegate to [`measure`](MonteCarlo::measure).
    #[cfg(feature = "mpi")]
    fn measure_with_comm(
        &mut self,
        ctx: &mut Context<Self::Rng>,
        _comm: &mpi::topology::SimpleCommunicator,
    ) {
        self.measure(ctx);
    }

    /// Optional: save state to HDF5 (default: empty).
    #[cfg(feature = "hdf5")]
    fn save(&self, _out: &mut hdf5::Group) {}

    /// Optional: load state from HDF5 (default: empty).
    #[cfg(feature = "hdf5")]
    fn load(&mut self, _in: &hdf5::Group) {}

    /// Optional: algorithm name (default: "UnnamedMC").
    fn name(&self) -> &'static str {
        "UnnamedMC"
    }
}

/// Trait for constructing models from parameters.
pub trait FromParams: MonteCarlo {
    /// Construct model from params.
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError>;

    /// Optional: validate parameters before construction.
    fn validate_params(_params: &Params) -> Result<(), CarloError> {
        Ok(())
    }
}

/// Extension trait for MonteCarlo with optional lifecycle methods.
pub trait MonteCarloExt: MonteCarlo {
    /// Initialize the simulation. Called once at start.
    fn init(&mut self, _ctx: &mut Context<Self::Rng>, _params: &Params) {}

    /// Register derived observables for post-processing.
    fn register_evaluables(
        _mc_type: std::marker::PhantomData<Self>,
        _evaluator: &mut crate::evaluable::Evaluator,
        _params: &Params,
    ) {
    }
}

// Blanket implementation for all MonteCarlo
impl<MC: MonteCarlo> MonteCarloExt for MC {}

/// Checkpoint support for MonteCarlo implementations.
#[cfg(feature = "hdf5")]
pub trait MonteCarloCheckpoint: MonteCarlo {
    /// Write simulation state to HDF5 group.
    fn write_checkpoint(&self, _group: &mut hdf5::Group) -> Result<(), CarloError> {
        // Default: no state to save
        Ok(())
    }

    /// Read simulation state from HDF5 group.
    fn read_checkpoint(&mut self, _group: &hdf5::Group) -> Result<(), CarloError> {
        // Default: no state to load
        Ok(())
    }
}

#[cfg(feature = "hdf5")]
impl<MC: MonteCarlo> MonteCarloCheckpoint for MC {}
