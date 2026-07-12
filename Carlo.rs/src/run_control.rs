//! Optional algorithm-driven run control for adaptive Monte Carlo workflows.
//!
//! The fixed-count [`crate::Scheduler::run_one`] path remains the default.
//! This module adds a small control protocol for algorithms whose warmup ends
//! on a convergence condition rather than a predetermined sweep count.

use crate::{Context, MonteCarlo, RunPhase};

/// Decision returned after one completed sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunDecision {
    /// Remain in thermalization/adaptation.
    ContinueAdaptation,
    /// Leave thermalization and enter frozen production before the next sweep.
    BeginProduction,
    /// Remain in production/measurement.
    ContinueProduction,
    /// Finish the run after the completed sweep.
    Stop,
}

/// Algorithm-specific controller for adaptive run length and phase changes.
///
/// The controller observes the accepted Monte Carlo state and context after
/// every completed sweep. It does not perform updates itself, so physics stays
/// in the `MonteCarlo` implementation and execution stays in Carlo.rs.
pub trait AdaptiveRunControl<MC: MonteCarlo> {
    /// Initial active phase. Only thermalization and measurement are valid.
    fn initial_phase(&self) -> RunPhase {
        RunPhase::Thermalization
    }

    /// Decide what should happen before the next sweep.
    fn after_sweep(&mut self, mc: &MC, ctx: &Context<MC::Rng>) -> RunDecision;
}
