//! Algorithm trait and shared utilities.

use crate::lattice::interaction::Hamiltonian;
use crate::lattice::state::System;
use rand::Rng;

/// CMC update phase retained for source compatibility.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimulationPhase {
    Thermalization,
    Measurement,
}

impl SimulationPhase {
    #[inline]
    pub const fn allows_adaptation(self) -> bool {
        matches!(self, Self::Thermalization)
    }

    /// Map Carlo.rs's full lifecycle onto the two phases relevant to updates.
    /// Initialization and finished states use the frozen production kernel if
    /// a custom driver invokes a sweep outside the normal scheduler contract.
    #[inline]
    pub const fn from_run_phase(phase: carlo_rs::RunPhase) -> Self {
        match phase {
            carlo_rs::RunPhase::Thermalization => Self::Thermalization,
            carlo_rs::RunPhase::Initialization
            | carlo_rs::RunPhase::Measurement
            | carlo_rs::RunPhase::Finished => Self::Measurement,
        }
    }
}

/// One update policy. Carlo.rs owns scheduling; CMC.rs owns state transitions.
pub trait Algorithm<H: Hamiltonian>: Send {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    );

    /// Direct/manual sweeps default to the frozen production kernel.
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        self.sweep_with_phase(system, model, rng, SimulationPhase::Measurement);
    }

    fn name(&self) -> &'static str {
        "Unknown"
    }

    /// Lifecycle hook used by adaptive kernels to mark a completed run terminal.
    fn finish_run(&mut self) {}
}

/// Marker for fixed-parameter lattice kernels whose replica-exchange weight is `-βE`.
///
/// Generalized-ensemble kernels deliberately do not implement this marker:
/// changing `β` alone does not describe their target distribution.
pub trait CanonicalLatticeKernel {}

/// Validate and clamp a bond activation probability.
pub fn checked_probability(value: f64, algorithm: &str) -> f64 {
    assert!(
        value.is_finite() && (0.0..=1.0).contains(&value),
        "{algorithm} model returned invalid bond probability {value}"
    );
    value
}
