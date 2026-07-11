//! Reusable foundations for quantum Monte Carlo engines.
//!
//! Carlo.rs owns run orchestration, random-number seeding, measurement
//! accumulation, parallel execution, and result analysis. This module defines
//! the smaller physics-side contract shared by QMC algorithms: advance a
//! representation-specific configuration, validate it, and expose update
//! diagnostics.

use rand::Rng;

/// A physics update kernel independent of Carlo.rs scheduling.
pub trait QmcKernel<C, R>
where
    R: Rng + ?Sized,
{
    /// Runtime or configuration error.
    type Error;
    /// Per-kernel diagnostics accumulated across updates.
    type Diagnostics;

    /// Advance `configuration` by one algorithmic sweep.
    fn sweep(&mut self, configuration: &mut C, rng: &mut R) -> Result<(), Self::Error>;

    /// Validate representation-specific invariants.
    fn validate(&self, configuration: &C) -> Result<(), Self::Error>;

    /// Return immutable update diagnostics.
    fn diagnostics(&self) -> &Self::Diagnostics;
}

/// Fixed work assigned to one continuous-time QMC sweep.
///
/// Keeping the work count fixed during measurement is important: adapting the
/// number of elementary moves to the current expansion order changes the
/// amount of Markov time assigned to each state and biases sweep-sampled
/// observables. Adaptation may be performed during thermalization, then this
/// schedule must be frozen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UpdateSchedule {
    /// Number of diagonal add/remove proposals.
    pub diagonal_proposals: usize,
    /// Number of directed loops.
    pub directed_loops: usize,
    /// Safety factor multiplying the number of vertex legs before declaring a
    /// loop non-closing.
    pub max_loop_steps_factor: usize,
}

impl UpdateSchedule {
    /// Construct a validated schedule with all counts clamped to at least one.
    pub fn new(
        diagonal_proposals: usize,
        directed_loops: usize,
        max_loop_steps_factor: usize,
    ) -> Self {
        Self {
            diagonal_proposals: diagonal_proposals.max(1),
            directed_loops: directed_loops.max(1),
            max_loop_steps_factor: max_loop_steps_factor.max(1),
        }
    }
}

impl Default for UpdateSchedule {
    fn default() -> Self {
        Self::new(1, 1, 16)
    }
}
