//! Explicit Monte Carlo run phases shared by schedulers and update kernels.

use serde::{Deserialize, Serialize};

/// Current lifecycle phase of a Monte Carlo run.
///
/// The phase is stored in [`crate::Context`] so algorithms do not need to infer
/// warmup state from sweep counters.  This is especially important for
/// adaptive kernels: adaptation is permitted only in [`RunPhase::Thermalization`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum RunPhase {
    /// Model construction or restored state before the first scheduled sweep.
    #[default]
    Initialization,
    /// Warmup sweeps. Measurements are normally not accumulated here.
    Thermalization,
    /// Production sweeps with a frozen transition kernel.
    Measurement,
    /// Terminal phase after the scheduler has finalized the run.
    Finished,
}

impl RunPhase {
    /// Whether adaptive transition-kernel parameters may be changed.
    #[inline]
    pub const fn allows_adaptation(self) -> bool {
        matches!(self, Self::Thermalization)
    }

    /// Whether measurements should be accumulated in this phase.
    #[inline]
    pub const fn collects_measurements(self) -> bool {
        matches!(self, Self::Measurement)
    }
}
