//! Explicit Monte Carlo and event-time clocks.
//!
//! Sweep counts are scheduler bookkeeping.  Dynamic Monte Carlo kernels may
//! additionally advance attempted-update, accepted-move, or physical event
//! time clocks.  Keeping the units distinct prevents kinetic time from being
//! confused with an arbitrary sweep definition.

use serde::{Deserialize, Serialize};

/// One typed simulation-clock reading.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SimulationClock {
    Sweeps(u64),
    Attempts(u64),
    AcceptedMoves(u64),
    EventTime(f64),
}

impl SimulationClock {
    /// Numeric value of this clock reading.
    #[inline]
    pub const fn value(self) -> f64 {
        match self {
            Self::Sweeps(value) | Self::Attempts(value) | Self::AcceptedMoves(value) => {
                value as f64
            }
            Self::EventTime(value) => value,
        }
    }
}
