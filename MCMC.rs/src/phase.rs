use serde::{Deserialize, Serialize};

/// Statistical sampling phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SamplingPhase {
    /// Tuning is permitted and draws are not posterior samples.
    Warmup,
    /// Tuning is frozen and draws are retained.
    Sampling,
}
