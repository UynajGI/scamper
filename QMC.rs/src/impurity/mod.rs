//! Quantum-impurity solvers built on QMC.rs and Carlo.rs.
//!
//! The production implementations currently provided here are the spin-boson
//! retarded-interaction wormhole solver and the longitudinal continuous-time
//! cluster solver. The module boundary is intentionally wider than either
//! representation so fermionic, bosonic, and Bose-Fermi impurity backends can
//! be added without forcing them into the same configuration or update type.

use thiserror::Error;

pub mod core;
pub mod spin_boson;

/// Impurity-solver construction and runtime errors.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ImpurityError {
    /// A physical or algorithmic parameter is outside its valid domain.
    #[error("invalid parameter `{field}`: {reason}")]
    InvalidParameter {
        /// Parameter name.
        field: String,
        /// Human-readable reason.
        reason: String,
    },

    /// A sampled operator configuration violates worldline invariants.
    #[error("invalid impurity configuration: {0}")]
    InvalidConfiguration(String),

    /// The directed loop exceeded its safety limit without closing.
    #[error("directed loop did not close after {steps} steps (limit {limit})")]
    LoopDidNotClose {
        /// Steps executed.
        steps: usize,
        /// Configured safety limit.
        limit: usize,
    },

    /// A tabulated bath was malformed.
    #[error("invalid tabulated bath: {0}")]
    InvalidBathTable(String),
}

impl ImpurityError {
    /// Convenience constructor for invalid parameters.
    pub fn parameter(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidParameter {
            field: field.into(),
            reason: reason.into(),
        }
    }
}

pub use core::estimator::{
    connected_susceptibility, register_connected_susceptibility, TransverseCorrelationSample,
};
pub use core::kernel::{KernelDirection, PairFlipGauge, SignFreeMetadata, SignFreeReport};
pub use core::local_hilbert::Spin;
pub use core::operators::{
    BasisTransform, PhysicalAxis, SignedAxis, VertexKind, A_IN, A_OUT, B_IN, B_OUT, LEGS_PER_VERTEX,
};
pub use spin_boson::bath::{Bath, BathSample, PowerLawBath, SingleModeBath, TabulatedBath};
pub use spin_boson::cluster::{
    build_segments, measure_cluster_observables, register_cluster_evaluables, ClusterDiagnostics,
    ClusterUpdateReport, ContinuousTimeClusterEngine, LongitudinalClusterObservables,
    LongitudinalSpinBosonClusterQmc, LongitudinalSpinBosonModel, LongitudinalWorldline,
    RetardedKernel, TimeInterval, WorldlineSegment,
};
pub use spin_boson::model::{
    CouplingNormalization, ImpurityModel, ImpurityModelKind, InteractionChannel,
};
pub use spin_boson::observables::{
    correlation_sigma_z, integrated_sigma_z, measure_observables, register_impurity_evaluables,
    ImpurityObservables,
};
pub use spin_boson::occupation::{
    CavityMode, OccupationBasis, OccupationModelKind, OccupationObservables,
    OccupationSpinBosonModel, OccupationWorldlineQmc, OccupationWorldlineSampler, SpinState,
};
pub use spin_boson::wormhole::configuration::{
    EndpointId, Event, LegId, LegSide, Vertex, VertexId, WormholeConfiguration,
};
pub use spin_boson::wormhole::mc::ImpurityQmc;
pub use spin_boson::wormhole::scattering::{
    kind_after_flips, ScatteringChoice, ScatteringDiagnostics, ScatteringPolicy, ScatteringTable,
};
pub use spin_boson::wormhole::updates::{LoopStartPolicy, WormholeEngine, WormholeUpdateStats};
