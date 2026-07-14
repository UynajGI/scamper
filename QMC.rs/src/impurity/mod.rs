//! Continuous-time wormhole QMC for quantum impurity models.
//!
//! # Representation
//!
//! Quadratic bosons are integrated out exactly.  The remaining partition
//! function is sampled as a stochastic expansion in retarded four-leg spin
//! vertices `nu = (interaction, kind, omega, tau, tau')`.  A diagonal update
//! changes the expansion order and samples the bath variables.  A directed
//! loop changes diagonal and off-diagonal vertex kinds; crossing from one
//! endpoint to the other is the nonlocal wormhole move.
//!
//! # Included model catalogs
//!
//! - Jaynes-Cummings (directed `S_+ D S_-` vertex)
//! - directed rotating/counter-rotating (RW-CRW) impurity
//! - U(1)-symmetric XXZ impurity
//! - fully anisotropic XYZ impurity
//! - original impurity / single-mode Rabi after a spin-axis rotation
//!
//! # Included bath samplers
//!
//! - single mode
//! - sharp-cutoff power law
//! - arbitrary positive discrete/tabulated spectrum
//!
//! [`ImpurityQmc`] is the Carlo.rs-facing entry point.  Lower-level users can
//! combine [`ImpurityModel`], [`WormholeConfiguration`], and
//! [`WormholeEngine`] directly.

pub mod bath;
pub mod configuration;
pub mod error;
pub mod mc;
pub mod model;
pub mod observables;
pub mod scattering;
pub mod updates;
pub mod vertex;

pub use bath::{Bath, BathSample, KernelDirection, PowerLawBath, SingleModeBath, TabulatedBath};
pub use configuration::WormholeConfiguration;
pub use error::ImpurityError;
pub use mc::ImpurityQmc;
pub use model::{CouplingNormalization, ImpurityModel, ImpurityModelKind, InteractionChannel};
pub use observables::{
    correlation_sigma_z, integrated_sigma_z, measure_observables, ImpurityObservables,
};
pub use scattering::{
    kind_after_flips, ScatteringChoice, ScatteringDiagnostics, ScatteringPolicy, ScatteringTable,
};
pub use updates::{LoopStartPolicy, WormholeEngine, WormholeUpdateStats};
pub use vertex::{
    EndpointId, Event, LegId, LegSide, Spin, Vertex, VertexId, VertexKind, A_IN, A_OUT, B_IN, B_OUT,
};
