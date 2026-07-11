//! Continuous-time lattice QMC on arbitrary adjacency graphs.
//!
//! The runtime pipeline is:
//!
//! `CsrGraph -> LocalHilbertSpace -> sparse OperatorTerm catalog ->`
//! `continuous-time configuration -> diagonal updates + directed loops`.
//!
//! The production implementation supports arbitrary quantum spin `S`,
//! Heisenberg/XY/XXZ/XYZ exchange, transverse and longitudinal fields, and
//! single-ion anisotropy. Fermionic local spaces are reserved at the trait
//! boundary but require a future sign/determinant backend.

pub mod configuration;
pub mod error;
pub mod mc;
pub mod model;
pub mod observables;
pub mod scattering;
pub mod updates;
pub mod vertex;

pub use configuration::{LatticeConfiguration, WorldlineIndex};
pub use error::LatticeQmcError;
pub use mc::LatticeSpinQmc;
pub use model::{
    EdgeCoupling, GaugePolicy, OperatorTerm, PositiveOperatorModel, SiteCoupling, SpinLatticeModel,
    SpinModelBuilder, TermLocation,
};
pub use observables::{measure_observables, site_time_averaged_magnetizations, LatticeObservables};
pub use scattering::{ScatteringChoice, ScatteringDiagnostics, ScatteringPolicy, ScatteringTable};
pub use updates::{ContinuousLatticeEngine, LatticeUpdateStats};
pub use vertex::{Event, Vertex, VertexKind};
