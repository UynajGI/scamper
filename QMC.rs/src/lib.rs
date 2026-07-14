//! # QMC.rs — reusable quantum Monte Carlo algorithms
//!
//! QMC.rs is the quantum-physics layer built on Carlo.rs. Carlo.rs owns run
//! scheduling, random seeds, accumulation, error analysis, checkpoint
//! orchestration, and parallel execution. QMC.rs owns representations, sparse
//! operator catalogs, update kernels, and estimators.
//!
//! ## Current production backends
//!
//! - [`lattice`] — continuous-time interaction-expansion directed-loop QMC on
//!   arbitrary CSR adjacency graphs, with arbitrary quantum spin `S`.
//! - [`impurity`] — continuous-time retarded-interaction wormhole QMC for
//!   quantum impurity models.
//!
//! Discrete-time prototypes and the old chain-specific Heisenberg adapter have
//! been removed. Lattice geometry is now data, not an algorithm type.

pub mod algorithm;
pub mod graph;
pub mod impurity;
pub mod lattice;
pub mod local_space;

pub use algorithm::{QmcKernel, UpdateSchedule};
pub use graph::{CsrGraph, Edge, EdgeSpec, GraphError, Neighbor};
pub use impurity::{
    Bath, ImpurityModel, ImpurityModelKind, ImpurityQmc, PowerLawBath, SingleModeBath,
    TabulatedBath, WormholeConfiguration, WormholeEngine,
};
pub use lattice::{
    ContinuousLatticeEngine, EdgeCoupling, GaugePolicy, LatticeConfiguration, LatticeObservables,
    LatticeQmcError, LatticeSpinQmc, LatticeUpdateStats, OperatorTerm, PositiveOperatorModel,
    ScatteringPolicy, SiteCoupling, SpinLatticeModel, SpinModelBuilder, TermLocation,
    Vertex as LatticeVertex, VertexKind as LatticeVertexKind, WorldlineIndex,
};
pub use local_space::{
    BasisState, LocalHilbertSpace, LocalSpaceError, ParticleStatistics, SpinSpace,
};
