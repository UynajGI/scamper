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
//! - [`impurity`] — retarded-interaction wormhole QMC and a continuous-time
//!   cluster solver for longitudinal spin-boson impurities.
//!
//! The [`variational`] family (continuum VMC, layer L0) shares the Carlo.rs
//! hosting conventions and grows toward optimizers, DMC and NQS behind the
//! `WaveFunction` trait.
//!
//! Discrete-time prototypes and the old chain-specific Heisenberg adapter have
//! been removed. Lattice geometry is now data, not an algorithm type.

pub mod algorithm;
pub mod graph;
pub mod impurity;
pub mod lattice;
pub mod local_space;
pub mod variational;

pub use algorithm::{QmcKernel, UpdateSchedule};
pub use graph::{CsrGraph, Edge, EdgeSpec, GraphError, Neighbor};
pub use impurity::{
    // Observable helpers
    connected_susceptibility,
    correlation_sigma_z,
    integrated_sigma_z,
    measure_cluster_observables,
    measure_observables,
    register_cluster_evaluables,
    register_connected_susceptibility,
    register_impurity_evaluables,
    // Core types
    BasisTransform,
    Bath,
    BathSample,
    // Occupation solver
    CavityMode,
    ClusterDiagnostics,
    ClusterUpdateReport,
    ContinuousTimeClusterEngine,
    CouplingNormalization,
    ImpurityError,
    ImpurityModel,
    ImpurityModelKind,
    ImpurityObservables,
    ImpurityQmc,
    InteractionChannel,
    KernelDirection,
    LongitudinalClusterObservables,
    LongitudinalSpinBosonClusterQmc,
    LongitudinalSpinBosonModel,
    LongitudinalWorldline,
    LoopStartPolicy,
    OccupationBasis,
    OccupationModelKind,
    OccupationObservables,
    OccupationSpinBosonModel,
    OccupationWorldlineQmc,
    OccupationWorldlineSampler,
    PairFlipGauge,
    PhysicalAxis,
    PowerLawBath,
    RetardedKernel,
    SignFreeMetadata,
    SignFreeReport,
    SingleModeBath,
    Spin,
    SpinState,
    TabulatedBath,
    TransverseCorrelationSample,
    VertexKind,
    WormholeConfiguration,
    WormholeEngine,
    WormholeUpdateStats,
    A_IN,
    A_OUT,
    B_IN,
    B_OUT,
    LEGS_PER_VERTEX,
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
pub use variational::{
    harmonic_closed_shell_electrons, harmonic_closed_shell_energy, harmonic_trap_orbitals,
    local_energy, Backflow, BlockStats, ContinuumHamiltonian, DeltaLog, DmcKernel, DmcStats,
    GaussianTrap, GradBuffer, GtoOrbital, HarmonicJastrow, HarmonicTrap, LinearMethod, LocalEnergy,
    McMillanJastrow, Optimizer, PairPotential, ParamGradBuffer, Point, Positions, Product,
    ReferenceSample, SlaterDeterminant, StochasticReconfiguration, VarianceMinimization,
    VarianceMinimizationResult, VarianceObjective, VariationalError, VmcKernel, VmcStats, Walker,
    WaveFunction, WaveFunctionParams, DIM, DMC_CHECKPOINT_FORMAT, VMC_CHECKPOINT_FORMAT,
};
