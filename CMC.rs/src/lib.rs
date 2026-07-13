//! CMC.rs — reusable classical lattice and particle Monte Carlo kernels on top of Carlo.rs.
//!
//! Carlo.rs owns execution concerns (RNG contexts, thermalization/measurement
//! scheduling, parallel backends, accumulation and parallel tempering). CMC.rs
//! owns graph/particle physical models, state transitions and observables.
//!
//! # Quick start
//!
//! ```ignore
//! use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
//! use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};
//!
//! type Simulation = ClassicalMC<IsingModel, MetropolisCore>;
//! let mut params = Params::new();
//! params.set("Lx", 16);
//! params.set("Ly", 16);
//! params.set("beta", 0.44);
//! let results = Scheduler::new(RayonBackend::new(1), RunConfig::default())
//!     .run_one::<Simulation>(&params);
//! ```
//!
//! # Module structure
//!
//! | Directory | Purpose |
//! |-----------|---------|
//! | [`core`] | Move types, trial evaluation, ensemble, acceptance rules, visit schedules |
//! | [`lattice`] | Graph topology, physical state, Hamiltonian traits, models, proposals |
//! | [`algorithms`] | Update kernels (Metropolis, Wolff, SW, heat bath, microcanonical, hybrid) |
//! | [`observables`] | Measurement traits and built-in lattice observables |
//! | [`particle`] | Periodic particles, Lennard-Jones potentials, molecular moves and NVT/NPT/μVT kernels |
//!
//! # Extensibility
//!
//! * [`CsrLattice`] is an arbitrary weighted undirected multigraph.
//! * [`Hamiltonian`] defines physical onsite/bond energy.
//! * Capability traits opt models into Metropolis, cluster, heat-bath or
//!   over-relaxation kernels without a monolithic model interface.
//! * [`ClassicalMC`] composes model + kernel + observable set into Carlo.rs's
//!   [`carlo_rs::MonteCarlo`] trait.

// ── Module tree ──────────────────────────────────────────────

// Top-level adapter modules
pub mod classical_mc;
pub mod multi_spin;
pub mod postprocess;
pub mod statistics;

// Hierarchical modules

pub mod algorithms;
pub mod audit;
pub mod core;
pub mod generalized;
pub mod lattice;
pub mod observables;
pub mod particle;

// ── Flat public re-exports (backward-compatible) ─────────────

pub use algorithms::{
    Algorithm, ContinuousHeatBathCore, HeatBathCore, HybridCore, MetropolisCore,
    MicrocanonicalCore, SWCore, SimulationPhase, WolffCore,
};
pub use audit::{
    audit_lattice_cache, audit_macrostate_bin, audit_particle_cache, automatic_cache_audit_enabled,
    effective_cache_audit_interval, should_audit_cache, DEFAULT_CACHE_AUDIT_INTERVAL,
};
pub use classical_mc::{build_lattice_from_params, ClassicalMC, FromHamiltonianParams};
pub use core::acceptance::{AcceptanceRule, MetropolisHastingsAcceptance};
pub use core::cache::{BatchEnergyPatch, BatchEnergyWorkspace, EnergyPatch};
pub use core::ensemble::{
    CanonicalEnsemble, Ensemble, GrandCanonical, IsothermalIsobaric, ThermodynamicDelta,
};
pub use core::r#move::{BatchSpinMove, SiteSpinMove, Spin};
pub use core::trial::{metropolis_hastings_step, ProposedMove, TrialEvaluator, TrialOutcome};
pub use core::visit::{SiteOrder, VisitSchedule};
pub use generalized::{
    canonical_reweight, enumerate_ising_density_of_states, BinnedAxis, CanonicalReweighting,
    DiscreteAxis, EnergyBiasCore, EnergyMacrostate, ExactIsingDensityOfStates, FixedBias,
    GeneralizedError, HarmonicUmbrellaBias, Histogram, IsingWangLandau, LogBias,
    LogDensityOfStates, Macrostate, MacrostateAxis, MagnetizationMacrostate, MulticanonicalBias,
    ParticleNumberMacrostate, WangLandauConfig, WangLandauCore, WangLandauPhase,
    WangLandauRefinement, WangLandauRunControl, WangLandauState, WangLandauTermination,
};
pub use lattice::graph::{
    build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square, build_triangular,
    Bond, BondType, CsrLattice,
};
pub use lattice::interaction::{
    ClusterAuxiliary, ClusterModel, ContinuousHeatBathable, Hamiltonian, HeatBathable,
    Initializable, LocalFieldModel, Measurable, PairInteraction, Proposable,
};
pub use lattice::models::{HeisenbergModel, IsingModel, ONModel, PottsModel, XYModel};
pub use lattice::proposal::{OPSSStrategy, ProposalStrategy, ProposedSpin, StandardStrategy};
pub use lattice::state::{SiteChange, System};
pub use multi_spin::{MultiSpinIsing, N_REPLICAS};
pub use observables::{
    compute_correlation_1d, DefaultObservableSet, EmptyObservableSet, EnergyPerSite, Magnetization,
    MomentSpec, Observable, ObservableSet, TotalEnergy,
};
pub use particle::{
    compute_total_energy as compute_particle_energy, CanonicalParticleKernel, CellList,
    CutoffTreatment, GrandCanonicalMove, GrandCanonicalPatch, InsertDeleteParticle,
    IsotropicVolumeChange, LennardJones, LennardJonesMuVt, LennardJonesNpt, LennardJonesNvt,
    LennardJonesSpecies, LogVolumeScale, MolecularMetropolisCore, MolecularMoveKind,
    MoleculeTopology, MoveMixture, OrthorhombicCell, PairPotential, ParticleAlgorithm,
    ParticleBatchMove, ParticleBatchPatch, ParticleConfiguration, ParticleDeletion,
    ParticleEnergyPatch, ParticleError, ParticleGrandCanonicalCore, ParticleInsertion, ParticleMC,
    ParticleMetropolisCore, ParticleNptMetropolisCore, ParticleSystem, ParticleTranslation,
    RigidMoleculeRotation, RigidMoleculeTranslation, SimulationCell, TorsionDefinition,
    TorsionRotation, TranslateParticle, VolumeChangePatch, WeightedMove,
};
pub use postprocess::{binder_cumulant, specific_heat, susceptibility};
pub use statistics::{statistical_efficiency, StatisticalEfficiency};
