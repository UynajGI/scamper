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
pub mod topology;
pub mod worm;

// Hierarchical modules

pub mod algorithms;
pub mod audit;
pub mod core;
pub mod dynamics;
pub mod generalized;
pub mod lattice;
pub mod observables;
pub mod particle;
pub mod percolation;

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
pub use dynamics::{
    BklEvent, BklIsingKernel, DynamicsError, EventChainOutcome, GillespieEvent, GillespieKernel,
    HardSphereEventChain, HardSphereEventChainMC, KawasakiCore, KawasakiIsingMC, KineticIsingBklMC,
    KineticIsingModel, KineticRateLaw, RejectionFreeModel,
};
pub use generalized::{
    canonical_reweight, enumerate_ising_density_of_states, BinnedAxis, CanonicalReweighting,
    DiscreteAxis, EnergyBiasCore, EnergyMacrostate, ExactIsingDensityOfStates, FixedBias,
    GeneralizedError, HarmonicUmbrellaBias, Histogram, IsingWangLandau, LogBias,
    LogDensityOfStates, Macrostate, MacrostateAxis, MagnetizationMacrostate, MulticanonicalBias,
    ParticleNumberMacrostate, WangLandauConfig, WangLandauCore, WangLandauPhase,
    WangLandauRefinement, WangLandauRunControl, WangLandauState, WangLandauTermination,
};
pub use lattice::graph::{
    build_chain, build_chain_with_embedding, build_honeycomb, build_honeycomb_with_embedding,
    build_hypercubic, build_hypercubic_with_embedding, build_kagome, build_kagome_with_embedding,
    build_square, build_square_with_embedding, build_triangular, build_triangular_with_embedding,
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
    CutoffTreatment, DipolarExternalField, GrandCanonicalMove, GrandCanonicalPatch,
    InsertDeleteParticle, IsotropicVolumeChange, LennardJones, LennardJonesMuVt, LennardJonesNpt,
    LennardJonesNvt, LennardJonesSpecies, LogVolumeScale, MolecularMetropolisCore,
    MolecularMoveKind, MoleculeTopology, MoveMixture, OrthorhombicCell, PairPotential,
    ParticleAlgorithm, ParticleBatchMove, ParticleBatchPatch, ParticleConfiguration,
    ParticleDeletion, ParticleEnergyPatch, ParticleError, ParticleGrandCanonicalCore,
    ParticleInsertion, ParticleMC, ParticleMetropolisCore, ParticleNptMetropolisCore,
    ParticleSystem, ParticleTranslation, RigidMoleculeRotation, RigidMoleculeTranslation,
    SimulationCell, TorsionDefinition, TorsionRotation, TranslateParticle, VolumeChangePatch,
    WeightedMove,
};
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use percolation::WorkspaceCapacityAudit;
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use percolation::WorkspaceWindingCapacityAudit;
pub use percolation::{
    analyze, analyze_bonds, analyze_bonds_with_labels, analyze_strong_components,
    analyze_weak_components, analyze_winding, analyze_with_labels,
    candidate_pairs_from_configuration, candidate_pairs_with_cell_list, chemical_distance_run,
    csr_from_edge_list, first_passage_run, geodesic_aging, geodesic_count_run, geodesic_path,
    invasion_run, kcore_decompose, long_range_run, measure_diameter_ladder, overlap_pairs,
    overlap_pairs_from_configuration, reachable, reweight, reweight_heterogeneous,
    reweight_heterogeneous_with, reweight_with, scan_bond, scan_site, AchlioptasError,
    AchlioptasGrid, AchlioptasProcessMC, AchlioptasRule, AchlioptasStepRecord, AllActiveVertices,
    AnalysisError, AnalysisResult, AttackGrid, AttackOrder, AttackOrderContext, AttackPolicy,
    BinderCumulant, BinderEstimator, BondBernoulli, BondConfigurationError, BondSample,
    BooleanError, BooleanPercolationMC, BootstrapDynamics, BootstrapError, BootstrapInterval,
    BootstrapOutcome, BootstrapPercolationMC, BootstrapSeeds, BootstrapThreshold,
    BootstrapWorkspace, BoundaryQuery, BoundaryQueryError, CandidatePairGraph, CandidateSampling,
    ChemicalDistanceError, ChemicalDistanceEstimate, ChemicalDistanceEstimator,
    ChemicalDistanceOutcome, ChemicalDistanceSampler, ChemicalDistanceWorkspace, ClusterQ,
    CollapsePoint, CollapseQuality, ComponentLabel, ComponentSummary, ComponentWinding,
    ComponentWorkspace, ContinuumEmbedding, ContinuumError, ContinuumWindingError, CoreThresholds,
    CorrelationPoint, CrossingEstimate, DegreeAttackMC, DegreeScratch, DiameterLadder,
    DiameterLadderError, DiameterRegime, DiameterRow, DiameterTracking, DirectedAnalysisError,
    Direction, DpProcessError, DpSeedPolicy, DuplicateEdgePolicy, EdgeActivity, EdgeListSpec,
    EdgeUniverse, EdgeWeightLaw, EdwardsSokalError, EdwardsSokalMC, ExternalBondConfiguration,
    FirstPassageError, FirstPassageOutcome, FirstPassageWorkspace, FixedSeries, FssCurve, FssError,
    FssPoint, GeodesicAging, GeodesicCountOutcome, GeodesicCountWorkspace, GeodesicPathOutcome,
    GeodesicTracking, HetReweightResult, HetReweightWorkspace, IncrementalError,
    IncrementalSummary, IncrementalUnionWorkspace, InvasionError, InvasionMode, InvasionOutcome,
    InvasionPercolationMC, InvasionStop, InvasionStopReason, InvasionTrapping, InvasionWorkspace,
    KCoreError, KCoreMC, KCoreSummary, KCoreWorkspace, LayeredDirectedPercolationMC, LongRangeBox,
    LongRangeError, LongRangeOutcome, LongRangePercolationMC, LongRangeWorkspace, MixedBernoulli,
    Moments, MultiplexError, MultiplexNetwork, MutualComponent, NetworkError, NetworkInputError,
    NetworkPercolationError, NetworkProfile, NewmanZiffWorkspace, ObservablePlan,
    ObservablePlanError, ObservableSpec, OccupationLaw, OverlapGraph, PairGraph, PassageEstimate,
    PassageEstimator, PassageGeometry, PassageObservablesError, PassageSampler, PercolationAdapter,
    PointProcess, Probability, ProbabilityField, ProcessTrajectory, Radii, RandomClusterError,
    RandomClusterMC, RandomFailureMC, ReachabilityWorkspace, ReweightError, ReweightWorkspace,
    ReweightedEstimate, SampleView, SamplingError, ScalarObservable, ScalingWindow, ScanCurve,
    ScanError, ScanEstimator, ScanMode, ScanObservable, ScanObservableSpec, ScanPoint,
    SelfLoopPolicy, SeriesPad, ShapeError, ShapeEstimate, ShapeEstimator, ShapeGeometry,
    ShapeSampler, ShapeTarget, SiteBernoulli, StaticConfiguration, StaticLaw, StaticObservable,
    StaticPercolationError, StaticPercolationMC, StrongComponentResult, StrongComponentSummary,
    StrongComponentWorkspace, TiePolicy, VertexActivity, WeakComponentResult, WeakComponentSummary,
    WeakComponentWorkspace, WeightLaw, WindingError, WindingResult, WindingSummary,
    WindingWorkspace, DIAMETER_VERTEX_CAP, NO_PREDECESSOR, SHAPE_SIDE_CAP,
};
pub use percolation::{
    bootstrap, bracket, collapse_residual, critical_probability, crossing, flip_bracket,
    inverse_normal_cdf, multiplex_from_edges, phase_scan, second_moment_correlation_length,
    ChayesMachtaError, ChayesMachtaMC, CorrelatedError, CorrelatedFieldSpec,
    CorrelatedGaussianField, CorrelatedPercolationMC, CorrelatedSiteLaw, CorrelatedThresholdError,
    FlipBracket, FractalDimensionEstimate, HarrisBracket, Hypergraph, HypergraphError,
    HypergraphSummary, IntervalEstimate, LinearFit, NuEstimate, PcExtrapolation, PcShiftNuEstimate,
    PcShiftRung, PottsPhaseError, PottsPhaseEstimate, PottsPhaseSpec, RungObservables, RungSlope,
    SweenyAcceptance, SweenyError, SweenyMC, SweenyProposal, SweenyRates, ThresholdLadder,
    ThresholdLadderSpec, WeinribHalperinError, WeinribHalperinMap, BOOTSTRAP_RESAMPLES,
    PEAK_FIT_RADIUS,
};
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use percolation::{
    HetReweightCapacityAudit, ReweightCapacityAudit, WorkspaceScanCapacityAudit,
};
#[cfg(feature = "io")]
pub use percolation::{
    NetworkDocument, NetworkEdge, NetworkIoError, EDGE_LIST_HEADER_V1, NETWORK_JSON_FORMAT,
    NETWORK_JSON_VERSION,
};
pub use postprocess::{
    binder_cumulant, connected_order_parameter_fluctuation, specific_heat, susceptibility,
    zero_field_ising_susceptibility,
};
pub use statistics::{statistical_efficiency, StatisticalEfficiency};
pub use topology::{
    Arc, ArcId, ArcIncidence, BorrowedDirectedCsr, BorrowedDirectedOutCsr, BorrowedUndirectedCsr,
    DirectedCsr, DirectedDisplacement, DirectedGraphView, DirectedInGraphView,
    DirectedOutGraphView, EdgeId, EmbeddingError, GraphView, Incidence, LatticeEmbedding,
    PeriodicEmbedding, TopologyError, UndirectedGraphView, VertexId,
};
pub use worm::{
    enumerate_ising_graph_expansion, EndpointPairHistogram, ExactIsingGraphExpansion,
    IsingComponentWorm, IsingGraphConfiguration, IsingGraphPatch, IsingGraphWormEnsemble,
    IsingGraphWormMC, IsingGraphWormModel, IsingWormStep, WormConfig, WormError, WormKernel,
    WormModel, WormSector, WormState, WormStepDelta, WormStepProposal, WormTransition,
    WormTransitionStatistics,
};
