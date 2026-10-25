//! Site and bond percolation on arbitrary [`crate::CsrLattice`] graphs.
//!
//! Percolation is sampled as independent configurations rather than a Markov
//! chain: every sweep redraws occupancy (sites or bonds open independently
//! with probability `p`) and every measurement analyses the occupied subgraph
//! through union-find. Samples are i.i.d., so thermalization is meaningless —
//! set `thermalization_sweeps = 0` and read every sweep as one independent
//! sample.
//!
//! The scheduler-ready adapter is [`StaticPercolationMC`], composed from an
//! owned [`crate::CsrLattice`], a uniform or heterogeneous (per-vertex and
//! per-edge [`Probability`] fields) [`StaticLaw`], a
//! [`StaticConfiguration`], a reusable [`ComponentWorkspace`], and an explicit
//! [`ObservablePlan`]. [`analyze`] and [`ComponentSummary`] are public for
//! direct, RNG-free analysis of fixed configurations (exact-enumeration
//! validation builds on this); borrowed graphs and custom laws call them
//! directly, optionally through `carlo_rs::Run::from_parts`.
//!
//! Z1 adds a separate core parameter-scan layer: [`scan_site`]/[`scan_bond`]
//! activate one random permutation incrementally, [`reweight`] applies exact
//! binomial weights to the resulting curve, and [`ScanEstimator`] aggregates
//! independent scans. [`reweight_heterogeneous`] extends the same unchanged
//! curves to per-element probability fields (T2). This is not a scheduler
//! adapter and has no vector `Results` output until C0.

mod activity;
mod analyzer;
mod bonds;
mod boolean;
mod carlo;
mod chayes_machta;
mod configuration;
mod continuum;
mod continuum_winding;
mod correlated;
mod correlated_thresholds;
mod directed;
mod edwards_sokal;
mod first_passage;
mod first_passage_observables;
mod first_passage_shape;
mod law;
mod long_range;
mod long_range_diameter;
mod long_range_distances;
mod network;
mod observable;
mod potts_phase;
mod probability;
mod process;
mod query;
mod random_cluster;
mod scan;
mod summary;
mod sweeny;
mod weinrib_halperin;
mod workspace;

pub use activity::{EdgeActivity, VertexActivity};
pub use analyzer::{
    analyze, analyze_winding, analyze_with_labels, AnalysisError, AnalysisResult, ComponentWinding,
    WindingError, WindingResult, WindingSummary,
};
pub use bonds::{
    analyze_bonds, analyze_bonds_with_labels, BondConfigurationError, BondSample,
    ExternalBondConfiguration,
};
pub use boolean::{
    overlap_pairs, overlap_pairs_from_configuration, BooleanError, BooleanPercolationMC,
    OverlapGraph, PointProcess, Radii,
};
pub use carlo::{StaticLaw, StaticPercolationError, StaticPercolationMC};
pub use chayes_machta::{ChayesMachtaError, ChayesMachtaMC};
pub use configuration::StaticConfiguration;
pub use continuum::{
    candidate_pairs_from_configuration, candidate_pairs_with_cell_list, CandidatePairGraph,
    ContinuumError, PairGraph,
};
pub use continuum_winding::{ContinuumEmbedding, ContinuumWindingError};
pub use correlated::{
    inverse_normal_cdf, CorrelatedError, CorrelatedFieldSpec, CorrelatedGaussianField,
    CorrelatedPercolationMC, CorrelatedSiteLaw,
};
pub use correlated_thresholds::{
    CorrelatedThresholdError, FractalDimensionEstimate, IntervalEstimate, LinearFit, NuEstimate,
    PcExtrapolation, PcShiftNuEstimate, PcShiftRung, RungObservables, RungSlope, ThresholdLadder,
    ThresholdLadderSpec, BOOTSTRAP_RESAMPLES, PEAK_FIT_RADIUS,
};
pub use directed::{
    analyze_strong_components, analyze_weak_components, reachable, DirectedAnalysisError,
    DpProcessError, DpSeedPolicy, LayeredDirectedPercolationMC, ReachabilityWorkspace,
    StrongComponentResult, StrongComponentSummary, StrongComponentWorkspace, WeakComponentResult,
    WeakComponentSummary, WeakComponentWorkspace,
};
pub use edwards_sokal::{AllActiveVertices, EdwardsSokalError, EdwardsSokalMC};
pub use first_passage::{
    first_passage_run, FirstPassageError, FirstPassageOutcome, FirstPassageWorkspace,
    GeodesicTracking, NO_PREDECESSOR,
};
pub use first_passage_observables::{
    EdgeWeightLaw, PassageEstimate, PassageEstimator, PassageGeometry, PassageObservablesError,
    PassageSampler,
};
pub use first_passage_shape::{
    geodesic_aging, geodesic_count_run, geodesic_path, Direction, GeodesicAging,
    GeodesicCountOutcome, GeodesicCountWorkspace, GeodesicPathOutcome, ShapeError, ShapeEstimate,
    ShapeEstimator, ShapeGeometry, ShapeSampler, ShapeTarget, SHAPE_SIDE_CAP,
};
pub use law::{BondBernoulli, MixedBernoulli, SamplingError, SiteBernoulli};
pub use long_range::{
    long_range_run, LongRangeBox, LongRangeError, LongRangeOutcome, LongRangePercolationMC,
    LongRangeWorkspace,
};
pub use long_range_diameter::{
    measure_diameter_ladder, DiameterLadder, DiameterLadderError, DiameterRegime, DiameterRow,
};
pub use long_range_distances::{
    chemical_distance_run, ChemicalDistanceError, ChemicalDistanceEstimate,
    ChemicalDistanceEstimator, ChemicalDistanceOutcome, ChemicalDistanceSampler,
    ChemicalDistanceWorkspace, DiameterTracking, DIAMETER_VERTEX_CAP,
};
pub use network::{
    csr_from_edge_list, AttackGrid, AttackPolicy, DegreeAttackMC, DegreeScratch,
    DuplicateEdgePolicy, EdgeListSpec, NetworkError, NetworkInputError, NetworkPercolationError,
    NetworkProfile, RandomFailureMC, SelfLoopPolicy,
};
pub use observable::{
    ObservablePlan, ObservablePlanError, ObservableSpec, SampleView, ScalarObservable,
    StaticObservable,
};
pub use potts_phase::{
    critical_probability, phase_scan, PottsPhaseError, PottsPhaseEstimate, PottsPhaseSpec,
};
pub use probability::{Probability, ProbabilityError, ProbabilityField};
pub use process::{
    invasion_run, kcore_decompose, AchlioptasError, AchlioptasGrid, AchlioptasProcessMC,
    AchlioptasRule, AchlioptasStepRecord, BootstrapDynamics, BootstrapError, BootstrapOutcome,
    BootstrapPercolationMC, BootstrapSeeds, BootstrapThreshold, BootstrapWorkspace,
    CandidateSampling, CoreThresholds, EdgeUniverse, IncrementalError, IncrementalSummary,
    IncrementalUnionWorkspace, InvasionError, InvasionMode, InvasionOutcome, InvasionPercolationMC,
    InvasionStop, InvasionStopReason, InvasionTrapping, InvasionWorkspace, KCoreError, KCoreMC,
    KCoreSummary, KCoreWorkspace, NeighborCounters, ProcessError, ProcessQueue, ProcessTrajectory,
    StampFlags, TiePolicy,
};
pub use query::{BoundaryQuery, BoundaryQueryError};
pub use random_cluster::{ClusterQ, RandomClusterError, RandomClusterMC};
pub use scan::{
    reweight, reweight_heterogeneous, scan_bond, scan_site, HetReweightResult,
    HetReweightWorkspace, NewmanZiffWorkspace, ReweightError, ReweightWorkspace,
    ReweightedEstimate, ScanCurve, ScanError, ScanEstimator, ScanMode, ScanObservable, ScanPoint,
};
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use scan::{HetReweightCapacityAudit, ReweightCapacityAudit, WorkspaceScanCapacityAudit};
pub use summary::ComponentSummary;
pub use sweeny::{SweenyAcceptance, SweenyError, SweenyMC, SweenyProposal, SweenyRates};
pub use weinrib_halperin::{
    bracket, flip_bracket, FlipBracket, HarrisBracket, WeinribHalperinError, WeinribHalperinMap,
};
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use workspace::WorkspaceCapacityAudit;
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use workspace::WorkspaceWindingCapacityAudit;
pub use workspace::{ComponentLabel, ComponentWorkspace, WindingWorkspace};
