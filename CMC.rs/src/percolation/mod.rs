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
mod configuration;
mod continuum;
mod directed;
mod edwards_sokal;
mod law;
mod network;
mod observable;
mod probability;
mod process;
mod query;
mod random_cluster;
mod scan;
mod summary;
mod workspace;

pub use activity::{EdgeActivity, VertexActivity};
pub use analyzer::{
    analyze, analyze_winding, analyze_with_labels, AnalysisError, AnalysisResult, ComponentWinding,
    WindingError, WindingResult, WindingSummary,
};
pub use bonds::{analyze_bonds, BondConfigurationError, BondSample, ExternalBondConfiguration};
pub use boolean::{
    overlap_pairs, overlap_pairs_from_configuration, BooleanError, BooleanPercolationMC,
    OverlapGraph, PointProcess, Radii,
};
pub use carlo::{StaticLaw, StaticPercolationError, StaticPercolationMC};
pub use configuration::StaticConfiguration;
pub use continuum::{
    candidate_pairs_from_configuration, candidate_pairs_with_cell_list, CandidatePairGraph,
    ContinuumError, PairGraph,
};
pub use directed::{
    analyze_strong_components, analyze_weak_components, reachable, DirectedAnalysisError,
    DpProcessError, DpSeedPolicy, LayeredDirectedPercolationMC, ReachabilityWorkspace,
    StrongComponentResult, StrongComponentSummary, StrongComponentWorkspace, WeakComponentResult,
    WeakComponentSummary, WeakComponentWorkspace,
};
pub use edwards_sokal::{AllActiveVertices, EdwardsSokalError, EdwardsSokalMC};
pub use law::{BondBernoulli, MixedBernoulli, SamplingError, SiteBernoulli};
pub use network::{
    csr_from_edge_list, AttackGrid, AttackPolicy, DegreeAttackMC, DegreeScratch,
    DuplicateEdgePolicy, EdgeListSpec, NetworkError, NetworkInputError, NetworkPercolationError,
    NetworkProfile, RandomFailureMC, SelfLoopPolicy,
};
pub use observable::{ObservablePlan, ObservablePlanError, StaticObservable};
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
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use workspace::WorkspaceCapacityAudit;
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use workspace::WorkspaceWindingCapacityAudit;
pub use workspace::{ComponentLabel, ComponentWorkspace, WindingWorkspace};
