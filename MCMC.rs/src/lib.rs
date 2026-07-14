//! Statistical Markov-chain Monte Carlo kernels for Scuttle.
//!
//! `mcmc-rs` owns target-density evaluation, transition kernels, warmup
//! adaptation, constrained-parameter transforms, posterior traces, replica
//! exchange and multi-chain diagnostics. `carlo-rs` remains responsible for
//! generic run lifecycle and scalar online measurements.

pub mod adaptation;
pub mod carlo_adapter;
pub mod checkpoint;
pub mod diagnostics;
pub mod error;
pub mod gradient;
pub mod integrator;
pub mod kernel;
pub mod metric;
pub mod multichain;
pub mod phase;
pub mod proposal;
pub mod sampler;
pub mod state;
pub mod target;
pub mod tempering;
pub mod trace;
pub mod transform;

pub use adaptation::{
    DenseCovarianceAdaptation, DiagonalCovarianceAdaptation, DualAveraging, HmcWarmup,
    MetricAdaptation, RobbinsMonroScale, StepSizeSearch, WarmupWindowConfig,
};
pub use carlo_adapter::McmcSampler;
pub use checkpoint::{ChainCheckpoint, TargetFingerprint, CHECKPOINT_FORMAT};
pub use diagnostics::{diagnose, energy_bfmi, MultiChainDiagnostics, ParameterDiagnostics};
pub use error::McmcError;
pub use gradient::{
    check_gradient, GradientCheckConfig, GradientCheckReport, GradientComponentCheck,
};
pub use integrator::{IntegrationReport, LeapfrogIntegrator, PhasePoint};
pub use kernel::{
    ComponentWiseMetropolis, GibbsKernel, GibbsUpdate, GibbsUpdateResult, Mixture, Nuts,
    RandomWalkMetropolis, Repeat, SliceSampler, StaticHmc, Then, TransitionKernel,
    TransitionReport,
};
pub use metric::{DenseMetric, DiagonalMetric, Metric, MetricKind, UnitMetric};
pub use multichain::{run_multichain, ChainOutput, McmcConfig, McmcOutput};
pub use phase::SamplingPhase;
pub use proposal::GaussianScale;
pub use sampler::ChainRunner;
pub use state::{ChainState, EuclideanCache, EuclideanState};
pub use target::{DifferentiableLogDensity, FnDifferentiableLogDensity, FnLogDensity, LogDensity};
pub use tempering::{
    run_parallel_tempering, ExchangeEdgeDiagnostics, TemperedChainOutput, TemperingConfig,
    TemperingOutput,
};
pub use trace::{MemoryTrace, TraceStore, TraceView};
pub use transform::{
    Bijector, DifferentiableBijector, Identity, Interval, Ordered, Positive, Product, Simplex,
    TransformError, TransformedTarget,
};
