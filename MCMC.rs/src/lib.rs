//! Statistical Markov-chain Monte Carlo kernels for Scuttle.
//!
//! `mcmc-rs` owns target-density evaluation, transition kernels, warmup
//! adaptation, posterior traces and multi-chain diagnostics. `carlo-rs` remains
//! responsible for generic run lifecycle and scalar online measurements.

pub mod adaptation;
pub mod carlo_adapter;
pub mod checkpoint;
pub mod diagnostics;
pub mod error;
pub mod kernel;
pub mod multichain;
pub mod phase;
pub mod proposal;
pub mod sampler;
pub mod state;
pub mod target;
pub mod trace;

pub use adaptation::{DiagonalCovarianceAdaptation, RobbinsMonroScale};
pub use carlo_adapter::McmcSampler;
pub use checkpoint::{ChainCheckpoint, TargetFingerprint, CHECKPOINT_FORMAT};
pub use diagnostics::{diagnose, MultiChainDiagnostics, ParameterDiagnostics};
pub use error::McmcError;
pub use kernel::{
    ComponentWiseMetropolis, RandomWalkMetropolis, SliceSampler, TransitionKernel, TransitionReport,
};
pub use multichain::{run_multichain, ChainOutput, McmcConfig, McmcOutput};
pub use phase::SamplingPhase;
pub use proposal::GaussianScale;
pub use sampler::ChainRunner;
pub use state::{ChainState, EuclideanCache, EuclideanState};
pub use target::{DifferentiableLogDensity, FnLogDensity, LogDensity};
pub use trace::{MemoryTrace, TraceStore, TraceView};
