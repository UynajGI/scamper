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
//! owned [`crate::CsrLattice`], a uniform [`StaticLaw`], a
//! [`StaticConfiguration`], a reusable [`ComponentWorkspace`], and an explicit
//! [`ObservablePlan`]. [`analyze`] and [`ComponentSummary`] are public for
//! direct, RNG-free analysis of fixed configurations (exact-enumeration
//! validation builds on this); borrowed graphs and custom laws call them
//! directly, optionally through `carlo_rs::Run::from_parts`.

mod activity;
mod analyzer;
mod carlo;
mod configuration;
mod law;
mod observable;
mod probability;
mod query;
mod summary;
mod workspace;

pub use activity::{EdgeActivity, VertexActivity};
pub use analyzer::{analyze, analyze_with_labels, AnalysisError, AnalysisResult};
pub use carlo::{StaticLaw, StaticPercolationError, StaticPercolationMC};
pub use configuration::StaticConfiguration;
pub use law::{BondBernoulli, MixedBernoulli, SamplingError, SiteBernoulli};
pub use observable::{ObservablePlan, ObservablePlanError, StaticObservable};
pub use probability::{Probability, ProbabilityError, ProbabilityField};
pub use query::{BoundaryQuery, BoundaryQueryError};
pub use summary::ComponentSummary;
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use workspace::WorkspaceCapacityAudit;
pub use workspace::{ComponentLabel, ComponentWorkspace};
