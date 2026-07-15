//! Continuous-imaginary-time cluster solver for the longitudinal spin-boson
//! model.
//!
//! This backend samples a periodic piecewise-constant `sigma_z(tau)`
//! worldline. Existing real kinks are augmented by Poisson-distributed
//! potential cuts with rate `Delta/2`. Equal-spin segments are connected by
//! the exact retarded bond probability
//!
//! `1 - exp[-2 int_I d tau int_J d tau' K_beta(tau-tau')]`,
//!
//! after which every connected component is assigned a new orientation. At
//! zero longitudinal bias the two orientations are equiprobable; finite bias
//! is handled by an exact cluster heat bath.

pub mod cluster_builder;
pub mod mc;
pub mod retarded_bonds;
pub mod segments;

pub use cluster_builder::{ClusterDiagnostics, ClusterUpdateReport, ContinuousTimeClusterEngine};
pub use mc::{
    measure_cluster_observables, register_cluster_evaluables, LongitudinalClusterObservables,
    LongitudinalSpinBosonClusterQmc,
};
pub use retarded_bonds::{LongitudinalSpinBosonModel, RetardedKernel};
pub use segments::{build_segments, LongitudinalWorldline, TimeInterval, WorldlineSegment};
