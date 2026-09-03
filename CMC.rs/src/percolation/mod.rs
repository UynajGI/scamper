//! Site and bond percolation on arbitrary [`crate::CsrLattice`] graphs.
//!
//! Percolation is sampled as independent configurations rather than a Markov
//! chain: every sweep redraws occupancy (sites or bonds open independently
//! with probability `p`) and every measurement analyses the occupied subgraph
//! through union-find. Samples are i.i.d., so thermalization is meaningless —
//! set `thermalization_sweeps = 0` and read every sweep as one independent
//! sample.
//!
//! The scheduler-ready adapter is [`PercolationMC`]. [`cluster_stats`] and
//! [`UnionFind`] are public for direct, RNG-free analysis of fixed
//! configurations (exact-enumeration validation builds on this).

mod activity;
mod cluster;
mod configuration;
mod law;
mod mc;
mod probability;
mod state;

pub use activity::{EdgeActivity, VertexActivity};
pub use cluster::{cluster_stats, ClusterStats, UnionFind};
pub use configuration::StaticConfiguration;
pub use law::{BondBernoulli, MixedBernoulli, SamplingError, SiteBernoulli};
pub use mc::PercolationMC;
pub use probability::{Probability, ProbabilityError, ProbabilityField};
pub use state::{OccupancyState, PercolationMode};
