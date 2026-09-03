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

mod cluster;
mod mc;
mod state;

pub use cluster::{cluster_stats, ClusterStats, UnionFind};
pub use mc::PercolationMC;
pub use state::{OccupancyState, PercolationMode};
