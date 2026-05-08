//! Lattice topology module.

mod bond;
pub mod builders;

pub use bond::BondType;
pub use builders::{build_chain, build_square};

/// Neighbor entry in adjacency list.
#[derive(Clone, Debug)]
pub struct Neighbor {
    /// Target site index
    pub target: usize,
    /// Bond type for direction-dependent weights
    pub bond_type: BondType,
}

/// Lattice topology represented as adjacency list.
#[derive(Clone, Debug)]
pub struct Lattice {
    /// Adjacency list: sites[i] = neighbors of site i
    pub sites: Vec<Vec<Neighbor>>,
    /// Total number of sites
    pub n_sites: usize,
    /// Total number of bonds (counting each bond once)
    pub n_bonds: usize,
}