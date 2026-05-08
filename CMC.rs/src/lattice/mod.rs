//! Lattice topology module.

mod bond;
mod builders;

pub use bond::BondType;
pub use builders::{build_chain, build_square};

/// Neighbor entry in lattice adjacency list.
#[derive(Clone, Debug)]
pub struct Neighbor {
    pub target: usize,
    pub bond_type: BondType,
}

/// Lattice represented as adjacency list.
#[derive(Clone, Debug)]
pub struct Lattice {
    pub sites: Vec<Vec<Neighbor>>,
    pub n_sites: usize,
    pub n_bonds: usize,
}

/// Domain trait for lattice Monte Carlo methods.
pub trait LatticeMC {
    fn lattice(&self) -> &Lattice;

    fn n_sites(&self) -> usize {
        self.lattice().n_sites
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_empty() {
        let lattice = Lattice {
            sites: vec![],
            n_sites: 0,
            n_bonds: 0,
        };
        assert_eq!(lattice.n_sites, 0);
    }

    #[test]
    fn test_lattice_single_site() {
        let lattice = Lattice {
            sites: vec![vec![]],
            n_sites: 1,
            n_bonds: 0,
        };
        assert_eq!(lattice.n_sites, 1);
        assert!(lattice.sites[0].is_empty());
    }
}
