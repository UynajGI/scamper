//! SSE engine and operator sequence.

use crate::hilbert::{HilbertSpace, LocalState, OpType};
use crate::lattice::{BondType, Lattice};
use std::collections::HashMap;

use super::vertex_data::VertexData;

/// Vertex in SSE operator sequence.
#[derive(Clone, Debug)]
pub struct Vertex {
    /// Bond index in lattice
    pub bond_idx: usize,
    /// Operator type at this position
    pub op: OpType,
    /// Vertex sub-index encoding specific spin configuration.
    /// 0=Identity, 1-4=Diagonal(↑↑,↑↓,↓↑,↓↓), 5-6=OffDiagonal(↑↓→↓↑,↓↑→↑↓)
    pub vertex_idx: u8,
}

impl Vertex {
    /// Get operator type from vertex_idx.
    #[inline]
    pub fn op_type(&self) -> OpType {
        VertexData::op_type(self.vertex_idx)
    }
}

impl Default for Vertex {
    fn default() -> Self {
        Vertex {
            bond_idx: 0,
            op: OpType::Identity,
            vertex_idx: 0,
        }
    }
}

/// SSE operator sequence (operator string in imaginary time).
#[derive(Clone, Debug)]
pub struct OperatorSequence {
    /// Fixed-length array of vertices, filled with Identity
    pub vertices: Vec<Vertex>,
    /// Count of non-Identity operators
    pub n_operators: usize,
    /// Maximum sequence length (M = N_sites × β × factor)
    pub max_length: usize,
}

impl OperatorSequence {
    /// Create new empty operator sequence.
    pub fn new(max_length: usize) -> Self {
        let vertices = vec![Vertex::default(); max_length];
        OperatorSequence {
            vertices,
            n_operators: 0,
            max_length,
        }
    }
}

/// SSE engine with generic HilbertSpace for zero-cost abstraction.
pub struct SSEEngine<H: HilbertSpace> {
    /// Lattice topology
    pub lattice: Lattice,
    /// Current spin/particle configuration
    pub spins: Vec<LocalState>,
    /// Operator sequence
    pub op_seq: OperatorSequence,
    /// HilbertSpace implementation (static dispatch)
    pub hs: H,
    /// Coupling constants from bond_operators()
    pub weights: HashMap<BondType, f64>,
    /// Flattened bond list: [(site_i, site_j, bond_type), ...]
    pub bond_list: Vec<(usize, usize, BondType)>,
    /// Inverse temperature
    pub beta: f64,
    /// Total diagonal shift constant (makes matrix elements positive).
    pub diagonal_shift: f64,
    /// Loop parameter (e.g., XXZ anisotropy Δ).
    /// Used by directed loop algorithm for scatter probabilities.
    /// Δ = 1.0 (default) corresponds to Heisenberg.
    pub loop_param: f64,
}

impl<H: HilbertSpace> SSEEngine<H> {
    /// Build flattened bond list from lattice adjacency structure.
    ///
    /// Returns a list of (site_i, site_j, bond_type) tuples for each unique bond.
    fn build_bond_list(lattice: &Lattice) -> Vec<(usize, usize, BondType)> {
        let mut bonds = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for (site_i, neighbors) in lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                // Create a canonical bond key (smaller index first)
                let bond_key = if site_i < neighbor.target {
                    (site_i, neighbor.target)
                } else {
                    (neighbor.target, site_i)
                };
                if seen.insert(bond_key) {
                    bonds.push((site_i, neighbor.target, neighbor.bond_type));
                }
            }
        }
        bonds
    }

    /// Create new SSE engine.
    pub fn new(
        lattice: Lattice,
        hs: H,
        max_length: usize,
        weights: HashMap<BondType, f64>,
        beta: f64,
        diagonal_shift: f64,
        loop_param: f64,
    ) -> Self {
        let n_sites = lattice.n_sites;
        // Initialize spins in alternating (AFM) pattern to allow diagonal operator insertion.
        // With the shifted Heisenberg, diagonal operators only exist for anti-aligned spins.
        let spins: Vec<LocalState> = (0..n_sites).map(|i| (i % 2) as LocalState).collect();
        let op_seq = OperatorSequence::new(max_length);
        let bond_list = Self::build_bond_list(&lattice);

        // Validate that all bond types in bond_list have entries in weights
        #[cfg(debug_assertions)]
        {
            for (_, _, bond_type) in &bond_list {
                assert!(
                    weights.contains_key(bond_type),
                    "SSEEngine::new(): Missing weight for bond type {:?} in weights map. \
                     All bond types in the lattice must have corresponding weights.",
                    bond_type
                );
            }
        }

        SSEEngine {
            lattice,
            spins,
            op_seq,
            hs,
            weights,
            bond_list,
            beta,
            diagonal_shift,
            loop_param,
        }
    }
}