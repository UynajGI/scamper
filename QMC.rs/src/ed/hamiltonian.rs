//! Sparse Hamiltonian matrix in CSR format for spin-1/2 Heisenberg model.
//!
//! H = J Σ_{<i,j>} S_i · S_j
//!
//! The Hilbert space dimension is 2^N. Each spin configuration |s_0 s_1 ... s_{N-1}>
//! maps to integer index Σ s_i * 2^i where s_i ∈ {0,1} (0=↑, 1=↓).

use crate::hilbert::LocalState;
use crate::lattice::Lattice;

/// CSR-format sparse Hamiltonian for spin-1/2 Heisenberg.
pub struct SparseHamiltonian {
    /// Hilbert space dimension (2^N)
    dim: usize,
    /// CSR row pointers (length dim + 1)
    row_ptr: Vec<usize>,
    /// CSR column indices
    col_idx: Vec<usize>,
    /// CSR values
    values: Vec<f64>,
}

impl SparseHamiltonian {
    /// Build Heisenberg Hamiltonian from lattice.
    ///
    /// H = J Σ S_i · S_j = J Σ [ Sz_i Sz_j + ½(S⁺_i S⁻_j + S⁻_i S⁺_j) ]
    ///
    /// Only practical for N ≤ 16 (dim ≤ 65536).
    pub fn from_heisenberg(lattice: &Lattice, j: f64) -> Self {
        let n_sites = lattice.n_sites;
        let dim = 1 << n_sites;

        // Build unique bond list
        let mut bonds = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for (site_i, neighbors) in lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let key = if site_i < neighbor.target {
                    (site_i, neighbor.target)
                } else {
                    (neighbor.target, site_i)
                };
                if seen.insert(key) {
                    bonds.push((site_i, neighbor.target));
                }
            }
        }

        // First pass: count non-zeros per row
        let mut row_ptr = vec![0usize; dim + 1];

        for state in 0..dim {
            let mut nnz = 0;
            for &(si, sj) in &bonds {
                let spin_i = ((state >> si) & 1) as LocalState;
                let spin_j = ((state >> sj) & 1) as LocalState;
                // Diagonal Sz_i Sz_j term — always present
                nnz += 1;
                // Off-diagonal S⁺S⁻ + S⁻S⁺ — only for anti-aligned spins
                if spin_i != spin_j {
                    nnz += 1;
                }
            }
            row_ptr[state + 1] = row_ptr[state] + nnz;
        }

        let total_nnz = row_ptr[dim];
        let mut col_idx = Vec::with_capacity(total_nnz);
        let mut values = Vec::with_capacity(total_nnz);

        // Second pass: fill values
        // Pre-allocate next_col to track insertion position per row
        let mut next_col = row_ptr[..dim].to_vec();

        for state in 0..dim {
            for &(si, sj) in &bonds {
                let spin_i = ((state >> si) & 1) as LocalState;
                let spin_j = ((state >> sj) & 1) as LocalState;

                // Diagonal: J * Sz_i * Sz_j
                // Sz = +½ for spin=0 (↑), -½ for spin=1 (↓)
                let sz_i = if spin_i == 0 { 0.5 } else { -0.5 };
                let sz_j = if spin_j == 0 { 0.5 } else { -0.5 };
                let pos = next_col[state];
                col_idx.push(state);
                values.push(j * sz_i * sz_j);
                next_col[state] = pos + 1;

                // Off-diagonal: J/2 * (S⁺_i S⁻_j + S⁻_i S⁺_j)
                // S⁺|↓> = |↑>, S⁻|↑> = |↓>, S⁺|↑> = 0, S⁻|↓> = 0
                // For anti-aligned spins: S⁺_i S⁻_j |↑↓> = |↓↑>, etc.
                // Matrix element = J/2
                if spin_i != spin_j {
                    let flipped = state ^ (1 << si) ^ (1 << sj);
                    let pos = next_col[state];
                    col_idx.push(flipped);
                    values.push(j * 0.5);
                    next_col[state] = pos + 1;
                }
            }
        }

        SparseHamiltonian {
            dim,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Matrix-vector product: y = H · x
    pub fn mat_vec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.dim);
        debug_assert_eq!(y.len(), self.dim);
        y.fill(0.0);

        for i in 0..self.dim {
            let xi = x[i];
            if xi.abs() < 1e-15 {
                continue;
            }
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                y[self.col_idx[j]] += self.values[j] * xi;
            }
        }
    }

    /// Hilbert space dimension.
    pub fn dim(&self) -> usize {
        self.dim
    }
}
