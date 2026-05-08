//! Measurements for SSE algorithm.

use super::SSEEngine;
use crate::hilbert::HilbertSpace;

impl<H: HilbertSpace> SSEEngine<H> {
    /// Compute energy from operator sequence.
    ///
    /// For the shifted Hamiltonian where we expand in -H' = -H + C:
    ///   ⟨n⟩ = β(-E + C_total)
    ///   E = C_total / N - ⟨n⟩ / (β * N)
    /// where C_total is the total diagonal shift.
    pub fn compute_energy(&self) -> f64 {
        if self.beta <= 0.0 || self.lattice.n_sites == 0 {
            return 0.0;
        }

        let n = self.op_seq.n_operators as f64;
        let beta = self.beta;
        let n_sites = self.lattice.n_sites as f64;

        -n / (beta * n_sites) + self.diagonal_shift / n_sites
    }

    /// Compute magnetization from spin configuration.
    pub fn compute_magnetization(&self) -> f64 {
        let n_sites = self.lattice.n_sites;

        let m: i32 = self
            .spins
            .iter()
            .map(|s| if *s == 0 { 1 } else { -1 })
            .sum();

        m.abs() as f64 / n_sites as f64
    }
}
