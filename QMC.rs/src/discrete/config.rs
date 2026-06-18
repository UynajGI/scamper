//! Discrete-time Suzuki-Trotter space-time configuration.
//!
//! `SpaceTimeConfig` is the full many-body spin configuration on an
//! `N_sites × N_slices` periodic space-time torus: spatial periodicity comes
//! from the lattice (chain PBC), temporal periodicity from the trace
//! (slice `M` wraps to slice `0`). This is what the worm walks through.

use crate::hamiltonian::QuantumHamiltonian;
use crate::lattice::ChainLattice;

/// Spin index: `0 = ↓`, `1 = ↑` (S = 1/2).
pub type Spin = u8;

/// Discrete-time space-time configuration.
///
/// Storage is **slice-major**: `spins[slice * n_sites + site]`. This puts
/// each Trotter slice in a contiguous block, which is the access pattern the
/// worm uses when it moves along the time axis.
#[derive(Debug, Clone)]
pub struct SpaceTimeConfig {
    /// Flattened spins, length `n_sites * n_slices`. Each entry is `0` or `1`.
    pub spins: Vec<Spin>,
    /// Number of lattice sites.
    pub n_sites: usize,
    /// Number of Trotter slices `M`.
    pub n_slices: usize,
    /// Inverse temperature β.
    pub beta: f64,
    /// Slice width `Δτ = β / M`.
    pub dtau: f64,
    /// Lattice topology.
    pub lattice: ChainLattice,
}

impl SpaceTimeConfig {
    /// Create a configuration with all slices initialized to a random spin.
    pub fn new_random(
        lattice: ChainLattice,
        beta: f64,
        n_slices: usize,
        rng: &mut impl rand::Rng,
    ) -> Self {
        use rand::RngExt;
        let n_sites = lattice.n_sites;
        let total = n_sites * n_slices;
        let spins: Vec<Spin> = (0..total).map(|_| rng.random::<bool>() as u8).collect();
        Self {
            spins,
            n_sites,
            n_slices,
            beta,
            dtau: beta / n_slices as f64,
            lattice,
        }
    }

    /// Create a configuration with every spin set to `s` (useful for tests).
    pub fn new_uniform(lattice: ChainLattice, beta: f64, n_slices: usize, s: Spin) -> Self {
        let n_sites = lattice.n_sites;
        let total = n_sites * n_slices;
        Self {
            spins: vec![s; total],
            n_sites,
            n_slices,
            beta,
            dtau: beta / n_slices as f64,
            lattice,
        }
    }

    /// Spin at `(site, slice)`.
    #[inline]
    pub fn spin(&self, site: usize, slice: usize) -> Spin {
        self.spins[slice * self.n_sites + site]
    }

    /// Mutable spin at `(site, slice)`.
    #[inline]
    pub fn spin_mut(&mut self, site: usize, slice: usize) -> &mut Spin {
        &mut self.spins[slice * self.n_sites + site]
    }

    /// Flip spin at `(site, slice)`.
    #[inline]
    pub fn flip(&mut self, site: usize, slice: usize) {
        let idx = slice * self.n_sites + site;
        self.spins[idx] ^= 1;
    }

    /// Total energy of the configuration (diagonal estimator — classical/Néel
    /// contribution only).
    ///
    /// For each space-time site `(i, τ)`, sum the per-bond physical energy
    /// of bonds connected to `i`, averaged over the slice's extent, then
    /// divide by `M`. This is the diagonal-only estimator; for quantum
    /// models with off-diagonal sectors (Heisenberg), use
    /// [`energy_quantum`](Self::energy_quantum) instead.
    pub fn energy<H: QuantumHamiltonian>(&self, ham: &H) -> f64 {
        let m = self.n_slices as f64;
        let mut e = 0.0;
        for slice in 0..self.n_slices {
            for (i, j) in self.lattice.bonds() {
                let s_i = self.spin(i, slice);
                let s_j = self.spin(j, slice);
                e += ham.energy_per_bond(s_i, s_j);
            }
        }
        e / m
    }

    /// Full quantum path-integral energy estimator.
    ///
    /// Sums `bond_energy_estimator` over every space-time bond, classifying
    /// each as kink-bearing (off-diagonal sector) or not:
    /// - **Spatial** bond `(i,j)` at slice `τ`: kink iff `s_i(τ) ≠ s_j(τ)`
    ///   (the spin-exchange operator is the active matrix element).
    /// - **Temporal** bond `(i, τ)↔(i, τ±1)`: kink iff `s_i(τ) ≠ s_i(τ±1)`.
    ///
    /// Result divided by `M`. For Heisenberg this reproduces the quantum
    /// ⟨H⟩ including spin-exchange fluctuations that lower the AF chain
    /// ground state below the classical −J/4 Néel value.
    pub fn energy_quantum<H: QuantumHamiltonian>(&self, ham: &H) -> f64 {
        let m = self.n_slices;
        let dtau = self.dtau;
        let mut e = 0.0;

        for slice in 0..m {
            for (i, j) in self.lattice.bonds() {
                let s_i = self.spin(i, slice);
                let s_j = self.spin(j, slice);
                e += ham.bond_energy_estimator(s_i, s_j, dtau, s_i != s_j);
            }
        }
        for slice in 0..m {
            let next = (slice + 1) % m;
            for site in 0..self.n_sites {
                let s = self.spin(site, slice);
                let s_next = self.spin(site, next);
                e += ham.bond_energy_estimator(s, s_next, dtau, s != s_next);
            }
        }
        e / m as f64
    }

    /// Magnetization per site, averaged over space-time:
    /// `m = (1/NM) Σ_{i,τ} (2 s_{i,τ} - 1)` ∈ [-1, 1].
    pub fn magnetization(&self) -> f64 {
        let total = self.spins.len() as f64;
        let sum: i64 = self.spins.iter().map(|&s| 2 * s as i64 - 1).sum();
        sum as f64 / total
    }

    /// Number of kinks = number of (site, slice-bond) pairs where the spin
    /// differs between adjacent slices. Equal to the number of off-diagonal
    /// operators (spin exchanges) in the configuration.
    pub fn num_kinks(&self) -> usize {
        let mut count = 0;
        for slice in 0..self.n_slices {
            let next = (slice + 1) % self.n_slices;
            for site in 0..self.n_sites {
                if self.spin(site, slice) != self.spin(site, next) {
                    count += 1;
                }
            }
        }
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::HeisenbergChain;
    use crate::lattice::ChainLattice;
    use rand::SeedableRng;

    fn make_rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn uniform_config_has_zero_kinks() {
        let lat = ChainLattice::new(4);
        let cfg = SpaceTimeConfig::new_uniform(lat, 4.0, 8, 1);
        assert_eq!(cfg.num_kinks(), 0);
    }

    #[test]
    fn kink_count_between_alternating_slices() {
        // Construct by hand: slice 0 all ↑, slice 1 all ↓, M=2.
        let lat = ChainLattice::new(4);
        let mut cfg = SpaceTimeConfig::new_uniform(lat, 2.0, 2, 1);
        for site in 0..cfg.n_sites {
            *cfg.spin_mut(site, 1) = 0; // slice 1 = ↓
        }
        // Each of the 4 sites flips between slice 0↔1 and slice 1↔0 (periodic).
        // Kinks counted per (site, slice→next): 4 sites × 2 transitions = 8.
        assert_eq!(cfg.num_kinks(), 8);
    }

    #[test]
    fn energy_uniform_ferromagnetic_ground_state() {
        // All ↑, ferromagnetic Heisenberg: every bond E = -J/4.
        // N=4 chain has 4 bonds; E = 4 × (-J/4) = -J = -1.0
        let lat = ChainLattice::new(4);
        let cfg = SpaceTimeConfig::new_uniform(lat, 4.0, 8, 1);
        let ham = HeisenbergChain::new(1.0);
        let e = cfg.energy(&ham);
        assert!((e - (-1.0)).abs() < 1e-12, "ferro ground state E = {e}");
    }

    #[test]
    fn magnetization_uniform_up() {
        let lat = ChainLattice::new(4);
        let cfg = SpaceTimeConfig::new_uniform(lat, 4.0, 8, 1);
        assert!((cfg.magnetization() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn magnetization_random_near_zero() {
        let lat = ChainLattice::new(64);
        let mut rng = make_rng();
        let cfg = SpaceTimeConfig::new_random(lat, 1.0, 16, &mut rng);
        // With 64×16 = 1024 random spins, |m| should be small.
        assert!(cfg.magnetization().abs() < 0.2);
    }
}
