//! Trait definitions for classical Monte Carlo models.
//!
//! This module defines the core traits that separate concerns:
//! - [`Hamiltonian`]: Energy computation and physical constants
//! - [`ClusterModel`]: Cluster algorithm support (Wolff, Swendsen-Wang)
//! - [`Proposable`]: Spin proposal for Metropolis
//! - [`Measurable`]: Magnetization computation

use crate::lattice::Lattice;
use rand::Rng;

/// Core physics: energy computation and coupling constants.
///
/// Implementations define the Hamiltonian H = -J Σ s_i · s_j and provide
/// methods to compute local energy contributions. Models are stateless —
/// temperature (β) lives in [`System`](crate::System).
pub trait Hamiltonian: Send + Sync {
    /// Number of spin components per site (1 = Ising/Potts, 2 = XY, 3 = Heisenberg).
    fn spin_dim(&self) -> usize;

    /// Coupling constant J.
    fn coupling(&self) -> f64;

    /// Energy contribution from site `i` interacting with its neighbors,
    /// given a **proposed** new spin and inverse temperature β.
    ///
    /// Returns the contribution to total energy from site `i`, summing over
    /// all bonds connected to `i`. The returned value replaces the old
    /// contribution directly (not a delta).
    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &Lattice,
        site: usize,
        beta: f64,
        proposed: &[f64],
    ) -> f64;

    /// Compute initial energy for the full system.
    fn compute_total_energy(&self, spins: &[f64], lattice: &Lattice, beta: f64) -> f64 {
        let mut total = 0.0;
        let sd = self.spin_dim();
        for site in 0..lattice.n_sites {
            let proposed = spins[site * sd..(site + 1) * sd].to_vec();
            total += self.local_energy(spins, lattice, site, beta, &proposed);
        }
        total / 2.0 // double-counted bonds
    }
}

/// Cluster algorithm support for Wolff and Swendsen-Wang.
///
/// Scalar models (Ising, Potts) implement [`flip_in_place`](ClusterModel::flip_in_place)
/// and [`opposite_spin`](ClusterModel::opposite_spin).
///
/// Vector models (XY, Heisenberg) implement [`reflect`](ClusterModel::reflect) and
/// [`embedding_direction`](ClusterModel::embedding_direction).
pub trait ClusterModel: Hamiltonian {
    /// FK bond percolation probability.
    /// Default: `1 - exp(-2βJ)` for Ising-like models.
    fn fk_bond_probability(&self, beta: f64) -> f64 {
        1.0 - (-2.0 * self.coupling() * beta).exp()
    }

    /// Flip spin in-place (scalar models: Ising, Potts).
    ///
    /// For Ising: negates the spin.
    /// For Potts: assigns a random different state.
    fn flip_in_place(&self, _spin: &mut [f64], _rng: &mut impl Rng) {
        panic!("flip_in_place not implemented for vector models");
    }

    /// Opposite spin for Wolff cluster flip (scalar models only).
    fn opposite_spin(&self, _spin: f64, _rng: &mut impl Rng) -> f64 {
        panic!("opposite_spin only for scalar models");
    }

    /// Reflect spin across plane perpendicular to direction (vector models: XY, Heisenberg).
    fn reflect(&self, _spin: &mut [f64], _direction: &[f64]) {
        panic!("reflect only for vector models");
    }

    /// Random direction for embedding (vector models only).
    fn embedding_direction(&self, _rng: &mut impl Rng) -> Vec<f64> {
        panic!("embedding_direction only for vector models");
    }

    /// Random spin value for SW cluster assignment (scalar models only).
    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        panic!("random_cluster_spin only for scalar models");
    }
}

/// Spin proposal for Metropolis algorithm.
pub trait Proposable: Hamiltonian {
    /// Propose a random new spin (all components). Returns a Vec of length `spin_dim()`.
    fn propose(&self, rng: &mut impl Rng) -> Vec<f64>;

    /// Normalize a spin vector to unit length (for XY/Heisenberg). Ising/Potts skip this.
    fn normalize_spin(&self, _spin: &mut [f64]) {}
}

/// Magnetization computation for measurements.
pub trait Measurable: Hamiltonian {
    /// Total magnetization |M|/N from the current spin configuration.
    fn magnetization(&self, spins: &[f64]) -> f64;
}
