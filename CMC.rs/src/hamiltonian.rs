//! Trait definitions for classical Monte Carlo models.
//!
//! This module defines the core traits that separate concerns:
//! - [`Hamiltonian`]: Energy computation and physical constants
//! - [`ClusterModel`]: Cluster algorithm support (Wolff, Swendsen-Wang)
//! - [`Proposable`]: Spin proposal for Metropolis
//! - [`Measurable`]: Magnetization computation

use crate::lattice::CsrLattice;
use rand::Rng;
use smallvec::SmallVec;

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
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        proposed: &[f64],
    ) -> f64;

    /// Compute initial energy for the full system.
    fn compute_total_energy(&self, spins: &[f64], lattice: &CsrLattice, beta: f64) -> f64 {
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
    fn embedding_direction(&self, _rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        panic!("embedding_direction only for vector models");
    }

    /// Random spin value for SW cluster assignment (scalar models only).
    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        panic!("random_cluster_spin only for scalar models");
    }
}

/// Spin proposal for Metropolis algorithm.
pub trait Proposable: Hamiltonian {
    /// Propose a random new spin (all components). Returns a SmallVec of length `spin_dim()`.
    fn propose(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]>;

    /// Normalize a spin vector to unit length (for XY/Heisenberg). Ising/Potts skip this.
    fn normalize_spin(&self, _spin: &mut [f64]) {}
}

/// Magnetization computation for measurements.
pub trait Measurable: Hamiltonian {
    /// Total magnetization |M|/N from the current spin configuration.
    fn magnetization(&self, spins: &[f64]) -> f64;
}

/// Heat-bath (Glauber dynamics) support.
///
/// Models that implement this can be used with [`HeatBathCore`](crate::HeatBathCore),
/// which directly samples each spin from its equilibrium distribution given
/// the local neighbor field — no Metropolis rejection step needed.
pub trait HeatBathable: Hamiltonian {
    /// Number of discrete spin states (2 for Ising, q for Potts).
    fn n_states(&self) -> usize;

    /// Boltzmann weight of each state given the local neighbor field.
    ///
    /// `neighbors` contains the spin values of adjacent sites.
    /// Returns a `Vec<f64>` of length `n_states()`, where `w[k]` is the
    /// un-normalized probability of placing the site in state `k`.
    fn boltzmann_weights(&self, neighbors: &[f64], beta: f64) -> Vec<f64>;

    /// Sample a spin value from the given probability weights.
    ///
    /// `weights` has length `n_states()`. Returns the spin value (e.g. ±1 for Ising).
    fn sample_spin(&self, weights: &[f64], rng: &mut impl Rng) -> f64;
}

/// Continuous-spin heat-bath (Glauber dynamics) support.
///
/// For XY (S¹) and Heisenberg (S²) models, the conditional distribution
/// P(s_i | neighbors) ∝ exp(βJ s_i · h_i) is von Mises / von Mises-Fisher.
/// Sampling is dimension-agnostic — each model implements the appropriate
/// algorithm (Best-Fisher rejection for XY, inverse-CDF for Heisenberg).
pub trait ContinuousHeatBathable: Hamiltonian {
    /// Sample a new spin from P(s_i | neighbors) ∝ exp(βJ s_i · h_i).
    ///
    /// `neighbors` is a flat slice of neighbor spin components (sd per neighbor).
    /// Returns a unit vector of length `spin_dim()`.
    fn heat_bath_sample(
        &self,
        neighbors: &[f64],
        beta: f64,
        rng: &mut impl Rng,
    ) -> SmallVec<[f64; 3]>;
}
