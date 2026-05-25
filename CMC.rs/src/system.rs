//! Mutable simulation state.

use crate::lattice::Lattice;

/// Mutable state of a classical Monte Carlo simulation.
///
/// All fields are `pub` — algorithms read and write them directly.
/// Temperature (β) lives here as runtime state, not in the model,
/// which enables parallel tempering via simple `swap(system.beta)`.
#[derive(Debug, Clone)]
pub struct System {
    /// Lattice topology (immutable after construction).
    pub lattice: Lattice,

    /// Spin configuration, flattened: `spins[site * spin_dim + component]`.
    /// Length = `n_sites × spin_dim`.
    pub spins: Vec<f64>,

    /// Running total energy of the current configuration.
    pub energy: f64,

    /// Inverse temperature β = 1/(k_B T).
    pub beta: f64,
}

impl System {
    /// Create a system with all spins set to `init_value`.
    pub fn new(lattice: Lattice, spin_dim: usize, init_value: f64, beta: f64) -> Self {
        let n = lattice.n_sites * spin_dim;
        Self {
            spins: vec![init_value; n],
            energy: 0.0,
            lattice,
            beta,
        }
    }

    /// Number of lattice sites.
    pub fn n_sites(&self) -> usize {
        self.lattice.n_sites
    }

    /// Slice of spin components at a given site.
    pub fn spin_at(&self, site: usize, spin_dim: usize) -> &[f64] {
        let base = site * spin_dim;
        &self.spins[base..base + spin_dim]
    }

    /// Mutable slice of spin components at a given site.
    pub fn spin_at_mut(&mut self, site: usize, spin_dim: usize) -> &mut [f64] {
        let base = site * spin_dim;
        &mut self.spins[base..base + spin_dim]
    }
}
