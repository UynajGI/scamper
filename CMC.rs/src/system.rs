//! Mutable configuration and canonical-ensemble state.

use crate::hamiltonian::Hamiltonian;
use crate::lattice::CsrLattice;

/// Explicit one-site change, useful for custom update kernels.
#[derive(Debug, Clone)]
pub struct SiteChange {
    pub site: usize,
    pub new_spin: Vec<f64>,
}

/// Mutable state of one classical Monte Carlo chain.
#[derive(Debug, Clone)]
pub struct System {
    pub lattice: CsrLattice,
    /// Flat site-major spin storage.
    pub spins: Vec<f64>,
    /// Cached physical energy, never multiplied by beta.
    pub energy: f64,
    /// Canonical inverse temperature.
    pub beta: f64,
}

impl System {
    pub fn new(lattice: CsrLattice, spin_dim: usize, init_value: f64, beta: f64) -> Self {
        assert!(spin_dim > 0, "spin dimension must be positive");
        assert!(
            beta.is_finite() && beta >= 0.0,
            "beta must be finite and non-negative"
        );
        let spins = vec![init_value; lattice.n_sites * spin_dim];
        Self {
            lattice,
            spins,
            energy: 0.0,
            beta,
        }
    }

    #[inline]
    pub fn n_sites(&self) -> usize {
        self.lattice.n_sites
    }

    #[inline]
    pub fn spin_at(&self, site: usize, spin_dim: usize) -> &[f64] {
        let base = site * spin_dim;
        &self.spins[base..base + spin_dim]
    }

    #[inline]
    pub fn spin_at_mut(&mut self, site: usize, spin_dim: usize) -> &mut [f64] {
        let base = site * spin_dim;
        &mut self.spins[base..base + spin_dim]
    }

    pub fn set_beta(&mut self, beta: f64) -> Result<(), String> {
        if !beta.is_finite() || beta < 0.0 {
            return Err("beta must be finite and non-negative".to_string());
        }
        self.beta = beta;
        Ok(())
    }

    pub fn recompute_energy<H: Hamiltonian>(&mut self, model: &H) -> f64 {
        self.energy = model.compute_total_energy(&self.spins, &self.lattice, self.beta);
        self.energy
    }

    pub fn energy_error<H: Hamiltonian>(&self, model: &H) -> f64 {
        self.energy - model.compute_total_energy(&self.spins, &self.lattice, self.beta)
    }

    pub fn validate<H: Hamiltonian>(&self, model: &H) -> Result<(), String> {
        self.lattice.validate()?;
        let expected = self.n_sites() * model.spin_dim();
        if self.spins.len() != expected {
            return Err(format!(
                "spin buffer length mismatch: expected {expected}, got {}",
                self.spins.len()
            ));
        }
        if !self.beta.is_finite() || self.beta < 0.0 {
            return Err("beta must be finite and non-negative".to_string());
        }
        if self.spins.iter().any(|value| !value.is_finite()) {
            return Err("spin buffer contains non-finite values".to_string());
        }
        if !self.energy.is_finite() {
            return Err("cached energy is non-finite".to_string());
        }
        Ok(())
    }
}
