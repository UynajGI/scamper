//! Pluggable observable system for measurements.
//!
//! Users can define custom observables by implementing the [`Observable`] trait
//! and registering them with [`DefaultObservableSet`].

use crate::hamiltonian::{Hamiltonian, Measurable};
use crate::system::System;

/// A single observable that produces a scalar measurement from the current state.
pub trait Observable<H: Hamiltonian>: Send {
    /// Name of the observable (used as key in results).
    fn name(&self) -> &str;

    /// Compute the observable value from the current system state.
    fn measure(&self, system: &System, model: &H) -> f64;
}

// ── Built-in observables ────────────────────────────────────

/// Total energy E.
pub struct TotalEnergy;

impl<H: Hamiltonian> Observable<H> for TotalEnergy {
    fn name(&self) -> &str {
        "Energy"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy
    }
}

/// Energy per site: E/N.
pub struct EnergyPerSite;

impl<H: Hamiltonian> Observable<H> for EnergyPerSite {
    fn name(&self) -> &str {
        "EnergyPerSite"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy / system.n_sites() as f64
    }
}

/// Magnetization per site: |M|/N.
pub struct Magnetization;

impl<H: Hamiltonian + Measurable> Observable<H> for Magnetization {
    fn name(&self) -> &str {
        "Magnetization"
    }

    fn measure(&self, system: &System, model: &H) -> f64 {
        model.magnetization(&system.spins)
    }
}

// ── Observable set ───────────────────────────────────────────

/// A collection of observables that are measured at each measurement sweep.
pub struct DefaultObservableSet<H: Hamiltonian> {
    observables: Vec<Box<dyn Observable<H>>>,
}

impl<H: Hamiltonian + Measurable> DefaultObservableSet<H> {
    /// Create a new set with the default observables (Energy and Magnetization).
    pub fn new() -> Self {
        Self {
            observables: vec![Box::new(TotalEnergy), Box::new(Magnetization)],
        }
    }

    /// Add a custom observable to the set.
    pub fn add(&mut self, obs: Box<dyn Observable<H>>) {
        self.observables.push(obs);
    }

    /// Iterate over all observables.
    pub fn iter(&self) -> impl Iterator<Item = &dyn Observable<H>> {
        self.observables.iter().map(|o| o.as_ref())
    }
}

impl<H: Hamiltonian + Measurable> Default for DefaultObservableSet<H> {
    fn default() -> Self {
        Self::new()
    }
}
