//! Pluggable observables integrated with Carlo.rs measurement accumulation.

use crate::hamiltonian::{Hamiltonian, Measurable};
use crate::system::System;
use carlo_rs::Context;
use rand::{Rng, SeedableRng};

/// One derived raw moment of an observable.
#[derive(Debug, Clone, Copy)]
pub struct MomentSpec {
    pub name: &'static str,
    pub order: u32,
}

/// Scalar observable.
pub trait Observable<H: Hamiltonian>: Send {
    fn name(&self) -> &str;
    fn measure(&self, system: &System, model: &H) -> f64;

    /// Raw moments recorded from the same sample.  Names are explicit, so no
    /// string matching is needed in `ClassicalMC`.
    fn moments(&self) -> &[MomentSpec] {
        &[]
    }
}

pub struct TotalEnergy;

impl<H: Hamiltonian> Observable<H> for TotalEnergy {
    fn name(&self) -> &str {
        "Energy"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy
    }

    fn moments(&self) -> &[MomentSpec] {
        static MOMENTS: [MomentSpec; 1] = [MomentSpec {
            name: "E2",
            order: 2,
        }];
        &MOMENTS
    }
}

pub struct EnergyPerSite;

impl<H: Hamiltonian> Observable<H> for EnergyPerSite {
    fn name(&self) -> &str {
        "EnergyPerSite"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy / system.n_sites() as f64
    }
}

pub struct Magnetization;

impl<H: Hamiltonian + Measurable> Observable<H> for Magnetization {
    fn name(&self) -> &str {
        "Magnetization"
    }

    fn measure(&self, system: &System, model: &H) -> f64 {
        model.magnetization(&system.spins)
    }

    fn moments(&self) -> &[MomentSpec] {
        static MOMENTS: [MomentSpec; 2] = [
            MomentSpec {
                name: "M2",
                order: 2,
            },
            MomentSpec {
                name: "M4",
                order: 4,
            },
        ];
        &MOMENTS
    }
}

/// Collection abstraction used by the generic `ClassicalMC<H,A,O>` wrapper.
pub trait ObservableSet<H: Hamiltonian>: Send {
    fn measure_all<R: Rng + SeedableRng>(
        &self,
        system: &System,
        model: &H,
        context: &mut Context<R>,
    );
}

/// Runtime-extensible scalar observable collection.
pub struct DefaultObservableSet<H: Hamiltonian> {
    observables: Vec<Box<dyn Observable<H>>>,
}

impl<H: Hamiltonian + Measurable> DefaultObservableSet<H> {
    pub fn new() -> Self {
        Self {
            observables: vec![Box::new(TotalEnergy), Box::new(Magnetization)],
        }
    }
}

impl<H: Hamiltonian> DefaultObservableSet<H> {
    pub fn energy_only() -> Self {
        Self {
            observables: vec![Box::new(TotalEnergy)],
        }
    }

    pub fn empty() -> Self {
        Self {
            observables: Vec::new(),
        }
    }

    pub fn add(&mut self, observable: Box<dyn Observable<H>>) {
        self.observables.push(observable);
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn Observable<H>> {
        self.observables.iter().map(Box::as_ref)
    }
}

impl<H: Hamiltonian + Measurable> Default for DefaultObservableSet<H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<H: Hamiltonian> ObservableSet<H> for DefaultObservableSet<H> {
    fn measure_all<R: Rng + SeedableRng>(
        &self,
        system: &System,
        model: &H,
        context: &mut Context<R>,
    ) {
        for observable in &self.observables {
            let value = observable.measure(system, model);
            context.measure(observable.name(), value);
            for moment in observable.moments() {
                context.measure(moment.name, value.powi(moment.order as i32));
            }
        }
    }
}

/// Useful for simulations that measure inside a custom algorithm.
#[derive(Debug, Clone, Copy, Default)]
pub struct EmptyObservableSet;

impl<H: Hamiltonian> ObservableSet<H> for EmptyObservableSet {
    fn measure_all<R: Rng + SeedableRng>(
        &self,
        _system: &System,
        _model: &H,
        _context: &mut Context<R>,
    ) {
    }
}
