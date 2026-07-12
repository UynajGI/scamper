//! Common observables infrastructure: traits, sets, and moment specs.

use crate::lattice::interaction::{Hamiltonian, Measurable};
use crate::lattice::state::System;
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
            observables: vec![
                Box::new(crate::observables::energy::TotalEnergy),
                Box::new(crate::observables::magnetization::Magnetization),
            ],
        }
    }
}

impl<H: Hamiltonian> DefaultObservableSet<H> {
    pub fn energy_only() -> Self {
        Self {
            observables: vec![Box::new(crate::observables::energy::TotalEnergy)],
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
