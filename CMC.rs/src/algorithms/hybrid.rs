//! Statically composed hybrid update kernel.

use crate::algorithms::common::{Algorithm, SimulationPhase};
use crate::lattice::interaction::Hamiltonian;
use crate::lattice::state::System;
use rand::Rng;

/// Statically composed hybrid update without trait-object overhead.
#[derive(Debug, Clone)]
pub struct HybridCore<A, B> {
    pub first: A,
    pub second: B,
    pub first_repetitions: usize,
    pub second_repetitions: usize,
}

impl<A, B> HybridCore<A, B> {
    pub fn new(first: A, second: B) -> Self {
        Self {
            first,
            second,
            first_repetitions: 1,
            second_repetitions: 1,
        }
    }

    pub fn repetitions(mut self, first: usize, second: usize) -> Self {
        self.first_repetitions = first;
        self.second_repetitions = second;
        self
    }
}

impl<H, A, B> Algorithm<H> for HybridCore<A, B>
where
    H: Hamiltonian,
    A: Algorithm<H>,
    B: Algorithm<H>,
{
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        for _ in 0..self.first_repetitions {
            self.first.sweep_with_phase(system, model, rng, phase);
        }
        for _ in 0..self.second_repetitions {
            self.second.sweep_with_phase(system, model, rng, phase);
        }
    }

    fn name(&self) -> &'static str {
        "Hybrid"
    }
}
