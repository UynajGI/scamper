//! Exact small-system Ising density of states for validation.

use crate::generalized::{DiscreteAxis, GeneralizedError, LogDensityOfStates};
use crate::lattice::graph::CsrLattice;
use crate::lattice::interaction::Hamiltonian;
use crate::lattice::models::IsingModel;

/// Exact energy levels and integer degeneracies of a small Ising graph.
#[derive(Debug, Clone, PartialEq)]
pub struct ExactIsingDensityOfStates {
    energies: Vec<f64>,
    degeneracies: Vec<u64>,
}

impl ExactIsingDensityOfStates {
    #[inline]
    pub fn energies(&self) -> &[f64] {
        &self.energies
    }

    #[inline]
    pub fn degeneracies(&self) -> &[u64] {
        &self.degeneracies
    }

    #[inline]
    pub fn states(&self) -> u64 {
        self.degeneracies.iter().sum()
    }

    pub fn axis(&self) -> Result<DiscreteAxis, GeneralizedError> {
        DiscreteAxis::new(self.energies.clone())
    }

    pub fn log_density(&self) -> Result<LogDensityOfStates, GeneralizedError> {
        LogDensityOfStates::from_values(
            self.degeneracies
                .iter()
                .map(|&degeneracy| (degeneracy as f64).ln())
                .collect(),
            vec![true; self.degeneracies.len()],
        )
    }
}

/// Enumerate all `2^N` configurations of an Ising graph.
///
/// The explicit limit keeps accidental exponential jobs out of production
/// paths. The intended use is exact regression testing and small reference
/// calculations, not simulation of large systems.
pub fn enumerate_ising_density_of_states(
    lattice: &CsrLattice,
    model: &IsingModel,
) -> Result<ExactIsingDensityOfStates, GeneralizedError> {
    lattice
        .validate()
        .map_err(|detail| GeneralizedError::new(format!("invalid lattice: {detail}")))?;
    if lattice.n_sites > 24 {
        return Err(GeneralizedError::new(
            "exact Ising enumeration is limited to 24 sites",
        ));
    }
    let states = 1_u64
        .checked_shl(lattice.n_sites as u32)
        .ok_or_else(|| GeneralizedError::new("Ising state-count overflow"))?;
    let mut energies = Vec::with_capacity(states as usize);
    let mut spins = vec![-1.0; lattice.n_sites];
    for state in 0..states {
        for (site, spin) in spins.iter_mut().enumerate() {
            *spin = if state & (1_u64 << site) == 0 {
                -1.0
            } else {
                1.0
            };
        }
        energies.push(model.compute_total_energy(&spins, lattice, 0.0));
    }
    energies.sort_by(f64::total_cmp);

    // Bound forward-summation roundoff instead of using a loose
    // relative tolerance that could merge physically distinct weighted
    // energy levels. Each Ising bond term has magnitude |J w_e|.
    let terms = lattice.n_edges().max(1) as f64;
    let sum_abs = lattice
        .edges
        .iter()
        .map(|edge| (model.j * edge.weight).abs())
        .sum::<f64>();
    let gamma = terms * f64::EPSILON / (1.0 - terms * f64::EPSILON);
    let grouping_tolerance = 8.0 * gamma * sum_abs.max(1.0);

    let mut levels: Vec<f64> = Vec::new();
    let mut degeneracies: Vec<u64> = Vec::new();
    for energy in energies {
        if levels
            .last()
            .is_some_and(|&previous| (previous - energy).abs() <= grouping_tolerance)
        {
            let last = degeneracies
                .last_mut()
                .expect("an existing energy level has a degeneracy");
            *last += 1;
        } else {
            levels.push(energy);
            degeneracies.push(1_u64);
        }
    }
    Ok(ExactIsingDensityOfStates {
        energies: levels,
        degeneracies,
    })
}
