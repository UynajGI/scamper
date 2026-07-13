//! Wolff cluster-update kernel.

use crate::algorithms::common::{checked_probability, Algorithm, SimulationPhase};
use crate::core::cache::BatchEnergyPatch;
use crate::core::r#move::BatchSpinMove;
use crate::core::trial::TrialEvaluator;
use crate::lattice::interaction::{ClusterModel, Hamiltonian};
use crate::lattice::state::System;
use rand::{Rng, RngExt};

#[derive(Debug, Clone, Default)]
pub struct WolffCore {
    visit_stamp: Vec<u32>,
    generation: u32,
    stack: Vec<usize>,
    members: Vec<usize>,
    movement: BatchSpinMove,
    patch: BatchEnergyPatch,
    last_cluster_size: usize,
}

impl WolffCore {
    pub fn new() -> Self {
        Self {
            visit_stamp: Vec::new(),
            generation: 0,
            stack: Vec::new(),
            members: Vec::new(),
            movement: BatchSpinMove::default(),
            patch: BatchEnergyPatch {
                delta_energy: 0.0,
                workspace: crate::core::cache::BatchEnergyWorkspace::new(),
            },
            last_cluster_size: 0,
        }
    }

    fn begin_cluster(&mut self, n_sites: usize) {
        if self.visit_stamp.len() != n_sites {
            self.visit_stamp.resize(n_sites, 0);
        }
        if self.generation == u32::MAX {
            self.visit_stamp.fill(0);
            self.generation = 1;
        } else {
            self.generation += 1;
            if self.generation == 0 {
                self.generation = 1;
            }
        }
        self.stack.clear();
        self.members.clear();
    }

    #[inline]
    fn contains(&self, site: usize) -> bool {
        self.visit_stamp[site] == self.generation
    }

    /// Number of sites in the most recently completed cluster.
    #[inline]
    pub const fn last_cluster_size(&self) -> usize {
        self.last_cluster_size
    }

    #[inline]
    fn insert(&mut self, site: usize) {
        self.visit_stamp[site] = self.generation;
        self.stack.push(site);
    }
}

impl<H: Hamiltonian + ClusterModel> Algorithm<H> for WolffCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        if n_sites == 0 {
            return;
        }
        let spin_dim = model.spin_dim();
        self.begin_cluster(n_sites);

        let seed = rng.random_range(0..n_sites);
        let auxiliary = model.wolff_auxiliary(system.spin_at(seed, spin_dim), rng);
        self.insert(seed);

        while let Some(site) = self.stack.pop() {
            self.members.push(site);
            let left_base = site * spin_dim;
            let left = &system.spins[left_base..left_base + spin_dim];
            for (neighbor, edge_id) in system.lattice.incidences(site) {
                if self.contains(neighbor) {
                    continue;
                }
                let right_base = neighbor * spin_dim;
                let right = &system.spins[right_base..right_base + spin_dim];
                let probability = checked_probability(
                    model.cluster_bond_probability(
                        left,
                        right,
                        &system.lattice.edges[edge_id],
                        &auxiliary,
                        system.beta,
                    ),
                    "Wolff",
                );
                if rng.random::<f64>() < probability {
                    self.insert(neighbor);
                }
            }
        }

        self.last_cluster_size = self.members.len();
        self.movement.reset(spin_dim);
        for &site in &self.members {
            let transformed =
                model.transform_cluster_spin(system.spin_at(site, spin_dim), &auxiliary);
            assert_eq!(
                transformed.len(),
                spin_dim,
                "cluster transform dimension mismatch"
            );
            self.movement.push(site, &transformed);
        }
        system.evaluate_trial(model, &self.movement, &mut self.patch);
        <System as TrialEvaluator<H, BatchSpinMove>>::commit_trial(
            system,
            &self.movement,
            &self.patch,
        );
    }

    fn name(&self) -> &'static str {
        "Wolff"
    }
}
