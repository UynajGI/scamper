//! Swendsen-Wang cluster-update kernel.

use crate::algorithms::common::{checked_probability, Algorithm, SimulationPhase};
use crate::core::cache::BatchEnergyPatch;
use crate::core::r#move::BatchSpinMove;
use crate::core::trial::TrialEvaluator;
use crate::lattice::interaction::{ClusterAuxiliary, ClusterModel, Hamiltonian};
use crate::lattice::state::System;
use rand::{Rng, RngExt};

#[derive(Debug, Clone, Default)]
pub struct SWCore {
    parent: Vec<usize>,
    rank: Vec<u8>,
    root_auxiliary: Vec<Option<ClusterAuxiliary>>,
    movement: BatchSpinMove,
    patch: BatchEnergyPatch,
}

impl SWCore {
    pub fn new() -> Self {
        Self {
            parent: Vec::new(),
            rank: Vec::new(),
            root_auxiliary: Vec::new(),
            movement: BatchSpinMove::default(),
            patch: BatchEnergyPatch {
                delta_energy: 0.0,
                workspace: crate::core::cache::BatchEnergyWorkspace::new(),
            },
        }
    }

    fn reset_union_find(&mut self, n_sites: usize) {
        self.parent.clear();
        self.parent.extend(0..n_sites);
        self.rank.clear();
        self.rank.resize(n_sites, 0);
        self.root_auxiliary.clear();
        self.root_auxiliary.resize(n_sites, None);
    }

    fn find(&mut self, site: usize) -> usize {
        let mut root = site;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut current = site;
        while self.parent[current] != root {
            let next = self.parent[current];
            self.parent[current] = root;
            current = next;
        }
        root
    }

    fn union(&mut self, left: usize, right: usize) {
        let left_root = self.find(left);
        let right_root = self.find(right);
        if left_root == right_root {
            return;
        }
        match self.rank[left_root].cmp(&self.rank[right_root]) {
            std::cmp::Ordering::Less => self.parent[left_root] = right_root,
            std::cmp::Ordering::Greater => self.parent[right_root] = left_root,
            std::cmp::Ordering::Equal => {
                self.parent[right_root] = left_root;
                self.rank[left_root] += 1;
            }
        }
    }
}

impl<H: Hamiltonian + ClusterModel> Algorithm<H> for SWCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        self.reset_union_find(n_sites);
        let bond_auxiliary = model.sw_bond_auxiliary(rng);

        // Physical edges are visited exactly once.
        for edge in &system.lattice.edges {
            let left_base = edge.source * spin_dim;
            let right_base = edge.target * spin_dim;
            let probability = checked_probability(
                model.cluster_bond_probability(
                    &system.spins[left_base..left_base + spin_dim],
                    &system.spins[right_base..right_base + spin_dim],
                    edge,
                    &bond_auxiliary,
                    system.beta,
                ),
                "Swendsen-Wang",
            );
            if rng.random::<f64>() < probability {
                self.union(edge.source, edge.target);
            }
        }

        // Every root receives an independent target/reflection decision.
        for site in 0..n_sites {
            let root = self.find(site);
            if self.root_auxiliary[root].is_none() {
                self.root_auxiliary[root] = Some(model.sw_cluster_auxiliary(
                    system.spin_at(site, spin_dim),
                    &bond_auxiliary,
                    rng,
                ));
            }
        }

        self.movement.reset(spin_dim);
        for site in 0..n_sites {
            let root = self.find(site);
            let auxiliary = self.root_auxiliary[root]
                .as_ref()
                .expect("SW root transformation must be initialized");
            let transformed =
                model.transform_cluster_spin(system.spin_at(site, spin_dim), auxiliary);
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
        "Swendsen-Wang"
    }
}
