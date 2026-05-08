//! Swendsen-Wang cluster algorithm.

use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::{HashMap, HashSet};

/// Union-Find data structure for cluster identification.
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx != ry {
            self.parent[rx] = ry;
        }
    }
}

/// Swendsen-Wang cluster algorithm.
pub struct SWCore<MC: ModelMC> {
    model: MC,
    snapshot_interval: Option<u64>,
}

impl<MC: ModelMC> SWCore<MC> {
    pub fn new(model: MC) -> Self {
        SWCore {
            model,
            snapshot_interval: None,
        }
    }

    /// Set the snapshot recording interval (in sweeps).
    /// When set, every `interval` sweeps during measurement, the current
    /// spin configuration is recorded via `ctx.measure_array("Snapshot", ...)`.
    pub fn with_snapshot_interval(mut self, interval: u64) -> Self {
        self.snapshot_interval = Some(interval);
        self
    }

    pub fn model(&self) -> &MC {
        &self.model
    }

    pub fn model_mut(&mut self) -> &mut MC {
        &mut self.model
    }
}

impl<MC: ModelMC> MonteCarlo for SWCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        let p_add = self.model.fk_bond_probability();

        // Step 1: Build FK clusters
        // With bidirectional bonds, each physical bond appears twice (a→b and b→a).
        // We only consider each bond once by checking site < neighbor.target.
        let mut uf = UnionFind::new(n);

        for site in 0..n {
            let site_spin = self.model.spins()[site];
            for neighbor in &self.model.lattice().sites[site] {
                let target = neighbor.target;
                if site < target
                    && site_spin == self.model.spins()[target]
                    && ctx.rng.random::<f64>() < p_add
                {
                    uf.union(site, target);
                }
            }
        }

        // Step 2: Assign random spin to each cluster
        let mut cluster_roots = HashSet::new();
        for i in 0..n {
            cluster_roots.insert(uf.find(i));
        }

        let cluster_spins: HashMap<usize, f64> = cluster_roots
            .into_iter()
            .map(|root| {
                let spin = self.model.random_cluster_spin(&mut ctx.rng);
                (root, spin)
            })
            .collect();

        // Step 3: Apply cluster spins
        for i in 0..n {
            let root = uf.find(i);
            self.model.spins_mut()[i] = cluster_spins[&root];
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        let magnetization = self.model.magnetization();
        ctx.measure("Energy", energy);
        ctx.measure("Energy_Squared", energy * energy);
        ctx.measure("Magnetization", magnetization);
        ctx.measure("Magnetization_Squared", magnetization * magnetization);

        if let Some(interval) = self.snapshot_interval {
            if ctx.sweep_count() % interval == 0 {
                ctx.measure_array("Snapshot", &self.model.snapshot());
            }
        }
    }

    fn name(&self) -> &'static str {
        "SWCore"
    }
}

impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for SWCore<MC> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        Ok(SWCore {
            model: mc,
            snapshot_interval: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::models::IsingModel;
    use rand::SeedableRng;

    #[test]
    fn test_sw_sweep_does_not_crash() {
        let lattice = build_chain(8, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        let mut core = SWCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        core.sweep(&mut ctx);
    }

    #[test]
    fn test_sw_energy_decreases() {
        let lattice = build_chain(16, true);
        // Start with alternating spins (higher energy)
        let mut model = IsingModel::new(lattice, 0.5, 1.0);
        for i in 0..16 {
            if i % 2 == 1 {
                model.spins_mut()[i] = -1.0;
            }
        }
        let mut core = SWCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        let initial_energy = core.model.total_energy();

        for _ in 0..50 {
            core.sweep(&mut ctx);
        }

        let final_energy = core.model.total_energy();
        assert!(final_energy <= initial_energy + 1e-10);
    }

    #[test]
    fn test_sw_cluster_structure() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        let mut core = SWCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        // All spins up initially
        assert!(core.model.spins().iter().all(|&s| s > 0.0));

        core.sweep(&mut ctx);

        // SW can produce any cluster configuration at finite temperature
        // Just verify it doesn't crash and spins are still +/-1
        for &s in core.model.spins() {
            assert!((s - 1.0).abs() < 1e-10 || (s + 1.0).abs() < 1e-10);
        }
    }
}
