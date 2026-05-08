//! Wolff cluster algorithm.

use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::VecDeque;

/// Wolff cluster algorithm.
pub struct WolffCore<MC: ModelMC> {
    model: MC,
    snapshot_interval: Option<u64>,
}

impl<MC: ModelMC> WolffCore<MC> {
    pub fn new(model: MC) -> Self {
        WolffCore {
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

impl<MC: ModelMC> MonteCarlo for WolffCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        let p_add = self.model.fk_bond_probability();

        // Wolff: grow one cluster and flip it
        let seed = ctx.rng.random_range(0..n);
        let seed_spin = self.model.spins()[seed];
        let mut cluster = Vec::new();
        let mut stack = VecDeque::new();
        // Track considered bonds to avoid double-counting bidirectional bonds
        let mut considered = std::collections::HashSet::new();

        cluster.push(seed);
        stack.push_back(seed);

        while let Some(site) = stack.pop_front() {
            for neighbor in &self.model.lattice().sites[site] {
                let target = neighbor.target;
                let nb_spin = self.model.spins()[target];
                if nb_spin != seed_spin || cluster.contains(&target) {
                    continue;
                }
                // Each physical bond appears twice (a→b and b→a).
                // Use canonical ordering to consider each bond only once.
                let bond = if site < target { (site, target) } else { (target, site) };
                if considered.insert(bond) && ctx.rng.random::<f64>() < p_add {
                    cluster.push(target);
                    stack.push_back(target);
                }
            }
        }

        // Flip entire cluster
        let new_spin = self.model.opposite_spin(seed_spin, &mut ctx.rng);
        for &site in &cluster {
            self.model.spins_mut()[site] = new_spin;
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
        "WolffCore"
    }
}

impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for WolffCore<MC> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        Ok(WolffCore {
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
    fn test_wolff_sweep_does_not_crash() {
        let lattice = build_chain(8, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        let mut core = WolffCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        core.sweep(&mut ctx);
    }

    #[test]
    fn test_wolff_energy_decreases() {
        let lattice = build_chain(16, true);
        // Start with alternating spins (higher energy)
        let mut model = IsingModel::new(lattice, 0.5, 1.0);
        for i in 0..16 {
            if i % 2 == 1 {
                model.spins_mut()[i] = -1.0;
            }
        }
        let mut core = WolffCore::new(model);
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
    fn test_wolff_cluster_flip() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        let mut core = WolffCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        // Initial state: all spins up
        assert!(core.model.spins().iter().all(|&s| s > 0.0));

        core.sweep(&mut ctx);

        // After Wolff sweep, some spins should have flipped
        let n_up = core.model.spins().iter().filter(|&&s| s > 0.0).count();
        assert!(n_up < 4);
    }
}
