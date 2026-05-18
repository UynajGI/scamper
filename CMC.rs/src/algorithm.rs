//! Monte Carlo algorithms — sweep strategies.

use crate::model::Model;
use crate::proposal::ProposalStrategy;
use crate::system::System;
use rand::Rng;
use rand::RngExt;

/// Algorithm trait. One `sweep` = one full pass over the system.
///
/// The algorithm is responsible for updating both `system.spins` and `system.energy`.
pub trait Algorithm<M: Model>: Send {
    fn sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng);

    fn name(&self) -> &'static str {
        "Unknown"
    }
}

// ── Metropolis ──────────────────────────────────────────────

/// Metropolis algorithm with pluggable proposal strategy.
///
/// In one sweep, visits every site once in random order, proposes a new spin,
/// computes ΔE, and accepts with probability `min(1, exp(-β ΔE))`.
#[derive(Debug, Clone)]
pub struct MetropolisCore<S = crate::proposal::StandardStrategy> {
    strategy: S,
}

impl MetropolisCore {
    pub fn new() -> Self {
        Self {
            strategy: Default::default(),
        }
    }
}

impl Default for MetropolisCore {
    fn default() -> Self {
        Self::new()
    }
}

impl<S: Default> MetropolisCore<S> {
    pub fn with_strategy(strategy: S) -> Self {
        Self { strategy }
    }
}

impl<M, S> Algorithm<M> for MetropolisCore<S>
where
    M: Model,
    S: ProposalStrategy<M> + Clone,
{
    fn sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let beta = model.beta();

        // Random visit order
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for &site in &order {
            let old_energy = model.local_energy(
                &system.spins,
                &system.lattice,
                site,
                system.spin_at(site, sd),
            );

            let proposed = self.strategy.propose(model, system, site, rng);
            let new_energy = model.local_energy(&system.spins, &system.lattice, site, &proposed);

            let delta_e = new_energy - old_energy;
            let accepted = delta_e <= 0.0 || rng.random::<f64>() < (-beta * delta_e).exp();

            if accepted {
                system.spin_at_mut(site, sd).copy_from_slice(&proposed);
                system.energy += delta_e;
            }
        }

        self.strategy.adapt_after_sweep(model);
    }

    fn name(&self) -> &'static str {
        "Metropolis"
    }
}

// ── Wolff Cluster ───────────────────────────────────────────

/// Wolff single-cluster algorithm.
///
/// Grows a cluster from a random seed site, connecting aligned neighbors
/// with FK bond probability, then flips the entire cluster.
#[derive(Debug, Clone, Default)]
pub struct WolffCore;

impl WolffCore {
    pub fn new() -> Self {
        Self
    }
}

impl<M: Model> Algorithm<M> for WolffCore {
    fn sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let p_add = model.fk_bond_probability();

        // Pick random seed
        let seed = rng.random_range(0..n);
        let seed_spin = system.spin_at(seed, sd).to_vec();
        let opposite = model.opposite_spin(seed_spin[0], rng);

        // BFS cluster growth
        let mut cluster = vec![false; n];
        let mut stack = vec![seed];
        cluster[seed] = true;

        while let Some(site) = stack.pop() {
            for nb in &system.lattice.sites[site] {
                if cluster[nb.target] {
                    continue;
                }
                // Check spin alignment
                let nb_spin = system.spins[nb.target];
                if (nb_spin - seed_spin[0]).abs() > 1e-10 {
                    continue;
                }
                // FK bond
                if rng.random::<f64>() < p_add {
                    cluster[nb.target] = true;
                    stack.push(nb.target);
                }
            }
        }

        // Flip cluster and update energy
        for (site, &in_cluster) in cluster.iter().enumerate() {
            if in_cluster {
                let old_local = model.local_energy(
                    &system.spins,
                    &system.lattice,
                    site,
                    system.spin_at(site, sd),
                );
                let new_spin = vec![opposite];
                let new_local = model.local_energy(&system.spins, &system.lattice, site, &new_spin);
                system.energy += new_local - old_local;
                system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
            }
        }
    }

    fn name(&self) -> &'static str {
        "Wolff"
    }
}

// ── Swendsen-Wang ───────────────────────────────────────────

/// Swendsen-Wang cluster algorithm.
///
/// Visits all bonds, connects aligned neighbors with FK bond probability,
/// identifies clusters via union-find, then flips each cluster with 50% probability.
#[derive(Debug, Clone, Default)]
pub struct SWCore;

impl SWCore {
    pub fn new() -> Self {
        Self
    }
}

impl<M: Model> Algorithm<M> for SWCore {
    fn sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let p_add = model.fk_bond_probability();

        // Union-find
        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank = vec![0usize; n];

        fn find(parent: &mut [usize], x: usize) -> usize {
            let mut px = parent[x];
            while parent[px] != px {
                px = parent[px];
            }
            let root = px;
            let mut x = x;
            while parent[x] != root {
                let nx = parent[x];
                parent[x] = root;
                x = nx;
            }
            root
        }

        fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra == rb {
                return;
            }
            match rank[ra].cmp(&rank[rb]) {
                std::cmp::Ordering::Less => parent[ra] = rb,
                std::cmp::Ordering::Greater => parent[rb] = ra,
                std::cmp::Ordering::Equal => {
                    parent[rb] = ra;
                    rank[ra] += 1;
                }
            }
        }

        // Visit all bonds
        for i in 0..n {
            for nb in &system.lattice.sites[i] {
                if i >= nb.target {
                    continue; // process each bond once
                }
                if (system.spins[i] - system.spins[nb.target]).abs() > 1e-10 {
                    continue;
                }
                if rng.random::<f64>() < p_add {
                    union(&mut parent, &mut rank, i, nb.target);
                }
            }
        }

        // Decide flip for each cluster root
        let mut flip_root = vec![false; n];
        for i in 0..n {
            if parent[i] == i {
                flip_root[i] = rng.random::<bool>();
            }
        }

        // Apply flips
        let opposite = model.random_cluster_spin(rng);
        for site in 0..n {
            let root = find(&mut parent, site);
            if flip_root[root] {
                let old_local = model.local_energy(
                    &system.spins,
                    &system.lattice,
                    site,
                    system.spin_at(site, sd),
                );
                let new_spin = vec![opposite];
                let new_local = model.local_energy(&system.spins, &system.lattice, site, &new_spin);
                system.energy += new_local - old_local;
                system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
            }
        }
    }

    fn name(&self) -> &'static str {
        "Swendsen-Wang"
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::model::IsingModel;
    use rand::SeedableRng;

    fn make_rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn test_metropolis_single_site() {
        let lattice = build_chain(1, false);
        let mut system = System::new(lattice, 1, 1.0);
        system.energy = 0.0;

        let model = IsingModel::new(1.0, 1.0);
        let mut algo = MetropolisCore::new();
        let mut rng = make_rng();

        for _ in 0..100 {
            algo.sweep(&mut system, &model, &mut rng);
        }
        assert!((system.energy - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_metropolis_cools_to_ground() {
        let lattice = build_chain(4, true);
        let mut system = System::new(lattice.clone(), 1, 1.0);
        let model = IsingModel::new(1.0, 10.0); // beta = 10, very cold
        system.energy = model.compute_total_energy(&system.spins, &system.lattice);

        let mut algo = MetropolisCore::new();
        let mut rng = make_rng();

        for _ in 0..500 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        assert!(system.energy < -3.5);
        let all_same = system.spins.iter().all(|&s| s == system.spins[0]);
        assert!(all_same);
    }

    #[test]
    fn test_wolff_preserves_energy_sign() {
        // At low T, Wolff should maintain ferromagnetic order
        let lattice = build_chain(8, true);
        let mut system = System::new(lattice.clone(), 1, 1.0);
        let model = IsingModel::new(1.0, 5.0); // beta=5, cold
        system.energy = model.compute_total_energy(&system.spins, &system.lattice);

        let mut algo = WolffCore::new();
        let mut rng = make_rng();

        for _ in 0..100 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // After Wolff at low T, all spins should be aligned (cluster covers system)
        let all_same = system.spins.iter().all(|&s| s == system.spins[0]);
        assert!(all_same);
        assert!(system.energy < -7.0); // all aligned, 8 bonds × -J = -8
    }

    #[test]
    fn test_sw_preserves_energy_sign() {
        let lattice = build_chain(8, true);
        let mut system = System::new(lattice.clone(), 1, 1.0);
        let model = IsingModel::new(1.0, 5.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice);

        let mut algo = SWCore::new();
        let mut rng = make_rng();

        for _ in 0..100 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // SW at high β should converge to all aligned
        let all_same = system.spins.iter().all(|&s| s == system.spins[0]);
        assert!(all_same);
        assert!(system.energy < -7.0);
    }
}
