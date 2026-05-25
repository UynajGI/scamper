//! Monte Carlo algorithms — sweep strategies.

use crate::hamiltonian::{ClusterModel, ContinuousHeatBathable, Hamiltonian, HeatBathable, Proposable};
use crate::proposal::ProposalStrategy;
use crate::system::System;
use rand::Rng;
use rand::RngExt;
use smallvec::{smallvec, SmallVec};

/// Algorithm trait. One `sweep` = one full pass over the system.
///
/// The algorithm is responsible for updating both `system.spins` and `system.energy`.
pub trait Algorithm<H: Hamiltonian>: Send {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng);

    fn name(&self) -> &'static str {
        "Unknown"
    }
}

/// Dot product of two equal-length slices.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
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

impl<H, S> Algorithm<H> for MetropolisCore<S>
where
    H: Hamiltonian + Proposable,
    S: ProposalStrategy<H> + Clone,
{
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let beta = system.beta;

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
                beta,
                system.spin_at(site, sd),
            );

            let proposed = self.strategy.propose(model, system, site, rng);
            let new_energy =
                model.local_energy(&system.spins, &system.lattice, site, beta, &proposed);

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
/// For scalar spins (Ising, Potts): grows a cluster from a random seed,
/// connects aligned neighbors with FK bond probability, flips the cluster.
///
/// For vector spins (XY, Heisenberg): uses Wolff embedding — projects spins
/// onto a random direction, builds FK clusters in the projected Ising variables,
/// then reflects each spin across the plane perpendicular to the embedding direction.
#[derive(Debug, Clone, Default)]
pub struct WolffCore;

impl WolffCore {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + ClusterModel> Algorithm<H> for WolffCore {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let beta = system.beta;
        let p_add = model.fk_bond_probability(beta);

        if sd == 1 {
            // ── scalar fast path (Ising, Potts) ──
            let seed = rng.random_range(0..n);
            let seed_spin = system.spins[seed];
            let opposite = model.opposite_spin(seed_spin, rng);

            let mut cluster = vec![false; n];
            let mut stack = vec![seed];
            cluster[seed] = true;

            while let Some(site) = stack.pop() {
                for &nb in system.lattice.neighbors(site) {
                    if cluster[nb] {
                        continue;
                    }
                    let nb_spin = system.spins[nb];
                    if (nb_spin - seed_spin).abs() > 1e-10 {
                        continue;
                    }
                    if rng.random::<f64>() < p_add {
                        cluster[nb] = true;
                        stack.push(nb);
                    }
                }
            }

            for (site, &in_cluster) in cluster.iter().enumerate() {
                if in_cluster {
                    let old_local = model.local_energy(
                        &system.spins,
                        &system.lattice,
                        site,
                        beta,
                        system.spin_at(site, sd),
                    );
                    let new_spin: SmallVec<[f64; 3]> = smallvec![opposite];
                    let new_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, &new_spin);
                    system.energy += new_local - old_local;
                    system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
                }
            }
        } else {
            // ── embedding path (XY, Heisenberg) ──
            let direction = model.embedding_direction(rng);

            let seed = rng.random_range(0..n);
            let seed_spin = system.spin_at(seed, sd);
            let seed_proj: f64 = dot(seed_spin, &direction);

            // BFS cluster based on projected spins
            let mut cluster = vec![false; n];
            let mut stack = vec![seed];
            cluster[seed] = true;

            while let Some(site) = stack.pop() {
                for &nb in system.lattice.neighbors(site) {
                    if cluster[nb] {
                        continue;
                    }
                    let nb_spin = system.spin_at(nb, sd);
                    let nb_proj: f64 = dot(nb_spin, &direction);
                    if seed_proj * nb_proj <= 0.0 {
                        continue;
                    }
                    if rng.random::<f64>() < p_add {
                        cluster[nb] = true;
                        stack.push(nb);
                    }
                }
            }

            // Reflect spins across plane ⟂ direction
            for (site, &in_cluster) in cluster.iter().enumerate() {
                if in_cluster {
                    let spin = system.spin_at(site, sd);
                    let old_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, spin);
                    let mut new_spin = SmallVec::<[f64; 3]>::from_slice(spin);
                    model.reflect(&mut new_spin, &direction);
                    let new_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, &new_spin);
                    system.energy += new_local - old_local;
                    system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
                }
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
/// For scalar spins (Ising, Potts): visits all bonds, connects aligned neighbors
/// with FK bond probability, identifies clusters via union-find, flips each
/// cluster with 50% probability.
///
/// For vector spins (XY, Heisenberg): uses SW embedding — projects spins onto
/// a random direction, builds FK clusters, flips clusters via spin reflection.
#[derive(Debug, Clone, Default)]
pub struct SWCore;

impl SWCore {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + ClusterModel> Algorithm<H> for SWCore {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let beta = system.beta;
        let p_add = model.fk_bond_probability(beta);

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

        if sd == 1 {
            // ── scalar fast path (Ising, Potts) ──
            for i in 0..n {
                for &nb in system.lattice.neighbors(i) {
                    if i >= nb {
                        continue;
                    }
                    if (system.spins[i] - system.spins[nb]).abs() > 1e-10 {
                        continue;
                    }
                    if rng.random::<f64>() < p_add {
                        union(&mut parent, &mut rank, i, nb);
                    }
                }
            }

            let mut flip_root = vec![false; n];
            for i in 0..n {
                if parent[i] == i {
                    flip_root[i] = rng.random::<bool>();
                }
            }

            let new_cluster_spin = model.random_cluster_spin(rng);
            for site in 0..n {
                let root = find(&mut parent, site);
                if flip_root[root] {
                    let old_local = model.local_energy(
                        &system.spins,
                        &system.lattice,
                        site,
                        beta,
                        system.spin_at(site, sd),
                    );
                    let new_spin: SmallVec<[f64; 3]> = smallvec![new_cluster_spin];
                    let new_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, &new_spin);
                    system.energy += new_local - old_local;
                    system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
                }
            }
        } else {
            // ── embedding path (XY, Heisenberg) ──
            let direction = model.embedding_direction(rng);

            for i in 0..n {
                for &nb in system.lattice.neighbors(i) {
                    if i >= nb {
                        continue;
                    }
                    let proj_i: f64 = dot(system.spin_at(i, sd), &direction);
                    let proj_nb: f64 = dot(system.spin_at(nb, sd), &direction);
                    if proj_i * proj_nb <= 0.0 {
                        continue;
                    }
                    if rng.random::<f64>() < p_add {
                        union(&mut parent, &mut rank, i, nb);
                    }
                }
            }

            let mut flip_root = vec![false; n];
            for i in 0..n {
                if parent[i] == i {
                    flip_root[i] = rng.random::<bool>();
                }
            }

            for site in 0..n {
                let root = find(&mut parent, site);
                if flip_root[root] {
                    let spin = system.spin_at(site, sd);
                    let old_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, spin);
                    let mut new_spin = SmallVec::<[f64; 3]>::from_slice(spin);
                    model.reflect(&mut new_spin, &direction);
                    let new_local =
                        model.local_energy(&system.spins, &system.lattice, site, beta, &new_spin);
                    system.energy += new_local - old_local;
                    system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
                }
            }
        }
    }

    fn name(&self) -> &'static str {
        "Swendsen-Wang"
    }
}

// ── Microcanonical Over-Relaxation ───────────────────────────

/// Reflect a spin across the unit vector of the local field.
///
/// `s_new = 2 (s·ĥ) ĥ - s`. Energy is exactly preserved (ΔE = 0).
/// Works for any spin dimension (XY: sd=2, Heisenberg: sd=3).
#[inline]
fn reflect_spin(spin: &[f64], local_field: &[f64], sd: usize) -> SmallVec<[f64; 3]> {
    let h_norm: f64 = local_field.iter().map(|&x| x * x).sum::<f64>().sqrt();
    if h_norm < 1e-12 {
        return SmallVec::from_slice(spin);
    }
    let h_hat: SmallVec<[f64; 3]> = local_field.iter().map(|&x| x / h_norm).collect();
    let s_dot_h: f64 = spin.iter().zip(&h_hat).map(|(&s, &h)| s * h).sum();
    let mut new = SmallVec::from_elem(0.0, sd);
    for k in 0..sd {
        new[k] = 2.0 * s_dot_h * h_hat[k] - spin[k];
    }
    new
}

/// Microcanonical over-relaxation algorithm.
///
/// Visits every site once in random order and reflectes each spin across the
/// local field direction. Energy is exactly preserved — no acceptance step
/// needed. For XY and Heisenberg models. Mix with Metropolis/HeatBath sweeps
/// (typically 1 ergodic + 4-10 OR sweeps) to reduce critical slowing down.
#[derive(Debug, Clone, Default)]
pub struct MicrocanonicalCore;

impl MicrocanonicalCore {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + Proposable> Algorithm<H> for MicrocanonicalCore {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();

        // Random visit order
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for &site in &order {
            // Compute local field from neighbors
            let mut h = vec![0.0; sd];
            for &nb in system.lattice.neighbors(site) {
                let base = nb * sd;
                for (k, hk) in h.iter_mut().enumerate() {
                    *hk += system.spins[base + k];
                }
            }

            let old = system.spin_at(site, sd).to_vec();
            let reflected = reflect_spin(&old, &h, sd);
            let reflected_norm: f64 =
                reflected.iter().map(|&x| x * x).sum::<f64>().sqrt();

            // Normalize and write back
            let inv_norm = 1.0 / reflected_norm.max(1e-15);
            let spin = system.spin_at_mut(site, sd);
            for (k, s) in spin.iter_mut().enumerate() {
                *s = reflected.get(k).copied().unwrap_or(0.0) * inv_norm;
            }
        }
        // Energy is exactly preserved — no update needed
    }

    fn name(&self) -> &'static str {
        "Microcanonical"
    }
}

// ── Continuous Heat-Bath ──────────────────────────────────────

/// Continuous-spin heat-bath algorithm.
///
/// For XY and Heisenberg models. Visits every site once in random order and
/// samples a new spin from the equilibrium distribution P(s_i | neighbors)
/// using exact inverse-CDF (Heisenberg) or Best-Fisher rejection (XY).
#[derive(Debug, Clone, Default)]
pub struct ContinuousHeatBathCore;

impl ContinuousHeatBathCore {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + ContinuousHeatBathable> Algorithm<H> for ContinuousHeatBathCore {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let sd = model.spin_dim();
        let beta = system.beta;

        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for &site in &order {
            let old_spin = system.spin_at(site, sd).to_vec();
            let old_energy = model.local_energy(
                &system.spins,
                &system.lattice,
                site,
                beta,
                &old_spin,
            );

            // Collect neighbor spins as flat slice
            let nbs: Vec<f64> = system
                .lattice
                .neighbors(site)
                .iter()
                .flat_map(|&nb| {
                    let base = nb * sd;
                    system.spins[base..base + sd].to_vec()
                })
                .collect();

            let new_spin = model.heat_bath_sample(&nbs, beta, rng);
            let new_energy = model.local_energy(
                &system.spins,
                &system.lattice,
                site,
                beta,
                &new_spin,
            );

            system.energy += new_energy - old_energy;
            system.spin_at_mut(site, sd).copy_from_slice(&new_spin);
        }
    }

    fn name(&self) -> &'static str {
        "ContinuousHeatBath"
    }
}

// ── Heat-Bath (Glauber Dynamics) ─────────────────────────────

/// Heat-bath (Glauber dynamics) algorithm.
///
/// Visits every site once in random order and directly samples a new spin
/// from the equilibrium distribution P(s_i) ∝ exp(-β E(s_i | neighbors)),
/// with no Metropolis rejection step. Exact for scalar discrete models
/// (Ising, Potts).
#[derive(Debug, Clone, Default)]
pub struct HeatBathCore;

impl HeatBathCore {
    pub fn new() -> Self {
        Self
    }
}

impl<H: Hamiltonian + HeatBathable> Algorithm<H> for HeatBathCore {
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        let n = system.n_sites();
        let beta = system.beta;

        // Random visit order
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for &site in &order {
            // Collect neighbor spins
            let nbs: Vec<f64> = system
                .lattice
                .neighbors(site)
                .iter()
                .map(|&nb| system.spins[nb])
                .collect();

            let old_energy = model.local_energy(
                &system.spins,
                &system.lattice,
                site,
                beta,
                system.spin_at(site, 1),
            );

            let weights = model.boltzmann_weights(&nbs, beta);
            let new_val = model.sample_spin(&weights, rng);
            let new_spin_arr = [new_val];
            let new_energy =
                model.local_energy(&system.spins, &system.lattice, site, beta, &new_spin_arr);

            system.energy += new_energy - old_energy;
            system.spins[site] = new_val;
        }
    }

    fn name(&self) -> &'static str {
        "HeatBath"
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::models::{HeisenbergModel, IsingModel, PottsModel, XYModel};
    use rand::SeedableRng;

    fn make_rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    /// Check that all XY/Heisenberg spins are approximately parallel.
    fn all_aligned(spins: &[f64], sd: usize) -> bool {
        let ref_spin = &spins[..sd];
        for chunk in spins.chunks(sd) {
            let d: f64 = dot(ref_spin, chunk);
            if d < 0.99 {
                return false;
            }
        }
        true
    }

    #[test]
    fn test_metropolis_single_site() {
        let lattice = build_chain(1, false);
        let mut system = System::new(lattice, 1, 1.0, 1.0);
        system.energy = 0.0;

        let model = IsingModel::new(1.0);
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
        let model = IsingModel::new(1.0);
        let mut system = System::new(lattice.clone(), 1, 1.0, 10.0); // beta = 10, very cold
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

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
        let model = IsingModel::new(1.0);
        let mut system = System::new(lattice.clone(), 1, 1.0, 5.0); // beta=5, cold
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

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
        let model = IsingModel::new(1.0);
        let mut system = System::new(lattice.clone(), 1, 1.0, 5.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

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

    // ── Embedding tests (XY / Heisenberg) ──────────────────

    #[test]
    fn test_wolff_xy_single_site() {
        let lattice = build_chain(1, false);
        let mut system = System::new(lattice, 2, 1.0, 1.0);
        system.energy = 0.0;

        let model = XYModel::new(1.0);
        let mut algo = WolffCore::new();
        let mut rng = make_rng();

        for _ in 0..50 {
            algo.sweep(&mut system, &model, &mut rng);
        }
        assert!((system.energy - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_wolff_xy_cools_to_ground() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(1.0);
        let mut system = System::new(lattice.clone(), 2, 1.0, 10.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = WolffCore::new();
        let mut rng = make_rng();

        for _ in 0..800 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // At high beta, XY spins should align: 4 bonds × -J = -4
        assert!(system.energy < -3.5);
        assert!(all_aligned(&system.spins, 2));
    }

    #[test]
    fn test_wolff_heisenberg_cools_to_ground() {
        let lattice = build_chain(4, true);
        let model = HeisenbergModel::new(1.0);
        let mut system = System::new(lattice.clone(), 3, 1.0, 10.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = WolffCore::new();
        let mut rng = make_rng();

        for _ in 0..800 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        assert!(system.energy < -3.5);
        assert!(all_aligned(&system.spins, 3));
    }

    #[test]
    fn test_sw_xy_cools_to_ground() {
        let lattice = build_chain(8, true);
        let model = XYModel::new(1.0);
        let mut system = System::new(lattice.clone(), 2, 1.0, 5.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = SWCore::new();
        let mut rng = make_rng();

        for _ in 0..500 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        assert!(system.energy < -7.0);
        assert!(all_aligned(&system.spins, 2));
    }

    #[test]
    fn test_sw_heisenberg_cools_to_ground() {
        let lattice = build_chain(8, true);
        let model = HeisenbergModel::new(1.0);
        let mut system = System::new(lattice.clone(), 3, 1.0, 5.0);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = SWCore::new();
        let mut rng = make_rng();

        for _ in 0..500 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        assert!(system.energy < -7.0);
        assert!(all_aligned(&system.spins, 3));
    }

    #[test]
    fn test_heat_bath_ising_cools_to_ground() {
        let lattice = build_chain(8, true);
        let model = IsingModel::new(1.0);
        // Start random, beta=5 (cold)
        let mut system = System::new(lattice, 1, 1.0, 5.0);
        let mut rng = make_rng();
        // Randomize initial spins
        for i in 0..system.n_sites() {
            system.spins[i] = if rng.random::<bool>() { 1.0 } else { -1.0 };
        }
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = HeatBathCore::new();
        for _ in 0..200 {
            algo.sweep(&mut system, &model, &mut rng);
        }
        // At beta=5, should converge to ground state (all aligned)
        assert!(system.energy < -7.0, "energy = {}", system.energy);
    }

    // ── Microcanonical over-relaxation tests ───────────────────

    #[test]
    fn test_microcanonical_xy_energy_preserved() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(1.0);
        let mut system = System::new(lattice.clone(), 2, 0.0, 1.0);

        // Manual spin config: alternating (0°, 180°), not ground state
        // site 0: (1,0), site 1: (-1,0), site 2: (1,0), site 3: (-1,0)
        let spins_init: Vec<f64> = vec![1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
        system.spins.copy_from_slice(&spins_init);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let energy_before = system.energy;
        let mut algo = MicrocanonicalCore::new();
        let mut rng = make_rng();

        for _ in 0..10 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        let energy_after = system.energy;
        assert!(
            (energy_after - energy_before).abs() < 1e-10,
            "energy should be preserved: before={}, after={}",
            energy_before,
            energy_after
        );
    }

    #[test]
    fn test_microcanonical_heisenberg_energy_preserved() {
        let lattice = build_chain(4, true);
        let model = HeisenbergModel::new(1.0);
        let mut system = System::new(lattice.clone(), 3, 0.0, 1.0);

        // Manual 3D config
        #[rustfmt::skip]
        let spins_init: Vec<f64> = vec![
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            -1.0, 0.0, 0.0,
        ];
        system.spins.copy_from_slice(&spins_init);
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let energy_before = system.energy;
        let mut algo = MicrocanonicalCore::new();
        let mut rng = make_rng();

        for _ in 0..10 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        let energy_after = system.energy;
        assert!(
            (energy_after - energy_before).abs() < 1e-10,
            "energy should be preserved: before={}, after={}",
            energy_before,
            energy_after
        );
    }

    #[test]
    fn test_microcanonical_spins_change() {
        // Verify that microcanonical sweep actually changes spins
        // (i.e., it's not a no-op)
        let lattice = build_chain(8, true);
        let model = XYModel::new(1.0);
        let mut system = System::new(lattice.clone(), 2, 0.0, 1.0);

        // Random initial spins
        let mut rng = make_rng();
        for i in 0..system.n_sites() {
            let angle: f64 = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
            system.spins[2 * i] = angle.cos();
            system.spins[2 * i + 1] = angle.sin();
        }
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let spins_before = system.spins.clone();
        let energy_before = system.energy;
        let mut algo = MicrocanonicalCore::new();

        for _ in 0..5 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // Spins should have changed from the initial random config
        let changed = spins_before
            .iter()
            .zip(&system.spins)
            .any(|(a, b)| (a - b).abs() > 1e-10);
        assert!(changed, "spins should change during microcanonical sweep");

        // Energy must still be preserved
        assert!(
            (system.energy - energy_before).abs() < 1e-10,
            "energy preserved even with spin changes"
        );
    }

    // ── Continuous heat-bath tests ──────────────────────────

    #[test]
    fn test_continuous_heat_bath_xy_cools() {
        let lattice = build_chain(8, true);
        let model = XYModel::new(1.0);
        let mut system = System::new(lattice.clone(), 2, 0.0, 5.0); // beta=5 cold
        let mut rng = make_rng();
        // Random initial spins on S¹
        for i in 0..system.n_sites() {
            let angle: f64 = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
            system.spins[2 * i] = angle.cos();
            system.spins[2 * i + 1] = angle.sin();
        }
        system.energy =
            model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let energy_before = system.energy;
        let mut algo = ContinuousHeatBathCore::new();

        for _ in 0..200 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // At beta=5, should order: e/site < -0.7
        let e_per_site = system.energy / system.n_sites() as f64;
        assert!(
            e_per_site < -0.7,
            "XY heat-bath should cool: e/site = {:.4}",
            e_per_site
        );
        assert!(
            system.energy < energy_before,
            "energy should decrease from random state"
        );
    }

    #[test]
    fn test_continuous_heat_bath_heisenberg_cools() {
        let lattice = build_chain(8, true);
        let model = HeisenbergModel::new(1.0);
        let mut system = System::new(lattice.clone(), 3, 0.0, 5.0);
        let mut rng = make_rng();
        // Random initial spins on S²
        for i in 0..system.n_sites() {
            let z: f64 = rng.random::<f64>() * 2.0 - 1.0;
            let sin_theta = (1.0 - z * z).sqrt();
            let phi: f64 = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
            system.spins[3 * i] = sin_theta * phi.cos();
            system.spins[3 * i + 1] = sin_theta * phi.sin();
            system.spins[3 * i + 2] = z;
        }
        system.energy =
            model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let energy_before = system.energy;
        let mut algo = ContinuousHeatBathCore::new();

        for _ in 0..200 {
            algo.sweep(&mut system, &model, &mut rng);
        }

        // At beta=5, should order: 8 bonds × -J = -8, e/site < -0.7
        let e_per_site = system.energy / system.n_sites() as f64;
        assert!(
            e_per_site < -0.7,
            "Heisenberg heat-bath should cool: e/site = {:.4}",
            e_per_site
        );
        assert!(
            system.energy < energy_before,
            "energy should decrease from random state"
        );
    }

    #[test]
    fn test_heat_bath_potts_cools_to_ground() {
        let lattice = build_chain(8, true);
        let model = PottsModel::new(1.0, 3);
        let mut system = System::new(lattice, 1, 0.0, 5.0);
        let mut rng = make_rng();
        // Randomize initial spins
        for i in 0..system.n_sites() {
            system.spins[i] = rng.random_range(0..3) as f64;
        }
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, system.beta);

        let mut algo = HeatBathCore::new();
        for _ in 0..200 {
            algo.sweep(&mut system, &model, &mut rng);
        }
        // At beta=5, should converge (all same state, energy = -8)
        assert!(system.energy < -6.0, "energy = {}", system.energy);
    }
}
