//! Reusable classical Monte Carlo update kernels.

use crate::hamiltonian::{
    ClusterAuxiliary, ClusterModel, ContinuousHeatBathable, Hamiltonian, HeatBathable,
    LocalFieldModel, Proposable, Spin,
};
use crate::proposal::{ProposalStrategy, StandardStrategy};
use crate::system::System;
use rand::{Rng, RngExt};

/// Carlo.rs lifecycle phase visible to adaptive update kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimulationPhase {
    Thermalization,
    Measurement,
}

impl SimulationPhase {
    #[inline]
    pub const fn allows_adaptation(self) -> bool {
        matches!(self, Self::Thermalization)
    }
}

/// One update policy.  Carlo.rs owns scheduling; CMC.rs owns state transitions.
pub trait Algorithm<H: Hamiltonian>: Send {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    );

    /// Direct/manual sweeps default to the frozen measurement kernel.
    fn sweep(&mut self, system: &mut System, model: &H, rng: &mut impl Rng) {
        self.sweep_with_phase(system, model, rng, SimulationPhase::Measurement);
    }

    fn name(&self) -> &'static str {
        "Unknown"
    }
}

fn prepare_random_order(order: &mut Vec<usize>, n_sites: usize, rng: &mut impl Rng) {
    if order.len() != n_sites {
        order.clear();
        order.extend(0..n_sites);
    }
    for index in (1..n_sites).rev() {
        let swap_with = rng.random_range(0..=index);
        order.swap(index, swap_with);
    }
}

fn checked_probability(value: f64, algorithm: &str) -> f64 {
    assert!(
        value.is_finite() && (0.0..=1.0).contains(&value),
        "{algorithm} model returned invalid bond probability {value}"
    );
    value
}

// ── Metropolis-Hastings ─────────────────────────────────────

#[derive(Debug, Clone)]
pub struct MetropolisCore<S = StandardStrategy> {
    pub strategy: S,
    order: Vec<usize>,
    energy_check_interval: u64,
    sweeps: u64,
}

impl MetropolisCore<StandardStrategy> {
    pub fn new() -> Self {
        Self::with_strategy(StandardStrategy::new())
    }
}

impl<S: Default> Default for MetropolisCore<S> {
    fn default() -> Self {
        Self::with_strategy(S::default())
    }
}

impl<S> MetropolisCore<S> {
    pub fn with_strategy(strategy: S) -> Self {
        Self {
            strategy,
            order: Vec::new(),
            energy_check_interval: 0,
            sweeps: 0,
        }
    }

    /// Periodically replace the cached energy with an exact recomputation.
    /// Zero (default) disables the check.
    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }
}

impl<H, S> Algorithm<H> for MetropolisCore<S>
where
    H: Hamiltonian + Proposable,
    S: ProposalStrategy<H>,
{
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        prepare_random_order(&mut self.order, n_sites, rng);

        for &site in &self.order {
            let proposal = self.strategy.propose(model, system, site, rng);
            assert_eq!(
                proposal.spin.len(),
                spin_dim,
                "proposal dimension does not match the model"
            );
            let delta_energy =
                model.delta_energy(&system.spins, &system.lattice, site, &proposal.spin);
            assert!(
                delta_energy.is_finite(),
                "model returned non-finite delta energy"
            );
            assert!(
                proposal.log_reverse_over_forward.is_finite(),
                "proposal returned a non-finite Hastings correction"
            );
            let log_acceptance = -system.beta * delta_energy + proposal.log_reverse_over_forward;
            let accepted = log_acceptance >= 0.0
                || rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance;

            self.strategy.record_result(accepted);
            if accepted {
                system
                    .spin_at_mut(site, spin_dim)
                    .copy_from_slice(&proposal.spin);
                system.energy += delta_energy;
            }
        }

        self.strategy.finish_sweep(phase.allows_adaptation());
        self.sweeps = self.sweeps.wrapping_add(1);
        if self.energy_check_interval > 0 && self.sweeps.is_multiple_of(self.energy_check_interval)
        {
            system.recompute_energy(model);
        }
    }

    fn name(&self) -> &'static str {
        "Metropolis-Hastings"
    }
}

// ── Wolff ───────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct WolffCore {
    membership: Vec<bool>,
    stack: Vec<usize>,
    members: Vec<usize>,
}

impl WolffCore {
    pub const fn new() -> Self {
        Self {
            membership: Vec::new(),
            stack: Vec::new(),
            members: Vec::new(),
        }
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
        if self.membership.len() != n_sites {
            self.membership.resize(n_sites, false);
        }
        self.membership.fill(false);
        self.stack.clear();
        self.members.clear();

        let seed = rng.random_range(0..n_sites);
        let seed_spin = system.spin_at(seed, spin_dim).to_vec();
        let auxiliary = model.wolff_auxiliary(&seed_spin, rng);
        self.membership[seed] = true;
        self.stack.push(seed);

        while let Some(site) = self.stack.pop() {
            self.members.push(site);
            let left_base = site * spin_dim;
            let left = &system.spins[left_base..left_base + spin_dim];
            for (neighbor, edge_id) in system.lattice.incidences(site) {
                if self.membership[neighbor] {
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
                    self.membership[neighbor] = true;
                    self.stack.push(neighbor);
                }
            }
        }

        for &site in &self.members {
            let transformed =
                model.transform_cluster_spin(system.spin_at(site, spin_dim), &auxiliary);
            assert_eq!(
                transformed.len(),
                spin_dim,
                "cluster transform dimension mismatch"
            );
            system
                .spin_at_mut(site, spin_dim)
                .copy_from_slice(&transformed);
        }
        // Cluster moves touch many sites.  Exact recomputation prevents all
        // order-dependent local-energy accounting and accumulated drift.
        system.recompute_energy(model);
    }

    fn name(&self) -> &'static str {
        "Wolff"
    }
}

// ── Swendsen-Wang ───────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct SWCore {
    parent: Vec<usize>,
    rank: Vec<u8>,
    root_auxiliary: Vec<Option<ClusterAuxiliary>>,
}

impl SWCore {
    pub const fn new() -> Self {
        Self {
            parent: Vec::new(),
            rank: Vec::new(),
            root_auxiliary: Vec::new(),
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
                let representative = system.spin_at(site, spin_dim).to_vec();
                self.root_auxiliary[root] =
                    Some(model.sw_cluster_auxiliary(&representative, &bond_auxiliary, rng));
            }
        }

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
            system
                .spin_at_mut(site, spin_dim)
                .copy_from_slice(&transformed);
        }
        system.recompute_energy(model);
    }

    fn name(&self) -> &'static str {
        "Swendsen-Wang"
    }
}

// ── Exact microcanonical over-relaxation ────────────────────

#[derive(Debug, Clone, Default)]
pub struct MicrocanonicalCore {
    order: Vec<usize>,
    field: Vec<f64>,
}

impl MicrocanonicalCore {
    pub const fn new() -> Self {
        Self {
            order: Vec::new(),
            field: Vec::new(),
        }
    }
}

impl<H: Hamiltonian + LocalFieldModel> Algorithm<H> for MicrocanonicalCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        prepare_random_order(&mut self.order, n_sites, rng);
        self.field.resize(spin_dim, 0.0);

        for &site in &self.order {
            model.local_field(&system.spins, &system.lattice, site, &mut self.field);
            let norm_squared = self.field.iter().map(|value| value * value).sum::<f64>();
            if norm_squared < 1e-28 {
                continue;
            }
            let old = system.spin_at(site, spin_dim).to_vec();
            let projection = old
                .iter()
                .zip(&self.field)
                .map(|(spin, field)| spin * field)
                .sum::<f64>()
                / norm_squared;
            let mut reflected = Spin::from_slice(&old);
            for component in 0..spin_dim {
                reflected[component] = 2.0 * projection * self.field[component] - old[component];
            }
            system
                .spin_at_mut(site, spin_dim)
                .copy_from_slice(&reflected);
        }
        system.recompute_energy(model);
    }

    fn name(&self) -> &'static str {
        "Microcanonical"
    }
}

// ── Heat bath ───────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct HeatBathCore {
    order: Vec<usize>,
}

impl HeatBathCore {
    pub const fn new() -> Self {
        Self { order: Vec::new() }
    }
}

impl<H: Hamiltonian + HeatBathable> Algorithm<H> for HeatBathCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        prepare_random_order(&mut self.order, n_sites, rng);
        for &site in &self.order {
            let proposed =
                model.heat_bath_sample_site(&system.spins, &system.lattice, site, system.beta, rng);
            assert_eq!(
                proposed.len(),
                spin_dim,
                "heat-bath sample dimension mismatch"
            );
            let delta = model.delta_energy(&system.spins, &system.lattice, site, &proposed);
            assert!(delta.is_finite(), "model returned non-finite delta energy");
            system
                .spin_at_mut(site, spin_dim)
                .copy_from_slice(&proposed);
            system.energy += delta;
        }
    }

    fn name(&self) -> &'static str {
        "HeatBath"
    }
}

#[derive(Debug, Clone, Default)]
pub struct ContinuousHeatBathCore {
    order: Vec<usize>,
}

impl ContinuousHeatBathCore {
    pub const fn new() -> Self {
        Self { order: Vec::new() }
    }
}

impl<H: Hamiltonian + ContinuousHeatBathable> Algorithm<H> for ContinuousHeatBathCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let n_sites = system.n_sites();
        let spin_dim = model.spin_dim();
        prepare_random_order(&mut self.order, n_sites, rng);
        for &site in &self.order {
            let proposed =
                model.heat_bath_sample_site(&system.spins, &system.lattice, site, system.beta, rng);
            assert_eq!(
                proposed.len(),
                spin_dim,
                "heat-bath sample dimension mismatch"
            );
            let delta = model.delta_energy(&system.spins, &system.lattice, site, &proposed);
            assert!(delta.is_finite(), "model returned non-finite delta energy");
            system
                .spin_at_mut(site, spin_dim)
                .copy_from_slice(&proposed);
            system.energy += delta;
        }
    }

    fn name(&self) -> &'static str {
        "ContinuousHeatBath"
    }
}

// ── Hybrid composition ──────────────────────────────────────

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::models::{IsingModel, PottsModel, XYModel};
    use rand::SeedableRng;

    fn rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn metropolis_cache_matches_exact_energy() {
        let model = IsingModel::new(1.0);
        let mut system = System::new(build_chain(8, true), 1, 1.0, 0.8);
        system.recompute_energy(&model);
        let mut algorithm = MetropolisCore::new();
        let mut random = rng();
        for _ in 0..100 {
            algorithm.sweep(&mut system, &model, &mut random);
        }
        assert!(system.energy_error(&model).abs() < 1e-10);
    }

    #[test]
    fn wolff_on_energy_is_exact_after_batch_move() {
        let model = XYModel::new(1.0);
        let mut system = System::new(build_chain(8, true), 2, 0.0, 1.0);
        for spin in system.spins.chunks_exact_mut(2) {
            spin[0] = 1.0;
        }
        system.recompute_energy(&model);
        let mut algorithm = WolffCore::new();
        let mut random = rng();
        for _ in 0..20 {
            algorithm.sweep(&mut system, &model, &mut random);
            assert!(system.energy_error(&model).abs() < 1e-10);
        }
    }

    #[test]
    fn sw_potts_assigns_valid_independent_states() {
        let model = PottsModel::new(1.0, 5);
        let mut system = System::new(build_chain(16, false), 1, 0.0, 0.0);
        system.recompute_energy(&model);
        let mut algorithm = SWCore::new();
        let mut random = rng();
        algorithm.sweep(&mut system, &model, &mut random);
        assert!(system.spins.iter().all(|spin| (0.0..5.0).contains(spin)));
        // beta=0 forms singleton clusters; independent assignments should
        // almost surely produce more than one state with this fixed seed.
        assert!(system.spins.windows(2).any(|pair| pair[0] != pair[1]));
    }

    #[test]
    fn microcanonical_preserves_on_energy() {
        let model = XYModel::new(1.0);
        let mut system = System::new(build_chain(4, true), 2, 0.0, 1.0);
        system
            .spins
            .copy_from_slice(&[1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0]);
        system.recompute_energy(&model);
        let before = system.energy;
        let mut algorithm = MicrocanonicalCore::new();
        algorithm.sweep(&mut system, &model, &mut rng());
        assert!((system.energy - before).abs() < 1e-10);
    }
}
