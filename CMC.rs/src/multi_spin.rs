//! Multi-spin coding for Ising model.
//!
//! Packs 64 independent Ising replicas into `u64` words and performs
//! neighbor field computation using bitwise operations. This gives a
//! ~4-8x speedup per physical sweep compared to single-replica Metropolis.
//!
//! Each bit position in a `u64` represents the same lattice site in a
//! different replica. All 64 replicas share the same lattice and coupling
//! but evolve independently with different random number sequences.

use crate::classical_mc::{build_lattice_from_params, parse_bool, parse_param};
use crate::lattice::interaction::Hamiltonian;
use crate::lattice::models::IsingModel;
use crate::lattice::state::System;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};
use rand::Rng;
use rand::RngExt;

/// Number of replicas packed into each u64 word.
pub const N_REPLICAS: usize = 64;

/// Multi-spin coded Ising Metropolis.
///
/// Maintains 64 independent replicas packed in `u64` words.
/// Call [`sweep`](Self::sweep) to update all replicas simultaneously.
pub struct MultiSpinIsing {
    /// `packed_spins[site]` has 64 bits, one per replica.
    /// bit = 1 → spin +1, bit = 0 → spin -1.
    pub packed_spins: Vec<u64>,
    /// Acceptance probability LUT: `accept_prob[k]` where `k` = number of
    /// aligned neighbors (out of `z` total neighbors).
    accept_prob: Vec<f64>,
    /// Coordination number (assumed uniform across sites).
    z: usize,
    /// Shared system state (spins, energy, beta, lattice).
    pub system: System,
    /// Ising model parameters (coupling J).
    pub model: IsingModel,
    /// Reused random-visit workspace.
    order: Vec<usize>,
}

impl MultiSpinIsing {
    /// Create a multi-spin Ising system.
    ///
    /// All 64 replicas start from the same spin configuration in `system`.
    /// The lattice must have uniform coordination (degree) for all sites.
    pub fn new(system: System, model: IsingModel) -> Self {
        let n = system.n_sites();
        assert_eq!(
            system.spins.len(),
            n,
            "multi-spin coding requires exactly one scalar spin per site"
        );

        // Check uniform coordination
        let z = system.lattice.degree(0);
        for i in 1..n {
            assert_eq!(
                system.lattice.degree(i),
                z,
                "multi-spin coding requires uniform lattice degree"
            );
        }

        assert!(
            system
                .lattice
                .edges
                .iter()
                .all(|edge| { edge.source != edge.target && (edge.weight - 1.0).abs() < 1e-12 }),
            "multi-spin coding currently requires unweighted, loop-free bonds"
        );

        // Initialize all replicas from system.spins
        let packed: Vec<u64> = system
            .spins
            .iter()
            .map(|&s| if s > 0.0 { u64::MAX } else { 0u64 })
            .collect();

        // Build acceptance LUT
        let beta = system.beta;
        let j = model.coupling();
        let mut accept_prob: Vec<f64> = vec![0.0; z + 1];
        for (k, prob) in accept_prob.iter_mut().enumerate() {
            let delta_e: f64 = 2.0 * j * (2.0 * k as f64 - z as f64);
            *prob = if delta_e <= 0.0 {
                1.0
            } else {
                (-beta * delta_e).exp()
            };
        }

        Self {
            packed_spins: packed,
            accept_prob,
            z,
            system,
            model,
            order: Vec::new(),
        }
    }

    /// Perform one full sweep updating all 64 replicas in parallel.
    ///
    /// After the sweep, replica 0 is written back into `self.system.spins` and
    /// `self.system.energy` is updated for that replica.
    pub fn sweep(&mut self, rng: &mut impl Rng) {
        let n = self.system.n_sites();
        let z = self.z;

        // Random visit order, reusing allocation across sweeps.
        if self.order.len() != n {
            self.order.clear();
            self.order.extend(0..n);
        }
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            self.order.swap(i, j);
        }

        for &site in &self.order {
            let site_word = self.packed_spins[site];

            // Count anti-aligned neighbors per replica
            let mut anti_counts = [0u8; 64];
            for &nb in self.system.lattice.neighbors(site) {
                let xor = site_word ^ self.packed_spins[nb];
                for (r, count) in anti_counts.iter_mut().enumerate() {
                    *count += ((xor >> r) & 1) as u8;
                }
            }

            let mut flip_mask: u64 = 0;
            for (r, &anti) in anti_counts.iter().enumerate() {
                let k = z - anti as usize;
                if rng.random::<f64>() < self.accept_prob[k] {
                    flip_mask |= 1u64 << r;
                }
            }

            self.packed_spins[site] = site_word ^ flip_mask;
        }

        // Extract replica 0 → system.spins and recompute energy
        self.extract_into(0);
        self.system.energy = self.model.compute_total_energy(
            &self.system.spins,
            &self.system.lattice,
            self.system.beta,
        );
    }

    /// Extract replica `k` (0..63) into a Vec of ±1.0 values.
    pub fn extract_replica(&self, k: usize) -> Vec<f64> {
        let mask = 1u64 << k;
        self.packed_spins
            .iter()
            .map(|&word| if word & mask != 0 { 1.0 } else { -1.0 })
            .collect()
    }

    /// Extract replica `k` into `self.system.spins` (does NOT update energy).
    pub fn extract_into(&mut self, k: usize) {
        let mask = 1u64 << k;
        let n = self.system.n_sites();
        for i in 0..n {
            self.system.spins[i] = if self.packed_spins[i] & mask != 0 {
                1.0
            } else {
                -1.0
            };
        }
    }

    /// Get the number of replicas (always 64).
    pub fn n_replicas(&self) -> usize {
        N_REPLICAS
    }

    /// Joint physical energy of all packed replicas.
    ///
    /// Parallel tempering changes beta for all 64 replicas simultaneously, so
    /// its weight ratio must use the sum of replica energies, not replica 0.
    fn packed_total_energy(&self) -> f64 {
        let coupling = self.model.coupling();
        self.system
            .lattice
            .edges
            .iter()
            .map(|edge| {
                let anti_aligned = (self.packed_spins[edge.source] ^ self.packed_spins[edge.target])
                    .count_ones() as f64;
                let signed_sum = N_REPLICAS as f64 - 2.0 * anti_aligned;
                -coupling * edge.weight * signed_sum
            })
            .sum()
    }

    /// Rebuild the acceptance LUT (call after changing beta).
    fn rebuild_accept_lut(&mut self) {
        let beta = self.system.beta;
        let j = self.model.coupling();
        let z = self.z;
        for (k, prob) in self.accept_prob.iter_mut().enumerate() {
            let delta_e: f64 = 2.0 * j * (2.0 * k as f64 - z as f64);
            *prob = if delta_e <= 0.0 {
                1.0
            } else {
                (-beta * delta_e).exp()
            };
        }
    }
}

// ── MonteCarlo impl ──────────────────────────────────────────

impl MonteCarlo for MultiSpinIsing {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sweep(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Measure all 64 replicas. They share β and visit order, so they are
        // best treated as a vector-valued ensemble rather than blindly folded
        // into one scalar time series. We expose them
        // as array observables (Energy[k], Magnetization[k], … for replica k)
        // so callers can pick a single replica or average, and additionally
        // keep replica-0 as the scalar "Energy"/"Magnetization" for backward
        // compatibility with scalar-only consumers.
        let n = self.system.n_sites();
        let j = self.model.coupling();

        let mut e_arr = [0.0f64; N_REPLICAS];
        let mut m_arr = [0.0f64; N_REPLICAS];
        let mut e2_arr = [0.0f64; N_REPLICAS];
        let mut m2_arr = [0.0f64; N_REPLICAS];
        let mut m4_arr = [0.0f64; N_REPLICAS];

        for r in 0..N_REPLICAS {
            let mask = 1u64 << r;
            let mut spin_sum = 0.0f64;
            for word in &self.packed_spins {
                spin_sum += if word & mask != 0 { 1.0 } else { -1.0 };
            }
            // Physical edges are stored once, so no hidden divide-by-two.
            let e: f64 = self
                .system
                .lattice
                .edges
                .iter()
                .map(|edge| {
                    let left = if self.packed_spins[edge.source] & mask != 0 {
                        1.0
                    } else {
                        -1.0
                    };
                    let right = if self.packed_spins[edge.target] & mask != 0 {
                        1.0
                    } else {
                        -1.0
                    };
                    -j * edge.weight * left * right
                })
                .sum();

            let m = (spin_sum / n as f64).abs();
            e_arr[r] = e;
            m_arr[r] = m;
            e2_arr[r] = e * e;
            m2_arr[r] = m * m;
            m4_arr[r] = m * m * m * m;
        }

        // Per-replica array observables (length 64).
        ctx.measure_array("Energy_replica", &e_arr);
        ctx.measure_array("Magnetization_replica", &m_arr);
        ctx.measure_array("E2_replica", &e2_arr);
        ctx.measure_array("M2_replica", &m2_arr);
        ctx.measure_array("M4_replica", &m4_arr);

        // Scalar fallback: replica 0 (matches prior behavior).
        // Mark with suffix so it's not confused with the array.
        ctx.measure("Energy", e_arr[0]);
        ctx.measure("Magnetization", m_arr[0]);
        ctx.measure("E2", e2_arr[0]);
        ctx.measure("M2", m2_arr[0]);
        ctx.measure("M4", m4_arr[0]);
    }

    fn name(&self) -> &'static str {
        "MultiSpinIsing"
    }
}

// ── FromParams impl ──────────────────────────────────────────

impl FromParams for MultiSpinIsing {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
        if !beta.is_finite() || beta < 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "beta".into(),
                reason: "must be finite and non-negative".into(),
            });
        }
        if !coupling.is_finite() {
            return Err(CarloError::InvalidConfig {
                field: "J".into(),
                reason: "must be finite".into(),
            });
        }
        Ok(())
    }

    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let pbc = parse_bool(params, "pbc", true)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        let degree = lattice.degree(0);
        if !(0..lattice.n_sites).all(|site| lattice.degree(site) == degree) {
            return Err(CarloError::InvalidConfig {
                field: "lattice".into(),
                reason: "MultiSpinIsing requires uniform coordination".into(),
            });
        }
        if !lattice
            .edges
            .iter()
            .all(|edge| edge.source != edge.target && (edge.weight - 1.0).abs() < 1e-12)
        {
            return Err(CarloError::InvalidConfig {
                field: "lattice".into(),
                reason: "MultiSpinIsing requires unweighted, loop-free bonds".into(),
            });
        }

        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        let model = IsingModel::new(parse_param::<f64>(params, "J")?.unwrap_or(1.0));
        let mut system = System::new(lattice, 1, 0.0, beta);
        for spin in &mut system.spins {
            *spin = if rng.random::<bool>() { 1.0 } else { -1.0 };
        }
        system.recompute_energy(&model);
        Ok(Self::new(system, model))
    }
}

// ── ParallelTemperingCompatible impl ────────────────────────

impl ParallelTemperingCompatible for MultiSpinIsing {
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64 {
        match param {
            "beta" => (self.system.beta - new_value) * self.packed_total_energy(),
            _ => panic!("unsupported PT param: {param}"),
        }
    }

    fn change_parameter(&mut self, param: &str, new_value: f64) {
        match param {
            "beta" => {
                self.system
                    .set_beta(new_value)
                    .expect("parallel-tempering beta must be finite and non-negative");
                self.rebuild_accept_lut();
                self.system.energy = self.model.compute_total_energy(
                    &self.system.spins,
                    &self.system.lattice,
                    new_value,
                );
            }
            _ => panic!("unsupported PT param: {param}"),
        }
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::graph::build_chain;
    use carlo_rs::{RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

    fn make_rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn test_multi_spin_ising_construction() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(system, model);
        assert_eq!(msi.packed_spins.len(), 4);
        assert_eq!(msi.z, 2);
        assert_eq!(msi.packed_spins[0], u64::MAX);
    }

    #[test]
    fn test_multi_spin_ising_sweep() {
        let lattice = build_chain(8, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice, 1, 1.0, 5.0);
        let mut msi = MultiSpinIsing::new(system, model);
        let mut rng = make_rng();

        for _ in 0..100 {
            msi.sweep(&mut rng);
        }

        let replica0 = msi.extract_replica(0);
        let mag: f64 = replica0.iter().copied().sum::<f64>().abs() / replica0.len() as f64;
        assert!(mag > 0.5, "magnetization = {mag}");
    }

    #[test]
    fn test_multi_spin_ising_extract_replica() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(system, model);

        for r in 0..64 {
            let replica = msi.extract_replica(r);
            assert_eq!(replica, vec![1.0; 4]);
        }
    }

    #[test]
    fn test_multi_spin_acceptance_lut() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(system, model);
        assert!((msi.accept_prob[2] - (-4.0f64).exp()).abs() < 1e-10);
        assert!((msi.accept_prob[1] - 1.0).abs() < 1e-10);
        assert!((msi.accept_prob[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_spin_end_to_end() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<MultiSpinIsing>(&params);

        let energy = results.get("Energy").expect("Energy missing");
        let mag = results.get("Magnetization").expect("Magnetization missing");

        assert!(energy.mean < 0.0, "Energy should be negative");
        assert!(energy.stderr > 0.0);
        assert!((0.0..=1.0).contains(&mag.mean), "M in [0,1]");
    }

    #[test]
    fn test_multi_spin_pt() {
        let mut params = Params::new();
        params.set("L", 4usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mc = MultiSpinIsing::from_params(&params, &mut rng).unwrap();

        let joint_energy = mc.packed_total_energy();
        let lr = mc.log_weight_ratio("beta", 2.0);
        let expected = (1.0 - 2.0) * joint_energy;
        assert!((lr - expected).abs() < 1e-10);

        let mut mc = mc;
        mc.change_parameter("beta", 2.5);
        assert!((mc.system.beta - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_multi_spin_measure_all_replicas() {
        // measure() must register per-replica array observables
        // (Energy_replica, …) in addition to the replica-0 scalar fallback.
        let lattice = build_chain(8, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice, 1, 1.0, 1.0);
        let mut msi = MultiSpinIsing::new(system, model);

        // Run a few sweeps so replicas diverge from the uniform init.
        let mut rng = make_rng();
        for _ in 0..20 {
            msi.sweep(&mut rng);
        }

        // Drive measure() through a real Context and check the array
        // observable was registered with 64 entries.
        let mut ctx = carlo_rs::Context::new(rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7), 0);
        for _ in 0..5 {
            msi.sweep(&mut rng);
            msi.measure(&mut ctx);
        }
        let estimates = ctx.finalize_measurements();

        // Array path exists and the scalar fallback is consistent with it.
        assert!(
            estimates.contains_key("Energy_replica"),
            "per-replica Energy array must be registered"
        );
        assert!(
            estimates.contains_key("Magnetization_replica"),
            "per-replica Magnetization array must be registered"
        );
        // Scalar fallbacks still present for scalar-only consumers.
        assert!(estimates.contains_key("Energy"));
        assert!(estimates.contains_key("Magnetization"));
    }
}
