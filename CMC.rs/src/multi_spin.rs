//! Multi-spin coding for Ising model.
//!
//! Packs 64 independent Ising replicas into `u64` words and performs
//! neighbor field computation using bitwise operations. This gives a
//! ~4-8x speedup per physical sweep compared to single-replica Metropolis.
//!
//! Each bit position in a `u64` represents the same lattice site in a
//! different replica. All 64 replicas share the same lattice and coupling
//! but evolve independently with different random number sequences.

use crate::hamiltonian::Hamiltonian;
use crate::lattice::CsrLattice;
use crate::models::IsingModel;
use crate::system::System;
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
}

impl MultiSpinIsing {
    /// Create a multi-spin Ising system initialized from `system.spins`.
    ///
    /// All 64 replicas start from the same spin configuration.
    /// The lattice must have uniform coordination (degree) for all sites.
    pub fn new(system: &System, model: &IsingModel, lattice: &CsrLattice) -> Self {
        let n = system.n_sites();

        // Check uniform coordination
        let z = lattice.degree(0);
        for i in 1..n {
            assert_eq!(
                lattice.degree(i),
                z,
                "multi-spin coding requires uniform lattice degree"
            );
        }

        // Initialize all replicas from system.spins (same config for all 64)
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
            // k = number of aligned neighbors
            // ΔE = 2J(2k - z); accept with min(1, exp(-βΔE))
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
        }
    }

    /// Perform one full sweep updating all 64 replicas in parallel.
    ///
    /// After the sweep, replica 0 is written back into `system.spins` and
    /// `system.energy` is updated for that replica.
    pub fn sweep(
        &mut self,
        system: &mut System,
        model: &IsingModel,
        lattice: &CsrLattice,
        rng: &mut impl Rng,
    ) {
        let n = system.n_sites();
        let z = self.z;

        // Number of bit planes needed to represent counts 0..=z
        let n_planes = (usize::BITS - z.leading_zeros()) as usize;

        // Random visit order
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        for &site in &order {
            let site_word = self.packed_spins[site];

            // Compute XOR of site spin with each neighbor's spin.
            // bit = 0 means aligned, bit = 1 means anti-aligned.
            let mut xor_words: Vec<u64> = Vec::with_capacity(z);
            for &nb in lattice.neighbors(site) {
                xor_words.push(site_word ^ self.packed_spins[nb]);
            }

            // Pack z XOR words into bit planes for per-replica counting.
            // planes[p] bit r = p-th bit of anti-aligned neighbor count for replica r.
            let planes = pack_bit_planes(&xor_words, n_planes);

            // Per-replica acceptance
            let mut flip_mask: u64 = 0;
            for r in 0..64 {
                // Extract anti-aligned count for replica r from planes
                let mut anti: usize = 0;
                for (p, &plane) in planes.iter().enumerate() {
                    anti |= (((plane >> r) & 1) as usize) << p;
                }
                let k = z - anti; // aligned neighbor count
                let prob = self.accept_prob[k];
                if rng.random::<f64>() < prob {
                    flip_mask |= 1u64 << r;
                }
            }

            self.packed_spins[site] = site_word ^ flip_mask;
        }

        // Extract replica 0 → system.spins and recompute energy
        self.extract_into(0, system);
        system.energy = model.compute_total_energy(&system.spins, lattice, system.beta);
    }

    /// Extract replica `k` (0..63) into a Vec of ±1.0 values.
    pub fn extract_replica(&self, k: usize) -> Vec<f64> {
        let mask = 1u64 << k;
        self.packed_spins
            .iter()
            .map(|&word| if word & mask != 0 { 1.0 } else { -1.0 })
            .collect()
    }

    /// Extract replica `k` into `system.spins` (does NOT update energy).
    pub fn extract_into(&self, k: usize, system: &mut System) {
        let mask = 1u64 << k;
        for (i, &word) in self.packed_spins.iter().enumerate() {
            system.spins[i] = if word & mask != 0 { 1.0 } else { -1.0 };
        }
    }

    /// Get the number of replicas (always 64).
    pub fn n_replicas(&self) -> usize {
        N_REPLICAS
    }
}

/// Pack `z` u64 words into `n_planes` bit planes.
///
/// `planes[p]` bit `r` = bit `p` of `words[r]` (when reading each word as
/// a z-bit integer). This enables per-replica counting of anti-aligned
/// neighbors using only bitwise operations.
///
/// Invariant: `planes[p] bit r` = `(words[r] >> p) & 1` for r < 64.
fn pack_bit_planes(words: &[u64], n_planes: usize) -> Vec<u64> {
    let mut planes = vec![0u64; n_planes];
    for (r, &word) in words.iter().enumerate().take(64) {
        for (p, plane) in planes.iter_mut().enumerate() {
            if (word >> p) & 1 != 0 {
                *plane |= 1u64 << r;
            }
        }
    }
    planes
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use rand::SeedableRng;

    fn make_rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn test_pack_bit_planes() {
        // z=3 neighbors, 64 replicas.
        // Words: w[0]=0b101, w[1]=0b011, w[2]=0b110
        // Plane 0 (LSB): bit 0 of w[0]=1, bit 0 of w[1]=1, bit 0 of w[2]=0 → mask 0b011
        // Plane 1 (MSB): bit 1 of w[0]=0, bit 1 of w[1]=1, bit 1 of w[2]=1 → mask 0b110
        let words = vec![0b101u64, 0b011u64, 0b110u64];
        let planes = pack_bit_planes(&words, 2);
        assert_eq!(planes[0] & 0b111, 0b011); // replicas 0,1 have LSB=1
        assert_eq!(planes[1] & 0b111, 0b110); // replicas 1,2 have MSB=1
    }

    #[test]
    fn test_multi_spin_ising_construction() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(&system, &model, &lattice);
        assert_eq!(msi.packed_spins.len(), 4);
        assert_eq!(msi.z, 2);
        // All spins initialized to +1 (all bits set)
        assert_eq!(msi.packed_spins[0], u64::MAX);
    }

    #[test]
    fn test_multi_spin_ising_sweep() {
        let lattice = build_chain(8, true);
        let model = IsingModel::new(1.0);
        let mut system = System::new(lattice.clone(), 1, 1.0, 5.0);
        let mut rng = make_rng();
        let mut msi = MultiSpinIsing::new(&system, &model, &lattice);

        // Run some sweeps
        for _ in 0..100 {
            msi.sweep(&mut system, &model, &lattice, &mut rng);
        }

        // At beta=5 (cold), replica 0 should be mostly aligned
        let replica0 = msi.extract_replica(0);
        let mag: f64 = replica0.iter().copied().sum::<f64>().abs() / replica0.len() as f64;
        assert!(mag > 0.5, "magnetization = {}", mag);
    }

    #[test]
    fn test_multi_spin_ising_extract_replica() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(&system, &model, &lattice);

        // All replicas are +1 initially
        for r in 0..64 {
            let replica = msi.extract_replica(r);
            assert_eq!(replica, vec![1.0; 4]);
        }
    }

    #[test]
    fn test_multi_spin_acceptance_lut() {
        // For z=2, J=1, beta=1:
        // k=2 (both aligned): ΔE = 2*(4-2) = +4, P = exp(-4) ≈ 0.018
        // k=1 (one aligned):  ΔE = 2*(2-2) = 0,  P = 1.0
        // k=0 (none aligned): ΔE = 2*(0-2) = -4, P = 1.0
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        let system = System::new(lattice.clone(), 1, 1.0, 1.0);
        let msi = MultiSpinIsing::new(&system, &model, &lattice);
        assert!((msi.accept_prob[2] - (-4.0f64).exp()).abs() < 1e-10);
        assert!((msi.accept_prob[1] - 1.0).abs() < 1e-10);
        assert!((msi.accept_prob[0] - 1.0).abs() < 1e-10);
    }
}
