//! Classical q-state Potts model: H = -J Σ_{<ij>} δ(s_i, s_j)

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Classical q-state Potts model: H = -J Σ_{<ij>} δ(s_i, s_j)
///
/// Spins s_i ∈ {0, 1, ..., q-1} are encoded as f64 values 0.0, 1.0, ..., (q-1).0.
pub struct PottsModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    q: usize,
    spins: Vec<f64>,
}

impl PottsModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64, q: usize) -> Self {
        let n_sites = lattice.n_sites;
        let spins = vec![0.0; n_sites];
        PottsModel {
            lattice,
            beta,
            j,
            q,
            spins,
        }
    }

    /// Number of Potts states.
    pub fn q(&self) -> usize {
        self.q
    }
}

impl LatticeMC for PottsModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}

impl ModelMC for PottsModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64) {
        let old = self.spins[site];
        let old_int = old as usize;

        // Pick a random new state from [0, q-1] that differs from old
        let new_int = if self.q == 1 {
            old_int // Only one state possible, must stay the same
        } else {
            let pick = rng.random_range(0..self.q - 1);
            if pick < old_int {
                pick
            } else {
                pick + 1
            }
        };

        (old, new_int as f64)
    }

    fn local_energy_change(&self, site: usize, old: f64, new: f64) -> f64 {
        let mut delta_e = 0.0;

        // Outgoing bonds: bonds where `site` is the source
        for neighbor in &self.lattice.sites[site] {
            let neighbor_spin = self.spins[neighbor.target];
            // If neighbor spin equals old value, we lose a favorable bond: +J
            if (neighbor_spin - old).abs() < 1e-10 {
                delta_e += self.j;
            }
            // If neighbor spin equals new value, we gain a favorable bond: -J
            if (neighbor_spin - new).abs() < 1e-10 {
                delta_e -= self.j;
            }
        }

        // Incoming bonds: bonds where `site` is the target
        for (source_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            if source_idx == site {
                continue;
            }
            for neighbor in neighbors {
                if neighbor.target == site {
                    let source_spin = self.spins[source_idx];
                    if (source_spin - old).abs() < 1e-10 {
                        delta_e += self.j;
                    }
                    if (source_spin - new).abs() < 1e-10 {
                        delta_e -= self.j;
                    }
                }
            }
        }

        // Bidirectional bonds count each physical bond twice; divide by 2
        delta_e / 2.0
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                if (self.spins[site_idx] - self.spins[neighbor.target]).abs() < 1e-10 {
                    energy -= self.j;
                }
            }
        }
        // Directed bonds count each physical bond twice; divide by 2
        energy / 2.0
    }

    fn spins(&self) -> &[f64] {
        &self.spins
    }

    fn spins_mut(&mut self) -> &mut [f64] {
        &mut self.spins
    }

    fn magnetization(&self) -> f64 {
        // For q=1, only one state exists → always fully ordered
        if self.q <= 1 {
            return 1.0;
        }
        let n = self.spins.len();
        let mut counts = vec![0usize; self.q];
        for &s in self.spins.iter() {
            let idx = s as usize;
            if idx < self.q {
                counts[idx] += 1;
            }
        }
        let max_count = counts.into_iter().max().unwrap_or(0);
        (self.q as f64 * max_count as f64 - n as f64) / (n as f64 * (self.q as f64 - 1.0))
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        rng.random_range(0..self.q) as f64
    }

    fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64 {
        let current = spin as usize;
        if self.q <= 1 {
            return spin;
        }
        let pick = rng.random_range(0..self.q - 1);
        if pick < current {
            pick as f64
        } else {
            (pick + 1) as f64
        }
    }

    fn fk_bond_probability(&self) -> f64 {
        // Potts model: H = -J Σ δ(s_i, s_j) → p_FK = 1 - exp(-βJ)
        // (No factor of 2, unlike Ising H = -J Σ s_i s_j → p_FK = 1 - exp(-2βJ))
        1.0 - (-self.coupling() * self.beta()).exp()
    }
}

// MonteCarlo implementation — placeholder, actual simulation done by *Core engines
impl MonteCarlo for PottsModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — engines handle the actual sweep
    }
}

impl FromParams for PottsModel {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let n_sites = params
            .get::<usize>("L")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "L".into(),
                reason: "System size L is required".into(),
            })?;

        let beta = params
            .get::<f64>("beta")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "beta".into(),
                reason: "Inverse temperature beta is required".into(),
            })?;

        let j = params.get::<f64>("J").unwrap_or(1.0);
        let q = params.get::<usize>("q").unwrap_or(3);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);

        let lattice = build_chain(n_sites, pbc);

        Ok(PottsModel::new(lattice, beta, j, q))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use rand_xoshiro::rand_core::SeedableRng;

    #[test]
    fn test_potts_ground_state() {
        // All spins same (0.0) → each undirected bond contributes -J
        // PBC chain of 4 sites has 4 bonds
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 3);
        // Expected: -J * n_bonds = -1.0 * 4 = -4.0
        assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_potts_local_energy_change() {
        // All spins 0.0, flip one spin to 1.0
        // PBC chain of 4: site 0 has 2 neighbors (sites 1 and 3), both spin 0.0
        // Losing 2 favorable bonds: delta_e = +2*J (outgoing) +2*J (incoming) / 2 = +2*J
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 3);

        let de = model.local_energy_change(0, 0.0, 1.0);
        // Each physical bond to neighbor contributes +J when we lose a match
        // Site 0 has 2 physical bonds (left and right), both currently matched
        // Outgoing: 2 neighbors both match old → +2J
        // Incoming: 2 bonds targeting site 0, both sources match old → +2J
        // Total: (2+2)/2 * J = 2.0
        assert!((de - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_from_params() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 2.0f64);
        params.set("J", 1.0f64);
        params.set("q", 4usize);

        let model = PottsModel::from_params(&params, &mut rng).unwrap();
        assert_eq!(model.n_sites(), 8);
        assert!((model.beta() - 2.0).abs() < 1e-10);
        assert_eq!(model.q(), 4);
    }

    #[test]
    fn test_potts_propose_flip() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 4);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);

        // Run many proposals, verify new spin always differs from old
        for _ in 0..100 {
            let (old, new) = model.propose_flip(0, &mut rng);
            assert!(
                (old - new).abs() > 1e-10,
                "proposed new spin must differ from old"
            );
            assert!(new >= 0.0 && new < 4.0);
            assert!(new == new.floor()); // new is an integer
        }
    }

    #[test]
    fn test_potts_magnetization_all_same() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 3);
        // All spins 0 → max_count = 4 → (3*4 - 4) / (4*2) = 8/8 = 1.0
        assert!((model.magnetization() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_uniform() {
        // q=4, 4 sites, each state once → M = 0
        let lattice = build_chain(4, true);
        let mut model = PottsModel::new(lattice, 1.0, 1.0, 4);
        for (i, s) in model.spins_mut().iter_mut().enumerate() {
            *s = i as f64;
        }
        assert!(model.magnetization() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_q1() {
        // q=1 → always fully ordered → M = 1
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 1);
        assert!((model.magnetization() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_random_cluster_spin() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 3);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
        for _ in 0..20 {
            let spin = model.random_cluster_spin(&mut rng);
            assert!(spin >= 0.0 && spin < 3.0);
            assert!(spin == spin.floor()); // Must be integer
        }
    }

    #[test]
    fn test_potts_opposite_spin() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 1.0, 1.0, 4);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
        for _ in 0..20 {
            let opp = model.opposite_spin(0.0, &mut rng);
            assert!(opp != 0.0, "opposite of 0 must not be 0");
            assert!(opp >= 1.0 && opp <= 3.0);
        }
    }
}
