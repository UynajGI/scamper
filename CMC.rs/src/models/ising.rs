//! Classical Ising model: H = -J Σ S_i S_j

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Classical Ising model: H = -J Σ S_i S_j
pub struct IsingModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    spins: Vec<f64>,
}

impl IsingModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64) -> Self {
        let n_sites = lattice.n_sites;
        let spins = vec![1.0; n_sites];
        IsingModel {
            lattice,
            beta,
            j,
            spins,
        }
    }

    pub fn j(&self) -> f64 {
        self.j
    }
}

impl LatticeMC for IsingModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}

impl ModelMC for IsingModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, _site: usize, _rng: &mut impl Rng) -> (f64, f64) {
        let old = self.spins[_site];
        (old, -old)
    }

    fn local_energy_change(&self, site: usize, old: f64, _new: f64) -> f64 {
        let mut delta_e = 0.0;

        // Outgoing bonds: bonds where `site` is the source
        for neighbor in &self.lattice.sites[site] {
            let neighbor_spin = self.spins[neighbor.target];
            // H = -J Σ_{undirected} S_i S_j, ΔE for flipping S_i: 2*J*S_i_old*S_j
            delta_e += 2.0 * self.j * old * neighbor_spin;
        }

        // Incoming bonds: bonds where `site` is the target
        for (source_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            if source_idx == site {
                continue;
            }
            for neighbor in neighbors {
                if neighbor.target == site {
                    let source_spin = self.spins[source_idx];
                    delta_e += 2.0 * self.j * source_spin * old;
                }
            }
        }

        // Each physical bond contributes once. Since directed bonds count each
        // physical bond twice (once in each direction), we divide by 2.
        delta_e / 2.0
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                energy -= self.j * self.spins[site_idx] * self.spins[neighbor.target];
            }
        }
        // Directed bonds count each physical bond twice; divide by 2 for
        // the standard undirected-bond Hamiltonian H = -J Σ_{<ij>} S_i S_j
        energy / 2.0
    }

    fn spins(&self) -> &[f64] {
        &self.spins
    }

    fn spins_mut(&mut self) -> &mut [f64] {
        &mut self.spins
    }

    fn magnetization(&self) -> f64 {
        let sum: f64 = self.spins.iter().sum();
        sum.abs() / self.spins.len() as f64
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        if rng.random::<f64>() < 0.5 {
            1.0
        } else {
            -1.0
        }
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin
    }
}

// MonteCarlo implementation — placeholder, actual simulation done by *Core
impl MonteCarlo for IsingModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — engines handle the actual sweep
    }
}

impl FromParams for IsingModel {
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
        let pbc = params.get::<bool>("pbc").unwrap_or(true);

        let lattice = build_chain(n_sites, pbc);

        Ok(IsingModel::new(lattice, beta, j))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;

    #[test]
    fn test_ising_all_up_energy() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        // All spins up: each bond contributes -J * 1 * 1 = -1
        // 4 bonds in PBC chain of 4 sites
        assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_ising_local_energy_change_flip() {
        let lattice = build_chain(4, true);
        let mut model = IsingModel::new(lattice, 1.0, 1.0);

        // Bidirectional bonds: each site has 2 neighbors (left and right)
        // Site 0: neighbors are 1 (right) and 3 (left)
        // Flipping site 0: old=+1, new=-1
        // Outgoing ΔE = 2*J*1*1 + 2*J*1*1 = 4
        // Incoming ΔE = 2*J*1*1 + 2*J*1*1 = 4
        // Total (divided by 2) = 4.0
        let de = model.local_energy_change(0, 1.0, -1.0);
        assert!((de - 4.0).abs() < 1e-10);

        // Actually flip the spin
        model.spins_mut()[0] = -1.0;
        let new_energy = model.total_energy();
        // 8 directed bonds: 4 involve site 0 (contribute +1 each), 4 don't (contribute -1 each)
        // E = (4 - 4) / 2 = 0.0
        assert!((new_energy - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_from_params() {
        use rand_xoshiro::rand_core::SeedableRng;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 2.0f64);
        params.set("J", 1.0f64);

        let model = IsingModel::from_params(&params, &mut rng).unwrap();
        assert_eq!(model.n_sites(), 8);
        assert!((model.beta() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_magnetization_all_up() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(lattice, 1.0, 1.0);
        assert!((model.magnetization() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_magnetization_alternating() {
        // Alternating spins: [1, -1, 1, -1] → M = 0
        let lattice = build_chain(4, true);
        let mut model = IsingModel::new(lattice, 1.0, 1.0);
        for (i, s) in model.spins_mut().iter_mut().enumerate() {
            *s = if i % 2 == 0 { 1.0 } else { -1.0 };
        }
        assert!(model.magnetization() < 1e-10);
    }
}
