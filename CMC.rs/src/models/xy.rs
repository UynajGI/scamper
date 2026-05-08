//! Classical XY model: H = -J Σ S⃗_i · S⃗_j
//!
//! S⃗_i = (cos θ_i, sin θ_i) — unit vector on S¹.
//! Spins stored as 2 consecutive f64 values per site: [x0, y0, x1, y1, ...]

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::f64::consts::PI;

/// Classical XY model: H = -J Σ S⃗_i · S⃗_j
///
/// Each spin is a unit 2D vector (cos θ, sin θ).
pub struct XYModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    spins: Vec<f64>,          // 2 * n_sites: [x0, y0, x1, y1, ...]
    proposal_width: f64,      // Angular perturbation range [-δ, δ]
}

impl XYModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64, proposal_width: f64) -> Self {
        let n_sites = lattice.n_sites;
        let mut spins = vec![0.0; 2 * n_sites];
        for i in 0..n_sites {
            spins[2 * i] = 1.0;     // x = cos(0) = 1.0
            spins[2 * i + 1] = 0.0; // y = sin(0) = 0.0
        }
        XYModel {
            lattice, beta, j, spins, proposal_width,
        }
    }

    pub fn proposal_width(&self) -> f64 {
        self.proposal_width
    }
}

impl LatticeMC for XYModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}

impl ModelMC for XYModel {
    fn spin_dim(&self) -> usize {
        2
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64) {
        let (old, new) = self.propose_flip_spin(site, rng);
        (old[0], new[0]) // Fallback: return x-component
    }

    fn local_energy_change(&self, _site: usize, _old: f64, _new: f64) -> f64 {
        unimplemented!("Use local_energy_change_spin for XYModel")
    }

    fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let x = self.spins[2 * site];
        let y = self.spins[2 * site + 1];

        // 2D rotation by random angle in [-proposal_width, proposal_width]
        let angle = rng.random_range(-self.proposal_width..self.proposal_width);
        let c = angle.cos();
        let s = angle.sin();

        let nx = c * x - s * y;
        let ny = s * x + c * y;

        // Normalize for numerical stability
        let norm = (nx * nx + ny * ny).sqrt();
        (vec![x, y], vec![nx / norm, ny / norm])
    }

    fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
        let mut delta_e = 0.0;
        for neighbor in &self.lattice.sites[site] {
            let t = neighbor.target;
            let nx = self.spins[2 * t];
            let ny = self.spins[2 * t + 1];
            delta_e -= self.j * (
                new[0] * nx + new[1] * ny
                - old[0] * nx - old[1] * ny
            );
        }
        delta_e
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let sx1 = self.spins[2 * site_idx];
                let sy1 = self.spins[2 * site_idx + 1];
                let sx2 = self.spins[2 * neighbor.target];
                let sy2 = self.spins[2 * neighbor.target + 1];
                energy -= self.j * (sx1 * sx2 + sy1 * sy2);
            }
        }
        energy / 2.0
    }

    fn spins(&self) -> &[f64] {
        &self.spins
    }

    fn spins_mut(&mut self) -> &mut [f64] {
        &mut self.spins
    }

    fn magnetization(&self) -> f64 {
        let n = self.lattice.n_sites as f64;
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..self.lattice.n_sites {
            mx += self.spins[2 * i];
            my += self.spins[2 * i + 1];
        }
        (mx * mx + my * my).sqrt() / n
    }

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        0.0 // Not applicable for continuous spin models
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin // Negation: x → -x
    }
}

impl MonteCarlo for XYModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — engines handle the actual sweep
    }
}

impl FromParams for XYModel {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let n_sites = params.get::<usize>("L").ok_or_else(|| CarloError::InvalidConfig {
            field: "L".into(),
            reason: "System size L is required".into(),
        })?;
        let beta = params.get::<f64>("beta").ok_or_else(|| CarloError::InvalidConfig {
            field: "beta".into(),
            reason: "Inverse temperature beta is required".into(),
        })?;
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);
        let proposal_width = params.get::<f64>("proposal_width").unwrap_or(PI / 4.0);
        let lattice = build_chain(n_sites, pbc);
        Ok(XYModel::new(lattice, beta, j, proposal_width))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use rand_xoshiro::rand_core::SeedableRng;

    #[test]
    fn test_xy_ground_state() {
        // All spins (1,0) → dot product = 1, energy = -J * n_physical_bonds
        // PBC chain of 4: 4 physical bonds
        let lattice = build_chain(4, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_xy_spin_norm() {
        let lattice = build_chain(8, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        for i in 0..8 {
            let sx = model.spins()[2 * i];
            let sy = model.spins()[2 * i + 1];
            let norm = (sx * sx + sy * sy).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "Spin {} has norm {}", i, norm);
        }
    }

    #[test]
    fn test_xy_total_energy() {
        // Set specific spin directions and verify total energy
        let lattice = build_chain(4, true);
        let mut model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        // Set angles: 0, PI/2, PI, 3*PI/2 → (1,0), (0,1), (-1,0), (0,-1)
        let angles = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0];
        for (i, &a) in angles.iter().enumerate() {
            model.spins_mut()[2 * i] = a.cos();
            model.spins_mut()[2 * i + 1] = a.sin();
        }

        // Bonds (undirected): (0,1), (1,2), (2,3), (3,0)
        // dot products: 1*0+0*1=0, 0*(-1)+1*0=0, (-1)*0+0*(-1)=0, 0*1+(-1)*0=0
        // energy = -1.0 * (0 + 0 + 0 + 0) = 0.0
        assert!((model.total_energy() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_propose_flip() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);

        // Run many proposals, verify proposed spin differs from old and has unit norm
        for _ in 0..100 {
            let (old, new) = model.propose_flip_spin(0, &mut rng);
            assert!(
                (old[0] - new[0]).abs() > 1e-10 || (old[1] - new[1]).abs() > 1e-10,
                "proposed new spin must differ from old"
            );
            let norm = (new[0] * new[0] + new[1] * new[1]).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "proposed spin must have unit norm, got {}",
                norm
            );
        }
    }

    #[test]
    fn test_xy_from_params() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 2.0f64);
        params.set("J", 1.0f64);

        let model = XYModel::from_params(&params, &mut rng).unwrap();
        assert_eq!(model.n_sites(), 8);
        assert!((model.beta() - 2.0).abs() < 1e-10);
        assert!((model.proposal_width() - PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_all_aligned() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        // All spins (1,0) → perfect alignment
        assert!((model.magnetization() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_uniform() {
        // Uniformly spread angles: 0, PI/2, PI, 3PI/2 → M = 0
        let lattice = build_chain(4, true);
        let mut model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        let angles = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0];
        for (i, &a) in angles.iter().enumerate() {
            model.spins_mut()[2 * i] = a.cos();
            model.spins_mut()[2 * i + 1] = a.sin();
        }
        assert!(model.magnetization() < 1e-10);
    }
}
