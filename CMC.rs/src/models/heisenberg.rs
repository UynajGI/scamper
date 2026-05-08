//! Classical Heisenberg model: H = -J Σ S⃗_i · S⃗_j
//!
//! S⃗_i = (sinθ cosφ, sinθ sinφ, cosθ) — unit vector on S².
//! Spins stored as 3 consecutive f64 values per site.

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::f64::consts::PI;

/// Classical Heisenberg model with naive Metropolis proposal.
pub struct HeisenbergModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    spins: Vec<f64>,          // 3 * n_sites: [x0, y0, z0, x1, y1, z1, ...]
    proposal_width: f64,      // Angular perturbation range [-δ, δ]
}

impl HeisenbergModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64, proposal_width: f64) -> Self {
        let n_sites = lattice.n_sites;
        let mut spins = vec![0.0; 3 * n_sites];
        for i in 0..n_sites {
            spins[3 * i + 2] = 1.0; // z = 1.0
        }
        HeisenbergModel {
            lattice, beta, j, spins, proposal_width,
        }
    }

    pub fn proposal_width(&self) -> f64 {
        self.proposal_width
    }
}

impl LatticeMC for HeisenbergModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}

impl ModelMC for HeisenbergModel {
    fn spin_dim(&self) -> usize {
        3
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, _rng: &mut impl Rng) -> (f64, f64) {
        // Fallback: return z-component (not used by MetropolisCore)
        (self.spins[3 * site + 2], 1.0)
    }

    fn local_energy_change(&self, _site: usize, _old: f64, _new: f64) -> f64 {
        unimplemented!("Use local_energy_change_spin for Heisenberg")
    }

    fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let sx = self.spins[3 * site];
        let sy = self.spins[3 * site + 1];
        let sz = self.spins[3 * site + 2];
        let old = vec![sx, sy, sz];

        // Random rotation axis (uniform on S²)
        let cos_theta = rng.random_range(-1.0..1.0);
        let sin_theta = f64::sqrt(1.0 - cos_theta * cos_theta);
        let phi = rng.random_range(0.0..2.0 * PI);
        let ax = sin_theta * phi.cos();
        let ay = sin_theta * phi.sin();
        let az = cos_theta;

        // Random rotation angle in [-proposal_width, proposal_width]
        let angle = rng.random_range(-self.proposal_width..self.proposal_width);
        let c = angle.cos();
        let s = angle.sin();
        let omc = 1.0 - c;

        // Rodrigues' rotation formula
        let dot = ax * sx + ay * sy + az * sz;
        let nx = c * sx + s * (ay * sz - az * sy) + omc * dot * ax;
        let ny = c * sy + s * (az * sx - ax * sz) + omc * dot * ay;
        let nz = c * sz + s * (ax * sy - ay * sx) + omc * dot * az;

        // Normalize
        let norm = (nx * nx + ny * ny + nz * nz).sqrt();
        let new = vec![nx / norm, ny / norm, nz / norm];
        (old, new)
    }

    fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
        // For bidirectional lattices, outgoing neighbors already cover all bonds
        // touching this site (each physical bond appears as both a→b and b→a).
        // The total_energy divides by 2, but the outgoing-only sum already gives
        // the correct local energy change (no division needed).
        let mut delta_e = 0.0;
        for neighbor in &self.lattice.sites[site] {
            let t = neighbor.target;
            let nx = self.spins[3 * t];
            let ny = self.spins[3 * t + 1];
            let nz = self.spins[3 * t + 2];
            delta_e -= self.j * (
                new[0] * nx + new[1] * ny + new[2] * nz
                - old[0] * nx - old[1] * ny - old[2] * nz
            );
        }
        delta_e
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let sx1 = self.spins[3 * site_idx];
                let sy1 = self.spins[3 * site_idx + 1];
                let sz1 = self.spins[3 * site_idx + 2];
                let sx2 = self.spins[3 * neighbor.target];
                let sy2 = self.spins[3 * neighbor.target + 1];
                let sz2 = self.spins[3 * neighbor.target + 2];
                energy -= self.j * (sx1 * sx2 + sy1 * sy2 + sz1 * sz2);
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
        let mut mz = 0.0;
        for i in 0..self.lattice.n_sites {
            mx += self.spins[3 * i];
            my += self.spins[3 * i + 1];
            mz += self.spins[3 * i + 2];
        }
        (mx * mx + my * my + mz * mz).sqrt() / n
    }

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        0.0 // Not applicable
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin
    }
}

// MonteCarlo implementation — placeholder, actual simulation done by *Core
impl MonteCarlo for HeisenbergModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — engines handle the actual sweep
    }
}

impl FromParams for HeisenbergModel {
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
        let proposal_width = params.get::<f64>("proposal_width").unwrap_or(PI / 8.0);
        let lattice = build_chain(n_sites, pbc);
        Ok(HeisenbergModel::new(lattice, beta, j, proposal_width))
    }
}
