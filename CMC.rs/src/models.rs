//! Physics model implementations.
//!
//! Each model implements the appropriate subset of traits:
//! - Hamiltonian: Energy computation
//! - ClusterModel: Cluster algorithm support
//! - Proposable: Spin proposal
//! - Measurable: Magnetization

use crate::hamiltonian::{
    ClusterModel, ContinuousHeatBathable, Hamiltonian, HeatBathable, Measurable, Proposable,
};
use crate::lattice::CsrLattice;
use rand::Rng;
use rand::RngExt;
use smallvec::{smallvec, SmallVec};

// ── Ising Model ─────────────────────────────────────────────

/// Ising model: H = -J Σ_{⟨i,j⟩} σ_i σ_j, σ_i ∈ {±1}.
///
/// Parameters: `j` (coupling). Temperature (β) is stored in [`System`](crate::System).
#[derive(Debug, Clone)]
pub struct IsingModel {
    pub j: f64,
}

impl IsingModel {
    pub fn new(j: f64) -> Self {
        Self { j }
    }
}

impl Hamiltonian for IsingModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        let s = proposed[0];
        let mut e = 0.0;
        for &nb in lattice.neighbors(site) {
            e += -self.j * s * spins[nb];
        }
        e
    }
}

impl ClusterModel for IsingModel {
    fn fk_bond_probability(&self, beta: f64) -> f64 {
        1.0 - (-2.0 * self.j * beta).exp()
    }

    fn flip_in_place(&self, spin: &mut [f64], _rng: &mut impl Rng) {
        spin[0] = -spin[0];
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        if rng.random::<bool>() {
            1.0
        } else {
            -1.0
        }
    }
}

impl Proposable for IsingModel {
    fn propose(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        if rng.random::<bool>() {
            smallvec![1.0]
        } else {
            smallvec![-1.0]
        }
    }
}

impl Measurable for IsingModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let sum: f64 = spins.iter().sum();
        (sum / spins.len() as f64).abs()
    }
}

impl HeatBathable for IsingModel {
    fn n_states(&self) -> usize {
        2
    }

    fn boltzmann_weights(&self, neighbors: &[f64], beta: f64) -> Vec<f64> {
        let h: f64 = neighbors.iter().sum::<f64>() * self.j;
        // w[0] = P(+1) ∝ exp(βh), w[1] = P(-1) ∝ exp(-βh)
        vec![(beta * h).exp(), (-beta * h).exp()]
    }

    fn sample_spin(&self, weights: &[f64], rng: &mut impl Rng) -> f64 {
        let total = weights[0] + weights[1];
        if rng.random::<f64>() < weights[0] / total {
            1.0
        } else {
            -1.0
        }
    }
}

// ── Potts Model ──────────────────────────────────────────────

/// q-state Potts model: H = -J Σ δ(s_i, s_j), s_i ∈ {0, 1, ..., q-1}.
///
/// FK bond probability: `1 - exp(-βJ)` (no factor 2, unlike Ising).
#[derive(Debug, Clone)]
pub struct PottsModel {
    pub j: f64,
    pub q: usize,
}

impl PottsModel {
    pub fn new(j: f64, q: usize) -> Self {
        assert!(q >= 2, "Potts q must be >= 2");
        Self { j, q }
    }

    fn state_counts(&self, spins: &[f64]) -> Vec<usize> {
        let mut counts = vec![0usize; self.q];
        for &s in spins {
            let k = s as usize;
            if k < self.q {
                counts[k] += 1;
            }
        }
        counts
    }
}

impl Hamiltonian for PottsModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        let s = proposed[0] as usize;
        let mut e = 0.0;
        for &nb in lattice.neighbors(site) {
            let nb_state = spins[nb] as usize;
            if s == nb_state {
                e += -self.j;
            }
        }
        e
    }
}

impl ClusterModel for PottsModel {
    fn fk_bond_probability(&self, beta: f64) -> f64 {
        1.0 - (-self.j * beta).exp()
    }

    fn flip_in_place(&self, spin: &mut [f64], rng: &mut impl Rng) {
        let current = spin[0] as usize;
        let mut new = rng.random_range(0..self.q - 1);
        if new >= current {
            new += 1;
        }
        spin[0] = new as f64;
    }

    fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64 {
        let current = spin as usize;
        let mut new = rng.random_range(0..self.q - 1);
        if new >= current {
            new += 1;
        }
        new as f64
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        rng.random_range(0..self.q) as f64
    }
}

impl Proposable for PottsModel {
    fn propose(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        smallvec![rng.random_range(0..self.q) as f64]
    }
}

impl Measurable for PottsModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let n = spins.len();
        if n == 0 {
            return 0.0;
        }
        let counts = self.state_counts(spins);
        let max_n = counts.iter().max().copied().unwrap_or(0);
        (self.q as f64 * max_n as f64 - n as f64) / (n as f64 * (self.q - 1) as f64)
    }
}

impl HeatBathable for PottsModel {
    fn n_states(&self) -> usize {
        self.q
    }

    fn boltzmann_weights(&self, neighbors: &[f64], beta: f64) -> Vec<f64> {
        let mut counts = vec![0usize; self.q];
        for &s in neighbors {
            let k = s as usize;
            if k < self.q {
                counts[k] += 1;
            }
        }
        // w[k] = exp(βJ * n_k) where n_k = number of neighbors in state k
        counts
            .iter()
            .map(|&n| (beta * self.j * n as f64).exp())
            .collect()
    }

    fn sample_spin(&self, weights: &[f64], rng: &mut impl Rng) -> f64 {
        let total: f64 = weights.iter().sum();
        let mut u = rng.random::<f64>() * total;
        for (k, &w) in weights.iter().enumerate() {
            u -= w;
            if u <= 0.0 {
                return k as f64;
            }
        }
        (weights.len() - 1) as f64
    }
}

// ── XY Model ─────────────────────────────────────────────────

/// XY model: H = -J Σ cos(θ_i - θ_j) = -J Σ s_i · s_j, |s_i| = 1.
///
/// Each spin is a unit vector in 2D: `(cos θ, sin θ)`.
#[derive(Debug, Clone)]
pub struct XYModel {
    pub j: f64,
}

impl XYModel {
    pub fn new(j: f64) -> Self {
        Self { j }
    }
}

impl Hamiltonian for XYModel {
    fn spin_dim(&self) -> usize {
        2
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        let (sx, sy) = (proposed[0], proposed[1]);
        let mut e = 0.0;
        for &nb in lattice.neighbors(site) {
            let base = nb * 2;
            e += -self.j * (sx * spins[base] + sy * spins[base + 1]);
        }
        e
    }
}

impl ClusterModel for XYModel {
    fn fk_bond_probability(&self, beta: f64) -> f64 {
        1.0 - (-self.j * beta).exp()
    }

    fn reflect(&self, spin: &mut [f64], direction: &[f64]) {
        let proj = spin[0] * direction[0] + spin[1] * direction[1];
        spin[0] -= 2.0 * proj * direction[0];
        spin[1] -= 2.0 * proj * direction[1];
        self.normalize_spin(spin);
    }

    fn embedding_direction(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        let theta: f64 = rng.random_range(0.0..std::f64::consts::TAU);
        smallvec![theta.cos(), theta.sin()]
    }
}

impl Proposable for XYModel {
    fn propose(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        let theta: f64 = rng.random_range(0.0..std::f64::consts::TAU);
        smallvec![theta.cos(), theta.sin()]
    }

    fn normalize_spin(&self, spin: &mut [f64]) {
        let r = (spin[0] * spin[0] + spin[1] * spin[1]).sqrt();
        if r > 1e-12 {
            spin[0] /= r;
            spin[1] /= r;
        }
    }
}

impl Measurable for XYModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let (mut sx, mut sy) = (0.0, 0.0);
        for chunk in spins.chunks(2) {
            sx += chunk[0];
            sy += chunk[1];
        }
        let n = (spins.len() / 2) as f64;
        (sx * sx + sy * sy).sqrt() / n
    }
}

// ── Heisenberg Model ─────────────────────────────────────────

/// Heisenberg model: H = -J Σ s_i · s_j, |s_i| = 1.
///
/// Each spin is a unit vector in 3D. Uses Marsaglia's method for uniform sampling on S².
#[derive(Debug, Clone)]
pub struct HeisenbergModel {
    pub j: f64,
}

impl HeisenbergModel {
    pub fn new(j: f64) -> Self {
        Self { j }
    }
}

impl Hamiltonian for HeisenbergModel {
    fn spin_dim(&self) -> usize {
        3
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        let (sx, sy, sz) = (proposed[0], proposed[1], proposed[2]);
        let mut e = 0.0;
        for &nb in lattice.neighbors(site) {
            let base = nb * 3;
            e += -self.j * (sx * spins[base] + sy * spins[base + 1] + sz * spins[base + 2]);
        }
        e
    }
}

impl ClusterModel for HeisenbergModel {
    fn fk_bond_probability(&self, beta: f64) -> f64 {
        1.0 - (-self.j * beta).exp()
    }

    fn reflect(&self, spin: &mut [f64], direction: &[f64]) {
        let proj = spin[0] * direction[0] + spin[1] * direction[1] + spin[2] * direction[2];
        spin[0] -= 2.0 * proj * direction[0];
        spin[1] -= 2.0 * proj * direction[1];
        spin[2] -= 2.0 * proj * direction[2];
        self.normalize_spin(spin);
    }

    fn embedding_direction(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        let theta: f64 = rng.random_range(0.0..std::f64::consts::TAU);
        let phi: f64 = rng.random_range(0.0..std::f64::consts::PI);
        smallvec![phi.sin() * theta.cos(), phi.sin() * theta.sin(), phi.cos()]
    }
}

impl Proposable for HeisenbergModel {
    fn propose(&self, rng: &mut impl Rng) -> SmallVec<[f64; 3]> {
        // Marsaglia: sample (x,y) in unit disk, reject if outside
        let (x, y) = loop {
            let x: f64 = rng.random_range(-1.0..1.0);
            let y: f64 = rng.random_range(-1.0..1.0);
            if x * x + y * y <= 1.0 {
                break (x, y);
            }
        };
        let r = (x * x + y * y).sqrt();
        smallvec![
            2.0 * x * (1.0 - r * r).sqrt(),
            2.0 * y * (1.0 - r * r).sqrt(),
            1.0 - 2.0 * (x * x + y * y),
        ]
    }

    fn normalize_spin(&self, spin: &mut [f64]) {
        let r = (spin[0] * spin[0] + spin[1] * spin[1] + spin[2] * spin[2]).sqrt();
        if r > 1e-12 {
            spin[0] /= r;
            spin[1] /= r;
            spin[2] /= r;
        }
    }
}

impl Measurable for HeisenbergModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
        for chunk in spins.chunks(3) {
            sx += chunk[0];
            sy += chunk[1];
            sz += chunk[2];
        }
        let n = (spins.len() / 3) as f64;
        (sx * sx + sy * sy + sz * sz).sqrt() / n
    }
}

impl ContinuousHeatBathable for HeisenbergModel {
    fn heat_bath_sample(
        &self,
        neighbors: &[f64],
        beta: f64,
        rng: &mut impl Rng,
    ) -> SmallVec<[f64; 3]> {
        // Local field h = Σ s_j (sum over neighbors)
        let hx: f64 = neighbors.chunks(3).map(|c| c[0]).sum();
        let hy: f64 = neighbors.chunks(3).map(|c| c[1]).sum();
        let hz: f64 = neighbors.chunks(3).map(|c| c[2]).sum();
        let h_norm = (hx * hx + hy * hy + hz * hz).sqrt();

        // κ = βJ|h|
        let kappa = beta * self.j * h_norm;

        if kappa < 1e-12 {
            // Isotropic: sample uniformly on S²
            let z: f64 = rng.random::<f64>() * 2.0 - 1.0;
            let sin_theta = (1.0 - z * z).sqrt();
            let phi: f64 = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
            return smallvec![
                sin_theta * phi.cos(),
                sin_theta * phi.sin(),
                z,
            ];
        }

        // cosθ via inverse CDF: t = ln(u * 2sinh(κ) + e^{-κ}) / κ
        let u: f64 = rng.random::<f64>();
        let two_sinh_k = kappa.exp() - (-kappa).exp(); // 2sinh(κ) = e^κ - e^{-κ}
        let exp_neg_k = (-kappa).exp();
        let t = (u * two_sinh_k + exp_neg_k).ln() / kappa;

        let sin_theta = (1.0 - t * t).sqrt().max(0.0);
        let phi: f64 = rng.random::<f64>() * 2.0 * std::f64::consts::PI;

        // Spin in local frame where h aligns with z
        let (sx_local, sy_local) = (sin_theta * phi.cos(), sin_theta * phi.sin());

        // Rotate (0,0,1) → ĥ
        let (ux, uy, uz) = (hx / h_norm, hy / h_norm, hz / h_norm);
        let r = (ux * ux + uy * uy).sqrt(); // sqrt(1 - uz²)
        let (sx, sy, sz) = if r < 1e-12 {
            // μ ≈ ±z, trivial
            if uz > 0.0 {
                (sx_local, sy_local, t)
            } else {
                (sx_local, sy_local, -t)
            }
        } else {
            let inv_r = 1.0 / r;
            (
                -uy * inv_r * sx_local - ux * uz * inv_r * sy_local + ux * t,
                ux * inv_r * sx_local - uy * uz * inv_r * sy_local + uy * t,
                r * sy_local + uz * t,
            )
        };

        smallvec![sx, sy, sz]
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;

    #[test]
    fn test_ising_local_energy_ferro() {
        let lattice = build_chain(2, false);
        let model = IsingModel::new(1.0);
        // Two aligned spins
        let spins = vec![1.0, 1.0];
        // site 0 energy = -J * 1.0 * 1.0 = -1.0
        let e0 = model.local_energy(&spins, &lattice, 0, 1.0, &[1.0]);
        assert!((e0 + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_local_energy_antiferro() {
        let lattice = build_chain(2, false);
        let model = IsingModel::new(1.0);
        let spins = vec![1.0, -1.0];
        // site 0 energy = -J * 1.0 * (-1.0) = 1.0
        let e0 = model.local_energy(&spins, &lattice, 0, 1.0, &[1.0]);
        assert!((e0 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_total_energy() {
        let lattice = build_chain(4, true); // ring
        let model = IsingModel::new(1.0);
        // All up: 4 bonds, each -1, total -4
        let spins = vec![1.0, 1.0, 1.0, 1.0];
        let total = model.compute_total_energy(&spins, &lattice, 1.0);
        assert!((total + 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_magnetization() {
        let model = IsingModel::new(1.0);
        assert!((model.magnetization(&[1.0, -1.0, 1.0, -1.0]) - 0.0).abs() < 1e-10);
        assert!((model.magnetization(&[1.0, 1.0, 1.0, 1.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_fk_bond_prob() {
        let model = IsingModel::new(1.0);
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((model.fk_bond_probability(0.5) - expected).abs() < 1e-10);
    }

    // ── Potts tests ───────────────────────────────────────

    #[test]
    fn test_potts_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = PottsModel::new(1.0, 3);
        // Both in state 0
        let spins = vec![0.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[0.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_local_energy_misaligned() {
        let lattice = build_chain(2, false);
        let model = PottsModel::new(1.0, 3);
        let spins = vec![0.0, 1.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[0.0]);
        assert!((e - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_fk_bond_prob() {
        let model = PottsModel::new(1.0, 3);
        let expected = 1.0 - (-0.5_f64).exp();
        assert!((model.fk_bond_probability(0.5) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_ordered() {
        let model = PottsModel::new(1.0, 4);
        // All in state 2 → m = (4*4-4)/(4*3) = 12/12 = 1.0
        let spins = vec![2.0, 2.0, 2.0, 2.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_random() {
        let model = PottsModel::new(1.0, 4);
        // One in each state → m = (4*1-4)/(4*3) = 0
        let spins = vec![0.0, 1.0, 2.0, 3.0];
        assert!((model.magnetization(&spins) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_total_energy() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(1.0, 3);
        // All in state 0: 4 bonds, each -1, total -4
        let spins = vec![0.0, 0.0, 0.0, 0.0];
        let total = model.compute_total_energy(&spins, &lattice, 1.0);
        assert!((total + 4.0).abs() < 1e-10);
    }

    // ── XY tests ──────────────────────────────────────────

    #[test]
    fn test_xy_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = XYModel::new(1.0);
        // Both pointing along x: (1,0)
        let spins = vec![1.0, 0.0, 1.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[1.0, 0.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_local_energy_anti() {
        let lattice = build_chain(2, false);
        let model = XYModel::new(1.0);
        // Opposite: (1,0) and (-1,0)
        let spins = vec![1.0, 0.0, -1.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[1.0, 0.0]);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_aligned() {
        let model = XYModel::new(1.0);
        // Two aligned spins
        let spins = vec![1.0, 0.0, 1.0, 0.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_opposite() {
        let model = XYModel::new(1.0);
        let spins = vec![1.0, 0.0, -1.0, 0.0];
        assert!((model.magnetization(&spins) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_normalize_preserves_direction() {
        let model = XYModel::new(1.0);
        let mut v = vec![2.0, 0.0];
        model.normalize_spin(&mut v);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!(v[1].abs() < 1e-10);
    }

    #[test]
    fn test_xy_energy_total_ring() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(1.0);
        // All x-aligned, 4 bonds, each -J = -1, total = -4
        let spins = vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let total = model.compute_total_energy(&spins, &lattice, 1.0);
        assert!((total + 4.0).abs() < 1e-10);
    }

    // ── Heisenberg tests ──────────────────────────────────

    #[test]
    fn test_heisenberg_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = HeisenbergModel::new(1.0);
        // Both along z
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[0.0, 0.0, 1.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_local_energy_anti() {
        let lattice = build_chain(2, false);
        let model = HeisenbergModel::new(1.0);
        // Opposite along z
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0];
        let e = model.local_energy(&spins, &lattice, 0, 1.0, &[0.0, 0.0, 1.0]);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_magnetization_aligned() {
        let model = HeisenbergModel::new(1.0);
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_normalize() {
        let model = HeisenbergModel::new(1.0);
        let mut v = vec![0.0, 3.0, 0.0];
        model.normalize_spin(&mut v);
        assert!((v[1] - 1.0).abs() < 1e-10);
        assert!((v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_energy_total_ring() {
        let lattice = build_chain(4, true);
        let model = HeisenbergModel::new(1.0);
        // All z-aligned, 4 bonds, each -J = -1, total = -4
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let total = model.compute_total_energy(&spins, &lattice, 1.0);
        assert!((total + 4.0).abs() < 1e-10);
    }
}
