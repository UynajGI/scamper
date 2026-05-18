//! Physics models — stateless formulas for spin systems.

use crate::lattice::Lattice;
use rand::Rng;
use rand::RngExt;

/// Physics model trait. Implementations define spin dimensionality, coupling constants,
/// and energy/proposal formulas. Models are stateless — all state lives in [`System`](crate::System).
pub trait Model: Send + Sync {
    /// Number of spin components per site (1 = Ising/Potts, 2 = XY, 3 = Heisenberg).
    fn spin_dim(&self) -> usize;

    /// Coupling constant J.
    fn coupling(&self) -> f64;

    /// Inverse temperature β = 1/(k_B T).
    fn beta(&self) -> f64;

    /// Energy contribution from site `i` interacting with its neighbors,
    /// given a **proposed** new spin. Hamiltonian convention: H = -J Σ_{⟨i,j⟩} s_i · s_j.
    ///
    /// `local_energy` should return the contribution to total energy from site `i`,
    /// summing over all bonds connected to `i`. The returned value replaces the old
    /// contribution directly (not a delta).
    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64;

    /// Propose a random new spin (all components). Returns a Vec of length `spin_dim()`.
    fn propose(&self, rng: &mut impl Rng) -> Vec<f64>;

    /// FK bond percolation probability for cluster algorithms.
    /// Default: `1 - exp(-2βJ)` for Ising-like models (H = -J Σ s_i s_j).
    fn fk_bond_probability(&self) -> f64 {
        1.0 - (-2.0 * self.coupling() * self.beta()).exp()
    }

    /// Total magnetization |M|/N from the current spin configuration.
    fn magnetization(&self, spins: &[f64]) -> f64;

    /// Random spin value for SW cluster assignment.
    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64;

    /// Opposite of a given spin (Wolff cluster flip).
    fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64;

    /// Compute initial energy for the full system.
    fn compute_total_energy(&self, spins: &[f64], lattice: &Lattice) -> f64 {
        let mut total = 0.0;
        let sd = self.spin_dim();
        for site in 0..lattice.n_sites {
            let proposed = spins[site * sd..(site + 1) * sd].to_vec();
            total += self.local_energy(spins, lattice, site, &proposed);
        }
        total / 2.0 // double-counted bonds
    }

    /// Normalize a spin vector to unit length (for XY/Heisenberg). Ising/Potts skip this.
    fn normalize_spin(&self, _spin: &mut [f64]) {}

    /// Random unit-length spin vector (for XY/Heisenberg). Ising/Potts override.
    fn random_spin(&self, rng: &mut impl Rng) -> Vec<f64> {
        vec![self.random_cluster_spin(rng)]
    }
}

// ── Ising Model ─────────────────────────────────────────────

/// Ising model: H = -J Σ_{⟨i,j⟩} σ_i σ_j, σ_i ∈ {±1}.
///
/// Parameters: `j` (coupling), `beta` (inverse temperature).
#[derive(Debug, Clone)]
pub struct IsingModel {
    pub j: f64,
    pub beta: f64,
}

impl IsingModel {
    pub fn new(j: f64, beta: f64) -> Self {
        Self { j, beta }
    }
}

impl Model for IsingModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64 {
        let s = proposed[0];
        let mut e = 0.0;
        for nb in &lattice.sites[site] {
            e += -self.j * s * spins[nb.target];
        }
        e
    }

    fn propose(&self, rng: &mut impl Rng) -> Vec<f64> {
        if rng.random::<bool>() {
            vec![1.0]
        } else {
            vec![-1.0]
        }
    }

    fn magnetization(&self, spins: &[f64]) -> f64 {
        let sum: f64 = spins.iter().sum();
        (sum / spins.len() as f64).abs()
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        if rng.random::<bool>() {
            1.0
        } else {
            -1.0
        }
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin
    }
}

// ── Potts Model ──────────────────────────────────────────────

/// q-state Potts model: H = -J Σ δ(s_i, s_j), s_i ∈ {0, 1, ..., q-1}.
///
/// FK bond probability: `1 - exp(-βJ)` (no factor 2, unlike Ising).
#[derive(Debug, Clone)]
pub struct PottsModel {
    pub j: f64,
    pub beta: f64,
    pub q: usize,
}

impl PottsModel {
    pub fn new(j: f64, beta: f64, q: usize) -> Self {
        assert!(q >= 2, "Potts q must be >= 2");
        Self { j, beta, q }
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

impl Model for PottsModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64 {
        let s = proposed[0] as usize;
        let mut e = 0.0;
        for nb in &lattice.sites[site] {
            let nb_state = spins[nb.target] as usize;
            if s == nb_state {
                e += -self.j;
            }
        }
        e
    }

    fn propose(&self, rng: &mut impl Rng) -> Vec<f64> {
        vec![rng.random_range(0..self.q) as f64]
    }

    /// Potts: p_FK = 1 - exp(-βJ)
    fn fk_bond_probability(&self) -> f64 {
        1.0 - (-self.j * self.beta).exp()
    }

    /// Potts magnetization: (q·max(n_k) - N) / (N·(q-1))
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let n = spins.len();
        if n == 0 {
            return 0.0;
        }
        let counts = self.state_counts(spins);
        let max_n = counts.iter().max().copied().unwrap_or(0);
        (self.q as f64 * max_n as f64 - n as f64) / (n as f64 * (self.q - 1) as f64)
    }

    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
        rng.random_range(0..self.q) as f64
    }

    fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64 {
        let current = spin as usize;
        let mut new = rng.random_range(0..self.q - 1);
        if new >= current {
            new += 1;
        }
        new as f64
    }
}

// ── XY Model ─────────────────────────────────────────────────

/// XY model: H = -J Σ cos(θ_i - θ_j) = -J Σ s_i · s_j, |s_i| = 1.
///
/// Each spin is a unit vector in 2D: `(cos θ, sin θ)`.
#[derive(Debug, Clone)]
pub struct XYModel {
    pub j: f64,
    pub beta: f64,
}

impl XYModel {
    pub fn new(j: f64, beta: f64) -> Self {
        Self { j, beta }
    }
}

impl Model for XYModel {
    fn spin_dim(&self) -> usize {
        2
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64 {
        let (sx, sy) = (proposed[0], proposed[1]);
        let mut e = 0.0;
        for nb in &lattice.sites[site] {
            let base = nb.target * 2;
            e += -self.j * (sx * spins[base] + sy * spins[base + 1]);
        }
        e
    }

    fn propose(&self, rng: &mut impl Rng) -> Vec<f64> {
        let theta: f64 = rng.random_range(0.0..std::f64::consts::TAU);
        vec![theta.cos(), theta.sin()]
    }

    fn magnetization(&self, spins: &[f64]) -> f64 {
        let (mut sx, mut sy) = (0.0, 0.0);
        for chunk in spins.chunks(2) {
            sx += chunk[0];
            sy += chunk[1];
        }
        let n = (spins.len() / 2) as f64;
        (sx * sx + sy * sy).sqrt() / n
    }

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        // Not used for continuous spins (cluster algorithms not applicable)
        0.0
    }

    fn opposite_spin(&self, _spin: f64, _rng: &mut impl Rng) -> f64 {
        0.0
    }

    fn normalize_spin(&self, spin: &mut [f64]) {
        let r = (spin[0] * spin[0] + spin[1] * spin[1]).sqrt();
        if r > 1e-12 {
            spin[0] /= r;
            spin[1] /= r;
        }
    }

    fn random_spin(&self, rng: &mut impl Rng) -> Vec<f64> {
        self.propose(rng)
    }
}

// ── Heisenberg Model ─────────────────────────────────────────

/// Heisenberg model: H = -J Σ s_i · s_j, |s_i| = 1.
///
/// Each spin is a unit vector in 3D. Uses Marsaglia's method for uniform sampling on S².
#[derive(Debug, Clone)]
pub struct HeisenbergModel {
    pub j: f64,
    pub beta: f64,
}

impl HeisenbergModel {
    pub fn new(j: f64, beta: f64) -> Self {
        Self { j, beta }
    }
}

impl Model for HeisenbergModel {
    fn spin_dim(&self) -> usize {
        3
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64 {
        let (sx, sy, sz) = (proposed[0], proposed[1], proposed[2]);
        let mut e = 0.0;
        for nb in &lattice.sites[site] {
            let base = nb.target * 3;
            e += -self.j * (sx * spins[base] + sy * spins[base + 1] + sz * spins[base + 2]);
        }
        e
    }

    fn propose(&self, rng: &mut impl Rng) -> Vec<f64> {
        // Marsaglia: sample (x,y) in unit disk, reject if outside
        let (x, y) = loop {
            let x: f64 = rng.random_range(-1.0..1.0);
            let y: f64 = rng.random_range(-1.0..1.0);
            if x * x + y * y <= 1.0 {
                break (x, y);
            }
        };
        let r = (x * x + y * y).sqrt();
        vec![
            2.0 * x * (1.0 - r * r).sqrt(),
            2.0 * y * (1.0 - r * r).sqrt(),
            1.0 - 2.0 * (x * x + y * y),
        ]
    }

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

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        0.0
    }

    fn opposite_spin(&self, _spin: f64, _rng: &mut impl Rng) -> f64 {
        0.0
    }

    fn normalize_spin(&self, spin: &mut [f64]) {
        let r = (spin[0] * spin[0] + spin[1] * spin[1] + spin[2] * spin[2]).sqrt();
        if r > 1e-12 {
            spin[0] /= r;
            spin[1] /= r;
            spin[2] /= r;
        }
    }

    fn random_spin(&self, rng: &mut impl Rng) -> Vec<f64> {
        self.propose(rng)
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
        let model = IsingModel::new(1.0, 1.0);
        // Two aligned spins
        let spins = vec![1.0, 1.0];
        // site 0 energy = -J * 1.0 * 1.0 = -1.0
        let e0 = model.local_energy(&spins, &lattice, 0, &[1.0]);
        assert!((e0 + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_local_energy_antiferro() {
        let lattice = build_chain(2, false);
        let model = IsingModel::new(1.0, 1.0);
        let spins = vec![1.0, -1.0];
        // site 0 energy = -J * 1.0 * (-1.0) = 1.0
        let e0 = model.local_energy(&spins, &lattice, 0, &[1.0]);
        assert!((e0 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_total_energy() {
        let lattice = build_chain(4, true); // ring
        let model = IsingModel::new(1.0, 1.0);
        // All up: 4 bonds, each -1, total -4
        let spins = vec![1.0, 1.0, 1.0, 1.0];
        let total = model.compute_total_energy(&spins, &lattice);
        assert!((total + 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_magnetization() {
        let model = IsingModel::new(1.0, 1.0);
        assert!((model.magnetization(&[1.0, -1.0, 1.0, -1.0]) - 0.0).abs() < 1e-10);
        assert!((model.magnetization(&[1.0, 1.0, 1.0, 1.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_fk_bond_prob() {
        let model = IsingModel::new(1.0, 0.5);
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((model.fk_bond_probability() - expected).abs() < 1e-10);
    }

    // ── Potts tests ───────────────────────────────────────

    #[test]
    fn test_potts_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = PottsModel::new(1.0, 1.0, 3);
        // Both in state 0
        let spins = vec![0.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, &[0.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_local_energy_misaligned() {
        let lattice = build_chain(2, false);
        let model = PottsModel::new(1.0, 1.0, 3);
        let spins = vec![0.0, 1.0];
        let e = model.local_energy(&spins, &lattice, 0, &[0.0]);
        assert!((e - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_fk_bond_prob() {
        let model = PottsModel::new(1.0, 0.5, 3);
        let expected = 1.0 - (-0.5_f64).exp();
        assert!((model.fk_bond_probability() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_ordered() {
        let model = PottsModel::new(1.0, 1.0, 4);
        // All in state 2 → m = (4*4-4)/(4*3) = 12/12 = 1.0
        let spins = vec![2.0, 2.0, 2.0, 2.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_magnetization_random() {
        let model = PottsModel::new(1.0, 1.0, 4);
        // One in each state → m = (4*1-4)/(4*3) = 0
        let spins = vec![0.0, 1.0, 2.0, 3.0];
        assert!((model.magnetization(&spins) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_potts_total_energy() {
        let lattice = build_chain(4, true);
        let model = PottsModel::new(1.0, 1.0, 3);
        // All in state 0: 4 bonds, each -1, total -4
        let spins = vec![0.0, 0.0, 0.0, 0.0];
        let total = model.compute_total_energy(&spins, &lattice);
        assert!((total + 4.0).abs() < 1e-10);
    }

    // ── XY tests ──────────────────────────────────────────

    #[test]
    fn test_xy_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = XYModel::new(1.0, 1.0);
        // Both pointing along x: (1,0)
        let spins = vec![1.0, 0.0, 1.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, &[1.0, 0.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_local_energy_anti() {
        let lattice = build_chain(2, false);
        let model = XYModel::new(1.0, 1.0);
        // Opposite: (1,0) and (-1,0)
        let spins = vec![1.0, 0.0, -1.0, 0.0];
        let e = model.local_energy(&spins, &lattice, 0, &[1.0, 0.0]);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_aligned() {
        let model = XYModel::new(1.0, 1.0);
        // Two aligned spins
        let spins = vec![1.0, 0.0, 1.0, 0.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_magnetization_opposite() {
        let model = XYModel::new(1.0, 1.0);
        let spins = vec![1.0, 0.0, -1.0, 0.0];
        assert!((model.magnetization(&spins) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_normalize_preserves_direction() {
        let model = XYModel::new(1.0, 1.0);
        let mut v = vec![2.0, 0.0];
        model.normalize_spin(&mut v);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!(v[1].abs() < 1e-10);
    }

    #[test]
    fn test_xy_energy_total_ring() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(1.0, 1.0);
        // All x-aligned, 4 bonds, each -J = -1, total = -4
        let spins = vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let total = model.compute_total_energy(&spins, &lattice);
        assert!((total + 4.0).abs() < 1e-10);
    }

    // ── Heisenberg tests ──────────────────────────────────

    #[test]
    fn test_heisenberg_local_energy_aligned() {
        let lattice = build_chain(2, false);
        let model = HeisenbergModel::new(1.0, 1.0);
        // Both along z
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let e = model.local_energy(&spins, &lattice, 0, &[0.0, 0.0, 1.0]);
        assert!((e + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_local_energy_anti() {
        let lattice = build_chain(2, false);
        let model = HeisenbergModel::new(1.0, 1.0);
        // Opposite along z
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, -1.0];
        let e = model.local_energy(&spins, &lattice, 0, &[0.0, 0.0, 1.0]);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_magnetization_aligned() {
        let model = HeisenbergModel::new(1.0, 1.0);
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        assert!((model.magnetization(&spins) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_normalize() {
        let model = HeisenbergModel::new(1.0, 1.0);
        let mut v = vec![0.0, 3.0, 0.0];
        model.normalize_spin(&mut v);
        assert!((v[1] - 1.0).abs() < 1e-10);
        assert!((v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heisenberg_energy_total_ring() {
        let lattice = build_chain(4, true);
        let model = HeisenbergModel::new(1.0, 1.0);
        // All z-aligned, 4 bonds, each -J = -1, total = -4
        let spins = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let total = model.compute_total_energy(&spins, &lattice);
        assert!((total + 4.0).abs() < 1e-10);
    }
}
