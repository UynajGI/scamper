//! Built-in lattice-spin models.
//!
//! Ising and Potts use discrete scalar storage.  XY and Heisenberg are aliases
//! of the const-generic [`ONModel`], which also enables arbitrary O(N) spins.

use crate::core::r#move::Spin;
use crate::lattice::graph::{Bond, CsrLattice};
use crate::lattice::interaction::{
    ClusterAuxiliary, ClusterModel, ContinuousHeatBathable, HeatBathable, Initializable,
    LocalFieldModel, Measurable, PairInteraction, Proposable,
};
use rand::{Rng, RngExt};
use smallvec::{smallvec, SmallVec};

#[inline]
fn clamp_probability(value: f64) -> f64 {
    if value.is_finite() {
        value.clamp(0.0, 1.0)
    } else if value.is_sign_positive() {
        1.0
    } else {
        0.0
    }
}

// ── Ising ───────────────────────────────────────────────────

/// Ising model `H = -J Σ_e w_e σ_i σ_j`, `σ_i ∈ {±1}`.
#[derive(Debug, Clone)]
pub struct IsingModel {
    pub j: f64,
}

impl IsingModel {
    pub fn new(j: f64) -> Self {
        assert!(j.is_finite(), "Ising coupling must be finite");
        Self { j }
    }

    pub const fn spin_dim(&self) -> usize {
        1
    }

    pub const fn coupling(&self) -> f64 {
        self.j
    }
}

impl PairInteraction for IsingModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        -self.j * bond.weight * left[0] * right[0]
    }

    fn validate_spin(&self, spin: &[f64]) -> Result<(), String> {
        if spin.len() != 1 {
            return Err(format!(
                "Ising spin dimension must be 1, got {}",
                spin.len()
            ));
        }
        if spin[0] != -1.0 && spin[0] != 1.0 {
            return Err(format!(
                "Ising spin must be exactly -1 or +1, got {}",
                spin[0]
            ));
        }
        Ok(())
    }
}

impl Initializable for IsingModel {
    fn random_spin(&self, rng: &mut impl Rng) -> Spin {
        if rng.random::<bool>() {
            smallvec![1.0]
        } else {
            smallvec![-1.0]
        }
    }

    fn ordered_spin(&self) -> Spin {
        smallvec![1.0]
    }
}

impl Proposable for IsingModel {
    fn propose(&self, rng: &mut impl Rng) -> Spin {
        self.random_spin(rng)
    }

    fn propose_from(&self, current: &[f64], _rng: &mut impl Rng) -> Spin {
        smallvec![-current[0]]
    }
}

impl Measurable for IsingModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        if spins.is_empty() {
            return 0.0;
        }
        spins.iter().sum::<f64>().abs() / spins.len() as f64
    }
}

impl ClusterModel for IsingModel {
    fn wolff_auxiliary(&self, seed_spin: &[f64], _rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::DiscreteTarget(-seed_spin[0])
    }

    fn sw_bond_auxiliary(&self, _rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::None
    }

    fn sw_cluster_auxiliary(
        &self,
        _representative_spin: &[f64],
        _bond_auxiliary: &ClusterAuxiliary,
        rng: &mut impl Rng,
    ) -> ClusterAuxiliary {
        ClusterAuxiliary::DiscreteTarget(if rng.random::<bool>() { 1.0 } else { -1.0 })
    }

    fn cluster_bond_probability(
        &self,
        left: &[f64],
        right: &[f64],
        bond: &Bond,
        _auxiliary: &ClusterAuxiliary,
        beta: f64,
    ) -> f64 {
        let coupling = self.j * bond.weight;
        assert!(
            coupling >= 0.0,
            "ferromagnetic Ising cluster updates require J * bond.weight >= 0"
        );
        if coupling == 0.0 || (left[0] - right[0]).abs() > 1e-12 {
            return 0.0;
        }
        clamp_probability(1.0 - (-2.0 * beta * coupling).exp())
    }

    fn transform_cluster_spin(&self, spin: &[f64], auxiliary: &ClusterAuxiliary) -> Spin {
        match auxiliary {
            ClusterAuxiliary::DiscreteTarget(target) => smallvec![*target],
            ClusterAuxiliary::Identity | ClusterAuxiliary::None => SmallVec::from_slice(spin),
            ClusterAuxiliary::Reflection(_) => SmallVec::from_slice(spin),
        }
    }
}

impl HeatBathable for IsingModel {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin {
        let field: f64 = lattice
            .incidences(site)
            .filter_map(|(neighbor, edge_id)| {
                let edge = lattice.edges[edge_id];
                (edge.source != edge.target).then_some(self.j * edge.weight * spins[neighbor])
            })
            .sum();
        let x = 2.0 * beta * field;
        let p_plus = if x >= 0.0 {
            1.0 / (1.0 + (-x).exp())
        } else {
            x.exp() / (1.0 + x.exp())
        };
        if rng.random::<f64>() < p_plus {
            smallvec![1.0]
        } else {
            smallvec![-1.0]
        }
    }
}

// ── Potts ───────────────────────────────────────────────────

/// q-state Potts model `H = -J Σ_e w_e δ(s_i,s_j)`.
#[derive(Debug, Clone)]
pub struct PottsModel {
    pub j: f64,
    pub q: usize,
}

impl PottsModel {
    pub fn new(j: f64, q: usize) -> Self {
        assert!(j.is_finite(), "Potts coupling must be finite");
        assert!(q >= 2, "Potts q must be >= 2");
        Self { j, q }
    }

    pub const fn spin_dim(&self) -> usize {
        1
    }

    pub const fn coupling(&self) -> f64 {
        self.j
    }

    fn random_other_state(&self, current: usize, rng: &mut impl Rng) -> usize {
        let mut state = rng.random_range(0..self.q - 1);
        if state >= current {
            state += 1;
        }
        state
    }
}

impl PairInteraction for PottsModel {
    fn spin_dim(&self) -> usize {
        1
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        if left[0] as usize == right[0] as usize {
            -self.j * bond.weight
        } else {
            0.0
        }
    }

    fn validate_spin(&self, spin: &[f64]) -> Result<(), String> {
        if spin.len() != 1 {
            return Err(format!(
                "Potts spin dimension must be 1, got {}",
                spin.len()
            ));
        }
        let value = spin[0];
        if !value.is_finite() || value.fract() != 0.0 || value < 0.0 || value >= self.q as f64 {
            return Err(format!(
                "Potts spin must be an integer in [0, {}), got {value}",
                self.q
            ));
        }
        Ok(())
    }
}

impl Initializable for PottsModel {
    fn random_spin(&self, rng: &mut impl Rng) -> Spin {
        smallvec![rng.random_range(0..self.q) as f64]
    }

    fn ordered_spin(&self) -> Spin {
        smallvec![0.0]
    }
}

impl Proposable for PottsModel {
    fn propose(&self, rng: &mut impl Rng) -> Spin {
        self.random_spin(rng)
    }

    fn propose_from(&self, current: &[f64], rng: &mut impl Rng) -> Spin {
        smallvec![self.random_other_state(current[0] as usize, rng) as f64]
    }
}

impl Measurable for PottsModel {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        if spins.is_empty() {
            return 0.0;
        }
        let mut counts = vec![0usize; self.q];
        for &spin in spins {
            let state = spin as usize;
            if state < self.q {
                counts[state] += 1;
            }
        }
        let largest = counts.into_iter().max().unwrap_or(0);
        let n = spins.len() as f64;
        (self.q as f64 * largest as f64 - n) / (n * (self.q - 1) as f64)
    }
}

impl ClusterModel for PottsModel {
    fn wolff_auxiliary(&self, seed_spin: &[f64], rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::DiscreteTarget(self.random_other_state(seed_spin[0] as usize, rng) as f64)
    }

    fn sw_bond_auxiliary(&self, _rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::None
    }

    fn sw_cluster_auxiliary(
        &self,
        _representative_spin: &[f64],
        _bond_auxiliary: &ClusterAuxiliary,
        rng: &mut impl Rng,
    ) -> ClusterAuxiliary {
        ClusterAuxiliary::DiscreteTarget(rng.random_range(0..self.q) as f64)
    }

    fn cluster_bond_probability(
        &self,
        left: &[f64],
        right: &[f64],
        bond: &Bond,
        _auxiliary: &ClusterAuxiliary,
        beta: f64,
    ) -> f64 {
        let coupling = self.j * bond.weight;
        assert!(
            coupling >= 0.0,
            "ferromagnetic Potts cluster updates require J * bond.weight >= 0"
        );
        if coupling == 0.0 || left[0] as usize != right[0] as usize {
            return 0.0;
        }
        clamp_probability(1.0 - (-beta * coupling).exp())
    }

    fn transform_cluster_spin(&self, spin: &[f64], auxiliary: &ClusterAuxiliary) -> Spin {
        match auxiliary {
            ClusterAuxiliary::DiscreteTarget(target) => smallvec![*target],
            ClusterAuxiliary::Identity | ClusterAuxiliary::None => SmallVec::from_slice(spin),
            ClusterAuxiliary::Reflection(_) => SmallVec::from_slice(spin),
        }
    }
}

impl HeatBathable for PottsModel {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin {
        let mut log_weights = vec![0.0f64; self.q];
        for (neighbor, edge_id) in lattice.incidences(site) {
            let edge = lattice.edges[edge_id];
            if edge.source == edge.target {
                continue;
            }
            let state = spins[neighbor] as usize;
            if state < self.q {
                log_weights[state] += beta * self.j * edge.weight;
            }
        }

        let max_log = log_weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let mut total = 0.0;
        for value in &mut log_weights {
            *value = (*value - max_log).exp();
            total += *value;
        }

        let mut threshold = rng.random::<f64>() * total;
        for (state, weight) in log_weights.into_iter().enumerate() {
            threshold -= weight;
            if threshold <= 0.0 {
                return smallvec![state as f64];
            }
        }
        smallvec![(self.q - 1) as f64]
    }
}

// ── O(N) ────────────────────────────────────────────────────

/// Ferromagnetic O(N) model `H = -J Σ_e w_e s_i·s_j`, `|s_i|=1`.
#[derive(Debug, Clone)]
pub struct ONModel<const D: usize> {
    pub j: f64,
}

impl<const D: usize> ONModel<D> {
    pub fn new(j: f64) -> Self {
        assert!(D >= 2, "O(N) spin dimension must be >= 2");
        assert!(j.is_finite(), "O(N) coupling must be finite");
        Self { j }
    }

    pub const fn spin_dim(&self) -> usize {
        D
    }

    pub const fn coupling(&self) -> f64 {
        self.j
    }

    fn random_unit_vector(&self, rng: &mut impl Rng) -> Spin {
        // Normalize independent standard-normal components generated in pairs
        // with Box-Muller.  This is uniform on S^(D-1) for every D >= 2.
        let mut vector = SmallVec::<[f64; 3]>::with_capacity(D);
        while vector.len() < D {
            let u1 = rng.random::<f64>().max(f64::MIN_POSITIVE);
            let u2 = rng.random::<f64>();
            let radius = (-2.0 * u1.ln()).sqrt();
            let angle = std::f64::consts::TAU * u2;
            vector.push(radius * angle.cos());
            if vector.len() < D {
                vector.push(radius * angle.sin());
            }
        }
        normalize(&mut vector);
        vector
    }
}

pub type XYModel = ONModel<2>;
pub type HeisenbergModel = ONModel<3>;

#[inline]
fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn normalize(spin: &mut [f64]) {
    let norm = spin.iter().map(|value| value * value).sum::<f64>().sqrt();
    if norm > 1e-15 {
        for value in spin {
            *value /= norm;
        }
    }
}

impl<const D: usize> PairInteraction for ONModel<D> {
    fn spin_dim(&self) -> usize {
        D
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        -self.j * bond.weight * dot(left, right)
    }

    fn validate_spin(&self, spin: &[f64]) -> Result<(), String> {
        if spin.len() != D {
            return Err(format!(
                "O(N) spin dimension must be {D}, got {}",
                spin.len()
            ));
        }
        if spin.iter().any(|component| !component.is_finite()) {
            return Err("O(N) spin contains a non-finite component".to_string());
        }
        let norm_squared = spin
            .iter()
            .map(|component| component * component)
            .sum::<f64>();
        let tolerance = 1e-10 * D as f64;
        if (norm_squared - 1.0).abs() > tolerance {
            return Err(format!(
                "O(N) spin must have unit norm: norm^2={norm_squared}"
            ));
        }
        Ok(())
    }
}

impl<const D: usize> Initializable for ONModel<D> {
    fn random_spin(&self, rng: &mut impl Rng) -> Spin {
        self.random_unit_vector(rng)
    }

    fn ordered_spin(&self) -> Spin {
        let mut spin = SmallVec::from_elem(0.0, D);
        spin[0] = 1.0;
        spin
    }
}

impl<const D: usize> Proposable for ONModel<D> {
    fn propose(&self, rng: &mut impl Rng) -> Spin {
        self.random_unit_vector(rng)
    }

    fn normalize_spin(&self, spin: &mut [f64]) {
        normalize(spin);
    }
}

impl<const D: usize> Measurable for ONModel<D> {
    fn magnetization(&self, spins: &[f64]) -> f64 {
        let n_sites = spins.len() / D;
        if n_sites == 0 {
            return 0.0;
        }
        let mut sum = vec![0.0; D];
        for spin in spins.chunks_exact(D) {
            for (component, value) in sum.iter_mut().zip(spin) {
                *component += *value;
            }
        }
        sum.iter().map(|value| value * value).sum::<f64>().sqrt() / n_sites as f64
    }
}

impl<const D: usize> ClusterModel for ONModel<D> {
    fn wolff_auxiliary(&self, _seed_spin: &[f64], rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::Reflection(self.random_unit_vector(rng))
    }

    fn sw_bond_auxiliary(&self, rng: &mut impl Rng) -> ClusterAuxiliary {
        ClusterAuxiliary::Reflection(self.random_unit_vector(rng))
    }

    fn sw_cluster_auxiliary(
        &self,
        _representative_spin: &[f64],
        bond_auxiliary: &ClusterAuxiliary,
        rng: &mut impl Rng,
    ) -> ClusterAuxiliary {
        if rng.random::<bool>() {
            bond_auxiliary.clone()
        } else {
            ClusterAuxiliary::Identity
        }
    }

    fn cluster_bond_probability(
        &self,
        left: &[f64],
        right: &[f64],
        bond: &Bond,
        auxiliary: &ClusterAuxiliary,
        beta: f64,
    ) -> f64 {
        let ClusterAuxiliary::Reflection(direction) = auxiliary else {
            return 0.0;
        };
        let coupling = self.j * bond.weight;
        assert!(
            coupling >= 0.0,
            "ferromagnetic O(N) cluster updates require J * bond.weight >= 0"
        );
        let product = dot(left, direction) * dot(right, direction);
        if coupling == 0.0 || product <= 0.0 {
            return 0.0;
        }
        // Correct Wolff embedded-Ising activation probability.
        clamp_probability(1.0 - (-2.0 * beta * coupling * product).exp())
    }

    fn transform_cluster_spin(&self, spin: &[f64], auxiliary: &ClusterAuxiliary) -> Spin {
        match auxiliary {
            ClusterAuxiliary::Reflection(direction) => {
                let projection = dot(spin, direction);
                let mut reflected = SmallVec::from_slice(spin);
                for (value, normal) in reflected.iter_mut().zip(direction) {
                    *value -= 2.0 * projection * normal;
                }
                normalize(&mut reflected);
                reflected
            }
            ClusterAuxiliary::Identity | ClusterAuxiliary::None => SmallVec::from_slice(spin),
            ClusterAuxiliary::DiscreteTarget(_) => SmallVec::from_slice(spin),
        }
    }
}

impl<const D: usize> LocalFieldModel for ONModel<D> {
    fn local_field(&self, spins: &[f64], lattice: &CsrLattice, site: usize, output: &mut [f64]) {
        assert_eq!(output.len(), D);
        output.fill(0.0);
        for (neighbor, edge_id) in lattice.incidences(site) {
            let edge = lattice.edges[edge_id];
            if edge.source == edge.target {
                continue;
            }
            let scale = self.j * edge.weight;
            let base = neighbor * D;
            for component in 0..D {
                output[component] += scale * spins[base + component];
            }
        }
    }
}

impl ContinuousHeatBathable for ONModel<2> {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin {
        let mut field = [0.0; 2];
        self.local_field(spins, lattice, site, &mut field);
        let norm = field.iter().map(|value| value * value).sum::<f64>().sqrt();
        if norm < 1e-15 {
            return self.random_unit_vector(rng);
        }
        let theta = sample_von_mises(rng, beta * norm) + field[1].atan2(field[0]);
        smallvec![theta.cos(), theta.sin()]
    }
}

impl ContinuousHeatBathable for ONModel<3> {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin {
        let mut field = [0.0; 3];
        self.local_field(spins, lattice, site, &mut field);
        let norm = field.iter().map(|value| value * value).sum::<f64>().sqrt();
        let kappa = beta * norm;
        if kappa < 1e-12 {
            return self.random_unit_vector(rng);
        }

        // Stable inverse CDF for exp(kappa * cos(theta)).
        let u = rng.random::<f64>().max(f64::MIN_POSITIVE);
        let cos_theta = 1.0 + (u + (1.0 - u) * (-2.0 * kappa).exp()).ln() / kappa;
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let phi = std::f64::consts::TAU * rng.random::<f64>();
        let local_x = sin_theta * phi.cos();
        let local_y = sin_theta * phi.sin();

        let ux = field[0] / norm;
        let uy = field[1] / norm;
        let uz = field[2] / norm;
        let transverse = (ux * ux + uy * uy).sqrt();
        let (x, y, z) = if transverse < 1e-12 {
            if uz >= 0.0 {
                (local_x, local_y, cos_theta)
            } else {
                (local_x, local_y, -cos_theta)
            }
        } else {
            let inverse = 1.0 / transverse;
            (
                -uy * inverse * local_x - ux * uz * inverse * local_y + ux * cos_theta,
                ux * inverse * local_x - uy * uz * inverse * local_y + uy * cos_theta,
                transverse * local_y + uz * cos_theta,
            )
        };
        smallvec![x, y, z]
    }
}

/// Best-Fisher rejection sampler for `VM(0, kappa)`.
fn sample_von_mises(rng: &mut impl Rng, kappa: f64) -> f64 {
    if kappa < 1e-8 {
        return rng.random::<f64>() * std::f64::consts::TAU;
    }
    let tau = 1.0 + (1.0 + 4.0 * kappa * kappa).sqrt();
    let rho = (tau - (2.0 * tau).sqrt()) / (2.0 * kappa);
    let r = (1.0 + rho * rho) / (2.0 * rho);

    let cosine = loop {
        let z = (std::f64::consts::PI * rng.random::<f64>()).cos();
        let value = (1.0 + r * z) / (r + z);
        let c = kappa * (r - value);
        let u = rng.random::<f64>().max(f64::MIN_POSITIVE);
        if c * (2.0 - c) > u || (c / u).ln() + 1.0 >= c {
            break value.clamp(-1.0, 1.0);
        }
    };

    if rng.random::<bool>() {
        cosine.acos()
    } else {
        -cosine.acos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::graph::{build_chain, BondType};
    use crate::lattice::interaction::Hamiltonian;
    use rand::SeedableRng;

    #[test]
    fn physical_edges_are_counted_once() {
        let lattice = build_chain(4, true);
        let model = IsingModel::new(1.0);
        assert_eq!(model.compute_total_energy(&[1.0; 4], &lattice, 1.0), -4.0);
    }

    #[test]
    fn weighted_bond_energy() {
        let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 2.5)]);
        let model = IsingModel::new(2.0);
        assert_eq!(model.compute_total_energy(&[1.0, 1.0], &lattice, 1.0), -5.0);
    }

    #[test]
    fn arbitrary_on_dimension_is_normalized() {
        let model = ONModel::<7>::new(1.0);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        let spin = model.random_spin(&mut rng);
        let norm = spin.iter().map(|value| value * value).sum::<f64>();
        assert_eq!(spin.len(), 7);
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn on_bond_probability_depends_on_projection() {
        let model = XYModel::new(1.0);
        let bond = Bond::new(0, 1, BondType::Generic, 1.0);
        let aux = ClusterAuxiliary::Reflection(smallvec![1.0, 0.0]);
        let parallel = model.cluster_bond_probability(&[1.0, 0.0], &[1.0, 0.0], &bond, &aux, 1.0);
        let perpendicular =
            model.cluster_bond_probability(&[0.0, 1.0], &[0.0, 1.0], &bond, &aux, 1.0);
        assert!(parallel > perpendicular);
        assert_eq!(perpendicular, 0.0);
    }
}
