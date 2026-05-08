//! Models module — defines physics for classical Monte Carlo.

mod heisenberg;
mod ising;
mod ising_2d;
mod potts;
mod xy;

pub use heisenberg::HeisenbergModel;
pub use ising::IsingModel;
pub use ising_2d::IsingModel2D;
pub use potts::PottsModel;
pub use xy::XYModel;

use crate::lattice::LatticeMC;
use rand::Rng;

/// Method trait for Classical Monte Carlo.
/// Models implement this trait to define their physics.
pub trait ModelMC: LatticeMC {
    /// Spin dimension: Ising=1, XY=2, Heisenberg=3
    fn spin_dim(&self) -> usize;

    /// Coupling constant J (for bond probability in cluster algorithms)
    fn coupling(&self) -> f64;

    /// Simulation inverse temperature
    fn beta(&self) -> f64;

    /// Propose a spin flip at the given site.
    /// Returns `(old_spin_value, proposed_new_spin_value)`.
    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64);

    /// Energy change when flipping site from old to new state
    fn local_energy_change(&self, site: usize, old: f64, new: f64) -> f64;

    /// Propose a full spin flip. Returns (old_spin_vec, new_spin_vec).
    /// Default: wraps propose_flip for scalar spin models.
    fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let (old, new) = self.propose_flip(site, rng);
        (vec![old], vec![new])
    }

    /// Energy change for a vector spin flip.
    /// Default: wraps local_energy_change for scalar spin models.
    fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
        self.local_energy_change(site, old[0], new[0])
    }

    /// Total energy of current configuration
    fn total_energy(&self) -> f64;

    /// Access spin configuration
    fn spins(&self) -> &[f64];

    /// Mutable access to spin configuration
    fn spins_mut(&mut self) -> &mut [f64];

    /// Magnetization of the current configuration.
    /// Model-specific: Ising=|Σs_i|/N, Potts=(q·max(n_k)-N)/(N·(q-1)), XY=|Σ(cosθ,sinθ)|/N
    fn magnetization(&self) -> f64;

    /// Random spin value for cluster assignment (SW algorithm).
    /// Only meaningful for discrete spin models (Ising, Potts).
    fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64;

    /// Opposite of a given spin (Wolff cluster flip).
    /// Only meaningful for discrete spin models with reflection symmetry.
    fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64;

    /// FK bond percolation probability for cluster algorithms.
    /// Default: 1 - exp(-2*beta*J) for Ising-like models (H = -J Σ s_i s_j).
    /// Potts overrides to 1 - exp(-beta*J) (H = -J Σ δ(s_i, s_j)).
    fn fk_bond_probability(&self) -> f64 {
        1.0 - (-2.0 * self.coupling() * self.beta()).exp()
    }

    /// Raw spin configuration snapshot as Vec<f64>.
    fn snapshot(&self) -> Vec<f64> {
        self.spins().to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::models::ising::IsingModel;
    use crate::models::potts::PottsModel;

    fn _requires_model_mc<M: ModelMC>(_m: &M) {}

    #[test]
    fn test_model_mc_trait_bound() {
        fn _compile_check() {}
    }

    #[test]
    fn test_ising_fk_bond_probability() {
        // Ising: H = -J Σ s_i s_j → p_FK = 1 - exp(-2βJ)
        let lattice = build_chain(4, true);
        let model = IsingModel::new(lattice, 0.5_f64, 1.0_f64);
        let expected = 1.0_f64 - (-2.0_f64 * 1.0_f64 * 0.5_f64).exp();
        assert!((model.fk_bond_probability() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_potts_fk_bond_probability() {
        // Potts: H = -J Σ δ(s_i, s_j) → p_FK = 1 - exp(-βJ)
        let lattice = build_chain(4, true);
        let model = PottsModel::new(lattice, 0.5_f64, 1.0_f64, 3);
        let expected = 1.0_f64 - (-1.0_f64 * 0.5_f64).exp();
        assert!((model.fk_bond_probability() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_ising_vs_potts_fk_probability_differs() {
        // Verify that Ising and Potts give different FK probabilities
        // for the same β and J parameters
        let lattice = build_chain(4, true);
        let ising = IsingModel::new(lattice.clone(), 1.0, 1.0);
        let potts = PottsModel::new(lattice, 1.0, 1.0, 3);

        let p_ising = ising.fk_bond_probability();
        let p_potts = potts.fk_bond_probability();

        // Ising: 1 - exp(-2) ≈ 0.8647
        // Potts: 1 - exp(-1) ≈ 0.6321
        assert!((p_ising - p_potts).abs() > 0.1,
            "Ising and Potts FK probabilities should differ significantly");
        assert!(p_ising > p_potts,
            "Ising FK probability should be larger than Potts (factor of 2 in exponent)");
    }
}
