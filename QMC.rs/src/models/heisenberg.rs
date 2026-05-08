//! Heisenberg model implementation.

use crate::hilbert::{OpType, SpinHalfHS};
use crate::lattice::{BondType, Lattice};
use crate::sse::{LatticeQMC, SSEMonteCarlo};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Heisenberg model H = J Σ S_i · S_j
pub struct HeisenbergModel {
    /// Lattice topology
    lattice: Lattice,
    /// Inverse temperature
    beta: f64,
    /// Coupling constant
    j: f64,
}

impl HeisenbergModel {
    /// Create new Heisenberg model.
    pub fn new(lattice: Lattice, beta: f64, j: f64) -> Self {
        HeisenbergModel { lattice, beta, j }
    }

    /// Access coupling constant J.
    pub fn j(&self) -> f64 {
        self.j
    }
}

// LatticeQMC implementation
impl LatticeQMC for HeisenbergModel {
    fn lattice(&self) -> &Lattice {
        &self.lattice
    }
}

// SSEMonteCarlo implementation - defines physics
impl SSEMonteCarlo for HeisenbergModel {
    type HilbertSpace = SpinHalfHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        // Shifted Heisenberg: H_b = J/2 (S^z_i S^z_j + 1/4) + J/2 (S^+_i S^-_j + S^-_i S^+_j)
        // Both diagonal and off-diagonal have weight J/2 for ALL spin configurations.
        // The diagonal shift per bond is J/8 (constant added to Hamiltonian).
        vec![
            (OpType::Diagonal, self.j * 0.5),
            (OpType::OffDiagonal, self.j * 0.5),
        ]
    }

    fn hilbert_space(&self) -> &SpinHalfHS {
        static HS: SpinHalfHS = SpinHalfHS;
        &HS
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    /// Total diagonal shift for the shifted Hamiltonian.
    /// We expand in -H' = -H + C, so E = C_total/N - ⟨n⟩/(β*N).
    /// The shift per bond is J/4, ensuring diagonal operators exist for anti-aligned spins.
    fn diagonal_shift(&self) -> f64 {
        let n_bonds = self.lattice.n_bonds as f64;
        self.j * n_bonds / 4.0
    }
}

// MonteCarlo implementation - required by FromParams
// Note: The actual simulation is done by SSECore, this is just a placeholder
impl MonteCarlo for HeisenbergModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder - SSECore handles the actual sweep
        // This should not be called directly
    }
}

// FromParams implementation - allows creation from Params
impl FromParams for HeisenbergModel {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        // Required parameters
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

        // Optional parameters with defaults
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);

        // Build lattice (1D chain for now)
        let lattice = crate::lattice::builders::build_chain(n_sites, pbc);

        Ok(HeisenbergModel::new(lattice, beta, j))
    }
}