//! XXZ model implementation.
//!
//! H = J Σ [Δ S^z_i S^z_j + ½(S^+_i S^-_j + S^-_i S^+_j)]
//!
//! For |Δ| ≤ 1: diagonal and off-diagonal weights are both positive
//! with shift = J/4 per bond. The directed loop uses deterministic
//! diagonal ↔ off-diagonal switching (same as Heisenberg).
//!
//! For Δ > 1: requires bounce probabilities in the directed loop.
//! The bounce probability is (Δ - 1)/(Δ + 1) for anti-aligned spins.
//!
//! Supported limits:
//! - Δ = 1: Heisenberg XXX model
//! - Δ = 0: XY model (pure off-diagonal)
//! - Δ → ∞: Ising model (pure diagonal)

use crate::hilbert::{OpType, SpinHalfHS};
use crate::lattice::{BondType, Lattice};
use crate::sse::{LatticeQMC, SSEMonteCarlo};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand_xoshiro::Xoshiro256PlusPlus;

/// XXZ model H = J Σ [Δ S^z_i S^z_j + ½(S^+_i S^-_j + S^-_i S^+_j)]
pub struct XxzModel {
    /// Lattice topology
    lattice: Lattice,
    /// Inverse temperature
    beta: f64,
    /// Coupling constant
    j: f64,
    /// Anisotropy parameter Δ
    delta: f64,
}

impl XxzModel {
    /// Create new XXZ model.
    ///
    /// # Arguments
    /// * `lattice` - Lattice topology
    /// * `beta` - Inverse temperature
    /// * `j` - Coupling constant (J > 0 for antiferromagnetic)
    /// * `delta` - Anisotropy parameter
    ///   - Δ = 0: XY model
    ///   - Δ = 1: Heisenberg model
    ///   - Δ > 1: Ising-like (requires bounce in directed loop)
    ///   - |Δ| ≤ 1: bounce-free directed loop
    pub fn new(lattice: Lattice, beta: f64, j: f64, delta: f64) -> Self {
        XxzModel {
            lattice,
            beta,
            j,
            delta,
        }
    }

    /// Coupling constant J.
    pub fn j(&self) -> f64 {
        self.j
    }

    /// Anisotropy parameter Δ.
    pub fn delta(&self) -> f64 {
        self.delta
    }
}

impl LatticeQMC for XxzModel {
    fn lattice(&self) -> &Lattice {
        &self.lattice
    }
}

impl SSEMonteCarlo for XxzModel {
    type HilbertSpace = SpinHalfHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        // XXZ model: H = J Σ [Δ SzSz + ½(S+S- + S-S+)]
        //
        // Shifted Hamiltonian: H' = H + J*Δ/4 per bond
        // SSE expands in -H' + C:
        //   Aligned diagonal:     0 (no diagonal on aligned bonds)
        //   Anti-aligned diagonal: J*Δ/2
        //   Off-diagonal:          J/2
        //
        // The diagonal weight returned here is the actual matrix element
        // for anti-aligned spins (J*Δ/2), not the average (J*Δ/4).
        let w_diag = self.j * self.delta / 2.0;
        let w_offdiag = self.j * 0.5;
        vec![
            (OpType::Diagonal, w_diag),
            (OpType::OffDiagonal, w_offdiag),
        ]
    }

    fn hilbert_space(&self) -> &SpinHalfHS {
        static HS: SpinHalfHS = SpinHalfHS;
        &HS
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    /// Total diagonal shift: J*Δ/4 per bond.
    /// H' = H + J*Δ/4 * N_bonds, so E = J*Δ/4 * N_bonds/N - <n>/(β*N).
    fn diagonal_shift(&self) -> f64 {
        let n_bonds = self.lattice.n_bonds as f64;
        self.j * self.delta / 4.0 * n_bonds
    }

    /// XXZ anisotropy parameter — used by the directed loop algorithm.
    fn loop_parameter(&self) -> f64 {
        self.delta
    }
}

impl MonteCarlo for XxzModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — SSECore handles the actual sweep
    }
}

impl FromParams for XxzModel {
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
        let delta = params.get::<f64>("Delta").unwrap_or(1.0);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);

        let lattice = crate::lattice::builders::build_chain(n_sites, pbc);

        Ok(XxzModel::new(lattice, beta, j, delta))
    }
}
