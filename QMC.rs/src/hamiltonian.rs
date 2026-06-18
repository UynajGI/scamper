//! Quantum Hamiltonian abstraction for lattice spin models.
//!
//! Designed to be reused by both discrete-time and (future) continuous-time
//! worm algorithms. The worm needs, per bond, the **diagonal** and
//! **off-diagonal** matrix elements of the Trotter-decomposed Boltzmann
//! operator `B(Δτ) = exp(-Δτ H_bond)`.
//!
//! For S = 1/2 models, the spin state at a site is encoded as a `u8`
//! (`0 = ↓`, `1 = ↑`).

use std::f64::consts;

/// Quantum spin Hamiltonian for the worm algorithm.
///
/// All matrix elements are of the **Trotter bond operator**
/// `B(Δτ) = exp(-Δτ H_bond)` for a single pair of neighboring sites.
///
/// `Δτ` is passed in (not stored) so a single model struct can serve
/// parallel-tempering runs across β, and so the worm can probe weights at
/// the actual slice width of its `SpaceTimeConfig`.
pub trait QuantumHamiltonian: Send + Sync {
    /// Diagonal matrix element ⟨sᵢ sⱼ| B(Δτ) |sᵢ sⱼ⟩.
    ///
    /// For aligned-spins configurations. Always ≥ 0.
    fn bond_diagonal(&self, s_i: u8, s_j: u8, dtau: f64) -> f64;

    /// Off-diagonal matrix element ⟨flip(sᵢ,sⱼ)| B(Δτ) |sᵢ sⱼ⟩.
    ///
    /// For S = 1/2 Heisenberg/XY this is the spin-exchange term
    /// `½(S⁺S⁻ + S⁻S⁺)`. Returns 0 if the pair can't be connected by a
    /// single bond operator (e.g. ⟨↑↑|B|↓↓⟩ = 0 for Heisenberg).
    fn bond_offdiag(&self, s_i: u8, s_j: u8, dtau: f64) -> f64;

    /// Physical energy of one bond in the state |sᵢ sⱼ⟩ — the **diagonal**
    /// energy estimator (classical/Néel contribution only). Sign convention:
    /// matches H.
    ///
    /// For the full quantum energy estimator including spin-exchange
    /// fluctuations, use [`QuantumHamiltonian::bond_energy_estimator`].
    fn energy_per_bond(&self, s_i: u8, s_j: u8) -> f64;

    /// Quantum path-integral energy estimator for one bond.
    ///
    /// `is_kink` is true when this bond carries an off-diagonal matrix
    /// element (i.e. for an antiparallel spin pair, or a temporal bond
    /// across slices with differing spins). Returns
    /// `−(d/dΔτ) ln⟨sᵢsⱼ| B(Δτ) |sᵢsⱼ⟩` evaluated in the relevant sector.
    ///
    /// Default implementation uses the diagonal-only value — correct for
    /// classical/Ising-type models but **missing quantum fluctuations** for
    /// Heisenberg. Models with non-trivial off-diagonal sectors should
    /// override this.
    fn bond_energy_estimator(&self, s_i: u8, s_j: u8, dtau: f64, is_kink: bool) -> f64 {
        let _ = (dtau, is_kink);
        self.energy_per_bond(s_i, s_j)
    }

    /// Coupling constant J (informational; PT and reweighting may use it).
    fn coupling(&self) -> f64;
}

// ── Heisenberg S = 1/2 chain ─────────────────────────────────

/// Heisenberg model: H = -J Σ_{⟨i,j⟩} Sᵢ·Sⱼ, S = 1/2.
///
/// The Trotter bond operator is `B(Δτ) = exp(Δτ J Sᵢ·Sⱼ)` (note: weight is
/// `exp(-Δτ H)`, and `H_bond = -J Sᵢ·Sⱼ`, so the exponent is `+Δτ J Sᵢ·Sⱼ`).
///
/// `Sᵢ·Sⱼ` is a 4×4 operator on the two-spin Hilbert space. Its eigenstructure:
/// - Triplet `|↑↑⟩`, `|↓↓⟩` (and the symmetric combo of `|↑↓⟩+|↓↑⟩`):
///   eigenvalue `+1/4` → `B = exp(Δτ J/4)`.
/// - Singlet `(|↑↓⟩ − |↓↑⟩)/√2`: eigenvalue `−3/4` → `B = exp(−3Δτ J/4)`.
///
/// In the product basis `{|↑↑⟩, |↓↓⟩, |↑↓⟩, |↓↑⟩}` the nontrivial 2×2 block on
/// the antiparallel subspace `{|↑↓⟩, |↓↑⟩}` has matrix elements:
/// - `⟨↑↓| B |↑↓⟩ = ⟨↓↑| B |↓↑⟩ = ½(exp(ΔτJ/4) + exp(−3ΔτJ/4))`
/// - `⟨↑↓| B |↓↑⟩ = ⟨↓↑| B |↑↓⟩ = ½(exp(ΔτJ/4) − exp(−3ΔτJ/4))`
///
/// The worm / single-flip kernel reads these as:
/// - **Parallel** pair `↑↑`/`↓↓`: bond weight = `exp(Δτ J/4)` (diagonal).
/// - **Antiparallel** pair `↑↓`: the bond contributes the diagonal element
///   `½(exp(ΔτJ/4) + exp(−3ΔτJ/4))` if we're *not* tracking the exchange
///   kink, or the off-diagonal `½(exp(ΔτJ/4) − exp(−3ΔτJ/4))` if a kink
///   (spin exchange) is present at this bond.
///
/// For the **path-integral representation** the temporal bond between slices
/// carries a kink whenever the spins differ across slices — that's the
/// off-diagonal matrix element. The single-spin-flip Metropolis kernel flips
/// spins and recomputes matrix elements; both branches below give the correct
/// weight for the corresponding (diagonal / kink) sector.
#[derive(Debug, Clone)]
pub struct HeisenbergChain {
    /// Coupling J. J > 0 = ferromagnetic, J < 0 = antiferromagnetic.
    pub j: f64,
}

impl HeisenbergChain {
    pub fn new(j: f64) -> Self {
        Self { j }
    }

    /// Eigenvalue-factorized matrix elements, to avoid cancellation.
    /// `diag_ap(dtau) = ½(exp(ΔτJ/4) + exp(−3ΔτJ/4))`.
    #[inline]
    fn diag_antiparallel(&self, dtau: f64) -> f64 {
        0.5 * ((dtau * self.j * 0.25).exp() + (-3.0 * dtau * self.j * 0.25).exp())
    }

    /// `offdiag_ap(dtau) = ½(exp(ΔτJ/4) − exp(−3ΔτJ/4))`.
    #[inline]
    fn offdiag_antiparallel(&self, dtau: f64) -> f64 {
        0.5 * ((dtau * self.j * 0.25).exp() - (-3.0 * dtau * self.j * 0.25).exp())
    }
}

impl QuantumHamiltonian for HeisenbergChain {
    #[inline]
    fn bond_diagonal(&self, s_i: u8, s_j: u8, dtau: f64) -> f64 {
        if s_i == s_j {
            // Parallel ↑↑ or ↓↓: triplet, eigenvalue +1/4.
            (dtau * self.j * 0.25).exp()
        } else {
            // Antiparallel, no kink: diagonal of the singlet/triplet block.
            self.diag_antiparallel(dtau)
        }
    }

    #[inline]
    fn bond_offdiag(&self, s_i: u8, s_j: u8, dtau: f64) -> f64 {
        if s_i == s_j {
            // ⟨↑↑|B|↓↓⟩ = 0 — different total Sz, no coupling.
            0.0
        } else {
            // Antiparallel, kink present (spin exchange): off-diagonal.
            self.offdiag_antiparallel(dtau)
        }
    }

    #[inline]
    fn energy_per_bond(&self, s_i: u8, s_j: u8) -> f64 {
        // ⟨sᵢsⱼ| H_bond |sᵢsⱼ⟩ = -J ⟨Sᵢ·Sⱼ⟩. Diagonal matrix element of Sᵢ·Sⱼ:
        // parallel ↑↑ → +1/4, antiparallel ↑↓ → -1/4.
        // (Note: this is the *diagonal* energy estimator; for a closed PI
        // config the actual ⟨H⟩ includes both diagonal and off-diagonal
        // contributions, but the diagonal-only estimator is what the
        // standard Metropolis samples in the spin basis.)
        let sdot = if s_i == s_j { 0.25 } else { -0.25 };
        -self.j * sdot
    }

    /// Full quantum PI estimator: `−(d/dΔτ) ln B` evaluated in the sector
    /// the bond actually carries.
    ///
    /// For the Heisenberg chain this is what reproduces the Bethe-ansatz
    /// ground-state energy: antiparallel bonds *with kinks* (the
    /// spin-exchange sector) contribute below the classical −J/4.
    #[inline]
    fn bond_energy_estimator(&self, s_i: u8, s_j: u8, dtau: f64, is_kink: bool) -> f64 {
        if s_i == s_j {
            // Parallel: triplet, B = exp(ΔτJ/4). −d/dΔτ = -J/4.
            -0.25 * self.j
        } else {
            // Antiparallel. e₊ = exp(ΔτJ/4), e₋ = exp(-3ΔτJ/4).
            let ep = (dtau * self.j * 0.25).exp();
            let em = (-3.0 * dtau * self.j * 0.25).exp();
            if is_kink {
                // Off-diagonal: B = ½(e₊ − e₋). −d/dΔτ ln B = -J(e₊ + 3e₋)/(e₊ − e₋).
                -self.j * (ep + 3.0 * em) / (ep - em)
            } else {
                // Diagonal: B = ½(e₊ + e₋). −d/dΔτ ln B = -J(e₊ − 3e₋)/(e₊ + e₋).
                -self.j * (ep - 3.0 * em) / (ep + em)
            }
        }
    }

    #[inline]
    fn coupling(&self) -> f64 {
        self.j
    }
}

/// Exact ground-state energy per site of the S = 1/2 Heisenberg chain
/// (Bethe ansatz): `e₀ = J(1/4 − ln 2)`. Antiferromagnetic (J > 0) → negative.
///
/// Used by the validation test as the acceptance target.
pub fn heisenberg_chain_ground_energy_per_site(j: f64) -> f64 {
    j * (0.25 - consts::LN_2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_parallel_value() {
        let h = HeisenbergChain::new(1.0);
        let dtau = 0.1;
        // ↑↑ : exp(Δτ J/4)
        let para = h.bond_diagonal(1, 1, dtau);
        assert!((para - (dtau * 0.25).exp()).abs() < 1e-12);
    }

    #[test]
    fn diagonal_antiparallel_value() {
        let h = HeisenbergChain::new(1.0);
        let dtau = 0.1;
        // ↑↓ : ½(exp(ΔτJ/4) + exp(-3ΔτJ/4))
        let anti = h.bond_diagonal(0, 1, dtau);
        let expected = 0.5 * ((dtau * 0.25).exp() + (-3.0 * dtau * 0.25).exp());
        assert!((anti - expected).abs() < 1e-12);
    }

    #[test]
    fn offdiag_only_for_antiparallel() {
        let h = HeisenbergChain::new(1.0);
        let dtau = 0.1;
        assert_eq!(h.bond_offdiag(1, 1, dtau), 0.0);
        assert_eq!(h.bond_offdiag(0, 0, dtau), 0.0);
        // ↑↓ : ½(exp(ΔτJ/4) − exp(-3ΔτJ/4))
        let off = h.bond_offdiag(1, 0, dtau);
        let expected = 0.5 * ((dtau * 0.25).exp() - (-3.0 * dtau * 0.25).exp());
        assert!((off - expected).abs() < 1e-12);
    }

    #[test]
    fn energy_per_bond_signs() {
        let h = HeisenbergChain::new(1.0);
        // Parallel: S_i·S_j diagonal = +1/4, H_bond = -J(1/4) = -0.25
        assert!((h.energy_per_bond(1, 1) - (-0.25)).abs() < 1e-12);
        // Antiparallel: diagonal S_i·S_j = -1/4, H_bond = -J(-1/4) = +0.25
        assert!((h.energy_per_bond(0, 1) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn bethe_ansatz_value() {
        // e₀/J = 1/4 - ln(2) ≈ -0.4431472
        let e = heisenberg_chain_ground_energy_per_site(1.0);
        assert!((e - (0.25 - std::f64::consts::LN_2)).abs() < 1e-12);
        assert!(e < 0.0);
        assert!((e - (-0.4431)).abs() < 1e-3);
    }
}
