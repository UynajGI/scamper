//! Local-energy estimator.
//!
//! `E_L(R) = ψ_T⁻¹ H ψ_T (R) = −½ Σ_i (∇_i² ln|ψ| + |∇_i ln|ψ||²) + V(R)`
//! (units ħ = m = 1). All derivative physics comes from the
//! [`WaveFunction`] implementation; all potential physics from the
//! [`ContinuumHamiltonian`]; this module only combines them. The gradient
//! buffer is caller-owned so repeated measurements never allocate.

use super::super::hamiltonian::ContinuumHamiltonian;
use super::super::wavefunction::{GradBuffer, WaveFunction};

/// One local-energy sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalEnergy {
    /// `E_L(R)`.
    pub value: f64,
    /// `Σ_i |∇_i ln|ψ|(R)|²` — the drift-velocity norm squared, measured
    /// for DMC-facing diagnostics and reusable by variance analysis.
    pub log_grad_squared: f64,
}

/// Evaluate the local energy of `cfg` under `wave_function` and
/// `hamiltonian`, using `grad` as reusable `∇ ln|ψ|` scratch space.
pub fn local_energy<W: WaveFunction>(
    wave_function: &W,
    hamiltonian: &ContinuumHamiltonian,
    cfg: &W::Config,
    grad: &mut GradBuffer,
) -> LocalEnergy {
    // log_grad accumulates (composite ansatz support); start from zero.
    grad.clear();
    wave_function.log_grad(cfg, grad);
    let log_grad_squared = grad.as_slice().iter().map(|&x| x * x).sum();
    let kinetic = -0.5 * (wave_function.log_laplacian(cfg) + log_grad_squared);
    LocalEnergy {
        value: kinetic + hamiltonian.potential_energy(cfg),
        log_grad_squared,
    }
}
