//! Slater determinants of cartesian Gaussian orbitals (L1).
//!
//! `psi_T = det^up * det^down` over contracted cartesian GTOs
//! `phi_m(x) = sum_k c_k x^{p_k} y^{q_k} z^{s_k} exp(-alpha |x|^2)`,
//! optionally evaluated at backflow quasiparticle coordinates
//! `x_i = r_i + sum_j eta(r_ij)(r_i - r_j)` (Kwon–Ceperley–Martin,
//! [`Backflow`](super::backflow::Backflow)). The exponent/coefficients of
//! every orbital and the backflow scale `lambda` are variational.
//!
//! # Conventions
//!
//! Spin blocks: particles `[0, n_up)` are spin-up (columns of the up
//! determinant), particles `[n_up, n_up + n_down)` spin-down. Up and down
//! blocks carry independent orbital lists; a block may be empty (a
//! fully-polarized single-species determinant), but not both. Spin
//! degeneracy is **not** required — any `(n_up, n_down)` pair with
//! matching orbital counts is accepted (open shells lose the closed-shell
//! zero-variance guarantee but are physically legitimate).
//!
//! Determinant matrix per block: `D[(m, i)] = phi_m(x_i)` — rows orbitals,
//! columns particles (so rows of `D^{-1}` are particle indices).
//! `ln|psi| = sum_blocks ln|det D|` via nalgebra LU (mature-crate policy:
//! no hand-rolled factorizations). Gradients and Laplacians are
//! hand-derived chains through `(D^{-1}, grad phi, Hess phi)` and — when
//! backflow is active — the displacement Jacobians:
//!
//! ```text
//! grad_k ln|det| = sum_{i in block} (A^{ik})^T O_i,
//!     O_i    = sum_m (D^{-1})_{im} grad phi_m(x_i)
//! lap_k ln|det|  = sum_{i in block} [ Tr(Q_i A^{ik} (A^{ik})^T) + O_i . lap_k x_i ]
//!                  - sum_{i,l in block} (A^{ik} w^{il}) . (A^{lk} w^{li}),
//!     Q_i    = sum_m (D^{-1})_{im} Hess phi_m(x_i),
//!     w^{il} = sum_m (D^{-1})_{lm} grad phi_m(x_i),   A^{ik} = dx_i / dr_k
//! ```
//!
//! Derivation: the plain-Slater case is the standard textbook chain
//! (Becca–Sorella 2017 ch. 3.3); the backflow extension applies the
//! second-derivative chain rule to `ln det D(x(r))` — the Hessian of
//! `ln det` w.r.t. the matrix entries is
//! `∂²ln det/∂D[m,i]∂D[m',l] = −(D^{-1})_{im'}(D^{-1})_{lm}`, which after
//! chaining through the orbitals leaves the rank-structured cross term
//! over `w^{il} ⊗ w^{li}`. (Kwon–Ceperley–Martin 1998 §II states the
//! transformation; Holzmann et al. 2019 note that the derivatives follow
//! by "applying the chain rule in the usual way" — neither prints the
//! expanded equations, so they are derived here and enforced against
//! central finite differences in CI.) One load-bearing fact: the fully
//! general chain rule contracts `w^{il}` against
//! `A^{ik} (A^{lk})^T`, while the form above contracts against
//! `(A^{ik})^T A^{lk}` — identical here because every `A^{ik}` is
//! **symmetric** (the pairwise displacement Jacobian
//! `J = eta' t t^T / r + eta I` is symmetric). A future backflow with
//! non-symmetric Jacobians must revisit the `log_grad`/`log_laplacian`
//! contractions (the finite-difference gate would catch the slip).
//! Without backflow,
//! `A^{ik} = δ_{ik} I` collapses the cross term to `−|O_k|²`, recovering
//! the elementary identity
//! `∇² ln|ψ| = ψ^{-1}∇²ψ − |∇ ln|ψ||²` — which is why the estimator's
//! `−½(∇² ln|ψ| + |∇ ln|ψ||²)` reproduces the standard
//! `−½ ψ^{-1}∇²ψ` determinant kinetic energy.
//!
//! # Single-particle fast path (plain Slater)
//!
//! Moving particle `i` changes only column `i` of its own block's `D`, so
//! the Metropolis ratio is the Sherman–Morrison column identity
//!
//! ```text
//! det(D') / det(D) = (D^{-1} d)_i,     d_m = phi_m(r_new),
//! ```
//!
//! an O(N) dot product against the cached inverse, with the O(N²) inverse
//! row update `D'^{-1} = D^{-1} - v (e_i^T D^{-1}) / (D^{-1} d)_i`,
//! `v = D^{-1} d - e_i`, applied on accept. Floating-point drift of the
//! incrementally updated inverse is repaired by periodic full rebuilds
//! (the `rebuild` hook; [`VmcKernel`](crate::variational::VmcKernel)
//! rebuilds and re-anchors every walker at its sweep entry — the
//! K-rebuild policy). With **backflow active** every
//! quasiparticle coordinate changes when one particle moves, so no rank-1
//! identity survives and `delta_log` is a full recompute — the standard
//! backflow cost noted in Kwon et al. 1998 §II.
//!
//! # Allocation contract
//!
//! The plain-Slater proposal path (`delta_log` against a current cache) is
//! allocation-free. Accepted plain-Slater moves (`commit_move`) apply the
//! O(N²) Sherman–Morrison update through nalgebra and allocate its O(N²)
//! result matrix; full determinant evaluations (rebuilds, backflow moves,
//! and `log_grad`/`log_laplacian`/`log_grad_params`, which
//! recompute the LU from scratch for drift-free measurements) go through
//! nalgebra and allocate a bounded number of N×N matrices per call — a
//! documented deviation from the stateless-L0 zero-alloc wording,
//! amortized over O(N³) work and off the proposal path.
//!
//! # Fermion signs
//!
//! VMC samples `|psi|^2`; the determinant sign is dropped everywhere
//! (`ln|det|`), which is exact for observables quadratic in `psi`. The
//! sign problem / fixed-node constraint is L3 business.
//!
//! # Configuration-length mismatches
//!
//! The ansatz fixes `n_particles = n_up + n_down`. Evaluation methods
//! handed a configuration of any other length return `NaN` (`log_psi`,
//! `log_laplacian`, `delta_log`) and contribute nothing (`log_grad`)
//! rather than panicking — the `VmcKernel` constructor rejects such
//! populations up front via its non-finite `log|psi|` check (criterion G).

use std::cell::RefCell;

use nalgebra::DMatrix;

use super::super::error::VariationalError;
use super::backflow::Backflow;
use super::{
    read_particle, write_particle, DeltaLog, GradBuffer, ParamGradBuffer, Point, Positions,
    WaveFunction, WaveFunctionParams, DIM,
};

/// Highest monomial power per Cartesian direction accepted in a GTO
/// primitive (the L1 orbital catalog — s/p/d shells — needs at most 2).
const MAX_POWER: u8 = 4;

/// Integer power for tiny exponents (avoids `powi` overhead and keeps
/// `0^0 = 1` well-defined).
#[inline]
fn ipow(x: f64, n: u8) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => x * x * x,
        4 => {
            let s = x * x;
            s * s
        }
        _ => x.powi(n as i32),
    }
}

/// Value, gradient and Hessian of an orbital at one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitalValue {
    /// `phi(t)`.
    pub value: f64,
    /// `grad phi(t)`.
    pub grad: [f64; DIM],
    /// `Hess phi(t)` (symmetric).
    pub hessian: [[f64; DIM]; DIM],
}

/// Value/gradient/Hessian of one monomial `t_0^p t_1^q t_2^s` — elementary
/// power-rule derivatives (no divisions: terms with a negative effective
/// power are zero, so the origin is safe).
fn monomial_derivs(powers: &[u8; DIM], t: &Point) -> OrbitalValue {
    let mut powv = [1.0_f64; DIM];
    let mut value = 1.0;
    for k in 0..DIM {
        powv[k] = ipow(t[k], powers[k]);
        value *= powv[k];
    }
    // Product over directions excluding `skip` (and `skip2`).
    let rest = |skip: usize, skip2: usize| {
        let mut out = 1.0;
        for (l, &value) in powv.iter().enumerate() {
            if l != skip && l != skip2 {
                out *= value;
            }
        }
        out
    };
    let mut grad = [0.0; DIM];
    let mut hessian = [[0.0; DIM]; DIM];
    for k in 0..DIM {
        if powers[k] > 0 {
            grad[k] = powers[k] as f64 * ipow(t[k], powers[k] - 1) * rest(k, k);
        }
        for l in 0..DIM {
            if k == l {
                if powers[k] >= 2 {
                    hessian[k][k] = powers[k] as f64
                        * (powers[k] - 1) as f64
                        * ipow(t[k], powers[k] - 2)
                        * rest(k, k);
                }
            } else if powers[k] > 0 && powers[l] > 0 {
                hessian[k][l] = powers[k] as f64
                    * powers[l] as f64
                    * ipow(t[k], powers[k] - 1)
                    * ipow(t[l], powers[l] - 1)
                    * rest(k, l);
            }
        }
    }
    OrbitalValue {
        value,
        grad,
        hessian,
    }
}

#[inline]
fn gaussian(exponent: f64, t: &Point) -> f64 {
    let r2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    (-exponent * r2).exp()
}

/// A contracted cartesian Gaussian orbital
/// `phi(t) = e^{-alpha |t|^2} sum_k c_k t_0^{p_k} t_1^{q_k} t_2^{s_k}`.
///
/// Analytic derivatives (hand-derived, elementary product rule): with
/// `S(t) = sum_k c_k m_k(t)`, `G(t) = e^{-alpha |t|^2}` and
/// `grad G = -2 alpha t G`, `Hess G = (4 alpha^2 t t^T - 2 alpha I) G`:
///
/// ```text
/// phi         = G S
/// grad_i phi  = (grad G)_i S + G S_i
/// Hess phi    = (Hess G) S + (grad G)(grad S)^T + (grad S)(grad G)^T + G Hess S
/// ```
///
/// Normalization constants are deliberately omitted: a per-orbital scale
/// multiplies the determinant by a constant, which cancels in every
/// physical quantity (`E_L`, `|psi|^2` sampling, move ratios).
#[derive(Debug, Clone, PartialEq)]
pub struct GtoOrbital {
    exponent: f64,
    primitives: Vec<(f64, [u8; DIM])>,
}

impl GtoOrbital {
    /// Construct from a Gaussian exponent (finite, > 0) and a non-empty
    /// primitive list of (coefficient, monomial powers); coefficients must
    /// be finite and powers at most [`MAX_POWER`].
    pub fn new(exponent: f64, primitives: Vec<(f64, [u8; DIM])>) -> Result<Self, VariationalError> {
        VariationalError::require_positive("exponent", exponent)?;
        if primitives.is_empty() {
            return Err(VariationalError::invalid(
                "primitives",
                "an orbital needs at least one cartesian primitive",
            ));
        }
        for (coefficient, powers) in &primitives {
            if !coefficient.is_finite() {
                return Err(VariationalError::invalid(
                    "coefficient",
                    format!("must be finite, got {coefficient}"),
                ));
            }
            if powers.iter().any(|&p| p > MAX_POWER) {
                return Err(VariationalError::invalid(
                    "powers",
                    format!("monomial powers must be <= {MAX_POWER}"),
                ));
            }
        }
        Ok(Self {
            exponent,
            primitives,
        })
    }

    /// The Gaussian exponent.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The primitive list (coefficient, monomial powers) in construction
    /// order — the parameter-vector order for coefficients.
    pub fn primitives(&self) -> &[(f64, [u8; DIM])] {
        &self.primitives
    }

    /// `phi(t)`.
    pub fn value(&self, t: &Point) -> f64 {
        let mut s = 0.0;
        for &(coefficient, powers) in &self.primitives {
            s += coefficient * monomial_derivs(&powers, t).value;
        }
        s * gaussian(self.exponent, t)
    }

    /// Value, gradient and Hessian in one pass.
    pub fn evaluate(&self, t: &Point) -> OrbitalValue {
        let mut s = 0.0;
        let mut s_grad = [0.0; DIM];
        let mut s_hess = [[0.0; DIM]; DIM];
        for &(coefficient, powers) in &self.primitives {
            let monomial = monomial_derivs(&powers, t);
            s += coefficient * monomial.value;
            for k in 0..DIM {
                s_grad[k] += coefficient * monomial.grad[k];
                for (l, hess_kl) in s_hess[k].iter_mut().enumerate() {
                    *hess_kl += coefficient * monomial.hessian[k][l];
                }
            }
        }
        let g = gaussian(self.exponent, t);
        // grad G = -2 alpha t G; Hess G = (4 alpha^2 t t^T - 2 alpha I) G.
        let mut g_grad = [0.0; DIM];
        let mut g_hess = [[0.0; DIM]; DIM];
        for k in 0..DIM {
            g_grad[k] = -2.0 * self.exponent * t[k] * g;
            for l in 0..DIM {
                g_hess[k][l] = 4.0 * self.exponent * self.exponent * t[k] * t[l] * g;
            }
            g_hess[k][k] -= 2.0 * self.exponent * g;
        }
        let mut grad = [0.0; DIM];
        let mut hessian = [[0.0; DIM]; DIM];
        for k in 0..DIM {
            grad[k] = g_grad[k] * s + g * s_grad[k];
            for l in 0..DIM {
                hessian[k][l] = g_hess[k][l] * s
                    + g_grad[k] * s_grad[l]
                    + g_grad[l] * s_grad[k]
                    + g * s_hess[k][l];
            }
        }
        OrbitalValue {
            value: g * s,
            grad,
            hessian,
        }
    }

    /// Primitive values `g_k(t) = m_k(t) e^{-alpha |t|^2}` — the
    /// `partial phi / partial c_k` of the parameter-gradient chain
    /// (optimizer-facing path, allocation allowed).
    pub fn primitive_values(&self, t: &Point) -> Vec<f64> {
        let g = gaussian(self.exponent, t);
        self.primitives
            .iter()
            .map(|&(_, powers)| monomial_derivs(&powers, t).value * g)
            .collect()
    }
}

/// Physicists' Hermite polynomial `H_n(sqrt(omega) t)` as monomial
/// coefficients over `t` (index = power): `H_0 = 1`,
/// `H_1 = 2 sqrt(w) t`, `H_{n+1} = 2 sqrt(w) t H_n - 2 n H_{n-1}`.
fn hermite_coeffs(n: usize, omega: f64) -> Vec<f64> {
    let u = omega.sqrt();
    if n == 0 {
        return vec![1.0];
    }
    let mut prev = vec![1.0_f64];
    let mut current = vec![0.0, 2.0 * u];
    for k in 1..n {
        let mut next = vec![0.0; current.len() + 1];
        for (i, &c) in current.iter().enumerate() {
            next[i + 1] += 2.0 * u * c;
        }
        for (i, &c) in prev.iter().enumerate() {
            next[i] -= 2.0 * k as f64 * c;
        }
        prev = current;
        current = next;
    }
    current
}

/// Exact harmonic-oscillator eigen-orbitals (per spin) filling the closed
/// shells `0..n_shells` — the zero-variance reference of the L1 tests.
///
/// The 3-D HO eigenfunctions are products of 1-D Hermite functions:
/// `phi_{n_x n_y n_z}(t) = H_{n_x}(sqrt(w) t_x) H_{n_y}(sqrt(w) t_y)
/// H_{n_z}(sqrt(w) t_z) e^{-w |t|^2 / 2}` with energy
/// `(n_x + n_y + n_z + 3/2) w`. With the GTO exponent `alpha = w/2` these
/// are contracted cartesian Gaussians; every polynomial coefficient is
/// kept exactly (the `4 w t^2 - 2`-type lower-power parts are what make a
/// filled shell the exact eigenstate — raw `t^2`-type orbitals carry
/// shell-0 admixtures and would *not* give a constant local energy).
/// Per-orbital normalization is omitted (it cancels; see [`GtoOrbital`]).
/// Shell degeneracies: 1, 3, 6 orbitals (shells 0, 1, 2) — `n_shells` of
/// 1, 2, 3 gives 1, 4, 10 orbitals per spin, i.e. 2, 8, 20 electrons for
/// the two-spin closed shell.
pub fn harmonic_trap_orbitals(
    omega: f64,
    n_shells: usize,
) -> Result<Vec<GtoOrbital>, VariationalError> {
    VariationalError::require_positive("omega", omega)?;
    if n_shells == 0 || n_shells > 3 {
        return Err(VariationalError::invalid(
            "n_shells",
            "L1 provides the exact HO shells 0..=2 (2, 8 or 20 electrons)",
        ));
    }
    let alpha = omega / 2.0;
    let mut orbitals = Vec::new();
    for shell in 0..n_shells {
        for nx in 0..=shell {
            for ny in 0..=(shell - nx) {
                let nz = shell - nx - ny;
                let px = hermite_coeffs(nx, omega);
                let py = hermite_coeffs(ny, omega);
                let pz = hermite_coeffs(nz, omega);
                let mut primitives = Vec::new();
                for (a, &cx) in px.iter().enumerate() {
                    for (b, &cy) in py.iter().enumerate() {
                        for (c, &cz) in pz.iter().enumerate() {
                            let coefficient = cx * cy * cz;
                            if coefficient != 0.0 {
                                primitives.push((coefficient, [a as u8, b as u8, c as u8]));
                            }
                        }
                    }
                }
                debug_assert!(!primitives.is_empty());
                orbitals.push(GtoOrbital::new(alpha, primitives)?);
            }
        }
    }
    Ok(orbitals)
}

/// Electron count of the two-spin closed shell `n_shells`:
/// `N = sum_{n=0}^{K} 2 [(n+1)(n+2)/2] = sum_{n=0}^{K} (n+1)(n+2)
///    = (K+1)(K+2)(K+3)/3` with `K = n_shells - 1`
/// (2, 8, 20 for `n_shells` 1, 2, 3).
pub fn harmonic_closed_shell_electrons(n_shells: usize) -> Result<usize, VariationalError> {
    if n_shells == 0 || n_shells > 3 {
        return Err(VariationalError::invalid(
            "n_shells",
            "L1 provides the exact HO shells 0..=2 (2, 8 or 20 electrons)",
        ));
    }
    Ok(n_shells * (n_shells + 1) * (n_shells + 2) / 3)
}

/// Exact non-interacting energy of the two-spin closed shell `n_shells`:
///
/// `E_0 = w sum_{n=0}^{K} 2 [(n+1)(n+2)/2] (n + 3/2)
///      = (w/2) sum_{n=0}^{K} (n+1)(n+2)(2n+3)`
///
/// — every one of the `2 (n+1)(n+2)/2` fermions of shell `n` carries
/// `(n + 3/2) w`. Values: `3w`, `18w`, `60w` for 2, 8, 20 electrons.
pub fn harmonic_closed_shell_energy(omega: f64, n_shells: usize) -> Result<f64, VariationalError> {
    VariationalError::require_positive("omega", omega)?;
    if n_shells == 0 || n_shells > 3 {
        return Err(VariationalError::invalid(
            "n_shells",
            "L1 provides the exact HO shells 0..=2 (2, 8 or 20 electrons)",
        ));
    }
    let mut total = 0.0;
    for n in 0..n_shells {
        let shell = n as f64;
        total += (shell + 1.0) * (shell + 2.0) * (2.0 * shell + 3.0);
    }
    Ok(0.5 * omega * total)
}

/// Per-block cached state (the L1 determinant machinery's incremental
/// inverse) plus preallocated proposal scratch. Interior mutability: the
/// `WaveFunction` evaluation methods take `&self` (kernel contract), so
/// the lazily rebuilt cache lives behind a `RefCell` — single-threaded
/// kernel use, never shared across the sweep loop.
#[derive(Debug, Clone)]
struct DetCache {
    /// Configuration coordinates the cached inverse/log-dets belong to.
    fingerprint: Vec<f64>,
    up_to_date: bool,
    log_det: [f64; 2],
    inverse: [DMatrix<f64>; 2],
    /// Reused proposal column (`phi_m(r_new)`).
    column: Vec<f64>,
    /// Reused moved-configuration buffer (backflow `delta_log` path).
    moved: Vec<f64>,
}

impl DetCache {
    fn matches(&self, coords: &[f64]) -> bool {
        self.up_to_date && self.fingerprint == coords
    }
}

/// A two-spin-block Slater determinant of contracted cartesian GTOs with
/// optional backflow (see the module docs for conventions and chains).
#[derive(Debug, Clone)]
pub struct SlaterDeterminant {
    orbitals: [Vec<GtoOrbital>; 2],
    n_up: usize,
    n_down: usize,
    backflow: Option<Backflow>,
    cache: RefCell<DetCache>,
}

impl SlaterDeterminant {
    /// Plain Slater determinant (no backflow) from validated orbital
    /// lists; at least one block must be non-empty.
    pub fn new(
        orbitals_up: Vec<GtoOrbital>,
        orbitals_down: Vec<GtoOrbital>,
    ) -> Result<Self, VariationalError> {
        Self::assemble(orbitals_up, orbitals_down, None)
    }

    /// Slater determinant of backflow quasiparticle coordinates.
    pub fn with_backflow(
        orbitals_up: Vec<GtoOrbital>,
        orbitals_down: Vec<GtoOrbital>,
        backflow: Backflow,
    ) -> Result<Self, VariationalError> {
        Self::assemble(orbitals_up, orbitals_down, Some(backflow))
    }

    /// Closed-shell harmonic-trap determinant (exact HO orbitals, both
    /// spins) — the zero-variance L1 reference.
    pub fn harmonic_trap(omega: f64, n_shells: usize) -> Result<Self, VariationalError> {
        let orbitals = harmonic_trap_orbitals(omega, n_shells)?;
        Self::new(orbitals.clone(), orbitals)
    }

    /// Closed-shell harmonic-trap determinant with backflow.
    pub fn harmonic_trap_with_backflow(
        omega: f64,
        n_shells: usize,
        backflow: Backflow,
    ) -> Result<Self, VariationalError> {
        let orbitals = harmonic_trap_orbitals(omega, n_shells)?;
        Self::with_backflow(orbitals.clone(), orbitals, backflow)
    }

    fn assemble(
        orbitals_up: Vec<GtoOrbital>,
        orbitals_down: Vec<GtoOrbital>,
        backflow: Option<Backflow>,
    ) -> Result<Self, VariationalError> {
        if orbitals_up.is_empty() && orbitals_down.is_empty() {
            return Err(VariationalError::invalid(
                "orbitals",
                "at least one spin block must carry an orbital",
            ));
        }
        // Orbital lists are built through GtoOrbital::new (validated);
        // re-assert cheaply for defense in depth.
        for (block, orbitals) in [&orbitals_up, &orbitals_down].iter().enumerate() {
            for (m, orbital) in orbitals.iter().enumerate() {
                if !(orbital.exponent().is_finite() && orbital.exponent() > 0.0)
                    || orbital.primitives().is_empty()
                {
                    return Err(VariationalError::invalid(
                        "orbitals",
                        format!("block {block} orbital {m} is not a valid GTO"),
                    ));
                }
            }
        }
        let n_up = orbitals_up.len();
        let n_down = orbitals_down.len();
        let max_block = n_up.max(n_down);
        Ok(Self {
            orbitals: [orbitals_up, orbitals_down],
            n_up,
            n_down,
            backflow,
            cache: RefCell::new(DetCache {
                fingerprint: Vec::new(),
                up_to_date: false,
                log_det: [f64::NAN, f64::NAN],
                inverse: [DMatrix::zeros(0, 0), DMatrix::zeros(0, 0)],
                column: vec![0.0; max_block],
                moved: vec![0.0; DIM * (n_up + n_down)],
            }),
        })
    }

    /// Number of spin-up (block-0) particles.
    #[inline]
    pub const fn n_up(&self) -> usize {
        self.n_up
    }

    /// Number of spin-down (block-1) particles.
    #[inline]
    pub const fn n_down(&self) -> usize {
        self.n_down
    }

    /// The configuration length this ansatz evaluates (`n_up + n_down`).
    #[inline]
    pub const fn expected_particles(&self) -> usize {
        self.n_up + self.n_down
    }

    /// The backflow transformation, if constructed with one.
    #[inline]
    pub fn backflow(&self) -> Option<&Backflow> {
        self.backflow.as_ref()
    }

    /// The backflow, but only when it actually displaces coordinates
    /// (`lambda != 0`). At `lambda == 0` every chain routes through the
    /// plain-Slater code path, which makes the backflow-off reduction
    /// **bit-exact** rather than merely exact in exact arithmetic.
    #[inline]
    fn active_backflow(&self) -> Option<&Backflow> {
        match &self.backflow {
            Some(backflow) if backflow.lambda() != 0.0 => Some(backflow),
            _ => None,
        }
    }

    /// `(offset, len)` of spin block `block`.
    #[inline]
    const fn block_bounds(&self, block: usize) -> (usize, usize) {
        if block == 0 {
            (0, self.n_up)
        } else {
            (self.n_up, self.n_down)
        }
    }

    /// `(ln|det D|, D^{-1})` for one spin block at `coords` (pure
    /// evaluation, no cache); `None` marks a singular or non-finite
    /// determinant.
    fn block_lu(
        &self,
        coords: &[f64],
        block: usize,
        backflow: Option<&Backflow>,
    ) -> Option<(f64, DMatrix<f64>)> {
        let matrix = self.block_matrix(coords, block, backflow);
        lu_logdet_inverse(&matrix)
    }

    /// The determinant matrix `D[(m, i)] = phi_m(x_i)` of one block
    /// (allocating; full-evaluation paths only).
    fn block_matrix(
        &self,
        coords: &[f64],
        block: usize,
        backflow: Option<&Backflow>,
    ) -> DMatrix<f64> {
        let (offset, n_blk) = self.block_bounds(block);
        let orbitals = &self.orbitals[block];
        let mut matrix = DMatrix::zeros(n_blk, n_blk);
        for i in 0..n_blk {
            let x = quasiparticle(coords, offset + i, backflow);
            for (m, orbital) in orbitals.iter().enumerate() {
                matrix[(m, i)] = orbital.value(&x);
            }
        }
        matrix
    }

    /// Full recompute of the cached inverse/log-dets for `coords`.
    fn rebuild_coords(&self, coords: &[f64]) {
        let backflow = self.active_backflow();
        let mut cache = self.cache.borrow_mut();
        // clear + extend (not copy_from_slice): the first rebuild starts
        // from an empty fingerprint buffer.
        cache.fingerprint.clear();
        cache.fingerprint.extend_from_slice(coords);
        cache.up_to_date = true;
        for block in 0..2 {
            match self.block_lu(coords, block, backflow) {
                Some((log_det, inverse)) => {
                    cache.log_det[block] = log_det;
                    cache.inverse[block] = inverse;
                }
                None => {
                    let (_, n_blk) = self.block_bounds(block);
                    cache.log_det[block] = f64::NAN;
                    cache.inverse[block] = DMatrix::repeat(n_blk, n_blk, f64::NAN);
                }
            }
        }
    }

    /// Make the cache current for `coords` (O(N³) when stale — e.g. the
    /// walker switch inside a kernel sweep; a no-op when current).
    fn ensure_cache(&self, coords: &[f64]) {
        if !self.cache.borrow().matches(coords) {
            self.rebuild_coords(coords);
        }
    }

    /// Total `ln|psi|` at `coords` (pure, both blocks; the backflow
    /// `delta_log` full recompute).
    fn total_logdet(&self, coords: &[f64], backflow: Option<&Backflow>) -> f64 {
        let mut total = 0.0;
        for block in 0..2 {
            total += match self.block_lu(coords, block, backflow) {
                Some((log_det, _)) => log_det,
                None => f64::NAN,
            };
        }
        total
    }
}

/// Quasiparticle coordinate of particle `index` (backflow displacement
/// added when active).
#[inline]
fn quasiparticle(coords: &[f64], index: usize, backflow: Option<&Backflow>) -> Point {
    let mut x = read_particle(coords, index);
    if let Some(backflow) = backflow {
        let displacement = backflow.displacement(coords, index);
        for k in 0..DIM {
            x[k] += displacement[k];
        }
    }
    x
}

/// `(ln|det M|, M^{-1})` via nalgebra LU with partial pivoting (the
/// mature-crate policy: never hand-roll factorizations). `None` on
/// singular or non-finite input.
fn lu_logdet_inverse(matrix: &DMatrix<f64>) -> Option<(f64, DMatrix<f64>)> {
    let n = matrix.nrows();
    if n == 0 {
        return Some((0.0, DMatrix::zeros(0, 0)));
    }
    let lu = matrix.clone().lu();
    let det = lu.determinant();
    let inverse = lu.solve(&DMatrix::identity(n, n))?;
    if !det.is_finite() || det == 0.0 || !inverse.iter().all(|x| x.is_finite()) {
        return None;
    }
    Some((det.abs().ln(), inverse))
}

impl WaveFunction for SlaterDeterminant {
    type Config = Positions;

    fn log_psi(&self, cfg: &Self::Config) -> f64 {
        let coords = cfg.as_ref();
        if coords.len() != DIM * self.expected_particles() {
            return f64::NAN;
        }
        self.ensure_cache(coords);
        let cache = self.cache.borrow();
        cache.log_det[0] + cache.log_det[1]
    }

    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer) {
        let coords = cfg.as_ref();
        if coords.len() != DIM * self.expected_particles() {
            return; // documented NaN policy; the kernel gates the length
        }
        let n = self.expected_particles();
        let backflow = self.active_backflow();
        let out = out.as_mut_slice();
        for block in 0..2 {
            let (offset, n_blk) = self.block_bounds(block);
            if n_blk == 0 {
                continue;
            }
            let chain = BlockChain::build(self, coords, block, backflow);
            // Identity-Jacobian part: grad_k = O_k for own-block particles.
            for i in 0..n_blk {
                let base = DIM * (offset + i);
                for k in 0..DIM {
                    out[base + k] += chain.o[i][k];
                }
            }
            // Backflow part: grad_k = sum_{i in block} O_i . (dx_i / dr_k).
            if let Some(backflow) = backflow {
                for k in 0..n {
                    let base = DIM * k;
                    if k >= offset && k < offset + n_blk {
                        let o_k = chain.o[k - offset];
                        for j in 0..n {
                            if j == k {
                                continue;
                            }
                            add_jacobian_action(
                                &mut out[base..base + DIM],
                                &backflow.pair_jacobian(coords, k, j),
                                &o_k,
                                1.0,
                            );
                        }
                    }
                    for i in 0..n_blk {
                        let global = offset + i;
                        if global == k {
                            continue;
                        }
                        add_jacobian_action(
                            &mut out[base..base + DIM],
                            &backflow.pair_jacobian(coords, global, k),
                            &chain.o[i],
                            -1.0,
                        );
                    }
                }
            }
        }
    }

    fn log_laplacian(&self, cfg: &Self::Config) -> f64 {
        let coords = cfg.as_ref();
        if coords.len() != DIM * self.expected_particles() {
            return f64::NAN;
        }
        let n = self.expected_particles();
        let backflow = self.active_backflow();
        let mut total = 0.0;
        for block in 0..2 {
            let (offset, n_blk) = self.block_bounds(block);
            if n_blk == 0 {
                continue;
            }
            let chain = BlockChain::build(self, coords, block, backflow);
            for k in 0..n {
                // A^{ik} = dx_i/dr_k = δ_ik I + Σ_{j≠i} (δ_ki − δ_kj) J^{ij}
                // (just δ_ik I when backflow is inactive — the loops below
                // then contribute exact zeros, so a λ = 0 backflow
                // determinant computes bit-identically to a plain one).
                let mut a = vec![[[0.0; DIM]; DIM]; n_blk];
                // ∇²_k x_i = Σ_{j≠i} (δ_ki + δ_kj) L^{ij}
                let mut lap = vec![[0.0; DIM]; n_blk];
                if k >= offset && k < offset + n_blk {
                    a[k - offset] = identity3();
                }
                if let Some(backflow) = backflow {
                    for i in 0..n_blk {
                        let global = offset + i;
                        for j in 0..n {
                            if j == global {
                                continue;
                            }
                            let jacobian = backflow.pair_jacobian(coords, global, j);
                            if k == global {
                                add_matrix_scaled(&mut a[i], &jacobian, 1.0);
                            } else if k == j {
                                add_matrix_scaled(&mut a[i], &jacobian, -1.0);
                            }
                            if k == global || k == j {
                                let laplacian = backflow.pair_laplacian(coords, global, j);
                                for c in 0..DIM {
                                    lap[i][c] += laplacian[c];
                                }
                            }
                        }
                    }
                }
                for i in 0..n_blk {
                    total += contract_q_aat(&chain.q[i], &a[i]);
                    total += dot3(&chain.o[i], &lap[i]);
                }
                // Second-derivative cross term (module docs):
                // −Σ_{i,l in block} (A^{ik} w^{il}) . (A^{lk} w^{li}).
                // Without backflow only i = l = k survives: −|O_k|².
                for i in 0..n_blk {
                    for l in 0..n_blk {
                        let left = mat_vec3(&a[i], &chain.w[i][l]);
                        let right = mat_vec3(&a[l], &chain.w[l][i]);
                        total -= dot3(&left, &right);
                    }
                }
            }
        }
        total
    }

    #[inline]
    fn n_params(&self) -> usize {
        let mut count = 0;
        for block in 0..2 {
            for orbital in &self.orbitals[block] {
                count += 1 + orbital.primitives().len();
            }
        }
        if self.backflow.is_some() {
            count += 1;
        }
        count
    }

    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer) {
        let coords = cfg.as_ref();
        if coords.len() != DIM * self.expected_particles() {
            return;
        }
        let backflow = self.active_backflow();
        let n = self.expected_particles();
        let out = out.as_mut_slice();
        let mut cursor = 0;
        for block in 0..2 {
            let (_, n_blk) = self.block_bounds(block);
            if n_blk == 0 {
                continue;
            }
            let chain = BlockChain::build(self, coords, block, backflow);
            for (m, orbital) in self.orbitals[block].iter().enumerate() {
                // d ln|det| / d alpha_m = sum_i (D^{-1})_{im} (-|x_i|^2 phi_m(x_i)):
                // d/d alpha of e^{-alpha r^2} is -r^2 times itself.
                let mut exponent_grad = 0.0;
                for i in 0..n_blk {
                    let x = &chain.x[i];
                    let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
                    exponent_grad += -r2 * chain.inverse[(i, m)] * chain.values[(m, i)];
                }
                out[cursor] += exponent_grad;
                cursor += 1;
                // d ln|det| / d c_{mk} = sum_i (D^{-1})_{im} g_k(x_i).
                for i in 0..n_blk {
                    let weight = chain.inverse[(i, m)];
                    for (k, g) in orbital
                        .primitive_values(&chain.x[i])
                        .into_iter()
                        .enumerate()
                    {
                        out[cursor + k] += weight * g;
                    }
                }
                cursor += orbital.primitives().len();
            }
        }
        // d ln|psi| / d lambda = sum_blocks sum_{i in block} O_i . (d x_i / d lambda)
        // with d x_i / d lambda = sum_{j != i} f(r_ij) (r_i - r_j)  (eta = lambda f).
        // Evaluated with the currently active chains; at lambda = 0 the
        // plain-Slater chain is exactly the derivative at zero.
        if let Some(backflow) = &self.backflow {
            let mut lambda_grad = 0.0;
            for block in 0..2 {
                let (offset, n_blk) = self.block_bounds(block);
                if n_blk == 0 {
                    continue;
                }
                let chain = BlockChain::build(self, coords, block, Some(backflow));
                for i in 0..n_blk {
                    let global = offset + i;
                    let origin = read_particle(coords, global);
                    let mut direction = [0.0; DIM];
                    for j in 0..n {
                        if j == global {
                            continue;
                        }
                        let other = read_particle(coords, j);
                        let f = backflow.shape_value(point_distance(&origin, &other));
                        for k in 0..DIM {
                            direction[k] += f * (origin[k] - other[k]);
                        }
                    }
                    lambda_grad += dot3(&chain.o[i], &direction);
                }
            }
            out[cursor] += lambda_grad;
        }
    }

    #[inline]
    fn update_params(&mut self, delta: &[f64]) {
        debug_assert_eq!(delta.len(), self.n_params());
        let mut cursor = 0;
        for block in 0..2 {
            for orbital in &mut self.orbitals[block] {
                orbital.exponent += delta[cursor];
                cursor += 1;
                for (coefficient, _) in &mut orbital.primitives {
                    *coefficient += delta[cursor];
                    cursor += 1;
                }
            }
        }
        if let Some(backflow) = &mut self.backflow {
            backflow.add_lambda(delta[cursor]);
        }
        self.cache.borrow_mut().up_to_date = false;
    }

    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog {
        let coords = cfg.as_ref();
        if coords.len() != DIM * self.expected_particles() || particle >= self.expected_particles()
        {
            return DeltaLog {
                log_ratio: f64::NAN,
            };
        }
        let backflow = self.active_backflow();
        if let Some(backflow) = backflow {
            // Backflow: every quasiparticle coordinate changes with one
            // particle — full recompute of the moved log-determinant
            // (Kwon et al. 1998 §II: no rank-1 update survives).
            self.ensure_cache(coords);
            let old = {
                let cache = self.cache.borrow();
                cache.log_det[0] + cache.log_det[1]
            };
            let new = {
                let mut cache = self.cache.borrow_mut();
                let moved = &mut cache.moved;
                moved.copy_from_slice(coords);
                moved[DIM * particle..DIM * particle + DIM].copy_from_slice(new_pos);
                self.total_logdet(moved, Some(backflow))
            };
            DeltaLog {
                log_ratio: new - old,
            }
        } else {
            // Plain Slater: O(N) Sherman–Morrison column ratio against the
            // cached inverse — allocation-free.
            self.ensure_cache(coords);
            let block = if particle < self.n_up { 0 } else { 1 };
            let (offset, n_blk) = self.block_bounds(block);
            let local = particle - offset;
            let log_ratio = {
                let mut cache = self.cache.borrow_mut();
                for (m, orbital) in self.orbitals[block].iter().enumerate() {
                    cache.column[m] = orbital.value(new_pos);
                }
                let mut ratio = 0.0;
                for m in 0..n_blk {
                    ratio += cache.inverse[block][(local, m)] * cache.column[m];
                }
                ratio.abs().ln()
            };
            DeltaLog { log_ratio }
        }
    }

    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point) {
        let length_ok = cfg.as_ref().len() == DIM * self.expected_particles();
        let in_range = particle < cfg.as_ref().len() / DIM;
        let cache_was_current = length_ok && in_range && self.cache.borrow().matches(cfg.as_ref());
        if in_range {
            write_particle(cfg, particle, *new_pos);
        }
        if !length_ok || !in_range {
            return; // defensive: never touch determinant caches for a
                    // mismatched configuration (log_psi reports NaN)
        }
        let coords = cfg.as_ref();
        if self.active_backflow().is_some() || !cache_was_current {
            // Backflow always rebuilds; a stale plain-Slater cache cannot
            // be rank-1 updated, so refresh it from the moved config.
            self.rebuild_coords(coords);
            return;
        }
        // Plain Slater: Sherman–Morrison O(N²) inverse row update.
        let block = if particle < self.n_up { 0 } else { 1 };
        let (offset, n_blk) = self.block_bounds(block);
        let local = particle - offset;
        let column = DMatrix::from_fn(n_blk, 1, |m, _| self.orbitals[block][m].value(new_pos));
        let mut cache = self.cache.borrow_mut();
        let product = &cache.inverse[block] * &column;
        let ratio = product[(local, 0)];
        if !ratio.is_finite() || ratio == 0.0 {
            drop(cache);
            self.rebuild_coords(coords);
            return;
        }
        let mut v = product;
        v[(local, 0)] -= 1.0;
        let updated = &cache.inverse[block] - &v * cache.inverse[block].row(local) / ratio;
        cache.inverse[block] = updated;
        cache.log_det[block] += ratio.abs().ln();
        let base = DIM * particle;
        cache.fingerprint[base..base + DIM].copy_from_slice(new_pos);
    }

    fn rebuild(&mut self, cfg: &Self::Config) {
        let coords = cfg.as_ref();
        if coords.len() == DIM * self.expected_particles() {
            self.rebuild_coords(coords);
        }
    }
}

/// Per-block derivative-chain data: quasiparticle coordinates, orbital
/// values `values[(m, i)] = phi_m(x_i)` (orbital rows, so rows of the LU
/// inverse are particle indices), the contracted
/// `O_i = sum_m (D^{-1})_{im} grad phi_m(x_i)` /
/// `Q_i = sum_m (D^{-1})_{im} Hess phi_m(x_i)` and the cross weights
/// `w^{il} = sum_m (D^{-1})_{lm} grad phi_m(x_i)` of the Laplacian's
/// second-derivative term.
struct BlockChain {
    x: Vec<Point>,
    values: DMatrix<f64>,
    inverse: DMatrix<f64>,
    o: Vec<[f64; DIM]>,
    q: Vec<[[f64; DIM]; DIM]>,
    w: Vec<Vec<[f64; DIM]>>,
}

impl BlockChain {
    /// Build for one block; `inverse`/`o`/`q`/`w` are all-NaN when the
    /// block determinant is singular (the measurement path of a NaN
    /// `log_psi`).
    fn build(
        slater: &SlaterDeterminant,
        coords: &[f64],
        block: usize,
        backflow: Option<&Backflow>,
    ) -> Self {
        let (offset, n_blk) = slater.block_bounds(block);
        let orbitals = &slater.orbitals[block];
        let mut x = Vec::with_capacity(n_blk);
        let mut values = DMatrix::zeros(n_blk, n_blk);
        let mut grads = vec![[0.0; DIM]; n_blk * n_blk];
        let mut hessians = vec![[[0.0; DIM]; DIM]; n_blk * n_blk];
        for i in 0..n_blk {
            let point = quasiparticle(coords, offset + i, backflow);
            x.push(point);
            for (m, orbital) in orbitals.iter().enumerate() {
                let evaluated = orbital.evaluate(&point);
                values[(m, i)] = evaluated.value;
                grads[i * n_blk + m] = evaluated.grad;
                hessians[i * n_blk + m] = evaluated.hessian;
            }
        }
        // The LU runs on the orbital-row matrix; the resulting inverse is
        // indexed [(particle, orbital)] — exactly the weights the chain
        // rules need. (Inverting the transposed, particle-row matrix would
        // silently transpose every weight; the finite-difference tests
        // guard this indexing.)
        let Some((_, inverse)) = lu_logdet_inverse(&values) else {
            let nan_matrix = [[f64::NAN; DIM]; DIM];
            return Self {
                x,
                values,
                inverse: DMatrix::repeat(n_blk, n_blk, f64::NAN),
                o: vec![[f64::NAN; DIM]; n_blk],
                q: vec![nan_matrix; n_blk],
                w: vec![vec![[f64::NAN; DIM]; n_blk]; n_blk],
            };
        };
        let mut o = vec![[0.0; DIM]; n_blk];
        let mut q = vec![[[0.0; DIM]; DIM]; n_blk];
        for i in 0..n_blk {
            for m in 0..n_blk {
                let weight = inverse[(i, m)];
                let grad = &grads[i * n_blk + m];
                let hessian = &hessians[i * n_blk + m];
                for k in 0..DIM {
                    o[i][k] += weight * grad[k];
                    for l in 0..DIM {
                        q[i][k][l] += weight * hessian[k][l];
                    }
                }
            }
        }
        let mut w = vec![vec![[0.0; DIM]; n_blk]; n_blk];
        for i in 0..n_blk {
            for l in 0..n_blk {
                for m in 0..n_blk {
                    let weight = inverse[(l, m)];
                    let grad = &grads[i * n_blk + m];
                    for k in 0..DIM {
                        w[i][l][k] += weight * grad[k];
                    }
                }
            }
        }
        Self {
            x,
            values,
            inverse,
            o,
            q,
            w,
        }
    }
}

#[inline]
fn identity3() -> [[f64; DIM]; DIM] {
    let mut a = [[0.0; DIM]; DIM];
    for (k, row) in a.iter_mut().enumerate() {
        row[k] = 1.0;
    }
    a
}

/// `target += sign * matrix` elementwise.
#[inline]
fn add_matrix_scaled(target: &mut [[f64; DIM]; DIM], matrix: &[[f64; DIM]; DIM], sign: f64) {
    for k in 0..DIM {
        for l in 0..DIM {
            target[k][l] += sign * matrix[k][l];
        }
    }
}

/// `(M v)_a = sum_b M[a][b] v[b]`.
#[inline]
fn mat_vec3(matrix: &[[f64; DIM]; DIM], v: &[f64; DIM]) -> [f64; DIM] {
    let mut out = [0.0; DIM];
    for (slot, row) in out.iter_mut().zip(matrix) {
        *slot = dot3(row, v);
    }
    out
}

/// `out += sign * J v` with the row convention `(J v)_a = sum_b J_ab v_b`.
#[inline]
fn add_jacobian_action(out: &mut [f64], jacobian: &[[f64; DIM]; DIM], v: &[f64; DIM], sign: f64) {
    for (alpha, slot) in out.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (beta, &v_beta) in v.iter().enumerate() {
            acc += jacobian[alpha][beta] * v_beta;
        }
        *slot += sign * acc;
    }
}

/// `sum_{beta,gamma} Q[beta][gamma] (A A^T)[beta][gamma]` — the quadratic
/// Jacobian contraction of the Laplacian chain (Q symmetric, so the
/// transpose bookkeeping is trivial).
#[inline]
fn contract_q_aat(q: &[[f64; DIM]; DIM], a: &[[f64; DIM]; DIM]) -> f64 {
    let mut total = 0.0;
    for (beta, q_row) in q.iter().enumerate() {
        for (gamma, &q_bg) in q_row.iter().enumerate() {
            let mut aat = 0.0;
            for (a_beta, a_gamma) in a[beta].iter().zip(a[gamma]) {
                aat += a_beta * a_gamma;
            }
            total += q_bg * aat;
        }
    }
    total
}

#[inline]
fn dot3(a: &[f64; DIM], b: &[f64; DIM]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn point_distance(a: &Point, b: &Point) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

impl WaveFunctionParams for SlaterDeterminant {
    fn param_values(&self) -> Vec<f64> {
        let mut values = Vec::with_capacity(self.n_params());
        for block in 0..2 {
            for orbital in &self.orbitals[block] {
                values.push(orbital.exponent);
                for (coefficient, _) in &orbital.primitives {
                    values.push(*coefficient);
                }
            }
        }
        if let Some(backflow) = &self.backflow {
            values.push(backflow.lambda());
        }
        values
    }

    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError> {
        if values.len() != self.n_params() {
            return Err(VariationalError::invalid(
                "params",
                format!("expected {} values, got {}", self.n_params(), values.len()),
            ));
        }
        // Validate everything first: a rejection leaves self untouched.
        let mut cursor = 0;
        for block in 0..2 {
            for orbital in &self.orbitals[block] {
                VariationalError::require_positive("exponent", values[cursor])?;
                for _ in orbital.primitives() {
                    if !values[cursor + 1].is_finite() {
                        return Err(VariationalError::invalid(
                            "coefficient",
                            format!("must be finite, got {}", values[cursor + 1]),
                        ));
                    }
                    cursor += 1;
                }
                cursor += 1;
            }
        }
        if self.backflow.is_some() && !values[cursor].is_finite() {
            return Err(VariationalError::invalid(
                "lambda",
                format!("must be finite, got {}", values[cursor]),
            ));
        }
        // Assign (validated above).
        let mut cursor = 0;
        for block in 0..2 {
            for orbital in &mut self.orbitals[block] {
                orbital.exponent = values[cursor];
                cursor += 1;
                for (coefficient, _) in &mut orbital.primitives {
                    *coefficient = values[cursor];
                    cursor += 1;
                }
            }
        }
        if let Some(backflow) = &mut self.backflow {
            backflow.set_lambda(values[cursor])?;
        }
        self.cache.borrow_mut().up_to_date = false;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gto_validation_rejects_garbage() {
        assert!(GtoOrbital::new(0.5, vec![(1.0, [0, 0, 0])]).is_ok());
        for bad_alpha in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(GtoOrbital::new(bad_alpha, vec![(1.0, [0, 0, 0])]).is_err());
        }
        assert!(GtoOrbital::new(0.5, vec![]).is_err());
        assert!(GtoOrbital::new(0.5, vec![(f64::NAN, [0, 0, 0])]).is_err());
        assert!(GtoOrbital::new(0.5, vec![(1.0, [9, 0, 0])]).is_err());
    }

    #[test]
    fn gaussian_orbital_closed_forms() {
        // phi = e^{-alpha r^2}: grad = -2 alpha t phi,
        // Hess = (-2 alpha I + 4 alpha^2 t t^T) phi.
        let alpha = 0.7;
        let s_orbital = GtoOrbital::new(alpha, vec![(1.0, [0, 0, 0])]).unwrap();
        let t = [0.3, -0.4, 1.2];
        let r2 = 0.09 + 0.16 + 1.44;
        let evaluated = s_orbital.evaluate(&t);
        let phi = (-alpha * r2).exp();
        assert!((evaluated.value - phi).abs() < 1e-15);
        for (&grad_k, t_k) in evaluated.grad.iter().zip(t) {
            assert!((grad_k + 2.0 * alpha * t_k * phi).abs() < 1e-14);
        }
        assert!((evaluated.hessian[1][2] - 4.0 * alpha * alpha * t[1] * t[2] * phi).abs() < 1e-14);
        assert!(
            (evaluated.hessian[0][0] - (-2.0 * alpha + 4.0 * alpha * alpha * t[0] * t[0]) * phi)
                .abs()
                < 1e-14
        );

        // p orbital phi = y e^{-alpha r^2} on the y axis: by-hand
        // product-rule values.
        let p_orbital = GtoOrbital::new(alpha, vec![(1.0, [0, 1, 0])]).unwrap();
        let t = [0.0, 0.9, 0.0];
        let g = (-alpha * 0.81).exp();
        let evaluated = p_orbital.evaluate(&t);
        assert!((evaluated.value - 0.9 * g).abs() < 1e-15);
        assert!((evaluated.grad[1] - g * (1.0 - 2.0 * alpha * 0.81)).abs() < 1e-14);
        assert!(
            (evaluated.hessian[1][1] - g * (-6.0 * alpha * 0.9 + 4.0 * alpha * alpha * 0.729))
                .abs()
                < 1e-14
        );
        assert!(evaluated.hessian[0][1].abs() < 1e-15);
    }

    #[test]
    fn hermite_polynomials_match_closed_forms() {
        // H_2(sqrt(w) t) = 4 w t^2 - 2, H_3 = 8 w^{3/2} t^3 - 12 sqrt(w) t.
        let h2 = hermite_coeffs(2, 1.3);
        assert!((h2[0] + 2.0).abs() < 1e-14);
        assert!(h2[1].abs() < 1e-14);
        assert!((h2[2] - 4.0 * 1.3).abs() < 1e-14);
        let h3 = hermite_coeffs(3, 1.3);
        assert!((h3[3] - 8.0 * 1.3 * 1.3_f64.sqrt()).abs() < 1e-14);
        assert!((h3[1] + 12.0 * 1.3_f64.sqrt()).abs() < 1e-14);
    }

    #[test]
    fn closed_shell_counts_and_energies() {
        assert_eq!(harmonic_closed_shell_electrons(1).unwrap(), 2);
        assert_eq!(harmonic_closed_shell_electrons(2).unwrap(), 8);
        assert_eq!(harmonic_closed_shell_electrons(3).unwrap(), 20);
        assert!((harmonic_closed_shell_energy(1.7, 1).unwrap() - 3.0 * 1.7).abs() < 1e-12);
        assert!((harmonic_closed_shell_energy(1.7, 2).unwrap() - 18.0 * 1.7).abs() < 1e-12);
        assert!((harmonic_closed_shell_energy(1.7, 3).unwrap() - 60.0 * 1.7).abs() < 1e-12);
        assert!(harmonic_closed_shell_electrons(0).is_err());
        assert!(harmonic_closed_shell_electrons(4).is_err());
        assert!(harmonic_trap_orbitals(-1.0, 1).is_err());
        assert_eq!(harmonic_trap_orbitals(1.0, 1).unwrap().len(), 1);
        assert_eq!(harmonic_trap_orbitals(1.0, 2).unwrap().len(), 4);
        assert_eq!(harmonic_trap_orbitals(1.0, 3).unwrap().len(), 10);
        // The first shell-2 orbital H_2(sqrt(w) x) = 4 w x^2 - 2 carries
        // exactly the two-primitive contraction.
        let shells = harmonic_trap_orbitals(1.0, 3).unwrap();
        assert_eq!(shells[4].primitives().len(), 2);
    }

    #[test]
    fn single_particle_determinant_matches_orbital() {
        // One up electron, empty down block: psi = phi, so
        // ln|psi| = ln phi, grad ln|psi| = grad phi / phi, etc.
        let slater = SlaterDeterminant::new(
            vec![GtoOrbital::new(0.5, vec![(1.0, [0, 0, 0])]).unwrap()],
            vec![],
        )
        .unwrap();
        assert_eq!(slater.n_up(), 1);
        assert_eq!(slater.n_down(), 0);
        let cfg = Positions::from_flat(vec![0.4, -0.2, 0.7]).unwrap();
        let alpha = 0.5;
        let r2 = 0.16 + 0.04 + 0.49;
        assert!((slater.log_psi(&cfg) - (-alpha * r2)).abs() < 1e-15);

        let mut grad = GradBuffer::new(1);
        slater.log_grad(&cfg, &mut grad);
        for k in 0..DIM {
            assert!((grad.as_slice()[k] + 2.0 * alpha * cfg.as_ref()[k]).abs() < 1e-14);
        }
        assert!((slater.log_laplacian(&cfg) + 6.0 * alpha).abs() < 1e-14);
    }

    #[test]
    fn configuration_length_mismatch_is_nan_not_panic() {
        let slater = SlaterDeterminant::harmonic_trap(1.0, 1).unwrap();
        // `Positions` guarantees a multiple of DIM; a particle-count
        // mismatch (4 particles vs the expected 2) is the reachable case.
        let wrong = Positions::from_flat(vec![0.0; 12]).unwrap();
        assert!(slater.log_psi(&wrong).is_nan());
        assert!(slater.log_laplacian(&wrong).is_nan());
        assert!(slater.delta_log(&wrong, 0, &[0.0; DIM]).log_ratio.is_nan());
        let mut grad = GradBuffer::new(1);
        slater.log_grad(&wrong, &mut grad);
        assert!(grad.as_slice().iter().all(|&x| x == 0.0));
        // Out-of-range particle on a right-sized config is also NaN.
        let right = Positions::from_flat(vec![0.3; 6]).unwrap();
        assert!(slater.delta_log(&right, 7, &[0.0; DIM]).log_ratio.is_nan());
    }

    #[test]
    fn same_spin_exchange_flips_sign_of_det() {
        // Swapping two same-spin particles swaps two columns:
        // det -> -det, so ln|det| is invariant (LU pivoting differs, so
        // compare at roundoff, not bit-exact).
        let slater = SlaterDeterminant::harmonic_trap(1.0, 2).unwrap();
        let mut cfg = Positions::from_flat(vec![
            0.6, -0.2, 0.4, //
            -0.9, 0.3, 0.1, //
            0.2, 1.1, -0.5, //
            0.0, -0.4, 0.8, //
            1.3, -0.6, 0.2, //
            -1.1, 0.5, 0.0, //
            0.4, 0.7, -0.3, //
            0.1, 0.2, 0.9,
        ])
        .unwrap();
        let direct = slater.log_psi(&cfg);
        let a = cfg.particle(0);
        let b = cfg.particle(2);
        cfg.set_particle(0, b);
        cfg.set_particle(2, a);
        let swapped = slater.log_psi(&cfg);
        assert!((direct - swapped).abs() < 1e-12);
    }

    #[test]
    fn parameter_layout_round_trip_and_validation() {
        let slater = SlaterDeterminant::harmonic_trap_with_backflow(
            1.0,
            1,
            Backflow::new_electron_gas_shape(0.1).unwrap(),
        )
        .unwrap();
        // Shell 0 for both spins: one orbital each with one primitive ->
        // (exponent, coefficient) x 2 + lambda = 5 parameters.
        assert_eq!(slater.n_params(), 5);
        let values = slater.param_values();
        assert_eq!(values.len(), 5);
        assert_eq!(values[0], 0.5); // exponent
        assert_eq!(values[1], 1.0); // coefficient
        assert_eq!(*values.last().unwrap(), 0.1); // lambda

        let mut updated = slater.clone();
        assert!(updated.try_set_params(&[0.5, 1.0, 0.5, 1.0, 0.1]).is_ok());
        assert_eq!(updated.param_values(), vec![0.5, 1.0, 0.5, 1.0, 0.1]);
        assert!(updated.try_set_params(&values).is_ok());
        // Rejections leave the state untouched.
        assert!(updated.try_set_params(&[0.0, 1.0, 0.5, 1.0, 0.1]).is_err());
        assert!(updated
            .try_set_params(&[f64::NAN, 1.0, 0.5, 1.0, 0.1])
            .is_err());
        assert!(updated
            .try_set_params(&[0.5, 1.0, 0.5, f64::NAN, 0.1])
            .is_err());
        assert!(updated
            .try_set_params(&[0.5, 1.0, 0.5, 1.0, f64::INFINITY])
            .is_err());
        assert!(updated.try_set_params(&[0.5]).is_err());
        assert_eq!(updated.param_values(), values);

        // update_params invalidates the cache: the next log_psi reflects
        // the shifted exponent.
        let cfg = Positions::from_flat(vec![0.4, -0.1, 0.2, 0.5, 0.3, -0.2]).unwrap();
        let before = updated.log_psi(&cfg);
        updated.update_params(&[0.1, 0.0, 0.1, 0.0, 0.0]);
        let after = updated.log_psi(&cfg);
        assert!(before.is_finite() && after.is_finite() && (before - after).abs() > 1e-3);
    }

    #[test]
    fn constructor_rejects_empty_determinants() {
        assert!(SlaterDeterminant::new(vec![], vec![]).is_err());
        // A single block (fully polarized) is fine.
        let polarized = SlaterDeterminant::new(harmonic_trap_orbitals(1.0, 1).unwrap(), vec![]);
        assert!(polarized.is_ok());
        assert_eq!(polarized.unwrap().expected_particles(), 1);
    }

    #[test]
    fn shermann_morrison_updates_match_fresh_lu() {
        // Two complementary views of the DESIGN.md L1 gate
        // ("Sherman-Morrison rank-1 update == full inverse, 1e-12"):
        //
        // 1. The identity itself on conditioning-controlled matrices:
        //    replacing one column of a well-conditioned O(1) matrix and
        //    applying the rank-1 inverse update reproduces the fresh
        //    inverse to <= 1e-12 entrywise.
        // 2. The deployed path: long accept-only move sequences through
        //    delta_log/commit_move with the kernel's per-pass `rebuild`
        //    re-anchoring. The quantities the kernel actually consumes —
        //    the Metropolis ratio and the accumulated log-determinant —
        //    agree with fresh LU recomputes to <= 1e-12. (Absolute
        //    inverse ENTRIES drift to ~1e-9..1e-6 at force-accepted
        //    near-cancellation ratios: a single update multiplies by
        //    1/ratio ~ e^{|log ratio|}, so the entry error is the
        //    condition-number floor eps * cond / |ratio|, not an
        //    implementation error — measured against an exact-arithmetic
        //    Python reference during development. Ratios and log-dets —
        //    all the walk consumes — remain at the 1e-13 level.)
        use rand::{RngExt, SeedableRng};
        use rand_xoshiro::Xoshiro256PlusPlus;

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5EED1);

        // -- Part 1: the identity on conditioning-controlled draws. -------
        let mut worst_identity = 0.0_f64;
        let mut worst_ratio = 0.0_f64;
        for _ in 0..200 {
            let n = rng.random_range(2..=6);
            // Well-conditioned O(1) matrices, guaranteed nonsingular by
            // construction (no unbounded rejection loop): strictly
            // diagonally dominant with diag 1.2n vs off-diagonal row sums
            // <= 0.85(n-1), so the Varah bound gives
            // |det| >= (0.35n + 0.85)^n > 0.05 in a single draw.
            let matrix = DMatrix::from_fn(n, n, |m, i| {
                let spread = 1.0 + 0.6 * (m as f64 - i as f64);
                let base = (-((m as f64 - i as f64 + 3.0 * (m % 3) as f64) / spread).powi(2)).exp();
                if m == i {
                    1.2 * n as f64 + 0.05 * rng.random_range(-1.0..1.0)
                } else {
                    0.8 * base + 0.05 * rng.random_range(-1.0..1.0)
                }
            });
            let Some((log_det, inverse)) = lu_logdet_inverse(&matrix) else {
                panic!("LU of a conditioning-controlled matrix failed");
            };
            let column = DMatrix::from_fn(n, 1, |m, _| {
                (-((m as f64 - 2.0) * 0.4).powi(2)).exp() + 0.05 * rng.random_range(-1.0..1.0)
            });
            // Moved matrix: column `local` replaced.
            let local = rng.random_range(0..n);
            let mut moved = matrix.clone();
            moved.set_column(local, &column.column(0));
            let Some((moved_log_det, moved_inverse)) = lu_logdet_inverse(&moved) else {
                panic!("LU of the moved matrix failed");
            };
            // The delta_log convention: ratio = (D^{-1} d)_local.
            let product = &inverse * &column;
            let ratio = product[(local, 0)];
            worst_ratio = worst_ratio.max((ratio.abs().ln() - (moved_log_det - log_det)).abs());
            if !(0.05..20.0).contains(&ratio.abs()) {
                continue; // keep only moderate ratios for the entry check
            }
            // The commit_move convention: rank-1 row update.
            let mut v = product;
            v[(local, 0)] -= 1.0;
            let updated = &inverse - &v * inverse.row(local) / ratio;
            for i in 0..n {
                for m in 0..n {
                    worst_identity =
                        worst_identity.max((updated[(i, m)] - moved_inverse[(i, m)]).abs());
                }
            }
        }
        assert!(
            worst_identity <= 1e-12,
            "SM identity drift: {worst_identity}"
        );
        assert!(worst_ratio <= 1e-12, "SM log-ratio drift: {worst_ratio}");

        // -- Part 2: the deployed path through delta_log/commit_move. ----
        let mut slater = SlaterDeterminant::harmonic_trap(1.2, 2).unwrap(); // 8 electrons
        let n = slater.expected_particles();
        let mut cfg =
            Positions::from_flat((0..DIM * n).map(|_| rng.random_range(-1.5..1.5)).collect())
                .unwrap();
        let mut worst_walk_ratio = 0.0_f64;
        let mut worst_log_det = 0.0_f64;
        slater.log_psi(&cfg); // populate the cache
                              // 50 kernel-like walker passes: n accepted single-particle moves,
                              // then the per-pass `rebuild` re-anchor (the kernel's K-rebuild
                              // policy). Without the re-anchor, 400 raw Sherman-Morrison
                              // accumulations drift to ~2e-9 — the drift the policy bounds.
        for _ in 0..50 {
            for _ in 0..n {
                let particle = rng.random_range(0..n);
                let old = cfg.particle(particle);
                let new = [
                    old[0] + rng.random_range(-0.4..0.4),
                    old[1] + rng.random_range(-0.4..0.4),
                    old[2] + rng.random_range(-0.4..0.4),
                ];
                let delta = slater.delta_log(&cfg, particle, &new);
                let mut moved = cfg.clone();
                moved.set_particle(particle, new);
                // Fresh recomputes via a witness clone whose cache is force-
                // rebuilt (the test must NOT re-anchor the ansatz under
                // test).
                let mut witness = slater.clone();
                witness.rebuild(&moved);
                let fresh_plus = witness.log_psi(&moved);
                witness.rebuild(&cfg);
                let fresh_old = witness.log_psi(&cfg);
                worst_walk_ratio =
                    worst_walk_ratio.max((delta.log_ratio - (fresh_plus - fresh_old)).abs());
                slater.commit_move(&mut cfg, particle, &new);
                // Accumulated log-determinant (SM chain) vs fresh LU.
                let cached = slater.log_psi(&cfg); // cache is current
                let mut fresh_total = 0.0;
                for block in 0..2 {
                    let Some((log_det, _)) = slater.block_lu(cfg.as_ref(), block, None) else {
                        panic!("fresh LU of a nonsingular configuration failed");
                    };
                    fresh_total += log_det;
                }
                worst_log_det =
                    worst_log_det.max((cached - fresh_total).abs() / (1.0 + fresh_total.abs()));
            }
            slater.rebuild(&cfg);
        }
        assert!(
            worst_walk_ratio <= 1e-12,
            "walk ratio vs full recompute: {worst_walk_ratio}"
        );
        assert!(worst_log_det <= 1e-12, "SM log-det drift: {worst_log_det}");
    }
}
