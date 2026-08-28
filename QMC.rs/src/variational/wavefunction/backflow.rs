//! Backflow coordinate transformation (L1).
//!
//! Implements the quasiparticle displacement of Kwon, Ceperley & Martin,
//! Phys. Rev. B 58, 6800 (1998), eq. (6) (the same form appears as eq. (6)
//! of Holzmann et al. 2019):
//!
//! ```text
//! x_i = r_i + sum_{j != i} eta(r_ij) (r_i - r_j),
//! ```
//!
//! with the radial backflow function of eq. (10) of the same paper,
//!
//! ```text
//! eta(r) = lambda_B (1 + s_B r) / (r_B + w_B r + r^4),
//! ```
//!
//! rewritten here as `eta(r) = lambda * f(r)` with the fixed shape
//! `f(r) = (1 + s_B r)/(r_B + w_B r + r^4)` and the single variational
//! scale `lambda`. Keeping the factorization explicit is what makes the
//! `lambda = 0` reduction to the plain Slater determinant **bit-exact**:
//! every quantity below is computed as `lambda * (finite bracket)`, so at
//! `lambda = 0` all displacements, Jacobians and Laplacians are exactly
//! zero (and [`SlaterDeterminant`](super::determinant::SlaterDeterminant)
//! additionally routes around the backflow chains entirely when
//! `lambda == 0`, so no reliance on zero-arithmetic remains).
//!
//! The shape parameters follow the Kwon et al. Table I semantics: `r_B`
//! controls the short-distance scale (their optimized values sit near
//! 0.2–0.3), `w_B` the intermediate range and `s_B` the `~1/r^3` tail
//! weight. They are construction-time constants at L1; only `lambda` is
//! variational (one parameter, owned by the host determinant's parameter
//! vector).
//!
//! # Hand-derived derivative chain
//!
//! Elementary derivations (verified against central finite differences in
//! CI). Write the pair displacement `u_ij = eta(r_ij) t_ij` with
//! `t_ij = r_i - r_j` and `r = |t|`:
//!
//! - Jacobian w.r.t. either endpoint (a symmetric 3×3 matrix):
//!
//! ```text
//! du_ij / dr_i = J^{ij} = eta'(r) t t^T / r + eta(r) I,    du_ij / dr_j = -J^{ij}
//! ```
//!
//!   (from `dr/dt_alpha = t_alpha / r` and `du_beta/dt_alpha =
//!   eta'(r) (t_alpha/r) t_beta + eta(r) delta_{alpha beta}`).
//!
//! - Laplacian w.r.t. either endpoint (same sign: `d/dr_j = -d/dr_t`):
//!
//! ```text
//! lap_{r_i} u_ij = lap_{r_j} u_ij = (eta''(r) + 4 eta'(r)/r) t_ij =: L^{ij}
//! ```
//!
//!   Trace of `d^2 u_beta / dt_alpha dt_gamma` (the full matrix is in the
//!   derivation comment on [`Backflow::pair_laplacian`]): the `eta''` term
//!   contributes `eta'' t_beta`, the three `eta'/r` terms contribute
//!   `3 eta' t_beta / r` and `eta' t_beta / r` respectively, the
//!   `eta'/r^3` term cancels one of them — net `(eta'' + 4 eta'/r) t`.
//!
//! - Derivative of the scale (`eta = lambda f`):
//!
//! ```text
//! d x_i / d lambda = sum_{j != i} f(r_ij) t_ij
//! ```
//!
//! # Small-distance guard
//!
//! At exactly `r_ij = 0` the pair direction `t/r` is undefined (measure
//! zero for continuous sampling). Pairs closer than [`SEPARATION_FLOOR`]
//! contribute the smooth limits `J = eta(0) I` and `L = 0` instead of
//! `0/0 = NaN` — never panic, never poison the walk.

use super::super::error::VariationalError;
use super::{read_particle, Point, DIM};

/// Below this pair distance the backflow pair terms use their smooth
/// `r -> 0` limits (direction undefined at exactly zero; probability-zero
/// event in continuous sampling, guarded only against NaN).
const SEPARATION_FLOOR: f64 = 1e-12;

/// Kwon–Ceperley–Martin backflow transformation (see module docs).
///
/// `lambda` is the single variational scale; `(s_b, r_b, w_b)` are fixed
/// shape parameters (Kwon et al. 1998 eq. (10)). The denominator
/// `r_b + w_b r + r^4` is strictly positive for `r >= 0` given the
/// constructor validation (`r_b > 0`, `w_b >= 0`), so `f`, `f'`, `f''` are
/// finite everywhere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Backflow {
    lambda: f64,
    s_b: f64,
    r_b: f64,
    w_b: f64,
}

impl Backflow {
    /// Construct with a variational scale `lambda` (any finite value —
    /// zero is the legitimate "off" value) and the Kwon et al. eq. (10)
    /// shape parameters: `s_b` finite, `r_b` finite and strictly positive,
    /// `w_b` finite and non-negative.
    pub fn new(lambda: f64, s_b: f64, r_b: f64, w_b: f64) -> Result<Self, VariationalError> {
        if !lambda.is_finite() {
            return Err(VariationalError::invalid(
                "lambda",
                format!("must be finite, got {lambda}"),
            ));
        }
        if !s_b.is_finite() {
            return Err(VariationalError::invalid(
                "s_b",
                format!("must be finite, got {s_b}"),
            ));
        }
        VariationalError::require_positive("r_b", r_b)?;
        if !w_b.is_finite() || w_b < 0.0 {
            return Err(VariationalError::invalid(
                "w_b",
                format!("must be finite and non-negative, got {w_b}"),
            ));
        }
        Ok(Self {
            lambda,
            s_b,
            r_b,
            w_b,
        })
    }

    /// Construct with the electron-gas-like shape of Kwon et al. Table I
    /// (`s_b = 0.4`, `r_b = 0.25`, `w_b = 0.7`, near their `r_s = 1`
    /// optimum) and only the scale variational.
    pub fn new_electron_gas_shape(lambda: f64) -> Result<Self, VariationalError> {
        Self::new(lambda, 0.4, 0.25, 0.7)
    }

    /// The variational scale `lambda`.
    #[inline]
    pub const fn lambda(&self) -> f64 {
        self.lambda
    }

    /// Set the scale (validated; used by the host determinant's parameter
    /// paths).
    pub fn set_lambda(&mut self, lambda: f64) -> Result<(), VariationalError> {
        if !lambda.is_finite() {
            return Err(VariationalError::invalid(
                "lambda",
                format!("must be finite, got {lambda}"),
            ));
        }
        self.lambda = lambda;
        Ok(())
    }

    /// Additive scale update — the optimizer-facing `update_params` path of
    /// the host determinant (same no-validation convention as the L0
    /// ansätze: the kernel constructor and `try_set_params` enforce
    /// finiteness).
    pub(crate) fn add_lambda(&mut self, delta: f64) {
        self.lambda += delta;
    }

    /// Shape value `f(r) = (1 + s_b r)/(r_b + w_b r + r^4)` (the `lambda`
    /// bracket).
    #[inline]
    fn shape(&self, r: f64) -> f64 {
        (1.0 + self.s_b * r) / self.denominator(r)
    }

    /// Shape derivatives `(f, f', f'')` in one pass. With
    /// `D = r_b + w_b r + r^4`, `D' = w_b + 4 r^3`, `N = s_b D - (1 + s_b r) D'`
    /// (so `f' = N / D^2`), and — because the `s_b D'` terms cancel —
    /// `N' = -(1 + s_b r) D'' = -(1 + s_b r) 12 r^2`, hence
    /// `f'' = (N' D - 2 N D') / D^3`.
    #[inline]
    fn shape_derivs(&self, r: f64) -> (f64, f64, f64) {
        let d = self.denominator(r);
        let d_prime = self.w_b + 4.0 * r * r * r;
        let numerator = self.s_b * d - (1.0 + self.s_b * r) * d_prime;
        let d_double_prime = 12.0 * r * r;
        let numerator_prime = -(1.0 + self.s_b * r) * d_double_prime;
        let f = (1.0 + self.s_b * r) / d;
        let f_prime = numerator / (d * d);
        let f_double_prime = (numerator_prime * d - 2.0 * numerator * d_prime) / (d * d * d);
        (f, f_prime, f_double_prime)
    }

    #[inline]
    fn denominator(&self, r: f64) -> f64 {
        self.r_b + self.w_b * r + r * r * r * r
    }

    /// `eta(r) = lambda f(r)`.
    #[inline]
    pub fn eta(&self, r: f64) -> f64 {
        self.lambda * self.shape(r)
    }

    /// `(eta, eta', eta'') = lambda (f, f', f'')`.
    #[inline]
    pub fn eta_derivs(&self, r: f64) -> (f64, f64, f64) {
        let (f, f_prime, f_double_prime) = self.shape_derivs(r);
        (
            self.lambda * f,
            self.lambda * f_prime,
            self.lambda * f_double_prime,
        )
    }

    /// The `lambda` bracket `f(r)` for parameter-gradient chains
    /// (`d x_i / d lambda = sum_j f(r_ij) t_ij`).
    #[inline]
    pub fn shape_value(&self, r: f64) -> f64 {
        self.shape(r)
    }

    /// Backflow displacement of particle `index`: `x_index - r_index =
    /// sum_{j != index} eta(r_indexj) (r_index - r_j)`.
    pub fn displacement<C: AsRef<[f64]> + ?Sized>(&self, cfg: &C, index: usize) -> Point {
        let n = cfg.as_ref().len() / DIM;
        let origin = read_particle(cfg, index);
        let mut out = [0.0; DIM];
        for j in 0..n {
            if j == index {
                continue;
            }
            let other = read_particle(cfg, j);
            let r = point_distance(&origin, &other);
            let eta = self.eta(r);
            for k in 0..DIM {
                out[k] += eta * (origin[k] - other[k]);
            }
        }
        out
    }

    /// Jacobian `J^{ij}` of the pair displacement `u_ij` w.r.t. `r_i`
    /// (symmetric; the `r_j` derivative is its negative):
    /// `J = eta'(r) t t^T / r + eta(r) I`.
    pub fn pair_jacobian<C: AsRef<[f64]> + ?Sized>(
        &self,
        cfg: &C,
        i: usize,
        j: usize,
    ) -> [[f64; 3]; 3] {
        let a = read_particle(cfg, i);
        let b = read_particle(cfg, j);
        let t = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let r2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
        let (eta, eta_prime, _) = self.eta_derivs(r2.sqrt());
        if r2 < SEPARATION_FLOOR * SEPARATION_FLOOR {
            // Smooth r -> 0 limit: the t t^T / r term vanishes, J -> eta(0) I.
            return [[eta, 0.0, 0.0], [0.0, eta, 0.0], [0.0, 0.0, eta]];
        }
        let scale = eta_prime / r2.sqrt();
        let mut out = [[0.0; 3]; 3];
        for alpha in 0..DIM {
            for beta in 0..DIM {
                out[alpha][beta] = scale * t[alpha] * t[beta];
            }
            out[alpha][alpha] += eta;
        }
        out
    }

    /// Laplacian vector `L^{ij} = (eta''(r) + 4 eta'(r)/r) t_ij` of the
    /// pair displacement w.r.t. either endpoint.
    ///
    /// Derivation: `u_beta = eta(r) t_beta` gives
    /// `d^2 u_beta / dt_alpha dt_gamma = eta'' t_alpha t_beta t_gamma / r^2
    /// + eta'[(delta_{alpha gamma} t_beta + t_alpha delta_{beta gamma})/r
    /// - t_alpha t_beta t_gamma/r^3] + eta' (t_gamma/r) delta_{alpha beta}`;
    /// tracing `alpha = gamma` collapses the bracket to `4 t_beta / r`
    /// (the `1/r^3` term cancels against `r^2/r^3` from the `eta''` trace),
    /// leaving `(eta'' + 4 eta'/r) t_beta`.
    pub fn pair_laplacian<C: AsRef<[f64]> + ?Sized>(&self, cfg: &C, i: usize, j: usize) -> Point {
        let a = read_particle(cfg, i);
        let b = read_particle(cfg, j);
        let t = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let r2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
        if r2 < SEPARATION_FLOOR * SEPARATION_FLOOR {
            // Direction undefined at r = 0; the isotropic eta(0) I part of
            // the displacement has zero Laplacian, and the (measure-zero)
            // eta'(0)/r direction term is set to zero rather than NaN.
            return [0.0; DIM];
        }
        let r = r2.sqrt();
        let (_, eta_prime, eta_double_prime) = self.eta_derivs(r);
        let scale = eta_double_prime + 4.0 * eta_prime / r;
        [scale * t[0], scale * t[1], scale * t[2]]
    }
}

#[inline]
fn point_distance(a: &Point, b: &Point) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> Vec<f64> {
        vec![
            0.4, -0.1, 0.2, // 0
            1.2, 0.5, -0.3, // 1
            -0.7, 0.8, 0.1, // 2
        ]
    }

    #[test]
    fn constructor_validates_shape_parameters() {
        assert!(Backflow::new(0.1, 0.4, 0.25, 0.7).is_ok());
        assert!(Backflow::new(0.0, 0.4, 0.25, 0.7).is_ok());
        for bad_lambda in [f64::NAN, f64::INFINITY] {
            assert!(Backflow::new(bad_lambda, 0.4, 0.25, 0.7).is_err());
        }
        for bad_s in [f64::NAN, f64::INFINITY] {
            assert!(Backflow::new(0.1, bad_s, 0.25, 0.7).is_err());
        }
        for bad_r in [0.0, -1.0, f64::NAN] {
            assert!(Backflow::new(0.1, 0.4, bad_r, 0.7).is_err());
        }
        for bad_w in [-1.0, f64::NAN] {
            assert!(Backflow::new(0.1, 0.4, 0.25, bad_w).is_err());
        }
        assert!(Backflow::new(0.1, 0.4, 0.25, 0.0).is_ok());
        assert!(Backflow::new_electron_gas_shape(0.2).is_ok());
    }

    #[test]
    fn eta_closed_form_and_derivatives() {
        let backflow = Backflow::new(0.3, 0.4, 0.25, 0.7).unwrap();
        let r = 1.1_f64;
        let d = 0.25 + 0.7 * r + r.powi(4);
        let expected = 0.3 * (1.0 + 0.4 * r) / d;
        let (eta, _, _) = backflow.eta_derivs(r);
        assert!((eta - expected).abs() < 1e-15);

        // Central differences on eta, eta', eta''.
        let h = 1e-6;
        let (fp, fm, _) = (backflow.eta(r + h), backflow.eta(r - h), 0.0);
        let (_, eta_prime, _) = backflow.eta_derivs(r);
        assert!((eta_prime - (fp - fm) / (2.0 * h)).abs() < 1e-6);
        let (_, _, eta_double_prime) = backflow.eta_derivs(r);
        let second = (backflow.eta(r + h) - 2.0 * eta + backflow.eta(r - h)) / (h * h);
        assert!((eta_double_prime - second).abs() < 1e-4);

        // lambda factorization: eta derivatives scale linearly with lambda.
        let doubled = Backflow::new(0.6, 0.4, 0.25, 0.7).unwrap();
        let (e1, p1, s1) = backflow.eta_derivs(r);
        let (e2, p2, s2) = doubled.eta_derivs(r);
        assert!((2.0 * e1 - e2).abs() < 1e-15);
        assert!((2.0 * p1 - p2).abs() < 1e-15);
        assert!((2.0 * s1 - s2).abs() < 1e-15);
    }

    #[test]
    fn zero_lambda_displacements_are_exactly_zero() {
        let backflow = Backflow::new(0.0, 0.4, 0.25, 0.7).unwrap();
        let cfg = config();
        let displacement = backflow.displacement(&cfg, 1);
        assert_eq!(displacement, [0.0; DIM]);
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert_eq!(backflow.pair_jacobian(&cfg, i, j), [[0.0; 3]; 3]);
                    assert_eq!(backflow.pair_laplacian(&cfg, i, j), [0.0; DIM]);
                }
            }
        }
    }

    #[test]
    fn coincident_pairs_produce_finite_limits_not_nan() {
        let backflow = Backflow::new(0.5, 0.4, 0.25, 0.7).unwrap();
        let cfg = vec![0.3, 0.2, -0.1, 0.3, 0.2, -0.1, 1.0, 0.0, 0.0];
        let jacobian = backflow.pair_jacobian(&cfg, 0, 1);
        let laplacian = backflow.pair_laplacian(&cfg, 0, 1);
        for row in &jacobian {
            assert!(row.iter().all(|x| x.is_finite()));
        }
        assert!(laplacian.iter().all(|x| x.is_finite()));
        assert_eq!(jacobian[0][0], backflow.eta(0.0));
        assert_eq!(laplacian, [0.0; DIM]);
    }

    #[test]
    fn jacobian_and_laplacian_match_finite_differences() {
        // The pair displacement u_ij(cfg) differentiated w.r.t. r_i and r_j
        // components by central differences.
        let backflow = Backflow::new(0.3, 0.4, 0.25, 0.7).unwrap();
        let base = config();
        let h = 1e-6;

        let displacement_at = |cfg: &[f64], i: usize, j: usize| -> [f64; 3] {
            let a = [cfg[3 * i], cfg[3 * i + 1], cfg[3 * i + 2]];
            let b = [cfg[3 * j], cfg[3 * j + 1], cfg[3 * j + 2]];
            let t = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let r = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let eta = backflow.eta(r);
            [eta * t[0], eta * t[1], eta * t[2]]
        };

        // du/d(r_0[alpha]): finite differences vs J^{01} rows.
        let jacobian = backflow.pair_jacobian(&base, 0, 1);
        for alpha in 0..3 {
            let mut plus = base.clone();
            let mut minus = base.clone();
            plus[alpha] += h;
            minus[alpha] -= h;
            for (beta, jac_entry) in jacobian[alpha].iter().enumerate() {
                let fd = (displacement_at(&plus, 0, 1)[beta] - displacement_at(&minus, 0, 1)[beta])
                    / (2.0 * h);
                assert!(
                    (fd - jac_entry).abs() < 1e-7,
                    "J[{alpha}][{beta}] FD {fd} vs analytic {jac_entry}"
                );
            }
        }

        // Laplacian w.r.t. r_0 (sum of second differences over components).
        let laplacian = backflow.pair_laplacian(&base, 0, 1);
        for (beta, &lap_beta) in laplacian.iter().enumerate() {
            let mut total = 0.0;
            for alpha in 0..3 {
                let mut plus = base.clone();
                let mut minus = base.clone();
                plus[alpha] += h;
                minus[alpha] -= h;
                total += (displacement_at(&plus, 0, 1)[beta]
                    - 2.0 * displacement_at(&base, 0, 1)[beta]
                    + displacement_at(&minus, 0, 1)[beta])
                    / (h * h);
            }
            assert!(
                (total - lap_beta).abs() < 1e-3,
                "L[{beta}] FD {total} vs analytic {lap_beta}"
            );
        }
    }
}
