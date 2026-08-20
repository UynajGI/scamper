//! One-body Gaussian trial state for a harmonic trap.

use super::super::error::VariationalError;
use super::{
    read_particle, write_particle, DeltaLog, GradBuffer, ParamGradBuffer, Point, Positions,
    WaveFunction, WaveFunctionParams,
};

/// Gaussian trap ansatz `ψ_T = Π_i exp(−α |r_i − r₀|²)`.
///
/// The exact ground state of the harmonic trap
/// `V = ½ ω² Σ_i |r_i − r₀|²` (units ħ = m = 1) at `α = ω/2`, with total
/// energy `E₀ = 3Nω/2`. This makes it the L0 zero-variance reference: at the
/// exact parameter the local energy is the constant `3Nα·1 + 0` (see the
/// derivation below), so every sampled configuration reproduces `E₀` to
/// machine precision and the Metropolis pipeline itself is validated.
///
/// Analytic derivatives (hand-derived, `d = r_i − r₀`):
///
/// - `ln|ψ| = −α Σ_i |d|²`
/// - `∇_i ln|ψ| = −2α d`
/// - `∇_i² ln|ψ| = −6α` (three dimensions), so `Σ_i = −6αN`
/// - `∂ ln|ψ|/∂α = −Σ_i |d|²`
///
/// Local energy against the trap:
/// `E_L = 3Nα + Σ_i (½ω² − 2α²) |d|²`, constant `= 3Nω/2` at `α = ω/2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianTrap {
    alpha: f64,
    center: Point,
}

impl GaussianTrap {
    /// Construct with variational width `alpha` (finite, > 0) and trap
    /// center `r₀`.
    pub fn new(alpha: f64, center: Point) -> Result<Self, VariationalError> {
        VariationalError::require_positive("alpha", alpha)?;
        if !center.iter().all(|x| x.is_finite()) {
            return Err(VariationalError::invalid(
                "center",
                "trap center must be finite in every coordinate",
            ));
        }
        Ok(Self { alpha, center })
    }

    /// Variational width `α`.
    #[inline]
    pub const fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Trap center `r₀`.
    #[inline]
    pub const fn center(&self) -> Point {
        self.center
    }
}

impl WaveFunction for GaussianTrap {
    type Config = Positions;

    #[inline]
    fn log_psi(&self, cfg: &Self::Config) -> f64 {
        let mut sum = 0.0;
        let coords = cfg.as_ref();
        for index in 0..coords.len() / super::DIM {
            let d = read_particle(cfg, index);
            sum += squared_distance(&d, &self.center);
        }
        -self.alpha * sum
    }

    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer) {
        let coords = cfg.as_ref();
        for index in 0..coords.len() / super::DIM {
            let d = read_particle(cfg, index);
            let base = super::DIM * index;
            for (k, slot) in out.as_mut_slice()[base..base + super::DIM]
                .iter_mut()
                .enumerate()
            {
                *slot += -2.0 * self.alpha * (d[k] - self.center[k]);
            }
        }
    }

    #[inline]
    fn log_laplacian(&self, cfg: &Self::Config) -> f64 {
        -6.0 * self.alpha * (cfg.as_ref().len() / super::DIM) as f64
    }

    #[inline]
    fn n_params(&self) -> usize {
        1
    }

    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer) {
        let mut sum = 0.0;
        let coords = cfg.as_ref();
        for index in 0..coords.len() / super::DIM {
            let d = read_particle(cfg, index);
            sum += squared_distance(&d, &self.center);
        }
        out.as_mut_slice()[0] += -sum;
    }

    #[inline]
    fn update_params(&mut self, delta: &[f64]) {
        self.alpha += delta[0];
    }

    #[inline]
    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog {
        let old = read_particle(cfg, particle);
        DeltaLog {
            log_ratio: -self.alpha
                * (squared_distance(new_pos, &self.center) - squared_distance(&old, &self.center)),
        }
    }

    #[inline]
    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point) {
        write_particle(cfg, particle, *new_pos);
    }

    #[inline]
    fn rebuild(&mut self, _cfg: &Self::Config) {
        // Stateless ansatz: pure function of the configuration, nothing cached.
    }
}

impl WaveFunctionParams for GaussianTrap {
    fn param_values(&self) -> Vec<f64> {
        vec![self.alpha]
    }

    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError> {
        if values.len() != self.n_params() {
            return Err(VariationalError::invalid(
                "params",
                format!("expected {} values, got {}", self.n_params(), values.len()),
            ));
        }
        VariationalError::require_positive("alpha", values[0])?;
        self.alpha = values[0];
        Ok(())
    }
}

#[inline]
fn squared_distance(a: &Point, b: &Point) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x - y) * (x - y))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> Positions {
        Positions::from_flat(vec![
            0.1, -0.2, 0.3, //
            1.0, 0.5, -0.4, //
            -0.7, 0.2, 0.9,
        ])
        .unwrap()
    }

    #[test]
    fn constructor_rejects_invalid_parameters() {
        assert!(GaussianTrap::new(0.5, [0.0; 3]).is_ok());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(GaussianTrap::new(bad, [0.0; 3]).is_err());
        }
        assert!(GaussianTrap::new(0.5, [f64::NAN, 0.0, 0.0]).is_err());
    }

    #[test]
    fn closed_form_log_and_derivatives() {
        let wave_function = GaussianTrap::new(0.7, [0.3, -0.2, 0.1]).unwrap();
        let cfg = config();
        // Sum of squared center displacements: 0.08 + 1.23 + 1.80.
        assert!((wave_function.log_psi(&cfg) - (-0.7 * 3.11)).abs() < 1e-15);

        let mut grad = GradBuffer::new(3);
        wave_function.log_grad(&cfg, &mut grad);
        assert_eq!(grad.as_slice()[0], -2.0 * 0.7 * (0.1 - 0.3));
        assert_eq!(wave_function.log_laplacian(&cfg), -6.0 * 0.7 * 3.0);

        let mut params = ParamGradBuffer::new(1);
        wave_function.log_grad_params(&cfg, &mut params);
        assert!((params.as_slice()[0] - (-3.11)).abs() < 1e-14);

        let mut restored = wave_function;
        restored.try_set_params(&[0.9]).unwrap();
        assert_eq!(restored.param_values(), vec![0.9]);
        restored.update_params(&[0.1]);
        assert_eq!(restored.alpha(), 1.0);
        assert!(restored.try_set_params(&[0.5, 0.5]).is_err());
        assert!(restored.try_set_params(&[f64::NAN]).is_err());
    }
}
