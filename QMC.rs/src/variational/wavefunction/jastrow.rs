//! Translation-invariant pair Jastrow trial states.

use super::super::error::VariationalError;
use super::{
    read_particle, write_particle, DeltaLog, GradBuffer, ParamGradBuffer, Point, Positions,
    WaveFunction, WaveFunctionParams, DIM,
};

/// Squared distance between two flat-config particles.
#[inline]
fn displacement(cfg: &impl AsRef<[f64]>, i: usize, j: usize) -> Point {
    let a = read_particle(cfg, i);
    let b = read_particle(cfg, j);
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// McMillan repulsive Jastrow `ψ_T = Π_{i<j} exp(−½ (b/r_ij)⁵)`.
///
/// The classic liquid/solid He-4 correlation factor (McMillan 1965): a
/// short-range repulsive core punch-out whose exponent `n = 5` matches the
/// `r⁻⁶` tail of the Lennard-Jones pair.
///
/// Analytic derivatives with `u(r) = −½ (b/r)⁵`, `u'(r) = (5/2) b⁵ r⁻⁶`,
/// `u''(r) = −15 b⁵ r⁻⁷` and the 3-D identity
/// `∇² u(r) = u''(r) + 2 u'(r)/r` (hand-derived, verified against finite
/// differences in CI):
///
/// - `ln|ψ| = Σ_{i<j} u(r_ij)`
/// - `∇_i ln|ψ| = Σ_{j≠i} (5/2) b⁵ r_ij⁻⁷ (r_i − r_j)`
/// - `∇_i² ln|ψ| = Σ_{j≠i} (u'' + 2u'/r)(r_ij) = Σ_{j≠i} (−10 b⁵ r_ij⁻⁷)`
/// - `∂ ln|ψ|/∂b = −(5/2) b⁴ Σ_{i<j} r_ij⁻⁵`
///
/// A proposal that lands exactly on another particle (`r = 0`) produces
/// `log_ratio = −∞`, which the Metropolis test rejects — the short-distance
/// divergence is self-regularizing, no special casing needed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct McMillanJastrow {
    b: f64,
}

impl McMillanJastrow {
    /// Construct with the repulsive length `b` (finite, > 0).
    pub fn new(b: f64) -> Result<Self, VariationalError> {
        VariationalError::require_positive("b", b)?;
        Ok(Self { b })
    }

    /// Repulsive length `b`.
    #[inline]
    pub const fn b(&self) -> f64 {
        self.b
    }

    /// `u(r) = −½ (b/r)⁵`.
    #[inline]
    fn pair_log(&self, r: f64) -> f64 {
        -0.5 * (self.b / r).powi(5)
    }
}

impl WaveFunction for McMillanJastrow {
    type Config = Positions;

    fn log_psi(&self, cfg: &Self::Config) -> f64 {
        let n = cfg.as_ref().len() / DIM;
        let mut sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                sum += self.pair_log(distance(cfg, i, j));
            }
        }
        sum
    }

    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer) {
        let n = cfg.as_ref().len() / DIM;
        let coefficient = 2.5 * self.b.powi(5);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = distance(cfg, i, j);
                let scale = coefficient / r.powi(7);
                let d = displacement(cfg, i, j);
                let base = DIM * i;
                for (slot, &dk) in out.as_mut_slice()[base..base + DIM]
                    .iter_mut()
                    .zip(d.iter())
                {
                    *slot += scale * dk;
                }
            }
        }
    }

    fn log_laplacian(&self, cfg: &Self::Config) -> f64 {
        let n = cfg.as_ref().len() / DIM;
        // Each pair contributes -10 b^5 / r^7 to both endpoints:
        // sum_i sum_{j != i} (-10 b^5 r_ij^-7) = -20 b^5 sum_{i<j} r_ij^-7.
        let mut sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                sum += 1.0 / distance(cfg, i, j).powi(7);
            }
        }
        -20.0 * self.b.powi(5) * sum
    }

    #[inline]
    fn n_params(&self) -> usize {
        1
    }

    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer) {
        let n = cfg.as_ref().len() / DIM;
        let mut sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                sum += 1.0 / distance(cfg, i, j).powi(5);
            }
        }
        out.as_mut_slice()[0] += -2.5 * self.b.powi(4) * sum;
    }

    #[inline]
    fn update_params(&mut self, delta: &[f64]) {
        self.b += delta[0];
    }

    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog {
        let n = cfg.as_ref().len() / DIM;
        let mut log_ratio = 0.0;
        for j in 0..n {
            if j == particle {
                continue;
            }
            let other = read_particle(cfg, j);
            let old_r = point_distance(&read_particle(cfg, particle), &other);
            let new_r = point_distance(new_pos, &other);
            log_ratio += self.pair_log(new_r) - self.pair_log(old_r);
        }
        DeltaLog { log_ratio }
    }

    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point) {
        write_particle(cfg, particle, *new_pos);
    }

    #[inline]
    fn rebuild(&mut self, _cfg: &Self::Config) {
        // Stateless ansatz: pure function of the configuration, nothing cached.
    }
}

impl WaveFunctionParams for McMillanJastrow {
    fn param_values(&self) -> Vec<f64> {
        vec![self.b]
    }

    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError> {
        if values.len() != self.n_params() {
            return Err(VariationalError::invalid(
                "params",
                format!("expected {} values, got {}", self.n_params(), values.len()),
            ));
        }
        VariationalError::require_positive("b", values[0])?;
        self.b = values[0];
        Ok(())
    }
}

/// Harmonic pair Jastrow `ψ_T = Π_{i<j} exp(−a r_ij²)`.
///
/// The L0 exact-eigenstate reference: it is the nodeless (hence ground-state,
/// by Perron–Frobenius) exact eigenfunction of the pair-harmonic Hamiltonian
/// `H = −½ Σ_i ∇_i² + Σ_{i<j} ½ k r_ij²` when `k = 4 a² N`, with total
/// energy `E₀ = 3aN(N−1)`. Derivation (hand-derived, CI-enforced):
///
/// - `Q = Σ_{i<j} r_ij² = N Σ_i |r_i − R|²` with `R` the center of mass,
///   because `Σ_{i<j} |r_i − r_j|² = N Σ_i |r_i|² − |Σ_i r_i|²`.
/// - `∇_i ln|ψ| = −2a Σ_{j≠i} (r_i − r_j) = −2aN (r_i − R)`, so
///   `Σ_i |∇_i ln|ψ||² = 4a²N² · Q/N = 4a² N Q`.
/// - `∇_i² ln|ψ| = −2a Σ_{j≠i} ∇_i·(r_i − r_j) = −6a (N−1)`, so
///   `Σ_i ∇_i² ln|ψ| = −6aN(N−1)`.
/// - Kinetic `T = −½(∇²ln|ψ| + |∇ln|ψ||²) = 3aN(N−1) − 2a²NQ`.
/// - `E_L = 3aN(N−1) + (½k − 2a²N) Q`, constant `≡ 3aN(N−1)` iff `k = 4a²N`.
///
/// `∂ ln|ψ|/∂a = −Q`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicJastrow {
    a: f64,
}

impl HarmonicJastrow {
    /// Construct with the pair width `a` (finite, > 0).
    pub fn new(a: f64) -> Result<Self, VariationalError> {
        VariationalError::require_positive("a", a)?;
        Ok(Self { a })
    }

    /// Pair width `a`.
    #[inline]
    pub const fn a(&self) -> f64 {
        self.a
    }

    /// The pair-harmonic spring constant `k = 4a²N` making this ansatz the
    /// exact ground state of `−½Σ∇² + Σ_{i<j} ½k r_ij²` for `N` particles
    /// (see the type-level derivation).
    #[inline]
    pub fn exact_pair_spring_constant(&self, n_particles: usize) -> f64 {
        4.0 * self.a * self.a * n_particles as f64
    }

    /// The exact ground-state energy `3aN(N−1)` under
    /// [`Self::exact_pair_spring_constant`].
    #[inline]
    pub fn exact_ground_state_energy(&self, n_particles: usize) -> f64 {
        3.0 * self.a * n_particles as f64 * (n_particles - 1) as f64
    }
}

impl WaveFunction for HarmonicJastrow {
    type Config = Positions;

    fn log_psi(&self, cfg: &Self::Config) -> f64 {
        let n = cfg.as_ref().len() / DIM;
        let mut sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = displacement(cfg, i, j);
                sum += d.iter().map(|&x| x * x).sum::<f64>();
            }
        }
        -self.a * sum
    }

    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer) {
        let n = cfg.as_ref().len() / DIM;
        for i in 0..n {
            let base = DIM * i;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d = displacement(cfg, i, j);
                for (slot, &dk) in out.as_mut_slice()[base..base + DIM]
                    .iter_mut()
                    .zip(d.iter())
                {
                    *slot += -2.0 * self.a * dk;
                }
            }
        }
    }

    #[inline]
    fn log_laplacian(&self, cfg: &Self::Config) -> f64 {
        let n = (cfg.as_ref().len() / DIM) as f64;
        -6.0 * self.a * n * (n - 1.0)
    }

    #[inline]
    fn n_params(&self) -> usize {
        1
    }

    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer) {
        // dln|psi|/da of -a*Q is -Q; log_psi = -a*Q so log_psi/a = -Q.
        out.as_mut_slice()[0] += self.log_psi(cfg) / self.a;
    }

    #[inline]
    fn update_params(&mut self, delta: &[f64]) {
        self.a += delta[0];
    }

    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog {
        let n = cfg.as_ref().len() / DIM;
        let mut old_q = 0.0;
        let mut new_q = 0.0;
        for j in 0..n {
            if j == particle {
                continue;
            }
            let other = read_particle(cfg, j);
            let old_d = read_particle(cfg, particle);
            for k in 0..DIM {
                old_q += (old_d[k] - other[k]) * (old_d[k] - other[k]);
                new_q += (new_pos[k] - other[k]) * (new_pos[k] - other[k]);
            }
        }
        DeltaLog {
            log_ratio: -self.a * (new_q - old_q),
        }
    }

    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point) {
        write_particle(cfg, particle, *new_pos);
    }

    #[inline]
    fn rebuild(&mut self, _cfg: &Self::Config) {
        // Stateless ansatz: pure function of the configuration, nothing cached.
    }
}

impl WaveFunctionParams for HarmonicJastrow {
    fn param_values(&self) -> Vec<f64> {
        vec![self.a]
    }

    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError> {
        if values.len() != self.n_params() {
            return Err(VariationalError::invalid(
                "params",
                format!("expected {} values, got {}", self.n_params(), values.len()),
            ));
        }
        VariationalError::require_positive("a", values[0])?;
        self.a = values[0];
        Ok(())
    }
}

#[inline]
fn distance(cfg: &impl AsRef<[f64]>, i: usize, j: usize) -> f64 {
    point_distance(&read_particle(cfg, i), &read_particle(cfg, j))
}

#[inline]
fn point_distance(a: &Point, b: &Point) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> Positions {
        Positions::from_flat(vec![
            0.5, -0.2, 0.3, //
            1.1, 0.4, -0.4, //
            -0.6, 0.9, 0.2,
        ])
        .unwrap()
    }

    #[test]
    fn constructors_reject_invalid_parameters() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(McMillanJastrow::new(bad).is_err());
            assert!(HarmonicJastrow::new(bad).is_err());
        }
        assert!(McMillanJastrow::new(1.2).is_ok());
        assert!(HarmonicJastrow::new(0.4).is_ok());
    }

    #[test]
    fn mcmillan_closed_forms() {
        let wave_function = McMillanJastrow::new(1.1).unwrap();
        let cfg = config();
        let r12 = distance(&cfg, 0, 1);
        let r13 = distance(&cfg, 0, 2);
        let r23 = distance(&cfg, 1, 2);
        let expected = -0.5 * ((1.1 / r12).powi(5) + (1.1 / r13).powi(5) + (1.1 / r23).powi(5));
        assert!((wave_function.log_psi(&cfg) - expected).abs() < 1e-15);

        let expected_laplacian =
            -20.0 * 1.1f64.powi(5) * (r12.powi(-7) + r13.powi(-7) + r23.powi(-7));
        assert!((wave_function.log_laplacian(&cfg) - expected_laplacian).abs() < 1e-12);

        let mut params = ParamGradBuffer::new(1);
        wave_function.log_grad_params(&cfg, &mut params);
        let expected_grad_b = -2.5 * 1.1f64.powi(4) * (r12.powi(-5) + r13.powi(-5) + r23.powi(-5));
        assert!((params.as_slice()[0] - expected_grad_b).abs() < 1e-12);

        let mut restored = wave_function;
        restored.try_set_params(&[1.3]).unwrap();
        assert_eq!(restored.param_values(), vec![1.3]);
        assert!(restored.try_set_params(&[0.0]).is_err());
        assert!(restored.try_set_params(&[]).is_err());
    }

    #[test]
    fn harmonic_jastrow_closed_forms_and_exact_constants() {
        let wave_function = HarmonicJastrow::new(0.45).unwrap();
        let cfg = config();
        let q = distance(&cfg, 0, 1).powi(2)
            + distance(&cfg, 0, 2).powi(2)
            + distance(&cfg, 1, 2).powi(2);
        assert!((wave_function.log_psi(&cfg) - (-0.45 * q)).abs() < 1e-14);
        assert!((wave_function.log_laplacian(&cfg) - (-6.0 * 0.45 * 3.0 * 2.0)).abs() < 1e-12);

        let mut params = ParamGradBuffer::new(1);
        wave_function.log_grad_params(&cfg, &mut params);
        assert!((params.as_slice()[0] - (-q)).abs() < 1e-12);

        assert_eq!(
            wave_function.exact_pair_spring_constant(6),
            4.0 * 0.45 * 0.45 * 6.0
        );
        assert_eq!(
            wave_function.exact_ground_state_energy(6),
            3.0 * 0.45 * 6.0 * 5.0
        );
    }
}
