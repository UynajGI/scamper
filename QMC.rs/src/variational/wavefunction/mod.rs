//! The `WaveFunction` trait — the anti-debt hinge of the variational family.
//!
//! All wave-function physics lives in implementations of [`WaveFunction`];
//! all sampling statistics live in the kernels; all parameter adaptation
//! lives in the (L2) optimizers. No layer reaches across. L0 ships
//! translation-invariant pair Jastrow ansätze and the one-body Gaussian
//! ([`GaussianTrap`] and [`McMillanJastrow`]/[`HarmonicJastrow`]); L1 will
//! add determinants and backflow behind the same trait.
//!
//! # Configuration layout
//!
//! `DIM` is fixed at 3 for L0: every validated ansatz and Hamiltonian is
//! written for three-dimensional continuum systems (harmonic trap, McMillan
//! He-4 Jastrow), and hand-derived `∇`/`∇²` formulas keep their 3-D
//! `∇² f(r) = f'' + 2f'/r` factors only for `D = 3`. The canonical
//! [`WaveFunction::Config`] is [`Positions`]: one flat, single-allocation
//! buffer of `DIM * n_particles` doubles with particle-interleaved layout
//! `(x_0, y_0, z_0, x_1, y_1, z_1, …)`. Interleaving keeps a moved
//! particle's three coordinates cache-contiguous at offset `DIM * i`, which
//! is what the single-particle hot path touches; the SoA guarantee that
//! matters (no per-particle `Vec`, no pointer chasing) is preserved.
//! Open-box coordinates only at L0: no periodic boundaries, no minimum-image
//! convention (documented debt for the future homogeneous-bulk layer).
//!
//! # Hot-path contract
//!
//! `delta_log`, `commit_move`, `log_grad` and `log_laplacian` must not
//! allocate: kernels hand them reused [`GradBuffer`]s and stack [`Point`]s,
//! and every L0 implementation is a pure function of its arguments plus
//! owned parameters.

pub mod gaussian;
pub mod jastrow;

use super::error::VariationalError;
pub use gaussian::GaussianTrap;
pub use jastrow::{HarmonicJastrow, McMillanJastrow};

/// Spatial dimension of L0 continuum configurations.
///
/// Fixed (not a generic const) by design: the hand-derived 3-D Laplacian
/// factors are load-bearing, and nothing at L0 needs 2-D. Generalizing to
/// `const D` later is a mechanical change confined to this module family.
pub const DIM: usize = 3;

/// A single particle position (one `[f64; DIM]`).
pub type Point = [f64; DIM];

/// Owned many-particle configuration: flat, particle-interleaved positions.
///
/// The only constructor paths enforce the invariants every downstream
/// consumer relies on (length a multiple of [`DIM`], all coordinates
/// finite), so ansatz implementations may index without re-validating.
#[derive(Debug, Clone, PartialEq)]
pub struct Positions {
    coords: Vec<f64>,
}

impl Positions {
    /// Build from a flat buffer; rejects lengths that are not a multiple of
    /// [`DIM`] and non-finite coordinates.
    pub fn from_flat(coords: Vec<f64>) -> Result<Self, VariationalError> {
        if !coords.len().is_multiple_of(DIM) {
            return Err(VariationalError::invalid(
                "positions",
                format!("length {} is not a multiple of DIM = {DIM}", coords.len()),
            ));
        }
        if let Some(bad) = coords.iter().position(|x| !x.is_finite()) {
            return Err(VariationalError::invalid(
                "positions",
                format!("coordinate {bad} is non-finite"),
            ));
        }
        Ok(Self { coords })
    }

    /// Number of particles.
    #[inline]
    pub fn n_particles(&self) -> usize {
        self.coords.len() / DIM
    }

    /// Copy out particle `i`'s position.
    #[inline]
    pub fn particle(&self, index: usize) -> Point {
        read_particle(self, index)
    }

    /// Overwrite particle `i`'s position.
    #[inline]
    pub fn set_particle(&mut self, index: usize, position: Point) {
        write_particle(self, index, position);
    }

    /// The flat coordinate buffer.
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.coords
    }

    /// The flat coordinate buffer, mutably.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.coords
    }
}

impl AsRef<[f64]> for Positions {
    #[inline]
    fn as_ref(&self) -> &[f64] {
        &self.coords
    }
}

impl AsMut<[f64]> for Positions {
    #[inline]
    fn as_mut(&mut self) -> &mut [f64] {
        &mut self.coords
    }
}

/// Read particle `index` from any flat particle-interleaved configuration.
#[inline]
pub(crate) fn read_particle(cfg: &impl AsRef<[f64]>, index: usize) -> Point {
    let base = DIM * index;
    let coords = cfg.as_ref();
    [coords[base], coords[base + 1], coords[base + 2]]
}

/// Write particle `index` into any flat particle-interleaved configuration.
#[inline]
pub(crate) fn write_particle(cfg: &mut impl AsMut<[f64]>, index: usize, position: Point) {
    let base = DIM * index;
    cfg.as_mut()[base..base + DIM].copy_from_slice(&position);
}

/// Reusable `∇_i ln|ψ|` buffer: flat `DIM * n_particles` doubles.
///
/// Component `k` of particle `i` lives at index `DIM * i + k`. Owned by the
/// kernel (or test) and reused sweep after sweep so the hot path never
/// allocates.
#[derive(Debug, Clone, PartialEq)]
pub struct GradBuffer {
    data: Vec<f64>,
}

impl GradBuffer {
    /// A zeroed buffer for `n_particles` particles.
    pub fn new(n_particles: usize) -> Self {
        Self {
            data: vec![0.0; n_particles * DIM],
        }
    }

    /// Resize (and zero) for a different particle count.
    pub fn resize(&mut self, n_particles: usize) {
        self.data.clear();
        self.data.resize(n_particles * DIM, 0.0);
    }

    /// Zero all entries. [`WaveFunction::log_grad`] implementations
    /// accumulate, so callers that want the plain gradient clear first.
    pub fn clear(&mut self) {
        self.data.fill(0.0);
    }

    /// Flat view.
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Flat mutable view.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Total number of stored components (`DIM * n_particles`).
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the buffer is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// Reusable `∂ ln|ψ| / ∂ p_k` buffer, one entry per variational parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct ParamGradBuffer {
    data: Vec<f64>,
}

impl ParamGradBuffer {
    /// A zeroed buffer for `n_params` parameters.
    pub fn new(n_params: usize) -> Self {
        Self {
            data: vec![0.0; n_params],
        }
    }

    /// Resize (and zero) for a different parameter count.
    pub fn resize(&mut self, n_params: usize) {
        self.data.clear();
        self.data.resize(n_params, 0.0);
    }

    /// Zero all entries.
    pub fn clear(&mut self) {
        self.data.fill(0.0);
    }

    /// Flat view.
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Flat mutable view.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Number of stored entries.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the buffer is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// Result of an incremental single-particle evaluation:
/// `log_ratio = ln|ψ(cfg with particle moved)| − ln|ψ(cfg)|`.
///
/// The kernel multiplies by 2 for the `|ψ|²` Metropolis weight and adds it
/// to the walker's cached `ln|ψ|` on acceptance. It is a plain `Copy` value
/// at L0; the L1 determinant machinery will grow richer payloads (e.g.
/// Sherman–Morrison row updates) without changing the call sites.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeltaLog {
    /// Log wave-function change of the proposed single-particle move.
    pub log_ratio: f64,
}

/// Trial wave functions with hand-derived analytic derivatives.
///
/// The shape follows `research/vmc/DESIGN.md` §3. `delta_log` must equal the
/// full `log_psi` recompute difference at (near) machine precision — that
/// equivalence is enforced by a day-1 test and is what keeps the fast paths
/// honest.
///
/// Units: `ħ = m = 1`, so `E_L = −½ Σ_i (∇²_i ln|ψ| + |∇_i ln|ψ||²) + V`.
pub trait WaveFunction {
    /// Configuration type; the canonical L0 layout is [`Positions`]. The
    /// `AsRef`/`AsMut` bound is the crate's flat-position contract.
    type Config: AsRef<[f64]> + AsMut<[f64]>;

    /// `ln|ψ_T(cfg)|`.
    fn log_psi(&self, cfg: &Self::Config) -> f64;

    /// Per-particle gradient `∇_i ln|ψ_T|` **added** into `out`
    /// (flattened as `DIM * i + k`).
    ///
    /// Accumulate contract: implementations add their contribution without
    /// clearing, so composite ansätze ([`Product`]) evaluate factor by
    /// factor into one caller-owned buffer with no temporaries. Callers
    /// that want the plain gradient (e.g. the local-energy estimator)
    /// zero `out` first.
    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer);

    /// `Σ_i ∇_i² ln|ψ_T|(cfg)`.
    fn log_laplacian(&self, cfg: &Self::Config) -> f64;

    /// Number of variational parameters.
    fn n_params(&self) -> usize;

    /// `∂ ln|ψ_T| / ∂ p_k` **added** into `out` (same accumulate contract as
    /// [`WaveFunction::log_grad`]).
    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer);

    /// Apply additive parameter deltas `p_k += delta[k]`.
    ///
    /// Optimizer-facing (L2); L0 callers only perturb by validated finite
    /// steps. Implementations may assume `delta.len() == n_params()`.
    fn update_params(&mut self, delta: &[f64]);

    /// Incremental log-ratio of moving one particle — O(1) for one-body,
    /// O(N) for pair Jastrow terms, O(N²) for rank-1 determinant updates.
    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog;

    /// Commit a previously evaluated move into the configuration.
    ///
    /// L0 ansätze are stateless (pure functions of `cfg`), so this only
    /// writes the coordinates; L1 determinant ansätze will additionally
    /// update inverse-matrix caches here.
    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point);

    /// Periodic O(N³)-scale refresh of internal caches.
    ///
    /// A no-op for the stateless L0 ansätze: they cache nothing, so there is
    /// no floating-point drift to repair between full recomputes. The hook
    /// exists so the L1 determinant layer can rebuild inverses every K
    /// accepted moves without changing kernel call sites.
    fn rebuild(&mut self, cfg: &Self::Config);
}

/// Read/write access to an ansatz's parameter vector for checkpointing.
///
/// A small companion to [`WaveFunction`] (which, per DESIGN.md, stays free
/// of serialization concerns). Kernels use it to serialize `params` into
/// versioned snapshots and to restore them with full validation on load;
/// nothing on the hot path touches it.
pub trait WaveFunctionParams: WaveFunction {
    /// Current parameter values.
    fn param_values(&self) -> Vec<f64>;

    /// Validate and assign a full parameter vector (must be finite and
    /// individually valid; length must match [`WaveFunction::n_params`]).
    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError>;
}

/// Product ansatz `ψ_T = ψ_A · ψ_B` (log-additive).
///
/// The combinator behind every Jastrow-times-confined-state construction —
/// at L0 the He-droplet state `GaussianTrap × McMillanJastrow`; at L1 the
/// Slater–Jastrow state will be exactly this type with a determinant
/// factor. Gradients and Laplacians add because
/// `ln|ψ_A ψ_B| = ln|ψ_A| + ln|ψ_B|`; the parameter vector is the
/// concatenation `(p_A, p_B)`.
///
/// Hot path (`delta_log`, `log_grad`, `commit_move`): allocation-free —
/// both factors accumulate into the caller-owned buffers. `log_grad_params`
/// builds two small factor-sized temporaries (optimizer-facing path only).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Product<A, B> {
    left: A,
    right: B,
}

impl<A, B> Product<A, B> {
    /// Combine two factors over the same configuration type.
    pub fn new(left: A, right: B) -> Self {
        Self { left, right }
    }

    /// The first factor.
    #[inline]
    pub const fn left(&self) -> &A {
        &self.left
    }

    /// The second factor.
    #[inline]
    pub const fn right(&self) -> &B {
        &self.right
    }
}

impl<A, B> WaveFunction for Product<A, B>
where
    A: WaveFunction,
    B: WaveFunction<Config = A::Config>,
{
    type Config = A::Config;

    #[inline]
    fn log_psi(&self, cfg: &Self::Config) -> f64 {
        self.left.log_psi(cfg) + self.right.log_psi(cfg)
    }

    #[inline]
    fn log_grad(&self, cfg: &Self::Config, out: &mut GradBuffer) {
        self.left.log_grad(cfg, out);
        self.right.log_grad(cfg, out);
    }

    #[inline]
    fn log_laplacian(&self, cfg: &Self::Config) -> f64 {
        self.left.log_laplacian(cfg) + self.right.log_laplacian(cfg)
    }

    #[inline]
    fn n_params(&self) -> usize {
        self.left.n_params() + self.right.n_params()
    }

    fn log_grad_params(&self, cfg: &Self::Config, out: &mut ParamGradBuffer) {
        let n_left = self.left.n_params();
        let mut left_buffer = ParamGradBuffer::new(n_left);
        self.left.log_grad_params(cfg, &mut left_buffer);
        out.as_mut_slice()[..n_left].copy_from_slice(left_buffer.as_slice());
        let mut right_buffer = ParamGradBuffer::new(self.right.n_params());
        self.right.log_grad_params(cfg, &mut right_buffer);
        out.as_mut_slice()[n_left..].copy_from_slice(right_buffer.as_slice());
    }

    #[inline]
    fn update_params(&mut self, delta: &[f64]) {
        let (left_delta, right_delta) = delta.split_at(self.left.n_params());
        self.left.update_params(left_delta);
        self.right.update_params(right_delta);
    }

    #[inline]
    fn delta_log(&self, cfg: &Self::Config, particle: usize, new_pos: &Point) -> DeltaLog {
        DeltaLog {
            log_ratio: self.left.delta_log(cfg, particle, new_pos).log_ratio
                + self.right.delta_log(cfg, particle, new_pos).log_ratio,
        }
    }

    #[inline]
    fn commit_move(&mut self, cfg: &mut Self::Config, particle: usize, new_pos: &Point) {
        self.left.commit_move(cfg, particle, new_pos);
        self.right.commit_move(cfg, particle, new_pos);
    }

    #[inline]
    fn rebuild(&mut self, cfg: &Self::Config) {
        self.left.rebuild(cfg);
        self.right.rebuild(cfg);
    }
}

impl<A, B> WaveFunctionParams for Product<A, B>
where
    A: WaveFunctionParams,
    B: WaveFunctionParams<Config = A::Config>,
{
    fn param_values(&self) -> Vec<f64> {
        let mut values = self.left.param_values();
        values.extend(self.right.param_values());
        values
    }

    fn try_set_params(&mut self, values: &[f64]) -> Result<(), VariationalError> {
        if values.len() != self.n_params() {
            return Err(VariationalError::invalid(
                "params",
                format!("expected {} values, got {}", self.n_params(), values.len()),
            ));
        }
        let (left_values, right_values) = values.split_at(self.left.n_params());
        self.left.try_set_params(left_values)?;
        self.right.try_set_params(right_values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positions_enforces_layout_and_finiteness() {
        assert!(Positions::from_flat(vec![0.0; 6]).is_ok());
        assert_eq!(Positions::from_flat(vec![0.0; 6]).unwrap().n_particles(), 2);
        assert!(Positions::from_flat(vec![0.0; 7]).is_err());
        assert!(Positions::from_flat(vec![0.0, f64::NAN, 0.0]).is_err());
        assert!(Positions::from_flat(vec![f64::INFINITY, 0.0, 0.0]).is_err());
    }

    #[test]
    fn positions_particle_accessors_round_trip() {
        let mut positions = Positions::from_flat(vec![0.0; 9]).unwrap();
        positions.set_particle(2, [1.0, 2.0, 3.0]);
        assert_eq!(positions.particle(2), [1.0, 2.0, 3.0]);
        assert_eq!(positions.as_ref()[6..9], [1.0, 2.0, 3.0]);
    }

    #[test]
    fn buffers_start_zeroed_and_resize() {
        let mut grad = GradBuffer::new(4);
        assert_eq!(grad.len(), 4 * DIM);
        assert!(grad.as_slice().iter().all(|&x| x == 0.0));
        grad.as_mut_slice()[0] = 1.0;
        grad.resize(2);
        assert_eq!(grad.len(), 2 * DIM);
        assert!(grad.as_slice().iter().all(|&x| x == 0.0));

        let mut params = ParamGradBuffer::new(3);
        assert_eq!(params.len(), 3);
        params.resize(1);
        assert_eq!(params.len(), 1);
        assert!(params.as_slice().iter().all(|&x| x == 0.0));
    }

    #[test]
    fn product_is_log_additive_with_concatenated_parameters() {
        let gaussian = GaussianTrap::new(0.4, [0.1, 0.0, -0.1]).unwrap();
        let mcmillan = McMillanJastrow::new(1.0).unwrap();
        let product = Product::new(gaussian, mcmillan);
        let cfg = Positions::from_flat(vec![
            0.6, -0.2, 0.4, //
            -0.9, 0.3, 0.1, //
            0.2, 1.1, -0.5,
        ])
        .unwrap();

        assert_eq!(product.n_params(), 2);
        assert_eq!(product.param_values(), vec![0.4, 1.0]);
        assert!(
            (product.log_psi(&cfg) - (gaussian.log_psi(&cfg) + mcmillan.log_psi(&cfg))).abs()
                < 1e-15
        );
        assert!(
            (product.log_laplacian(&cfg)
                - (gaussian.log_laplacian(&cfg) + mcmillan.log_laplacian(&cfg)))
            .abs()
                < 1e-12
        );

        // Gradients accumulate: sum of factor gradients == product gradient.
        let mut factors = GradBuffer::new(3);
        gaussian.log_grad(&cfg, &mut factors);
        mcmillan.log_grad(&cfg, &mut factors);
        let mut combined = GradBuffer::new(3);
        product.log_grad(&cfg, &mut combined);
        for (a, b) in factors.as_slice().iter().zip(combined.as_slice()) {
            assert!((a - b).abs() < 1e-14);
        }

        // Parameter gradients concatenate.
        let mut params = ParamGradBuffer::new(2);
        product.log_grad_params(&cfg, &mut params);
        let mut expected_left = ParamGradBuffer::new(1);
        gaussian.log_grad_params(&cfg, &mut expected_left);
        assert_eq!(params.as_slice()[0], expected_left.as_slice()[0]);

        // update_params splits; try_set_params validates both halves.
        let mut updated = product;
        updated.update_params(&[0.1, -0.1]);
        assert_eq!(updated.param_values(), vec![0.5, 0.9]);
        assert!(updated.try_set_params(&[0.5, -1.0]).is_err());
        assert!(updated.try_set_params(&[0.5]).is_err());
        assert!(updated.try_set_params(&[0.5, 1.0]).is_ok());

        // delta_log adds the factor ratios.
        let new_pos = [0.7, 0.0, 0.0];
        let incremental = product.delta_log(&cfg, 1, &new_pos).log_ratio;
        let expected = gaussian.delta_log(&cfg, 1, &new_pos).log_ratio
            + mcmillan.delta_log(&cfg, 1, &new_pos).log_ratio;
        assert!((incremental - expected).abs() < 1e-15);
    }
}
