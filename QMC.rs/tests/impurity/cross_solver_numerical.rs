//! Cross-solver validation: wormhole↔occupation.
//!
//! Both solvers are compared against the exact analytic result for the
//! free two-level system (g=0): |⟨σ⟩| = tanh(βΔ/2). The wormhole and
//! occupation solvers use different basis conventions (see
//! `cross_solver.rs` for the full catalogue), so the sign of the
//! measured observable differs, but the magnitude must agree and each
//! must match the exact result individually.
//!
//! For the interacting case (g>0), the two solvers measure different
//! physical quantities because the wormhole uses a rotated basis and a
//! different coupling normalisation (g_eff = g/2 in the pre-rotation
//! Hamiltonian).  This is a fundamental limitation, not a bug.
//!
//! Instead of forcing a direct solver-to-solver comparison, both solvers
//! are validated against a shared exact-diagonalisation (ED) reference
//! computed via the scaling-and-squaring matrix exponential (same pattern
//! as `lattice_ed.rs`).

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use qmc_rs::impurity::ImpurityQmc;
use qmc_rs::OccupationWorldlineQmc;

// ─── P1.1: Wormhole ↔ Occupation (free two-level system) ────────────────

#[test]
fn occupation_and_wormhole_match_exact_free_two_level() {
    // Free two-level system: g=0, only splitting Δ.
    //
    // The two solvers use different basis conventions (see cross_solver.rs):
    //   • Wormhole (rotated basis, σz_sampled = σx_physical):
    //       H = -(Δ/2)σx  →  MagnetizationSigmaZ = +tanh(βΔ/2)
    //   • Occupation (occupation basis):
    //       H = +(Δ/2)σz  →  OccupationSigmaZ = -tanh(βΔ/2)
    //
    // The sign flip is convention difference #4. Both solvers measure the
    // SAME physical splitting Δ, so their magnitudes must agree and each
    // must match the exact result individually.
    let beta: f64 = 4.0;
    let delta: f64 = 1.0;
    let exact_tanh: f64 = (beta * delta / 2.0).tanh(); // tanh(2.0) ≈ 0.964

    // ── Occupation solver ──────────────────────────────────────────
    let mut params_occ = Params::new();
    params_occ.set("beta", beta);
    params_occ.set("kind", "rabi");
    params_occ.set("spin_splitting", delta);
    params_occ.set("g", 0.0); // free
    params_occ.set("omega0", 1.0);
    params_occ.set("cutoff", 5);
    let config_occ = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_occ = Scheduler::new(RayonBackend::new(1), config_occ)
        .run_one::<OccupationWorldlineQmc>(&params_occ);
    let sz_occ = results_occ
        .get("OccupationSigmaZ")
        .expect("OccupationSigmaZ");

    // ── Wormhole solver ───────────────────────────────────────────
    // NOTE: the wormhole uses `tunnelling` (not `spin_splitting`).
    let mut params_wh = Params::new();
    params_wh.set("beta", beta);
    params_wh.set("model", "rabi");
    params_wh.set("bath", "single");
    params_wh.set("omega0", 1.0);
    params_wh.set("g", 0.0); // free
    params_wh.set("tunnelling", delta);
    params_wh.set("h_z", 0.0);
    let config_wh = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 10000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_wh =
        Scheduler::new(RayonBackend::new(1), config_wh).run_one::<ImpurityQmc>(&params_wh);
    let sz_wh = results_wh
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ");

    // ── 1. Occupation ⟨σz⟩ vs exact ──────────────────────────────
    // Occupation convention: ⟨σz⟩ = -tanh(βΔ/2)
    let exact_occ = -exact_tanh;
    assert!(
        (sz_occ.mean - exact_occ).abs() < 4.0 * sz_occ.stderr.max(0.02),
        "Occupation ⟨σz⟩={:.4}±{:.4}, exact={:.4}",
        sz_occ.mean,
        sz_occ.stderr,
        exact_occ
    );

    // ── 2. Wormhole ⟨σz⟩ (= physical ⟨σx⟩) vs exact ──────────────
    // Wormhole convention: MagnetizationSigmaZ = +tanh(βΔ/2)
    let exact_wh = exact_tanh;
    assert!(
        (sz_wh.mean - exact_wh).abs() < 4.0 * sz_wh.stderr.max(0.02),
        "Wormhole ⟨σz⟩={:.4}±{:.4}, exact={:.4}",
        sz_wh.mean,
        sz_wh.stderr,
        exact_wh
    );

    // ── 3. Cross-solver: occupation vs wormhole ──────────────────
    // Sign-corrected comparison: -⟨σz⟩_occ should agree with ⟨σz⟩_wh
    // within combined 4σ (sqrt of sum of squared stderrs).
    let combined_sigma = (sz_occ.stderr.powi(2) + sz_wh.stderr.powi(2))
        .sqrt()
        .max(0.02);
    assert!(
        ((-sz_occ.mean) - sz_wh.mean).abs() < 4.0 * combined_sigma,
        "Cross-solver mismatch: -occ⟨σz⟩={:.4}±{:.4}, wh⟨σz⟩={:.4}±{:.4}, combined 4σ={:.4}",
        -sz_occ.mean,
        sz_occ.stderr,
        sz_wh.mean,
        sz_wh.stderr,
        4.0 * combined_sigma
    );
}

// ── Exact diagonalisation helper ─────────────────────────────────────────

/// Minimal dense matrix and matrix-exponential tools, following the
/// same scaling-and-squaring pattern as [`super::lattice_ed::DenseMatrix`].
mod dense {
    /// Row-major dense matrix.
    pub(crate) struct Matrix {
        pub(crate) dim: usize,
        pub(crate) elements: Vec<f64>,
    }

    impl Matrix {
        pub(crate) fn zero(dim: usize) -> Self {
            Self {
                dim,
                elements: vec![0.0; dim * dim],
            }
        }

        fn identity(dim: usize) -> Self {
            let mut m = Self::zero(dim);
            for i in 0..dim {
                m.elements[i * dim + i] = 1.0;
            }
            m
        }

        pub(crate) fn get(&self, i: usize, j: usize) -> f64 {
            self.elements[i * self.dim + j]
        }

        fn set(&mut self, i: usize, j: usize, val: f64) {
            self.elements[i * self.dim + j] = val;
        }

        pub(crate) fn add(&mut self, i: usize, j: usize, val: f64) {
            self.elements[i * self.dim + j] += val;
        }

        fn multiply(&self, other: &Self) -> Self {
            let dim = self.dim;
            let mut result = Self::zero(dim);
            for i in 0..dim {
                for j in 0..dim {
                    let mut sum = 0.0;
                    for k in 0..dim {
                        sum += self.get(i, k) * other.get(k, j);
                    }
                    result.set(i, j, sum);
                }
            }
            result
        }

        fn scale(&self, s: f64) -> Self {
            Self {
                dim: self.dim,
                elements: self.elements.iter().map(|&x| x * s).collect(),
            }
        }

        pub(crate) fn trace(&self) -> f64 {
            (0..self.dim).map(|i| self.get(i, i)).sum()
        }

        /// Matrix exponential exp(-beta * H) via scaling and squaring.
        pub(crate) fn expm_negative(&self, beta: f64) -> Self {
            let mut a = self.scale(-beta);
            // Scale down so max |element| ≤ 0.5
            let mut max_el = 0.0_f64;
            for &v in &a.elements {
                max_el = max_el.max(v.abs());
            }
            let mut n_scales = 0;
            while max_el > 0.5 {
                a = a.scale(0.5);
                max_el *= 0.5;
                n_scales += 1;
            }
            // Taylor series: I + A + A²/2! + A³/3! + ...
            let dim = self.dim;
            let mut result = Self::identity(dim);
            let mut term = Self::identity(dim);
            for k in 1..=40 {
                term = term.multiply(&a).scale(1.0 / k as f64);
                result
                    .elements
                    .iter_mut()
                    .zip(&term.elements)
                    .for_each(|(r, &t)| *r += t);
            }
            // Square n_scales times
            for _ in 0..n_scales {
                result = result.multiply(&result);
            }
            result
        }
    }
}

use dense::Matrix;

/// Build the single-mode Rabi Hamiltonian in the **occupation**
/// (quantum-optics) convention and return exact thermal expectations
/// computed via the matrix-exponential density matrix:
///
/// ```text
///   H_occ = ω a†a + (Δ/2) σz + g σx (a + a†)
/// ```
///
/// Basis: `|n, s⟩ → 2n + s`, `s = 0` for ↓ (`σz = −1`), `s = 1` for ↑
/// (`σz = +1`).  Hilbert-space dimension = `2 × cutoff`.
///
/// Returns `(exact_energy, exact_sigma_z, exact_sigma_x)`.  Note that
/// `exact_sigma_x` = 0 for this Hamiltonian because σx is odd under the
/// Rabi parity P = σz exp(iπ a†a) and the thermal state is even.
fn ed_occupation_rabi(
    beta: f64,
    omega: f64,
    spin_splitting: f64,
    g: f64,
    cutoff: usize,
) -> (f64, f64, f64) {
    let dim = 2 * cutoff;
    let mut h = Matrix::zero(dim);

    for n in 0..cutoff {
        for s in 0..2usize {
            let i = 2 * n + s;
            let sz = if s == 0 { -1.0 } else { 1.0 };

            h.add(i, i, omega * n as f64);
            h.add(i, i, 0.5 * spin_splitting * sz);

            // g σx (a + a†): connects |n, s⟩ ↔ |n+1, 1−s⟩
            if n + 1 < cutoff {
                let amplitude = (n as f64 + 1.0).sqrt();
                let j = 2 * (n + 1) + (1 - s);
                let val = g * amplitude;
                h.add(i, j, val);
                h.add(j, i, val);
            }
        }
    }

    let rho = h.expm_negative(beta);
    let z = rho.trace();

    // ⟨σz⟩ — diagonal in the occupation basis
    let mut sigma_z_exp = 0.0;
    for n in 0..cutoff {
        for s in 0..2usize {
            let i = 2 * n + s;
            let sz = if s == 0 { -1.0 } else { 1.0 };
            sigma_z_exp += sz * rho.get(i, i);
        }
    }
    sigma_z_exp /= z;

    // ⟨σx⟩ — off-diagonal, spin-flip at same n
    // (zero by parity for this Hamiltonian; computed for completeness)
    let mut sigma_x_exp = 0.0;
    for n in 0..cutoff {
        let i0 = 2 * n;
        let i1 = 2 * n + 1;
        sigma_x_exp += rho.get(i1, i0) + rho.get(i0, i1);
    }
    sigma_x_exp /= z;

    let mut energy = 0.0;
    for i in 0..dim {
        for j in 0..dim {
            energy += h.get(i, j) * rho.get(j, i);
        }
    }
    energy /= z;

    (energy, sigma_z_exp, sigma_x_exp)
}

/// Build the wormhole solver's **physical** (pre-rotation) Hamiltonian
/// and return exact thermal expectations:
///
/// ```text
///   H_wh = ω a†a − (Δ/2) σx + (g/2) σz (a + a†)
/// ```
///
/// where Δ = `tunnelling`.  The wormhole samples in a rotated basis
/// where its `MagnetizationSigmaZ` observable measures the physical
/// ⟨σx⟩ of this Hamiltonian.
///
/// Returns `(exact_energy, exact_sigma_z, exact_sigma_x)`.
fn ed_wormhole_rabi(
    beta: f64,
    omega: f64,
    tunnelling: f64,
    g: f64,
    cutoff: usize,
) -> (f64, f64, f64) {
    let dim = 2 * cutoff;
    let mut h = Matrix::zero(dim);

    for n in 0..cutoff {
        for s in 0..2usize {
            let i = 2 * n + s;
            let sz = if s == 0 { -0.5 } else { 0.5 };
            // Boson kinetic
            h.add(i, i, omega * n as f64);

            // −(Δ/2) σx: spin-flip at same n
            let j = 2 * n + (1 - s);
            h.add(i, j, -0.5 * tunnelling);

            // (g/2) σz (a + a†): diagonal in spin, n↔n±1
            if n + 1 < cutoff {
                let amplitude = (n as f64 + 1.0).sqrt();
                let k = 2 * (n + 1) + s; // same spin, one boson higher
                let val = g * sz * amplitude; // g/2 * (2*sz) * amplitude = g * sz * amplitude
                h.add(i, k, val);
                h.add(k, i, val);
            }
        }
    }

    let rho = h.expm_negative(beta);
    let z = rho.trace();

    // ⟨σx⟩ — spin-flip at same n
    let mut sigma_x_exp = 0.0;
    for n in 0..cutoff {
        let i0 = 2 * n;
        let i1 = 2 * n + 1;
        sigma_x_exp += rho.get(i1, i0) + rho.get(i0, i1);
    }
    sigma_x_exp /= z;

    // ⟨σz⟩ — diagonal
    let mut sigma_z_exp = 0.0;
    for n in 0..cutoff {
        for s in 0..2usize {
            let i = 2 * n + s;
            let sz = if s == 0 { -1.0 } else { 1.0 };
            sigma_z_exp += sz * rho.get(i, i);
        }
    }
    sigma_z_exp /= z;

    let mut energy = 0.0;
    for i in 0..dim {
        for j in 0..dim {
            energy += h.get(i, j) * rho.get(j, i);
        }
    }
    energy /= z;

    (energy, sigma_z_exp, sigma_x_exp)
}

// ─── P1.1b: Wormhole ↔ Occupation (interacting, ED reference) ────────────

#[test]
fn occupation_and_wormhole_both_match_ed_reference_interacting() {
    // Each solver is validated against an ED reference built for its
    // own Hamiltonian convention, using the same DenseMatrix expm
    // infrastructure.
    //
    // The occupation solver's Hamiltonian:
    //   H_occ = ω a†a + (Δ/2) σz + g σx (a + a†)
    // → OccupationSigmaZ = ⟨σz⟩_occ
    //
    // The wormhole solver's pre-rotation Hamiltonian:
    //   H_wh  = ω a†a − (Δ/2) σx + (g/2) σz (a + a†)
    // → MagnetizationSigmaZ = ⟨σx⟩_wh   (sampled σz after rotation)
    //
    // H_occ and H_wh differ by a spin rotation AND a coupling
    // normalisation (g vs g/2).  This is a fundamental convention
    // difference, not a bug — see cross_solver.rs for the full
    // catalogue.  Each solver is therefore checked against its
    // own ED reference.

    let beta = 3.0;
    let omega = 1.0;
    let spin_splitting = 1.0;
    let g = 0.3;
    let cutoff = 8;

    // ── Occupation ED reference ───────────────────────────────────
    let (ed_occ_e, ed_occ_sz, ed_occ_sx) =
        ed_occupation_rabi(beta, omega, spin_splitting, g, cutoff);
    eprintln!(
        "ED (occupation, cutoff={cutoff}): ⟨E⟩={ed_occ_e:.8}, ⟨σz⟩={ed_occ_sz:.8}, ⟨σx⟩={ed_occ_sx:.8}"
    );

    // ── Wormhole ED reference ────────────────────────────────────
    let (ed_wh_e, ed_wh_sz, ed_wh_sx) = ed_wormhole_rabi(beta, omega, spin_splitting, g, cutoff);
    eprintln!(
        "ED (wormhole, cutoff={cutoff}):   ⟨E⟩={ed_wh_e:.8}, ⟨σz⟩={ed_wh_sz:.8}, ⟨σx⟩={ed_wh_sx:.8}"
    );

    // ── Occupation solver ──────────────────────────────────────────
    let mut params_occ = Params::new();
    params_occ.set("beta", beta);
    params_occ.set("kind", "rabi");
    params_occ.set("spin_splitting", spin_splitting);
    params_occ.set("g", g);
    params_occ.set("omega0", omega);
    params_occ.set("cutoff", cutoff);
    let config_occ = RunConfig {
        thermalization_sweeps: 4000,
        measurement_sweeps: 20000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_occ = Scheduler::new(RayonBackend::new(1), config_occ)
        .run_one::<OccupationWorldlineQmc>(&params_occ);
    let sz_occ = results_occ
        .get("OccupationSigmaZ")
        .expect("OccupationSigmaZ");

    eprintln!(
        "Occupation  ⟨σz⟩ = {:.6} ± {:.6}",
        sz_occ.mean, sz_occ.stderr
    );

    let tol_occ = (4.0 * sz_occ.stderr).max(0.04);
    assert!(
        (sz_occ.mean - ed_occ_sz).abs() < tol_occ,
        "OccupationSigmaZ {:.6} ± {:.6} vs ED ⟨σz⟩_occ {:.6} (Δ={:.6}, tol={:.6})",
        sz_occ.mean,
        sz_occ.stderr,
        ed_occ_sz,
        (sz_occ.mean - ed_occ_sz).abs(),
        tol_occ,
    );

    // ── Wormhole solver ───────────────────────────────────────────
    let mut params_wh = Params::new();
    params_wh.set("beta", beta);
    params_wh.set("model", "rabi");
    params_wh.set("bath", "single");
    params_wh.set("omega0", omega);
    params_wh.set("g", g);
    params_wh.set("tunnelling", spin_splitting);
    params_wh.set("h_z", 0.0);
    let config_wh = RunConfig {
        thermalization_sweeps: 4000,
        measurement_sweeps: 20000,
        binsize: 100,
        base_seed: 42,
        ..Default::default()
    };
    let results_wh =
        Scheduler::new(RayonBackend::new(1), config_wh).run_one::<ImpurityQmc>(&params_wh);
    let sz_wh = results_wh
        .get("MagnetizationSigmaZ")
        .expect("MagnetizationSigmaZ");

    eprintln!(
        "Wormhole    MagnetizationSigmaZ = {:.6} ± {:.6}",
        sz_wh.mean, sz_wh.stderr
    );

    // MagnetizationSigmaZ = sampled ⟨σz⟩ = physical ⟨σx⟩_wh
    let tol_wh = (4.0 * sz_wh.stderr).max(0.04);
    assert!(
        (sz_wh.mean - ed_wh_sx).abs() < tol_wh,
        "MagnetizationSigmaZ {:.6} ± {:.6} vs ED ⟨σx⟩_wh {:.6} (Δ={:.6}, tol={:.6})",
        sz_wh.mean,
        sz_wh.stderr,
        ed_wh_sx,
        (sz_wh.mean - ed_wh_sx).abs(),
        tol_wh,
    );

    // Physics sanity: at g=0.3, the boson occupation should be > 0.
    let n_occ = results_occ
        .get("OccupationBosonNumber")
        .expect("OccupationBosonNumber");
    assert!(
        n_occ.mean > 0.0 && n_occ.mean.is_finite(),
        "Boson number should be positive, got {:.4}",
        n_occ.mean
    );
    let order_wh = results_wh.get("ExpansionOrder").expect("ExpansionOrder");
    assert!(
        order_wh.mean > 0.0 && order_wh.mean.is_finite(),
        "Expansion order should be positive, got {:.4}",
        order_wh.mean
    );
}
