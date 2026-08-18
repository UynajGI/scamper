//! Cross-solver validation of local Metropolis for continuous spins.
//!
//! Gap closed: "Continuous spins (XY/Heisenberg) via Metropolis, frustrated
//! systems" — previously the continuous-spin Metropolis path (plane-rotation
//! proposals through [`cmc_rs::OPSSStrategy`]) was never cross-checked against
//! an independently validated solver, and no frustrated-geometry continuous
//! test existed.
//!
//! ## Part 1 (mandatory): XY/O(2) and Heisenberg/O(3) Metropolis vs Wolff
//!
//! Both solvers run the same 8×8 square lattice, PBC, J=1, at two
//! temperatures per symmetry (one test each, mirroring the runtime budget of
//! the default suite) chosen where local updates mix well:
//!   - XY: β = 0.70 and β = 0.90, both below the Kosterlitz–Thouless point
//!     (β_KT ≈ 1.1199 for J=1); verified empirically to keep local-update
//!     autocorrelation inside the 100-sweep bin size at L=8.
//!   - Heisenberg: β = 0.50 and β = 0.90 (2D O(3) has T_c = 0, so every β
//!     mixes rapidly — the upper value just probes the strongly correlated
//!     regime).
//!
//! Default statistics per test: 8 independent seeds per solver (raised by
//! `SCUTTLE_ZSCORE_SEEDS` for nightly monitoring), 3000 thermalization +
//! 2000 measurement sweeps, bin size 100. All four cross-solver pooled-z
//! values measured at these settings stay below |z| = 2.
//!
//! Wolff for O(N) is independently validated against exact Langevin (O(3))
//! and Bessel-ratio (O(2)) results on the 2-site ring in
//! `tests/physics/usage_exact.rs`, so agreement here transfers that
//! credibility to the Metropolis-continuous kernel.
//!
//! Statistic (repo cross-solver convention): per seed the Carlo.rs scheduler
//! returns binned (mean, stderr) for ⟨E⟩ and ⟨m²⟩. Each solver pools its n
//! independent seeds as mean = average of seed means and SEM = √(Σ se_i²)/n;
//! the cross-solver pooled-z is
//!   z = (mean_metro − mean_wolff) / √(SEM²_metro + SEM²_wolff),
//! required |z| < 4 per observable per temperature.
//!
//! Time-homogeneity note: [`cmc_rs::OPSSStrategy`] adapts its rotation width
//! σ only while `SimulationPhase::allows_adaptation()` — the Carlo.rs
//! scheduler freezes it for the entire measurement phase, so every sample in
//! the comparison is drawn from a fixed (time-homogeneous) transition kernel.
//! The dedicated acceptance-rate checks below additionally pin σ constant
//! from sweep one by driving the kernel manually through `sweep()` (which
//! defaults to the frozen measurement phase).
//!
//! ## Part 2 (frustrated): antiferromagnetic XY triangle vs exact quadrature
//!
//! No closed form exists for frustrated continuous-spin systems in general,
//! but the 3-site AFM XY triangle (J = −1, `build_chain(3, true)`) has a
//! 3-degree-of-freedom partition function whose global-rotation zero mode can
//! be factored out exactly, leaving a smooth, strictly periodic 2D integral.
//! Periodic trapezoidal quadrature converges on such integrands with spectral
//! accuracy, giving an effectively exact ⟨E⟩(β) and chirality moment ⟨κ²⟩(β)
//! at any β. Metropolis (fixed σ, multi-seed) is then checked against this
//! reference with per-seed exact-value z-scores — a full equilibrium check,
//! strictly stronger than the T→0 directional limit (E → E₀ = −3J/2 at the
//! 120° states, κ → ±1). Metropolis is sign-agnostic about J, so this also
//! validates the frustrated case; Wolff is not applicable (cluster updates
//! assert J ≥ 0).

use super::common::zscore_seed_count;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    build_chain, build_square, metropolis_hastings_step, Algorithm, ClassicalMC, EnergyPatch,
    Hamiltonian, Initializable, MetropolisCore, MetropolisHastingsAcceptance, ONModel,
    OPSSStrategy, ProposalStrategy, ProposedMove, SiteSpinMove, System, WolffCore,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Binned point estimate: (mean, stderr).
type Estimate = (f64, f64);
/// One solver run's observables: ((E, SE_E), (m², SE_m²)).
type SolverRun = (Estimate, Estimate);

const L: usize = 8;
const J_FM: f64 = 1.0;
const N_SEEDS: usize = 8;
const THERM_SWEEPS: u64 = 3_000;
const MEAS_SWEEPS: u64 = 2_000;
const BINSIZE: usize = 100;

/// Fixed Metropolis rotation width used by the manual (non-adaptive) runs.
const SIGMA: f64 = 0.5;

// ── Scheduler-driven runs (mirror zscore_validation.rs / zscore_extended.rs) ──

/// One scheduler run on the 8×8 FM O(D) model: returns ((E, SE_E), (m², SE_m²)).
fn run_8x8<const D: usize, A>(beta: f64, seed: u64) -> SolverRun
where
    ClassicalMC<ONModel<D>, A>: carlo_rs::FromParams,
    A: Algorithm<ONModel<D>> + Default,
{
    let mut params = Params::new();
    params.set("Lx", L);
    params.set("Ly", L);
    params.set("J", J_FM);
    params.set("beta", beta);
    let config = RunConfig {
        thermalization_sweeps: THERM_SWEEPS,
        measurement_sweeps: MEAS_SWEEPS,
        binsize: BINSIZE,
        base_seed: seed,
        ..Default::default()
    };
    let results =
        Scheduler::new(RayonBackend::new(1), config).run_one::<ClassicalMC<ONModel<D>, A>>(&params);
    let e = results.get("Energy").expect("Energy observable");
    let m2 = results.get("M2").expect("M2 observable");
    ((e.mean, e.stderr), (m2.mean, m2.stderr))
}

/// Pool independent-seed (mean, stderr) pairs: (pooled mean, SEM of the pool).
fn pooled(results: &[(f64, f64)]) -> (f64, f64) {
    let n = results.len() as f64;
    let mean = results.iter().map(|(m, _)| m).sum::<f64>() / n;
    let sem = results.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
    (mean, sem)
}

/// Cross-solver pooled-z agreement on both ⟨E⟩ and ⟨m²⟩ (|z| < 4 each).
fn assert_cross_solver_agreement(metro: &[SolverRun], wolff: &[SolverRun], label: &str) {
    let (e_m, e_m_se) = pooled(&metro.iter().map(|(e, _)| *e).collect::<Vec<_>>());
    let (e_w, e_w_se) = pooled(&wolff.iter().map(|(e, _)| *e).collect::<Vec<_>>());
    let z_e = (e_m - e_w) / (e_m_se * e_m_se + e_w_se * e_w_se).sqrt();

    let (m_m, m_m_se) = pooled(&metro.iter().map(|(_, m)| *m).collect::<Vec<_>>());
    let (m_w, m_w_se) = pooled(&wolff.iter().map(|(_, m)| *m).collect::<Vec<_>>());
    let z_m2 = (m_m - m_w) / (m_m_se * m_m_se + m_w_se * m_w_se).sqrt();

    eprintln!(
        "[cross-solver {label}] E: metro {e_m:.4}±{e_m_se:.4} wolff {e_w:.4}±{e_w_se:.4} z={z_e:.2} | \
         m²: metro {m_m:.4}±{m_m_se:.4} wolff {m_w:.4}±{m_w_se:.4} z={z_m2:.2}"
    );
    assert!(
        z_e.abs() < 4.0,
        "{label}: ⟨E⟩ pooled-z = {z_e:.2} (metro {e_m:.4}±{e_m_se:.4}, wolff {e_w:.4}±{e_w_se:.4})"
    );
    assert!(
        z_m2.abs() < 4.0,
        "{label}: ⟨m²⟩ pooled-z = {z_m2:.2} (metro {m_m:.4}±{m_m_se:.4}, wolff {m_w:.4}±{m_w_se:.4})"
    );
}

fn cross_solver_at<const D: usize>(beta: f64, label: &str) {
    let n_seeds = zscore_seed_count(N_SEEDS);
    let metro: Vec<SolverRun> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, MetropolisCore<OPSSStrategy>>(beta, seed))
        .collect();
    let wolff: Vec<SolverRun> = (0..n_seeds as u64)
        .map(|seed| run_8x8::<D, WolffCore>(beta, seed))
        .collect();
    assert_cross_solver_agreement(&metro, &wolff, label);
}

#[test]
fn xy_metropolis_matches_wolff_8x8_beta_070() {
    cross_solver_at::<2>(0.70, "XY β=0.70");
}

#[test]
fn xy_metropolis_matches_wolff_8x8_beta_090() {
    cross_solver_at::<2>(0.90, "XY β=0.90");
}

#[test]
fn heisenberg_metropolis_matches_wolff_8x8_beta_050() {
    cross_solver_at::<3>(0.50, "O(3) β=0.50");
}

#[test]
fn heisenberg_metropolis_matches_wolff_8x8_beta_090() {
    cross_solver_at::<3>(0.90, "O(3) β=0.90");
}

// ── Acceptance rates of the pinned-width Metropolis-continuous kernel ─────────

/// Acceptance rate of fixed-σ plane-rotation Metropolis on the 8×8 O(D) model.
///
/// `OPSSStrategy` resets its internal counters at every `finish_sweep`, so the
/// rate cannot be read off a driven kernel. Instead this replays the exact
/// `MetropolisCore::sweep_with_phase` body — strategy proposal, then
/// `metropolis_hastings_step` — and counts `TrialOutcome::accepted` directly.
/// σ is pinned by never invoking the adaptive path.
fn metropolis_acceptance<const D: usize>(beta: f64, seed: u64) -> f64
where
    MetropolisCore<OPSSStrategy>: Algorithm<ONModel<D>>,
    OPSSStrategy: ProposalStrategy<ONModel<D>>,
{
    let lattice = build_square(L, L, true);
    let model = ONModel::<D>::new(J_FM);
    let mut system = System::new(lattice, D, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, D).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut kernel = MetropolisCore::with_strategy(OPSSStrategy::new().with_sigma(SIGMA));
    for _ in 0..200 {
        kernel.sweep(&mut system, &model, &mut rng);
    }

    const RATE_SWEEPS: u64 = 2_000;
    let acceptance = MetropolisHastingsAcceptance;
    let ensemble = system.canonical_ensemble();
    let mut patch = EnergyPatch::default();
    let mut accepted: u64 = 0;
    let mut attempted: u64 = 0;
    for _ in 0..RATE_SWEEPS {
        for site in 0..system.n_sites() {
            let proposal = kernel.strategy.propose(&model, &system, site, &mut rng);
            let movement = SiteSpinMove::new(site, proposal.spin);
            let proposal = ProposedMove::new(movement, proposal.log_reverse_over_forward);
            let outcome = metropolis_hastings_step(
                &mut system,
                &model,
                &proposal,
                &ensemble,
                &acceptance,
                &mut patch,
                &mut rng,
            );
            attempted += 1;
            accepted += u64::from(outcome.accepted);
        }
    }
    accepted as f64 / attempted as f64
}

#[test]
fn continuous_metropolis_acceptance_stays_in_working_window() {
    // With σ = 0.5 rad the fixed-width kernel must stay in the efficient
    // window across both symmetries and all validated temperatures; too-small
    // acceptance would signal a random-walk-dominated (unvalidated-mixing)
    // regime, too-large acceptance a nearly-no-op proposal.
    let cases = [(2usize, 0.70_f64), (2, 0.90), (3, 0.50), (3, 0.90)];
    for (d, beta) in cases {
        let rate = match d {
            2 => metropolis_acceptance::<2>(beta, 0xACCE),
            _ => metropolis_acceptance::<3>(beta, 0xACCE),
        };
        eprintln!("[acceptance] O({d}) β={beta:.2} σ={SIGMA}: {rate:.3}");
        assert!(
            (0.25..0.98).contains(&rate),
            "O({d}) β={beta:.2}: acceptance {rate:.3} outside [0.25, 0.98]"
        );
    }
}

// ── Frustrated AFM XY triangle: exact quadrature reference ────────────────────

/// Triangle chirality κ = (2/(3√3)) Σ sin(θᵢ − θⱼ) over the directed edges.
/// κ = ±1 on the two 120° ground-state chiral sectors, 0 when disordered.
fn triangle_chirality(spins: &[f64]) -> f64 {
    let angle = |s: &[f64]| s[1].atan2(s[0]);
    let t1 = angle(&spins[0..2]);
    let t2 = angle(&spins[2..4]);
    let t3 = angle(&spins[4..6]);
    (2.0 / (3.0 * 3.0f64.sqrt())) * ((t1 - t2).sin() + (t2 - t3).sin() + (t3 - t1).sin())
}

/// Exact ⟨E⟩ and ⟨κ²⟩ of the AFM XY triangle (J = −1) at inverse temperature
/// β via 2D periodic trapezoidal quadrature.
///
/// With θ₃ = 0 (global-rotation zero mode factored out) the energy is
/// E(θ₁,θ₂) = cos(θ₁−θ₂) + cos θ₂ + cos θ₁. The integrand is smooth and
/// periodic, so the trapezoidal rule converges spectrally; n = 1024 per
/// dimension is orders of magnitude beyond what the β ≤ 3 features require.
fn afm_xy_triangle_exact(beta: f64) -> (f64, f64) {
    const N: usize = 1024;
    let mut w_sum = 0.0_f64;
    let mut we_sum = 0.0_f64;
    let mut wk2_sum = 0.0_f64;
    for i in 0..N {
        let t1 = std::f64::consts::TAU * i as f64 / N as f64;
        for j in 0..N {
            let t2 = std::f64::consts::TAU * j as f64 / N as f64;
            let energy = (t1 - t2).cos() + t2.cos() + t1.cos();
            let weight = (-beta * energy).exp();
            let kappa = (2.0 / (3.0 * 3.0f64.sqrt())) * ((t1 - t2).sin() + t2.sin() + (-t1).sin());
            w_sum += weight;
            we_sum += weight * energy;
            wk2_sum += weight * kappa * kappa;
        }
    }
    (we_sum / w_sum, wk2_sum / w_sum)
}

/// Binned (mean, stderr) of a sample series.
fn binned_stats(samples: &[f64], binsize: usize) -> (f64, f64) {
    assert!(
        samples.len().is_multiple_of(binsize),
        "sample count must bin exactly"
    );
    let bins: Vec<f64> = samples
        .chunks(binsize)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect();
    let n = bins.len() as f64;
    let mean = bins.iter().sum::<f64>() / n;
    let var = bins.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, (var / n).sqrt())
}

/// Fixed-σ Metropolis run on the 3-site AFM XY triangle:
/// returns ((E, SE_E), (κ², SE_κ²)).
fn run_afm_triangle(beta: f64, seed: u64) -> SolverRun {
    const THERM: usize = 2_000;
    const MEAS: usize = 60_000;
    const BIN: usize = 300;

    let lattice = build_chain(3, true);
    let model = ONModel::<2>::new(-1.0);
    let mut system = System::new(lattice, 2, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut rng);
        system.spin_at_mut(site, 2).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut kernel = MetropolisCore::with_strategy(OPSSStrategy::new().with_sigma(SIGMA));
    for _ in 0..THERM {
        kernel.sweep(&mut system, &model, &mut rng);
    }
    let mut energies = Vec::with_capacity(MEAS);
    let mut chirality_squared = Vec::with_capacity(MEAS);
    for _ in 0..MEAS {
        kernel.sweep(&mut system, &model, &mut rng);
        energies.push(system.energy);
        let kappa = triangle_chirality(&system.spins);
        chirality_squared.push(kappa * kappa);
    }
    (
        binned_stats(&energies, BIN),
        binned_stats(&chirality_squared, BIN),
    )
}

#[test]
fn afm_xy_triangle_ground_state_is_the_120_degree_manifold() {
    // Algebra sanity for the frustrated reference: the 120° configuration is
    // an exact ground state with E₀ = −3|J|/2 per triangle and |κ| = 1.
    let lattice = build_chain(3, true);
    let model = ONModel::<2>::new(-1.0);
    let mut spins = vec![0.0; 6];
    for (site, angle) in [
        0.0_f64,
        std::f64::consts::TAU / 3.0,
        2.0 * std::f64::consts::TAU / 3.0,
    ]
    .iter()
    .enumerate()
    {
        let (s, c) = angle.sin_cos();
        spins[2 * site] = c;
        spins[2 * site + 1] = s;
    }
    let ground_energy = model.compute_total_energy(&spins, &lattice, 1.0);
    assert!(
        (ground_energy + 1.5).abs() < 1e-12,
        "120° state energy {ground_energy} != −1.5"
    );
    assert!((triangle_chirality(&spins).abs() - 1.0).abs() < 1e-12);

    // Quadrature must reproduce the T→0 limit: ⟨E⟩ → −1.5 + kT·(2/2) harmonic
    // modes and ⟨κ²⟩ → 1 − O(1/β). At β = 50 the anharmonic remainder is
    // O(β⁻²) for E and O(β⁻¹) for 1 − ⟨κ²⟩.
    let (e, k2) = afm_xy_triangle_exact(50.0);
    let harmonic = -1.5 + 1.0 / 50.0;
    assert!(
        (e - harmonic).abs() < 5e-3,
        "quadrature low-T ⟨E⟩={e:.5}, harmonic prediction {harmonic:.5}"
    );
    assert!(k2 > 0.95, "quadrature low-T ⟨κ²⟩={k2:.5} should → 1");
}

#[test]
fn frustrated_afm_xy_triangle_metropolis_matches_exact_quadrature() {
    const TRI_SEEDS: usize = 12;
    let n_seeds = zscore_seed_count(TRI_SEEDS);
    for beta in [1.0_f64, 3.0] {
        let (e_exact, k2_exact) = afm_xy_triangle_exact(beta);
        let mut e_z_max = 0.0_f64;
        let mut k2_z_max = 0.0_f64;
        let mut e_z_sum = 0.0_f64;
        let mut k2_z_sum = 0.0_f64;
        for seed in 0..n_seeds as u64 {
            let ((e_mean, e_se), (k2_mean, k2_se)) = run_afm_triangle(beta, 0xF057 + seed);
            let z_e = (e_mean - e_exact) / e_se.max(1e-12);
            let z_k2 = (k2_mean - k2_exact) / k2_se.max(1e-12);
            e_z_max = e_z_max.max(z_e.abs());
            k2_z_max = k2_z_max.max(z_k2.abs());
            e_z_sum += z_e;
            k2_z_sum += z_k2;
        }
        eprintln!(
            "[frustrated β={beta:.2}] exact E={e_exact:.5} κ²={k2_exact:.5} | \
             max|z_E|={e_z_max:.2} max|z_κ²|={k2_z_max:.2}"
        );
        // Repo exact-value z convention: max |z| < 4 per seed, |z̄| < 2.
        assert!(
            e_z_max < 4.0,
            "β={beta}: ⟨E⟩ max |z| = {e_z_max:.2} vs exact {e_exact:.5}"
        );
        assert!(
            k2_z_max < 4.0,
            "β={beta}: ⟨κ²⟩ max |z| = {k2_z_max:.2} vs exact {k2_exact:.5}"
        );
        let n = n_seeds as f64;
        assert!(
            (e_z_sum / n).abs() < 2.0,
            "β={beta}: ⟨E⟩ mean z = {e_z_sum:.2}"
        );
        assert!(
            (k2_z_sum / n).abs() < 2.0,
            "β={beta}: ⟨κ²⟩ mean z = {k2_z_sum:.2}"
        );
    }
}
