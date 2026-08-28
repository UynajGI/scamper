//! L2 optimizer validation.
//!
//! Two complementary evidence layers:
//!
//! 1. **Deterministic quadrature statistics.** For the one-particle
//!    `GaussianTrap` family every block moment is analytically known:
//!    with `psi_alpha(r) = exp(-alpha r^2)` sampled at `|psi|^2`
//!    (weight `exp(-2 alpha r^2)`, each Cartesian component Gaussian
//!    `N(0, 1/(4 alpha))` in 3D):
//!    ```text
//!    ln psi  = -alpha r^2,            O      = d/dalpha ln psi = -r^2
//!    E_L(r)  = 3 alpha + c r^2,       c      = omega^2/2 - 2 alpha^2
//!    E(alpha)= 3 alpha/2 + 3 omega^2/(8 alpha),  minimum at alpha* = omega/2
//!    S(alpha)= Cov(O, E_L) = -c Var(r^2) = -c * 3/(8 alpha^2)
//!    G(alpha)= Var(r^2)              = 3/(8 alpha^2)
//!    ```
//!    A weighted quadrature grid (`push_weighted`) reproduces these to
//!    quadrature accuracy, so the optimizers can be iterated
//!    deterministically against exact statistics — convergence,
//!    step-quality comparisons and the linear-method eigenvalue bound
//!    become theorem-level assertions, not statistical ones.
//! 2. **Stochastic kernel statistics.** `collect_block_stats` through the
//!    real Metropolis kernel: the force vanishes for an exact state, and
//!    SR improves a genuinely two-parameter droplet.

use qmc_rs::{
    BlockStats, ContinuumHamiltonian, GaussianTrap, HarmonicTrap, LinearMethod, McMillanJastrow,
    Optimizer, PairPotential, Positions, Product, StochasticReconfiguration, VmcKernel,
    WaveFunctionParams, DIM,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Deterministic weighted quadrature of the one-particle Gaussian-trap
/// block statistics at parameter `alpha` (traps `omega`), over a
/// `points^3` grid on `[-half, half]^3` with weights `exp(-2 alpha r^2)`.
fn quadrature_stats(alpha: f64, omega: f64, points: usize, half: f64) -> BlockStats {
    let mut stats = BlockStats::new(1);
    let step = 2.0 * half / points as f64;
    for ix in 0..points {
        let x = -half + (ix as f64 + 0.5) * step;
        for iy in 0..points {
            let y = -half + (iy as f64 + 0.5) * step;
            for iz in 0..points {
                let z = -half + (iz as f64 + 0.5) * step;
                let r2 = x * x + y * y + z * z;
                let weight = (-2.0 * alpha * r2).exp();
                let local_energy = 3.0 * alpha + (0.5 * omega * omega - 2.0 * alpha * alpha) * r2;
                stats.push_weighted(local_energy, &[-r2], weight);
            }
        }
    }
    stats
}

/// Midpoint-grid resolution and box used below: `128^3` points on a
/// `[-5.5, 5.5]^3` box. The box dominates the error budget: at the
/// shallowest alpha = 0.25 the weight is `exp(-2 alpha r^2)`, so ±5.5 is
/// ±5.5 sigma per dimension and the truncated tail (~e^-15) is
/// negligible; a ±3.5 box leaves a ~3e-4 systematic in `<r^2>`-weighted
/// moments that no grid refinement removes (measured).
fn toy_stats(alpha: f64, omega: f64) -> BlockStats {
    quadrature_stats(alpha, omega, 128, 5.5)
}

#[test]
fn block_stats_reproduce_analytic_moments() {
    // The centering algebra of BlockStats (force, metric, energy) against
    // the closed forms in the module derivation, at two parameter values.
    let omega = 1.3;
    for alpha in [0.35_f64, 0.9] {
        let stats = toy_stats(alpha, omega);
        let expected_energy = 1.5 * alpha + 3.0 * omega * omega / (8.0 * alpha);
        let curvature = 0.5 * omega * omega - 2.0 * alpha * alpha;
        let var_r2 = 3.0 / (8.0 * alpha * alpha);
        let expected_force = -curvature * var_r2;
        let expected_metric = var_r2;
        let e = stats.energy();
        let s = stats.force()[0];
        let g = stats.metric()[(0, 0)];
        assert!(
            (e - expected_energy).abs() <= 1e-4 * expected_energy.abs(),
            "E: {e} vs {expected_energy}"
        );
        assert!(
            (s - expected_force).abs() <= 1e-4 * expected_force.abs(),
            "S: {s} vs {expected_force}"
        );
        assert!(
            (g - expected_metric).abs() <= 1e-4 * expected_metric.abs(),
            "G: {g} vs {expected_metric}"
        );
    }
}

#[test]
fn force_vanishes_for_the_exact_state() {
    // Zero-variance principle on the statistics layer: at alpha = omega/2
    // the local energy is configuration-independent, so its covariance
    // with ANY function of the configuration - the force included - is
    // zero. Quadrature-exact to the machine floor.
    let omega = 1.3;
    let stats = toy_stats(omega / 2.0, omega);
    assert!(
        stats.force()[0].abs() <= 1e-10 * stats.metric()[(0, 0)],
        "exact-state force {} vs metric {}",
        stats.force()[0],
        stats.metric()[(0, 0)]
    );
    // Zero variance up to the raw-moment cancellation floor: the variance
    // is computed as <E^2> - <E>^2 from accumulated sums (no second pass
    // over samples), so a configuration-independent E_L lands at
    // |Var| ~ eps * E^2 * O(1) ~ 1e-15..1e-11 - and may be negative.
    assert!(
        stats.energy_variance().abs() <= 1e-9,
        "exact-state variance {} exceeds the cancellation floor",
        stats.energy_variance()
    );
}

#[test]
fn linear_method_eigenvalue_never_exceeds_current_energy() {
    // The linear-method generalized eigenproblem contains the undisplaced
    // state as its first basis vector, so the lowest eigenvalue (the
    // predicted post-update energy) cannot exceed the current block
    // energy - a convexity property of the variational subspace.
    let omega = 1.3;
    for alpha in [0.3_f64, 0.45, 0.95] {
        let stats = toy_stats(alpha, omega);
        assert!(
            LinearMethod::predicted_energy(&stats) <= stats.energy() + 1e-9,
            "LM prediction exceeds E at alpha={alpha}"
        );
    }
}

/// Run one optimizer deterministically on quadrature statistics from
/// `alpha0` until the parameter is within `tol` of `omega/2` (or the
/// iteration cap hits). Returns the iteration count of convergence.
fn run_deterministic<O: Optimizer>(
    mut optimizer: O,
    alpha0: f64,
    omega: f64,
    tol: f64,
    max_iterations: usize,
) -> (usize, f64) {
    let mut alpha = alpha0;
    for iteration in 1..=max_iterations {
        let stats = toy_stats(alpha, omega);
        let delta = optimizer.propose(&stats).expect("valid statistics");
        optimizer.feedback(true);
        alpha += delta[0];
        if optimizer.converged() || (alpha - omega / 2.0).abs() < tol {
            return (iteration, alpha);
        }
    }
    (usize::MAX, alpha)
}

#[test]
fn stochastic_reconfiguration_converges_to_the_closed_form_optimum() {
    // Learning rate 0.3 keeps the fixed-point map contracting (the SR
    // step is eps' * c(alpha) with |1 - 4 eps' alpha*| < 1 iff
    // eps' < 1/(2 alpha*) ~ 0.77 here).
    let omega = 1.3;
    let sr = StochasticReconfiguration::new(0.3, 0.1, 1e-7, 3).unwrap();
    let (iterations, alpha) = run_deterministic(sr, 0.25, omega, 1e-4, 200);
    assert!(
        (alpha - omega / 2.0).abs() < 1e-3,
        "SR stalled at alpha={alpha}"
    );
    let energy = 1.5 * alpha + 3.0 * omega * omega / (8.0 * alpha);
    let exact = 1.5 * omega / 2.0 + 3.0 * omega * omega / (8.0 * (omega / 2.0));
    assert!(
        (energy - exact).abs() <= 1e-4,
        "E(alpha_final) = {energy} vs exact {exact} (iterations {iterations})"
    );
}

#[test]
fn linear_method_converges_on_the_deterministic_toy() {
    // L2 gate for the linear method on this landscape: convergence to the
    // closed-form optimum within a bounded iteration count. The step-
    // quality comparison against SR is deliberately NOT asserted here: on
    // a one-parameter smooth toy with a well-tuned learning rate SR wins
    // (8 vs ~17 iterations); the linear method's advantage - stable large
    // steps on stiff multi-parameter landscapes - has no room to appear,
    // and asserting iteration dominance would be a fake gate. The honest
    // model-level bound (predicted energy <= current energy, a convexity
    // property of the enlarged subspace) is covered separately below.
    let omega = 1.3;
    let lm = LinearMethod::new(0.5, 0.1, 1e-7, 3).unwrap();
    let (iterations, alpha) = run_deterministic(lm, 0.25, omega, 1e-4, 60);
    assert!(
        (alpha - omega / 2.0).abs() < 1e-3,
        "LM stalled at alpha={alpha} after {iterations} iterations"
    );
    let energy = 1.5 * alpha + 3.0 * omega * omega / (8.0 * alpha);
    let exact = 1.5 * omega;
    assert!(
        (energy - exact).abs() <= 1e-4,
        "E(alpha_final) = {energy} vs exact {exact} (iterations {iterations})"
    );
}

#[test]
fn trust_region_escalates_and_relaxes() {
    let mut sr = StochasticReconfiguration::new(0.5, 0.1, 1e-7, 3).unwrap();
    sr.feedback(false);
    assert!((sr.shift() - 0.2).abs() < 1e-15);
    sr.feedback(false);
    assert!((sr.shift() - 0.4).abs() < 1e-15);
    sr.feedback(true);
    assert!((sr.shift() - 0.2).abs() < 1e-15);

    let mut lm = LinearMethod::new(0.05, 0.1, 1e-7, 3).unwrap();
    lm.feedback(false);
    assert!((lm.metric_shift() - 0.2).abs() < 1e-15);
    lm.feedback(true);
    assert!((lm.metric_shift() - 0.1).abs() < 1e-15);
}

#[test]
fn optimizer_inputs_and_empty_blocks_are_rejected() {
    for bad in [0.0_f64, -1.0, f64::NAN] {
        assert!(StochasticReconfiguration::new(bad, 0.1, 1e-7, 3).is_err());
        assert!(LinearMethod::new(bad, 0.1, 1e-7, 3).is_err());
    }
    assert!(StochasticReconfiguration::new(0.5, 0.1, 1e-7, 0).is_err());
    assert!(LinearMethod::new(0.05, 0.1, 1e-7, 0).is_err());
    let mut sr = StochasticReconfiguration::new(0.5, 0.1, 1e-7, 3).unwrap();
    assert!(sr.propose(&BlockStats::new(1)).is_err());
    let mut lm = LinearMethod::new(0.05, 0.1, 1e-7, 3).unwrap();
    assert!(lm.propose(&BlockStats::new(1)).is_err());
}

#[test]
fn convergence_requires_patience_consecutive_quiet_blocks() {
    // Feeding exact-state statistics (zero force) must report convergence
    // only after `patience` consecutive blocks; one loud block resets.
    let omega = 1.3;
    let mut sr = StochasticReconfiguration::new(1.0, 0.1, 1e-7, 3).unwrap();
    let quiet = toy_stats(omega / 2.0, omega);
    let loud = toy_stats(0.25, omega);
    sr.propose(&quiet).unwrap();
    assert!(!sr.converged());
    sr.propose(&quiet).unwrap();
    assert!(!sr.converged());
    sr.propose(&loud).unwrap(); // resets the patience counter
    sr.propose(&quiet).unwrap();
    sr.propose(&quiet).unwrap();
    assert!(!sr.converged());
    sr.propose(&quiet).unwrap();
    assert!(sr.converged());
}

#[test]
fn kernel_block_stats_zero_force_for_exact_state_through_metropolis() {
    // The stochastic counterpart of the zero-force identity: sampling the
    // exact Gaussian at alpha = omega/2 through the kernel, the measured
    // force is zero to statistical precision (E_L is a constant, so its
    // sampled covariance with O is pure noise ~ sigma/sqrt(n) - asserted
    // at a generous 10-sigma-of-the-metric scale).
    let omega = 1.1;
    let n_particles = 5;
    let wave_function = GaussianTrap::new(omega / 2.0, [0.0; DIM]).unwrap();
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap();
    let mut rng = Rng::seed_from_u64(0x17A2);
    let mut kernel = VmcKernel::new(
        wave_function,
        hamiltonian,
        8,
        n_particles,
        1.5,
        0.6,
        &mut rng,
    )
    .unwrap();

    let mut stats = BlockStats::new(1);
    for _ in 0..400 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
        kernel.collect_block_stats(&mut stats);
    }
    let force = stats.force()[0];
    let metric = stats.metric()[(0, 0)];
    let scale = (metric / stats.n_samples() as f64).sqrt();
    assert!(
        force.abs() <= 10.0 * scale,
        "exact-state force {force} vs noise scale {scale}"
    );
    // ...and the local-energy variance is the zero-variance signature.
    assert!(stats.energy_variance() <= 1e-10);
}

#[test]
fn sr_improves_the_two_parameter_droplet() {
    // Stochastic SR end-to-end on a genuinely two-parameter problem:
    // Product<GaussianTrap(alpha), McMillanJastrow(b)> under
    // trap(omega) + LennardJones. From a deliberately poor start the
    // block energy must improve beyond noise. No exact optimum is
    // asserted - the honest claims are improvement and physical
    // parameters throughout.
    let omega = 0.05;
    let n_particles = 8;
    let build = |alpha: f64, b: f64| {
        Product::new(
            GaussianTrap::new(alpha, [0.0; DIM]).unwrap(),
            McMillanJastrow::new(b).unwrap(),
        )
    };
    let hamiltonian = ContinuumHamiltonian::new(
        Some(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()),
        Some(PairPotential::LennardJones {
            epsilon: 1.0,
            sigma: 1.0,
        }),
    )
    .unwrap();

    let mut rng = Rng::seed_from_u64(0x5211);
    let mut kernel = VmcKernel::new(
        build(0.10, 2.5),
        hamiltonian,
        8,
        n_particles,
        1.5,
        0.8,
        &mut rng,
    )
    .unwrap();

    // Thermalize, then measure the baseline block at the poor start.
    for _ in 0..100 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Thermalization);
    }
    let mut baseline = BlockStats::new(2);
    for _ in 0..60 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
        kernel.collect_block_stats(&mut baseline);
    }
    let e_before = baseline.energy();
    let stderr_before = baseline.energy_stderr_naive();

    // Optimization loop: propose from the current block, apply through
    // the kernel's parameter API (walker configurations kept - common
    // random numbers), re-block, keep the step only if the energy did not
    // worsen beyond noise; otherwise revert the parameters.
    let mut optimizer = StochasticReconfiguration::new(0.3, 0.3, 1e-4, 3).unwrap();
    let mut e_current = e_before;
    for _ in 0..12 {
        let mut probe = BlockStats::new(2);
        for _ in 0..60 {
            kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
            kernel.collect_block_stats(&mut probe);
        }
        let Ok(delta) = optimizer.propose(&probe) else {
            break;
        };
        if delta.iter().all(|&d| d == 0.0) {
            break; // converged by patience
        }
        if kernel.update_wave_function_params(&delta).is_err() {
            optimizer.feedback(false);
            continue;
        }
        let mut candidate = BlockStats::new(2);
        for _ in 0..60 {
            kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
            kernel.collect_block_stats(&mut candidate);
        }
        let accepted = candidate.energy() <= e_current + 2.0 * stderr_before;
        optimizer.feedback(accepted);
        if accepted {
            e_current = candidate.energy();
        } else {
            let revert: Vec<f64> = delta.iter().map(|&d| -d).collect();
            kernel
                .update_wave_function_params(&revert)
                .expect("reverting a previously valid parameter set");
        }
    }

    assert!(
        e_current <= e_before - 2.0 * stderr_before,
        "E before {e_before} +/- {stderr_before}, after {e_current}"
    );
    let params = kernel.wave_function().param_values();
    assert!(
        params.iter().all(|p| p.is_finite()),
        "parameters left the physical domain: {params:?}"
    );
}

// ---- L2-b: correlated-sampling variance minimization (argmin) ---------

use qmc_rs::{ReferenceSample, VarianceMinimization, VarianceObjective};

/// Uniform-grid reference samples for the one-particle Gaussian toy at
/// reference parameter `alpha0`: with the correlated reweighting
/// `w = exp(2(ln psi_alpha - ln psi_alpha0))` the uniform measure times
/// the weight is exactly `exp(-2 alpha r^2)` at any candidate alpha, so
/// the objective is deterministic quadrature of the variance.
fn grid_samples(alpha0: f64, points: usize, half: f64) -> Vec<ReferenceSample> {
    let reference = GaussianTrap::new(alpha0, [0.0; DIM]).unwrap();
    let step = 2.0 * half / points as f64;
    let mut samples = Vec::with_capacity(points * points * points);
    for ix in 0..points {
        let x = -half + (ix as f64 + 0.5) * step;
        for iy in 0..points {
            let y = -half + (iy as f64 + 0.5) * step;
            for iz in 0..points {
                let z = -half + (iz as f64 + 0.5) * step;
                let cfg = Positions::from_flat(vec![x, y, z]).unwrap();
                samples.push(ReferenceSample::new(&reference, cfg));
            }
        }
    }
    samples
}

#[test]
fn variance_objective_matches_closed_form_and_vanishes_at_exact_state() {
    // Var(alpha) = c(alpha)^2 Var_{psi_alpha}(r^2)
    //            = c(alpha)^2 * 3/(8 alpha^2),  c = omega^2/2 - 2 alpha^2,
    // which is zero exactly at alpha* = omega/2 (zero variance).
    let omega = 1.3;
    let alpha0 = 0.35;
    let objective = VarianceObjective::new(
        GaussianTrap::new(alpha0, [0.0; DIM]).unwrap(),
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
        grid_samples(alpha0, 40, 5.5),
    )
    .unwrap();
    for alpha in [0.35_f64, 0.5, 0.8] {
        let curvature = 0.5 * omega * omega - 2.0 * alpha * alpha;
        let expected = curvature * curvature * 3.0 / (8.0 * alpha * alpha);
        let sampled = objective.variance_at(&[alpha]);
        assert!(
            (sampled - expected).abs() <= 1e-3 * expected.max(1e-12),
            "Var({alpha}) = {sampled} vs closed form {expected}"
        );
    }
    // Two-pass estimator: the exact-state variance sits at the clean
    // quadrature floor, not at a cancellation floor.
    let exact = objective.variance_at(&[omega / 2.0]);
    assert!(exact <= 1e-9, "exact-state variance {exact}");
}

#[test]
fn variance_minimization_converges_deterministically() {
    let omega = 1.3;
    let alpha0 = 0.35;
    let minimizer = VarianceMinimization::new(200, 1e-10, 0.3).unwrap();
    let result = minimizer
        .minimize(
            GaussianTrap::new(alpha0, [0.0; DIM]).unwrap(),
            ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
            grid_samples(alpha0, 40, 5.5),
        )
        .unwrap();
    assert!(
        (result.params[0] - omega / 2.0).abs() <= 1e-3,
        "alpha = {} vs omega/2 = {}",
        result.params[0],
        omega / 2.0
    );
    assert!(result.variance <= 1e-8, "variance {}", result.variance);
    assert!(result.iterations < 200);
}

#[test]
fn out_of_domain_parameters_cost_a_penalty_not_a_panic() {
    // A negative Gaussian width leaves the ansatz's domain; the objective
    // must return the finite penalty so the simplex can walk around the
    // infeasible vertex (criterion G spirit on the optimizer surface).
    let omega = 1.3;
    let objective = VarianceObjective::new(
        GaussianTrap::new(0.35, [0.0; DIM]).unwrap(),
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
        grid_samples(0.35, 8, 2.0),
    )
    .unwrap();
    assert_eq!(objective.variance_at(&[-0.5]), 1e6);
    assert!(objective.variance_at(&[f64::NAN]).is_finite());
}

#[test]
fn variance_minimization_on_kernel_sampled_configurations() {
    // Stochastic end-to-end: configurations sampled by the Metropolis
    // kernel at a deliberately poor alpha0; the correlated-sampling
    // reweighting drives the variance down by more than an order of
    // magnitude and lands alpha near the exact omega/2.
    let omega = 1.1;
    let alpha0 = 0.3 * (omega / 2.0);
    let mut rng = Rng::seed_from_u64(0x7A17);
    let mut kernel = VmcKernel::new(
        GaussianTrap::new(alpha0, [0.0; DIM]).unwrap(),
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
        8,
        5,
        1.5,
        0.8,
        &mut rng,
    )
    .unwrap();
    for _ in 0..80 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Thermalization);
    }
    let reference = GaussianTrap::new(alpha0, [0.0; DIM]).unwrap();
    let mut samples = Vec::new();
    for _ in 0..150 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
        for walker in kernel.walkers() {
            samples.push(ReferenceSample::new(
                &reference,
                walker.configuration().clone(),
            ));
        }
    }

    let objective = VarianceObjective::new(
        reference,
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
        samples.clone(),
    )
    .unwrap();
    let variance_before = objective.variance_at(&[alpha0]);

    let minimizer = VarianceMinimization::new(150, 1e-8, 0.3).unwrap();
    let result = minimizer
        .minimize(
            GaussianTrap::new(alpha0, [0.0; DIM]).unwrap(),
            ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
            samples,
        )
        .unwrap();
    assert!(
        variance_before / result.variance.max(1e-12) >= 10.0,
        "variance before {variance_before}, after {}",
        result.variance
    );
    assert!(
        (result.params[0] - omega / 2.0).abs() <= 0.1 * (omega / 2.0),
        "alpha = {} vs omega/2 = {}",
        result.params[0],
        omega / 2.0
    );
}
