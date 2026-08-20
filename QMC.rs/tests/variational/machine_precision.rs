//! Machine-precision validation of the L0 variational pipeline.
//!
//! The signature property of a correct VMC implementation: for a trial
//! state that is the *exact* eigenstate of the sampled Hamiltonian, the
//! local energy is a configuration-independent constant ("zero variance"),
//! so every Metropolis sample reproduces the exact energy to machine
//! precision. Together with `delta_log` ≡ full-recompute equivalence and
//! finite-difference agreement of every hand-derived derivative, this pins
//! the whole pipeline (ansatz derivatives, estimator, kernel).

use qmc_rs::{
    local_energy, ContinuumHamiltonian, GaussianTrap, GradBuffer, HarmonicJastrow, HarmonicTrap,
    McMillanJastrow, PairPotential, ParamGradBuffer, Positions, Product, WaveFunction,
    WaveFunctionParams, DIM,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Deterministic random configuration: uniform in a cube of half-width
/// `scale` (kept ≳ 0.9 so pair distances stay away from the Jastrow core).
fn random_positions(rng: &mut Rng, n_particles: usize, scale: f64) -> Positions {
    let coords = (0..n_particles * DIM)
        .map(|_| rng.random_range(-scale..scale))
        .collect();
    Positions::from_flat(coords).expect("finite coordinates")
}

/// A copy of `cfg` with one particle replaced.
fn moved(cfg: &Positions, particle: usize, new_pos: [f64; DIM]) -> Positions {
    let mut moved = cfg.clone();
    moved.set_particle(particle, new_pos);
    moved
}

#[test]
fn gaussian_trap_zero_variance_at_exact_parameter() {
    // psi = prod exp(-alpha |r - r0|^2) is the exact ground state of
    // V = sum 1/2 omega^2 |r - r0|^2 at alpha = omega/2:
    //   E_L = 3N alpha + (1/2 omega^2 - 2 alpha^2) sum |r - r0|^2
    //       = 3N omega/2   (constant, since 1/2 w^2 = 2 (w/2)^2).
    let omega = 1.3;
    let alpha = omega / 2.0;
    let n_particles = 7;
    let center = [0.2, -0.4, 0.6];
    let wave_function = GaussianTrap::new(alpha, center).unwrap();
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, center).unwrap()).unwrap();
    let exact = 1.5 * n_particles as f64 * omega;

    let mut rng = Rng::seed_from_u64(0xD1CE);
    let mut grad = GradBuffer::new(n_particles);
    for _ in 0..50 {
        let cfg = random_positions(&mut rng, n_particles, 1.4);
        let sample = local_energy(&wave_function, &hamiltonian, &cfg, &mut grad);
        assert!(
            (sample.value - exact).abs() <= 1e-12,
            "E_L = {} vs exact {exact}",
            sample.value
        );
    }
}

#[test]
fn harmonic_jastrow_zero_variance_under_pair_harmonic_trap() {
    // psi = prod_{i<j} exp(-a r_ij^2) is the exact nodeless ground state of
    // H = -1/2 sum_i lap_i + sum_{i<j} 1/2 k r_ij^2 with k = 4 a^2 N, with
    // total energy E_0 = 3 a N (N-1). Derivation (HarmonicJastrow docs):
    //   Q = sum_{i<j} r_ij^2 = N sum_i |r_i - R_cm|^2
    //   sum_i |grad_i ln|psi||^2 = 4 a^2 N Q
    //   sum_i lap_i ln|psi|      = -6 a N (N-1)
    //   E_L = 3aN(N-1) + (k/2 - 2 a^2 N) Q -> constant 3aN(N-1) at k = 4a^2N.
    let a = 0.35;
    let n_particles = 6;
    let wave_function = HarmonicJastrow::new(a).unwrap();
    let k = wave_function.exact_pair_spring_constant(n_particles);
    assert!((k - 4.0 * a * a * n_particles as f64).abs() < 1e-15);
    let hamiltonian =
        ContinuumHamiltonian::pair_only(PairPotential::Harmonic { spring_constant: k }).unwrap();
    let exact = wave_function.exact_ground_state_energy(n_particles);

    let mut rng = Rng::seed_from_u64(0x5EED);
    let mut grad = GradBuffer::new(n_particles);
    for _ in 0..50 {
        let cfg = random_positions(&mut rng, n_particles, 0.9);
        let sample = local_energy(&wave_function, &hamiltonian, &cfg, &mut grad);
        assert!(
            (sample.value - exact).abs() <= 1e-11,
            "E_L = {} vs exact {exact}",
            sample.value
        );
    }
}

#[test]
fn zero_variance_survives_the_metropolis_pipeline() {
    // The same Gaussian exactness, sampled through the full kernel: after
    // every sweep of single-particle Metropolis the population mean local
    // energy must equal the exact value to machine precision.
    use qmc_rs::VmcKernel;

    let omega = 1.1;
    let (n_walkers, n_particles) = (6, 5);
    let center = [0.0; DIM];
    let wave_function = GaussianTrap::new(omega / 2.0, center).unwrap();
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, center).unwrap()).unwrap();
    let mut rng = Rng::seed_from_u64(0xABCD);
    let mut kernel = VmcKernel::new(
        wave_function,
        hamiltonian,
        n_walkers,
        n_particles,
        1.5,
        0.6,
        &mut rng,
    )
    .unwrap();

    let exact = 1.5 * n_particles as f64 * omega;
    let mut worst = 0.0_f64;
    for _ in 0..300 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
        worst = worst.max((kernel.population_mean_local_energy().value - exact).abs());
    }
    assert!(worst <= 1e-12, "worst sweep deviation {worst}");
}

/// Smallest pairwise distance in a configuration.
fn min_pair_distance(cfg: &Positions) -> f64 {
    let n = cfg.n_particles();
    (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| {
            let a = cfg.particle(i);
            let b = cfg.particle(j);
            let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

/// Fast-path honesty invariant for one ansatz: the incremental O(N)
/// single-particle log-ratio must equal the difference of full O(N^2)
/// recomputes. Configurations and trial positions are resampled to keep
/// every pair distance >= 0.8 so that |ln|psi|| = O(10) and the <= 1e-14
/// assertion is exactly the eps * |ln|psi|| machine-precision floor of the
/// two summation orders (the identity itself is exact in exact arithmetic).
fn check_delta_log<W: WaveFunction<Config = Positions>>(
    label: &str,
    wave_function: &W,
    rng: &mut Rng,
) {
    let n_particles = 4;
    let mut worst = 0.0_f64;
    for _ in 0..40 {
        let mut cfg;
        loop {
            cfg = random_positions(rng, n_particles, 1.3);
            if min_pair_distance(&cfg) >= 0.8 {
                break;
            }
        }
        for particle in 0..n_particles {
            let old = cfg.particle(particle);
            let mut new_pos;
            loop {
                new_pos = [
                    old[0] + rng.random_range(-0.3..0.3),
                    old[1] + rng.random_range(-0.3..0.3),
                    old[2] + rng.random_range(-0.3..0.3),
                ];
                let mut moved_cfg = cfg.clone();
                moved_cfg.set_particle(particle, new_pos);
                if min_pair_distance(&moved_cfg) >= 0.8 {
                    break;
                }
            }
            let incremental = wave_function.delta_log(&cfg, particle, &new_pos).log_ratio;
            let full = wave_function.log_psi(&moved(&cfg, particle, new_pos))
                - wave_function.log_psi(&cfg);
            worst = worst.max((incremental - full).abs());
        }
    }
    assert!(
        worst <= 1e-14,
        "{label}: delta_log vs full recompute worst error {worst}"
    );
}

#[test]
fn delta_log_matches_full_recompute_for_all_ansatze() {
    let mut rng = Rng::seed_from_u64(0xDE1A);
    check_delta_log(
        "GaussianTrap",
        &GaussianTrap::new(0.6, [0.1, -0.1, 0.2]).unwrap(),
        &mut rng,
    );
    check_delta_log(
        "McMillanJastrow",
        &McMillanJastrow::new(1.0).unwrap(),
        &mut rng,
    );
    check_delta_log(
        "HarmonicJastrow",
        &HarmonicJastrow::new(0.4).unwrap(),
        &mut rng,
    );
    check_delta_log(
        "Product",
        &Product::new(
            GaussianTrap::new(0.5, [0.1, -0.1, 0.2]).unwrap(),
            McMillanJastrow::new(1.0).unwrap(),
        ),
        &mut rng,
    );
}

/// Hand-derived gradient and Laplacian against central finite differences
/// for one ansatz. Gradient: first differences (h = 1e-6; truncation O(h^2),
/// roundoff ~eps|ln psi|/h -> observed ~1e-9). Total Laplacian: second
/// differences (h = 1e-4; roundoff ~4 eps |ln psi| / h^2 is the FD floor,
/// ~1e-6 relative at these scales).
fn check_derivatives<W: WaveFunction<Config = Positions>>(
    label: &str,
    wave_function: &W,
    cfg: &Positions,
) {
    let n_particles = cfg.n_particles();

    let mut grad = GradBuffer::new(n_particles);
    grad.clear();
    wave_function.log_grad(cfg, &mut grad);
    let h = 1e-6;
    for particle in 0..n_particles {
        for k in 0..DIM {
            let bump = |delta: f64| {
                let p = moved(cfg, particle, {
                    let mut p = cfg.particle(particle);
                    p[k] += delta;
                    p
                });
                wave_function.log_psi(&p)
            };
            let fd = (bump(h) - bump(-h)) / (2.0 * h);
            let analytic = grad.as_slice()[DIM * particle + k];
            assert!(
                (fd - analytic).abs() <= 1e-6,
                "{label}: grad[{particle}][{k}] FD {fd} vs analytic {analytic}"
            );
        }
    }

    let h = 1e-4;
    let mut fd_total = 0.0;
    for particle in 0..n_particles {
        for k in 0..DIM {
            let bump = |delta: f64| {
                let mut p = cfg.particle(particle);
                p[k] += delta;
                wave_function.log_psi(&moved(cfg, particle, p))
            };
            fd_total += (bump(h) - 2.0 * wave_function.log_psi(cfg) + bump(-h)) / (h * h);
        }
    }
    let analytic = wave_function.log_laplacian(cfg);
    assert!(
        (fd_total - analytic).abs() <= 1e-6 * (1.0 + analytic.abs()),
        "{label}: laplacian FD {fd_total} vs analytic {analytic}"
    );
}

#[test]
fn log_grad_and_laplacian_vs_central_finite_differences() {
    let mut rng = Rng::seed_from_u64(0xFADE);
    check_derivatives(
        "GaussianTrap",
        &GaussianTrap::new(0.6, [0.15, -0.25, 0.35]).unwrap(),
        &random_positions(&mut rng, 5, 1.1),
    );
    check_derivatives(
        "McMillanJastrow",
        &McMillanJastrow::new(1.0).unwrap(),
        &random_positions(&mut rng, 5, 1.1),
    );
    check_derivatives(
        "HarmonicJastrow",
        &HarmonicJastrow::new(0.45).unwrap(),
        &random_positions(&mut rng, 5, 1.1),
    );
    check_derivatives(
        "Product",
        &Product::new(
            GaussianTrap::new(0.5, [0.1, -0.1, 0.2]).unwrap(),
            McMillanJastrow::new(1.0).unwrap(),
        ),
        &random_positions(&mut rng, 5, 1.1),
    );
}

/// Parameter gradient against central parameter differences for one ansatz.
fn check_param_gradient<W>(label: &str, wave_function: &W, cfg: &Positions)
where
    W: WaveFunctionParams<Config = Positions> + Clone,
{
    let mut buffer = ParamGradBuffer::new(wave_function.n_params());
    buffer.clear();
    wave_function.log_grad_params(cfg, &mut buffer);
    let h = 1e-6;
    for (k, &analytic) in buffer.as_slice().iter().enumerate() {
        let mut delta = vec![0.0; wave_function.n_params()];
        let mut plus = wave_function.clone();
        delta[k] = h;
        plus.update_params(&delta);
        let mut minus = wave_function.clone();
        delta[k] = -h;
        minus.update_params(&delta);
        let fd = (plus.log_psi(cfg) - minus.log_psi(cfg)) / (2.0 * h);
        assert!(
            (fd - analytic).abs() <= 1e-6 * (1.0 + analytic.abs()),
            "{label}: d/dp{k} FD {fd} vs analytic {analytic}"
        );
    }
}

#[test]
fn log_grad_params_vs_finite_differences() {
    let mut rng = Rng::seed_from_u64(0xFEED);
    let gaussian = GaussianTrap::new(0.6, [0.1, -0.1, 0.2]).unwrap();
    check_param_gradient(
        "GaussianTrap",
        &gaussian,
        &random_positions(&mut rng, 5, 1.2),
    );
    check_param_gradient(
        "McMillanJastrow",
        &McMillanJastrow::new(1.0).unwrap(),
        &random_positions(&mut rng, 5, 1.2),
    );
    check_param_gradient(
        "HarmonicJastrow",
        &HarmonicJastrow::new(0.4).unwrap(),
        &random_positions(&mut rng, 5, 1.2),
    );
    check_param_gradient(
        "Product",
        &Product::new(gaussian, McMillanJastrow::new(1.0).unwrap()),
        &random_positions(&mut rng, 5, 1.2),
    );
}
