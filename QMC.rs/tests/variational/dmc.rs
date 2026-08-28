//! L3 DMC validation.
//!
//! # The exact anchor (derived, not cited)
//!
//! Two particles in an isotropic trap with a harmonic repulsive pair,
//! `H = -1/2 (lap_1^2 + lap_2^2) + 1/2 omega^2 (r_1^2 + r_2^2)
//! + 1/2 k |r_1 - r_2|^2`, separates in center-of-mass `R = (r_1+r_2)/2`
//! and relative `r = r_1 - r_2` coordinates:
//!
//! ```text
//! kinetic:  -1/2(lap_1^2+lap_2^2) = -1/4 lap_R^2 - lap_r^2
//! trap:     1/2 omega^2 (r_1^2+r_2^2) = omega^2 R^2 + 1/4 omega^2 r^2
//! CM:       mass 2, frequency omega            ->  E_R = 3/2 omega
//! rel:      mass 1/2, Omega^2 = omega^2 + 2k   ->  E_r = 3/2 Omega
//! E_0 = 3/2 (omega + sqrt(omega^2 + 2k))
//! ```
//!
//! The exact ground state is `exp(-omega R^2 - Omega r^2 / 4)`, i.e. the
//! `Product<GaussianTrap(alpha = omega/2), HarmonicJastrow(a = Omega/4)>`
//! ansatz. A nodeless trial state makes fixed-node exact, so DMC MUST
//! converge to `E_0` from any approximate nodeless psi_T — the flagship
//! gate. With the exact psi_T the local energy is constant, every
//! branching weight is identical, and the whole drift-diffusion +
//! branching + population-control machinery reduces to an identity.

use qmc_rs::{
    ContinuumHamiltonian, DmcKernel, GaussianTrap, HarmonicJastrow, HarmonicTrap, PairPotential,
    Positions, Product, VariationalError, VmcKernel, WaveFunctionParams, DIM,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

const OMEGA: f64 = 1.2;
const SPRING: f64 = 0.6;
/// `sqrt(omega^2 + 2k)` (const fn cannot call sqrt; value of
/// `sqrt(2.64)` to full precision).
const OMEGA_REL: f64 = 1.624807680927192;
/// `E_0 = 3/2 (omega + Omega)` — see the derivation above.
const EXACT_GROUND_STATE: f64 = 1.5 * (OMEGA + OMEGA_REL);

fn pair_trap_hamiltonian() -> ContinuumHamiltonian {
    ContinuumHamiltonian::new(
        Some(HarmonicTrap::new(OMEGA, [0.0; DIM]).unwrap()),
        Some(PairPotential::Harmonic {
            spring_constant: SPRING,
        }),
    )
    .unwrap()
}

/// Approximate nodeless trial state (both parameters deliberately off:
/// exact would be alpha = omega/2 ~= 0.6, a = Omega/4 ~= 0.406).
fn approximate_trial_state() -> Product<GaussianTrap, HarmonicJastrow> {
    Product::new(
        GaussianTrap::new(0.5, [0.0; DIM]).unwrap(),
        HarmonicJastrow::new(0.3).unwrap(),
    )
}

/// Exact trial state (zero variance under `pair_trap_hamiltonian`):
/// `psi = exp(-omega R^2 - Omega r^2 / 4)`, and since
/// `GaussianTrap(omega/2) = exp(-omega R^2 - omega r^2 / 4)`, the
/// HarmonicJastrow must carry only the REMAINDER of the relative
/// gaussian: `a = (Omega - omega)/4` (not `Omega/4` — with `a` too
/// large the `r^2` coefficient of E_L turns negative and the population
/// diffuses to the divergent tail; caught by the identity test).
fn exact_trial_state() -> Product<GaussianTrap, HarmonicJastrow> {
    Product::new(
        GaussianTrap::new(OMEGA / 2.0, [0.0; DIM]).unwrap(),
        HarmonicJastrow::new((OMEGA_REL - OMEGA) / 4.0).unwrap(),
    )
}

/// Batch-means standard error of a correlated step series.
fn batch_means_stderr(samples: &[f64], block: usize) -> f64 {
    let n_blocks = samples.len() / block;
    if n_blocks < 2 {
        return f64::INFINITY;
    }
    let means: Vec<f64> = (0..n_blocks)
        .map(|b| samples[b * block..(b + 1) * block].iter().sum::<f64>() / block as f64)
        .collect();
    let mean = means.iter().sum::<f64>() / means.len() as f64;
    let variance =
        means.iter().map(|m| (m - mean) * (m - mean)).sum::<f64>() / (means.len() - 1) as f64;
    (variance / means.len() as f64).sqrt()
}

fn build_dmc<W: WaveFunctionParams<Config = Positions> + Clone>(
    wave_function: W,
    n_walkers: usize,
    tau: f64,
    n_delay: usize,
    seed: u64,
) -> DmcKernel<W> {
    let mut rng = Rng::seed_from_u64(seed);
    DmcKernel::new(
        wave_function,
        pair_trap_hamiltonian(),
        n_walkers,
        2,
        tau,
        n_delay,
        EXACT_GROUND_STATE, // a sane initial energy keeps step 1 tame
        1.2,
        &mut rng,
    )
    .unwrap()
}

#[test]
fn dmc_converges_to_the_exact_separable_ground_state() {
    // Nodeless (fixed-node-exact) trial state, deliberately approximate:
    // DMC must land on E_0 within statistical precision, and strictly
    // closer than VMC at the same psi_T.
    let (vmc_mean, vmc_stderr) = {
        let mut rng = Rng::seed_from_u64(0x0DD1);
        let mut kernel = VmcKernel::new(
            approximate_trial_state(),
            pair_trap_hamiltonian(),
            16,
            2,
            1.5,
            0.7,
            &mut rng,
        )
        .unwrap();
        for _ in 0..200 {
            kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Thermalization);
        }
        let mut samples = Vec::new();
        for _ in 0..2000 {
            kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
            samples.push(kernel.population_mean_local_energy().value);
        }
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        (mean, batch_means_stderr(&samples, 100))
    };

    let mut kernel = build_dmc(approximate_trial_state(), 48, 0.02, 10, 0x0DD2);
    let mut rng = Rng::seed_from_u64(0x0DD2 ^ 0xA11CE);
    let equilibrate = 400;
    let steps = 3000;
    let mut prev = 0.0_f64;
    let mut prev_n = 0_u64;
    let mut per_step = Vec::with_capacity(steps);
    for step in 0..(equilibrate + steps) {
        kernel.step(&mut rng).unwrap();
        if step >= equilibrate {
            // Recover this step's population mean from the accumulators.
            let stats = kernel.stats();
            per_step.push((stats.energy_sum - prev) / (stats.n_energy_samples - prev_n) as f64);
            prev = stats.energy_sum;
            prev_n = stats.n_energy_samples;
        }
    }
    let dmc_mean = per_step.iter().sum::<f64>() / per_step.len() as f64;
    let dmc_stderr = batch_means_stderr(&per_step, 100);

    let z = (dmc_mean - EXACT_GROUND_STATE).abs() / dmc_stderr;
    assert!(
        z < 4.0,
        "E_DMC = {dmc_mean} +/- {dmc_stderr} vs exact {EXACT_GROUND_STATE} (z = {z})"
    );
    // VMC at the same psi_T sits above E_0 (variational) and farther from
    // it than DMC (projection improves a nodeless state).
    assert!(
        vmc_mean - EXACT_GROUND_STATE > -4.0 * vmc_stderr,
        "VMC violated the variational bound: {vmc_mean} vs {EXACT_GROUND_STATE}"
    );
    assert!(
        (dmc_mean - EXACT_GROUND_STATE).abs() < (vmc_mean - EXACT_GROUND_STATE).abs(),
        "DMC ({dmc_mean}) not closer to E_0 than VMC ({vmc_mean})"
    );
}

#[test]
fn dmc_reduces_to_identity_for_the_exact_trial_state() {
    // Zero variance => constant E_L => uniform branching weights and a
    // constant mixed estimator: the whole machinery collapses to the
    // exact energy regardless of drift, diffusion or population state.
    let mut kernel = build_dmc(exact_trial_state(), 24, 0.02, 8, 0x2ACE);
    let mut rng = Rng::seed_from_u64(0x2ACE ^ 0xA11CE);
    for _ in 0..600 {
        kernel.step(&mut rng).unwrap();
    }
    let mixed = kernel.stats().mixed_energy();
    assert!(
        (mixed - EXACT_GROUND_STATE).abs() <= 1e-10,
        "exact-state mixed estimator {mixed} vs {EXACT_GROUND_STATE}"
    );
    // The forward-walking pure estimator agrees (both are the constant).
    let pure = kernel.stats().pure_energy();
    assert!(
        (pure - EXACT_GROUND_STATE).abs() <= 1e-10,
        "exact-state pure estimator {pure} vs {EXACT_GROUND_STATE}"
    );
    // The population stayed healthy (no die-out, no explosion).
    let population = kernel.population();
    assert!((6..=200).contains(&population), "population {population}");
}

#[test]
fn population_control_bias_shrinks_with_the_population() {
    // Equal walker-step budgets: N = 8 with 4x the steps vs N = 64. The
    // population-control bias scales as O(1/N_w), so the larger
    // population's estimate is closer to (or as close as) the exact
    // energy. Point-estimate comparison; both runs are individually
    // z-consistent with E_0 at generous bounds.
    let (bias_small, z_small) = {
        let mut kernel = build_dmc(approximate_trial_state(), 8, 0.02, 8, 0x8177);
        let mut rng = Rng::seed_from_u64(0x8177 ^ 0xA11CE);
        let (equilibrate, steps) = (400, 8000);
        let mut prev = 0.0;
        let mut prev_n = 0;
        let mut per_step = Vec::with_capacity(steps);
        for step in 0..(equilibrate + steps) {
            kernel.step(&mut rng).unwrap();
            if step >= equilibrate {
                let stats = kernel.stats();
                per_step.push((stats.energy_sum - prev) / (stats.n_energy_samples - prev_n) as f64);
                prev = stats.energy_sum;
                prev_n = stats.n_energy_samples;
            }
        }
        let mean = per_step.iter().sum::<f64>() / per_step.len() as f64;
        let stderr = batch_means_stderr(&per_step, 200);
        (
            (mean - EXACT_GROUND_STATE).abs(),
            (mean - EXACT_GROUND_STATE).abs() / stderr,
        )
    };
    let (bias_large, z_large) = {
        let mut kernel = build_dmc(approximate_trial_state(), 64, 0.02, 8, 0x8188);
        let mut rng = Rng::seed_from_u64(0x8188 ^ 0xA11CE);
        let (equilibrate, steps) = (400, 1000);
        let mut prev = 0.0;
        let mut prev_n = 0;
        let mut per_step = Vec::with_capacity(steps);
        for step in 0..(equilibrate + steps) {
            kernel.step(&mut rng).unwrap();
            if step >= equilibrate {
                let stats = kernel.stats();
                per_step.push((stats.energy_sum - prev) / (stats.n_energy_samples - prev_n) as f64);
                prev = stats.energy_sum;
                prev_n = stats.n_energy_samples;
            }
        }
        let mean = per_step.iter().sum::<f64>() / per_step.len() as f64;
        let stderr = batch_means_stderr(&per_step, 50);
        (
            (mean - EXACT_GROUND_STATE).abs(),
            (mean - EXACT_GROUND_STATE).abs() / stderr,
        )
    };
    assert!(
        z_small.is_finite() && z_small < 6.0,
        "N=8 run: z = {z_small}"
    );
    assert!(
        z_large.is_finite() && z_large < 6.0,
        "N=64 run: z = {z_large}"
    );
    assert!(
        bias_large <= bias_small,
        "population-control bias did not shrink: |bias|_64 = {bias_large} vs |bias|_8 = {bias_small}"
    );
}

#[test]
fn dmc_software_gates_determinism_checkpoints_and_validation() {
    // Same seed -> bit-identical statistics.
    let run = |seed: u64| {
        let mut kernel = build_dmc(approximate_trial_state(), 12, 0.03, 6, seed);
        let mut rng = Rng::seed_from_u64(seed ^ 0x5EED);
        for _ in 0..200 {
            kernel.step(&mut rng).unwrap();
        }
        kernel
    };
    let left = run(0x1234);
    let right = run(0x1234);
    assert_eq!(
        format!("{:?}", left.stats()),
        format!("{:?}", right.stats())
    );
    assert_eq!(left.population(), right.population());

    // Checkpoint round-trip and bit-identical replay.
    let snapshot = left.save_snapshot();
    assert_eq!(
        snapshot["format"].as_str(),
        Some("qmc-rs-dmc-v1") // DMC_CHECKPOINT_FORMAT
    );
    let mut restored = run(0xBEEF);
    restored.load_snapshot(&snapshot).unwrap();
    let mut rng_a = Rng::seed_from_u64(0x0DD);
    let mut rng_b = Rng::seed_from_u64(0x0DD);
    let mut continued = left;
    for _ in 0..100 {
        continued.step(&mut rng_a).unwrap();
        restored.step(&mut rng_b).unwrap();
    }
    assert_eq!(
        format!("{:?}", continued.stats()),
        format!("{:?}", restored.stats())
    );

    // Corruption and mismatch rejections (loud, never a panic).
    let mut bad = snapshot.clone();
    bad["format"] = serde_json::json!("qmc-rs-dmc-v9");
    let mut target = run(0x7777);
    assert!(matches!(
        target.load_snapshot(&bad),
        Err(VariationalError::CheckpointCorrupted { .. })
    ));
    assert!(target.load_snapshot(&serde_json::json!({})).is_err());
    // Particle-count mismatch: a 3-particle kernel refuses a 2-particle
    // snapshot.
    let mut rng = Rng::seed_from_u64(3);
    let mut three_particle = DmcKernel::new(
        approximate_trial_state(),
        pair_trap_hamiltonian(),
        4,
        3,
        0.03,
        4,
        5.0,
        1.0,
        &mut rng,
    )
    .unwrap();
    assert!(matches!(
        three_particle.load_snapshot(&snapshot),
        Err(VariationalError::CheckpointCorrupted { .. })
    ));

    // Constructor validation (criterion G).
    let mut rng = Rng::seed_from_u64(5);
    let wf = approximate_trial_state();
    assert!(DmcKernel::new(
        wf,
        pair_trap_hamiltonian(),
        0,
        2,
        0.02,
        4,
        1.0,
        1.0,
        &mut rng
    )
    .is_err());
    assert!(DmcKernel::new(
        wf,
        pair_trap_hamiltonian(),
        4,
        0,
        0.02,
        4,
        1.0,
        1.0,
        &mut rng
    )
    .is_err());
    assert!(DmcKernel::new(
        wf,
        pair_trap_hamiltonian(),
        4,
        2,
        0.0,
        4,
        1.0,
        1.0,
        &mut rng
    )
    .is_err());
    assert!(DmcKernel::new(
        wf,
        pair_trap_hamiltonian(),
        4,
        2,
        -0.1,
        4,
        1.0,
        1.0,
        &mut rng
    )
    .is_err());
    assert!(DmcKernel::new(
        wf,
        pair_trap_hamiltonian(),
        4,
        2,
        0.02,
        0,
        1.0,
        1.0,
        &mut rng
    )
    .is_err());
    assert!(DmcKernel::new(
        approximate_trial_state(),
        pair_trap_hamiltonian(),
        4,
        2,
        0.02,
        4,
        f64::NAN,
        1.0,
        &mut rng
    )
    .is_err());
}
