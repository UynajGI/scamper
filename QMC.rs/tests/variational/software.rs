//! Software-behavior tests: determinism, Carlo.rs integration, input
//! validation (criterion G), and checkpoint round-trip / corruption
//! rejection.

use carlo_rs::{Context, MonteCarlo, Run, RunConfig, RunId, TaskId};
use qmc_rs::{
    ContinuumHamiltonian, GaussianTrap, HarmonicJastrow, HarmonicTrap, McMillanJastrow,
    PairPotential, Positions, SlaterDeterminant, VariationalError, VmcKernel, WaveFunction,
    WaveFunctionParams, DIM, VMC_CHECKPOINT_FORMAT,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use serde_json::json;

type Rng = Xoshiro256PlusPlus;
type TrapKernel = VmcKernel<GaussianTrap>;

fn trap_kernel(seed: u64, n_walkers: usize, n_particles: usize) -> TrapKernel {
    let wave_function = GaussianTrap::new(0.5, [0.1, 0.0, -0.1]).unwrap();
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(1.0, [0.1, 0.0, -0.1]).unwrap()).unwrap();
    let mut rng = Rng::seed_from_u64(seed);
    VmcKernel::new(
        wave_function,
        hamiltonian,
        n_walkers,
        n_particles,
        1.5,
        0.5,
        &mut rng,
    )
    .expect("valid kernel inputs")
}

#[test]
fn same_seed_bit_identical_through_the_carlo_adapter() {
    // Two identically seeded runs driven through the MonteCarlo adapter
    // (context RNG -> per-walker derived streams) must agree bit-for-bit;
    // a different seed must diverge.
    let drive = |seed: u64| -> Vec<Vec<f64>> {
        let mut kernel = trap_kernel(seed, 5, 4);
        let mut ctx = Context::new_with_binsize(Rng::seed_from_u64(seed), 10, 10);
        for _ in 0..100 {
            kernel.sweep(&mut ctx);
        }
        kernel
            .walkers()
            .iter()
            .map(|walker| walker.configuration().as_slice().to_vec())
            .collect()
    };
    let left = drive(0x5EED1);
    let right = drive(0x5EED1);
    assert_eq!(left, right);
    assert_ne!(left, drive(0x5EED2));
}

#[test]
fn run_from_parts_pipeline_reports_zero_variance_for_exact_state() {
    // The Carlo.rs Run integration (manual `Run::from_parts` route — see
    // adapters/carlo.rs for why FromParams is not implemented) must carry
    // the zero-variance property through to the final estimates: the
    // LocalEnergy mean equals the exact 3Nw/2 and its stderr is zero.
    let omega = 1.0_f64;
    let (n_walkers, n_particles) = (4, 5);
    let wave_function = GaussianTrap::new(omega / 2.0, [0.0; DIM]).unwrap();
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap();
    let mut rng = Rng::seed_from_u64(31);
    let kernel = VmcKernel::new(
        wave_function,
        hamiltonian,
        n_walkers,
        n_particles,
        1.5,
        0.6,
        &mut rng,
    )
    .unwrap();
    assert_eq!(kernel.name(), "VariationalMetropolis");

    let context = Context::new_with_binsize(Rng::seed_from_u64(77), 20, 10);
    let config = RunConfig {
        thermalization_sweeps: 20,
        measurement_sweeps: 200,
        binsize: 10,
        ..RunConfig::default()
    };
    let mut run = Run::from_parts(context, kernel, TaskId::new(0), RunId::new(0), config);
    run.run(220);

    let estimates = run.context_mut().finalize_measurements();
    let exact = 1.5 * n_particles as f64 * omega;
    let energy = &estimates["LocalEnergy"];
    assert!(
        (energy.mean - exact).abs() <= 1e-12,
        "LocalEnergy mean {} vs exact {exact}",
        energy.mean
    );
    assert!(
        energy.stderr <= 1e-12,
        "zero-variance state has stderr {}",
        energy.stderr
    );
    assert!(estimates.contains_key("LocalEnergySquared"));
    assert!(estimates.contains_key("LogGradSquared"));
    assert!(estimates.contains_key("EnergyPerParticle"));
}

#[test]
fn invalid_inputs_are_rejected_never_accepted_or_panicking() {
    // Criterion G across the whole public construction surface.
    let valid_hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(1.0, [0.0; DIM]).unwrap()).unwrap();

    // Ansatz parameters.
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(GaussianTrap::new(bad, [0.0; DIM]).is_err());
        assert!(McMillanJastrow::new(bad).is_err());
        assert!(HarmonicJastrow::new(bad).is_err());
    }
    // try_set_params length/content validation.
    let mut gaussian = GaussianTrap::new(0.5, [0.0; DIM]).unwrap();
    assert!(gaussian.try_set_params(&[f64::NAN]).is_err());
    assert!(gaussian.try_set_params(&[0.5, 0.5]).is_err());

    // Hamiltonians.
    assert!(HarmonicTrap::new(0.0, [0.0; DIM]).is_err());
    assert!(HarmonicTrap::new(1.0, [f64::NAN, 0.0, 0.0]).is_err());
    assert!(ContinuumHamiltonian::new(None, None).is_err());
    assert!(PairPotential::LennardJones {
        epsilon: 1.0,
        sigma: 0.0
    }
    .validate()
    .is_err());
    assert!(PairPotential::Harmonic {
        spring_constant: f64::NAN
    }
    .validate()
    .is_err());

    // Configurations.
    assert!(Positions::from_flat(vec![0.0; 7]).is_err());
    assert!(Positions::from_flat(vec![f64::NAN; 6]).is_err());

    // Kernel populations and widths.
    let build = |walkers: usize, particles: usize, spread: f64, width: f64| {
        let mut rng = Rng::seed_from_u64(1);
        VmcKernel::new(
            GaussianTrap::new(0.5, [0.0; DIM]).unwrap(),
            valid_hamiltonian,
            walkers,
            particles,
            spread,
            width,
            &mut rng,
        )
    };
    assert!(matches!(
        build(0, 4, 1.0, 0.5),
        Err(VariationalError::InvalidConfig { .. })
    ));
    assert!(build(8, 0, 1.0, 0.5).is_err());
    assert!(build(8, 4, 0.0, 0.5).is_err());
    assert!(build(8, 4, 1.0, f64::NAN).is_err());
    assert!(build(8, 4, 1.0, -0.5).is_err());

    // Non-finite wave-function parameters slip through `update_params`
    // (an optimizer-facing additive setter) but must still be caught by
    // the kernel constructor.
    let mut poisoned = GaussianTrap::new(0.5, [0.0; DIM]).unwrap();
    poisoned.update_params(&[f64::NAN]);
    assert!(build_kernel_with(poisoned).is_err());
}

fn build_kernel_with(wave_function: GaussianTrap) -> Result<TrapKernel, VariationalError> {
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(1.0, [0.0; DIM]).unwrap()).unwrap();
    let mut rng = Rng::seed_from_u64(2);
    VmcKernel::new(wave_function, hamiltonian, 4, 4, 1.5, 0.5, &mut rng)
}

#[test]
fn checkpoint_round_trip_replay_and_corruption_rejection() {
    let mut kernel = trap_kernel(0xBEEF, 4, 4);
    let mut rng = Rng::seed_from_u64(9);
    for _ in 0..50 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
    }
    let snapshot = kernel.save_snapshot();
    assert_eq!(snapshot["format"], VMC_CHECKPOINT_FORMAT);

    // Round trip: a freshly constructed kernel restores everything.
    let mut restored = trap_kernel(0x1234, 4, 4);
    restored.load_snapshot(&snapshot).unwrap();
    assert_eq!(restored.stats(), kernel.stats());
    for (a, b) in restored.walkers().iter().zip(kernel.walkers()) {
        assert_eq!(a.configuration().as_slice(), b.configuration().as_slice());
        assert_eq!(a.log_psi().to_bits(), b.log_psi().to_bits());
    }

    // Replay determinism: both continue identically under the same salts.
    let mut left_rng = Rng::seed_from_u64(0xC0);
    let mut right_rng = Rng::seed_from_u64(0xC0);
    for _ in 0..50 {
        kernel.sweep_with_phase(&mut left_rng, carlo_rs::RngPhase::Measurement);
        restored.sweep_with_phase(&mut right_rng, carlo_rs::RngPhase::Measurement);
    }
    for (a, b) in restored.walkers().iter().zip(kernel.walkers()) {
        assert_eq!(a.configuration().as_slice(), b.configuration().as_slice());
        assert_eq!(a.log_psi().to_bits(), b.log_psi().to_bits());
    }
    assert_eq!(restored.stats(), kernel.stats());

    // Loud rejection matrix (never a panic).
    let reject = |label: &str, snapshot: serde_json::Value| {
        let mut target = trap_kernel(0x1234, 4, 4);
        assert!(
            matches!(
                target.load_snapshot(&snapshot),
                Err(VariationalError::CheckpointCorrupted { .. })
            ),
            "{label}: expected loud CheckpointCorrupted"
        );
    };
    let mut bad;
    bad = snapshot.clone();
    bad["format"] = json!("qmc-rs-vmc-v2");
    reject("future format tag", bad);
    reject("empty object", json!({}));
    bad = snapshot.clone();
    bad["n_walkers"] = json!(5);
    reject("walker count mismatch", bad);
    bad = snapshot.clone();
    bad["n_particles"] = json!(3);
    reject("particle count mismatch", bad);
    bad = snapshot.clone();
    bad["walkers"] = json!(null);
    reject("walkers not an array", bad);
    bad = snapshot.clone();
    bad["walkers"].as_array_mut().unwrap().pop();
    reject("truncated walker list", bad);
    bad = snapshot.clone();
    bad["walkers"][1]["positions"][2] = json!("not a number");
    reject("non-numeric coordinate", bad);
    bad = snapshot.clone();
    bad["walkers"][2]["log_psi"] = json!(123.456);
    reject("tampered cached log_psi", bad);
    bad = snapshot.clone();
    bad["params"] = json!([0.7]);
    reject("parameter mismatch", bad);
    bad = snapshot.clone();
    bad["hamiltonian"]["trap"]["omega"] = json!(2.0);
    reject("hamiltonian mismatch", bad);
    bad = snapshot.clone();
    bad["sweeps"] = json!("many");
    reject("non-u64 counter", bad);
}

/// L1: checkpoint round-trip with a `SlaterDeterminant` kernel — the ansatz
/// state is fully described by its variational parameters, so a restored
/// kernel must replay bit-identically — plus the particle-count gate
/// (criterion G): a Slater ansatz fixes `n_particles = n_up + n_down`, and
/// the kernel constructor must reject any other count loudly (its
/// non-finite initial `log|psi|` check) instead of sampling a determinant
/// of the wrong shape.
#[test]
fn slater_kernel_checkpoint_round_trip_and_mismatch_rejection() {
    let omega = 1.3;
    let hamiltonian =
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap();
    let build_kernel = |n_particles: usize, seed: u64| {
        let mut rng = Rng::seed_from_u64(seed);
        VmcKernel::new(
            SlaterDeterminant::harmonic_trap(omega, 2).unwrap(),
            hamiltonian,
            4,
            n_particles,
            1.5,
            0.6,
            &mut rng,
        )
    };

    let mut kernel = build_kernel(8, 0x57A7).unwrap();
    let mut rng = Rng::seed_from_u64(0x57A7);
    for _ in 0..50 {
        kernel.sweep_with_phase(&mut rng, carlo_rs::RngPhase::Measurement);
    }
    let snapshot = kernel.save_snapshot();
    assert_eq!(snapshot["format"], VMC_CHECKPOINT_FORMAT);

    let mut restored = build_kernel(8, 0x9999).unwrap();
    restored.load_snapshot(&snapshot).unwrap();
    let mut left = Rng::seed_from_u64(0x0DD);
    let mut right = Rng::seed_from_u64(0x0DD);
    for _ in 0..50 {
        kernel.sweep_with_phase(&mut left, carlo_rs::RngPhase::Measurement);
        restored.sweep_with_phase(&mut right, carlo_rs::RngPhase::Measurement);
    }
    for (a, b) in restored.walkers().iter().zip(kernel.walkers()) {
        assert_eq!(a.configuration().as_slice(), b.configuration().as_slice());
        assert_eq!(a.log_psi().to_bits(), b.log_psi().to_bits());
    }
    assert_eq!(restored.stats(), kernel.stats());

    // Wrong particle count: loud rejection, never a wrong-shape walk.
    assert!(
        build_kernel(7, 1).is_err(),
        "particle-count mismatch must be rejected at construction"
    );
    // Corrupted snapshot on the Slater kernel too (same loud path).
    let mut bad = snapshot.clone();
    bad["format"] = json!("qmc-rs-vmc-v0");
    let mut target = build_kernel(8, 0x1234).unwrap();
    assert!(matches!(
        target.load_snapshot(&bad),
        Err(VariationalError::CheckpointCorrupted { .. })
    ));
}
