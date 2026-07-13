//! Small stochastic checks against analytically soluble spin-boson limits.

use qmc_rs::spin_boson::{
    integrated_sigma_z, Bath, CouplingNormalization, SingleModeBath, SpinBosonModel,
    WormholeConfiguration, WormholeEngine,
};
use qmc_rs::{QmcKernel, UpdateSchedule};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn single_mode() -> Bath {
    Bath::SingleMode(SingleModeBath::new(1.0).expect("single mode"))
}

fn sample_model(
    model: SpinBosonModel,
    beta: f64,
    seed: u64,
    warmup: usize,
    samples: usize,
) -> (f64, f64, f64) {
    let mut engine = WormholeEngine::new(model, UpdateSchedule::new(4, 4, 64));
    let mut configuration = WormholeConfiguration::new(beta, 1).expect("configuration");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    for _ in 0..warmup {
        engine
            .sweep(&mut configuration, &mut rng)
            .expect("warmup sweep");
    }

    let mut magnetization_sum = 0.0;
    let mut order_sum = 0.0;
    let mut order_squared_sum = 0.0;
    for sample in 0..samples {
        engine
            .sweep(&mut configuration, &mut rng)
            .expect("measurement sweep");
        if sample % 1_000 == 0 {
            configuration
                .validate(engine.model())
                .expect("valid worldline during micro simulation");
        }
        let magnetization =
            integrated_sigma_z(&configuration, engine.model()).expect("magnetization");
        let order = configuration.expansion_order() as f64;
        magnetization_sum += magnetization;
        order_sum += order;
        order_squared_sum += order * order;
    }
    configuration
        .validate(engine.model())
        .expect("valid final worldline");

    let count = samples as f64;
    let mean_magnetization = magnetization_sum / count;
    let mean_order = order_sum / count;
    let order_variance = order_squared_sum / count - mean_order * mean_order;
    (mean_magnetization, mean_order, order_variance)
}

#[test]
fn constant_diagonal_activity_has_poisson_expansion_order() {
    let beta = 3.0;
    let constant = 0.6;
    let model = SpinBosonModel::xxz(single_mode(), 0.0, 0.0, 0.0, Some(constant))
        .expect("constant diagonal model");
    let (magnetization, mean_order, order_variance) =
        sample_model(model, beta, 7_001, 3_000, 20_000);
    let expected_order = beta * constant;

    assert!(magnetization.abs() < 0.12, "m={magnetization}");
    assert!(
        (mean_order - expected_order).abs() < 0.15,
        "<n>={mean_order}, exact={expected_order}"
    );
    assert!(
        (order_variance - expected_order).abs() < 0.25,
        "var(n)={order_variance}, exact={expected_order}"
    );
}

#[test]
fn diagonal_field_matches_two_state_partition_function() {
    let beta = 2.0;
    let field = 0.5;
    let constant = 0.7;
    let model = SpinBosonModel::xxz(single_mode(), 0.0, 0.0, field, Some(constant))
        .expect("diagonal field model");
    let (magnetization, mean_order, _) = sample_model(model, beta, 8_002, 4_000, 25_000);

    let exact_magnetization = (0.5_f64 * beta * field).tanh();
    let exact_order = beta * (constant + 0.5 * field * exact_magnetization);
    assert!(
        (magnetization - exact_magnetization).abs() < 0.12,
        "m={magnetization}, exact={exact_magnetization}"
    );
    assert!(
        (mean_order - exact_order).abs() < 0.18,
        "<n>={mean_order}, exact={exact_order}"
    );
}

#[test]
fn equal_rw_crw_amplitudes_preserve_spin_inversion_symmetry() {
    let model = SpinBosonModel::rw_crw(
        single_mode(),
        0.25,
        1.0,
        0.0,
        CouplingNormalization::FixedTotal,
        Some(0.3),
    )
    .expect("symmetric RW-CRW model");
    let (magnetization, mean_order, _) = sample_model(model, 4.0, 9_003, 4_000, 20_000);

    assert!(magnetization.abs() < 0.18, "m={magnetization}");
    assert!(mean_order > 0.2, "the interacting sector was not sampled");
}
