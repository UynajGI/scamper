use mcmc_rs::{
    DenseCovarianceAdaptation, EuclideanState, FnLogDensity, GaussianScale, RandomWalkMetropolis,
    SamplingPhase, TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn dense_adaptation_cholesky_reconstructs_regularized_covariance() {
    let mut adaptation = DenseCovarianceAdaptation::new(2, 1.0e-6).unwrap();
    for point in [[-2.0, -1.0], [-1.0, -2.0], [1.0, 0.0], [2.0, 3.0]] {
        adaptation.observe(&point).unwrap();
    }
    let covariance = adaptation.covariance().unwrap();
    let lower = adaptation.finalize_cholesky().unwrap();
    assert!(adaptation.is_frozen());

    for row in 0..2 {
        for column in 0..2 {
            let reconstructed = (0..2)
                .map(|inner| lower[row * 2 + inner] * lower[column * 2 + inner])
                .sum::<f64>();
            assert!((reconstructed - covariance[row * 2 + column]).abs() < 1.0e-10);
        }
    }
}

#[test]
fn dense_random_walk_freezes_to_dense_geometry_after_warmup() {
    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * (position[0] * position[0] + position[1] * position[1])
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0]).unwrap();
    let mut kernel = RandomWalkMetropolis::isotropic(2, 0.5)
        .unwrap()
        .with_dense_covariance_adaptation(1.0e-4)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(17);

    for _ in 0..200 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    kernel
        .on_phase_end(&mut target, SamplingPhase::Warmup, &state)
        .unwrap();
    assert!(kernel.adaptation_is_frozen());
    assert!(matches!(kernel.scale(), GaussianScale::Dense { .. }));
}

#[test]
fn dense_scale_generates_correlated_displacements_without_allocation_api() {
    let scale = GaussianScale::dense_cholesky(2, vec![1.0, 0.0, 0.8, 0.6]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(19);
    let mut normal = vec![0.0; 2];
    let mut displacement = vec![0.0; 2];
    let mut cross = 0.0;
    let mut first_square = 0.0;
    let mut second_square = 0.0;
    for _ in 0..20_000 {
        scale
            .fill_displacement(&mut rng, &mut normal, &mut displacement)
            .unwrap();
        cross += displacement[0] * displacement[1];
        first_square += displacement[0] * displacement[0];
        second_square += displacement[1] * displacement[1];
    }
    let correlation = cross / (first_square * second_square).sqrt();
    assert!(
        (correlation - 0.8).abs() < 0.04,
        "correlation={correlation}"
    );
}

#[test]
fn gaussian_geometry_rejects_zero_dimension_and_bad_covariance_shape() {
    assert!(RandomWalkMetropolis::isotropic(0, 1.0).is_err());
    assert!(RandomWalkMetropolis::diagonal(Vec::new()).is_err());
    assert!(GaussianScale::dense_cholesky(0, Vec::new()).is_err());
    assert!(GaussianScale::dense_from_covariance(2, &[1.0, 0.0, 1.0], 1.0e-6).is_err());
    assert!(RandomWalkMetropolis::dense_covariance(2, &[1.0, 0.5, 0.5, 1.0], 1.0e-8,).is_ok());
}

#[test]
fn v01_random_walk_and_report_json_remain_readable() {
    let kernel = RandomWalkMetropolis::isotropic(2, 0.4).unwrap();
    let mut kernel_json = serde_json::to_value(kernel).unwrap();
    let object = kernel_json.as_object_mut().unwrap();
    object.remove("dense_covariance_adaptation");
    object.remove("normal_buffer");
    let mut restored: RandomWalkMetropolis = serde_json::from_value(kernel_json).unwrap();

    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * (position[0] * position[0] + position[1] * position[1])
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(23);
    restored
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();

    let mut report_json = serde_json::to_value(mcmc_rs::TransitionReport::default()).unwrap();
    report_json
        .as_object_mut()
        .unwrap()
        .remove("subtransitions");
    let report: mcmc_rs::TransitionReport = serde_json::from_value(report_json).unwrap();
    assert_eq!(report.subtransitions, 0);
}
