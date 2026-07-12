use mcmc_rs::{
    ComponentWiseMetropolis, EuclideanState, FnLogDensity, McmcError, RandomWalkMetropolis,
    SamplingPhase, SliceSampler, TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn fallible_component_transition_leaves_accepted_state_unchanged() {
    let mut calls = 0_u32;
    let mut target = FnLogDensity::new(move |_position: &[f64]| {
        calls += 1;
        if calls == 1 {
            0.0
        } else {
            f64::NAN
        }
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0]).unwrap();
    let original_position = state.position().clone();
    let original_log_density = state.log_density();
    let original_iteration = state.iteration();
    let mut kernel = ComponentWiseMetropolis::new(vec![1.0, 1.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(101);

    let error = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap_err();

    assert!(matches!(error, McmcError::InvalidLogDensity { .. }));
    assert_eq!(state.position(), &original_position);
    assert_eq!(state.log_density(), original_log_density);
    assert_eq!(state.iteration(), original_iteration);
}

#[test]
fn fallible_slice_transition_leaves_accepted_state_unchanged() {
    let mut calls = 0_u32;
    let mut target = FnLogDensity::new(move |_position: &[f64]| {
        calls += 1;
        if calls == 1 {
            0.0
        } else {
            f64::NAN
        }
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let original_position = state.position().clone();
    let original_log_density = state.log_density();
    let original_iteration = state.iteration();
    let mut kernel = SliceSampler::new(vec![1.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(102);

    let error = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap_err();

    assert!(matches!(error, McmcError::InvalidLogDensity { .. }));
    assert_eq!(state.position(), &original_position);
    assert_eq!(state.log_density(), original_log_density);
    assert_eq!(state.iteration(), original_iteration);
}

#[test]
fn every_kernel_advances_iteration_once_per_transition() {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(103);

    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0, 0.0]).unwrap();
    let mut random_walk = RandomWalkMetropolis::isotropic(3, 0.5).unwrap();
    random_walk
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(state.iteration(), 1);

    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0, 0.0]).unwrap();
    let mut component = ComponentWiseMetropolis::new(vec![0.5, 0.5, 0.5]).unwrap();
    component
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(state.iteration(), 1);

    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0, 0.0]).unwrap();
    let mut slice = SliceSampler::new(vec![1.0, 1.0, 1.0])
        .unwrap()
        .with_limits(4, 100);
    slice
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(state.iteration(), 1);
}

#[test]
fn invalid_slice_limits_fail_without_mutating_state() {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut kernel = SliceSampler::new(vec![1.0]).unwrap().with_limits(4, 0);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(104);

    let error = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap_err();

    assert!(matches!(error, McmcError::InvalidConfig(_)));
    assert_eq!(state.position(), &[0.0]);
    assert_eq!(state.log_density(), 0.0);
    assert_eq!(state.iteration(), 0);
}

#[test]
fn legacy_component_checkpoint_without_workspace_remains_usable() {
    let mut kernel: ComponentWiseMetropolis =
        serde_json::from_str(r#"{"scales":[0.5,0.5],"adaptations":[null,null]}"#).unwrap();
    let mut target = FnLogDensity::new(|position: &[f64]| {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.0, 0.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(105);

    kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();

    assert_eq!(state.iteration(), 1);
}

#[test]
fn legacy_slice_checkpoint_without_workspace_remains_usable() {
    let mut kernel: SliceSampler =
        serde_json::from_str(r#"{"widths":[1.0],"max_steps_out":4,"max_shrink_steps":100}"#)
            .unwrap();
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(106);

    kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();

    assert_eq!(state.iteration(), 1);
}
