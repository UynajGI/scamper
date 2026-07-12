use mcmc_rs::{
    EuclideanState, FnLogDensity, RandomWalkMetropolis, SamplingPhase, TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn proposal_scale_is_constant_in_sampling_phase() {
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2));
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut kernel = RandomWalkMetropolis::isotropic(1, 1.0)
        .unwrap()
        .with_scale_adaptation(0.44)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(9);
    for _ in 0..500 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    kernel.on_phase_end(SamplingPhase::Warmup, &state).unwrap();
    let frozen = kernel.effective_global_multiplier();
    for _ in 0..500 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        assert_eq!(kernel.effective_global_multiplier(), frozen);
    }
}
