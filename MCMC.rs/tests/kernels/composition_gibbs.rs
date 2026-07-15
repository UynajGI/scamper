use mcmc_rs::proposal::standard_normal;
use mcmc_rs::{
    EuclideanState, GibbsKernel, GibbsUpdate, GibbsUpdateResult, LogDensity, McmcError, Mixture,
    Repeat, SamplingPhase, Then, TransitionKernel, TransitionReport,
};
use rand::Rng;
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

struct ConditionalNormal {
    conditional_mean: f64,
}

impl LogDensity<[f64]> for ConditionalNormal {
    fn log_density(&mut self, state: &[f64]) -> f64 {
        -0.5 * (state[0] - self.conditional_mean).powi(2)
    }
}

#[derive(Clone, Copy)]
struct DrawConditionalNormal;

impl GibbsUpdate<ConditionalNormal> for DrawConditionalNormal {
    fn update<R>(
        &mut self,
        target: &mut ConditionalNormal,
        _current: &EuclideanState,
        proposed_position: &mut [f64],
        rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<GibbsUpdateResult, McmcError>
    where
        R: Rng + ?Sized,
    {
        proposed_position[0] = target.conditional_mean + standard_normal(rng);
        Ok(GibbsUpdateResult::requiring_target_evaluation())
    }

    fn name(&self, _target: &ConditionalNormal) -> &'static str {
        "DrawConditionalNormal"
    }
}

#[derive(Clone, Copy)]
struct EmptyReportKernel;

impl TransitionKernel<ConditionalNormal> for EmptyReportKernel {
    fn transition<R>(
        &mut self,
        _target: &mut ConditionalNormal,
        state: &mut EuclideanState,
        _rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        state.mark_rejected_transition();
        Ok(TransitionReport::default())
    }

    fn name(&self, _target: &ConditionalNormal) -> &'static str {
        "EmptyReportKernel"
    }
}

struct FailingUpdate;

impl GibbsUpdate<ConditionalNormal> for FailingUpdate {
    fn update<R>(
        &mut self,
        _target: &mut ConditionalNormal,
        _current: &EuclideanState,
        proposed_position: &mut [f64],
        _rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<GibbsUpdateResult, McmcError>
    where
        R: Rng + ?Sized,
    {
        proposed_position[0] = 99.0;
        Err(McmcError::InvalidConfig("intentional failure".to_string()))
    }
}

#[test]
fn target_specific_gibbs_update_is_atomic_and_exact() {
    let mut target = ConditionalNormal {
        conditional_mean: 2.5,
    };
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut kernel = GibbsKernel::new(DrawConditionalNormal);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert!(state.position()[0].is_finite());
    assert_eq!(
        state.log_density(),
        -0.5 * (state.position()[0] - target.conditional_mean).powi(2)
    );
    assert_eq!(state.iteration(), 1);
    assert_eq!(report.target_evaluations, 1);
    assert_eq!(report.subtransitions, 1);
}

#[test]
fn failing_gibbs_update_does_not_corrupt_accepted_state() {
    let mut target = ConditionalNormal {
        conditional_mean: 0.0,
    };
    let mut state = EuclideanState::initialize(&mut target, vec![1.0]).unwrap();
    let before = state.position().clone();
    let before_log_density = state.log_density();
    let mut kernel = GibbsKernel::new(FailingUpdate);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(9);

    assert!(kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .is_err());
    assert_eq!(state.position(), &before);
    assert_eq!(state.log_density(), before_log_density);
    assert_eq!(state.iteration(), 0);
}

#[test]
fn static_composition_aggregates_reports() {
    let mut target = ConditionalNormal {
        conditional_mean: -1.0,
    };
    let mut state = EuclideanState::initialize(&mut target, vec![4.0]).unwrap();
    let first = GibbsKernel::new(DrawConditionalNormal);
    let second = GibbsKernel::new(DrawConditionalNormal);
    let mut kernel = Then::new(first, second);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(report.proposals, 2);
    assert_eq!(report.acceptances, 2);
    assert_eq!(report.target_evaluations, 2);
    assert_eq!(report.subtransitions, 2);
    assert_eq!(report.accepted, None);
    assert_eq!(state.iteration(), 2);
}

#[test]
fn repeat_and_mixture_are_deterministic_at_probability_boundaries() {
    let mut target = ConditionalNormal {
        conditional_mean: 3.0,
    };
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let repeated = Repeat::new(GibbsKernel::new(DrawConditionalNormal), 3).unwrap();
    let fallback = GibbsKernel::new(FailingUpdate);
    let mut kernel = Mixture::new(repeated, fallback, 1.0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(13);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(report.subtransitions, 3);
    assert_eq!(report.proposals, 3);
    assert!(state.position()[0].is_finite());
}

#[test]
fn composition_normalizes_legacy_or_empty_child_reports() {
    let mut target = ConditionalNormal {
        conditional_mean: 0.0,
    };
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut kernel = Then::new(EmptyReportKernel, EmptyReportKernel);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(15);

    let report = kernel
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .unwrap();
    assert_eq!(report.subtransitions, 2);
    assert_eq!(state.iteration(), 2);
    assert_eq!(kernel.name(&target), "Then");
}
