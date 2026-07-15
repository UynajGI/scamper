use mcmc_rs::{
    DenseMetric, DiagonalMetric, DifferentiableLogDensity, EuclideanState, LogDensity, Metric,
    Nuts, SamplingPhase, StaticHmc, StepSizeSearch, TransitionKernel, UnitMetric,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Clone, Copy)]
struct StandardNormal;
impl LogDensity<[f64]> for StandardNormal {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * x.iter().map(|v| v * v).sum::<f64>()
    }
}
impl DifferentiableLogDensity for StandardNormal {
    fn log_density_and_gradient(&mut self, x: &[f64], g: &mut [f64]) -> f64 {
        for (g, x) in g.iter_mut().zip(x.iter().copied()) {
            *g = -x;
        }
        self.log_density(x)
    }
}

#[test]
fn hmc_and_nuts_search_rescue_bad_initial_scale() {
    for nuts in [false, true] {
        let mut target = StandardNormal;
        let mut state = EuclideanState::initialize(&mut target, vec![0.5]).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51E5 + u64::from(nuts));
        let report = if nuts {
            let mut k = Nuts::unit(1, 1.0e4, 7)
                .unwrap()
                .with_dual_averaging(40, 0.8)
                .unwrap()
                .with_step_size_search(StepSizeSearch::default())
                .unwrap();
            k.transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
                .unwrap()
        } else {
            let mut k = StaticHmc::unit(1, 1.0e4, 5)
                .unwrap()
                .with_dual_averaging(40, 0.8)
                .unwrap()
                .with_step_size_search(StepSizeSearch::default())
                .unwrap();
            k.transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
                .unwrap()
        };
        let used = report.proposal_scale.unwrap();
        assert!(used > 0.0 && used < 1.0e4);
        assert!(report.target_evaluations > report.leapfrog_steps);
    }
}

#[test]
fn search_requires_warmup_and_old_json_defaults_to_disabled() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut k = Nuts::unit(1, 0.1, 5)
        .unwrap()
        .with_step_size_search(StepSizeSearch::default())
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(12);
    assert!(k
        .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
        .is_err());

    let mut json = serde_json::to_value(Nuts::unit(1, 0.2, 5).unwrap()).unwrap();
    json.as_object_mut().unwrap().remove("step_size_search");
    json.as_object_mut()
        .unwrap()
        .remove("step_size_search_complete");
    let _: Nuts<UnitMetric> = serde_json::from_value(json).unwrap();

    let mut hmc_json = serde_json::to_value(StaticHmc::unit(1, 0.2, 3).unwrap()).unwrap();
    hmc_json.as_object_mut().unwrap().remove("step_size_search");
    hmc_json
        .as_object_mut()
        .unwrap()
        .remove("step_size_search_complete");
    let _: StaticHmc<UnitMetric> = serde_json::from_value(hmc_json).unwrap();
}

#[test]
fn built_in_metric_summed_products_match_explicit_velocity() {
    fn check<M: Metric>(m: &M) {
        let p = [0.4, -0.7];
        let a = [1.2, 0.3];
        let b = [-0.2, 0.9];
        let mut v = [0.0; 2];
        m.velocity(&p, &mut v).unwrap();
        let expected = v
            .iter()
            .zip(a.iter().zip(b.iter()))
            .map(|(v, (a, b))| v * (a + b))
            .sum::<f64>();
        let actual = m.velocity_dot_momentum_sum(&p, &a, &b).unwrap();
        assert!((actual - expected).abs() < 1.0e-12);
    }
    check(&UnitMetric::new(2).unwrap());
    check(&DiagonalMetric::new(vec![2.0, 0.5]).unwrap());
    check(&DenseMetric::from_inverse_mass(2, &[1.5, 0.3, 0.3, 0.8], 1e-12).unwrap());
}

#[test]
fn search_repeats_after_metric_updates() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![1.0]).unwrap();
    let mut k = StaticHmc::new(DiagonalMetric::unit(1).unwrap(), 20.0, 3)
        .unwrap()
        .with_diagonal_adaptation(80, 0.8, 1e-3)
        .unwrap()
        .with_step_size_search(StepSizeSearch::default())
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xA0A0);
    let mut searches = 0;
    for _ in 0..80 {
        let r = k
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
        if r.target_evaluations > r.leapfrog_steps {
            searches += 1;
        }
    }
    assert!(searches >= 2);
    assert!(k.adaptation_is_frozen());
}

#[derive(Clone, Copy)]
struct Correlated;
impl LogDensity<[f64]> for Correlated {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        let r = 0.8;
        -0.5 * (x[0] * x[0] - 2.0 * r * x[0] * x[1] + x[1] * x[1]) / (1.0 - r * r)
    }
}
impl DifferentiableLogDensity for Correlated {
    fn log_density_and_gradient(&mut self, x: &[f64], g: &mut [f64]) -> f64 {
        let r = 0.8;
        let d = 1.0 - r * r;
        g[0] = -(x[0] - r * x[1]) / d;
        g[1] = -(x[1] - r * x[0]) / d;
        self.log_density(x)
    }
}

#[test]
fn generalized_nuts_recovers_correlated_covariance() {
    let mut target = Correlated;
    let mut state = EuclideanState::initialize(&mut target, vec![1.5, -1.0]).unwrap();
    let metric = DenseMetric::from_inverse_mass(2, &[1.0, 0.8, 0.8, 1.0], 1e-12).unwrap();
    let mut k = Nuts::new(metric, 0.3, 8).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xD1A6);
    for _ in 0..300 {
        k.transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
    }
    let n = 6000.0;
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for _ in 0..6000 {
        k.transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        let x = state.position()[0];
        let y = state.position()[1];
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
    }
    let (mx, my) = (sx / n, sy / n);
    let (vx, vy) = (sxx / n - mx * mx, syy / n - my * my);
    let cov = sxy / n - mx * my;
    assert!(mx.abs() < 0.08 && my.abs() < 0.08);
    assert!((vx - 1.0).abs() < 0.15 && (vy - 1.0).abs() < 0.15);
    assert!((cov - 0.8).abs() < 0.15, "cov={cov}");
}

#[test]
fn step_size_search_checkpoint_preserves_future_trajectory() {
    let mut target = StandardNormal;
    let mut state = EuclideanState::initialize(&mut target, vec![1.0]).unwrap();
    let mut kernel = Nuts::unit(1, 50.0, 7)
        .unwrap()
        .with_dual_averaging(30, 0.8)
        .unwrap()
        .with_step_size_search(StepSizeSearch::default())
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xC0FFEE);
    for _ in 0..10 {
        kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }
    let encoded = serde_json::to_string(&(kernel, state, rng)).unwrap();
    let (mut left_kernel, mut left_state, mut left_rng): (
        Nuts<UnitMetric>,
        EuclideanState,
        Xoshiro256PlusPlus,
    ) = serde_json::from_str(&encoded).unwrap();
    let (mut right_kernel, mut right_state, mut right_rng): (
        Nuts<UnitMetric>,
        EuclideanState,
        Xoshiro256PlusPlus,
    ) = serde_json::from_str(&encoded).unwrap();
    let mut left_target = StandardNormal;
    let mut right_target = StandardNormal;
    for index in 0..40 {
        let phase = if index < 20 {
            SamplingPhase::Warmup
        } else {
            SamplingPhase::Sampling
        };
        let left = left_kernel
            .transition(&mut left_target, &mut left_state, &mut left_rng, phase)
            .unwrap();
        let right = right_kernel
            .transition(&mut right_target, &mut right_state, &mut right_rng, phase)
            .unwrap();
        assert_eq!(left, right);
        assert_eq!(left_state.position(), right_state.position());
        assert_eq!(left_state.log_density(), right_state.log_density());
    }
}
