use mcmc_rs::{run_multichain, FnLogDensity, McmcConfig, RandomWalkMetropolis};

#[test]
fn fixed_seed_is_reproducible_across_runs() {
    let config = McmcConfig {
        chains: 4,
        warmup: 200,
        samples: 500,
        base_seed: 123,
        parameter_names: vec!["x".to_string()],
        ..McmcConfig::default()
    };
    let initials = vec![vec![-2.0], vec![-1.0], vec![1.0], vec![2.0]];
    let run = || {
        run_multichain(
            |_| FnLogDensity::new(|position: &[f64]| -0.5 * position[0].powi(2)),
            |_| {
                RandomWalkMetropolis::isotropic(1, 0.8)
                    .unwrap()
                    .with_scale_adaptation(0.44)
                    .unwrap()
            },
            initials.clone(),
            config.clone(),
        )
        .unwrap()
    };
    let first = run();
    let second = run();
    for (left, right) in first.chains.iter().zip(&second.chains) {
        assert_eq!(left.trace, right.trace);
        assert_eq!(left.final_position, right.final_position);
        assert_eq!(left.final_log_density, right.final_log_density);
    }
}
