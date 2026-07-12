use mcmc_rs::{run_multichain, FnLogDensity, McmcConfig, RandomWalkMetropolis};

fn main() -> Result<(), mcmc_rs::McmcError> {
    let kernel = RandomWalkMetropolis::isotropic(2, 0.5)?
        .with_scale_adaptation(0.234)?
        .with_diagonal_covariance_adaptation(1e-3)?;
    let output = run_multichain(
        |_| {
            FnLogDensity::new(|position: &[f64]| {
                let x = position[0];
                let y = position[1];
                -0.5 * (x * x + (y - 0.8 * x).powi(2) / 0.36)
            })
        },
        |_| kernel.clone(),
        vec![
            vec![-2.0, -2.0],
            vec![-1.0, 2.0],
            vec![1.0, -2.0],
            vec![2.0, 2.0],
        ],
        McmcConfig {
            warmup: 2_000,
            samples: 5_000,
            parameter_names: vec!["x".to_string(), "y".to_string()],
            ..McmcConfig::default()
        },
    )?;
    for parameter in output.diagnostics.parameters {
        println!(
            "{}: mean={:.3}, sd={:.3}, rhat={:.4}, ess_bulk={:.0}",
            parameter.name, parameter.mean, parameter.std_dev, parameter.rhat, parameter.ess_bulk
        );
    }
    Ok(())
}
