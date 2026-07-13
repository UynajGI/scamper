use mcmc_rs::{run_multichain, DifferentiableLogDensity, LogDensity, McmcConfig, StaticHmc};

#[derive(Clone, Copy)]
struct CorrelatedGaussian;

impl LogDensity<[f64]> for CorrelatedGaussian {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        let x = position[0];
        let y = position[1];
        -(x * x - 1.6 * x * y + y * y) / 0.72
    }
}

impl DifferentiableLogDensity for CorrelatedGaussian {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        let x = position[0];
        let y = position[1];
        gradient[0] = -(2.0 * x - 1.6 * y) / 0.72;
        gradient[1] = -(-1.6 * x + 2.0 * y) / 0.72;
        self.log_density(position)
    }
}

fn main() -> Result<(), mcmc_rs::McmcError> {
    let warmup = 1_000;
    let config = McmcConfig {
        chains: 4,
        warmup,
        samples: 2_000,
        parameter_names: vec!["x".to_string(), "y".to_string()],
        ..McmcConfig::default()
    };
    let output = run_multichain(
        |_| CorrelatedGaussian,
        |_| {
            StaticHmc::diagonal(vec![1.0, 1.0], 0.15, 8)
                .and_then(|kernel| kernel.with_diagonal_adaptation(warmup, 0.8, 1.0e-3))
                .expect("valid HMC configuration")
        },
        vec![
            vec![-2.0, -1.0],
            vec![2.0, 1.0],
            vec![-1.0, 2.0],
            vec![1.0, -2.0],
        ],
        config,
    )?;

    for parameter in output.diagnostics.parameters {
        println!(
            "{}: mean={:.4}, rhat={:.4}, bulk_ess={:.1}",
            parameter.name, parameter.mean, parameter.rhat, parameter.ess_bulk
        );
    }
    println!("divergences={}", output.diagnostics.total_divergences);
    Ok(())
}
