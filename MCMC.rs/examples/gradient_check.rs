use mcmc_rs::{
    check_gradient, DifferentiableLogDensity, GradientCheckConfig, LogDensity, McmcError,
};

struct Gaussian;
impl LogDensity<[f64]> for Gaussian {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * x.iter().map(|value| value * value).sum::<f64>()
    }
}
impl DifferentiableLogDensity for Gaussian {
    fn log_density_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> f64 {
        for (gradient, value) in gradient.iter_mut().zip(x.iter().copied()) {
            *gradient = -value;
        }
        self.log_density(x)
    }
}
fn main() -> Result<(), McmcError> {
    let report = check_gradient(&mut Gaussian, &[0.3, -1.2], GradientCheckConfig::default())?;
    println!(
        "passed={}, max_error={}",
        report.passed, report.maximum_absolute_error
    );
    Ok(())
}
