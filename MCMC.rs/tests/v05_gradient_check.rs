use mcmc_rs::{check_gradient, DifferentiableLogDensity, GradientCheckConfig, LogDensity};

#[derive(Clone, Copy)]
struct CorrelatedGaussian;
impl LogDensity<[f64]> for CorrelatedGaussian {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * (2.0 * x[0] * x[0] + 2.0 * x[0] * x[1] + 3.0 * x[1] * x[1])
    }
}
impl DifferentiableLogDensity for CorrelatedGaussian {
    fn log_density_and_gradient(&mut self, x: &[f64], g: &mut [f64]) -> f64 {
        g[0] = -2.0 * x[0] - x[1];
        g[1] = -x[0] - 3.0 * x[1];
        self.log_density(x)
    }
}

#[test]
fn finite_difference_checker_accepts_correct_gradient() {
    let report = check_gradient(
        &mut CorrelatedGaussian,
        &[1.2, -0.7],
        GradientCheckConfig::default(),
    )
    .unwrap();
    assert!(report.passed, "{report:?}");
    assert_eq!(report.target_evaluations, 5);
    assert!(report.maximum_absolute_error < 1.0e-8);
}

struct WrongGradient;
impl LogDensity<[f64]> for WrongGradient {
    fn log_density(&mut self, x: &[f64]) -> f64 {
        -0.5 * x.iter().map(|v| v * v).sum::<f64>()
    }
}
impl DifferentiableLogDensity for WrongGradient {
    fn log_density_and_gradient(&mut self, x: &[f64], g: &mut [f64]) -> f64 {
        g.copy_from_slice(x);
        self.log_density(x)
    }
}

#[test]
fn finite_difference_checker_reports_wrong_components() {
    let report = check_gradient(
        &mut WrongGradient,
        &[0.75, -1.25],
        GradientCheckConfig::default(),
    )
    .unwrap();
    assert!(!report.passed);
    assert_eq!(report.components.iter().filter(|c| !c.passed).count(), 2);
}

#[test]
fn gradient_checker_rejects_invalid_input() {
    let invalid = GradientCheckConfig {
        step: 0.0,
        ..GradientCheckConfig::default()
    };
    assert!(check_gradient(&mut CorrelatedGaussian, &[0.0, 0.0], invalid).is_err());
    assert!(check_gradient(&mut CorrelatedGaussian, &[], GradientCheckConfig::default()).is_err());
}
