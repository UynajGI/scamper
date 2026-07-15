use mcmc_rs::{DifferentiableLogDensity, LogDensity, Positive, Simplex, TransformedTarget};

#[derive(Clone, Copy)]
struct PositiveGaussian;

impl LogDensity<[f64]> for PositiveGaussian {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -0.5 * position[0] * position[0]
    }
}

impl DifferentiableLogDensity for PositiveGaussian {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        gradient[0] = -position[0];
        self.log_density(position)
    }
}

#[derive(Clone, Copy)]
struct SimplexQuadratic;

impl LogDensity<[f64]> for SimplexQuadratic {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -0.5 * position
            .iter()
            .enumerate()
            .map(|(index, value)| (index + 1) as f64 * value * value)
            .sum::<f64>()
    }
}

impl DifferentiableLogDensity for SimplexQuadratic {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        for (index, (gradient, value)) in gradient
            .iter_mut()
            .zip(position.iter().copied())
            .enumerate()
        {
            *gradient = -((index + 1) as f64) * value;
        }
        self.log_density(position)
    }
}

#[test]
fn positive_transform_gradient_includes_pullback_and_jacobian() {
    let mut target = TransformedTarget::new(PositiveGaussian, Positive).unwrap();
    let mut gradient = vec![0.0];
    let z = [0.3];
    let log_density = target.log_density_and_gradient(&z, &mut gradient);
    let x = z[0].exp();
    assert!((log_density - (-0.5 * x * x + z[0])).abs() < 1.0e-12);
    assert!((gradient[0] - (1.0 - x * x)).abs() < 1.0e-12);
}

#[test]
fn simplex_transformed_gradient_matches_finite_difference() {
    let mut target = TransformedTarget::new(SimplexQuadratic, Simplex::new(4).unwrap()).unwrap();
    let position = vec![-0.7, 0.4, 1.1];
    let mut gradient = vec![0.0; 3];
    let value = target.log_density_and_gradient(&position, &mut gradient);
    assert!(value.is_finite());

    let step = 1.0e-6;
    for index in 0..position.len() {
        let mut left = position.clone();
        let mut right = position.clone();
        left[index] -= step;
        right[index] += step;
        let finite_difference =
            (target.log_density(&right) - target.log_density(&left)) / (2.0 * step);
        assert!(
            (gradient[index] - finite_difference).abs() < 2.0e-6,
            "index={index}, analytic={}, finite_difference={finite_difference}",
            gradient[index]
        );
    }
}
