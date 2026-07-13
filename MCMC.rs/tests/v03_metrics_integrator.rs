use mcmc_rs::{
    DenseMetric, DiagonalMetric, DifferentiableLogDensity, LeapfrogIntegrator, LogDensity, Metric,
    PhasePoint, UnitMetric,
};

#[derive(Clone, Copy)]
struct StandardNormal;

impl LogDensity<[f64]> for StandardNormal {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -0.5 * position.iter().map(|value| value * value).sum::<f64>()
    }
}

impl DifferentiableLogDensity for StandardNormal {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        for (gradient, position) in gradient.iter_mut().zip(position.iter().copied()) {
            *gradient = -position;
        }
        self.log_density(position)
    }
}

#[test]
fn diagonal_and_dense_metrics_apply_inverse_mass_consistently() {
    let diagonal = DiagonalMetric::new(vec![2.0, 0.5]).unwrap();
    let mut velocity = vec![0.0; 2];
    diagonal.velocity(&[1.0, 2.0], &mut velocity).unwrap();
    assert_eq!(velocity, vec![2.0, 1.0]);
    assert!((diagonal.kinetic_energy(&[1.0, 2.0]).unwrap() - 2.0).abs() < 1.0e-12);

    let dense = DenseMetric::from_inverse_mass(2, &[2.0, 0.5, 0.5, 1.0], 1.0e-12).unwrap();
    dense.velocity(&[1.0, 2.0], &mut velocity).unwrap();
    assert!((velocity[0] - 3.0).abs() < 1.0e-10);
    assert!((velocity[1] - 2.5).abs() < 1.0e-10);
    assert!((dense.kinetic_energy(&[1.0, 2.0]).unwrap() - 4.0).abs() < 1.0e-10);
}

#[test]
fn leapfrog_is_reversible_up_to_roundoff() {
    let metric = UnitMetric::new(1).unwrap();
    let mut target = StandardNormal;
    let mut integrator = LeapfrogIntegrator::with_dimension(1);
    let mut point = PhasePoint {
        position: vec![0.7],
        momentum: vec![-0.4],
        gradient: vec![-0.7],
        log_density: -0.5 * 0.7_f64.powi(2),
    };
    let initial = point.clone();
    integrator
        .integrate(&mut target, &metric, &mut point, 0.05, 20)
        .unwrap();
    for momentum in &mut point.momentum {
        *momentum = -*momentum;
    }
    integrator
        .integrate(&mut target, &metric, &mut point, 0.05, 20)
        .unwrap();
    assert!((point.position[0] - initial.position[0]).abs() < 1.0e-12);
    assert!((point.momentum[0] - -initial.momentum[0]).abs() < 1.0e-12);
}

#[test]
fn smaller_leapfrog_steps_reduce_energy_error() {
    let metric = UnitMetric::new(1).unwrap();
    let mut target = StandardNormal;
    let initial = PhasePoint {
        position: vec![1.1],
        momentum: vec![0.8],
        gradient: vec![-1.1],
        log_density: -0.5 * 1.1_f64.powi(2),
    };
    let energy =
        |point: &PhasePoint| -point.log_density + metric.kinetic_energy(&point.momentum).unwrap();

    let mut coarse = initial.clone();
    LeapfrogIntegrator::with_dimension(1)
        .integrate(&mut target, &metric, &mut coarse, 0.2, 5)
        .unwrap();
    let mut fine = initial.clone();
    LeapfrogIntegrator::with_dimension(1)
        .integrate(&mut target, &metric, &mut fine, 0.1, 10)
        .unwrap();
    let coarse_error = (energy(&coarse) - energy(&initial)).abs();
    let fine_error = (energy(&fine) - energy(&initial)).abs();
    assert!(
        fine_error < coarse_error,
        "fine={fine_error}, coarse={coarse_error}"
    );
}
