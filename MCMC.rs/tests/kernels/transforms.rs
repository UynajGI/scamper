use mcmc_rs::{
    Bijector, Interval, LogDensity, Ordered, Positive, Product, Simplex, TransformedTarget,
};

fn assert_close(left: f64, right: f64) {
    assert!((left - right).abs() < 1.0e-10, "left={left}, right={right}");
}

#[test]
fn scalar_and_product_transforms_round_trip_with_inverse_jacobians() {
    let mut transform = Product::new(Positive, Interval::new(-2.0, 5.0).unwrap());
    let input = [0.3, -0.7];
    let mut constrained = [0.0; 2];
    let forward = transform.forward(&input, &mut constrained).unwrap();
    assert!(constrained[0] > 0.0);
    assert!(constrained[1] > -2.0 && constrained[1] < 5.0);

    let mut recovered = [0.0; 2];
    let inverse = transform.inverse(&constrained, &mut recovered).unwrap();
    assert_close(recovered[0], input[0]);
    assert_close(recovered[1], input[1]);
    assert_close(forward + inverse, 0.0);
}

#[test]
fn ordered_and_simplex_transforms_round_trip() {
    let mut ordered = Ordered::new(3).unwrap();
    let input = [-1.0, 0.2, -0.4];
    let mut constrained = [0.0; 3];
    let forward = ordered.forward(&input, &mut constrained).unwrap();
    assert!(constrained[0] < constrained[1] && constrained[1] < constrained[2]);
    let mut recovered = [0.0; 3];
    let inverse = ordered.inverse(&constrained, &mut recovered).unwrap();
    for (recovered, expected) in recovered.iter().zip(input) {
        assert_close(*recovered, expected);
    }
    assert_close(forward + inverse, 0.0);

    let mut simplex = Simplex::new(3).unwrap();
    let input = [0.4, -0.8];
    let mut constrained = [0.0; 3];
    let forward = simplex.forward(&input, &mut constrained).unwrap();
    assert_close(constrained.iter().sum(), 1.0);
    assert!(constrained.iter().all(|value| *value > 0.0));
    let mut recovered = [0.0; 2];
    let inverse = simplex.inverse(&constrained, &mut recovered).unwrap();
    for (recovered, expected) in recovered.iter().zip(input) {
        assert_close(*recovered, expected);
    }
    assert_close(forward + inverse, 0.0);
}

struct ExponentialTarget;

impl LogDensity<[f64]> for ExponentialTarget {
    fn log_density(&mut self, state: &[f64]) -> f64 {
        if state[0] > 0.0 {
            -state[0]
        } else {
            f64::NEG_INFINITY
        }
    }
}

#[test]
fn transformed_target_adds_log_jacobian() {
    let mut target = TransformedTarget::new(ExponentialTarget, Positive).unwrap();
    let unconstrained = [2.0_f64.ln()];
    assert_close(target.log_density(&unconstrained), -2.0 + 2.0_f64.ln());
    assert_close(target.constrained_position()[0], 2.0);
}
