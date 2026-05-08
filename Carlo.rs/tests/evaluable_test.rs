use carlo_rs::evaluable::jackknife;
use carlo_rs::evaluable::{Evaluator, MultiplexEvaluator};
use carlo_rs::merge::ResultObservable;
use carlo_rs::run::timing;
use ndarray::Array1;
use std::collections::HashMap;

#[test]
fn test_jackknife_simple_mean() {
    let samples = Array1::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]).into_dyn();
    let sample_sets = vec![samples];

    let (mean, error, _cov) = jackknife(
        |args| args[0].clone(), // Identity function
        &sample_sets,
        false,
    )
    .unwrap();

    // Mean should be close to 3.0
    assert!((mean[0] - 3.0).abs() < 0.5);
    assert!(error[0] > 0.0);
}

#[test]
fn test_evaluator_creation() {
    let obs = ResultObservable::<f64> {
        internal_bin_length: 100,
        rebin_length: 500,
        mean: Array1::from_vec(vec![1.0]).into_dyn(),
        error: Array1::from_vec(vec![0.1]).into_dyn(),
        covariance: None,
        autocorrelation_time: Array1::from_vec(vec![5.0]).into_dyn(),
        rebin_means: Array1::from_vec(vec![1.0; 10]).into_dyn(),
    };

    let observables = HashMap::from([("Energy".to_string(), obs)]);
    let evaluator = Evaluator::new(observables, false);
    assert!(evaluator.observables().contains_key("Energy"));
}

#[test]
fn test_timing_constants() {
    // Verify timing observable names match Carlo.jl convention
    assert_eq!(timing::SWEEP_TIME, "_ll_sweep_time");
    assert_eq!(timing::MEASURE_TIME, "_ll_measure_time");
    assert_eq!(timing::CHECKPOINT_READ_TIME, "_ll_checkpoint_read_time");
    assert_eq!(timing::CHECKPOINT_WRITE_TIME, "_ll_checkpoint_write_time");
}

#[test]
fn test_multiplex_evaluator_creation() {
    let multi_eval = MultiplexEvaluator::new(4); // 4 PT chains
    assert_eq!(multi_eval.entry_count(), 4);
    assert_eq!(multi_eval.len(), 0);
    assert!(multi_eval.is_empty());
}

#[test]
fn test_multiplex_evaluator_register() {
    let mut multi_eval = MultiplexEvaluator::new(2);

    // Register evaluator for chain 0
    multi_eval.evaluate("energy", &["magnetization"], |args| args[0].clone() * 2.0);

    // Register evaluator for chain 1
    multi_eval.evaluate("energy", &["magnetization"], |args| args[0].clone() * 3.0);

    assert_eq!(multi_eval.len(), 1);
    assert!(!multi_eval.is_empty());
}

#[test]
fn test_evaluator_evaluate() {
    // Create mock observables
    let mut observables: HashMap<String, ResultObservable<f64>> = HashMap::new();

    // Create a ResultObservable with some samples
    let samples = ndarray::Array1::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]).into_dyn();
    let mean_val = samples.mean().unwrap_or(0.0);
    let obs = ResultObservable {
        internal_bin_length: 1,
        rebin_length: 1,
        mean: ndarray::Array0::from_elem((), mean_val).into_dyn(),
        error: ndarray::Array0::from_elem((), 0.0).into_dyn(),
        covariance: None,
        autocorrelation_time: ndarray::ArrayD::zeros(vec![]),
        rebin_means: samples.clone(),
    };

    observables.insert("magnetization".to_string(), obs);

    let mut evaluator = Evaluator::new(observables, false);

    // Define a simple evaluable: energy = 2 * magnetization
    evaluator
        .evaluate("energy", &["magnetization"], |args| args[0].clone() * 2.0)
        .unwrap();

    // Check that the evaluable was created
    assert!(evaluator.evaluables().contains_key("energy"));

    // The mean should be 2 * mean(magnetization) = 2 * 3.0 = 6.0
    let energy = evaluator.evaluables().get("energy").unwrap();
    let expected = 2.0 * (1.0 + 2.0 + 3.0 + 4.0 + 5.0) / 5.0;
    assert!(
        (energy[[0]] - expected).abs() < 1e-10,
        "Expected {}, got {}",
        expected,
        energy[[0]]
    );
}

#[test]
fn test_evaluator_missing_ingredient() {
    let observables: HashMap<String, ResultObservable<f64>> = HashMap::new();
    let mut evaluator = Evaluator::new(observables, false);

    // Try to define an evaluable with missing ingredient
    let result = evaluator.evaluate("energy", &["nonexistent"], |_args| {
        ndarray::Array0::from_elem((), 0.0).into_dyn()
    });

    // Should succeed but not add the evaluable
    assert!(result.is_ok());
    assert!(!evaluator.evaluables().contains_key("energy"));
}
