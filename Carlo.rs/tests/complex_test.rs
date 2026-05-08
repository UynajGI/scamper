use carlo_rs::{
    Accumulator, ComplexAccumulator, ComplexEstimate, ComplexResult, Context, Estimate,
    Measurements, Results,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_complex_accumulator_scalar() {
    let mut acc = ComplexAccumulator::new(1);

    // Add 4 complex samples: (1+2i), (3+4i), (5+6i), (7+8i)
    for (re, im) in [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (7.0, 8.0)] {
        acc.add(re, im);
    }

    assert_eq!(acc.num_bins(), 4);
    assert_eq!(acc.total_count(), 4);

    let est = acc.finalize();
    // Mean of re: (1+3+5+7)/4 = 4.0
    assert!(
        (est.re.mean - 4.0).abs() < 1e-10,
        "Expected re mean 4.0, got {}",
        est.re.mean
    );
    // Mean of im: (2+4+6+8)/4 = 5.0
    assert!(
        (est.im.mean - 5.0).abs() < 1e-10,
        "Expected im mean 5.0, got {}",
        est.im.mean
    );
}

#[test]
fn test_complex_accumulator_binning() {
    let mut acc = ComplexAccumulator::new(5); // bin size 5

    // Add 10 samples
    for i in 0..10 {
        acc.add(i as f64, (i * 2) as f64);
    }

    // Should produce 2 bins (10 / 5 = 2)
    assert_eq!(acc.num_bins(), 2);

    let est = acc.finalize();
    // Bin 0: mean of re [0,1,2,3,4] = 2.0, im [0,2,4,6,8] = 4.0
    // Bin 1: mean of re [5,6,7,8,9] = 7.0, im [10,12,14,16,18] = 14.0
    // Overall mean re: (2+7)/2 = 4.5
    // Overall mean im: (4+14)/2 = 9.0
    assert!((est.re.mean - 4.5).abs() < 1e-10);
    assert!((est.im.mean - 9.0).abs() < 1e-10);
}

#[test]
fn test_complex_estimate_format() {
    let est = ComplexEstimate::new(
        Estimate {
            mean: 1.5,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 10,
        },
        Estimate {
            mean: 2.5,
            stderr: 0.2,
            autocorr_time: 1.5,
            n_bins: 10,
        },
    );
    let formatted = est.format();
    assert!(formatted.contains("1.500000"));
    assert!(formatted.contains("2.500000"));
    assert!(formatted.contains("i"));
}

#[test]
fn test_measurements_complex() {
    let mut meas = Measurements::new(1);

    // Add complex samples
    meas.add_sample_complex("phase", 0.5, 0.866);
    meas.add_sample_complex("phase", 0.3, 0.954);
    meas.add_sample_complex("phase", 0.7, 0.714);

    let complex_results = meas.finalize_complex();
    assert!(complex_results.contains_key("phase"));
    assert_eq!(complex_results["phase"].re.n_bins, 3);
    assert_eq!(complex_results["phase"].im.n_bins, 3);
}

#[test]
fn test_context_measure_complex() {
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 10);

    ctx.measure_complex("complex_obs", 1.0, 2.0);
    ctx.measure_complex("complex_obs", 3.0, 4.0);
    ctx.advance_sweep();

    let complex_results = ctx.finalize_complex_measurements();
    assert!(complex_results.contains_key("complex_obs"));
}

#[test]
fn test_results_complex_json() {
    let mut results = Results::new();

    results.add("energy", Estimate::from_bins(&[1.0, 2.0, 3.0]));

    let complex_est = ComplexEstimate::new(
        Estimate::from_bins(&[0.1, 0.2, 0.3]),
        Estimate::from_bins(&[0.4, 0.5, 0.6]),
    );
    results.add_complex("order_parameter", complex_est);

    let json = results.to_json().unwrap();
    assert!(json.contains("complex_observables"));
    assert!(json.contains("order_parameter"));
    assert!(json.contains("\"re\""));
    assert!(json.contains("\"im\""));
}

#[test]
fn test_complex_result_roundtrip() {
    let original = ComplexResult {
        re: Estimate {
            mean: 1.5,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 10,
        },
        im: Estimate {
            mean: 2.5,
            stderr: 0.2,
            autocorr_time: 1.5,
            n_bins: 10,
        },
    };

    let json = serde_json::to_string(&original).unwrap();
    let restored: ComplexResult = serde_json::from_str(&json).unwrap();

    assert!((restored.re.mean - original.re.mean).abs() < 1e-10);
    assert!((restored.im.mean - original.im.mean).abs() < 1e-10);
    assert!((restored.re.stderr - original.re.stderr).abs() < 1e-10);
    assert!((restored.im.stderr - original.im.stderr).abs() < 1e-10);
}

#[test]
fn test_complex_estimate_magnitude() {
    let est = ComplexEstimate::new(
        Estimate {
            mean: 3.0,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 10,
        },
        Estimate {
            mean: 4.0,
            stderr: 0.2,
            autocorr_time: 1.5,
            n_bins: 10,
        },
    );
    // magnitude should be sqrt(3^2 + 4^2) = 5.0
    assert!((est.magnitude() - 5.0).abs() < 1e-10);
}

#[test]
fn test_accumulator_autocorr_time_from_bins() {
    // Correlated data: each value = previous + small perturbation
    let mut acc = Accumulator::new(1);
    let mut value = 0.0;
    for _ in 0..100 {
        acc.add(value);
        value += (value * 0.9).min(5.0); // autocorrelated random walk
    }
    let autocorr = acc.autocorr_time_from_bins();
    assert!(
        autocorr >= 0.0,
        "Autocorrelation should be non-negative, got {autocorr}"
    );
}
