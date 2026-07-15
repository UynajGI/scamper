use carlo_rs::Estimate;

#[test]
fn test_estimate_from_empty_bins() {
    let bins: Vec<f64> = vec![];
    let est = Estimate::from_bins(&bins);

    assert_eq!(est.mean, 0.0);
    assert_eq!(est.stderr, 0.0);
    assert_eq!(est.n_bins, 0);
}

#[test]
fn test_estimate_from_single_bin() {
    let bins = vec![1.0];
    let est = Estimate::from_bins(&bins);

    assert!((est.mean - 1.0).abs() < 1e-10);
    assert_eq!(est.n_bins, 1);
    // Single bin has no error estimate
    assert_eq!(est.stderr, 0.0);
}

#[test]
fn test_estimate_from_two_bins() {
    let bins = vec![1.0, 3.0];
    let est = Estimate::from_bins(&bins);

    // Mean = (1 + 3) / 2 = 2
    assert!((est.mean - 2.0).abs() < 1e-10);
    assert_eq!(est.n_bins, 2);
    // Std of [1, 3] = sqrt(2)
    // Stderr = sqrt(2) / sqrt(2) = 1
    assert!((est.stderr - 1.0).abs() < 1e-10);
}

#[test]
fn test_estimate_from_constant_bins() {
    let bins = vec![5.0, 5.0, 5.0, 5.0];
    let est = Estimate::from_bins(&bins);

    assert!((est.mean - 5.0).abs() < 1e-10);
    // Constant values should have zero stderr
    assert!(est.stderr < 1e-10);
}

#[test]
fn test_estimate_from_large_bins() {
    let bins: Vec<f64> = (1..=1000).map(|i| i as f64).collect();
    let est = Estimate::from_bins(&bins);

    // Mean of 1..=1000 is 500.5
    assert!((est.mean - 500.5).abs() < 1e-10);
    assert_eq!(est.n_bins, 1000);
}

#[test]
fn test_estimate_format() {
    let est = Estimate {
        mean: 1.234567,
        stderr: 0.012345,
        autocorr_time: 1.0,
        n_bins: 100,
    };

    let formatted = est.format();
    assert!(formatted.contains("1.234567"));
    assert!(formatted.contains("0.012345"));
}

#[test]
fn test_estimate_negative_values() {
    let bins = vec![-1.0, -2.0, -3.0];
    let est = Estimate::from_bins(&bins);

    assert!((est.mean - (-2.0)).abs() < 1e-10);
}

#[test]
fn test_estimate_mixed_values() {
    let bins = vec![-1.0, 0.0, 1.0];
    let est = Estimate::from_bins(&bins);

    assert!((est.mean).abs() < 1e-10);
}

#[test]
fn test_estimate_serialization() {
    let est = Estimate {
        mean: 1.5,
        stderr: 0.1,
        autocorr_time: 2.0,
        n_bins: 100,
    };

    // Test JSON serialization
    let json = serde_json::to_string(&est).unwrap();
    assert!(json.contains("1.5"));
    assert!(json.contains("0.1"));

    // Test deserialization
    let est2: Estimate = serde_json::from_str(&json).unwrap();
    assert!((est2.mean - est.mean).abs() < 1e-10);
}
