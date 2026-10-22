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

#[test]
fn test_array_estimate_new_validates_element_count() {
    use carlo_rs::{ArrayEstimate, CarloError};

    let elements: Vec<Estimate> = (0..4)
        .map(|i| Estimate {
            mean: i as f64,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 10,
        })
        .collect();

    let arr = ArrayEstimate::new(vec![2, 2], elements.clone()).unwrap();
    assert_eq!(arr.shape(), &[2, 2]);
    assert_eq!(arr.n_elements(), 4);
    assert_eq!(arr.element(3).unwrap().mean, 3.0);
    // Row-major indexing: at([1, 0]) is flat index 2.
    assert_eq!(arr.at(&[1, 0]).unwrap().mean, 2.0);
    assert!(arr.at(&[2, 0]).is_none());
    assert!(arr.at(&[0]).is_none());

    // Wrong element count for the shape is an explicit error.
    let short: Vec<Estimate> = elements[..3].to_vec();
    match ArrayEstimate::new(vec![2, 2], short) {
        Err(CarloError::ShapeMismatch {
            expected_len,
            got_len,
            ..
        }) => {
            assert_eq!(expected_len, 4);
            assert_eq!(got_len, 3);
        }
        other => panic!("expected ShapeMismatch, got {other:?}"),
    }

    // Scalar shapes are rejected: scalars belong in a plain Estimate.
    assert!(ArrayEstimate::new(vec![], vec![Estimate::from_bins(&[1.0])]).is_err());
}

#[test]
fn test_array_estimate_from_element_bins() {
    use carlo_rs::ArrayEstimate;

    // Two elements with distinct bin columns.
    let columns = vec![vec![1.0, 2.0, 3.0, 4.0], vec![10.0, 12.0, 14.0, 16.0]];
    let arr = ArrayEstimate::from_element_bins(&[2], &columns).unwrap();

    let reference_a = Estimate::from_bins(&columns[0]);
    let reference_b = Estimate::from_bins(&columns[1]);
    let a = arr.element(0).unwrap();
    let b = arr.element(1).unwrap();
    assert!((a.mean - reference_a.mean).abs() < 1e-12);
    assert!((a.stderr - reference_a.stderr).abs() < 1e-12);
    assert!((b.mean - reference_b.mean).abs() < 1e-12);
    assert!((b.stderr - reference_b.stderr).abs() < 1e-12);
    assert_eq!(a.n_bins, 4);
    assert_eq!(b.n_bins, 4);

    // Ragged columns are rejected.
    let ragged = vec![vec![1.0, 2.0], vec![1.0, 2.0, 3.0]];
    assert!(ArrayEstimate::from_element_bins(&[2], &ragged).is_err());

    // Column count must match the shape product.
    assert!(ArrayEstimate::from_element_bins(&[3], &columns).is_err());
}

#[test]
fn test_array_estimate_json_round_trip() {
    use carlo_rs::ArrayEstimate;

    let columns = vec![vec![1.0, 2.0, 3.0], vec![7.0, 8.0, 9.0]];
    let arr = ArrayEstimate::from_element_bins(&[2], &columns).unwrap();

    let json = serde_json::to_string(&arr).unwrap();
    assert!(json.contains("\"shape\""));
    assert!(json.contains("\"elements\""));

    let back: ArrayEstimate = serde_json::from_str(&json).unwrap();
    assert_eq!(back.shape(), &[2]);
    for (lhs, rhs) in back.elements().iter().zip(arr.elements()) {
        assert!((lhs.mean - rhs.mean).abs() < 1e-15);
        assert!((lhs.stderr - rhs.stderr).abs() < 1e-15);
        assert_eq!(lhs.n_bins, rhs.n_bins);
    }
}
