use carlo_rs::Accumulator;

#[test]
fn test_accumulator_binning() {
    let bin_length = 10;
    let mut acc = Accumulator::new(bin_length);

    // Fill first bin
    for i in 1..=bin_length {
        acc.add(i as f64);
    }

    assert!(acc.has_complete_bins());

    // Mean of 1..=10 is 5.5
    let est = acc.finalize();
    assert!((est.mean - 5.5).abs() < 1e-10);
}

#[test]
fn test_accumulator_multiple_bins() {
    let bin_length = 5;
    let mut acc = Accumulator::new(bin_length);

    // Add 15 samples = 3 full bins
    for i in 1..=15 {
        acc.add(i as f64);
    }

    let est = acc.finalize();
    assert!(est.mean > 0.0);
    assert!(est.stderr > 0.0);
    assert_eq!(est.n_bins, 3);
}

#[test]
fn test_accumulator_partial_bin() {
    let bin_length = 10;
    let mut acc = Accumulator::new(bin_length);

    // Add 15 samples = 1 full bin + 5 partial
    for i in 1..=15 {
        acc.add(i as f64);
    }

    let est = acc.finalize();
    // Should include partial bin
    assert_eq!(est.n_bins, 2);
}

#[test]
fn test_accumulator_empty() {
    let bin_length = 10;
    let acc = Accumulator::new(bin_length);

    assert!(!acc.has_complete_bins());

    let est = acc.finalize();
    assert_eq!(est.n_bins, 0);
    assert_eq!(est.mean, 0.0);
}

#[test]
fn test_accumulator_single_sample() {
    let bin_length = 10;
    let mut acc = Accumulator::new(bin_length);

    acc.add(1.0);

    assert!(!acc.has_complete_bins());

    let est = acc.finalize();
    assert_eq!(est.n_bins, 1);
    assert!((est.mean - 1.0).abs() < 1e-10);
}

#[test]
fn test_accumulator_bin_capacity() {
    let bin_length = 5;
    let acc = Accumulator::new(bin_length);

    assert_eq!(acc.bin_capacity(), 5);
}

#[test]
fn test_accumulator_bins_access() {
    let bin_length = 3;
    let mut acc = Accumulator::new(bin_length);

    for i in 1..=6 {
        acc.add(i as f64);
    }

    let bins = acc.bins();
    assert_eq!(bins.len(), 2);
}

#[test]
fn test_accumulator_zero_values() {
    let bin_length = 5;
    let mut acc = Accumulator::new(bin_length);

    for _ in 0..5 {
        acc.add(0.0);
    }

    let est = acc.finalize();
    assert!((est.mean).abs() < 1e-10);
}

#[test]
fn test_accumulator_negative_values() {
    let bin_length = 3;
    let mut acc = Accumulator::new(bin_length);

    for i in -3..=0 {
        acc.add(i as f64);
    }

    let est = acc.finalize();
    assert!(est.mean < 0.0);
}

#[test]
fn test_accumulator_large_values() {
    let bin_length = 2;
    let mut acc = Accumulator::new(bin_length);

    acc.add(1e10);
    acc.add(1e10);

    let est = acc.finalize();
    assert!((est.mean - 1e10).abs() < 1e5);
}

// ── Array observable tests ────────────────────────────────────────────────

#[test]
fn test_accumulator_with_shape() {
    let acc = Accumulator::with_shape(5, &[3]);
    assert_eq!(acc.shape(), &[3]);
    assert_eq!(acc.bin_capacity(), 5);
}

#[test]
fn test_accumulator_with_shape_2d() {
    let acc = Accumulator::with_shape(10, &[2, 3]);
    assert_eq!(acc.shape(), &[2, 3]);
}

#[test]
fn test_accumulator_add_array_basic() {
    let mut acc = Accumulator::new(2);
    acc.add_array(&[1.0, 2.0, 3.0]);

    assert_eq!(acc.shape(), &[3]);
    assert_eq!(acc.total_count(), 1);
}

#[test]
fn test_accumulator_add_array_completes_bins() {
    let bin_size = 2;
    let mut acc = Accumulator::new(bin_size);

    // Fill two bins (4 samples)
    for _ in 0..2 {
        acc.add_array(&[1.0, 2.0]);
    }
    assert_eq!(acc.num_bins(), 1);

    for _ in 0..2 {
        acc.add_array(&[3.0, 4.0]);
    }
    assert_eq!(acc.num_bins(), 2);
}

#[test]
fn test_accumulator_add_array_finalize() {
    let bin_size = 2;
    let mut acc = Accumulator::new(bin_size);

    // Fill one bin: each sample is [10.0, 20.0]
    acc.add_array(&[10.0, 20.0]);
    acc.add_array(&[10.0, 20.0]);

    let est = acc.finalize();
    // Each bin stores mean per component: [10, 20]
    // Finalize computes mean over components then mean over bins
    // For a single bin: mean of [10.0, 20.0] = 15.0
    assert!((est.mean - 15.0).abs() < 1e-10);
    assert_eq!(est.n_bins, 1);
}

#[test]
fn test_accumulator_add_array_empty_ignored() {
    let mut acc = Accumulator::new(2);
    acc.add_array(&[]);
    assert_eq!(acc.total_count(), 0);
    assert!(acc.shape().is_empty());
}

#[test]
fn test_accumulator_shape_scalar() {
    let acc: Accumulator = Accumulator::new(5);
    assert!(acc.shape().is_empty());
}

#[test]
fn test_accumulator_total_count() {
    let mut acc = Accumulator::new(10);
    acc.add(1.0);
    acc.add(2.0);
    acc.add(3.0);
    assert_eq!(acc.total_count(), 3);
    assert_eq!(acc.total_samples(), 3);
}

#[test]
fn test_accumulator_rebin_means_scalar() {
    let mut acc = Accumulator::new(2);
    // Fill 2 bins: bin0 = mean(1,2)=1.5, bin1 = mean(3,4)=3.5
    acc.add(1.0);
    acc.add(2.0);
    acc.add(3.0);
    acc.add(4.0);

    let means = acc.rebin_means();
    assert_eq!(means.len(), 2);
    assert!((means[0] - 1.5).abs() < 1e-10);
    assert!((means[1] - 3.5).abs() < 1e-10);
}

#[test]
fn test_accumulator_rebin_means_empty() {
    let acc = Accumulator::new(5);
    assert!(acc.rebin_means().is_empty());
}

#[test]
fn test_accumulator_bin_matrix() {
    let mut acc = Accumulator::new(2);
    acc.add_array(&[1.0, 2.0]);
    acc.add_array(&[3.0, 4.0]); // completes bin 0: [2, 3]
    acc.add_array(&[5.0, 6.0]);
    acc.add_array(&[7.0, 8.0]); // completes bin 1: [6, 7]

    let matrix = acc.bin_matrix();
    assert_eq!(matrix.shape(), &[2, 2]);
}

#[test]
fn test_accumulator_bin_matrix_empty() {
    let acc = Accumulator::new(5);
    let matrix = acc.bin_matrix();
    assert_eq!(matrix.shape(), &[0]);
}

#[test]
fn test_accumulator_autocorr_time_few_bins() {
    let acc = Accumulator::new(2);
    assert_eq!(acc.autocorr_time(), 1.0);
}

#[test]
fn test_accumulator_autocorr_time_uncorrelated() {
    let mut acc = Accumulator::new(1);
    // Alternating values → low autocorrelation
    for i in 0..20 {
        acc.add(if i % 2 == 0 { 1.0 } else { -1.0 });
    }
    let tau = acc.autocorr_time();
    assert!(tau >= 0.0);
}

#[test]
fn test_accumulator_covariance_returns_none_for_scalar() {
    let mut acc = Accumulator::new(1);
    acc.add(1.0);
    acc.add(2.0);
    acc.add(3.0);
    acc.add(4.0);
    assert!(acc.covariance().is_none());
}

#[test]
fn test_accumulator_covariance_array() {
    let mut acc = Accumulator::with_shape(1, &[2]);
    // 4 samples with 2 components each
    acc.add_array(&[1.0, 2.0]);
    acc.add_array(&[3.0, 4.0]);
    acc.add_array(&[5.0, 6.0]);
    acc.add_array(&[7.0, 8.0]);

    let cov = acc.covariance();
    assert!(cov.is_some());
    let cov = cov.unwrap();
    assert_eq!(cov.shape(), &[2, 2]);
}

#[test]
fn test_accumulator_covariance_insufficient_bins() {
    let mut acc = Accumulator::with_shape(10, &[2]);
    acc.add_array(&[1.0, 2.0]);
    acc.add_array(&[3.0, 4.0]);
    // Only 0 complete bins
    assert!(acc.covariance().is_none());
}

#[test]
fn test_accumulator_num_bins() {
    let mut acc = Accumulator::new(2);
    assert_eq!(acc.num_bins(), 0);
    acc.add(1.0);
    acc.add(2.0);
    assert_eq!(acc.num_bins(), 1);
    acc.add(3.0);
    acc.add(4.0);
    assert_eq!(acc.num_bins(), 2);
}
