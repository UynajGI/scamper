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
    // Mean of 1..=15 is exactly 8.0
    assert!((est.mean - 8.0).abs() < 1e-10);
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
    // binsize=3: bin0=[-3,-2,-1] mean=-2.0, partial bin1=[0] mean=0.0
    // finalize mean = (-2.0 + 0.0) / 2 = -1.0
    assert!((est.mean - (-1.0)).abs() < 1e-10);
}

#[test]
fn test_accumulator_large_values() {
    let bin_length = 2;
    let mut acc = Accumulator::new(bin_length);

    acc.add(1e10);
    acc.add(1e10);

    let est = acc.finalize();
    assert!((est.mean - 1e10).abs() < 1.0);
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
    acc.add_array(&[1.0, 2.0, 3.0]).unwrap();

    assert_eq!(acc.shape(), &[3]);
    assert_eq!(acc.total_count(), 1);
}

#[test]
fn test_accumulator_add_array_completes_bins() {
    let bin_size = 2;
    let mut acc = Accumulator::new(bin_size);

    // Fill two bins (4 samples)
    for _ in 0..2 {
        acc.add_array(&[1.0, 2.0]).unwrap();
    }
    assert_eq!(acc.num_bins(), 1);

    for _ in 0..2 {
        acc.add_array(&[3.0, 4.0]).unwrap();
    }
    assert_eq!(acc.num_bins(), 2);
}

#[test]
fn test_accumulator_add_array_finalize() {
    let bin_size = 2;
    let mut acc = Accumulator::new(bin_size);

    // Fill one bin: each sample is [10.0, 20.0]
    acc.add_array(&[10.0, 20.0]).unwrap();
    acc.add_array(&[10.0, 20.0]).unwrap();

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
    acc.add_array(&[]).unwrap();
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
    acc.add_array(&[1.0, 2.0]).unwrap();
    acc.add_array(&[3.0, 4.0]).unwrap(); // completes bin 0: [2, 3]
    acc.add_array(&[5.0, 6.0]).unwrap();
    acc.add_array(&[7.0, 8.0]).unwrap(); // completes bin 1: [6, 7]

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
    assert!(tau.is_finite(), "autocorr_time should be finite, got {tau}");
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
    acc.add_array(&[1.0, 2.0]).unwrap();
    acc.add_array(&[3.0, 4.0]).unwrap();
    acc.add_array(&[5.0, 6.0]).unwrap();
    acc.add_array(&[7.0, 8.0]).unwrap();

    let cov = acc.covariance();
    assert!(cov.is_some());
    let cov = cov.unwrap();
    assert_eq!(cov.shape(), &[2, 2]);
}

#[test]
fn test_accumulator_covariance_insufficient_bins() {
    let mut acc = Accumulator::with_shape(10, &[2]);
    acc.add_array(&[1.0, 2.0]).unwrap();
    acc.add_array(&[3.0, 4.0]).unwrap();
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

#[test]
fn test_add_array_shape_mismatch_is_explicit_error() {
    use carlo_rs::CarloError;

    let mut acc = Accumulator::new(10);
    acc.add_array(&[1.0, 2.0]).unwrap();

    // A shorter sample must be rejected, not silently zipped/truncated —
    // in every build profile (this check is not a debug_assert).
    let err = acc.add_array(&[3.0]).unwrap_err();
    match &err {
        CarloError::ShapeMismatch {
            expected_len,
            got_len,
            ..
        } => {
            assert_eq!(*expected_len, 2);
            assert_eq!(*got_len, 1);
        }
        other => panic!("expected ShapeMismatch, got {other:?}"),
    }

    // The rejected sample leaves no trace in the accumulator.
    assert_eq!(acc.total_count(), 1);
    assert_eq!(acc.shape(), &[2]);

    // Empty samples remain no-ops.
    assert!(acc.add_array(&[]).is_ok());
    assert_eq!(acc.total_count(), 1);
}

#[test]
fn test_finalize_array_none_for_scalar() {
    let mut acc = Accumulator::new(2);
    acc.add(1.0);
    acc.add(2.0);
    acc.add(3.0);
    assert!(acc.finalize_array().is_none());
}

#[test]
fn test_finalize_array_matches_independent_scalar_accumulators() {
    // Element-preserving finalize must agree exactly with running each
    // element through its own scalar accumulator fed the same samples.
    let mut acc = Accumulator::new(2);
    let samples: Vec<[f64; 2]> = vec![
        [1.0, 10.0],
        [3.0, 12.0],
        [5.0, 20.0],
        [7.0, 16.0],
        [9.0, 30.0],
    ];

    let mut scalar_a = Accumulator::new(2);
    let mut scalar_b = Accumulator::new(2);
    for [a, b] in &samples {
        acc.add_array(&[*a, *b]).unwrap();
        scalar_a.add(*a);
        scalar_b.add(*b);
    }

    let array = acc.finalize_array().expect("array accumulator");
    assert_eq!(array.shape(), &[2]);
    let est_a = array.element(0).unwrap();
    let est_b = array.element(1).unwrap();

    let ref_a = scalar_a.finalize();
    let ref_b = scalar_b.finalize();
    assert_eq!(est_a.n_bins, ref_a.n_bins);
    assert_eq!(est_b.n_bins, ref_b.n_bins);
    assert!((est_a.mean - ref_a.mean).abs() < 1e-12);
    assert!((est_a.stderr - ref_a.stderr).abs() < 1e-12);
    assert!((est_a.autocorr_time - ref_a.autocorr_time).abs() < 1e-12);
    assert!((est_b.mean - ref_b.mean).abs() < 1e-12);
    assert!((est_b.stderr - ref_b.stderr).abs() < 1e-12);
    assert!((est_b.autocorr_time - ref_b.autocorr_time).abs() < 1e-12);
}

#[test]
fn test_finalize_array_partial_bin_included_per_element() {
    // 5 samples with bin capacity 2 → 2 complete bins + 1 partial sample.
    // The partial bin must join each element's column as a normalized mean.
    let mut acc = Accumulator::new(2);
    acc.add_array(&[4.0, 8.0]).unwrap();
    acc.add_array(&[4.0, 8.0]).unwrap(); // bin 0 mean: [4, 8]
    acc.add_array(&[2.0, 6.0]).unwrap();
    acc.add_array(&[2.0, 6.0]).unwrap(); // bin 1 mean: [2, 6]
    acc.add_array(&[9.0, 3.0]).unwrap(); // partial: [9, 3]

    let array = acc.finalize_array().expect("array accumulator");
    let est = array.element(0).unwrap();
    // Element 0 column: [4, 2, 9] → mean 5, n_bins 3.
    assert_eq!(est.n_bins, 3);
    assert!((est.mean - 5.0).abs() < 1e-12);
    let est = array.element(1).unwrap();
    assert!((est.mean - 17.0 / 3.0).abs() < 1e-12);
}

#[test]
fn test_finalize_array_tau_ignores_partial_bin_like_scalar_path() {
    // Review follow-up (2026-10-22): the scalar path estimates tau from
    // completed bins only. Including the normalized partial bin in the
    // per-element tau column let one outlier sample inflate var_rebin and
    // spuriously raise tau (reviewer probe: 1.104 vs 1.0), so finalize_array
    // now matches the scalar convention exactly. The partial bin still
    // contributes to mean and stderr.
    let mut acc = Accumulator::new(2);
    let mut scalar = Accumulator::new(2);
    // Completed-bin means alternate 0/1; the partial-bin sample is an
    // outlier (500), which would dominate any variance that includes it.
    for round in 0..4_u32 {
        let value = f64::from(round % 2);
        for _ in 0..2 {
            acc.add_array(&[value, -value]).unwrap();
            scalar.add(value);
        }
    }
    acc.add_array(&[500.0, -500.0]).unwrap();
    scalar.add(500.0);

    let array = acc.finalize_array().expect("array accumulator");
    let reference = scalar.finalize();
    let element = array.element(0).unwrap();

    // Exact parity with the scalar statistics on identical samples.
    assert_eq!(element.autocorr_time, reference.autocorr_time);
    assert!((element.mean - reference.mean).abs() < 1e-9);
    assert!((element.stderr - reference.stderr).abs() < 1e-9);
    assert_eq!(element.n_bins, reference.n_bins);

    // The outlier partial sample is in the mean (not silently dropped)...
    assert!(element.mean > 50.0);
    // ...but does not leak into tau: the alternating completed bins leave
    // the raw estimator at its floor, stored as 1.0 exactly as the scalar
    // path reports it.
    assert_eq!(element.autocorr_time, 1.0);
}
