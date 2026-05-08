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
