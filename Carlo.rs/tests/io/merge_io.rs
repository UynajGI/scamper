// File: Carlo.rs/tests/merge_io_test.rs
use carlo_rs::merge::{calc_rebin_count, calc_rebin_length, list_meas_files, MergeOptions};
use std::path::PathBuf;

#[test]
fn test_merge_options_default() {
    let opts = MergeOptions::default();
    assert!(opts.rebin_length.is_none());
    assert_eq!(opts.sample_skip, 0);
    assert!(!opts.estimate_covariance);
}

#[test]
fn test_list_meas_files_empty() {
    let dir = PathBuf::from("/nonexistent/path");
    let result = list_meas_files(&dir);
    assert!(result.is_ok());
    assert!(result.unwrap().is_empty());
}

#[test]
fn test_calc_rebin_length() {
    // With explicit rebin_length
    assert_eq!(calc_rebin_length(1000, Some(100)), 100);

    // Auto calculation
    let auto = calc_rebin_length(1000, None);
    assert!(auto > 0);
    assert!(auto < 1000);

    // Edge case: zero samples
    assert_eq!(calc_rebin_length(0, None), 1);
}

#[test]
fn test_calc_rebin_count() {
    // Small sample count: no rebinning
    assert_eq!(calc_rebin_count(5, 10), 5);

    // Large sample count: rebinning applied
    let count = calc_rebin_count(1000, 10);
    assert!(count > 10);
    assert!(count < 1000);
}
