//! HDF5 result file merging test.
//!
//! Creates measurement files via Accumulator (raw bins), merges them,
//! and verifies the merged estimate matches manual calculation.

#![cfg(feature = "hdf5")]

use carlo_rs::merge::{merge_results_from_files, MergeOptions};
use carlo_rs::Accumulator;
use std::path::PathBuf;

fn make_temp_dir() -> PathBuf {
    let dir = std::env::temp_dir().join("carlo_merge_test");
    let _ = std::fs::create_dir_all(&dir);
    dir
}

/// Write a measurement file with one scalar observable and known bin values.
fn write_meas_file(path: &std::path::Path, obs_name: &str, values: &[f64], binsize: usize) {
    use hdf5::File as Hdf5File;
    let mut acc = Accumulator::new(binsize);
    for &v in values {
        acc.add(v);
    }
    let file = Hdf5File::create(path).unwrap();
    let mut obs_group = file.create_group("observables").unwrap();
    acc.write_hdf5(&mut obs_group, obs_name).unwrap();
}

#[test]
fn hdf5_merge_two_files_produces_combined_mean() {
    let dir = make_temp_dir();

    // File 1: Energy values centered around 1.0
    let vals1: Vec<f64> = (0..100).map(|i| 1.0 + 0.01 * (i as f64 - 50.0)).collect();
    let path1 = dir.join("task_001.meas.h5");
    write_meas_file(&path1, "Energy", &vals1, 10);

    // File 2: Energy values centered around 3.0
    let vals2: Vec<f64> = (0..100).map(|i| 3.0 + 0.01 * (i as f64 - 50.0)).collect();
    let path2 = dir.join("task_002.meas.h5");
    write_meas_file(&path2, "Energy", &vals2, 10);

    let files = vec![path1, path2];
    let opts = MergeOptions::default();
    let merged = merge_results_from_files(&files, &opts).expect("merge should succeed");

    let energy = merged.get("Energy").expect("Energy in merged results");
    let mean_val: f64 = energy.mean.iter().next().copied().unwrap_or(f64::NAN);

    // Combined mean should be ~2.0 (average of 1.0 and 3.0)
    assert!(
        (mean_val - 2.0).abs() < 0.05,
        "merged mean should be ~2.0, got {mean_val}"
    );
}

#[test]
fn hdf5_merge_single_file_returns_observables() {
    let dir = make_temp_dir();

    let vals: Vec<f64> = (0..50).map(|i| 0.5 + 0.001 * (i as f64)).collect();
    let path = dir.join("task_single.meas.h5");
    write_meas_file(&path, "Magnetization", &vals, 5);

    let files = vec![path];
    let opts = MergeOptions::default();
    let merged = merge_results_from_files(&files, &opts).expect("merge single file");

    let mag = merged
        .get("Magnetization")
        .expect("Magnetization should survive merge");
    let mean_val: f64 = mag.mean.iter().next().copied().unwrap_or(f64::NAN);
    assert!(
        (mean_val - 0.527).abs() < 0.05,
        "merged mean should be ~0.527, got {mean_val}"
    );
}

#[test]
fn hdf5_merge_empty_file_list_returns_empty() {
    let opts = MergeOptions::default();
    let merged = merge_results_from_files(&[], &opts).expect("empty merge");
    assert!(merged.is_empty());
}
