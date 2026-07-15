use carlo_rs::{Estimate, Metadata, Results};
use chrono::Utc;

#[test]
fn test_results_creation() {
    let mut results = Results::new();

    results.add(
        "Energy",
        Estimate {
            mean: 1.5,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let est = results.get("Energy").expect("Energy");
    assert!((est.mean - 1.5).abs() < 1e-10);
}

#[test]
fn test_results_json_output() {
    let mut results = Results::new();
    results.add(
        "Magnetization",
        Estimate {
            mean: 0.8,
            stderr: 0.01,
            autocorr_time: 1.0,
            n_bins: 50,
        },
    );

    let json = results.to_json().unwrap();
    assert!(json.contains("Magnetization"));
    assert!(json.contains("0.8"));
}

#[test]
fn test_results_merge_empty() {
    let results: Vec<Results> = vec![];
    let merged = Results::merge(&results);
    assert!(merged.estimates().is_empty());
}

#[test]
fn test_results_merge_single() {
    let mut results = Results::new();
    results.add(
        "Energy",
        Estimate {
            mean: 1.5,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let merged = Results::merge(&[results.clone()]);
    let est = merged.get("Energy").expect("Energy");
    assert!((est.mean - 1.5).abs() < 1e-10);
    assert_eq!(est.n_bins, 100);
}

#[test]
fn test_results_merge_two() {
    let mut r1 = Results::new();
    r1.add(
        "Energy",
        Estimate {
            mean: 1.5,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let mut r2 = Results::new();
    r2.add(
        "Energy",
        Estimate {
            mean: 1.6,
            stderr: 0.03,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let merged = Results::merge(&[r1, r2]);
    let est = merged.get("Energy").expect("Energy");

    // Weighted mean: (1.5*100 + 1.6*100) / 200 = 1.55
    assert!((est.mean - 1.55).abs() < 1e-10);
    assert_eq!(est.n_bins, 200);
}

#[test]
fn test_results_merge_multiple_observables() {
    let mut r1 = Results::new();
    r1.add(
        "Energy",
        Estimate {
            mean: 1.5,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );
    r1.add(
        "Magnetization",
        Estimate {
            mean: 0.5,
            stderr: 0.01,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let mut r2 = Results::new();
    r2.add(
        "Energy",
        Estimate {
            mean: 1.6,
            stderr: 0.03,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );
    r2.add(
        "Magnetization",
        Estimate {
            mean: 0.6,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );

    let merged = Results::merge(&[r1, r2]);

    assert!(merged.get("Energy").is_some());
    assert!(merged.get("Magnetization").is_some());

    let energy = merged.get("Energy").expect("Energy");
    assert!((energy.mean - 1.55).abs() < 1e-10);

    let mag = merged.get("Magnetization").expect("Magnetization");
    assert!((mag.mean - 0.55).abs() < 1e-10);
}

#[test]
fn test_results_metadata() {
    let mut results = Results::new();
    results.set_metadata(Metadata {
        version: "0.1.0".to_string(),
        timestamp: Utc::now(),
        base_seed: 42,
        thermalization_sweeps: 100,
        measurement_sweeps: 1000,
        n_tasks: 4,
    });

    let meta = results.metadata();
    assert_eq!(meta.base_seed, 42);
    assert_eq!(meta.measurement_sweeps, 1000);
}
