use carlo_rs::{ComplexEstimate, Estimate, Metadata, Results};
use chrono::Utc;
use std::collections::HashMap;

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

// ── Complex results ───────────────────────────────────────────────────────

#[test]
fn test_results_add_and_get_complex() {
    let mut results = Results::new();
    results.add_complex(
        "Susceptibility",
        ComplexEstimate::new(
            Estimate::from_bins(&[1.0, 2.0, 3.0]),
            Estimate::from_bins(&[4.0, 5.0, 6.0]),
        ),
    );

    let cr = results.get_complex("Susceptibility").expect("complex obs");
    assert!((cr.re.mean - 2.0).abs() < 1e-10);
    assert!((cr.im.mean - 5.0).abs() < 1e-10);
}

#[test]
fn test_results_from_measurements_with_complex() {
    let mut scalars = HashMap::new();
    scalars.insert(
        "Energy".into(),
        Estimate {
            mean: 1.0,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 5,
        },
    );

    let mut complex = HashMap::new();
    complex.insert(
        "Green".into(),
        ComplexEstimate::new(
            Estimate::from_bins(&[0.1, 0.2]),
            Estimate::from_bins(&[0.3, 0.4]),
        ),
    );

    let results = Results::from_measurements_with_complex(&scalars, &complex);

    assert!(results.get("Energy").is_some());
    assert!(results.get_complex("Green").is_some());

    let complex_map = results.complex_estimates();
    assert_eq!(complex_map.len(), 1);
    assert!(complex_map.contains_key("Green"));
}

#[test]
fn test_results_complex_json_no_complex_when_empty() {
    let mut results = Results::new();
    results.add("Energy", Estimate::from_bins(&[1.0, 2.0]));

    let json = results.to_json().unwrap();
    assert!(!json.contains("complex_observables"));
}

#[test]
fn test_results_merge_mismatched_observables() {
    let mut r1 = Results::new();
    r1.add(
        "Energy",
        Estimate {
            mean: 1.0,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 10,
        },
    );
    r1.add(
        "Mag",
        Estimate {
            mean: 0.5,
            stderr: 0.05,
            autocorr_time: 1.0,
            n_bins: 10,
        },
    );

    let mut r2 = Results::new();
    r2.add(
        "Energy",
        Estimate {
            mean: 2.0,
            stderr: 0.2,
            autocorr_time: 1.0,
            n_bins: 10,
        },
    );
    // r2 does NOT have "Mag"

    let merged = Results::merge(&[r1, r2]);
    assert!(merged.get("Energy").is_some());
    assert!(merged.get("Mag").is_some());

    let energy = merged.get("Energy").unwrap();
    assert!((energy.mean - 1.5).abs() < 1e-10);

    let mag = merged.get("Mag").unwrap();
    assert!((mag.mean - 0.5).abs() < 1e-10);
}

#[test]
fn test_results_default() {
    let results = Results::default();
    assert!(results.estimates().is_empty());
    assert!(results.complex_estimates().is_empty());
}

#[test]
fn test_results_estimates_accessor() {
    let mut results = Results::new();
    results.add("A", Estimate::from_bins(&[1.0]));
    results.add("B", Estimate::from_bins(&[2.0]));

    let ests = results.estimates();
    assert_eq!(ests.len(), 2);
    assert!(ests.contains_key("A"));
    assert!(ests.contains_key("B"));
}

#[test]
fn test_metadata_default() {
    let meta = Metadata::default();
    assert_eq!(meta.base_seed, 0);
    assert_eq!(meta.n_tasks, 1);
    assert!(!meta.version.is_empty());
}
