use carlo_rs::{save_json, Estimate, Metadata, Results};
use chrono::Utc;

#[test]
fn test_save_json_writes_file() {
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
    results.add(
        "Magnetization",
        Estimate {
            mean: 0.8,
            stderr: 0.01,
            autocorr_time: 2.0,
            n_bins: 100,
        },
    );

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("results.json");

    save_json(&results, &path).unwrap();

    assert!(path.exists());

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.contains("Energy"));
    assert!(content.contains("Magnetization"));
    assert!(content.contains("1.5"));
    assert!(content.contains("0.8"));
}

#[test]
fn test_save_json_creates_valid_json() {
    let mut results = Results::new();
    results.add("Obs", Estimate::from_bins(&[1.0, 2.0, 3.0]));

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("output.json");

    save_json(&results, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert!(parsed.get("observables").is_some());
    assert!(parsed.get("metadata").is_some());
}

#[test]
fn test_save_json_with_metadata() {
    let mut results = Results::new();
    results.add("E", Estimate::from_bins(&[1.0]));
    results.set_metadata(Metadata {
        version: "test-1.0".to_string(),
        timestamp: Utc::now(),
        base_seed: 999,
        thermalization_sweeps: 50,
        measurement_sweeps: 500,
        n_tasks: 3,
    });

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("meta.json");

    save_json(&results, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.contains("test-1.0"));
    assert!(content.contains("999"));
    assert!(content.contains("500"));
}

#[test]
fn test_save_json_empty_results() {
    let results = Results::new();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("empty.json");

    save_json(&results, &path).unwrap();
    assert!(path.exists());
}

#[test]
fn test_save_json_roundtrip() {
    let mut results = Results::new();
    results.add("Energy", Estimate::from_bins(&[1.0, 2.0, 3.0, 4.0]));
    results.add("Mag", Estimate::from_bins(&[0.5, 0.6, 0.7, 0.8]));

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("roundtrip.json");

    save_json(&results, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();

    let observables = parsed["observables"].as_object().unwrap();
    assert!(observables.contains_key("Energy"));
    assert!(observables.contains_key("Mag"));

    let energy_mean = observables["Energy"]["mean"].as_f64().unwrap();
    assert!((energy_mean - 2.5).abs() < 1e-10);
}
