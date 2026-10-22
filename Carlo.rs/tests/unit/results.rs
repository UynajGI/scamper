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
        attributes: std::collections::BTreeMap::new(),
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

// ── Array observables (C0: element-preserving vector estimates) ────────────

fn sample_array_observations() -> Vec<Vec<f64>> {
    // Four samples of a two-element observable, bin capacity 2.
    vec![
        vec![1.0, 10.0],
        vec![3.0, 12.0],
        vec![5.0, 20.0],
        vec![7.0, 16.0],
    ]
}

fn accumulated_array(name: &str) -> carlo_rs::ArrayEstimate {
    use carlo_rs::Measurements;
    let mut meas = Measurements::new(2);
    for sample in sample_array_observations() {
        meas.add_sample_array(name, &sample).unwrap();
    }
    meas.finalize_arrays().remove(name).expect("array estimate")
}

#[test]
fn test_results_array_accessors() {
    let arr = accumulated_array("Curve");
    let mut results = Results::new();
    results.add_array("Curve", arr.clone());

    assert!(results.get_array("Curve").is_some());
    assert_eq!(results.array_estimates().len(), 1);
    assert_eq!(results.get_array("Curve").unwrap().shape(), &[2]);
    assert!(results.get("Curve").is_none());
}

#[test]
fn test_results_array_json_round_trip() {
    let arr = accumulated_array("Curve");
    let mut results = Results::new();
    results.add(
        "Energy",
        Estimate {
            mean: 1.0,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 4,
        },
    );
    results.add_array("Curve", arr);

    let json = results.to_json().unwrap();
    assert!(json.contains("\"array_observables\""));
    assert!(json.contains("\"Curve\""));

    let back = Results::from_json(&json).unwrap();
    let curve = back.get_array("Curve").expect("Curve survives round trip");
    assert_eq!(curve.shape(), &[2]);
    for (lhs, rhs) in curve
        .elements()
        .iter()
        .zip(accumulated_array("Curve").elements())
    {
        assert!((lhs.mean - rhs.mean).abs() < 1e-15);
        assert!((lhs.stderr - rhs.stderr).abs() < 1e-15);
        assert_eq!(lhs.n_bins, rhs.n_bins);
        assert!((lhs.autocorr_time - rhs.autocorr_time).abs() < 1e-15);
    }
    assert!((back.get("Energy").unwrap().mean - 1.0).abs() < 1e-15);
}

#[test]
fn test_results_scalar_json_has_no_array_key() {
    // Legacy promise: scalar-only output is byte-compatible — the
    // array_observables key must not appear at all.
    let mut results = Results::new();
    results.add(
        "Energy",
        Estimate {
            mean: 1.0,
            stderr: 0.1,
            autocorr_time: 1.0,
            n_bins: 4,
        },
    );
    let json = results.to_json().unwrap();
    assert!(!json.contains("array_observables"));
    assert!(!json.contains("attributes"));
}

#[test]
fn test_results_from_json_legacy_fixture() {
    // A results file written before array observables and provenance
    // attributes existed must still load.
    let legacy = r#"{
  "observables": {
    "Energy": { "mean": -1.234, "stderr": 0.005, "autocorr_time": 1.5, "n_bins": 512 }
  },
  "metadata": {
    "version": "0.1.0",
    "timestamp": "2026-01-01T00:00:00Z",
    "base_seed": 42,
    "thermalization_sweeps": 100,
    "measurement_sweeps": 1000,
    "n_tasks": 1
  }
}"#;

    let results = Results::from_json(legacy).unwrap();
    let energy = results.get("Energy").expect("Energy");
    assert!((energy.mean - (-1.234)).abs() < 1e-15);
    assert_eq!(energy.n_bins, 512);
    assert!(results.array_estimates().is_empty());
    assert!(results.metadata().attributes.is_empty());
}

#[test]
fn test_merge_arrays_matches_scalar_merge_elementwise() {
    // Element-wise array merging must reproduce exactly the scalar merge
    // applied to the same per-element estimates.
    let build = |offset: f64| -> Results {
        let arr = carlo_rs::ArrayEstimate::from_element_bins(
            &[2],
            &[
                vec![1.0 + offset, 2.0 + offset, 3.0 + offset],
                vec![10.0, 11.0, 12.0],
            ],
        )
        .unwrap();
        let mut r = Results::new();
        r.add_array("Curve", arr);

        let e0 = Estimate::from_bins(&[1.0 + offset, 2.0 + offset, 3.0 + offset]);
        let e1 = Estimate::from_bins(&[10.0, 11.0, 12.0]);
        r.add("e0", e0);
        r.add("e1", e1);
        r
    };

    let merged = Results::merge(&[build(0.0), build(1.0)]);

    let curve = merged.get_array("Curve").expect("Curve merged");
    assert_eq!(curve.shape(), &[2]);
    for (flat, scalar_name) in [(0usize, "e0"), (1, "e1")] {
        let element = curve.element(flat).unwrap();
        let scalar = merged.get(scalar_name).expect(scalar_name);
        assert!(
            (element.mean - scalar.mean).abs() < 1e-12,
            "{scalar_name} mean"
        );
        assert!(
            (element.stderr - scalar.stderr).abs() < 1e-12,
            "{scalar_name} stderr"
        );
        assert_eq!(element.n_bins, scalar.n_bins);
        assert!((element.autocorr_time - scalar.autocorr_time).abs() < 1e-12);
    }
}

#[test]
fn test_merge_skips_shape_mismatched_arrays() {
    let mut r1 = Results::new();
    r1.add_array(
        "Curve",
        carlo_rs::ArrayEstimate::from_element_bins(&[2], &[vec![1.0, 2.0], vec![3.0, 4.0]])
            .unwrap(),
    );
    let mut r2 = Results::new();
    r2.add_array(
        "Curve",
        carlo_rs::ArrayEstimate::from_element_bins(
            &[3],
            &[vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]],
        )
        .unwrap(),
    );

    let merged = Results::merge(&[r1, r2]);
    assert!(
        merged.get_array("Curve").is_none(),
        "mismatched shapes must not be flattened"
    );
    assert!(merged.array_estimates().is_empty());
}

#[test]
fn test_merge_partial_array_observables() {
    // An array observable missing from one result merges the ones present.
    let mut r1 = Results::new();
    r1.add_array(
        "Curve",
        carlo_rs::ArrayEstimate::from_element_bins(&[2], &[vec![1.0, 2.0], vec![3.0, 4.0]])
            .unwrap(),
    );
    let r2 = Results::new();

    let merged = Results::merge(&[r1, r2]);
    assert!(merged.get_array("Curve").is_some());
    assert_eq!(merged.get_array("Curve").unwrap().n_elements(), 2);
}

#[test]
fn test_metadata_attributes_round_trip() {
    let mut metadata = Metadata::default();
    metadata.set_attribute("model_family", "directed_percolation");
    metadata.set_attribute("lattice_width", 64);

    let mut results = Results::new();
    results.set_metadata(metadata);

    let json = results.to_json().unwrap();
    assert!(json.contains("\"attributes\""));
    assert!(json.contains("directed_percolation"));

    let back = Results::from_json(&json).unwrap();
    let attributes = &back.metadata().attributes;
    assert_eq!(
        attributes.get("model_family").and_then(|v| v.as_str()),
        Some("directed_percolation")
    );
    assert_eq!(
        attributes.get("lattice_width").and_then(|v| v.as_u64()),
        Some(64)
    );
}

#[test]
fn test_metadata_attributes_merge_semantics() {
    let mut agree_a = Metadata::default();
    agree_a.set_attribute("policy", "degree_desc");
    let mut agree_b = Metadata::default();
    agree_b.set_attribute("policy", "degree_desc");
    let mut conflict = Metadata::default();
    conflict.set_attribute("policy", "random");

    let mut ra = Results::new();
    ra.set_metadata(agree_a);
    let mut rb = Results::new();
    rb.set_metadata(agree_b);
    let merged = Results::merge(&[ra, rb]);
    assert_eq!(
        merged
            .metadata()
            .attribute("policy")
            .and_then(|v| v.as_str()),
        Some("degree_desc")
    );

    let mut rc = Results::new();
    rc.set_metadata(conflict);
    let merged = Results::merge(&[merged, rc]);
    assert!(merged.metadata().attribute("policy").is_none());
}

#[test]
fn test_merge_stderr_is_pooled_bin_variance() {
    // Review follow-up (2026-10-22): merged stderr must equal the standard
    // error of the pooled-bin weighted mean, Var = sum_i n_i s_i^2 / N^2
    // with s_i^2 = stderr_i^2 * n_i. The pre-review formula divided by
    // n_runs^2 and underestimated by sqrt(n_runs) for equal runs.
    let run = |mean: f64, stderr: f64| {
        let mut r = Results::new();
        r.add(
            "Energy",
            Estimate {
                mean,
                stderr,
                autocorr_time: 1.0,
                n_bins: 100,
            },
        );
        r
    };

    // Two equal runs: s_1^2 = 0.02^2*100, s_2^2 = 0.03^2*100.
    // Var = (100*0.04 + 100*0.09) / 200^2 = 13/40000 -> stderr = sqrt(13)/200.
    let merged = Results::merge(&[run(1.5, 0.02), run(1.6, 0.03)]);
    let expected = (13.0_f64 / 40000.0).sqrt();
    let est = merged.get("Energy").expect("Energy");
    assert!(
        (est.stderr - expected).abs() < 1e-12,
        "merged stderr {} != pooled-bin {}",
        est.stderr,
        expected
    );
    // The old formula gave sqrt(0.13)/20 ~ 0.018028? No: it gave
    // sqrt(13/100)/200*sqrt(200)/sqrt(200)... its value was
    // sqrt(0.0013/4)/sqrt(2) — in any case strictly smaller; assert we
    // are not below the correct value (regression direction guard).
    assert!(est.stderr > 0.018 && est.stderr < 0.019);

    // Single run: merged stderr is the run's own stderr, unchanged.
    let merged = Results::merge(&[run(1.5, 0.02)]);
    assert!((merged.get("Energy").unwrap().stderr - 0.02).abs() < 1e-12);

    // Unequal bin counts weight correctly:
    // Var = (400*0.02^2*400... = sum n_i * (stderr_i^2 n_i) / (500^2).
    let mut a = Results::new();
    a.add(
        "Energy",
        Estimate {
            mean: 1.0,
            stderr: 0.02,
            autocorr_time: 1.0,
            n_bins: 400,
        },
    );
    let mut b = Results::new();
    b.add(
        "Energy",
        Estimate {
            mean: 2.0,
            stderr: 0.05,
            autocorr_time: 1.0,
            n_bins: 100,
        },
    );
    let merged = Results::merge(&[a, b]);
    let expected = ((400.0 * (0.02_f64 * 20.0).powi(2)) + (100.0 * (0.05_f64 * 10.0).powi(2)))
        / (500.0_f64).powi(2);
    let expected = expected.sqrt();
    assert!((merged.get("Energy").unwrap().stderr - expected).abs() < 1e-12);
}
