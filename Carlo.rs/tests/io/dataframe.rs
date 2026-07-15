use carlo_rs::dataframe;
use serde_json::json;

fn write_results_file(
    dir: &std::path::Path,
    filename: &str,
    data: serde_json::Value,
) -> std::path::PathBuf {
    let path = dir.join(filename);
    std::fs::write(&path, serde_json::to_string_pretty(&data).unwrap()).unwrap();
    path
}

#[test]
fn test_dataframe_reads_single_task() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 8, "beta": 0.5},
            "results": {
                "Energy": {
                    "mean": -1.23,
                    "error": 0.01,
                    "rebin_len": 100,
                    "autocorr_time": 2.5
                }
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].task, "task0001");
    assert_eq!(rows[0].observable, "Energy");
    assert_eq!(rows[0].mean, json!(-1.23));
    assert_eq!(rows[0].error, json!(0.01));
    assert_eq!(rows[0].rebin_len, 100);
    assert!(rows[0].covariance.is_none());
}

#[test]
fn test_dataframe_reads_multiple_observables() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 4},
            "results": {
                "Energy": {"mean": -1.0, "error": 0.1, "rebin_len": 50, "autocorr_time": 1.0},
                "Mag": {"mean": 0.5, "error": 0.05, "rebin_len": 50, "autocorr_time": 2.0}
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 2);

    let obs_names: Vec<&str> = rows.iter().map(|r| r.observable.as_str()).collect();
    assert!(obs_names.contains(&"Energy"));
    assert!(obs_names.contains(&"Mag"));
}

#[test]
fn test_dataframe_reads_multiple_tasks() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 4},
            "results": {
                "E": {"mean": -1.0, "error": 0.1, "rebin_len": 10, "autocorr_time": 1.0}
            }
        },
        {
            "task": "task0002",
            "parameters": {"L": 8},
            "results": {
                "E": {"mean": -2.0, "error": 0.2, "rebin_len": 10, "autocorr_time": 1.5}
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].task, "task0001");
    assert_eq!(rows[1].task, "task0002");
}

#[test]
fn test_dataframe_skips_null_observables() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {},
            "results": {
                "ValidObs": {"mean": 1.0, "error": 0.1, "rebin_len": 10, "autocorr_time": 1.0},
                "NullObs": null
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].observable, "ValidObs");
}

#[test]
fn test_dataframe_skips_missing_results_key() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 4}
        },
        {
            "task": "task0002",
            "parameters": {"L": 8},
            "results": {
                "E": {"mean": -1.0, "error": 0.1, "rebin_len": 10, "autocorr_time": 1.0}
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].task, "task0002");
}

#[test]
fn test_dataframe_empty_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = write_results_file(dir.path(), "empty.json", json!([]));

    let rows = dataframe(&path).unwrap();
    assert!(rows.is_empty());
}

#[test]
fn test_dataframe_nonexistent_file_errors() {
    let path = std::path::PathBuf::from("/nonexistent/results.json");
    let result = dataframe(&path);
    assert!(result.is_err());
}

#[test]
fn test_dataframe_malformed_json_errors() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("bad.json");
    std::fs::write(&path, "not valid json {{{").unwrap();

    let result = dataframe(&path);
    assert!(result.is_err());
}

#[test]
fn test_dataframe_preserves_parameters() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 16, "beta": 0.44, "J": 1.0},
            "results": {
                "E": {"mean": -1.0, "error": 0.1, "rebin_len": 10, "autocorr_time": 1.0}
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows[0].parameters.len(), 3);
    assert_eq!(rows[0].parameters["L"], json!(16));
    assert_eq!(rows[0].parameters["beta"], json!(0.44));
    assert_eq!(rows[0].parameters["J"], json!(1.0));
}

#[test]
fn test_dataframe_array_observable() {
    let dir = tempfile::tempdir().unwrap();
    let data = json!([
        {
            "task": "task0001",
            "parameters": {"L": 4},
            "results": {
                "Corr": {
                    "mean": [1.0, 0.5, 0.25],
                    "error": [0.1, 0.05, 0.025],
                    "rebin_len": 20,
                    "autocorr_time": [1.0, 1.5, 2.0]
                }
            }
        }
    ]);
    let path = write_results_file(dir.path(), "results.json", data);

    let rows = dataframe(&path).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(rows[0].mean.is_array());
    assert!(rows[0].error.is_array());
}
