//! ResultTools — load and analyze Carlo results.json files.
//!
//! Mirrors Carlo.jl's `resulttools/ResultTools.jl`. Provides utilities
//! for reading merged simulation results back into Rust for analysis.
//!
//! # Key Functions
//!
//! - [`dataframe()`]: Load a `*.results.json` file into a flat row structure
//! - [`measurement_from_obs()`]: Build a typed measurement from a JSON observable dict
//! - [`make_scalar()`]: Flatten single-element JSON arrays to scalars
//! - [`recursive_stack()`]: Reconstruct nested arrays from JSON nesting

use serde_json::Value as Json;
use std::path::Path;

/// A flat row of task results, suitable for DataFrame conversion.
///
/// Each row contains the task name, all parameters, and one observable's
/// value ± error (plus optional covariance). Rows are produced by [`dataframe()`].
#[derive(Debug, Clone)]
pub struct ResultRow {
    /// Task name (basename of task directory).
    pub task: String,
    /// Task parameters (all key/value pairs).
    pub parameters: serde_json::Map<String, Json>,
    /// Observable name.
    pub observable: String,
    /// Mean value (scalar or array).
    pub mean: Json,
    /// Standard error (scalar or array, matching mean shape).
    pub error: Json,
    /// Rebin length used in analysis.
    pub rebin_len: u64,
    /// Autocorrelation time estimate (scalar or array).
    pub autocorr_time: Json,
    /// Covariance tensor (optional, present only for array observables
    /// when covariance estimation was enabled).
    pub covariance: Option<Json>,
}

/// Convert a single-element JSON array to a scalar.
///
/// Mirrors Carlo.jl's `make_scalar`: if the value is an array of
/// length 1, returns the single element; otherwise returns the value as-is.
///
/// # Example
///
/// ```rust,ignore
/// let val = serde_json::json!([1.23]);
/// assert_eq!(make_scalar(&val), &serde_json::json!(1.23));
///
/// let val = serde_json::json!([1.0, 2.0]);
/// assert_eq!(make_scalar(&val), &val); // no flatten
/// ```
pub fn make_scalar(value: &Json) -> &Json {
    // We can't return a reference to a modified value, so this function
    // takes a reference and returns the scalar if possible.
    // For actual mutation, use `make_scalar_owned`.
    match value {
        Json::Array(arr) if arr.len() == 1 => &arr[0],
        _ => value,
    }
}

/// Owned version of [`make_scalar`].
pub fn make_scalar_owned(value: Json) -> Json {
    match value {
        Json::Array(ref arr) if arr.len() == 1 => arr[0].clone(),
        other => other,
    }
}

/// Recursively stack nested JSON arrays into a single ndarray-like structure.
///
/// When JSON serializes multi-dimensional arrays, they become nested
/// `[[...], [...]]`. This function attempts to recursively combine them
/// back into a flat structure, mirroring Carlo.jl's `recursive_stack`.
///
/// For JSON objects with `"re"` and `"im"` keys, returns the original
/// object (complex number format).
pub fn recursive_stack(value: Json) -> Json {
    match value {
        Json::Object(_) => value,
        Json::Array(ref arr) => {
            if arr.is_empty() || arr.iter().any(|v| v.is_null()) {
                return value;
            }
            // Check if all elements are arrays of same length → try to stack
            if arr.iter().all(|v| v.is_array()) {
                let inner_lens: Vec<usize> = arr
                    .iter()
                    .map(|v| v.as_array().map(|a| a.len()).unwrap_or(0))
                    .collect();
                let first_len = inner_lens[0];
                if first_len > 0 && inner_lens.iter().all(|&l| l == first_len) {
                    let stacked: Vec<Json> = arr
                        .iter()
                        .flat_map(|v| v.as_array().unwrap().clone())
                        .map(recursive_stack)
                        .collect();
                    Json::Array(stacked)
                } else {
                    Json::Array(arr.iter().cloned().map(recursive_stack).collect())
                }
            } else {
                Json::Array(arr.iter().cloned().map(recursive_stack).collect())
            }
        }
        _ => value,
    }
}

/// Build a typed measurement from a JSON observable dictionary.
///
/// Takes a JSON object with keys `mean`, `error`, `rebin_len`,
/// `autocorr_time`, and optionally `covariance`, and returns a
/// structured measurement with properly stacked arrays.
///
/// Mirrors Carlo.jl's `measurement_from_obs`.
///
/// # Returns
///
/// A tuple of (mean, error, rebin_len, autocorr_time, covariance)
/// where the array values have been through [`recursive_stack`].
pub fn measurement_from_obs(obs: &Json) -> Option<(Json, Json, u64, Json, Option<Json>)> {
    let obj = obs.as_object()?;

    let mean = recursive_stack(obj.get("mean")?.clone());
    let error = recursive_stack(obj.get("error")?.clone());
    let rebin_len = obj.get("rebin_len").and_then(|v| v.as_u64()).unwrap_or(0);
    let autocorr_time = recursive_stack(obj.get("autocorr_time").cloned().unwrap_or(Json::Null));
    let covariance = obj.get("covariance").cloned().map(recursive_stack);

    Some((mean, error, rebin_len, autocorr_time, covariance))
}

/// Load a Carlo `*.results.json` file and flatten into rows.
///
/// Reads the concatenated results file (an array of task result objects),
/// and produces one [`ResultRow`] per observable per task.
///
/// Each task result object has the format (matching Carlo.jl):
/// ```json
/// {
///   "task": "task_name",
///   "parameters": {"L": 10, "beta": 1.0},
///   "results": {
///     "Energy": {"mean": -1.23, "error": 0.01, "rebin_len": 100, "autocorr_time": 2.5}
///   }
/// }
/// ```
///
/// # Returns
///
/// A vector of [`ResultRow`]s, one per observable per task.
/// Suitable for conversion to a DataFrame (e.g. via `polars` or `DataFrame` crate).
pub fn dataframe(path: &Path) -> Result<Vec<ResultRow>, crate::CarloError> {
    let content = std::fs::read_to_string(path).map_err(|e| crate::CarloError::IoError {
        path: path.to_path_buf(),
        source: e,
    })?;

    let tasks: Vec<Json> =
        serde_json::from_str(&content).map_err(crate::CarloError::SerializationError)?;

    let mut rows = Vec::new();

    for task in tasks {
        let task_name = task["task"].as_str().unwrap_or("unknown").to_string();
        let parameters = task["parameters"].as_object().cloned().unwrap_or_default();
        let results = match task["results"].as_object() {
            Some(r) => r,
            None => continue,
        };

        for (obs_name, obs_value) in results {
            if obs_value.is_null() {
                continue;
            }

            let (mean, error, rebin_len, autocorr_time, covariance) =
                match measurement_from_obs(obs_value) {
                    Some(m) => m,
                    None => continue,
                };

            rows.push(ResultRow {
                task: task_name.clone(),
                parameters: parameters.clone(),
                observable: obs_name.clone(),
                mean,
                error,
                rebin_len,
                autocorr_time,
                covariance,
            });
        }
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_scalar_single_element() {
        let val = Json::Array(vec![Json::Number(
            serde_json::Number::from_f64(1.23).unwrap(),
        )]);
        let result = make_scalar_owned(val);
        assert_eq!(
            result,
            Json::Number(serde_json::Number::from_f64(1.23).unwrap())
        );
    }

    #[test]
    fn test_make_scalar_multi_element_unchanged() {
        let val = Json::Array(vec![
            Json::Number(serde_json::Number::from_f64(1.0).unwrap()),
            Json::Number(serde_json::Number::from_f64(2.0).unwrap()),
        ]);
        let result = make_scalar_owned(val.clone());
        assert_eq!(result, val);
    }

    #[test]
    fn test_recursive_stack_flat_unchanged() {
        let val = Json::Array(vec![
            Json::Number(serde_json::Number::from_f64(1.0).unwrap()),
            Json::Number(serde_json::Number::from_f64(2.0).unwrap()),
        ]);
        let result = recursive_stack(val.clone());
        assert_eq!(result, val);
    }

    #[test]
    fn test_recursive_stack_complex_object() {
        let val = serde_json::json!({"re": 1.0, "im": 0.5});
        let result = recursive_stack(val.clone());
        assert_eq!(result, val);
    }

    #[test]
    fn test_measurement_from_obs_scalar() {
        let obs = serde_json::json!({
            "mean": 1.23,
            "error": 0.05,
            "rebin_len": 100,
            "autocorr_time": 2.5
        });
        let (mean, error, rebin_len, autocorr_time, cov) = measurement_from_obs(&obs).unwrap();
        assert_eq!(mean, serde_json::json!(1.23));
        assert_eq!(error, serde_json::json!(0.05));
        assert_eq!(rebin_len, 100);
        assert_eq!(autocorr_time, serde_json::json!(2.5));
        assert!(cov.is_none());
    }

    #[test]
    fn test_measurement_from_obs_with_covariance() {
        let obs = serde_json::json!({
            "mean": [1.0, 2.0],
            "error": [0.1, 0.2],
            "rebin_len": 50,
            "autocorr_time": [1.5, 1.8],
            "covariance": [[1.0, 0.5], [0.5, 1.0]]
        });
        let (mean, _error, rebin_len, _autocorr_time, cov) = measurement_from_obs(&obs).unwrap();
        assert_eq!(mean, serde_json::json!([1.0, 2.0]));
        assert_eq!(rebin_len, 50);
        assert!(cov.is_some());
    }
}
