use carlo_rs::CarloError;
use std::error::Error;
use std::io;

#[test]
fn test_error_display() {
    let err = CarloError::InvalidConfig {
        field: "binsize".into(),
        reason: "must be positive".into(),
    };
    assert!(err.to_string().contains("binsize"));
    assert!(err.to_string().contains("must be positive"));
}

#[test]
fn test_io_error_display() {
    let err = CarloError::IoError {
        path: "/tmp/missing.json".into(),
        source: io::Error::new(io::ErrorKind::NotFound, "file missing"),
    };
    let msg = err.to_string();
    assert!(msg.contains("/tmp/missing.json"));
    assert!(msg.contains("file missing"));
}

#[test]
fn test_measurement_not_found_display() {
    let err = CarloError::MeasurementNotFound {
        name: "SpinCorrelation".into(),
    };
    assert!(err.to_string().contains("SpinCorrelation"));
}

#[test]
fn test_checkpoint_corrupted_display() {
    let err = CarloError::CheckpointCorrupted {
        detail: "invalid sweep_count field".into(),
    };
    assert!(err.to_string().contains("invalid sweep_count field"));
}

#[test]
fn test_convergence_timeout_display() {
    let err = CarloError::ConvergenceTimeout { sweeps: 50_000 };
    let msg = err.to_string();
    assert!(msg.contains("50000"));
    assert!(msg.contains("sweeps"));
}

#[test]
fn test_serialization_error_from_serde_json() {
    let json_err = serde_json::from_str::<serde_json::Value>("{bad}").unwrap_err();
    let err = CarloError::SerializationError(json_err);
    assert!(err.to_string().contains("Serialization error"));
}

#[test]
fn test_invalid_config_source_chain() {
    let err = CarloError::InvalidConfig {
        field: "beta".into(),
        reason: "must be positive".into(),
    };
    assert!(err.source().is_none());
}

#[test]
fn test_io_error_source_is_propagated() {
    let source = io::Error::new(io::ErrorKind::PermissionDenied, "access denied");
    let err = CarloError::IoError {
        path: "/root/secret".into(),
        source,
    };
    assert!(err.source().is_some());
    assert!(err.source().unwrap().to_string().contains("access denied"));
}

#[test]
fn test_serialization_error_source_is_propagated() {
    let json_err = serde_json::from_str::<serde_json::Value>("not json").unwrap_err();
    let err = CarloError::SerializationError(json_err);
    assert!(err.source().is_some());
}
