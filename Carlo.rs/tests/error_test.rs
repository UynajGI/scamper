use carlo_rs::CarloError;

#[test]
fn test_error_display() {
    let err = CarloError::InvalidConfig {
        field: "binsize".into(),
        reason: "must be positive".into(),
    };
    assert!(err.to_string().contains("binsize"));
}
