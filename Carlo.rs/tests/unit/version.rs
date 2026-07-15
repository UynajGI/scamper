use carlo_rs::Version;

#[test]
fn test_version_new_with_mc_version() {
    let v = Version::new(Some("ising-1.2"));
    assert!(!v.carlo_version.is_empty());
    assert_eq!(v.mc_version.as_deref(), Some("ising-1.2"));
    assert!(v.rng_version > 0);
}

#[test]
fn test_version_new_without_mc_version() {
    let v = Version::new(None);
    assert!(!v.carlo_version.is_empty());
    assert!(v.mc_version.is_none());
    assert!(v.rng_version > 0);
}

#[test]
fn test_version_current() {
    let v = Version::current();
    assert!(!v.carlo_version.is_empty());
    assert!(v.mc_version.is_none());
}

#[test]
fn test_version_carlo_matches_crate_version() {
    let v = Version::current();
    assert_eq!(v.carlo_version, env!("CARGO_PKG_VERSION"));
}

#[test]
fn test_version_serde_roundtrip() {
    let v = Version::new(Some("test-model-0.3"));
    let json = serde_json::to_string(&v).unwrap();
    let restored: Version = serde_json::from_str(&json).unwrap();

    assert_eq!(restored.carlo_version, v.carlo_version);
    assert_eq!(restored.mc_version, v.mc_version);
    assert_eq!(restored.rng_version, v.rng_version);
}

#[test]
fn test_version_serde_roundtrip_no_mc_version() {
    let v = Version::new(None);
    let json = serde_json::to_string(&v).unwrap();
    let restored: Version = serde_json::from_str(&json).unwrap();
    assert_eq!(restored.carlo_version, v.carlo_version);
    assert!(restored.mc_version.is_none());
}

#[test]
fn test_version_debug_format() {
    let v = Version::new(Some("ising"));
    let debug = format!("{:?}", v);
    assert!(debug.contains("Version"));
    assert!(debug.contains("carlo_version"));
}
