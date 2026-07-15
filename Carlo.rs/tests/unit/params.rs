use carlo_rs::lattice::LatticeParams;
use carlo_rs::{CarloError, Params};

#[test]
fn test_params_basic() {
    let mut params = Params::new();
    params.set("lattice_size", 32);
    params.set("temperature", 2.269);

    assert_eq!(params.get::<usize>("lattice_size"), Some(32));
    assert_eq!(params.get::<f64>("temperature"), Some(2.269));
    assert_eq!(params.get::<i32>("missing"), None);
}

#[test]
fn test_params_contains() {
    let mut params = Params::new();
    assert!(!params.contains("beta"));

    params.set("beta", 1.0);
    assert!(params.contains("beta"));
    assert!(!params.contains("gamma"));
}

#[test]
fn test_params_merge() {
    let mut base = Params::new();
    base.set("L", 16);
    base.set("beta", 0.5);

    let mut overlay = Params::new();
    overlay.set("beta", 1.0);
    overlay.set("J", 2.0);

    base.merge(&overlay);

    assert_eq!(base.get::<usize>("L"), Some(16));
    assert_eq!(base.get::<f64>("beta"), Some(1.0));
    assert_eq!(base.get::<f64>("J"), Some(2.0));
}

#[test]
fn test_params_merge_does_not_alter_source() {
    let mut overlay = Params::new();
    overlay.set("beta", 1.0);

    let mut target = Params::new();
    target.set("beta", 0.5);
    target.merge(&overlay);

    assert_eq!(overlay.get::<f64>("beta"), Some(1.0));
    assert_eq!(target.get::<f64>("beta"), Some(1.0));
}

#[test]
fn test_params_type_parsing() {
    let mut params = Params::new();
    params.set("count", 42u64);
    params.set("ratio", 0.75);
    params.set("flag", "true");

    assert_eq!(params.get::<u32>("count"), Some(42));
    assert_eq!(params.get::<f32>("ratio"), Some(0.75));
    assert_eq!(params.get::<String>("flag"), Some("true".into()));
    assert_eq!(params.get::<bool>("flag"), Some(true));
}

#[test]
fn test_params_unparseable_returns_none() {
    let mut params = Params::new();
    params.set("not_a_number", "abc");
    assert_eq!(params.get::<f64>("not_a_number"), None);
    assert_eq!(params.get::<usize>("not_a_number"), None);
}

#[test]
fn test_params_default() {
    let params = Params::default();
    assert!(!params.contains("anything"));
}

#[test]
fn test_params_serde_roundtrip() {
    let mut params = Params::new();
    params.set("L", 32);
    params.set("beta", 0.44);

    let json = serde_json::to_string(&params).unwrap();
    let restored: Params = serde_json::from_str(&json).unwrap();

    assert_eq!(restored.get::<usize>("L"), Some(32));
    assert_eq!(restored.get::<f64>("beta"), Some(0.44));
}

#[test]
fn test_params_equality() {
    let mut a = Params::new();
    a.set("L", 16);
    a.set("beta", 1.0);

    let mut b = Params::new();
    b.set("L", 16);
    b.set("beta", 1.0);

    assert_eq!(a, b);

    b.set("L", 32);
    assert_ne!(a, b);
}

#[test]
fn test_lattice_params() {
    let mut params = Params::new();
    params.set("lx", 16);
    params.set("ly", 16);

    let lattice = LatticeParams::from_params(&params).unwrap();
    assert_eq!(lattice.n_sites(), 256);
}

#[test]
fn test_lattice_params_ly_defaults_to_lx() {
    let mut params = Params::new();
    params.set("lx", 8);

    let lattice = LatticeParams::from_params(&params).unwrap();
    assert_eq!(lattice.lx, 8);
    assert_eq!(lattice.ly, 8);
    assert_eq!(lattice.n_sites(), 64);
}

#[test]
fn test_lattice_params_non_square() {
    let mut params = Params::new();
    params.set("lx", 4);
    params.set("ly", 6);

    let lattice = LatticeParams::from_params(&params).unwrap();
    assert_eq!(lattice.n_sites(), 24);
}

#[test]
fn test_lattice_params_missing_lx_errors() {
    let params = Params::new();
    let err = LatticeParams::from_params(&params).unwrap_err();
    assert!(matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "lx"));
}

#[test]
fn test_lattice_params_zero_dimension_errors() {
    let mut params = Params::new();
    params.set("lx", 0);

    let err = LatticeParams::from_params(&params).unwrap_err();
    assert!(matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "lattice"));

    params.set("lx", 4);
    params.set("ly", 0);
    let err = LatticeParams::from_params(&params).unwrap_err();
    assert!(matches!(err, CarloError::InvalidConfig { ref field, .. } if field == "lattice"));
}

#[test]
fn test_lattice_params_clone_debug() {
    let mut params = Params::new();
    params.set("lx", 3);
    let lattice = LatticeParams::from_params(&params).unwrap();
    let cloned = lattice.clone();
    assert_eq!(cloned.n_sites(), lattice.n_sites());
    let debug_str = format!("{:?}", lattice);
    assert!(debug_str.contains("LatticeParams"));
}
