use carlo_rs::lattice::LatticeParams;
use carlo_rs::Params;

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
fn test_lattice_params() {
    let mut params = Params::new();
    params.set("lx", 16);
    params.set("ly", 16);

    let lattice = LatticeParams::from_params(&params).unwrap();
    assert_eq!(lattice.n_sites(), 256);
}
