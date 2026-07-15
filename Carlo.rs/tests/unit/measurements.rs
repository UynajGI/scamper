use carlo_rs::Measurements;

#[test]
fn test_binning_accumulation() {
    let mut meas = Measurements::new(10);

    // Add 25 samples (should create 2 full bins + 5 partial)
    for i in 0..25 {
        meas.add_sample("Energy", i as f64);
    }

    // After finalize, should have 2 bins
    let results = meas.finalize();
    let estimate = results.get("Energy").expect("Energy observable");

    // Mean of bins: [0..10].mean()=4.5, [10..20].mean()=14.5 → overall ~9.5
    assert!(estimate.mean > 0.0);
    assert!(estimate.stderr > 0.0);
}

#[test]
fn test_measurements_register_scalar() {
    let mut meas = Measurements::new(5);
    meas.register("Magnetization", 3);

    meas.add_sample("Magnetization", 0.5);
    meas.add_sample("Magnetization", 0.7);
    meas.add_sample("Magnetization", 0.9);

    let results = meas.finalize();
    assert!(results.contains_key("Magnetization"));
    let est = &results["Magnetization"];
    assert!((est.mean - 0.7).abs() < 1e-10);
}

#[test]
fn test_measurements_add_sample_array() {
    let mut meas = Measurements::new(2);
    meas.add_sample_array("Correlation", &[1.0, 2.0, 3.0]);
    meas.add_sample_array("Correlation", &[4.0, 5.0, 6.0]);
    meas.add_sample_array("Correlation", &[7.0, 8.0, 9.0]);
    meas.add_sample_array("Correlation", &[10.0, 11.0, 12.0]);

    let results = meas.finalize();
    assert!(results.contains_key("Correlation"));
    let est = &results["Correlation"];
    assert_eq!(est.n_bins, 2);
}

#[test]
fn test_measurements_register_array() {
    let mut meas = Measurements::new(1);
    meas.register_array("SpinProfile", 1, &[4]);

    meas.add_sample_array("SpinProfile", &[1.0, 2.0, 3.0, 4.0]);
    meas.add_sample_array("SpinProfile", &[5.0, 6.0, 7.0, 8.0]);

    let results = meas.finalize();
    assert!(results.contains_key("SpinProfile"));
}

#[test]
fn test_measurements_register_complex() {
    let mut meas = Measurements::new(1);
    meas.register_complex("OrderParam", 2);

    meas.add_sample_complex("OrderParam", 1.0, 2.0);
    meas.add_sample_complex("OrderParam", 3.0, 4.0);

    let complex = meas.finalize_complex();
    assert!(complex.contains_key("OrderParam"));
    let est = &complex["OrderParam"];
    assert!((est.re.mean - 2.0).abs() < 1e-10);
    assert!((est.im.mean - 3.0).abs() < 1e-10);
}

#[test]
fn test_measurements_observables_accessor() {
    let mut meas = Measurements::new(5);
    meas.add_sample("Energy", 1.0);
    meas.add_sample("Mag", 0.5);

    let observables = meas.observables();
    assert!(observables.contains_key("Energy"));
    assert!(observables.contains_key("Mag"));
    assert_eq!(observables.len(), 2);
}

#[test]
fn test_measurements_complex_observables_accessor() {
    let mut meas = Measurements::new(5);
    meas.add_sample_complex("GreenFunc", 1.0, 0.5);

    let complex = meas.complex_observables();
    assert!(complex.contains_key("GreenFunc"));
    assert_eq!(complex.len(), 1);
}

#[test]
fn test_measurements_mixed_scalar_and_complex() {
    let mut meas = Measurements::new(1);

    meas.add_sample("Energy", -1.5);
    meas.add_sample("Energy", -2.0);
    meas.add_sample_complex("Susceptibility", 0.3, 0.4);
    meas.add_sample_complex("Susceptibility", 0.5, 0.6);

    let results = meas.finalize();
    assert!(results.contains_key("Energy"));

    let complex = meas.finalize_complex();
    assert!(complex.contains_key("Susceptibility"));
}

#[test]
fn test_measurements_multiple_observables_finalize() {
    let mut meas = Measurements::new(2);

    for i in 0..10 {
        meas.add_sample("A", i as f64);
        meas.add_sample("B", (i as f64) * 2.0);
        meas.add_sample("C", (i as f64) * 0.5);
    }

    let results = meas.finalize();
    assert_eq!(results.len(), 3);
    for key in &["A", "B", "C"] {
        assert!(results.contains_key(*key));
    }

    let b_est = &results["B"];
    let a_est = &results["A"];
    assert!((b_est.mean - 2.0 * a_est.mean).abs() < 1e-10);
}

#[test]
fn test_measurements_empty_finalize() {
    let meas = Measurements::new(5);
    let results = meas.finalize();
    assert!(results.is_empty());

    let complex = meas.finalize_complex();
    assert!(complex.is_empty());
}
