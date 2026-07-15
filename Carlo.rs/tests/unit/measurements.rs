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
