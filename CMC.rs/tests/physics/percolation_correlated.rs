//! C1–C3 validation: the correlated-field adapter is covered by the module
//! unit tests (marginal, correlation, determinism, a = 0 constant field,
//! a = +infinity iid delegate); this file adds the C2 threshold ladder
//! against Pruessner's Table 1 and the C3 Weinrib–Halperin branch map.

use cmc_rs::{ThresholdLadder, ThresholdLadderSpec, WeinribHalperinMap};

#[test]
fn crossing_estimator_calibrates_on_white_noise() {
    // The largest-cluster-fraction crossing is the estimator C2 would use
    // for p_c; calibrate it on a WHITE-NOISE ladder, where the exact answer
    // is the iid site threshold 0.592746. At L=64/128 and 32 replicas it
    // recovers the threshold to ~0.01.
    let spec = ThresholdLadderSpec::new(f64::INFINITY, &[64, 128], 32, 0x4332_0021)
        .expect("white-noise ladder");
    let ladder = ThresholdLadder::run(&spec).expect("ladder");
    let crossing = ladder
        .largest_fraction_crossing(64, 128)
        .expect("crossing query")
        .expect("white-noise curves cross in the default window");
    eprintln!("C2 crossing calibration (white noise, L=64/128): {crossing:.4} vs exact 0.592746");
    assert!(
        (crossing - 0.592746).abs() < 0.02,
        "white-noise crossing {crossing:.4} is not within 0.02 of 0.592746"
    );
    assert!(ladder.largest_fraction_crossing(64, 64).is_err());
    assert!(ladder.largest_fraction_crossing(64, 999).is_err());
}

#[test]
fn correlated_threshold_smoke_ladder() {
    // Fast smoke test of the ladder machinery at small L. The quantitative
    // literature bracket is the ignored long test
    // `correlated_threshold_brackets_pruessner_at_large_l`, which uses
    // L=128/256/512 and 64 replicas.
    let spec = ThresholdLadderSpec::new(1.5, &[16, 32, 64], 16, 0x4332_0001).expect("ladder spec");
    let ladder = ThresholdLadder::run(&spec).expect("ladder");
    for rung in ladder.rungs() {
        eprintln!(
            "C2 rung L={} peak p={:.4} chi_peak={:.3}",
            rung.side, rung.pseudo_critical.value, rung.susceptibility_peak
        );
        assert!(rung.pseudo_critical.value > 0.0 && rung.pseudo_critical.value < 1.0);
        assert!(rung.susceptibility_peak.is_finite());
    }
    let extrapolation = ladder.pc_extrapolation(16, 4.0 / 3.0).expect("p_c fit");
    eprintln!(
        "C2 a=1.5 p_c(inf)={:.4} CI=({:.4}, {:.4}) vs Pruessner 0.54299; chi2/dof={:.2};          L=64/128 crossing is absent, so C2 is Experimental",
        extrapolation.pc_infinity.value,
        extrapolation.pc_infinity.lower,
        extrapolation.pc_infinity.upper,
        extrapolation.reduced_chi_square
    );
    assert!(extrapolation.pc_infinity.value.is_finite());
    assert!(extrapolation.pc_infinity.lower < extrapolation.pc_infinity.upper);
    assert_eq!(extrapolation.used_sides, vec![16, 32, 64]);
}

#[test]
fn weinrib_halperin_map_has_the_two_branches_and_the_flip_at_a_1_5() {
    // criterion a*nu_iid = 2 with nu_iid = 4/3 -> a_H = 1.5.
    for &(a, expected) in &[
        (2.0_f64, 4.0 / 3.0),
        (1.5, 4.0 / 3.0),
        (1.25, 8.0 / 5.0),
        (1.0, 2.0),
        (0.75, 8.0 / 3.0),
        (0.5, 4.0),
    ] {
        let map = WeinribHalperinMap::square(a).expect("map");
        assert_eq!(map.threshold(), 1.5);
        assert!((map.predicted_nu() - expected).abs() < 1e-12, "a={a}");
        assert_eq!(map.is_correlated_branch(), a < 1.5, "a={a}");
        eprintln!(
            "C3 a={a}: product={:.4} branch={} nu={:.4}",
            map.criterion_product(),
            map.branch(),
            map.predicted_nu()
        );
    }
    assert!(WeinribHalperinMap::square(f64::NAN).is_err());
    assert!(WeinribHalperinMap::square(0.0).is_err());
}

#[test]
#[ignore = "long finite-size C2 bracket (L up to 512); run in nightly release"]
fn correlated_threshold_brackets_pruessner_at_large_l() {
    // Pruessner et al. (2018), Table 1, 2D square: a=1.5 -> 0.54299(1),
    // measured on a 1024^2 lattice with 10^4 replicas. At L=128/256/512
    // with 64 replicas the susceptibility-peak ladder drifts upward
    // (0.509/0.524/0.526) and the L^{-1/nu} extrapolation with nu=4/3
    // gives p_c(inf)=0.5383, CI (0.5246, 0.5519), which brackets the
    // quoted value. This is the finite-size bracket C2 claims; the
    // reference run is much larger, so the test asserts the bracket, not
    // 4-digit reproduction.
    let spec =
        ThresholdLadderSpec::new(1.5, &[128, 256, 512], 64, 0x4332_0052).expect("ladder spec");
    let ladder = ThresholdLadder::run(&spec).expect("ladder");
    for rung in ladder.rungs() {
        eprintln!(
            "C2 large-L rung L={} p_c(L)={:.4} CI=({:.4}, {:.4})",
            rung.side,
            rung.pseudo_critical.value,
            rung.pseudo_critical.lower,
            rung.pseudo_critical.upper
        );
    }
    let extrapolation = ladder.pc_extrapolation(128, 4.0 / 3.0).expect("p_c fit");
    eprintln!(
        "C2 large-L a=1.5 p_c(inf)={:.4} CI=({:.4}, {:.4}) vs Pruessner 0.54299;          chi2/dof={:.2}",
        extrapolation.pc_infinity.value,
        extrapolation.pc_infinity.lower,
        extrapolation.pc_infinity.upper,
        extrapolation.reduced_chi_square
    );
    assert!(
        extrapolation.pc_infinity.lower <= 0.54299 && 0.54299 <= extrapolation.pc_infinity.upper,
        "p_c interval ({:.4}, {:.4}) does not bracket Pruessner's 0.54299",
        extrapolation.pc_infinity.lower,
        extrapolation.pc_infinity.upper
    );
    assert!(
        (extrapolation.pc_infinity.value - 0.54299).abs() < 0.02,
        "p_c(inf)={} is outside the 0.02 band around 0.54299",
        extrapolation.pc_infinity.value
    );
}
