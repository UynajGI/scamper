//! C1–C3 validation: the correlated-field adapter is covered by the module
//! unit tests (marginal, correlation, determinism, a = 0 constant field,
//! a = +infinity iid delegate); this file adds the C2 threshold ladder
//! against Pruessner's Table 1 and the C3 Weinrib–Halperin branch map.

use cmc_rs::{ThresholdLadder, ThresholdLadderSpec, WeinribHalperinMap};

#[test]
fn correlated_threshold_ladder_reports_its_current_limitation() {
    // Pruessner et al. (2018), Table 1, 2D square: a=1.5 -> 0.54299(1).
    // The delivered estimator is NOT yet able to reproduce that value at
    // feasible sizes: the susceptibility-peak estimate is noisy and the
    // L^{-1/nu} fit has a large reduced chi-square. This gate therefore
    // validates the MACHINERY (finite estimates, ordered rungs, a usable
    // interval) and records the honest comparison rather than asserting a
    // four-digit reproduction. C2 remains Experimental pending a better
    // finite-size estimator.
    let spec = ThresholdLadderSpec::new(1.5, &[16, 32, 64], 16, 0x4332_0001).expect("ladder spec");
    let ladder = ThresholdLadder::run(&spec).expect("ladder");
    for rung in ladder.rungs() {
        eprintln!(
            "C2 rung L={} p_c(L)={:.4} CI=({:.4}, {:.4}) chi_peak={:.3}",
            rung.side,
            rung.pseudo_critical.value,
            rung.pseudo_critical.lower,
            rung.pseudo_critical.upper,
            rung.susceptibility_peak
        );
        assert!(rung.pseudo_critical.value > 0.0 && rung.pseudo_critical.value < 1.0);
        assert!(rung.susceptibility_peak.is_finite());
    }
    let extrapolation = ladder.pc_extrapolation(16, 4.0 / 3.0).expect("p_c fit");
    eprintln!(
        "C2 a=1.5: p_c(inf)={:.4} CI=({:.4}, {:.4}) vs Pruessner 0.54299;          chi2/dof={:.2} (the current estimator is Experimental)",
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
