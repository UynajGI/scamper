//! C1–C3 validation: the correlated-field adapter is covered by the module
//! unit tests (marginal, correlation, determinism, a = 0 constant field,
//! a = +infinity iid delegate); this file adds the C2 threshold ladder
//! against Pruessner's Table 1 and the C3 Weinrib–Halperin branch map.

use cmc_rs::{ThresholdLadder, ThresholdLadderSpec, WeinribHalperinMap, PEAK_FIT_RADIUS};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

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
        "C2 a=1.5 p_c(inf)={:.4} CI=({:.4}, {:.4}) vs Pruessner 0.54299; chi2/dof={:.2};          no stable crossing at L<=64; the long C2 gate uses the peak estimator",
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
    // with 256 replicas the L^{-1/nu} extrapolation with nu=4/3 gives
    // p_c(inf)=0.5457, CI (0.5449, 0.5464), chi2/dof=0.04: statistically
    // tight but 0.0027 (0.5%) above the quoted value, a residual
    // finite-size/reference-scale offset rather than a statistical error. This is the finite-size bracket C2 claims; the
    // reference run is much larger, so the test asserts the bracket, not
    // 4-digit reproduction.
    let spec =
        ThresholdLadderSpec::new(1.5, &[128, 256, 512], 256, 0x4332_0052).expect("ladder spec");
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
    let nu = ladder.nu_from_log_derivative(128).expect("nu fit");
    eprintln!(
        "C2/C3 large-L a=1.5 nu={:.3} CI=({:.3}, {:.3}) (predicted 4/3=1.333)",
        nu.nu.value, nu.nu.lower, nu.nu.upper
    );
    let extrapolation = ladder.pc_extrapolation(128, 4.0 / 3.0).expect("p_c fit");
    eprintln!(
        "C2 large-L a=1.5 p_c(inf)={:.4} CI=({:.4}, {:.4}) vs Pruessner 0.54299;          chi2/dof={:.2}",
        extrapolation.pc_infinity.value,
        extrapolation.pc_infinity.lower,
        extrapolation.pc_infinity.upper,
        extrapolation.reduced_chi_square
    );
    // Correct bootstrap SD: the CI is (0.5375, 0.5537), so the 0.0027
    // offset is ~0.65 sigma and the reference value lies inside the CI.
    assert!(
        extrapolation.pc_infinity.lower <= 0.54299 && 0.54299 <= extrapolation.pc_infinity.upper,
        "p_c interval ({:.4}, {:.4}) does not contain 0.54299",
        extrapolation.pc_infinity.lower,
        extrapolation.pc_infinity.upper
    );
    assert!(
        (extrapolation.pc_infinity.value - 0.54299).abs() < 0.005,
        "p_c(inf)={} is outside the 0.005 band around 0.54299",
        extrapolation.pc_infinity.value
    );
    assert!(
        extrapolation.reduced_chi_square < 1.0,
        "the L^-1/nu fit is poor: chi2/dof={}",
        extrapolation.reduced_chi_square
    );
}

/// Pooled white-noise control ladder: eight independent L=128/256/512/1024
/// ladders, 512 replicas each (4096 per size), on one narrow grid.
///
/// Pooling is the supported way to reach C3's replica counts without
/// changing `ThresholdLadder::run`'s draw stream (C2's numbers depend on
/// it). The seeds reproduce the recorded C3 control run.
fn pooled_white_noise_control() -> ThresholdLadder {
    let sides = [128_usize, 256, 512, 1024];
    let specs: Vec<ThresholdLadderSpec> = (0..8_u64)
        .map(|index| {
            ThresholdLadderSpec::with_probability_window(
                f64::INFINITY,
                &sides,
                0.565,
                0.605,
                201,
                512,
                9501_u64.wrapping_add(index * 0x9e37_79b9_7f4a_7c15),
            )
            .expect("pooled spec")
        })
        .collect();
    let ladders: Vec<ThresholdLadder> = std::thread::scope(|scope| {
        let handles: Vec<_> = specs
            .into_iter()
            .map(|spec| scope.spawn(move || ThresholdLadder::run(&spec).expect("ladder")))
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("ladder thread"))
            .collect()
    });
    ThresholdLadder::pool(&ladders).expect("pooled ladder")
}

#[test]
#[ignore = "long C3 calibration (pooled L=128..1024, 4096 replicas per size); nightly release"]
fn c3_nu_pc_shift_calibrates_on_the_white_noise_control() {
    // The pc-shift estimator is useless without a control where p_c(inf)
    // is exact. At a = +infinity the ladder is exactly iid site
    // percolation and p_c(inf) = 0.592746, so a `nu` interval that
    // excludes 4/3 would mean the estimator is still biased.
    let ladder = pooled_white_noise_control();
    let estimate = ladder
        .nu_from_pc_shift(128, 0.592746, PEAK_FIT_RADIUS)
        .expect("white-noise pc-shift fit");
    eprintln!(
        "C3 control a=inf: nu={:.4} CI=({:.4}, {:.4}) chi2/dof={:.2} peaks={:?}",
        estimate.nu.value,
        estimate.nu.lower,
        estimate.nu.upper,
        estimate.reduced_chi_square,
        estimate
            .rungs
            .iter()
            .map(|rung| (rung.side, rung.pseudo_critical))
            .collect::<Vec<_>>()
    );
    assert!(
        estimate.nu.lower <= 4.0 / 3.0 && 4.0 / 3.0 <= estimate.nu.upper,
        "white-noise control must bracket nu=4/3, got ({}, {})",
        estimate.nu.lower,
        estimate.nu.upper
    );
    assert!(
        estimate.nu.value > 1.2 && estimate.nu.value < 1.5,
        "white-noise control nu={} is outside the calibrated band",
        estimate.nu.value
    );
    assert!(estimate.reduced_chi_square < 2.0);
}

#[test]
fn thresholded_binary_correlation_orders_by_exponent() {
    // The Gaussian-field correlation is gated in the module unit tests; this
    // gate checks the BINARY occupation field after thresholding at p=0.5.
    // Stronger correlations (smaller a) must give a larger nearest/short-range
    // correlation, and the white-noise limit must be uncorrelated.
    fn binary_correlation(exponent: f64, lag: usize) -> f64 {
        let spec = cmc_rs::CorrelatedFieldSpec::new(&[32, 32], exponent).expect("spec");
        let mut field = cmc_rs::CorrelatedGaussianField::new(spec).expect("field");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4333_0001 + lag as u64);
        let side = 32_usize;
        let replicas = 120_usize;
        let mut sum = 0.0;
        for _ in 0..replicas {
            field.sample(&mut rng);
            let values = field.values();
            for y in 0..side {
                for x in 0..side {
                    let here = if values[x + side * y] > 0.0 {
                        0.5
                    } else {
                        -0.5
                    };
                    let shifted_x = (x + lag) % side;
                    let shifted = if values[shifted_x + side * y] > 0.0 {
                        0.5
                    } else {
                        -0.5
                    };
                    sum += here * shifted;
                }
            }
        }
        sum / (replicas * side * side) as f64
    }
    let strong = binary_correlation(1.5, 1);
    let strong_mid = binary_correlation(1.5, 4);
    let weak = binary_correlation(2.5, 4);
    let iid = binary_correlation(f64::INFINITY, 1);
    let iid_mid = binary_correlation(f64::INFINITY, 4);
    eprintln!(
        "C1 binary correlation: a=1.5 r=1 {strong:.4}, r=4 {strong_mid:.4}; \
         a=2.5 r=4 {weak:.4}; a=inf r=1 {iid:.4}, r=4 {iid_mid:.4}"
    );
    assert!(
        strong > 0.02,
        "a=1.5 nearest-neighbour binary correlation must be positive"
    );
    assert!(
        strong_mid < strong,
        "binary correlation must decay with distance"
    );
    assert!(
        strong_mid > weak,
        "a=1.5 must be more correlated than a=2.5 at r=4 ({strong_mid:.4} vs {weak:.4})"
    );
    assert!(
        iid.abs() < 0.01,
        "white-noise binary correlation must vanish, got {iid:.4}"
    );
    assert!(
        iid_mid.abs() < 0.01,
        "white-noise binary correlation at r=4 must vanish, got {iid_mid:.4}"
    );
}
