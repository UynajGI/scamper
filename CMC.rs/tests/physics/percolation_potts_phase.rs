//! Q3 validation: the exact 2D FK/Potts self-dual critical curve and an
//! exploratory finite-size phase scan at that curve.
//!
//! The critical probabilities are ANCHOR INPUTS (`sqrt(q)/(1+sqrt(q))`),
//! never fitted. The scan reports the bond-density mean/variance and the
//! descriptive bimodality coefficient; the q=2 vs q=8 ordering is recorded
//! as exploratory evidence, not a latent-heat claim.

use cmc_rs::{critical_probability, phase_scan, PottsPhaseSpec};

#[test]
fn exact_self_dual_critical_curve() {
    for &(q, expected) in &[
        (1.0_f64, 0.5_f64),
        (
            2.0,
            std::f64::consts::SQRT_2 / (1.0 + std::f64::consts::SQRT_2),
        ),
        (4.0, 2.0 / 3.0),
        (8.0, 8.0_f64.sqrt() / (1.0 + 8.0_f64.sqrt())),
    ] {
        let value = critical_probability(q).expect("critical probability");
        assert!((value - expected).abs() < 1e-15, "q={q}");
        eprintln!("Q3 p_c(q={q}) = {value:.6}");
    }
    assert!(critical_probability(0.0).is_err());
    assert!(critical_probability(f64::INFINITY).is_err());
}

#[test]
fn phase_scan_reports_bond_density_and_bimodality_at_the_exact_curve() {
    let mut estimates = Vec::new();
    for &q in &[2.0_f64, 8.0] {
        let spec = PottsPhaseSpec::new(q, 8, 200, 800, 0x5133_0001 + q as u64).expect("spec");
        let estimate = phase_scan(spec).expect("scan");
        eprintln!(
            "Q3 q={q} L=8 p_c={:.4}: density mean={:.4} var={:.5} skew={:.3} kurt={:.3} \
             bimodality={:.3} (exploratory)",
            estimate.critical_probability,
            estimate.density_mean,
            estimate.density_variance,
            estimate.skewness,
            estimate.kurtosis,
            estimate.bimodality
        );
        assert!(estimate.density_mean > 0.0 && estimate.density_mean < 1.0);
        assert!(estimate.density_variance > 0.0);
        estimates.push(estimate);
    }
    // The exact curve is increasing in q; the scan is anchored at it.
    assert!(
        estimates[1].critical_probability > estimates[0].critical_probability,
        "the exact p_c(q) curve must increase"
    );
    // Exploratory only: the finite-size bimodality coefficient is reported
    // for q=2 and q=8; no ordering or threshold is asserted at L=8.
}

#[test]
#[ignore = "long exploratory Q3 scan at L=16"]
fn q_eight_bimodality_ordering_at_l16() {
    let mut estimates = Vec::new();
    for &q in &[2.0_f64, 4.0, 6.0, 8.0] {
        let spec = PottsPhaseSpec::new(q, 16, 400, 3000, 0x5133_1001 + q as u64).expect("spec");
        let estimate = phase_scan(spec).expect("scan");
        eprintln!(
            "Q3 L=16 q={q} p_c={:.4}: mean={:.4} var={:.5} skew={:.3} kurt={:.3} bimodality={:.3}",
            estimate.critical_probability,
            estimate.density_mean,
            estimate.density_variance,
            estimate.skewness,
            estimate.kurtosis,
            estimate.bimodality
        );
        estimates.push(estimate);
    }
    // Descriptive classification: below the 5/9 bimodality threshold for
    // q <= 4 (continuous) and for q=6 at this L; above it only for q=8.
    const BIMODAL: f64 = 5.0 / 9.0;
    for estimate in &estimates[..2] {
        assert!(
            estimate.bimodality < BIMODAL,
            "q={} bimodality {:.3} should stay below 5/9 at L=16",
            estimate.q,
            estimate.bimodality
        );
    }
    // q=6 sits below 5/9 at L=16 (0.519): the first-order signal for the
    // intermediate case needs a larger lattice, so it is reported, not
    // asserted. q=8 is clearly bimodal at this size.
    assert!(
        estimates[3].bimodality > BIMODAL,
        "q=8 bimodality {:.3} should exceed 5/9 at L=16",
        estimates[3].bimodality
    );
    for pair in estimates.windows(2) {
        assert!(
            pair[1].bimodality > pair[0].bimodality,
            "bimodality must increase with q: {:.3} then {:.3}",
            pair[0].bimodality,
            pair[1].bimodality
        );
    }
    eprintln!(
        "Q3 q=6 bimodality {:.3} is below 5/9 at L=16 (intermediate case unresolved)",
        estimates[2].bimodality
    );
}

#[test]
fn two_mode_estimator_separates_first_order_from_continuous_cases() {
    // Exploratory q-ladder at L=12: q=1/2 are continuous, q=8 is first
    // order. The two-mode separation is a finite-size latent-heat proxy;
    // the gate asserts the expected ordering at this one size and prints
    // the L ladder, not a thermodynamic claim.
    let mut rows = Vec::new();
    for &q in &[1.0_f64, 2.0, 4.0, 8.0] {
        let spec = PottsPhaseSpec::new(q, 12, 400, 2000, 0x5133_2001 + q as u64).expect("spec");
        let scan = cmc_rs::phase_scan_samples(spec).expect("scan");
        let modes = cmc_rs::two_mode_estimate(&scan.densities).expect("two-mode split");
        eprintln!(
            "Q3 two-mode L=12 q={q}: gap={:.4} separation={:.3} weights=({:.2}, {:.2}) mean={:.4}",
            modes.gap,
            modes.separation,
            modes.lower_weight,
            1.0 - modes.lower_weight,
            scan.estimate.density_mean
        );
        rows.push(modes);
    }
    assert!(
        rows[3].separation > rows[1].separation,
        "q=8 separation {:.3} must exceed q=2 {:.3}",
        rows[3].separation,
        rows[1].separation
    );
    assert!(
        rows[3].gap > rows[1].gap,
        "q=8 gap {:.4} must exceed q=2 {:.4}",
        rows[3].gap,
        rows[1].gap
    );
    // q=1 (bond percolation at p=0.5) must not look more bimodal than q=8.
    assert!(rows[3].separation > rows[0].separation);
}
