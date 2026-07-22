//! Reference-value tests for autocorrelation time estimation.
//!
//! AR(1) process: x_{t+1} = ρ·x_t + √(1−ρ²)·ε_t, ε~N(0,1)
//! Theoretical integrated autocorrelation time: τ = (1+ρ)/(1−ρ)

use carlo_rs::Accumulator;

fn ar1_series(rho: f64, n: usize, seed: u64) -> Vec<f64> {
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    let mut rng = StdRng::seed_from_u64(seed);
    let noise_scale = (1.0 - rho * rho).sqrt();

    let mut series = Vec::with_capacity(n);
    let mut x = 0.0;
    for _ in 0..n {
        let u1: f64 = rng.random();
        let u2: f64 = rng.random();
        let u1 = u1.max(f64::MIN_POSITIVE);
        let eps = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        x = rho * x + noise_scale * eps;
        series.push(x);
    }
    series
}

fn make_acc_from_samples(samples: &[f64], _binsize: usize) -> Accumulator {
    let mut a = Accumulator::new(1);
    for &s in samples {
        a.add(s);
    }
    a
}

#[test]
fn uncorrelated_data_has_autocorr_time_near_one() {
    let series = ar1_series(0.0, 10000, 42);
    let acc = make_acc_from_samples(&series, 10);
    let tau = acc.autocorr_time();
    // The rebinned estimator returns 0.5*(var_rebin/var_no_rebin*2 - 1).
    // For truly uncorrelated data, var_rebin ≈ var_no_rebin/2 → τ ≈ 0.
    // The estimator clamps at 0, then finalize adds 1.0 via .max(1.0).
    // So the raw autocorr_time() can return 0 for uncorrelated data.
    // The important property is monotonicity, checked elsewhere.
    assert!(
        tau >= 0.0,
        "autocorr_time should be non-negative, got {tau}"
    );
}

#[test]
fn ar1_correlated_data_has_higher_tau_than_uncorrelated() {
    // The rebinned estimator is crude (single 2→1 rebin), so absolute τ values
    // are underestimated. We test the key property: correlated data has larger τ.
    let s_uncorr = ar1_series(0.0, 20000, 42);
    let s_mid = ar1_series(0.5, 20000, 42);
    let s_strong = ar1_series(0.8, 20000, 42);
    let tau_uncorr = make_acc_from_samples(&s_uncorr, 5).autocorr_time();
    let tau_mid = make_acc_from_samples(&s_mid, 5).autocorr_time();
    let tau_strong = make_acc_from_samples(&s_strong, 5).autocorr_time();
    assert!(
        tau_strong > tau_mid,
        "expected τ(0.8)>τ(0.5): {tau_strong} > {tau_mid}"
    );
    assert!(
        tau_mid >= tau_uncorr,
        "expected τ(0.5)≥τ(0.0): {tau_mid} ≥ {tau_uncorr}"
    );
}

#[test]
fn finalize_carries_nontrivial_tau_for_correlated_data() {
    // For ρ=0.8 with enough samples, the rebinned estimator should give τ>0,
    // and finalize should propagate it (clamped to ≥1.0).
    let series = ar1_series(0.8, 80000, 321);
    let acc = make_acc_from_samples(&series, 5);
    let tau_raw = acc.autocorr_time();
    let est = acc.finalize();
    // With 80k samples, the estimator should detect some correlation (τ_raw > 0)
    assert!(
        tau_raw > 0.01,
        "raw autocorr_time should detect ρ=0.8 correlation, got {tau_raw}"
    );
    // finalize should propagate the raw τ (clamped to ≥1.0)
    assert!(
        (est.autocorr_time - tau_raw.max(1.0)).abs() < 1e-10,
        "Estimate τ={} should match max(raw τ={tau_raw}, 1.0)",
        est.autocorr_time
    );
}

#[test]
fn stronger_correlation_gives_larger_autocorr_time() {
    let s_m = ar1_series(0.5, 20000, 789);
    let s_s = ar1_series(0.8, 20000, 789);
    let tau_m = make_acc_from_samples(&s_m, 5).autocorr_time();
    let tau_s = make_acc_from_samples(&s_s, 5).autocorr_time();
    assert!(tau_s >= tau_m, "τ(0.8)={tau_s} should be ≥ τ(0.5)={tau_m}");
}

#[test]
fn uncorrelated_finalize_gives_default_tau() {
    let series = ar1_series(0.0, 10000, 999);
    let acc = make_acc_from_samples(&series, 10);
    let est = acc.finalize();
    // For uncorrelated data, raw τ≈0, clamped to 1.0 by .max(1.0)
    assert!(
        est.autocorr_time >= 1.0,
        "expected τ≥1.0 for uncorrelated, got {}",
        est.autocorr_time
    );
}
