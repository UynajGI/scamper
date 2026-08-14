//! ESS calibration against AR(1) reference chains (tracker P0.5).
//!
//! The multi-chain ESS estimator behind `diagnose` was previously calibrated
//! only on IID draws. A stationary Gaussian AR(1) chain
//! `x_{t+1} = ρ·x_t + sqrt(1-ρ²)·ε_t` has an exact N(0,1) marginal,
//! integrated autocorrelation time τ = (1+ρ)/(1-ρ) and therefore a
//! theoretical ESS/N = (1-ρ)/(1+ρ). Because the marginal is exactly standard
//! normal, the rank normalization applied inside `diagnose` is an identity
//! transform in distribution, so the closed-form ESS carries over to
//! `ess_bulk`.

use mcmc_rs::proposal::standard_normal;
use mcmc_rs::{diagnose, EuclideanState, FnLogDensity, MemoryTrace, TraceStore, TransitionReport};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const CHAINS: usize = 4;
/// Burn-in decays the initial condition to insignificance even at ρ = 0.99
/// (0.99^1000 ≈ 4·10^-5).
const BURN_IN: usize = 1_000;

fn ar1_trace(chain_id: usize, rho: f64, draws: usize, seed: u64) -> MemoryTrace {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let innovation = (1.0 - rho * rho).sqrt();
    let mut value = 0.0;
    for _ in 0..BURN_IN {
        value = rho * value + innovation * standard_normal(&mut rng);
    }

    // dummy standard-normal target; only the recorded positions are used
    let mut target = FnLogDensity::new(|position: &[f64]| -0.5 * position[0] * position[0]);
    let mut state = EuclideanState::initialize(&mut target, vec![0.0]).unwrap();
    let mut trace = MemoryTrace::new(1, 1).unwrap();
    let report = TransitionReport {
        accepted: Some(true),
        proposals: 1,
        acceptances: 1,
        ..TransitionReport::default()
    };
    for _ in 0..draws {
        value = rho * value + innovation * standard_normal(&mut rng);
        state.replace(vec![value], -0.5 * value * value);
        trace.record(chain_id, &state, &report).unwrap();
    }
    trace
}

/// Estimated bulk ESS of AR(1) chains must track the closed-form
/// ESS = N·(1-ρ)/(1+ρ), with relative tolerances that widen with ρ.
///
/// Tolerance basis (measured across seeds with this estimator): the raw
/// Geyer initial-monotone estimator is unbiased to within ±2% at these
/// sizes, while the rank-normalize-and-split pipeline inside `diagnose`
/// adds a small conservative ESS underestimate for short-memory chains
/// that shrinks with N — hence 4×25k draws at ρ = 0.5 to hold a 5% band
/// (measured 0.4-2.6%; ρ = 0 returns ESS = N exactly). Long-memory chains
/// leave the truncation-lag choice of the paired-autocorrelation walk more
/// noisy: at ρ = 0.9 (τ = 19, measured 1-7% over 4×10k draws) a 15% band
/// applies, and at ρ = 0.99 (τ = 199, only ≈200 effective draws remain out
/// of 40k, measured 9-15%) a 30% band. Draw counts stay near the smallest
/// size meeting each band because the estimator costs O(draws²) per chain.
#[test]
fn ar1_chains_recover_theoretical_ess() {
    let cases = [
        (0.0f64, 10_000usize, 0.05f64),
        (0.5, 25_000, 0.05),
        (0.9, 10_000, 0.15),
        (0.99, 10_000, 0.30),
    ];
    for (case_index, &(rho, draws, tolerance)) in cases.iter().enumerate() {
        let traces = (0..CHAINS)
            .map(|chain| {
                ar1_trace(
                    chain,
                    rho,
                    draws,
                    9_241 + 1_031 * case_index as u64 + chain as u64,
                )
            })
            .collect::<Vec<_>>();
        let diagnostics = diagnose(&traces, &["x".to_string()]).unwrap();
        let ess = diagnostics.parameters[0].ess_bulk;
        let total = (CHAINS * draws) as f64;
        let theoretical = total * (1.0 - rho) / (1.0 + rho);
        let relative_error = ((ess - theoretical) / theoretical).abs();
        assert!(
            relative_error <= tolerance,
            "rho={rho}: ess={ess}, theoretical={theoretical}, \
             relative error={relative_error}"
        );
    }
}
