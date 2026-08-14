//! Physics validation of the parallel-tempering exchange acceptance rule.
//!
//! A replica exchange between two chains with combined log-weight ratio Δ is
//! accepted with Metropolis probability min(1, e^Δ). These tests pin down the
//! ingredients of that physics without MPI:
//!
//! 1. `accept_log_probability` is deterministic at the decision boundaries:
//!    Δ ≥ 0 always accepts (including Δ = 0 and Δ = f64::EPSILON), Δ = −∞
//!    never accepts, and NaN input is rejected (documented behavior).
//! 2. For a fixed negative log-ratio the empirical acceptance frequency over
//!    many seeded draws matches exp(Δ) within 5σ binomial error bars.
//! 3. `ParallelTemperingMC::set_chain_idx` propagates the new parameter value
//!    to the child MC through `change_parameter` (observable in child state),
//!    and `current_value()` follows `parameter_values[chain_idx]`.

use carlo_rs::accept_log_probability;
use carlo_rs::parallel_tempering::{ParallelTemperingConfig, ParallelTemperingMC};
use carlo_rs::{Context, MonteCarlo, ParallelTemperingCompatible};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ── Toy PT-compatible model ───────────────────────────────────────────────

/// Two-state model with weight W(x | β) = exp(−β·x), so π_β(x = 1) = 1/(1+e^β).
/// Used to observe `change_parameter` propagation through the PT wrapper.
struct TwoStateMc {
    beta: f64,
    x: u8,
}

impl MonteCarlo for TwoStateMc {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {}

    fn name(&self) -> &'static str {
        "TwoStateMc"
    }
}

impl ParallelTemperingCompatible for TwoStateMc {
    /// log[W(x|new)/W(x|cur)] = (β_cur − β_new)·x.
    fn log_weight_ratio(&self, _param: &str, new_value: f64) -> f64 {
        (self.beta - new_value) * f64::from(self.x)
    }

    fn change_parameter(&mut self, _param: &str, new_value: f64) {
        self.beta = new_value;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────

/// The log-Metropolis rule must be a step function at Δ = 0: any non-negative
/// log-ratio (up to and including +0.0 and the smallest positive subnormal
/// distance from zero) accepts unconditionally, and −∞ rejects unconditionally.
/// NaN is rejected per the documented contract.
#[test]
fn deterministic_acceptance_at_boundaries() {
    for seed in [0u64, 1, 7, 42, 2024, u64::MAX] {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for delta in [
            0.0,
            -0.0,
            f64::EPSILON,
            f64::MIN_POSITIVE,
            1.0,
            f64::INFINITY,
        ] {
            for _ in 0..1000 {
                assert!(
                    accept_log_probability(delta, &mut rng),
                    "delta = {delta} (seed {seed}) must always accept"
                );
            }
        }
        for _ in 0..1000 {
            assert!(
                !accept_log_probability(f64::NEG_INFINITY, &mut rng),
                "delta = -inf (seed {seed}) must never accept"
            );
        }
        for _ in 0..1000 {
            assert!(
                !accept_log_probability(f64::NAN, &mut rng),
                "NaN (seed {seed}) must be rejected"
            );
        }
    }
}

/// For Δ = ln(0.3) the acceptance probability is exactly e^Δ = 0.3. The
/// empirical rate over 200_000 seeded decisions must agree within 5 standard
/// deviations of the binomial error.
#[test]
fn acceptance_rate_matches_exp_delta() {
    let delta = 0.3f64.ln(); // ≈ -1.20397
    let n = 200_000u64;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(2024);

    let accepted = (0..n)
        .filter(|_| accept_log_probability(delta, &mut rng))
        .count();
    let empirical = accepted as f64 / n as f64;

    let sigma = (0.3 * 0.7 / n as f64).sqrt();
    assert!(
        (empirical - 0.3).abs() <= 5.0 * sigma,
        "empirical acceptance rate {empirical:.5} must match exp(delta) = 0.3 \
         within 5σ = {:.5} (accepted {accepted}/{n})",
        5.0 * sigma
    );
}

/// An accepted exchange calls `set_chain_idx`: the chain label moves, the
/// child MC must receive the new parameter value via `change_parameter`, and
/// `current_value()` must report `parameter_values[chain_idx]`.
#[test]
fn set_chain_idx_propagates_parameter_to_child() {
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.3, 0.7, 1.1],
        interval: 1,
    };
    let mc = TwoStateMc { beta: 0.3, x: 1 };
    let mut wrapper = ParallelTemperingMC::new(&config, 0, mc);

    // Initial state: chain 0 at beta = 0.3.
    assert_eq!(wrapper.chain_idx(), 0);
    assert_eq!(wrapper.current_value(), 0.3);
    assert_eq!(wrapper.child_mc.beta, 0.3);

    // Simulate an accepted exchange moving this replica to chain 2.
    wrapper.set_chain_idx(2);
    assert_eq!(wrapper.chain_idx(), 2);
    assert_eq!(wrapper.current_value(), 1.1);
    assert_eq!(
        wrapper.child_mc.beta, 1.1,
        "child MC must follow the new chain"
    );

    // ... and back down to chain 1 (exchanges move labels in both directions).
    wrapper.set_chain_idx(1);
    assert_eq!(wrapper.chain_idx(), 1);
    assert_eq!(wrapper.current_value(), 0.7);
    assert_eq!(wrapper.child_mc.beta, 0.7);
}
