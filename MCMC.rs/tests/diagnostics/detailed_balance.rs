//! Detailed-balance validation for Metropolis accept/reject transitions
//! (tracker P0.4).
//!
//! Two complementary layers:
//!
//! 1. Machine precision — the acceptance statistic reported by
//!    `RandomWalkMetropolis` is exactly `exp(min(0, Δlog p))`, the log-space
//!    Metropolis rule for a symmetric random-walk proposal, and that rule
//!    satisfies the detailed-balance identity `p(x)·A(x,y) == p(y)·A(y,x)`
//!    to within floating-point round-off.
//! 2. Statistical — long non-adaptive chains on a binned continuous target
//!    reproduce symmetric empirical flows `π̂(x)·P̂(x→y) ≈ π̂(y)·P̂(y→x)` for
//!    both `RandomWalkMetropolis` and `ComponentWiseMetropolis`, and the
//!    binned occupancy matches the analytically binned target density.

use std::sync::Mutex;

use mcmc_rs::{
    ComponentWiseMetropolis, EuclideanState, FnLogDensity, RandomWalkMetropolis, SamplingPhase,
    TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ── Machine precision ─────────────────────────────────────────────────────

/// The kernel's reported acceptance probability is the log-Metropolis rule
/// `A(x, y) = exp(min(0, log p(y) - log p(x)))`, bit-for-bit, and uphill
/// proposals are accepted deterministically.
#[test]
fn acceptance_statistic_follows_log_metropolis_formula() {
    let evaluations = Mutex::new(Vec::<f64>::new());
    let mut target = FnLogDensity::new(|position: &[f64]| {
        evaluations.lock().unwrap().push(position[0]);
        -0.5 * position[0] * position[0]
    });
    let mut state = EuclideanState::initialize(&mut target, vec![0.3]).unwrap();
    let mut kernel = RandomWalkMetropolis::isotropic(1, 1.2).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(2024);

    let log_pdf = |value: f64| -0.5 * value * value;
    let (mut uphill, mut rejections) = (0usize, 0usize);
    for _ in 0..4_000 {
        let position_before = state.position()[0];
        let log_density_before = state.log_density();
        let seen_before = evaluations.lock().unwrap().len();
        let report = kernel
            .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();

        // exactly one target evaluation per proposal
        assert_eq!(evaluations.lock().unwrap().len(), seen_before + 1);
        let proposed = *evaluations.lock().unwrap().last().unwrap();
        let delta = log_pdf(proposed) - log_density_before;
        assert_eq!(report.log_acceptance_probability, Some(delta.min(0.0)));
        assert_eq!(report.acceptance_statistic, Some(delta.min(0.0).exp()));
        if delta >= 0.0 {
            uphill += 1;
            // Δ ≥ 0 accepts without drawing a uniform variate
            assert_eq!(report.accepted, Some(true));
            assert_eq!(state.position()[0], proposed);
        } else if report.accepted == Some(false) {
            rejections += 1;
            assert_eq!(state.position()[0], position_before);
        }
    }
    assert!(uphill > 100, "only {uphill} uphill proposals observed");
    assert!(rejections > 100, "only {rejections} rejections observed");
}

/// For a symmetric proposal, detailed balance reduces to
/// `p(x)·A(x,y) == p(y)·A(y,x)` with
/// `A(x,y) = exp(min(0, log p(y) - log p(x)))` — the exact rule applied by the
/// kernels in log space. In log form the identity reads
/// `log p(x) + min(0, Δ) == log p(y) + min(0, -Δ)` with
/// `Δ = log p(y) - log p(x)`, whose exact value is `min(log p(x), log p(y))`.
/// We sweep constructed log-density pairs (including near-ties where Δ sits
/// at the ulp scale) and require agreement to floating-point round-off;
/// `acceptance_statistic_follows_log_metropolis_formula` ties this formula to
/// the kernels' real code path.
#[test]
fn log_metropolis_rule_satisfies_detailed_balance_at_machine_precision() {
    let log_p_values = [-800.0, -300.0, -42.0, -7.5, -1.0, -1e-3, -1e-9, 0.0];
    let delta_magnitudes = [
        0.0,
        f64::EPSILON,
        1e-12,
        1e-8,
        1e-3,
        0.5,
        3.0,
        40.0,
        300.0,
        700.0,
        1e5,
    ];
    for &log_p_x in &log_p_values {
        for &magnitude in &delta_magnitudes {
            for &sign in &[1.0, -1.0] {
                let delta = magnitude * sign;
                let log_p_y = log_p_x + delta;
                let forward = (log_p_y - log_p_x).min(0.0);
                let reverse = (log_p_x - log_p_y).min(0.0);

                // log-space flow: log[p(x) A(x,y)] vs log[p(y) A(y,x)]
                let flow_x = log_p_x + forward;
                let flow_y = log_p_y + reverse;
                let tolerance = 1e-14 * (1.0 + flow_x.abs() + flow_y.abs() + delta.abs());
                assert!(
                    (flow_x - flow_y).abs() <= tolerance,
                    "log flows differ at log p(x)={log_p_x}, Δ={delta}: \
                     {flow_x} vs {flow_y}"
                );
                assert!(
                    (flow_x - log_p_x.min(log_p_y)).abs() <= tolerance,
                    "flow does not reduce to min(log p) at \
                     log p(x)={log_p_x}, Δ={delta}: {flow_x}"
                );

                // linear-space identity where exp neither under- nor overflows
                if log_p_x.abs() < 50.0 && log_p_y.abs() < 50.0 {
                    let weight_x = log_p_x.exp() * forward.exp();
                    let weight_y = log_p_y.exp() * reverse.exp();
                    let reference = weight_x.abs().max(weight_y.abs()).max(f64::MIN_POSITIVE);
                    assert!(
                        (weight_x - weight_y).abs() <= 1e-12 * reference,
                        "linear flows differ at log p(x)={log_p_x}, Δ={delta}: \
                         {weight_x} vs {weight_y}"
                    );
                }
            }
        }
    }
}

// ── Statistical flow balance on a binned continuous target ────────────────

/// Number of equal-width bins covering the target support for flow counting.
const FLOW_BUCKETS: usize = 8;
/// Bin edges span [-FLOW_RANGE, FLOW_RANGE]; mass beyond the edges clamps
/// into the edge bins (negligible for the target below and still reversible).
const FLOW_RANGE: f64 = 3.5;
const BURN_IN: usize = 10_000;
const FLOW_DRAWS: usize = 300_000;

/// Bimodal mixture of N(±1.5, 0.8²): asymmetric within bins, mildly bimodal,
/// yet fast-mixing under unit-scale random-walk proposals.
fn mixture_log_density(x: f64) -> f64 {
    let left = -0.5 * ((x + 1.5) / 0.8).powi(2);
    let right = -0.5 * ((x - 1.5) / 0.8).powi(2);
    let peak = left.max(right);
    peak + ((left - peak).exp() + (right - peak).exp()).ln()
}

fn bucket_of(x: f64) -> usize {
    let width = 2.0 * FLOW_RANGE / FLOW_BUCKETS as f64;
    (((x + FLOW_RANGE) / width).floor()).clamp(0.0, (FLOW_BUCKETS - 1) as f64) as usize
}

/// Simpson quadrature of exp(mixture_log_density) over [lower, upper].
fn integrate_mixture(lower: f64, upper: f64, intervals: usize) -> f64 {
    let step = (upper - lower) / intervals as f64;
    let density = |x: f64| mixture_log_density(x).exp();
    let mut integral = density(lower) + density(upper);
    for index in 1..intervals {
        let weight = if index.is_multiple_of(2) { 2.0 } else { 4.0 };
        integral += weight * density(lower + index as f64 * step);
    }
    integral * step / 3.0
}

/// Binned target probabilities for the mixture, with the edge buckets
/// absorbing the tails (integrated out to ±8 where the mass is negligible).
fn mixture_bucket_probabilities() -> [f64; FLOW_BUCKETS] {
    let width = 2.0 * FLOW_RANGE / FLOW_BUCKETS as f64;
    let edge = |bucket: usize| -FLOW_RANGE + bucket as f64 * width;
    let weights = (0..FLOW_BUCKETS)
        .map(|bucket| {
            let lower = if bucket == 0 { -8.0 } else { edge(bucket) };
            let upper = if bucket == FLOW_BUCKETS - 1 {
                8.0
            } else {
                edge(bucket + 1)
            };
            integrate_mixture(lower, upper, 128)
        })
        .collect::<Vec<_>>();
    let total = weights.iter().sum::<f64>();
    let mut probabilities = [0.0f64; FLOW_BUCKETS];
    for (probability, weight) in probabilities.iter_mut().zip(&weights) {
        *probability = weight / total;
    }
    probabilities
}

/// Run one non-adaptive chain (burn-in discarded) and bin every retained
/// transition, producing an empirical bucket-to-bucket count matrix together
/// with the overall acceptance rate. A time-homogeneous kernel is essential
/// here, so no adaptation is configured.
fn binned_transition_counts<K, F>(
    kernel: &mut K,
    target: &mut FnLogDensity<F>,
    seed: u64,
) -> ([[usize; FLOW_BUCKETS]; FLOW_BUCKETS], f64)
where
    K: TransitionKernel<FnLogDensity<F>>,
    F: FnMut(&[f64]) -> f64 + Send,
{
    let mut state = EuclideanState::initialize(target, vec![0.0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..BURN_IN {
        kernel
            .transition(target, &mut state, &mut rng, SamplingPhase::Warmup)
            .unwrap();
    }

    let mut counts = [[0usize; FLOW_BUCKETS]; FLOW_BUCKETS];
    let mut proposals = 0usize;
    let mut acceptances = 0usize;
    let mut previous = bucket_of(state.position()[0]);
    for _ in 0..FLOW_DRAWS {
        let report = kernel
            .transition(target, &mut state, &mut rng, SamplingPhase::Sampling)
            .unwrap();
        let current = bucket_of(state.position()[0]);
        counts[previous][current] += 1;
        proposals += report.proposals as usize;
        acceptances += report.acceptances as usize;
        previous = current;
    }
    (counts, acceptances as f64 / proposals as f64)
}

fn assert_empirical_flow_balance(counts: &[[usize; FLOW_BUCKETS]; FLOW_BUCKETS], acceptance: f64) {
    // a non-degenerate accept/reject mix is required for the check to bite
    assert!(
        (0.15..0.9).contains(&acceptance),
        "degenerate acceptance rate {acceptance}"
    );
    let total: usize = counts.iter().flatten().sum();
    let row_sums: Vec<usize> = counts.iter().map(|row| row.iter().sum()).collect();
    for (bucket, &row) in row_sums.iter().enumerate() {
        assert!(row > 5_000, "bucket {bucket} saw only {row} visits");
    }

    // Empirical flow balance: π̂(i)·P̂(i→j) = n_ij / total, so detailed
    // balance demands n_ij ≈ n_ji. The count difference sits well inside a
    // Poisson envelope sqrt(n_ij + n_ji) for reversible chains (measured
    // worst case ≈ 2 envelope units over the 28 bin pairs), while a
    // deliberate ×1.5 acceptance asymmetry lands ≈ 12 envelope units out; the
    // 4-unit band below separates both by a wide margin.
    for (source, row) in counts.iter().enumerate() {
        for (destination, &forward) in row.iter().enumerate().skip(source + 1) {
            let reverse = counts[destination][source];
            let asymmetry = (forward as f64 - reverse as f64).abs();
            let bound = 4.0 * ((forward + reverse + 1) as f64).sqrt();
            assert!(
                asymmetry <= bound,
                "flow({source}→{destination}) breaks balance: \
                 {forward} vs {reverse} (bound {bound})"
            );
        }
    }

    // Stationarity under the binned chain: the occupancy π̂ must reproduce
    // the analytically binned target probabilities (Simpson quadrature of
    // the mixture density; edge bins accumulate the clamped tails). Sampling
    // noise keeps the deviation below ≈ 4·10^-3; a tilted acceptance rule
    // shifts bucket mass by ≈ 0.1, far outside the 8·10^-3 band.
    let analytic = mixture_bucket_probabilities();
    let mut worst = 0.0f64;
    for (bucket, &row) in row_sums.iter().enumerate() {
        let empirical = row as f64 / total as f64;
        worst = worst.max((empirical - analytic[bucket]).abs());
    }
    assert!(worst < 8e-3, "π̂ deviates from the binned target by {worst}");
}

#[test]
fn random_walk_flow_is_balanced_on_binned_mixture() {
    let mut kernel = RandomWalkMetropolis::isotropic(1, 1.0).unwrap();
    let mut target = FnLogDensity::new(|position: &[f64]| mixture_log_density(position[0]));
    let (counts, acceptance) = binned_transition_counts(&mut kernel, &mut target, 4242);
    assert_empirical_flow_balance(&counts, acceptance);
}

#[test]
fn component_wise_flow_is_balanced_on_binned_mixture() {
    let mut kernel = ComponentWiseMetropolis::new(vec![1.0]).unwrap();
    let mut target = FnLogDensity::new(|position: &[f64]| mixture_log_density(position[0]));
    let (counts, acceptance) = binned_transition_counts(&mut kernel, &mut target, 8281);
    assert_empirical_flow_balance(&counts, acceptance);
}
