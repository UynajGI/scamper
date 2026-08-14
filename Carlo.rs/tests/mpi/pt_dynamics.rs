//! MPI physics tests for parallel-tempering exchange dynamics.
//!
//! Validates the Metropolis exchange statistics of `PtExchange` against exact
//! analytic results for a two-state toy model with weight W(x | β) = exp(−β·x),
//! so π_β(x = 1) = 1/(1 + e^β):
//!
//! * size == 2 — the measured exchange acceptance rate matches the enumerated
//!   min(1, e^R) rate for β = [0.3, 0.7] (R = (β_a−β_b)(x_a−x_b), ≈ 0.9063).
//! * size >= 4 — replica-label ergodicity: every rank hosts every chain label
//!   (round trips need both even and odd pairings), and the final labels form
//!   a permutation of 0..size with `current_value() == parameter_values[label]`.
//! * size >= 2 — a rank/value count mismatch is rejected as a clean
//!   `InvalidConfig` error on every rank (`PtExchange::new` performs no
//!   collectives, so no rank can deadlock).
//! * size == 1 — a single-chain PT run is degenerate but valid: `try_exchange`
//!   is always `Ok(false)` and `try_step` completes.
//!
//! MPI can be initialized only once per process, so this file contains exactly
//! one `#[test]` that branches on the world size. Run it with:
//!
//! ```bash
//! mpirun -np 2 cargo test --features mpi --test suite -- --ignored --exact mpi_pt_dynamics::pt_exchange_dynamics_suite --nocapture
//! mpirun -np 4 cargo test --features mpi --test suite -- --ignored --exact mpi_pt_dynamics::pt_exchange_dynamics_suite --nocapture
//! ```

#![cfg(feature = "mpi")]

use carlo_rs::parallel_tempering::{ParallelTemperingConfig, PtExchange};
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};
use mpi::environment::Universe;
use mpi::traits::*;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

// ── Toy model ─────────────────────────────────────────────────────────────

/// Two-state model: x ∈ {0, 1} with weight W(x | β) = exp(−β·x).
/// Each sweep redraws x exactly from π_β, so the state is independent of the
/// exchange history and single-attempt acceptance rates are analytically
/// solvable.
struct TwoStateMc {
    beta: f64,
    x: u8,
}

/// Stationary probability of x = 1 at inverse temperature β.
fn pi_one(beta: f64) -> f64 {
    1.0 / (1.0 + beta.exp())
}

impl MonteCarlo for TwoStateMc {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Draw u ∈ [0, 1) once per sweep: x = 1 iff u < π_β(1).
        let u: f64 = ctx.rng.random();
        self.x = u8::from(u < pi_one(self.beta));
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("X", f64::from(self.x));
    }

    fn name(&self) -> &'static str {
        "TwoStateMc"
    }
}

impl FromParams for TwoStateMc {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let beta = params
            .get::<f64>("beta")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "beta".into(),
                reason: "TwoStateMc requires a floating-point beta parameter".into(),
            })?;
        // Draw the initial state from the stationary distribution.
        let u: f64 = rng.random();
        let x = u8::from(u < pi_one(beta));
        Ok(Self { beta, x })
    }
}

impl ParallelTemperingCompatible for TwoStateMc {
    /// log[W(x|new)/W(x|cur)] = (β_cur − β_new)·x.
    fn log_weight_ratio(&self, _param: &str, new_value: f64) -> f64 {
        (self.beta - new_value) * f64::from(self.x)
    }

    /// The parameter label travels; the configuration x stays on this rank.
    fn change_parameter(&mut self, _param: &str, new_value: f64) {
        self.beta = new_value;
    }
}

// ── Analytic reference ────────────────────────────────────────────────────

/// Exact single-attempt acceptance probability for exchanging the parameters
/// of two chains at (β_a, β_b) whose states are drawn independently from their
/// stationary distributions.
///
/// Swapping the two parameters changes the total weight by
/// e^R with R = (β_a − β_b)(x_a − x_b), and the Metropolis rule accepts with
/// probability min(1, e^R). Enumerating the four (x_a, x_b) pairs:
/// p_accept = Σ π_a(x_a)·π_b(x_b)·min(1, e^R).
fn analytic_acceptance(beta_a: f64, beta_b: f64) -> f64 {
    let pa = pi_one(beta_a);
    let pb = pi_one(beta_b);
    let mut total = 0.0;
    for (x_a, w_a) in [(0i32, 1.0 - pa), (1, pa)] {
        for (x_b, w_b) in [(0i32, 1.0 - pb), (1, pb)] {
            let log_r = (beta_a - beta_b) * f64::from(x_a - x_b);
            total += w_a * w_b * log_r.exp().min(1.0);
        }
    }
    total
}

/// Mirror of the even-odd pairing in `PtExchange::try_exchange`: at pairing
/// offset `pairing_offset`, chains whose index matches the offset pair
/// upward, the others pair downward; a boundary chain left without a partner
/// sits that exchange round out. For 2 chains this means the odd-offset
/// rounds attempt no exchange at all.
fn has_partner(chain_idx: usize, n_chains: usize, pairing_offset: u64) -> bool {
    if chain_idx % 2 == pairing_offset as usize {
        chain_idx + 1 < n_chains
    } else {
        chain_idx > 0
    }
}

// ── The suite (single #[test]: MPI initializes once per process) ──────────

#[test]
#[ignore = "requires mpirun"]
fn pt_exchange_dynamics_suite() {
    let universe = mpi::initialize().expect("test must run under mpirun/mpiexec");
    let size = universe.world().size();
    let rank = universe.world().rank();

    if size == 2 {
        scenario_acceptance_rate_matches_analytic(&universe);
    } else {
        eprintln!(
            "[rank {rank}] skip acceptance-rate scenario (needs exactly 2 ranks, got {size})"
        );
    }

    if size >= 4 {
        scenario_replica_round_trip(&universe);
    } else {
        eprintln!("[rank {rank}] skip round-trip scenario (needs >= 4 ranks, got {size})");
    }

    if size >= 2 {
        scenario_rank_value_mismatch_is_clean_error(&universe);
    } else {
        eprintln!("[rank {rank}] skip mismatch-error scenario (needs >= 2 ranks, got {size})");
    }

    if size == 1 {
        scenario_single_rank_degenerate(&universe);
    } else {
        eprintln!("[rank {rank}] skip single-rank scenario (needs exactly 1 rank, got {size})");
    }
}

// ── Scenario: size == 2, acceptance rate vs analytic ──────────────────────

/// `try_step` triggers `try_exchange` whenever sweep_count > 0 and
/// sweep_count % interval == 0; with interval = 1 that is every step. An
/// exchange is only actually attempted when the even-odd pairing leaves this
/// chain with a partner (for 2 chains, the odd-offset rounds pair nothing),
/// and an accepted exchange flips this rank's chain label, i.e. changes
/// `current_value()`. States are redrawn from π_β every sweep, so real
/// attempts are i.i.d. with the enumerated acceptance rate.
fn scenario_acceptance_rate_matches_analytic(universe: &Universe) {
    let world = universe.world();
    let rank = world.rank();
    let n_chains = world.size() as usize;
    eprintln!("[rank {rank}] scenario: exchange acceptance rate vs analytic (beta = [0.3, 0.7])");

    let beta_a = 0.3;
    let beta_b = 7.0 / 10.0;
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![beta_a, beta_b],
        interval: 1,
    };
    let mut params = Params::new();
    params.set("beta", "0.3");

    let mut exchange = PtExchange::<TwoStateMc, Xoshiro256PlusPlus>::new(
        universe.world(),
        &config,
        &params,
        2024,    // seed
        1,       // binsize
        100_000, // target sweeps == exchange attempts
    )
    .expect("2 ranks must match 2 PT values");

    let expected = analytic_acceptance(beta_a, beta_b);
    // Cross-check the enumeration against the known closed-form value.
    assert!(
        (expected - 0.9063).abs() < 5.0e-4,
        "enumerated p_accept = {expected}, expected ≈ 0.9063"
    );

    let mut steps = 0u64;
    let mut attempts = 0u64;
    let mut accepted = 0u64;
    while !exchange.is_complete() {
        let value_before = exchange.current_value();
        // Labels only move via exchanges, so the label during this step's
        // exchange decision is the one held before the step.
        let chain_before = exchange.chain_idx();
        exchange.try_step().expect("PT step must succeed");
        steps += 1;
        let exchange_step = exchange.context().sweep_count() / config.interval;
        let pairing_offset = exchange_step.saturating_sub(1) & 1;
        if has_partner(chain_before, n_chains, pairing_offset) {
            attempts += 1;
            if exchange.current_value() != value_before {
                accepted += 1;
            }
        }
    }

    // With 2 chains the odd-offset exchange rounds pair nothing, so exactly
    // every other step attempts an exchange.
    assert_eq!(
        attempts,
        steps / 2,
        "2 chains must attempt exchanges on exactly the even-offset rounds"
    );
    assert!(
        attempts >= 50_000,
        "expected >= 50k attempts, got {attempts}"
    );
    let n = attempts as f64;
    let empirical = accepted as f64 / n;
    let five_sigma = 5.0 * (expected * (1.0 - expected) / n).sqrt();
    assert!(
        (empirical - expected).abs() <= five_sigma,
        "empirical acceptance rate {empirical:.5} must match analytic {expected:.5} \
         within 5σ = {five_sigma:.5} ({accepted}/{attempts} accepted)"
    );

    // Both ranks observe the same broadcast exchange decisions, so their
    // counts must agree exactly.
    let mut counts = [0u64; 2];
    world.all_gather_into(&accepted, counts.as_mut_slice());
    assert_eq!(
        counts[0], counts[1],
        "every rank must count the same accepted exchanges: {counts:?}"
    );

    eprintln!(
        "[rank {rank}] acceptance rate: empirical = {empirical:.5}, analytic = {expected:.5}, \
         {attempts} attempts in {steps} steps"
    );
}

// ── Scenario: size >= 4, replica round-trip / label ergodicity ────────────

/// Chains at β = 0.0, 0.4, 0.8, 1.2, ... (spacing 0.4). Alternating even/odd
/// pairings are required for labels to random-walk across the whole ladder;
/// over 30k steps every rank must have hosted every label, and the final
/// labels must form a permutation with consistent parameter values.
fn scenario_replica_round_trip(universe: &Universe) {
    let world = universe.world();
    let rank = world.rank();
    let n_chains = world.size() as usize;
    eprintln!("[rank {rank}] scenario: replica round-trip / label ergodicity ({n_chains} chains)");

    let values: Vec<f64> = (0..n_chains).map(|i| 0.4 * i as f64).collect();
    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: values.clone(),
        interval: 1,
    };
    let mut params = Params::new();
    params.set("beta", "0.0");

    let mut exchange = PtExchange::<TwoStateMc, Xoshiro256PlusPlus>::new(
        universe.world(),
        &config,
        &params,
        2025,   // seed
        1,      // binsize
        30_000, // target sweeps
    )
    .expect("world size must match the number of PT values");

    let mut hosted = vec![false; n_chains];
    hosted[exchange.chain_idx()] = true;
    let mut steps = 0u64;
    while !exchange.is_complete() {
        exchange.try_step().expect("PT step must succeed");
        hosted[exchange.chain_idx()] = true;
        steps += 1;
    }

    let missing: Vec<usize> = hosted
        .iter()
        .enumerate()
        .filter_map(|(label, &seen)| (!seen).then_some(label))
        .collect();
    assert!(
        missing.is_empty(),
        "rank {rank} never hosted chain labels {missing:?} in {steps} steps — \
         labels are not ergodic across the ladder"
    );

    // Final labels across ranks must form a permutation of 0..n_chains, and
    // each rank must sit at the parameter value of its current label.
    let local_label = exchange.chain_idx() as u64;
    let mut labels = vec![0u64; n_chains];
    world.all_gather_into(&local_label, labels.as_mut_slice());
    let mut sorted = labels.clone();
    sorted.sort_unstable();
    let expected_labels: Vec<u64> = (0..n_chains).map(|i| i as u64).collect();
    assert_eq!(
        sorted, expected_labels,
        "final chain labels must be a permutation, got {labels:?}"
    );
    assert!(
        (exchange.current_value() - values[exchange.chain_idx()]).abs() < 1.0e-12,
        "current_value() must equal parameter_values[chain_idx]"
    );

    eprintln!(
        "[rank {rank}] hosted all {n_chains} labels after {steps} steps, final label = {local_label}"
    );
}

// ── Scenario: size >= 2, rank/value count mismatch ────────────────────────

/// `PtExchange::new` validates the rank/value count before any communication,
/// so a 3-value config on a larger world must fail identically on every rank
/// and the barrier afterwards proves nobody deadlocked.
fn scenario_rank_value_mismatch_is_clean_error(universe: &Universe) {
    let world = universe.world();
    let rank = world.rank();
    eprintln!("[rank {rank}] scenario: rank/value count mismatch is a clean error");

    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.3, 0.5, 0.7], // 3 values: mismatches any world size != 3
        interval: 1,
    };
    let mut params = Params::new();
    params.set("beta", "0.3");

    match PtExchange::<TwoStateMc, Xoshiro256PlusPlus>::new(
        universe.world(),
        &config,
        &params,
        2024,
        1,
        10,
    ) {
        Err(CarloError::InvalidConfig { field, reason }) => {
            assert_eq!(
                field, "pt_chains",
                "mismatch must be reported against pt_chains, got reason: {reason}"
            );
        }
        Err(other) => panic!("rank {rank}: expected InvalidConfig, got {other:?}"),
        Ok(_) => panic!(
            "rank {rank}: PtExchange::new must reject {} values on a {}-rank world",
            config.values.len(),
            world.size()
        ),
    }

    // Every rank reached the error handling identically: the barrier returns.
    world.barrier();
    eprintln!("[rank {rank}] mismatch rejected with InvalidConfig on every rank, no deadlock");
}

// ── Scenario: size == 1, degenerate single chain ──────────────────────────

fn scenario_single_rank_degenerate(universe: &Universe) {
    let world = universe.world();
    eprintln!("[rank 0] scenario: single-rank PT is degenerate but valid");

    let config = ParallelTemperingConfig {
        parameter: "beta".to_string(),
        values: vec![0.5],
        interval: 1,
    };
    let mut params = Params::new();
    params.set("beta", "0.5");

    let mut exchange =
        PtExchange::<TwoStateMc, Xoshiro256PlusPlus>::new(world, &config, &params, 2024, 1, 5)
            .expect("1 rank with 1 value must construct");

    // A lone chain has no exchange partner.
    let accepted = exchange
        .try_exchange()
        .expect("try_exchange must succeed on a single rank");
    assert!(!accepted, "a single chain can never accept an exchange");

    exchange.try_step().expect("try_step must succeed");
    assert_eq!(exchange.sweep_count(), 1);

    eprintln!("[rank 0] try_exchange returned Ok(false) and try_step completed");
}
