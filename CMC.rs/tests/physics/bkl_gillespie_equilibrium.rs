//! Long-time equilibrium validation for the continuous-time dynamics
//! kernels (`BklIsingKernel` / n-fold-way and the generic `GillespieKernel`).
//!
//! Closes two VALIDATION.md gaps:
//!
//! 1. **BKL long-time equilibrium distribution.**  The pre-existing coverage
//!    (integration/dynamics_stage6.rs) is a *single* β, N=4, ⟨E⟩-only
//!    fixed-time check with a loose absolute tolerance; the exact-trajectory
//!    reproducibility test says nothing about the sampled measure.  Here the
//!    N=4 and N=8 PBC chains are run to long times at three β values
//!    (disordered → ordered-ish) and both ⟨E⟩ and ⟨m²⟩ are compared against
//!    exact enumeration with multi-seed z-scores.
//!
//! 2. **Gillespie multi-state equilibrium distribution.**  Pre-existing
//!    coverage checks rate *selection* and the mean exponential wait.  Here a
//!    fully connected asymmetric 3-state CTMC is driven through the real
//!    `GillespieKernel` + `RejectionFreeModel` machinery and its time-weighted
//!    occupancies are compared against the exact stationary distribution
//!    πQ = 0 (solved in-code by Cramer's rule).  A second test drives the
//!    actual `KineticIsingModel` through `GillespieKernel` and pins ⟨E⟩
//!    against exact enumeration, tying the generic kernel to physics.
//!
//! ## Estimator: residence-time (per-time), not per-visit
//!
//! For a continuous-time Markov chain the canonical equilibrium average of an
//! observable O is the *time* average
//! ```text
//! ⟨O⟩_time = Σ_visits O(state) · Δt_visit / Σ_visits Δt_visit,
//! ```
//! where `Δt_visit` is the sojourn (residence) time of the configuration
//! before the next event.  Per-visit (event-count) averaging is a *different,
//! length-biased* estimator: states are visited with frequency ∝ π_i·q_i but
//! weighted by residence time ∝ 1/q_i, and the two weightings only agree when
//! the total exit rate q_i is state-independent.  Both kernels expose the
//! pre-event sojourn directly (`BklEvent::delta_time` / `GillespieEvent::
//! delta_time`, advanced before the commit), so the physically correct
//! residence-time estimator is used throughout.  The pre-existing fixed-time
//! BKL test samples snapshots at equal clock intervals — a discrete
//! approximation of the same time average — so this test pins the exact
//! continuous-time estimator, more β values, a second observable (⟨m²⟩), and
//! a second system size.
//!
//! ## z-score conventions
//!
//! Same framework as zscore_extended.rs: per-seed estimates with binning-based
//! standard errors, z_i = (mean_i − exact)/stderr_i, asserting max|z| < 4,
//! |z̄| < 2, and Σz > −2√n (one-sided bias floor; −2√8 ≈ −5.66 at the default
//! n = 8).  Bins are contiguous blocks of events; each bin yields the
//! residence-time-weighted mean over its own time span, so the bin-to-bin
//! scatter empirically captures the trajectory autocorrelation (a naive
//! binomial error √(π(1−π)/n_events) would be anticonservative here because
//! consecutive occupancies are correlated).  `SCUTTLE_ZSCORE_SEEDS` raises the
//! seed count for nightly monitoring, as elsewhere.

use super::common::{assert_close, exact_ising_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, BklIsingKernel, DynamicsError, GillespieKernel, Hamiltonian, IsingModel,
    KineticIsingModel, KineticRateLaw, RejectionFreeModel, System,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const N_SEEDS: usize = 8;

/// Cache-consistency audit cadence inside the long BKL runs (events).
const BKL_AUDIT_INTERVAL: u64 = 10_000;

/// Discarded equilibration events per BKL trajectory.  The N ≤ 8 ring
/// decorrelates within O(10²) events even at β = 1, so 20 000 is ≥ 100 τ.
const BKL_EQUILIBRATION_EVENTS: u64 = 20_000;

/// Measurement events per BKL trajectory, split into 5 000-event bins (each
/// bin spans ≥ 10 autocorrelation times, so the bin scatter is a valid
/// stderr estimate).  Sized to keep the default suite fast; raise
/// SCUTTLE_ZSCORE_SEEDS for higher-power nightly monitoring.
const BKL_MEASUREMENT_EVENTS: u64 = 120_000;
const BKL_BINS: usize = 24;

/// Toy-chain and Ising-via-Gillespie trajectory sizes.
const GILLESPIE_EQUILIBRATION_EVENTS: u64 = 10_000;
const GILLESPIE_MEASUREMENT_EVENTS: u64 = 160_000;
const GILLESPIE_BINS: usize = 32;

/// Splitmix64-style seed mixing: sequential seed indices must not produce
/// correlated Xoshiro streams.
fn mixed_seed(base: u64, index: u64) -> u64 {
    let mut z = base.wrapping_add(index).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z ^ (z >> 31)
}

/// Residence-time mean and binning-based standard error.
///
/// `sums[b]` = Σ O·Δt over bin b, `durations[b]` = Σ Δt over bin b.  The mean
/// is the pooled time average; the stderr is the bin-mean scatter / √n_bins.
fn time_weighted_stats(sums: &[f64], durations: &[f64]) -> (f64, f64) {
    assert_eq!(sums.len(), durations.len());
    let total: f64 = durations.iter().sum();
    let mean = sums.iter().sum::<f64>() / total;
    let n_bins = sums.len() as f64;
    let variance = sums
        .iter()
        .zip(durations)
        .map(|(sum, duration)| sum / duration - mean)
        .map(|deviation| deviation * deviation)
        .sum::<f64>()
        / (n_bins - 1.0);
    (mean, variance.sqrt() / n_bins.sqrt())
}

/// Assert the repo z-score criteria (zscore_extended.rs conventions).
fn assert_exact_value_z_scores(estimates: &[(f64, f64)], exact: f64, label: &str) {
    let n = estimates.len() as f64;
    let z_scores: Vec<f64> = estimates
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-12))
        .collect();
    let max_abs_z = z_scores.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let sum_z: f64 = z_scores.iter().sum();
    let mean_z = sum_z / n;
    let sum_z_floor = -2.0 * n.sqrt();
    assert!(
        max_abs_z < 4.0,
        "{label}: max |z| = {max_abs_z:.2} should be < 4 (exact = {exact:.6})"
    );
    assert!(
        mean_z.abs() < 2.0,
        "{label}: mean z = {mean_z:.2} should be |z̄| < 2"
    );
    assert!(
        sum_z > sum_z_floor,
        "{label}: sum z = {sum_z:.2} should be > {sum_z_floor:.2} (no one-sided bias)"
    );
}

// ════════════════════════════════════════════════════════════════════════
// A) BKL / n-fold-way: long-time equilibrium of the kinetic Ising chain
// ════════════════════════════════════════════════════════════════════════

/// Run one BKL trajectory; returns (⟨E⟩, stderr, ⟨m²⟩, stderr) using the
/// residence-time estimator.
fn bkl_time_weighted_moments(n_sites: usize, beta: f64, seed: u64) -> (f64, f64, f64, f64) {
    let lattice = build_chain(n_sites, true);
    let state = System::new(lattice, 1, 1.0, beta);
    let model = KineticIsingModel::new(1.0, KineticRateLaw::glauber(1.0).unwrap()).unwrap();
    let mut kernel = BklIsingKernel::new(model, state, BKL_AUDIT_INTERVAL).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    for _ in 0..BKL_EQUILIBRATION_EVENTS {
        kernel.step(&mut rng).unwrap().unwrap();
    }

    let events_per_bin = BKL_MEASUREMENT_EVENTS / BKL_BINS as u64;
    let mut energy_sums = vec![0.0; BKL_BINS];
    let mut m2_sums = vec![0.0; BKL_BINS];
    let mut durations = vec![0.0; BKL_BINS];
    for ((energy_sum, m2_sum), duration) in
        energy_sums.iter_mut().zip(&mut m2_sums).zip(&mut durations)
    {
        let mut previous_time = kernel.event_time();
        for _ in 0..events_per_bin {
            // Observe the pre-event configuration, then weight it by the
            // sojourn time `event.delta_time` the kernel reports for it.
            let energy = kernel.state().energy;
            let magnetization = kernel.state().spins.iter().sum::<f64>() / n_sites as f64;
            let event = kernel.step(&mut rng).unwrap().unwrap();
            assert!(event.delta_time > 0.0, "BKL sojourn times must be positive");
            let now = kernel.event_time();
            assert!(
                now > previous_time,
                "BKL event clock must strictly increase"
            );
            previous_time = now;
            *energy_sum += energy * event.delta_time;
            *m2_sum += magnetization * magnetization * event.delta_time;
            *duration += event.delta_time;
        }
    }
    kernel.validate().unwrap();
    // The incremental energy bookkeeping must still match a direct
    // recomputation after the full trajectory.
    let direct = IsingModel::new(1.0).compute_total_energy(
        &kernel.state().spins,
        &kernel.state().lattice,
        1.0,
    );
    assert_close(kernel.state().energy, direct, 1e-9);
    assert!(kernel.events() >= BKL_EQUILIBRATION_EVENTS + BKL_MEASUREMENT_EVENTS);

    let (energy_mean, energy_stderr) = time_weighted_stats(&energy_sums, &durations);
    let (m2_mean, m2_stderr) = time_weighted_stats(&m2_sums, &durations);
    (energy_mean, energy_stderr, m2_mean, m2_stderr)
}

fn check_bkl_equilibrium(n_sites: usize, betas: &[f64], seed_base: u64, label: &str) {
    let lattice = build_chain(n_sites, true);
    let n_seeds = zscore_seed_count(N_SEEDS);
    for &beta in betas {
        let (_, exact_energy, _, exact_m2) = exact_ising_moments(&lattice, 1.0, beta);
        let mut energy_estimates = Vec::with_capacity(n_seeds);
        let mut m2_estimates = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let (e, e_se, m2, m2_se) =
                bkl_time_weighted_moments(n_sites, beta, mixed_seed(seed_base, seed));
            energy_estimates.push((e, e_se));
            m2_estimates.push((m2, m2_se));
        }
        assert_exact_value_z_scores(
            &energy_estimates,
            exact_energy,
            &format!("{label} N={n_sites} β={beta} ⟨E⟩"),
        );
        assert_exact_value_z_scores(
            &m2_estimates,
            exact_m2,
            &format!("{label} N={n_sites} β={beta} ⟨m²⟩"),
        );
    }
}

#[test]
fn bkl_n4_time_weighted_equilibrium_matches_exact_ising_moments() {
    check_bkl_equilibrium(4, &[0.2, 0.6, 1.0], 0xB41_1201, "BKL");
}

#[test]
fn bkl_n8_time_weighted_equilibrium_matches_exact_ising_moments() {
    check_bkl_equilibrium(8, &[0.2, 0.6, 1.0], 0xB41_1202, "BKL");
}

// ════════════════════════════════════════════════════════════════════════
// B) Gillespie: multi-state equilibrium of an asymmetric toy CTMC
// ════════════════════════════════════════════════════════════════════════

/// Fully connected asymmetric 3-state chain.  Event `j` in state `i` is the
/// jump i → j with rate `rates[i][j]` (diagonal entries are unused/zero, so
/// the catalog contains every admissible transition and nothing else).
#[derive(Clone)]
struct AsymmetricChain {
    rates: [[f64; 3]; 3],
}

impl AsymmetricChain {
    fn new() -> Self {
        Self {
            rates: [[0.0, 1.0, 0.9], [0.5, 0.0, 0.7], [0.4, 1.3, 0.0]],
        }
    }
}

impl RejectionFreeModel for AsymmetricChain {
    type State = usize;
    type Patch = ();

    fn event_count(&self, _state: &Self::State) -> usize {
        3
    }

    fn event_rate(&self, state: &Self::State, event: usize) -> Result<f64, DynamicsError> {
        Ok(self.rates[*state][event])
    }

    fn prepare_event(
        &self,
        _state: &Self::State,
        _event: usize,
        _patch: &mut Self::Patch,
    ) -> Result<(), DynamicsError> {
        Ok(())
    }

    fn commit_event(&self, state: &mut Self::State, event: usize, _patch: &Self::Patch) {
        *state = event;
    }

    fn validate_state(&self, state: &Self::State) -> Result<(), DynamicsError> {
        if *state < 3 {
            Ok(())
        } else {
            Err(DynamicsError::new("asymmetric chain state out of range"))
        }
    }
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Exact stationary distribution π of the chain: solve πQ = 0 with Σπ = 1 by
/// Cramer's rule on Qᵀ with the third equation replaced by the normalization
/// (the bordered system is nonsingular for an irreducible chain).
fn exact_stationary_distribution(chain: &AsymmetricChain) -> [f64; 3] {
    let mut q = [[0.0_f64; 3]; 3];
    for (i, row) in q.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            if i != j {
                *entry = chain.rates[i][j];
            }
        }
        row[i] = -chain.rates[i].iter().sum::<f64>();
    }
    // Coefficient matrix of πQ = 0 is Qᵀ; replace the last row by ones.
    let mut bordered = [[0.0_f64; 3]; 3];
    for target in 0..3 {
        for source in 0..3 {
            bordered[target][source] = q[source][target];
        }
    }
    bordered[2] = [1.0; 3];
    let determinant = det3(&bordered);
    assert!(
        determinant.abs() > 1e-9,
        "irreducible chain must be solvable"
    );
    let mut pi = [0.0_f64; 3];
    for component in 0..3 {
        let mut replaced = bordered;
        for (row, entry) in replaced.iter_mut().enumerate() {
            entry[component] = if row == 2 { 1.0 } else { 0.0 };
        }
        pi[component] = det3(&replaced) / determinant;
    }
    // In-code sanity: π solves the stationary equations and is normalized.
    let mut residual = [0.0_f64; 3];
    for (source, q_row) in q.iter().enumerate() {
        for (target, &rate) in q_row.iter().enumerate() {
            residual[target] += pi[source] * rate;
        }
    }
    assert!(
        residual.iter().all(|value| value.abs() < 1e-12),
        "stationary residuals {residual:?} exceed 1e-12"
    );
    assert_close(pi.iter().sum::<f64>(), 1.0, 1e-12);
    pi
}

/// Time-weighted occupancy fractions of one Gillespie trajectory.
/// Returns per-state (fraction, stderr) plus sanity data
/// (event count, final event time, Σ_state occupied time).
fn gillespie_toy_occupancies(seed: u64) -> (Vec<(f64, f64)>, u64, f64, f64) {
    let mut kernel = GillespieKernel::new(AsymmetricChain::new(), 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..GILLESPIE_EQUILIBRATION_EVENTS {
        kernel.step(&mut rng).unwrap().unwrap();
    }

    let events_per_bin = GILLESPIE_MEASUREMENT_EVENTS / GILLESPIE_BINS as u64;
    let mut occupancy_sums = vec![[0.0_f64; 3]; GILLESPIE_BINS];
    let mut durations = vec![0.0_f64; GILLESPIE_BINS];
    let measurement_start = kernel.event_time();
    let mut previous_time = kernel.event_time();
    for (occupancy_sum, duration) in occupancy_sums.iter_mut().zip(&mut durations) {
        for _ in 0..events_per_bin {
            let state = *kernel.state();
            let event = kernel.step(&mut rng).unwrap().unwrap();
            assert!(event.delta_time > 0.0, "Gillespie waits must be positive");
            let now = kernel.event_time();
            assert!(
                now > previous_time,
                "Gillespie event clock must strictly increase"
            );
            previous_time = now;
            occupancy_sum[state] += event.delta_time;
            *duration += event.delta_time;
        }
    }
    kernel.validate().unwrap();

    let total_occupied: f64 = occupancy_sums.iter().flatten().sum();
    let total_time: f64 = durations.iter().sum();
    // Occupancy bookkeeping must exactly account for the measurement clock:
    // Σ_state occupied time = Σ sojourns = event-time span of the window.
    assert_close(total_occupied, total_time, 1e-9 * total_time);
    assert_close(
        kernel.event_time() - measurement_start,
        total_time,
        1e-9 * total_time,
    );

    let fractions = (0..3)
        .map(|state| {
            let sums: Vec<f64> = occupancy_sums.iter().map(|bin| bin[state]).collect();
            time_weighted_stats(&sums, &durations)
        })
        .collect();
    (
        fractions,
        kernel.events(),
        kernel.event_time(),
        total_occupied,
    )
}

#[test]
fn gillespie_asymmetric_three_state_occupancy_matches_exact_stationary_distribution() {
    let pi_exact = exact_stationary_distribution(&AsymmetricChain::new());
    let n_seeds = zscore_seed_count(N_SEEDS);
    let mut per_state: Vec<Vec<(f64, f64)>> = vec![Vec::new(); 3];
    for seed in 0..n_seeds as u64 {
        let (fractions, events, event_time, total_occupied) =
            gillespie_toy_occupancies(mixed_seed(0x61_5EED1, seed));
        assert!(events >= GILLESPIE_EQUILIBRATION_EVENTS + GILLESPIE_MEASUREMENT_EVENTS);
        assert!(event_time > 0.0);
        assert!(total_occupied > 0.0);
        for (state, estimate) in fractions.iter().enumerate() {
            per_state[state].push(*estimate);
        }
    }
    for (state, estimates) in per_state.iter().enumerate() {
        assert_exact_value_z_scores(
            estimates,
            pi_exact[state],
            &format!("Gillespie 3-state occupancy[{state}]"),
        );
    }
}

// ════════════════════════════════════════════════════════════════════════
// B') The same generic kernel on physical Ising rates
// ════════════════════════════════════════════════════════════════════════

/// Residence-time ⟨E⟩ of one `GillespieKernel<KineticIsingModel>` trajectory
/// on the N=4 PBC chain (rates re-derived from the actual Glauber rate law
/// before every event, unlike the Fenwick-cached BKL kernel).
fn gillespie_ising_energy(beta: f64, seed: u64) -> (f64, f64) {
    let n_sites = 4;
    let lattice = build_chain(n_sites, true);
    let mut state = System::new(lattice, 1, 1.0, beta);
    let model = KineticIsingModel::new(1.0, KineticRateLaw::glauber(1.0).unwrap()).unwrap();
    // Unlike BklIsingKernel::new, the generic GillespieKernel validates the
    // energy cache at construction, so it must be primed explicitly.
    state.recompute_energy(model.hamiltonian());
    let mut kernel = GillespieKernel::new(model, state).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    for _ in 0..GILLESPIE_EQUILIBRATION_EVENTS {
        kernel.step(&mut rng).unwrap().unwrap();
    }

    let events_per_bin = GILLESPIE_MEASUREMENT_EVENTS / GILLESPIE_BINS as u64;
    let mut energy_sums = vec![0.0_f64; GILLESPIE_BINS];
    let mut durations = vec![0.0_f64; GILLESPIE_BINS];
    for (energy_sum, duration) in energy_sums.iter_mut().zip(&mut durations) {
        for _ in 0..events_per_bin {
            let energy = kernel.state().energy;
            let event = kernel.step(&mut rng).unwrap().unwrap();
            assert!(event.delta_time > 0.0);
            *energy_sum += energy * event.delta_time;
            *duration += event.delta_time;
        }
    }
    kernel.validate().unwrap();
    assert!(kernel.events() >= GILLESPIE_EQUILIBRATION_EVENTS + GILLESPIE_MEASUREMENT_EVENTS);
    time_weighted_stats(&energy_sums, &durations)
}

#[test]
fn gillespie_kinetic_ising_time_weighted_energy_matches_exact() {
    let beta = 0.6;
    let lattice = build_chain(4, true);
    let (_, exact_energy, _, _) = exact_ising_moments(&lattice, 1.0, beta);
    let n_seeds = zscore_seed_count(N_SEEDS);
    let estimates: Vec<(f64, f64)> = (0..n_seeds as u64)
        .map(|seed| gillespie_ising_energy(beta, mixed_seed(0x61_5EED2, seed)))
        .collect();
    assert_exact_value_z_scores(&estimates, exact_energy, "Gillespie-Ising N=4 β=0.6 ⟨E⟩");
}
