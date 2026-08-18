//! Wang-Landau production on the continuous `BinnedAxis` energy axis.
//!
//! The already-validated route is `DiscreteAxis` with one bin per exact level
//! (see `integration/generalized_stage4.rs` and `physics/generalized_exact.rs`).
//! This file closes the remaining gap: a physical production run whose energy
//! axis is a `BinnedAxis`, i.e. bins of finite width that may integrate
//! several distinct microstate levels and reweight through bin centers.
//!
//! Two systems pin the two halves of that route:
//!
//! 1. **Weighted 6-ring, 14-bin unit-width axis.** Alternating bond weights
//!    `1.0/1.1` split the spectrum into eight levels
//!    `E ∈ {-6.3, -2.3, -2.1, -1.9, +1.9, +2.1, +2.3, +6.3}`; bins 4 and 9
//!    each integrate *two* levels (binned degeneracies {2, 24, 6, 6, 24, 2}),
//!    so per-bin agreement is a genuine test of bin integration, not a
//!    one-level-per-bin restatement of the discrete case.
//! 2. **Uniform 8-ring, 5-bin axis spanning [-10, 10]** whose bin centers
//!    coincide exactly with the levels `{-8, -4, 0, +4, +8}` (degeneracies
//!    `{2, 56, 140, 56, 2}`).  Reweighting through the centers is then
//!    geometry-exact, isolating the Wang-Landau statistical error when
//!    compared against `exact_ising_moments`.
//!
//! Both refinement routes of `WangLandauCore` are exercised through the
//! binned axis: geometric flat-histogram halving, and `1/t` refinement once
//! `ln f` drops below a threshold (which removes the flat-histogram error
//! saturation floor).
//!
//! Setting `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count for nightly
//! high-power monitoring (unset → the documented per-test default).

use super::common::{exact_ising_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, canonical_reweight, enumerate_ising_density_of_states, Algorithm, BinnedAxis,
    Bond, BondType, CsrLattice, IsingModel, MacrostateAxis, SimulationPhase, System,
    WangLandauConfig, WangLandauCore, WangLandauState, WangLandauTermination,
};
use rand::SeedableRng;

const DEFAULT_SEEDS: usize = 8;

/// Weighted periodic 6-ring with alternating bond weights.
fn weighted_ring() -> CsrLattice {
    let edges: Vec<Bond> = (0..6)
        .map(|site| {
            let weight = if site % 2 == 0 { 1.0 } else { 1.1 };
            let target = (site + 1) % 6;
            Bond::new(site, target, BondType::Generic, weight)
        })
        .collect();
    CsrLattice::from_edges(6, edges)
}

/// 1/t-refined production config: flat-histogram halving down to `ln f` =
/// 1/256, then `1/t` refinement to `ln f` = 1/65536.  Switching to 1/t only
/// after eight halvings keeps the error carried out of the flat-histogram
/// phase small; the 1/t phase then contracts it further as sqrt(t0/t).
fn one_over_t_config(minimum_visited_fraction: f64) -> WangLandauConfig {
    WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 1.0 / 65536.0,
        flatness: 0.8,
        flatness_check_interval: 25,
        discovery_sweeps: 0,
        one_over_t_threshold: 1.0 / 256.0,
        max_adaptation_sweeps: 2_000_000,
        minimum_visited_fraction,
    }
}

/// Pure flat-histogram production config (no 1/t refinement).
fn flat_histogram_config(
    max_adaptation_sweeps: u64,
    minimum_visited_fraction: f64,
) -> WangLandauConfig {
    WangLandauConfig {
        initial_log_f: 1.0,
        final_log_f: 1.0 / 4096.0,
        flatness: 0.8,
        flatness_check_interval: 25,
        discovery_sweeps: 0,
        one_over_t_threshold: 0.0,
        max_adaptation_sweeps,
        minimum_visited_fraction,
    }
}

/// Production Wang-Landau run: adapt to convergence through `axis`, then
/// return the frozen estimator state.
fn run_wang_landau(
    axis: &BinnedAxis,
    lattice: &CsrLattice,
    model: &IsingModel,
    config: &WangLandauConfig,
    seed: u64,
) -> WangLandauState {
    let mut system = System::new(lattice.clone(), 1, 1.0, 0.0);
    system.recompute_energy(model);
    assert!(
        axis.bin(system.energy).is_some(),
        "cold start energy must lie on the axis"
    );
    let mut kernel = WangLandauCore::new(*axis, config.clone()).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut guard = 0_u64;
    while kernel.estimator().is_adaptive() {
        kernel.sweep_with_phase(
            &mut system,
            model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
        guard += 1;
        assert!(guard < 1_000_000, "Wang-Landau adaptation did not converge");
    }
    assert_eq!(
        kernel.estimator().termination(),
        Some(WangLandauTermination::Converged),
        "run must end by convergence, not the sweep guard"
    );
    assert_eq!(kernel.out_of_range_proposals(), 0);
    kernel.estimator().clone()
}

fn sample_std(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    (values.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0)).sqrt()
}

#[test]
fn wang_landau_binned_dos_integrates_exact_levels_per_bin() {
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    assert_eq!(exact.states(), 64);
    let expected_levels = [-6.3, -2.3, -2.1, -1.9, 1.9, 2.1, 2.3, 6.3];
    assert_eq!(
        exact.energies().len(),
        expected_levels.len(),
        "alternating 1.0/1.1 weights split the 6-ring spectrum into 8 levels: {:?}",
        exact.energies()
    );
    for (&level, &expected) in exact.energies().iter().zip(&expected_levels) {
        super::common::assert_close(level, expected, 1e-9);
    }

    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    // Exact degeneracy integrated over each bin.
    let mut exact_counts = vec![0_u64; axis.bins()];
    for (&energy, &degeneracy) in exact.energies().iter().zip(exact.degeneracies()) {
        let bin = axis.bin(energy).expect("every level lies inside the axis");
        exact_counts[bin] += degeneracy;
    }
    // Bins {0, 4, 5, 8, 9, 13}: unit-width bins integrate two levels each in
    // bins 4 and 9; the remaining eight bins contain no physical level.
    let reachable: Vec<(usize, u64)> = exact_counts
        .iter()
        .enumerate()
        .filter(|&(_, &count)| count > 0)
        .map(|(bin, &count)| (bin, count))
        .collect();
    assert_eq!(
        reachable,
        vec![(0, 2), (4, 24), (5, 6), (8, 6), (9, 24), (13, 2)],
        "expected binned degeneracies from the 8 enumerated levels"
    );

    // Shift the exact log-DOS so the largest *reachable* value is zero.
    let max_log_count = reachable
        .iter()
        .map(|&(_, count)| (count as f64).ln())
        .fold(f64::NEG_INFINITY, f64::max);
    let exact_log_dos: Vec<Option<f64>> = exact_counts
        .iter()
        .map(|&count| (count > 0).then(|| (count as f64).ln() - max_log_count))
        .collect();

    // Only 6 of the 14 bins are physically reachable: ceil(0.4 * 14) = 6, so
    // the flatness gate demands exactly the reachable set — requiring more
    // would forbid convergence, fewer would accept an incomplete walk.
    let config = one_over_t_config(0.4);

    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let mut deviations: Vec<f64> = Vec::new();
    for seed in 0..n_seeds as u64 {
        let state = run_wang_landau(&axis, &lattice, &model, &config, seed);
        let visited: Vec<usize> = (0..axis.bins())
            .filter(|&bin| state.log_density().is_visited(bin))
            .collect();
        assert_eq!(
            visited,
            reachable.iter().map(|&(bin, _)| bin).collect::<Vec<_>>(),
            "seed {seed}: visited bins must equal the physically reachable bins"
        );
        // freeze_for_production() normalizes the estimate so the largest
        // visited log-DOS value is exactly zero; verify the gauge before use.
        let max_estimate = visited
            .iter()
            .map(|&bin| state.log_density().value(bin))
            .fold(f64::NEG_INFINITY, f64::max);
        super::common::assert_close(max_estimate, 0.0, 1e-12);
        for (bin, _) in &reachable {
            let estimate = state.log_density().value(*bin);
            let reference = exact_log_dos[*bin].unwrap();
            let deviation = (estimate - reference).abs();
            deviations.push(deviation);
            println!(
                "[wl binned seed={seed} bin={bin}] ln g = {estimate:+.4} \
                 exact = {reference:+.4} |d| = {deviation:.4}"
            );
            // Per-seed cap: a structural failure detector, not a statistical
            // gate — the worst deviation over seeds grows with the seed
            // count, so it is placed at a level (1) that is still ~7x below
            // the smallest resolvable gap between reachable bins,
            // |ln(24/6)| = 1.386, and would catch any wrong binned
            // degeneracy, missing level, or broken gauge.
            assert!(
                deviation <= 1.0,
                "seed {seed} bin {bin}: |Δ ln g| = {deviation:.4} exceeds 1.0"
            );
        }
    }
    let rms = (deviations.iter().map(|d| d * d).sum::<f64>() / deviations.len() as f64).sqrt();
    println!(
        "[wl binned] RMS |Δ ln g| = {rms:.4}, worst = {:.4}",
        deviations.iter().copied().fold(0.0, f64::max)
    );
    // Statistical gate on the seed-stable RMS.  Measured: 0.013 (n=8),
    // 0.013 (n=256); 0.05 is a ~4x margin.
    assert!(
        rms <= 0.05,
        "1/t route RMS |Δ ln g| = {rms:.4} exceeds 0.05"
    );
}

#[test]
fn wang_landau_binned_flat_histogram_route_bounds_dos_error() {
    // Pure geometric flat-histogram refinement (one_over_t_threshold = 0) is
    // the default WangLandauConfig route; its log-DOS error saturates near
    // the final ln f scale instead of vanishing.  Bound that floor explicitly
    // so the 1/t route above is known to be the accurate one, not the only
    // working one.
    let lattice = weighted_ring();
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = BinnedAxis::new(-7.0, 7.0, 14).unwrap();
    let mut exact_counts = vec![0_u64; axis.bins()];
    for (&energy, &degeneracy) in exact.energies().iter().zip(exact.degeneracies()) {
        exact_counts[axis.bin(energy).unwrap()] += degeneracy;
    }
    let max_log_count = exact_counts
        .iter()
        .copied()
        .filter(|&count| count > 0)
        .map(|count| (count as f64).ln())
        .fold(f64::NEG_INFINITY, f64::max);

    let config = flat_histogram_config(200_000, 0.4);
    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let mut deviations: Vec<f64> = Vec::new();
    for seed in 0..n_seeds as u64 {
        let state = run_wang_landau(&axis, &lattice, &model, &config, seed);
        for (bin, &count) in exact_counts.iter().enumerate() {
            if count == 0 || !state.log_density().is_visited(bin) {
                continue;
            }
            let reference = (count as f64).ln() - max_log_count;
            let deviation = (state.log_density().value(bin) - reference).abs();
            deviations.push(deviation);
            // The flat-histogram error distribution is heavy-tailed in the
            // seed count (worst over seeds: 0.16 at n=8, 0.30 at n=64, 0.54
            // at n=256), so the per-seed gate is only a structural cap well
            // below the smallest resolvable gap |ln(24/6)| = 1.386; the
            // accurate statistical gate is the RMS below.
            assert!(
                deviation <= 1.0,
                "seed {seed} bin {bin}: flat-histogram |Δ ln g| = {deviation:.4} exceeds 1.0"
            );
        }
    }
    let mean = deviations.iter().sum::<f64>() / deviations.len() as f64;
    let rms = (deviations.iter().map(|d| d * d).sum::<f64>() / deviations.len() as f64).sqrt();
    println!(
        "[wl binned flat-histogram] mean |Δ ln g| = {mean:.4}, rms = {rms:.4},          worst = {:.4}",
        deviations.iter().copied().fold(0.0, f64::max)
    );
    // Measured RMS: 0.09 (n=8) to 0.13 (n=256); 0.20 leaves margin in the
    // seed count while remaining below the 1/t route's tolerance scale.
    assert!(
        rms <= 0.20,
        "flat-histogram route RMS |Δ ln g| = {rms:.4} exceeds 0.20"
    );
}

#[test]
fn wang_landau_binned_reweights_canonical_energy_across_temperatures() {
    let lattice = build_chain(8, true);
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    // Even-cardinality antialigned-bond subsets of the 8-ring.
    assert_eq!(exact.energies(), &[-8.0, -4.0, 0.0, 4.0, 8.0]);
    assert_eq!(exact.degeneracies(), &[2, 56, 140, 56, 2]);
    // Bin centers {-8, -4, 0, +4, +8} coincide exactly with the levels, so
    // reweighting through the centers carries no binning bias.
    let axis = BinnedAxis::new(-10.0, 10.0, 5).unwrap();
    for (bin, &energy) in exact.energies().iter().enumerate() {
        assert!((axis.center(bin) - energy).abs() < 1e-12);
    }

    let config = one_over_t_config(1.0);
    let betas = [0.25, 0.6, 1.2];
    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    // One walk per seed; reweight the same frozen DOS at every temperature.
    let mut estimates: Vec<[f64; 3]> = Vec::with_capacity(n_seeds);
    for seed in 0..n_seeds as u64 {
        let state = run_wang_landau(&axis, &lattice, &model, &config, seed);
        assert_eq!(state.log_density().visited_bins(), 5);
        let mut per_seed = [0.0; 3];
        for (slot, &beta) in betas.iter().enumerate() {
            let reweighted = canonical_reweight(&axis, state.log_density(), beta).unwrap();
            per_seed[slot] = reweighted.mean_energy();
        }
        estimates.push(per_seed);
    }

    for (slot, &beta) in betas.iter().enumerate() {
        let (_, exact_energy, _, _) = exact_ising_moments(&lattice, 1.0, beta);
        let values: Vec<f64> = estimates.iter().map(|e| e[slot]).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let stderr = sample_std(&values) / (values.len() as f64).sqrt();
        let z = (mean - exact_energy) / stderr.max(1e-9);
        let max_deviation = values
            .iter()
            .map(|&v| (v - exact_energy).abs())
            .fold(0.0, f64::max);
        println!(
            "[wl reweight beta={beta}] E = {mean:.4} ± {stderr:.4} (exact {exact_energy:.4}, \
             z = {z:+.2}), worst per-seed |ΔE| = {max_deviation:.4}"
        );
        // Mean gate: a 3-sigma cross-seed band (Student-t with n-1 dof) with
        // an absolute floor of 0.02 — 4x the largest mean offset measured at
        // any seed count (0.005), and 0.25% of the energy scale — so that a
        // very large nightly seed count cannot shrink the standard error
        // below the physical resolution the test intends to assert.
        let mean_bound = (3.0 * stderr).max(0.02);
        assert!(
            (mean - exact_energy).abs() <= mean_bound,
            "beta={beta}: |⟨E⟩ - exact| = {:.4} exceeds {mean_bound:.4}",
            (mean - exact_energy).abs()
        );
        // Per-seed cap: single-seed spread is σ ≈ 0.01 with a mild heavy
        // tail (worst measured 0.086 over 256 seeds); 0.25 bounds the
        // worst-of-n for any seed count this suite supports.
        assert!(
            max_deviation <= 0.25,
            "beta={beta}: per-seed |ΔE| up to {max_deviation:.4} exceeds 0.25"
        );
    }
}

#[test]
#[ignore = "long Wang-Landau binned refinement run (nightly runs --ignored)"]
fn wang_landau_binned_long_refinement_run() {
    let lattice = build_chain(8, true);
    let model = IsingModel::new(1.0);
    let axis = BinnedAxis::new(-10.0, 10.0, 5).unwrap();
    let config = one_over_t_config(1.0);
    let betas = [0.25, 0.6, 1.2];
    let n_seeds = zscore_seed_count(32);
    let mut estimates: Vec<[f64; 3]> = Vec::with_capacity(n_seeds);
    for seed in 0..n_seeds as u64 {
        let state = run_wang_landau(&axis, &lattice, &model, &config, seed);
        let mut per_seed = [0.0; 3];
        for (slot, &beta) in betas.iter().enumerate() {
            per_seed[slot] = canonical_reweight(&axis, state.log_density(), beta)
                .unwrap()
                .mean_energy();
        }
        estimates.push(per_seed);
    }
    for (slot, &beta) in betas.iter().enumerate() {
        let (_, exact_energy, _, _) = exact_ising_moments(&lattice, 1.0, beta);
        let values: Vec<f64> = estimates.iter().map(|e| e[slot]).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let stderr = sample_std(&values) / (values.len() as f64).sqrt();
        let z = (mean - exact_energy) / stderr.max(1e-9);
        let mean_bound = (3.0 * stderr).max(0.02);
        println!(
            "[wl long reweight beta={beta}] E = {mean:.4} ± {stderr:.4} \
             (exact {exact_energy:.4}, z = {z:+.2})"
        );
        assert!(
            (mean - exact_energy).abs() <= mean_bound,
            "beta={beta}: long-run |⟨E⟩ - exact| exceeds {mean_bound:.4}"
        );
    }
}
