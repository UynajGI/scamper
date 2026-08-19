//! Per-update detailed balance for the occupation transfer-matrix sampler
//! (criterion D, statistical half).
//!
//! The deterministic half lives in the library test
//! `sweep_kernel_is_exact_heat_bath_on_closed_paths`
//! (`src/impurity/spin_boson/occupation/worldline.rs`): the sweep's bridge
//! recipe is the exact heat-bath kernel over closed worldlines at machine
//! precision, so `w(x) P(x'|x) = w(x') P(x|x')` holds as an identity.
//!
//! Here the *sampled* chain is audited empirically, in the style of the CMC
//! binned flow-balance tests: between occupation-state bins, the number of
//! observed transitions `A -> B` must match `B -> A` within Poisson errors,
//! and the stationary slice-state frequencies must match the exact
//! thermal marginals `rho_ss(beta)` from exact diagonalization.

use qmc_rs::impurity::spin_boson::occupation::transfer::SymmetricEigensystem;
use qmc_rs::{CavityMode, OccupationSpinBosonModel, OccupationWorldlineSampler};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn mode(omega: f64, coupling: f64, cutoff: usize) -> CavityMode {
    CavityMode::new(omega, coupling, cutoff).expect("mode")
}

/// Exact diagonal thermal marginal of one basis state, `rho_ss(beta)`.
fn thermal_marginal(model: &OccupationSpinBosonModel, state: usize, beta: f64) -> f64 {
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).expect("diagonalize");
    let ground = eigen.values[0];
    let mut numerator = 0.0;
    let mut partition = 0.0;
    for (level, &energy) in eigen.values.iter().enumerate() {
        let weight = (-beta * (energy - ground)).exp();
        numerator += weight * eigen.vectors[state][level] * eigen.vectors[state][level];
        partition += weight;
    }
    numerator / partition
}

/// Standard error of the mean over equal-sized blocks (jackknife-over-bins
/// house style: block means, standard deviation / sqrt(block count)).
fn blocked_stderr(values: &[f64], blocks: usize) -> f64 {
    let block_size = values.len() / blocks;
    assert!(block_size > 0, "not enough samples for {blocks} blocks");
    let mut block_means = Vec::with_capacity(blocks);
    for block in 0..blocks {
        let slice = &values[block * block_size..(block + 1) * block_size];
        block_means.push(slice.iter().sum::<f64>() / slice.len() as f64);
    }
    let mean = block_means.iter().sum::<f64>() / blocks as f64;
    let variance = block_means
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / (blocks - 1) as f64;
    variance.sqrt() / (blocks as f64).sqrt()
}

#[test]
fn occupation_empirical_flow_balance_and_stationary_marginal() {
    // Interacting Rabi model, dimension 10.
    let model = OccupationSpinBosonModel::rabi(0.9, vec![mode(1.2, 0.32, 5)]).expect("model");
    let dimension = model.basis().dimension();
    let beta = 2.0;
    let slices = 5;

    let mut sampler =
        OccupationWorldlineSampler::new(model.clone(), beta, slices, 0).expect("sampler");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0CCB_2026);
    for _ in 0..2_000 {
        sampler.sweep(&mut rng).expect("warmup");
    }

    let sweeps = 300_000_usize;
    let mut slice_zero = Vec::with_capacity(sweeps);
    for _ in 0..sweeps {
        sampler.sweep(&mut rng).expect("sweep");
        slice_zero.push(sampler.states()[0]);
    }

    // ── Empirical flow balance between distinct basis states ──────────────
    let mut transitions = vec![0_u64; dimension * dimension];
    for pair in slice_zero.windows(2) {
        transitions[pair[0] * dimension + pair[1]] += 1;
    }
    let mut checked_pairs = 0;
    for left in 0..dimension {
        for right in (left + 1)..dimension {
            let forward = transitions[left * dimension + right];
            let backward = transitions[right * dimension + left];
            let total = forward + backward;
            // Only audit pairs with enough statistics for a 4-sigma Poisson
            // bound to be meaningful.
            if total < 100 {
                continue;
            }
            checked_pairs += 1;
            let bound = 4.0 * (total as f64).sqrt();
            assert!(
                (forward as f64 - backward as f64).abs() < bound,
                "flow imbalance between states {left} and {right}: \
                 {forward} vs {backward} (4-sigma bound {bound})"
            );
        }
    }
    assert!(
        checked_pairs >= dimension,
        "only {checked_pairs} state pairs had enough transitions to audit"
    );

    // ── Stationary slice-state frequencies vs exact marginals ─────────────
    let indicators: Vec<Vec<f64>> = (0..dimension)
        .map(|target| {
            slice_zero
                .iter()
                .map(|&state| if state == target { 1.0 } else { 0.0 })
                .collect()
        })
        .collect();
    for (state, indicator) in indicators.iter().enumerate() {
        let exact = thermal_marginal(&model, state, beta);
        if exact < 0.01 {
            // Rare states give unstable blocked errors; the flow-balance
            // audit above already covers them statistically.
            continue;
        }
        let stderr = blocked_stderr(indicator, 20);
        let sampled = indicator.iter().sum::<f64>() / sweeps as f64;
        let z = (sampled - exact) / stderr.max(1.0e-12);
        assert!(
            z.abs() < 4.0,
            "state {state}: sampled frequency {sampled:.5} ± {stderr:.5} vs exact \
             marginal {exact:.5} (z = {z:.2})"
        );
    }
}
