//! Full multicanonical Monte Carlo distribution versus exact enumeration.
//!
//! The 6-site periodic Ising ring has four energy levels
//! `E ∈ {-6, -2, +2, +6}` with degeneracies `g(E) = {2, 30, 30, 2}` (one
//! even-cardinality subset of antialigned bonds per configuration, times the
//! global flip).  A bias `ln w(E) = -ln g(E)` built from that exact DOS turns
//! [`cmc_rs::EnergyBiasCore`] into a flat-histogram walk over the four
//! macrostates, which is the production route this file validates:
//!
//! 1. reweighted canonical moments `⟨E⟩(β)` and `⟨m²⟩(β)` at three
//!    temperatures versus `exact_ising_moments`, with per-seed z-scores
//!    against a delta-method standard error (repo z-score convention:
//!    `|z| < 4` per seed, `|z̄| < 1.5`, at least two seeds of each sign);
//! 2. the full reweighted energy distribution `P(E)` at one temperature
//!    against the exact Boltzmann distribution over enumerated states;
//! 3. flatness of the biased visit histogram — the defining property that
//!    distinguishes multicanonical sampling from any canonical run (at
//!    β = 0.5 the canonical bin ratio is `P(-6)/P(+6) = e^{6} ≈ 403`).
//!
//! Setting `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count for nightly
//! high-power monitoring (unset → the documented per-test default).

use super::common::{exact_ising_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, enumerate_ising_density_of_states, Algorithm, DiscreteAxis, EnergyBiasCore,
    IsingModel, MacrostateAxis, MulticanonicalBias, SimulationPhase, System,
};
use rand::SeedableRng;

const N_SITES: usize = 6;
const COUPLING: f64 = 1.0;
const DEFAULT_SEEDS: usize = 16;
const THERMALIZATION_SWEEPS: u64 = 1_000;
const MEASUREMENT_SWEEPS: usize = 24_000;
/// Block length for delta-method standard errors; far above the flat-histogram
/// autocorrelation time of the 4-level walk (a few sweeps).
const BLOCK: usize = 100;
const BETAS: [f64; 3] = [0.2, 0.5, 1.0];

/// One production run's per-sweep trajectory and biased visit histogram.
struct ProductionRun {
    energies: Vec<f64>,
    magnetization_squared: Vec<f64>,
    counts: Vec<u64>,
    out_of_range_proposals: u64,
}

fn run_production(
    seed: u64,
    thermalization: u64,
    measurement: usize,
    lattice: &cmc_rs::CsrLattice,
    model: &IsingModel,
    axis: &DiscreteAxis,
    bias: &MulticanonicalBias,
) -> ProductionRun {
    let mut system = System::new(lattice.clone(), 1, 1.0, 0.0);
    system.recompute_energy(model);
    assert!(
        axis.bin(system.energy).is_some(),
        "cold start energy must lie on the axis"
    );
    let mut kernel = EnergyBiasCore::new(axis.clone(), bias.clone());
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..thermalization {
        kernel.sweep_with_phase(
            &mut system,
            model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    kernel.clear_histogram();
    let mut energies = Vec::with_capacity(measurement);
    let mut magnetization_squared = Vec::with_capacity(measurement);
    for _ in 0..measurement {
        kernel.sweep_with_phase(&mut system, model, &mut rng, SimulationPhase::Measurement);
        energies.push(system.energy);
        let magnetization = system.spins.iter().sum::<f64>() / N_SITES as f64;
        magnetization_squared.push(magnetization * magnetization);
    }
    ProductionRun {
        energies,
        magnetization_squared,
        counts: kernel.histogram().counts().to_vec(),
        out_of_range_proposals: kernel.out_of_range_proposals(),
    }
}

fn sample_variance(samples: &[f64]) -> f64 {
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    samples
        .iter()
        .map(|&x| (x - mean) * (x - mean))
        .sum::<f64>()
        / (n - 1.0)
}

fn sample_covariance(left: &[f64], right: &[f64]) -> f64 {
    let n = left.len() as f64;
    let left_mean = left.iter().sum::<f64>() / n;
    let right_mean = right.iter().sum::<f64>() / n;
    left.iter()
        .zip(right)
        .map(|(&l, &r)| (l - left_mean) * (r - right_mean))
        .sum::<f64>()
        / (n - 1.0)
}

fn block_means(samples: &[f64], block: usize) -> Vec<f64> {
    assert!(samples.len().is_multiple_of(block));
    samples
        .chunks(block)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect()
}

/// Canonical estimate `Σ A_i g(E_i) e^{-βE_i} / Σ g(E_i) e^{-βE_i}` with a
/// delta-method standard error built from `BLOCK`-sweep block means of
/// numerator and denominator: `Var(R) ≈ (Var(Ā) − 2R·Cov(Ā,B̄) + R²·Var(B̄))/B̄²`.
///
/// The biased walk samples each level uniformly; unbiasing divides out the
/// applied weight `1/g(E)`, so every per-sweep sample carries the canonical
/// importance weight `g(E) e^{-βE}`.
fn reweighted_mean_and_stderr(
    energies: &[f64],
    values: &[f64],
    axis: &DiscreteAxis,
    log_density: &cmc_rs::LogDensityOfStates,
    beta: f64,
    block: usize,
) -> (f64, f64) {
    assert_eq!(energies.len(), values.len());
    let weights: Vec<f64> = energies
        .iter()
        .map(|&e| {
            let bin = axis.bin(e).expect("visited energy lies on the axis");
            (log_density.value(bin) - beta * e).exp()
        })
        .collect();
    let numerator: Vec<f64> = values
        .iter()
        .zip(&weights)
        .map(|(&value, &weight)| value * weight)
        .collect();
    let numerator_blocks = block_means(&numerator, block);
    let denominator_blocks = block_means(&weights, block);
    let n_blocks = numerator_blocks.len() as f64;
    let numerator_mean = numerator_blocks.iter().sum::<f64>() / n_blocks;
    let denominator_mean = denominator_blocks.iter().sum::<f64>() / n_blocks;
    let ratio = numerator_mean / denominator_mean;
    let var_numerator = sample_variance(&numerator_blocks) / n_blocks;
    let var_denominator = sample_variance(&denominator_blocks) / n_blocks;
    let covariance = sample_covariance(&numerator_blocks, &denominator_blocks) / n_blocks;
    let var_ratio = ((var_numerator - 2.0 * ratio * covariance + ratio * ratio * var_denominator)
        .max(0.0))
        / (denominator_mean * denominator_mean);
    (ratio, var_ratio.sqrt().max(1e-12))
}

/// Repo z-score convention: each seed within 4σ, no systematic offset, and
/// at least two seeds of each sign (a one-sided bias would fail the sign
/// check; `P(≤1 sign flip | null) = (1+n)/2^(n-1)`, i.e. 0.1% at n=16).
fn assert_zscores(z_scores: &[f64], label: &str) {
    let max_abs = z_scores.iter().map(|z| z.abs()).fold(0.0, f64::max);
    assert!(
        max_abs < 4.0,
        "{label}: max |z| = {max_abs:.2} should be < 4"
    );
    let mean = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    assert!(
        mean.abs() < 1.5,
        "{label}: mean z = {mean:.2} should be |z̄| < 1.5"
    );
    let n_positive = z_scores.iter().filter(|&&z| z > 0.0).count();
    assert!(
        (2..=z_scores.len() - 2).contains(&n_positive),
        "{label}: {n_positive}/{} seeds positive; at least two of each sign required",
        z_scores.len()
    );
}

#[test]
fn multicanonical_production_reweighted_moments_match_exact() {
    let lattice = build_chain(N_SITES, true);
    let model = IsingModel::new(COUPLING);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    // Anchor: even-cardinality antialigned-bond subsets of the 6-ring.
    assert_eq!(exact.energies(), &[-6.0, -2.0, 2.0, 6.0]);
    assert_eq!(exact.degeneracies(), &[2, 30, 30, 2]);
    assert_eq!(exact.states(), 1 << N_SITES);
    let axis = exact.axis().unwrap();
    let dos = exact.log_density().unwrap();
    let bias = MulticanonicalBias::from_log_density(&dos).unwrap();

    // One production run per seed; the same trajectory is reweighted at every
    // temperature.
    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let exact_moments: Vec<(f64, f64)> = BETAS
        .iter()
        .map(|&beta| {
            let (_, exact_energy, _, exact_m2) = exact_ising_moments(&lattice, COUPLING, beta);
            (exact_energy, exact_m2)
        })
        .collect();
    let mut z_energy: Vec<Vec<f64>> = BETAS.iter().map(|_| Vec::new()).collect();
    let mut z_m2: Vec<Vec<f64>> = BETAS.iter().map(|_| Vec::new()).collect();
    for seed in 0..n_seeds as u64 {
        let run = run_production(
            seed,
            THERMALIZATION_SWEEPS,
            MEASUREMENT_SWEEPS,
            &lattice,
            &model,
            &axis,
            &bias,
        );
        assert_eq!(run.out_of_range_proposals, 0);
        for (slot, &beta) in BETAS.iter().enumerate() {
            let (energy, energy_stderr) =
                reweighted_mean_and_stderr(&run.energies, &run.energies, &axis, &dos, beta, BLOCK);
            let (m2, m2_stderr) = reweighted_mean_and_stderr(
                &run.energies,
                &run.magnetization_squared,
                &axis,
                &dos,
                beta,
                BLOCK,
            );
            let (exact_energy, exact_m2) = exact_moments[slot];
            z_energy[slot].push((energy - exact_energy) / energy_stderr);
            z_m2[slot].push((m2 - exact_m2) / m2_stderr);
        }
    }
    for (slot, &beta) in BETAS.iter().enumerate() {
        let mean_z_energy = z_energy[slot].iter().sum::<f64>() / n_seeds as f64;
        let mean_z_m2 = z_m2[slot].iter().sum::<f64>() / n_seeds as f64;
        println!(
            "[muca beta={beta}] mean z(E) = {mean_z_energy:+.2}, mean z(m2) = {mean_z_m2:+.2}, \
             max |z| = {:.2}",
            z_energy[slot]
                .iter()
                .zip(&z_m2[slot])
                .map(|(a, b)| a.abs().max(b.abs()))
                .fold(0.0, f64::max),
        );
        assert_zscores(&z_energy[slot], &format!("⟨E⟩ at beta={beta}"));
        assert_zscores(&z_m2[slot], &format!("⟨m²⟩ at beta={beta}"));
    }
}

#[test]
fn multicanonical_production_energy_distribution_matches_boltzmann() {
    let lattice = build_chain(N_SITES, true);
    let model = IsingModel::new(COUPLING);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = exact.axis().unwrap();
    let dos = exact.log_density().unwrap();
    let bias = MulticanonicalBias::from_log_density(&dos).unwrap();
    let beta = 0.5;

    // Exact Boltzmann distribution over the enumerated levels.
    let weights: Vec<f64> = exact
        .energies()
        .iter()
        .zip(exact.degeneracies())
        .map(|(&energy, &degeneracy)| degeneracy as f64 * (-beta * energy).exp())
        .collect();
    let partition: f64 = weights.iter().sum();
    let exact_p: Vec<f64> = weights.iter().map(|&w| w / partition).collect();

    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let mut accumulated: Vec<f64> = vec![0.0; exact.energies().len()];
    let mut total_weight = 0.0;
    let mut samples = 0_usize;
    for seed in 0..n_seeds as u64 {
        let run = run_production(
            seed,
            THERMALIZATION_SWEEPS,
            MEASUREMENT_SWEEPS,
            &lattice,
            &model,
            &axis,
            &bias,
        );
        samples += run.energies.len();
        for &energy in &run.energies {
            let bin = axis.bin(energy).expect("visited energy lies on the axis");
            // Unbias: divide out the applied 1/g(E) weight, then apply the
            // canonical Boltzmann factor.
            let weight = (dos.value(bin) - beta * energy).exp();
            accumulated[bin] += weight;
            total_weight += weight;
        }
    }
    let estimated_p: Vec<f64> = accumulated.iter().map(|&w| w / total_weight).collect();

    // Conservative effective sample size: the flat walk over four levels
    // decorrelates within a few sweeps; counting one effective sample per 50
    // sweeps is a >= 10x margin against autocorrelation inflation.  The bound
    // is a 4-sigma binomial band for each reweighted bin probability.
    let effective_samples = (samples / 50).max(1) as f64;
    let mut total_variation = 0.0;
    for ((&estimated, &reference), &energy) in
        estimated_p.iter().zip(&exact_p).zip(exact.energies())
    {
        let deviation = (estimated - reference).abs();
        total_variation += deviation;
        let bound = 4.0 * (reference * (1.0 - reference) / effective_samples).sqrt() + 1e-4;
        println!(
            "[muca P(E) E={energy:+.0}] estimated={estimated:.5} exact={reference:.5} \
             |dp|={deviation:.5} bound={bound:.5}"
        );
        assert!(
            deviation <= bound,
            "P(E={energy:+.0}): |Δp| = {deviation:.5} exceeds binomial band {bound:.5}"
        );
    }
    let tv = total_variation / 2.0;
    println!("[muca P(E)] total variation distance = {tv:.5}");
    assert!(tv <= 0.02, "total variation distance {tv:.5} exceeds 0.02");
}

#[test]
fn multicanonical_production_histogram_is_flat() {
    let lattice = build_chain(N_SITES, true);
    let model = IsingModel::new(COUPLING);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = exact.axis().unwrap();
    let bias = MulticanonicalBias::from_log_density(&exact.log_density().unwrap()).unwrap();

    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let mut worst_ratio = f64::INFINITY;
    for seed in 0..n_seeds as u64 {
        let run = run_production(
            seed,
            THERMALIZATION_SWEEPS,
            MEASUREMENT_SWEEPS,
            &lattice,
            &model,
            &axis,
            &bias,
        );
        assert_eq!(run.out_of_range_proposals, 0);
        let visited = run.counts.iter().filter(|&&c| c > 0).count();
        assert_eq!(visited, exact.energies().len(), "all macrostates visited");
        let mean = run.counts.iter().sum::<u64>() as f64 / visited as f64;
        let minimum = *run.counts.iter().min().unwrap() as f64;
        let ratio = minimum / mean;
        worst_ratio = worst_ratio.min(ratio);
        println!("[muca flat seed={seed}] min/mean = {ratio:.4}");
        // Canonical sampling at any beta >= 0.2 spends < 3% of the time above
        // E = +2; a flat multicanonical histogram cannot dip below 90% of the
        // mean per bin (binomial 4σ with a ≥10x autocorrelation margin).
        assert!(
            ratio >= 0.90,
            "biased histogram not flat: min/mean = {ratio:.4}, counts = {:?}",
            run.counts
        );
    }
    println!("[muca flat] worst min/mean across seeds = {worst_ratio:.4}");
}

#[test]
#[ignore = "long multicanonical statistical run (nightly runs --ignored)"]
fn multicanonical_long_convergence_run() {
    let lattice = build_chain(N_SITES, true);
    let model = IsingModel::new(COUPLING);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = exact.axis().unwrap();
    let dos = exact.log_density().unwrap();
    let bias = MulticanonicalBias::from_log_density(&dos).unwrap();
    let n_seeds = zscore_seed_count(32);
    for &beta in &BETAS {
        let (_, exact_energy, _, exact_m2) = exact_ising_moments(&lattice, COUPLING, beta);
        let mut z_energy = Vec::with_capacity(n_seeds);
        let mut z_m2 = Vec::with_capacity(n_seeds);
        for seed in 0..n_seeds as u64 {
            let run = run_production(seed, 4_000, 96_000, &lattice, &model, &axis, &bias);
            let (energy, energy_stderr) =
                reweighted_mean_and_stderr(&run.energies, &run.energies, &axis, &dos, beta, 400);
            let (m2, m2_stderr) = reweighted_mean_and_stderr(
                &run.energies,
                &run.magnetization_squared,
                &axis,
                &dos,
                beta,
                400,
            );
            z_energy.push((energy - exact_energy) / energy_stderr);
            z_m2.push((m2 - exact_m2) / m2_stderr);
        }
        assert_zscores(&z_energy, &format!("long ⟨E⟩ at beta={beta}"));
        assert_zscores(&z_m2, &format!("long ⟨m²⟩ at beta={beta}"));
    }
}
