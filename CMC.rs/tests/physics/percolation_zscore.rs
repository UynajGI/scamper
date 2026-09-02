//! z-score statistical validation for percolation.
//!
//! A 4x4 open square has 2^16 = 65536 site configurations — small enough to
//! enumerate exactly. The exact P(spanning) and <MaxCluster> at p = 0.6
//! become the reference for 16 independent scheduler runs (default seeds,
//! raised via `SCUTTLE_ZSCORE_SEEDS` for nightly monitoring):
//!   - Each individual seed: |z| < 4
//!   - Mean z-score across seeds: |z̄| < 1.5 (no systematic bias)
//!   - z-scores are not all same sign (no one-sided bias)

use super::common::zscore_seed_count;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{build_square, cluster_stats, OccupancyState, PercolationMC, PercolationMode};

const L: usize = 4;
const P: f64 = 0.6;
const N_SEEDS: usize = 16;
const MEASUREMENT_SWEEPS: u64 = 100_000;

/// Exact moments from full configuration enumeration of the 4x4 square.
fn exact_moments() -> (f64, f64) {
    let lattice = build_square(L, L, false);
    let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Site);
    let (mut p_span, mut mean_max) = (0.0, 0.0);
    for mask in 0..(1usize << lattice.n_sites) {
        for (bit, open) in occupancy.site_open.iter_mut().enumerate() {
            *open = (mask >> bit) & 1 == 1;
        }
        let occupied = mask.count_ones() as usize;
        let weight = P.powi(occupied as i32) * (1.0 - P).powi((lattice.n_sites - occupied) as i32);
        let stats = cluster_stats(&lattice, &occupancy, &[0, 4, 8, 12], &[3, 7, 11, 15]);
        p_span += weight * f64::from(u8::from(stats.spanning));
        mean_max += weight * stats.max_size as f64;
    }
    (p_span, mean_max)
}

fn run_percolation(seed: u64) -> Vec<(f64, f64)> {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", L);
    params.set("Ly", L);
    params.set("p", P);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: MEASUREMENT_SWEEPS,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);
    ["Spanning", "MaxCluster"]
        .into_iter()
        .map(|name| {
            let estimate = results
                .get(name)
                .unwrap_or_else(|| panic!("missing {name}"));
            (estimate.mean, estimate.stderr)
        })
        .collect()
}

fn analyze(results: &[(f64, f64)], exact: f64) -> (f64, f64, f64) {
    let z_scores: Vec<f64> = results
        .iter()
        .map(|(mean, stderr)| (mean - exact) / stderr.max(1e-12))
        .collect();
    let max_abs_z = z_scores.iter().fold(0.0_f64, |acc, z| acc.max(z.abs()));
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    let frac_pos = z_scores.iter().filter(|z| **z > 0.0).count() as f64 / z_scores.len() as f64;
    (max_abs_z, mean_z, frac_pos)
}

#[test]
fn percolation_zscore_spanning_16_seeds() {
    let (exact_span, _) = exact_moments();
    let n_seeds = zscore_seed_count(N_SEEDS);
    let runs: Vec<Vec<(f64, f64)>> = (0..n_seeds as u64).map(run_percolation).collect();
    let spanning: Vec<(f64, f64)> = runs.iter().map(|run| run[0]).collect();
    let (max_z, mean_z, frac_pos) = analyze(&spanning, exact_span);
    assert!(max_z < 4.0, "Spanning max |z| = {max_z:.2} should be < 4");
    assert!(
        mean_z.abs() < 1.5,
        "Spanning mean z = {mean_z:.2} should be |z̄| < 1.5"
    );
    assert!(
        (0.15..=0.85).contains(&frac_pos),
        "Spanning fraction positive = {frac_pos:.2} indicates one-sided bias"
    );
}

#[test]
fn percolation_zscore_max_cluster_16_seeds() {
    let (_, exact_max) = exact_moments();
    let n_seeds = zscore_seed_count(N_SEEDS);
    let runs: Vec<Vec<(f64, f64)>> = (0..n_seeds as u64).map(run_percolation).collect();
    let max_cluster: Vec<(f64, f64)> = runs.iter().map(|run| run[1]).collect();
    let (max_z, mean_z, frac_pos) = analyze(&max_cluster, exact_max);
    assert!(max_z < 4.0, "MaxCluster max |z| = {max_z:.2} should be < 4");
    assert!(
        mean_z.abs() < 1.5,
        "MaxCluster mean z = {mean_z:.2} should be |z̄| < 1.5"
    );
    assert!(
        (0.15..=0.85).contains(&frac_pos),
        "MaxCluster fraction positive = {frac_pos:.2} indicates one-sided bias"
    );
}
