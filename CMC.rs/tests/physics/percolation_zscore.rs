//! z-score statistical validation for percolation.
//!
//! Both percolation modes are validated against exact full-configuration
//! enumeration (4x4 site = 2^16 configurations, 3x3 bond = 2^12 bond
//! configurations). Each reference feeds 16 independent scheduler runs
//! (default seed count, raised via `SCUTTLE_ZSCORE_SEEDS` for nightly
//! monitoring):
//!   - Each individual seed: |z| < 4
//!   - Mean z-score across seeds: |z̄| < 1.5 (no systematic bias)
//!   - z-scores are not all same sign (no one-sided bias)

use super::common::zscore_seed_count;
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    build_square, cluster_stats, CsrLattice, OccupancyState, PercolationMC, PercolationMode,
};

const N_SEEDS: usize = 16;
const MEASUREMENT_SWEEPS: u64 = 100_000;
/// Site-mode reference: 4x4 open square at p = 0.6, left/right columns.
const SITE: (usize, f64) = (4, 0.6);
/// Bond-mode reference: 3x3 open square at p = 0.55, near the bond p_c = 1/2
/// where spanning variance (and therefore statistical power) is largest.
const BOND: (usize, f64) = (3, 0.55);

/// Exact (P(spanning), <MaxCluster>) from full configuration enumeration.
/// `side_len` is the square side (row/column width) for the spanning sets.
fn exact_moments(
    lattice: &CsrLattice,
    mode: PercolationMode,
    side_len: usize,
    p: f64,
) -> (f64, f64) {
    let n_elements = match mode {
        PercolationMode::Site => lattice.n_sites,
        PercolationMode::Bond => lattice.n_edges(),
    };
    let from: Vec<usize> = (0..lattice.n_sites)
        .filter(|site| site % side_len == 0)
        .collect();
    let to: Vec<usize> = (0..lattice.n_sites)
        .filter(|site| site % side_len == side_len - 1)
        .collect();
    let mut occupancy = OccupancyState::new(lattice, mode);
    let (mut p_span, mut mean_max) = (0.0, 0.0);
    for mask in 0..(1usize << n_elements) {
        match mode {
            PercolationMode::Site => {
                for (bit, open) in occupancy.site_open.iter_mut().enumerate() {
                    *open = (mask >> bit) & 1 == 1;
                }
            }
            PercolationMode::Bond => {
                for (bit, open) in occupancy.bond_open.iter_mut().enumerate() {
                    *open = (mask >> bit) & 1 == 1;
                }
            }
        }
        let occupied = mask.count_ones() as usize;
        let weight = p.powi(occupied as i32) * (1.0 - p).powi((n_elements - occupied) as i32);
        let stats = cluster_stats(lattice, &occupancy, &from, &to);
        p_span += weight * f64::from(u8::from(stats.spanning));
        mean_max += weight * stats.max_size as f64;
    }
    (p_span, mean_max)
}

/// 16-seed scheduler runs; returns per-seed `(mean, stderr)` for `Spanning`
/// and `MaxCluster`.
fn run_seeds(mode: &str, side: usize, p: f64) -> Vec<[(f64, f64); 2]> {
    let n_seeds = zscore_seed_count(N_SEEDS);
    (0..n_seeds as u64)
        .map(|seed| {
            let mut params = Params::new();
            params.set("lattice_type", "square");
            params.set("Lx", side);
            params.set("Ly", side);
            params.set("mode", mode);
            params.set("p", p);
            let config = RunConfig {
                thermalization_sweeps: 0,
                measurement_sweeps: MEASUREMENT_SWEEPS,
                binsize: 100,
                base_seed: seed,
                ..Default::default()
            };
            let results =
                Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);
            ["Spanning", "MaxCluster"].map(|name| {
                let estimate = results
                    .get(name)
                    .unwrap_or_else(|| panic!("missing {name}"));
                (estimate.mean, estimate.stderr)
            })
        })
        .collect()
}

/// Assert the standard z-score gates: per-seed |z| < 4, |z̄| < 1.5, and no
/// one-sided bias.
fn assert_zscore_gates(runs: &[[(f64, f64); 2]], observable: usize, exact: f64, label: &str) {
    let z_scores: Vec<f64> = runs
        .iter()
        .map(|run| {
            let (mean, stderr) = run[observable];
            (mean - exact) / stderr.max(1e-12)
        })
        .collect();
    let max_abs_z = z_scores.iter().fold(0.0_f64, |acc, z| acc.max(z.abs()));
    let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
    let frac_pos = z_scores.iter().filter(|z| **z > 0.0).count() as f64 / z_scores.len() as f64;
    assert!(
        max_abs_z < 4.0,
        "{label} max |z| = {max_abs_z:.2} should be < 4"
    );
    assert!(
        mean_z.abs() < 1.5,
        "{label} mean z = {mean_z:.2} should be |z̄| < 1.5"
    );
    assert!(
        (0.15..=0.85).contains(&frac_pos),
        "{label} fraction positive = {frac_pos:.2} indicates one-sided bias"
    );
}

#[test]
fn percolation_zscore_site_spanning_16_seeds() {
    let lattice = build_square(SITE.0, SITE.0, false);
    let (exact, _) = exact_moments(&lattice, PercolationMode::Site, SITE.0, SITE.1);
    let runs = run_seeds("site", SITE.0, SITE.1);
    assert_zscore_gates(&runs, 0, exact, "site Spanning");
}

#[test]
fn percolation_zscore_site_max_cluster_16_seeds() {
    let lattice = build_square(SITE.0, SITE.0, false);
    let (_, exact) = exact_moments(&lattice, PercolationMode::Site, SITE.0, SITE.1);
    let runs = run_seeds("site", SITE.0, SITE.1);
    assert_zscore_gates(&runs, 1, exact, "site MaxCluster");
}

#[test]
fn percolation_zscore_bond_spanning_16_seeds() {
    let lattice = build_square(BOND.0, BOND.0, false);
    let (exact, _) = exact_moments(&lattice, PercolationMode::Bond, BOND.0, BOND.1);
    let runs = run_seeds("bond", BOND.0, BOND.1);
    assert_zscore_gates(&runs, 0, exact, "bond Spanning");
}

#[test]
fn percolation_zscore_bond_max_cluster_16_seeds() {
    let lattice = build_square(BOND.0, BOND.0, false);
    let (_, exact) = exact_moments(&lattice, PercolationMode::Bond, BOND.0, BOND.1);
    let runs = run_seeds("bond", BOND.0, BOND.1);
    assert_zscore_gates(&runs, 1, exact, "bond MaxCluster");
}
