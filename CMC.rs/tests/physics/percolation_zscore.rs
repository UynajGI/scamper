//! z-score statistical validation for percolation.
//!
//! All three percolation modes are validated against exact
//! full-configuration enumeration (4x4 site = 2^16 configurations, 3x3 bond
//! = 2^12 bond configurations, 3x3 site-bond = 2^21 mixed configurations —
//! computed once and shared across its two tests). Each reference feeds 16
//! independent scheduler runs (default seed count, raised via
//! `SCUTTLE_ZSCORE_SEEDS` for nightly monitoring):
//!   - Each individual seed: |z| < 4
//!   - Mean z-score across seeds: |z̄| < 1.5 (no systematic bias)
//!   - z-scores are not all same sign (no one-sided bias)

use std::sync::OnceLock;

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
/// Mixed reference: 3x3 open square at (p_site, p_bond) = (0.6, 0.7).
const SITE_BOND: (usize, f64, f64) = (3, 0.6, 0.7);

/// Exact (P(spanning), <MaxCluster>) from full configuration enumeration.
/// `side_len` is the square side (row/column width) for the spanning sets;
/// pure modes use only the probability of their sampled kind.
fn exact_moments(
    lattice: &CsrLattice,
    mode: PercolationMode,
    side_len: usize,
    p_site: f64,
    p_bond: f64,
) -> (f64, f64) {
    let n_bits = match mode {
        PercolationMode::Site => lattice.n_sites,
        PercolationMode::Bond => lattice.n_edges(),
        PercolationMode::SiteBond => lattice.n_sites + lattice.n_edges(),
    };
    let n_sites = lattice.n_sites;
    let from: Vec<usize> = (0..n_sites).filter(|site| site % side_len == 0).collect();
    let to: Vec<usize> = (0..n_sites)
        .filter(|site| site % side_len == side_len - 1)
        .collect();
    let mut occupancy = OccupancyState::new(lattice, mode);
    let (mut p_span, mut mean_max) = (0.0, 0.0);
    for mask in 0..(1usize << n_bits) {
        if mode.samples_sites() {
            for (bit, open) in occupancy.site_open.iter_mut().enumerate() {
                *open = (mask >> bit) & 1 == 1;
            }
        }
        if mode.samples_bonds() {
            let offset = match mode {
                PercolationMode::SiteBond => n_sites,
                _ => 0,
            };
            for (bit, open) in occupancy.bond_open.iter_mut().enumerate() {
                *open = (mask >> (offset + bit)) & 1 == 1;
            }
        }
        let (open_sites, open_bonds) = match mode {
            PercolationMode::Site => (mask.count_ones() as usize, 0),
            PercolationMode::Bond => (0, mask.count_ones() as usize),
            PercolationMode::SiteBond => {
                let site_mask = (1usize << n_sites) - 1;
                (
                    (mask & site_mask).count_ones() as usize,
                    (mask >> n_sites).count_ones() as usize,
                )
            }
        };
        let (ns, ne) = (lattice.n_sites, lattice.n_edges());
        let site_factor =
            p_site.powi(open_sites as i32) * (1.0 - p_site).powi((ns - open_sites) as i32);
        let bond_factor =
            p_bond.powi(open_bonds as i32) * (1.0 - p_bond).powi((ne - open_bonds) as i32);
        let weight = match mode {
            PercolationMode::Site => site_factor,
            PercolationMode::Bond => bond_factor,
            PercolationMode::SiteBond => site_factor * bond_factor,
        };
        let stats = cluster_stats(lattice, &occupancy, &from, &to);
        p_span += weight * f64::from(u8::from(stats.spanning));
        mean_max += weight * stats.max_size as f64;
    }
    (p_span, mean_max)
}

/// Exact mixed-mode reference for the 3x3 square; the 2^21-configuration
/// enumeration takes seconds, so share it across both mixed z-score tests.
fn mixed_exact() -> (f64, f64) {
    static MIXED: OnceLock<(f64, f64)> = OnceLock::new();
    *MIXED.get_or_init(|| {
        let lattice = build_square(SITE_BOND.0, SITE_BOND.0, false);
        exact_moments(
            &lattice,
            PercolationMode::SiteBond,
            SITE_BOND.0,
            SITE_BOND.1,
            SITE_BOND.2,
        )
    })
}

/// 16-seed scheduler runs; returns per-seed `(mean, stderr)` for `Spanning`
/// and `MaxCluster`. Pure modes read `p`; mixed mode reads `p_site`/`p_bond`.
fn run_seeds(mode: &str, side: usize, p_site: f64, p_bond: f64) -> Vec<[(f64, f64); 2]> {
    let n_seeds = zscore_seed_count(N_SEEDS);
    (0..n_seeds as u64)
        .map(|seed| {
            let mut params = Params::new();
            params.set("lattice_type", "square");
            params.set("Lx", side);
            params.set("Ly", side);
            params.set("mode", mode);
            if mode == "site-bond" {
                params.set("p_site", p_site);
                params.set("p_bond", p_bond);
            } else {
                params.set("p", p_site);
            }
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
    let (exact, _) = exact_moments(&lattice, PercolationMode::Site, SITE.0, SITE.1, SITE.1);
    let runs = run_seeds("site", SITE.0, SITE.1, SITE.1);
    assert_zscore_gates(&runs, 0, exact, "site Spanning");
}

#[test]
fn percolation_zscore_site_max_cluster_16_seeds() {
    let lattice = build_square(SITE.0, SITE.0, false);
    let (_, exact) = exact_moments(&lattice, PercolationMode::Site, SITE.0, SITE.1, SITE.1);
    let runs = run_seeds("site", SITE.0, SITE.1, SITE.1);
    assert_zscore_gates(&runs, 1, exact, "site MaxCluster");
}

#[test]
fn percolation_zscore_bond_spanning_16_seeds() {
    let lattice = build_square(BOND.0, BOND.0, false);
    let (exact, _) = exact_moments(&lattice, PercolationMode::Bond, BOND.0, BOND.1, BOND.1);
    let runs = run_seeds("bond", BOND.0, BOND.1, BOND.1);
    assert_zscore_gates(&runs, 0, exact, "bond Spanning");
}

#[test]
fn percolation_zscore_bond_max_cluster_16_seeds() {
    let lattice = build_square(BOND.0, BOND.0, false);
    let (_, exact) = exact_moments(&lattice, PercolationMode::Bond, BOND.0, BOND.1, BOND.1);
    let runs = run_seeds("bond", BOND.0, BOND.1, BOND.1);
    assert_zscore_gates(&runs, 1, exact, "bond MaxCluster");
}

#[test]
fn percolation_zscore_site_bond_spanning_16_seeds() {
    let (exact, _) = mixed_exact();
    let runs = run_seeds("site-bond", SITE_BOND.0, SITE_BOND.1, SITE_BOND.2);
    assert_zscore_gates(&runs, 0, exact, "site-bond Spanning");
}

#[test]
fn percolation_zscore_site_bond_max_cluster_16_seeds() {
    let (_, exact) = mixed_exact();
    let runs = run_seeds("site-bond", SITE_BOND.0, SITE_BOND.1, SITE_BOND.2);
    assert_zscore_gates(&runs, 1, exact, "site-bond MaxCluster");
}
