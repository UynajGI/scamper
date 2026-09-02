//! Exact-enumeration validation for site and bond percolation.
//!
//! The 2x2 open square has 4 sites and 4 bonds, so full enumeration covers
//! all 2^4 = 16 configurations of either mode. Expected moments are derived
//! by hand from the configuration table (independent of the implementation):
//!
//! Site mode, clusters of open sites (config counts out of 16):
//!   max: 4x1 + 4x2 + 2x1 + 4x3 + 1x4  → <max> = 30/16 = 1.875
//!   sum(s^2): 4x1 + 4x4 + 2x2 + 4x9 + 1x16 → <s2> = 76/16 = 4.75
//!   clusters: 4x1 + 4x1 + 2x2 + 4x1 + 1x1 → <n> = 17/16
//!   spanning configs {top row, bottom row, 4 triples, all} = 7/16 at p=1/2
//!
//! Bond mode, clusters over all sites (singletons included):
//!   max: 1x1 + 4x2 + 4x3 + 2x2 + 4x4 + 1x4 → <max> = 45/16 = 2.8125
//!   sum(s^2): 1x4 + 4x6 + 4x10 + 2x8 + 4x16 + 1x16 → <s2> = 164/16 = 10.25
//!   clusters: 1x4 + 4x3 + 4x2 + 2x2 + 4x1 + 1x1 → <n> = 33/16
//!   spanning configs = 12/16 at p=1/2

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{build_square, cluster_stats, OccupancyState, PercolationMC, PercolationMode};

const L: usize = 2;
const FROM: [usize; 2] = [0, 2];
const TO: [usize; 2] = [1, 3];

/// Enumerate all occupancy configurations, returning probability-weighted
/// `<max>`, `<sum(s^2)>`, `<n_clusters>` and P(spanning) at probability `p`.
fn enumerate_moments(mode: PercolationMode, p: f64) -> (f64, f64, f64, f64) {
    let lattice = build_square(L, L, false);
    let mut occupancy = OccupancyState::new(&lattice, mode);
    let n_elements = match mode {
        PercolationMode::Site => lattice.n_sites,
        PercolationMode::Bond => lattice.n_edges(),
    };
    let (mut mean_max, mut mean_s2, mut mean_n, mut p_span) = (0.0, 0.0, 0.0, 0.0);
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
        let stats = cluster_stats(&lattice, &occupancy, &FROM, &TO);
        mean_max += weight * stats.max_size as f64;
        mean_s2 += weight * stats.second_moment as f64;
        mean_n += weight * stats.n_clusters as f64;
        p_span += weight * f64::from(u8::from(stats.spanning));
    }
    (mean_max, mean_s2, mean_n, p_span)
}

#[test]
fn exact_enumeration_matches_hand_derived_moments() {
    // Site mode: (30/16, 76/16, 17/16, 7/16) at p = 1/2.
    let (max, s2, n, span) = enumerate_moments(PercolationMode::Site, 0.5);
    assert_close(max, 1.875, 1e-12, "site <max>");
    assert_close(s2, 4.75, 1e-12, "site <s^2>");
    assert_close(n, 17.0 / 16.0, 1e-12, "site <n_clusters>");
    assert_close(span, 0.4375, 1e-12, "site P(span)");

    // Bond mode: (45/16, 164/16, 33/16, 12/16) at p = 1/2.
    let (max, s2, n, span) = enumerate_moments(PercolationMode::Bond, 0.5);
    assert_close(max, 2.8125, 1e-12, "bond <max>");
    assert_close(s2, 10.25, 1e-12, "bond <s^2>");
    assert_close(n, 33.0 / 16.0, 1e-12, "bond <n_clusters>");
    assert_close(span, 0.75, 1e-12, "bond P(span)");
}

#[test]
fn exact_spanning_probability_follows_closed_form() {
    // Site mode: spanning configs are the two rows (2 of size 2), all four
    // triples and the full config: P = 2p^2(1-p)^2 + 4p^3(1-p) + p^4.
    for p in [0.2, 0.44, 0.5927, 0.8] {
        let (_, _, _, span) = enumerate_moments(PercolationMode::Site, p);
        let closed_form =
            2.0 * p * p * (1.0 - p) * (1.0 - p) + 4.0 * p.powi(3) * (1.0 - p) + p.powi(4);
        assert_close(span, closed_form, 1e-12, "site spanning polynomial");
    }
}

#[test]
fn single_site_lattice_is_degenerate_but_consistent() {
    let lattice = cmc_rs::build_square(1, 1, false);
    // Site mode with from == to = {0}: spanning iff the site is open.
    let mut site = OccupancyState::new(&lattice, PercolationMode::Site);
    site.site_open[0] = true;
    let stats = cluster_stats(&lattice, &site, &[0], &[0]);
    assert!(stats.spanning);
    assert_eq!(stats.max_size, 1);

    site.site_open[0] = false;
    let stats = cluster_stats(&lattice, &site, &[0], &[0]);
    assert!(!stats.spanning);
    assert_eq!(stats.max_size, 0);
}

/// End-to-end scheduler run: i.i.d. samples reproduce the exact 2x2 site
/// moments within 4 standard errors.
#[test]
fn scheduler_run_matches_exact_moments() {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", L);
    params.set("Ly", L);
    params.set("p", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: 200_000,
        binsize: 100,
        base_seed: 2026,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);

    let (exact_max, exact_s2, exact_n, exact_span) = enumerate_moments(PercolationMode::Site, 0.5);
    for (name, exact) in [
        ("MaxCluster", exact_max),
        ("SecondMoment", exact_s2),
        ("NClusters", exact_n),
        ("Spanning", exact_span),
    ] {
        let estimate = results
            .get(name)
            .unwrap_or_else(|| panic!("missing {name}"));
        let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
        assert!(
            z.abs() < 4.0,
            "{name}: mean {} vs exact {exact}, z = {z:.2}",
            estimate.mean
        );
    }
}

/// Long-run crossing check: at the square-lattice bond threshold p_c = 1/2
/// the crossing probability of a large L tends to 1/2.
#[test]
#[ignore = "long: bond percolation threshold crossing check (~20 s)"]
fn bond_crossing_at_critical_probability_tends_to_half() {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", 32);
    params.set("Ly", 32);
    params.set("mode", "bond");
    params.set("p", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: 200_000,
        binsize: 100,
        base_seed: 7,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);
    let spanning = results.get("Spanning").expect("Spanning measured");
    assert!(
        (spanning.mean - 0.5).abs() < 0.06,
        "P(cross) at p_c = {} deviates too far from 1/2",
        spanning.mean
    );
}

fn assert_close(value: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (value - expected).abs() < tolerance,
        "{label}: {value} vs expected {expected}"
    );
}
