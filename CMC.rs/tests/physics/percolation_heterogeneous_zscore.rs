//! Multi-seed z-score statistical validation for heterogeneous (N1) laws.
//!
//! Heterogeneous site, bond, and mixed percolation run through the full
//! scheduler stack (`StaticPercolationMC` + `Run::from_parts`, because the
//! `FromParams` schema stays uniform-only) on small open squares with
//! element-wise distinct probabilities — genuinely non-identically
//! distributed fields. Two exact reference families feed the gates:
//!
//! - **Full configuration enumeration** with per-element weights
//!   `prod_i p_i^{x_i} (1 - p_i)^{1 - x_i}`: 3x3 site = 2^9, 3x3 bond =
//!   2^12, and 3x3 mixed = 2^21 configurations (computed once and shared
//!   across its two tests). Gates `BoundaryCrossing0` and `LargestSize`.
//! - **Poisson-binomial activity means**: per sample,
//!   `<ActiveVertexCount> = sum p_v` and
//!   `<ActiveEdgeCount> = sum_e q_e p_u(e) p_v(e)` (mixed couples to both
//!   endpoints), so the adapter activity counts are checked against exact
//!   first moments without enumeration; z denominators are the per-run
//!   binned stderr.
//!
//! Each reference feeds 16 independent scheduler runs (default seed count,
//! raised via `SCUTTLE_ZSCORE_SEEDS` for nightly monitoring):
//!   - Each individual seed: |z| < 4
//!   - Mean z-score across seeds: |z̄| < 1.5 (no systematic bias)
//!   - z-scores are not all same sign (no one-sided bias)

use std::sync::OnceLock;

use super::common::zscore_seed_count;
use carlo_rs::{Context, Run, RunConfig, RunId, TaskId};
use cmc_rs::{
    analyze, build_square, BoundaryQuery, ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId,
    GraphView, ObservablePlan, Probability, StaticLaw, StaticPercolationMC, UndirectedGraphView,
    VertexActivity, VertexId,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const N_SEEDS: usize = 16;
/// Labeled seed base for the heterogeneous scheduler runs ("HETE"), so the
/// streams are distinct from the uniform z-score family and from the law-level
/// probes.
const SEED_BASE: u64 = 0x4845_5445_0000;
const MEASUREMENT_SWEEPS: u64 = 100_000;
/// Heterogeneous reference: 3x3 open square, element-wise distinct
/// probabilities cycling these tables (vertex 4-cycles, edge 3-cycles).
const SITE_P: [f64; 4] = [0.15, 0.35, 0.65, 0.85];
const BOND_P: [f64; 3] = [0.2, 0.45, 0.7];
const SIDE: usize = 3;

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("test probability")
}

fn vertex_field(lattice: &CsrLattice) -> Vec<f64> {
    (0..lattice.vertex_count())
        .map(|index| SITE_P[index % SITE_P.len()])
        .collect()
}

fn edge_field(lattice: &CsrLattice) -> Vec<f64> {
    (0..lattice.edge_count())
        .map(|index| BOND_P[index % BOND_P.len()])
        .collect()
}

struct EnumState {
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl EnumState {
    fn new(lattice: &CsrLattice) -> Self {
        Self {
            sites: vec![true; lattice.vertex_count()],
            bonds: vec![true; lattice.edge_count()],
        }
    }
}

impl VertexActivity for EnumState {
    fn vertex_count(&self) -> usize {
        self.sites.len()
    }
    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.sites[vertex.index()]
    }
}

impl EdgeActivity for EnumState {
    fn edge_count(&self) -> usize {
        self.bonds.len()
    }
    fn edge_active(&self, edge: EdgeId) -> bool {
        self.bonds[edge.index()]
    }
}

/// Exact `(P(spanning), <LargestSize>)` from full configuration enumeration
/// with per-element heterogeneous weights (sites low, bonds high).
fn exact_moments(lattice: &CsrLattice, site_p: &[f64], bond_p: &[f64]) -> (f64, f64) {
    let samples_sites = !site_p.is_empty();
    let samples_bonds = !bond_p.is_empty();
    let n_bits = site_p.len() + bond_p.len();
    let from: Vec<usize> = (0..lattice.vertex_count())
        .filter(|site| site % SIDE == 0)
        .collect();
    let to: Vec<usize> = (0..lattice.vertex_count())
        .filter(|site| site % SIDE == SIDE - 1)
        .collect();
    let queries = [BoundaryQuery::new(lattice, &from, &to).unwrap()];
    let mut workspace = ComponentWorkspace::new();
    let mut state = EnumState::new(lattice);
    let (mut p_span, mut mean_max) = (0.0, 0.0);
    for mask in 0..(1usize << n_bits) {
        if samples_sites {
            for (bit, open) in state.sites.iter_mut().enumerate() {
                *open = (mask >> bit) & 1 == 1;
            }
        }
        if samples_bonds {
            let offset = site_p.len();
            for (bit, open) in state.bonds.iter_mut().enumerate() {
                *open = (mask >> (offset + bit)) & 1 == 1;
            }
        }
        let weight = |values: &[f64], offset: usize| {
            values
                .iter()
                .enumerate()
                .map(|(bit, &p)| {
                    if (mask >> (offset + bit)) & 1 == 1 {
                        p
                    } else {
                        1.0 - p
                    }
                })
                .product::<f64>()
        };
        let total = weight(site_p, 0) * weight(bond_p, site_p.len());
        let result = analyze(lattice, &state, &state, &queries, &mut workspace).unwrap();
        p_span += total * f64::from(u8::from(result.query_outcomes()[0]));
        mean_max += total * result.summary().largest_component_size as f64;
    }
    (p_span, mean_max)
}

/// The 2^21 mixed enumeration takes seconds, so share it across both mixed
/// z-score tests.
fn mixed_exact() -> (f64, f64) {
    static MIXED: OnceLock<(f64, f64)> = OnceLock::new();
    *MIXED.get_or_init(|| {
        let lattice = build_square(SIDE, SIDE, false);
        exact_moments(&lattice, &vertex_field(&lattice), &edge_field(&lattice))
    })
}

/// One heterogeneous scheduler run through `Run::from_parts`.
fn run_once(law: &StaticLaw, seed: u64, plan: &ObservablePlan) -> carlo_rs::Results {
    let lattice = build_square(SIDE, SIDE, false);
    let adapter = StaticPercolationMC::new(lattice.clone(), law.clone(), plan.clone()).unwrap();
    let context =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(SEED_BASE + seed), 0, 100);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: MEASUREMENT_SWEEPS,
        binsize: 100,
        base_seed: SEED_BASE + seed,
        ..Default::default()
    };
    let mut run = Run::from_parts(context, adapter, TaskId::new(0), RunId::new(0), config);
    run.run(MEASUREMENT_SWEEPS);
    run.finalize(SEED_BASE + seed)
}

fn spanning_plan() -> ObservablePlan {
    let lattice = build_square(SIDE, SIDE, false);
    let from: Vec<usize> = (0..lattice.vertex_count())
        .filter(|site| site % SIDE == 0)
        .collect();
    let to: Vec<usize> = (0..lattice.vertex_count())
        .filter(|site| site % SIDE == SIDE - 1)
        .collect();
    ObservablePlan::all(vec![BoundaryQuery::new(&lattice, &from, &to).unwrap()])
}

/// 16-seed scheduler runs; returns per-seed `(mean, stderr)` for `Spanning`
/// and `MaxCluster`.
fn run_seeds(law: &StaticLaw) -> Vec<[(f64, f64); 2]> {
    let n_seeds = zscore_seed_count(N_SEEDS);
    let plan = spanning_plan();
    (0..n_seeds as u64)
        .map(|seed| {
            let results = run_once(law, seed, &plan);
            ["BoundaryCrossing0", "LargestSize"].map(|name| {
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
fn heterogeneous_zscore_site_spanning_16_seeds() {
    let lattice = build_square(SIDE, SIDE, false);
    let (exact, _) = exact_moments(&lattice, &vertex_field(&lattice), &[]);
    let law = StaticLaw::site_heterogeneous(
        vertex_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 0, exact, "heterogeneous site Spanning");
}

#[test]
fn heterogeneous_zscore_site_max_cluster_16_seeds() {
    let lattice = build_square(SIDE, SIDE, false);
    let (_, exact) = exact_moments(&lattice, &vertex_field(&lattice), &[]);
    let law = StaticLaw::site_heterogeneous(
        vertex_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 1, exact, "heterogeneous site MaxCluster");
}

#[test]
fn heterogeneous_zscore_bond_spanning_16_seeds() {
    let lattice = build_square(SIDE, SIDE, false);
    let (exact, _) = exact_moments(&lattice, &[], &edge_field(&lattice));
    let law = StaticLaw::bond_heterogeneous(
        edge_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 0, exact, "heterogeneous bond Spanning");
}

#[test]
fn heterogeneous_zscore_bond_max_cluster_16_seeds() {
    let lattice = build_square(SIDE, SIDE, false);
    let (_, exact) = exact_moments(&lattice, &[], &edge_field(&lattice));
    let law = StaticLaw::bond_heterogeneous(
        edge_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 1, exact, "heterogeneous bond MaxCluster");
}

#[test]
fn heterogeneous_zscore_mixed_spanning_16_seeds() {
    let (exact, _) = mixed_exact();
    let lattice = build_square(SIDE, SIDE, false);
    let law = StaticLaw::mixed_heterogeneous(
        vertex_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
        edge_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 0, exact, "heterogeneous mixed Spanning");
}

#[test]
fn heterogeneous_zscore_mixed_max_cluster_16_seeds() {
    let (_, exact) = mixed_exact();
    let lattice = build_square(SIDE, SIDE, false);
    let law = StaticLaw::mixed_heterogeneous(
        vertex_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
        edge_field(&lattice)
            .iter()
            .map(|&value| probability(value))
            .collect::<Vec<_>>(),
    );
    let runs = run_seeds(&law);
    assert_zscore_gates(&runs, 1, exact, "heterogeneous mixed MaxCluster");
}

/// Poisson-binomial activity gate: through the scheduler, mean occupied
/// counts of heterogeneous fields match exact first moments within |z| < 4
/// per seed (stderr comes from each run's own binned statistics).
///
/// - site law: `ActiveVertexCount` is a Poisson-binomial count with mean
///   `sum_i p_i`;
/// - bond law: `ActiveEdgeCount` is a Poisson-binomial count with mean
///   `sum_j q_j`;
/// - mixed law: `ActiveEdgeCount` counts an edge only when the edge and both
///   endpoint sites are open, so its exact mean is
///   `sum_j q_j p_u(j) p_v(j)` — independent elements, endpoint-coupled.
#[test]
fn heterogeneous_zscore_adapter_activity_matches_poisson_binomial() {
    let lattice = build_square(SIDE, SIDE, false);
    let site_values = vertex_field(&lattice);
    let bond_values = edge_field(&lattice);
    let site_mean: f64 = site_values.iter().sum();
    let bond_mean: f64 = bond_values.iter().sum();
    let coupled_bond_mean: f64 = lattice
        .edges
        .iter()
        .enumerate()
        .map(|(index, edge)| {
            bond_values[index] * site_values[edge.source] * site_values[edge.target]
        })
        .sum();
    let site_field = site_values
        .iter()
        .map(|&value| probability(value))
        .collect::<Vec<_>>();
    let bond_field = bond_values
        .iter()
        .map(|&value| probability(value))
        .collect::<Vec<_>>();
    let plan = ObservablePlan::all(vec![
        BoundaryQuery::new(&lattice, &[0, 3, 6], &[2, 5, 8]).unwrap()
    ]);
    let n_seeds = zscore_seed_count(N_SEEDS);
    for seed in 0..n_seeds as u64 {
        let site = run_once(
            &StaticLaw::site_heterogeneous(site_field.clone()),
            seed,
            &plan,
        );
        let bond = run_once(
            &StaticLaw::bond_heterogeneous(bond_field.clone()),
            seed,
            &plan,
        );
        let mixed = run_once(
            &StaticLaw::mixed_heterogeneous(site_field.clone(), bond_field.clone()),
            seed,
            &plan,
        );
        for (results, name, exact) in [
            (&site, "ActiveVertexCount", site_mean),
            (&bond, "ActiveEdgeCount", bond_mean),
            (&mixed, "ActiveVertexCount", site_mean),
            (&mixed, "ActiveEdgeCount", coupled_bond_mean),
        ] {
            let estimate = results
                .get(name)
                .unwrap_or_else(|| panic!("missing {name}"));
            let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
            assert!(
                z.abs() < 4.0,
                "{name} seed {seed}: {} vs {exact}, z = {z:.2}",
                estimate.mean
            );
        }
    }
}
