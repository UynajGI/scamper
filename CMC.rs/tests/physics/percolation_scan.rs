//! Z1 Newman-Ziff scan validation.
//!
//! Three layers:
//! 1. every prefix of the activation order is checked against a fresh
//!    `analyze` pass over the same occupancy;
//! 2. reweighted scan estimates are checked against full 3x3 enumeration with
//!    `SCUTTLE_ZSCORE_SEEDS`-scalable batches of independent scans;
//! 3. reweighted estimates are checked against direct Bernoulli sampling
//!    within a combined confidence interval.

use cmc_rs::{
    analyze, build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square,
    build_triangular, reweight, Bond, BondType, BorrowedUndirectedCsr, ComponentSummary,
    ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId, NewmanZiffWorkspace, Probability,
    ReweightWorkspace, ScanEstimator, ScanMode, ScanObservable, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::zscore_seed_count;

const METRICS: [ScanObservable; 5] = [
    ScanObservable::ComponentCount,
    ScanObservable::ActiveEdgeCount,
    ScanObservable::LargestComponentSize,
    ScanObservable::GiantFraction,
    ScanObservable::RawSecondMoment,
];

struct Occupancy {
    mode: ScanMode,
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl Occupancy {
    fn new(graph: &CsrLattice, mode: ScanMode) -> Self {
        Self {
            mode,
            sites: vec![matches!(mode, ScanMode::Bond); graph.n_sites],
            bonds: vec![matches!(mode, ScanMode::Site); graph.n_edges()],
        }
    }

    fn from_prefix(graph: &CsrLattice, mode: ScanMode, order: &[usize], count: usize) -> Self {
        let mut occupancy = Self::new(graph, mode);
        match mode {
            ScanMode::Site => {
                for &vertex in &order[..count] {
                    occupancy.sites[vertex] = true;
                }
            }
            ScanMode::Bond => {
                for &edge in &order[..count] {
                    occupancy.bonds[edge] = true;
                }
            }
        }
        occupancy
    }

    fn from_mask(graph: &CsrLattice, mode: ScanMode, mask: usize) -> Self {
        let mut occupancy = Self::new(graph, mode);
        match mode {
            ScanMode::Site => {
                for (index, open) in occupancy.sites.iter_mut().enumerate() {
                    *open = (mask >> index) & 1 == 1;
                }
            }
            ScanMode::Bond => {
                for (index, open) in occupancy.bonds.iter_mut().enumerate() {
                    *open = (mask >> index) & 1 == 1;
                }
            }
        }
        occupancy
    }

    fn sample(&mut self, probability: f64, rng: &mut Xoshiro256PlusPlus) {
        match self.mode {
            ScanMode::Site => {
                for open in &mut self.sites {
                    *open = rng.random::<f64>() < probability;
                }
            }
            ScanMode::Bond => {
                for open in &mut self.bonds {
                    *open = rng.random::<f64>() < probability;
                }
            }
        }
    }
}

impl VertexActivity for Occupancy {
    fn vertex_count(&self) -> usize {
        self.sites.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.sites[vertex.index()]
    }
}

impl EdgeActivity for Occupancy {
    fn edge_count(&self) -> usize {
        self.bonds.len()
    }

    fn edge_active(&self, edge: EdgeId) -> bool {
        self.bonds[edge.index()]
    }
}

fn summary(graph: &CsrLattice, occupancy: &Occupancy) -> ComponentSummary {
    let mut workspace = ComponentWorkspace::new();
    analyze(graph, occupancy, occupancy, &[], &mut workspace)
        .unwrap()
        .summary()
}

fn random_graph() -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
    let n = 10;
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            if rng.random::<f64>() < 0.35 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            }
        }
    }
    CsrLattice::try_from_edges(n, edges).expect("random graph is valid")
}

fn multigraph() -> CsrLattice {
    CsrLattice::try_from_edges(
        3,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
        ],
    )
    .expect("multigraph is valid")
}

fn assert_prefix_parity(graph: &CsrLattice, mode: ScanMode, seeds: u64) {
    let mut scan_workspace = NewmanZiffWorkspace::new();
    for seed in 0..seeds {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let curve = match mode {
            ScanMode::Site => scan_site_curve(graph, &mut scan_workspace, &mut rng),
            ScanMode::Bond => scan_bond_curve(graph, &mut scan_workspace, &mut rng),
        };
        for (occupation, point) in curve.points().iter().enumerate() {
            let occupancy =
                Occupancy::from_prefix(graph, mode, curve.activation_order(), occupation);
            let expected = summary(graph, &occupancy);
            assert_eq!(point.occupation, occupation);
            assert_eq!(
                point.active_edge_count, expected.active_edge_count,
                "{mode:?} prefix {occupation}: active edges"
            );
            assert_eq!(
                point.component_count, expected.component_count,
                "{mode:?} prefix {occupation}: components"
            );
            assert_eq!(
                point.largest_component_size, expected.largest_component_size,
                "{mode:?} prefix {occupation}: largest"
            );
            assert_eq!(
                point.raw_second_moment, expected.raw_second_moment,
                "{mode:?} prefix {occupation}: second moment"
            );
            if mode == ScanMode::Site {
                assert_eq!(expected.active_vertex_count, occupation);
            } else {
                assert_eq!(point.active_edge_count, occupation);
                assert_eq!(expected.active_vertex_count, graph.n_sites);
            }
        }
    }
}

fn scan_site_curve<'a>(
    graph: &CsrLattice,
    workspace: &'a mut NewmanZiffWorkspace,
    rng: &mut Xoshiro256PlusPlus,
) -> cmc_rs::ScanCurve<'a> {
    cmc_rs::scan_site(graph, workspace, rng).unwrap()
}

fn scan_bond_curve<'a>(
    graph: &CsrLattice,
    workspace: &'a mut NewmanZiffWorkspace,
    rng: &mut Xoshiro256PlusPlus,
) -> cmc_rs::ScanCurve<'a> {
    cmc_rs::scan_bond(graph, workspace, rng).unwrap()
}

#[test]
fn scan_prefixes_match_naive_analysis() {
    let cases = [
        ("chain", build_chain(8, false)),
        ("square", build_square(4, 4, false)),
        (
            "cubic",
            build_hypercubic(
                &[2, 2, 2],
                &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                false,
            ),
        ),
        ("triangular", build_triangular(3, 3)),
        ("honeycomb", build_honeycomb(2, 2)),
        ("kagome", build_kagome(2, 2)),
        ("random", random_graph()),
        ("multigraph", multigraph()),
    ];
    for (name, graph) in cases {
        for mode in [ScanMode::Site, ScanMode::Bond] {
            assert_prefix_parity(&graph, mode, 4);
        }
        eprintln!("scan prefix parity: {name} ok");
    }
}

#[test]
fn scan_handles_empty_borrowed_graph() {
    let graph = BorrowedUndirectedCsr::new(&[0], &[], &[], &[]).unwrap();
    let mut workspace = NewmanZiffWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);

    let site = cmc_rs::scan_site(&graph, &mut workspace, &mut rng).unwrap();
    assert_eq!(site.trials(), 0);
    assert_eq!(site.points().len(), 1);
    let bond = cmc_rs::scan_bond(&graph, &mut workspace, &mut rng).unwrap();
    assert_eq!(bond.trials(), 0);
    assert_eq!(bond.points().len(), 1);

    let mut reweight_workspace = ReweightWorkspace::new();
    let error = reweight(
        &bond,
        Probability::new(0.5).unwrap(),
        ScanObservable::GiantFraction,
        &mut reweight_workspace,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        cmc_rs::ReweightError::UndefinedObservable(ScanObservable::GiantFraction)
    ));
}

#[derive(Clone, Copy)]
struct ExactConfig {
    occupation: usize,
    active_edges: usize,
    component_count: usize,
    largest: usize,
    second_moment: u128,
}

fn enumerate_configs(graph: &CsrLattice, mode: ScanMode) -> Vec<ExactConfig> {
    let bits = match mode {
        ScanMode::Site => graph.n_sites,
        ScanMode::Bond => graph.n_edges(),
    };
    let mut workspace = ComponentWorkspace::new();
    (0..(1usize << bits))
        .map(|mask| {
            let occupancy = Occupancy::from_mask(graph, mode, mask);
            let summary = analyze(graph, &occupancy, &occupancy, &[], &mut workspace)
                .unwrap()
                .summary();
            ExactConfig {
                occupation: mask.count_ones() as usize,
                active_edges: summary.active_edge_count,
                component_count: summary.component_count,
                largest: summary.largest_component_size,
                second_moment: summary.raw_second_moment,
            }
        })
        .collect()
}

fn exact_means(
    configs: &[ExactConfig],
    bits: usize,
    total_vertices: usize,
    probability: f64,
) -> [f64; 5] {
    let mut means = [0.0; 5];
    let mut weight_sum = 0.0;
    for config in configs {
        let weight = probability.powi(config.occupation as i32)
            * (1.0 - probability).powi((bits - config.occupation) as i32);
        weight_sum += weight;
        means[0] += weight * config.component_count as f64;
        means[1] += weight * config.active_edges as f64;
        means[2] += weight * config.largest as f64;
        means[3] += weight * config.largest as f64 / total_vertices as f64;
        means[4] += weight * config.second_moment as f64;
    }
    assert!((weight_sum - 1.0).abs() < 1e-12, "weights sum to one");
    means
}

fn scan_estimators(
    graph: &CsrLattice,
    mode: ScanMode,
    probabilities: &[f64],
    scans: usize,
    seed_base: u64,
) -> Vec<[ScanEstimator; 5]> {
    let mut estimators = vec![[ScanEstimator::new(); 5]; probabilities.len()];
    let mut scan_workspace = NewmanZiffWorkspace::new();
    scan_workspace
        .prepare(graph.n_sites, graph.n_edges())
        .unwrap();
    let mut reweight_workspace = ReweightWorkspace::new();

    for index in 0..scans {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed_base + index as u64);
        let curve = match mode {
            ScanMode::Site => scan_site_curve(graph, &mut scan_workspace, &mut rng),
            ScanMode::Bond => scan_bond_curve(graph, &mut scan_workspace, &mut rng),
        };
        for (probability_index, &probability) in probabilities.iter().enumerate() {
            let probability = Probability::new(probability).unwrap();
            for (metric_index, metric) in METRICS.iter().enumerate() {
                let value =
                    reweight(&curve, probability, *metric, &mut reweight_workspace).unwrap();
                estimators[probability_index][metric_index].add(value);
            }
        }
    }
    estimators
}

#[test]
fn reweighted_scans_match_exact_enumeration() {
    let graph = build_square(3, 3, false);
    let probabilities = [0.3, 0.5, 0.7];
    let scans = 16 * zscore_seed_count(16);
    let mut z_scores = Vec::new();

    for (mode_index, mode) in [ScanMode::Site, ScanMode::Bond].into_iter().enumerate() {
        let configs = enumerate_configs(&graph, mode);
        let bits = match mode {
            ScanMode::Site => graph.n_sites,
            ScanMode::Bond => graph.n_edges(),
        };
        let estimates = scan_estimators(
            &graph,
            mode,
            &probabilities,
            scans,
            0x5CA9_0000 + mode_index as u64 * 0x1_0000,
        );

        for (probability_index, &probability) in probabilities.iter().enumerate() {
            let exact = exact_means(&configs, bits, graph.n_sites, probability);
            for (metric_index, metric) in METRICS.iter().enumerate() {
                let estimate = estimates[probability_index][metric_index]
                    .estimate()
                    .expect("at least two scans");
                assert_eq!(estimate.samples, scans);
                let expected = exact[metric_index];
                if estimate.standard_error <= 1e-14 {
                    assert!(
                        (estimate.mean - expected).abs() < 1e-9,
                        "{mode:?} p={probability} {metric:?}: deterministic mean {} vs {expected}",
                        estimate.mean
                    );
                } else {
                    let z = (estimate.mean - expected) / estimate.standard_error;
                    assert!(
                        z.abs() < 4.5,
                        "{mode:?} p={probability} {metric:?}: z={z}, mean={}, exact={expected}, se={}",
                        estimate.mean,
                        estimate.standard_error
                    );
                    z_scores.push(z);
                }
            }
        }
    }

    let positive = z_scores.iter().filter(|&&z| z > 0.0).count();
    let fraction = positive as f64 / z_scores.len() as f64;
    assert!(
        (0.15..=0.85).contains(&fraction),
        "one-sided z-score bias: {positive}/{} positive",
        z_scores.len()
    );
}

fn direct_estimates(
    graph: &CsrLattice,
    mode: ScanMode,
    probability: f64,
    samples: usize,
    seed: u64,
) -> ([f64; 2], [f64; 2]) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut occupancy = Occupancy::new(graph, mode);
    let mut sum = [0.0; 2];
    let mut sum_sq = [0.0; 2];
    for _ in 0..samples {
        occupancy.sample(probability, &mut rng);
        let summary = summary(graph, &occupancy);
        let values = [
            summary.component_count as f64,
            summary.largest_component_size as f64,
        ];
        for index in 0..2 {
            sum[index] += values[index];
            sum_sq[index] += values[index] * values[index];
        }
    }
    let count = samples as f64;
    let mut means = [0.0; 2];
    let mut errors = [0.0; 2];
    for index in 0..2 {
        means[index] = sum[index] / count;
        let variance = ((sum_sq[index] - sum[index] * sum[index] / count) / (count - 1.0)).max(0.0);
        errors[index] = (variance / count).sqrt();
    }
    (means, errors)
}

#[test]
fn reweighted_scans_match_direct_bernoulli_within_combined_ci() {
    let graph = build_square(3, 3, false);
    let probability = 0.55;
    let scans = 16 * zscore_seed_count(16);
    for (mode_index, mode) in [ScanMode::Site, ScanMode::Bond].into_iter().enumerate() {
        let estimates = scan_estimators(
            &graph,
            mode,
            &[probability],
            scans,
            0x6A5B_0000 + mode_index as u64 * 0x1_0000,
        );
        let (direct_mean, direct_error) =
            direct_estimates(&graph, mode, probability, 8_192, 0x9911);
        for (direct_index, metric_index) in [(0, 0), (1, 2)] {
            let scan = estimates[0][metric_index].estimate().unwrap();
            let combined =
                4.0 * (scan.standard_error.powi(2) + direct_error[direct_index].powi(2)).sqrt();
            assert!(
                (scan.mean - direct_mean[direct_index]).abs() <= combined,
                "{mode:?} metric {metric_index}: scan {} +/- {} vs direct {} +/- {}",
                scan.mean,
                scan.standard_error,
                direct_mean[direct_index],
                direct_error[direct_index]
            );
        }
    }
}

#[test]
fn reweight_endpoints_are_exact() {
    let graph = build_square(3, 3, false);
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut reweight_workspace = ReweightWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(31);
    let curve = scan_site_curve(&graph, &mut scan_workspace, &mut rng);

    for metric in METRICS {
        let zero = reweight(
            &curve,
            Probability::new(0.0).unwrap(),
            metric,
            &mut reweight_workspace,
        )
        .unwrap();
        let one = reweight(
            &curve,
            Probability::new(1.0).unwrap(),
            metric,
            &mut reweight_workspace,
        )
        .unwrap();
        let first = curve.points()[0];
        let last = curve.points()[curve.trials()];
        assert_eq!(zero, expected_value(first, metric, graph.n_sites));
        assert_eq!(one, expected_value(last, metric, graph.n_sites));
    }
}

fn expected_value(point: cmc_rs::ScanPoint, metric: ScanObservable, total_vertices: usize) -> f64 {
    match metric {
        ScanObservable::ComponentCount => point.component_count as f64,
        ScanObservable::ActiveEdgeCount => point.active_edge_count as f64,
        ScanObservable::LargestComponentSize => point.largest_component_size as f64,
        ScanObservable::GiantFraction => {
            point.largest_component_size as f64 / total_vertices as f64
        }
        ScanObservable::RawSecondMoment => point.raw_second_moment as f64,
    }
}

#[test]
fn fixed_seed_is_reproducible_and_workspace_reuse_is_clean() {
    let graph = build_chain(10, false);
    let mut workspace = NewmanZiffWorkspace::new();
    workspace.prepare(graph.n_sites, graph.n_edges()).unwrap();

    let mut first_rng = Xoshiro256PlusPlus::seed_from_u64(77);
    let first = scan_site_curve(&graph, &mut workspace, &mut first_rng);
    let first_order = first.activation_order().to_vec();
    let first_points = first.points().to_vec();

    let mut second_rng = Xoshiro256PlusPlus::seed_from_u64(77);
    let second = scan_site_curve(&graph, &mut workspace, &mut second_rng);
    assert_eq!(second.activation_order(), first_order.as_slice());
    assert_eq!(second.points(), first_points.as_slice());

    let mut shifted_rng = Xoshiro256PlusPlus::seed_from_u64(78);
    let shifted = scan_site_curve(&graph, &mut workspace, &mut shifted_rng);
    assert_ne!(shifted.activation_order(), first_order.as_slice());

    let mut bond_rng = Xoshiro256PlusPlus::seed_from_u64(79);
    let reused_bond = scan_bond_curve(&graph, &mut workspace, &mut bond_rng);
    let reused_order = reused_bond.activation_order().to_vec();
    let reused_points = reused_bond.points().to_vec();

    let mut fresh_workspace = NewmanZiffWorkspace::new();
    let mut fresh_rng = Xoshiro256PlusPlus::seed_from_u64(79);
    let fresh_bond = scan_bond_curve(&graph, &mut fresh_workspace, &mut fresh_rng);
    assert_eq!(fresh_bond.activation_order(), reused_order.as_slice());
    assert_eq!(fresh_bond.points(), reused_points.as_slice());
}
