//! T2 heterogeneous (Poisson-binomial) Newman-Ziff reweighting validation.
//!
//! Four layers, following the Z1 gate style:
//! 1. exact full-enumeration z-scores against a TEST-LOCAL BFS reference
//!    (never production union-find) for a moderate-spread U(0.3,0.8) field
//!    and a bimodal half-0.2/half-0.9 field, site and bond scans;
//! 2. the exact-unit-mass diagnostic (mean exactly 1 over scans) and
//!    finiteness of every returned weight sum at N = 10^4 mid-range fields;
//! 3. at-scale N1 cross-validation: reweighted site scans on L=64/128
//!    squares with a clipped linear ramp field centered at p_c vs direct
//!    `SiteBernoulli` sampling at the same field;
//! 4. uniform-field reduction to the v1 binomial `reweight` and typed
//!    rejections.
//!
//! The bimodal field is heavy-tailed at larger universes (the estimator's
//! log-weight excursion scale is sqrt((N/4) Var(logit p))), which is why the
//! exact gates use small graphs and the graded-field gates use narrow ramps;
//! see the reweight_heterogeneous rustdoc variance caveat.

use cmc_rs::{
    analyze, build_chain, build_square, reweight, reweight_heterogeneous, scan_bond, scan_site,
    Bond, BondType, ComponentWorkspace, CsrLattice, GraphView, HetReweightWorkspace,
    NewmanZiffWorkspace, Probability, ProbabilityField, ReweightError, ReweightWorkspace,
    ScanEstimator, ScanMode, ScanObservable, SiteBernoulli, StaticConfiguration,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::zscore_seed_count;

/// Labeled seed bases ("HET2"/"HEZ2") so these streams are distinct from the
/// Z1 and N1 z-score families.
const EXACT_SEED_BASE: u64 = 0x4845_5432_0000;
const SCALE_SEED_BASE: u64 = 0x4845_5A32_0000;
/// Heavy-tailed bimodal values need larger batches than the Z1 default to
/// keep the exact-enumeration z-gates reliable.
const EXACT_SCANS: usize = 200 * 16;
const SCALE_SCANS: usize = 25 * 16;
const DIRECT_SAMPLES: usize = 4096;
/// Ramp half-widths: the per-scan reweighting values are roughly lognormal
/// with parameter sqrt((N/4) Var(logit p)), so the graded-field gates use
/// ramps narrowed like 1/sqrt(N) to keep that parameter below ~0.5.
const P_C: f64 = 0.5927;
const RAMP_HALF_L64: f64 = 0.004;
const RAMP_HALF_L128: f64 = 0.002;

const OBSERVABLES: [ScanObservable; 2] = [
    ScanObservable::GiantFraction,
    ScanObservable::RawSecondMoment,
];

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("test probability")
}

fn field(values: &[f64]) -> ProbabilityField<'_> {
    ProbabilityField::Owned(values.iter().map(|&value| probability(value)).collect())
}

/// Moderate-spread deterministic U(0.3, 0.8) field.
fn spread_values(n: usize) -> Vec<f64> {
    (0..n)
        .map(|index| 0.3 + 0.5 * (((index * 37 + 11) % 97) as f64) / 96.0)
        .collect()
}

/// Bimodal half-0.2/half-0.9 field.
fn bimodal_values(n: usize) -> Vec<f64> {
    (0..n)
        .map(|index| if index % 2 == 0 { 0.2 } else { 0.9 })
        .collect()
}

/// Clipped linear ramp in the column coordinate, centered at p_c.
fn graded_values(lattice: &CsrLattice, side: usize, half: f64) -> Vec<f64> {
    (0..lattice.vertex_count())
        .map(|index| {
            let column = (index % side) as f64 / (side - 1) as f64 - 0.5;
            (P_C + half * column / 0.5).clamp(0.25, 0.85)
        })
        .collect()
}

fn edge_list(lattice: &CsrLattice) -> Vec<(usize, usize)> {
    lattice
        .edges
        .iter()
        .map(|bond| (bond.source, bond.target))
        .collect()
}

fn random_graph(n: usize) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(6);
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            if rng.random::<f64>() < 0.4 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            }
        }
    }
    CsrLattice::try_from_edges(n, edges).expect("random graph is valid")
}

/// Test-local BFS reference over an explicit site occupancy and edge list —
/// fully independent of the production union-find and analyzer.
fn bfs_largest_and_moment(
    vertex_count: usize,
    sites: &[bool],
    bonds: &[(usize, usize)],
) -> (usize, u128) {
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
    for &(left, right) in bonds {
        if sites[left] && sites[right] {
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
    }
    let mut seen = vec![false; vertex_count];
    let mut largest = 0;
    let mut moment = 0_u128;
    for seed in 0..vertex_count {
        if !sites[seed] || seen[seed] {
            continue;
        }
        seen[seed] = true;
        let mut queue = vec![seed];
        let mut size = 0;
        while let Some(vertex) = queue.pop() {
            size += 1;
            for &neighbor in &adjacency[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
        largest = largest.max(size);
        moment += (size as u128) * (size as u128);
    }
    (largest, moment)
}

/// Exact `(E[giant fraction], E[raw second moment])` under the heterogeneous
/// product measure, by full 2^N enumeration with the local BFS reference.
fn exact_heterogeneous_means(lattice: &CsrLattice, mode: ScanMode, values: &[f64]) -> (f64, f64) {
    let vertices = lattice.vertex_count();
    let bonds = edge_list(lattice);
    let bits = match mode {
        ScanMode::Site => vertices,
        ScanMode::Bond => bonds.len(),
    };
    let mut giant = 0.0;
    let mut moment = 0.0;
    let mut weight_sum = 0.0;
    for mask in 0_usize..(1_usize << bits) {
        let mut weight = 1.0;
        let mut sites = vec![true; vertices];
        let active_bonds: Vec<(usize, usize)> = match mode {
            ScanMode::Site => {
                for (index, open) in sites.iter_mut().enumerate() {
                    *open = (mask >> index) & 1 == 1;
                    weight *= if *open {
                        values[index]
                    } else {
                        1.0 - values[index]
                    };
                }
                bonds.clone()
            }
            ScanMode::Bond => bonds
                .iter()
                .enumerate()
                .filter_map(|(index, &edge)| {
                    let open = (mask >> index) & 1 == 1;
                    weight *= if open {
                        values[index]
                    } else {
                        1.0 - values[index]
                    };
                    open.then_some(edge)
                })
                .collect(),
        };
        let (largest, second) = bfs_largest_and_moment(vertices, &sites, &active_bonds);
        weight_sum += weight;
        giant += weight * largest as f64 / vertices as f64;
        moment += weight * second as f64;
    }
    assert!((weight_sum - 1.0).abs() < 1e-12, "weights sum to one");
    (giant, moment)
}

fn heterogeneous_estimators(
    lattice: &CsrLattice,
    mode: ScanMode,
    values: &[f64],
    scans: usize,
    seed_base: u64,
) -> ([ScanEstimator; 2], ScanEstimator) {
    let mut estimators = [ScanEstimator::new(); 2];
    let mut mass_estimator = ScanEstimator::new();
    let mut scan_workspace = NewmanZiffWorkspace::new();
    scan_workspace
        .prepare(lattice.vertex_count(), lattice.n_edges())
        .unwrap();
    let mut reweight_workspace = HetReweightWorkspace::new();
    let probability_field = field(values);

    for index in 0..scans {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed_base + index as u64);
        let curve = match mode {
            ScanMode::Site => scan_site(lattice, &mut scan_workspace, &mut rng).unwrap(),
            ScanMode::Bond => scan_bond(lattice, &mut scan_workspace, &mut rng).unwrap(),
        };
        for (metric_index, metric) in OBSERVABLES.iter().enumerate() {
            let result = reweight_heterogeneous(
                &curve,
                &probability_field,
                *metric,
                &mut reweight_workspace,
            )
            .unwrap();
            estimators[metric_index].add(result.value);
            if metric_index == 0 {
                mass_estimator.add(result.unit_mass);
            }
        }
    }
    (estimators, mass_estimator)
}

#[test]
fn heterogeneous_reweighting_matches_exact_enumeration() {
    let graphs = [
        ("random-8-site", random_graph(8), ScanMode::Site),
        ("chain-8-bond", build_chain(8, false), ScanMode::Bond),
    ];
    let scans = EXACT_SCANS * zscore_seed_count(16) / 16;
    let mut z_scores = Vec::new();

    for (graph_index, (name, graph, mode)) in graphs.iter().enumerate() {
        let spread = match mode {
            ScanMode::Site => spread_values(graph.vertex_count()),
            ScanMode::Bond => spread_values(graph.n_edges()),
        };
        let bimodal = bimodal_values(spread.len());
        for (field_index, values) in [spread, bimodal].into_iter().enumerate() {
            let (exact_giant, exact_moment) = exact_heterogeneous_means(graph, *mode, &values);
            let (estimators, mass) = heterogeneous_estimators(
                graph,
                *mode,
                &values,
                scans,
                EXACT_SEED_BASE + (graph_index * 4 + field_index) as u64 * 0x1_0000,
            );
            for (metric_index, metric) in OBSERVABLES.iter().enumerate() {
                let expected = if metric_index == 0 {
                    exact_giant
                } else {
                    exact_moment
                };
                let estimate = estimators[metric_index]
                    .estimate()
                    .expect("at least two scans");
                assert_eq!(estimate.samples, scans);
                let z = (estimate.mean - expected) / estimate.standard_error;
                assert!(
                    z.abs() <= 4.0,
                    "{name} field {field_index} {metric:?}: z = {z:.3}, mean = {}, exact = \
                     {expected}, sem = {}",
                    estimate.mean,
                    estimate.standard_error
                );
                eprintln!(
                    "het exact: {name} field{field_index} {}: z = {z:+.3} (mean {}, exact \
                     {expected}, sem {})",
                    metric.name(),
                    estimate.mean,
                    estimate.standard_error
                );
                z_scores.push(z);
            }
            let mass_estimate = mass.estimate().expect("at least two scans");
            let mass_z = (mass_estimate.mean - 1.0) / mass_estimate.standard_error;
            assert!(
                mass_z.abs() <= 4.0,
                "{name} field {field_index}: unit-mass z = {mass_z:.3} (mean {})",
                mass_estimate.mean
            );
            eprintln!(
                "het exact: {name} field{field_index} unit-mass z = {mass_z:+.3} (mean {}, sem \
                 {}, R = {scans})",
                mass_estimate.mean, mass_estimate.standard_error
            );
        }
    }
    eprintln!(
        "het exact: {} z-scores, max |z| = {:.3}",
        z_scores.len(),
        z_scores.iter().fold(0.0_f64, |acc, z| acc.max(z.abs()))
    );
}

#[test]
fn heterogeneous_unit_mass_gate_and_finiteness_at_scale() {
    // Exact-unit-mass gate on the L=128 ramp field of the cross-validation
    // test below: the scan average of the diagnostic must sit within 4 SEM
    // of exactly 1.
    let lattice = build_square(128, 128, false);
    let values = graded_values(&lattice, 128, RAMP_HALF_L128);
    let (_, mass) = heterogeneous_estimators(
        &lattice,
        ScanMode::Site,
        &values,
        SCALE_SCANS * zscore_seed_count(16) / 16,
        SCALE_SEED_BASE,
    );
    let estimate = mass.estimate().expect("at least two scans");
    let z = (estimate.mean - 1.0) / estimate.standard_error;
    assert!(
        z.abs() <= 4.0,
        "L=128 ramp unit-mass z = {z:.3} (mean {}, sem {})",
        estimate.mean,
        estimate.standard_error
    );
    eprintln!(
        "het unit-mass: L=128 ramp mean = {:.6}, sem = {:.6}, z = {z:+.3}",
        estimate.mean, estimate.standard_error
    );

    // Every internal and returned weight stays finite at N = 10^4 across
    // mid-range fields, including spreads far past the estimator's useful
    // envelope (finiteness is a numerics claim, not an accuracy claim).
    let large = build_square(100, 100, false);
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut reweight_workspace = HetReweightWorkspace::new();
    let fields = [
        spread_values(large.vertex_count()),
        bimodal_values(large.vertex_count()),
        vec![0.5; large.vertex_count()],
    ];
    for (field_index, values) in fields.iter().enumerate() {
        let probability_field = field(values);
        for seed in 0..4 {
            let mut rng =
                Xoshiro256PlusPlus::seed_from_u64(0x4649_4E00 + 16 * field_index as u64 + seed);
            let curve = scan_site(&large, &mut scan_workspace, &mut rng).unwrap();
            let result = reweight_heterogeneous(
                &curve,
                &probability_field,
                ScanObservable::GiantFraction,
                &mut reweight_workspace,
            )
            .unwrap();
            assert!(
                result.value.is_finite() && result.unit_mass.is_finite(),
                "field {field_index} seed {seed}: value {} mass {}",
                result.value,
                result.unit_mass
            );
        }
    }
    eprintln!("het finiteness: N=10^4 spread/bimodal/constant fields all finite (4 seeds each)");
}

/// Direct N1 reference: `SiteBernoulli` at the same field plus one F3
/// analysis per sample. Returns (means, sems) for giant fraction and raw
/// second moment.
fn direct_estimates(
    lattice: &CsrLattice,
    values: &[f64],
    samples: usize,
    seed: u64,
) -> ([f64; 2], [f64; 2]) {
    let probabilities = values
        .iter()
        .map(|&value| probability(value))
        .collect::<Vec<_>>();
    let law = SiteBernoulli::new(ProbabilityField::Borrowed(&probabilities));
    let mut configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.n_edges());
    let mut workspace = ComponentWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut sums = [0.0_f64; 2];
    let mut sums_sq = [0.0_f64; 2];
    let vertices = lattice.vertex_count() as f64;
    for _ in 0..samples {
        law.sample(lattice, &mut configuration, &mut rng).unwrap();
        let summary = analyze(lattice, &configuration, &configuration, &[], &mut workspace)
            .unwrap()
            .summary();
        let values = [
            summary.largest_component_size as f64 / vertices,
            summary.raw_second_moment as f64,
        ];
        for index in 0..2 {
            sums[index] += values[index];
            sums_sq[index] += values[index] * values[index];
        }
    }
    let count = samples as f64;
    let mut means = [0.0; 2];
    let mut sems = [0.0; 2];
    for index in 0..2 {
        means[index] = sums[index] / count;
        let variance =
            ((sums_sq[index] - sums[index] * sums[index] / count) / (count - 1.0)).max(0.0);
        sems[index] = (variance / count).sqrt();
    }
    (means, sems)
}

#[test]
fn heterogeneous_reweighting_matches_direct_sampling_at_scale() {
    let scans = SCALE_SCANS * zscore_seed_count(16) / 16;
    for (side, half) in [(64_usize, RAMP_HALF_L64), (128, RAMP_HALF_L128)] {
        let lattice = build_square(side, side, false);
        let values = graded_values(&lattice, side, half);
        let (estimators, _) = heterogeneous_estimators(
            &lattice,
            ScanMode::Site,
            &values,
            scans,
            SCALE_SEED_BASE + 0x1_0000 * u64::from(side == 128),
        );
        let (direct_means, direct_sems) =
            direct_estimates(&lattice, &values, DIRECT_SAMPLES, 0x4449_5200 + side as u64);
        for (metric_index, metric) in OBSERVABLES.iter().enumerate() {
            let estimate = estimators[metric_index]
                .estimate()
                .expect("at least two scans");
            let combined =
                (estimate.standard_error.powi(2) + direct_sems[metric_index].powi(2)).sqrt();
            let z = (estimate.mean - direct_means[metric_index]) / combined;
            assert!(
                z.abs() <= 4.0,
                "L={side} {}: z = {z:.3}, reweighted {} +/- {} vs direct {} +/- {}",
                metric.name(),
                estimate.mean,
                estimate.standard_error,
                direct_means[metric_index],
                direct_sems[metric_index]
            );
            eprintln!(
                "het at-scale L={side} {}: z = {z:+.3}, reweighted {:.6} +/- {:.6} vs direct \
                 {:.6} +/- {:.6} (R = {scans}, direct = {DIRECT_SAMPLES}, seed base {:#x})",
                metric.name(),
                estimate.mean,
                estimate.standard_error,
                direct_means[metric_index],
                direct_sems[metric_index],
                SCALE_SEED_BASE + 0x1_0000 * u64::from(side == 128)
            );
        }
    }
}

#[test]
fn heterogeneous_uniform_field_reduces_to_binomial_reweight() {
    let lattice = build_square(4, 4, false);
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut het_workspace = HetReweightWorkspace::new();
    let mut binomial_workspace = ReweightWorkspace::new();

    for mode in [ScanMode::Site, ScanMode::Bond] {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5E_4400 + mode as u64);
        let curve = match mode {
            ScanMode::Site => scan_site(&lattice, &mut scan_workspace, &mut rng).unwrap(),
            ScanMode::Bond => scan_bond(&lattice, &mut scan_workspace, &mut rng).unwrap(),
        };
        let trials = curve.trials();
        for probability_value in [0.1, 0.25, 0.5, 0.5927, 0.75, 0.9] {
            let constant = vec![probability_value; trials];
            let het = reweight_heterogeneous(
                &curve,
                &field(&constant),
                ScanObservable::LargestComponentSize,
                &mut het_workspace,
            )
            .unwrap();
            let binomial = reweight(
                &curve,
                probability(probability_value),
                ScanObservable::LargestComponentSize,
                &mut binomial_workspace,
            )
            .unwrap();
            let relative = ((het.value - binomial) / binomial).abs();
            assert!(
                relative < 1e-9,
                "{mode:?} p = {probability_value}: relative difference {relative:e}"
            );
            assert!((het.unit_mass - 1.0).abs() < 1e-9);

            // The Uniform variant delegates to the binomial path outright.
            let delegated = reweight_heterogeneous(
                &curve,
                &probability(probability_value).into(),
                ScanObservable::LargestComponentSize,
                &mut het_workspace,
            )
            .unwrap();
            assert_eq!(delegated.value, binomial);
            assert_eq!(delegated.unit_mass, 1.0);
        }
    }
}

#[test]
fn heterogeneous_rejects_invalid_fields_with_indices() {
    let lattice = build_chain(6, false);
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut reweight_workspace = HetReweightWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5E_5245);
    let curve = scan_site(&lattice, &mut scan_workspace, &mut rng).unwrap();

    let short = field(&[0.5, 0.5, 0.5]);
    assert!(matches!(
        reweight_heterogeneous(
            &curve,
            &short,
            ScanObservable::ComponentCount,
            &mut reweight_workspace
        )
        .unwrap_err(),
        ReweightError::ProbabilityFieldLength {
            expected: 6,
            actual: 3
        }
    ));

    let mut zero_values = vec![0.5; 6];
    zero_values[4] = 0.0;
    assert!(matches!(
        reweight_heterogeneous(
            &curve,
            &field(&zero_values),
            ScanObservable::ComponentCount,
            &mut reweight_workspace
        )
        .unwrap_err(),
        ReweightError::BoundaryProbability {
            index: 4,
            value: 0.0
        }
    ));

    let mut one_values = vec![0.5; 6];
    one_values[2] = 1.0;
    assert!(matches!(
        reweight_heterogeneous(
            &curve,
            &field(&one_values),
            ScanObservable::ComponentCount,
            &mut reweight_workspace
        )
        .unwrap_err(),
        ReweightError::BoundaryProbability {
            index: 2,
            value: 1.0
        }
    ));

    // NaN probabilities are unrepresentable: the field is built from
    // validated `Probability` values only.
    assert!(Probability::new(f64::NAN).is_err());
    assert!(Probability::new(f64::INFINITY).is_err());
}
