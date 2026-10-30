//! N2 production validation: random network failure and degree observables.
//!
//! Gates (roadmap §9 N2), each anchored to published results:
//!
//! 1. **Site-law reduction**: `RandomFailureMC` samples bitwise-identical
//!    activity to the F2 uniform site Bernoulli law and to
//!    `StaticPercolationMC` under the same RNG stream.
//! 2. **Hand-graph exact expectations**: path, cycle, star, and complete
//!    graphs at `keep in {0.3, 0.5, 0.7}` versus exhaustive enumeration of
//!    all vertex subsets, with per-subset quantities computed by an
//!    independent BFS reference; multi-seed |z| < 4 and |z̄| < 2 gates on
//!    every degree array element and scalar observable.
//! 3. **Relabeling invariance**: identical-in-distribution arrays under a
//!    vertex permutation (statistical gate; RNG draw order changes with
//!    labeling, so bitwise equality is not the claim).
//! 4. **Empty degree bins**: absent degree values keep the full
//!    `[0 ..= max_degree]` shape with exact zeros and explicit occupancy
//!    markers.
//! 5. **Literature anchors**:
//!    Callaway, Newman, Strogatz, Watts, PRL 85, 5468 (2000)
//!    (arXiv:cond-mat/0007300) — giant fraction under uniform random
//!    failure against the generating-function solution of their Eq. (5)
//!    using the graph's own degree distribution; and Cohen, Erez,
//!    ben-Avraham, Havlin, PRL 85, 4626 (2000) (arXiv:cond-mat/0007048) —
//!    the critical removal fraction `f_c = 1 - 1/(kappa - 1)` brackets the
//!    observed transition. Both PDFs are archived under `docs/ref/`.

use std::collections::HashMap;

use carlo_rs::{Context, Estimate, MonteCarlo};
use cmc_rs::{
    build_chain, GraphView, NetworkProfile, ObservablePlan, RandomFailureMC, SiteBernoulli,
    StaticConfiguration, StaticLaw, StaticObservable, StaticPercolationMC, VertexActivity,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::zscore_seed_count;

/// One seeded failure run: (scalar estimates, array estimates).
fn run_failure(
    mc: &mut RandomFailureMC,
    seed: u64,
    sweeps: u64,
) -> (
    HashMap<String, Estimate>,
    HashMap<String, carlo_rs::ArrayEstimate>,
) {
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 16);
    for _ in 0..sweeps {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    ctx.finalize_measurements_split()
}

fn star_lattice(leaves: usize) -> cmc_rs::CsrLattice {
    cmc_rs::CsrLattice::from_edges(
        leaves + 1,
        (1..=leaves)
            .map(|leaf| cmc_rs::Bond::new(0, leaf, cmc_rs::BondType::Generic, 1.0))
            .collect(),
    )
}

fn cycle_lattice(n: usize) -> cmc_rs::CsrLattice {
    cmc_rs::CsrLattice::from_edges(
        n,
        (0..n)
            .map(|site| cmc_rs::Bond::new(site, (site + 1) % n, cmc_rs::BondType::Generic, 1.0))
            .collect(),
    )
}

fn complete_lattice(n: usize) -> cmc_rs::CsrLattice {
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            edges.push(cmc_rs::Bond::new(
                left,
                right,
                cmc_rs::BondType::Generic,
                1.0,
            ));
        }
    }
    cmc_rs::CsrLattice::from_edges(n, edges)
}

fn permutation_lattice(base: &cmc_rs::CsrLattice, permutation: &[usize]) -> cmc_rs::CsrLattice {
    let edges = base
        .edges
        .iter()
        .map(|bond| {
            cmc_rs::Bond::new(
                permutation[bond.source],
                permutation[bond.target],
                cmc_rs::BondType::Generic,
                1.0,
            )
        })
        .collect();
    cmc_rs::CsrLattice::from_edges(base.n_sites, edges)
}

// ─────────────────────────────────────────────────────────────────────────
// Independent enumeration reference (no analyzer code).
// ─────────────────────────────────────────────────────────────────────────

/// Per-subset quantities from an independent adjacency BFS.
struct SubsetObservables {
    giant_fraction: f64,
    susceptibility: Option<f64>,
    active_fraction_by_degree: Vec<f64>,
    giant_membership_by_degree: Vec<f64>,
    first_moment: Vec<f64>,
    second_moment: Vec<f64>,
}

/// BFS components over the active subset `active` of `adjacency`;
/// returns `(size per vertex, canonical id per component)`.
fn bfs_components(adjacency: &[Vec<usize>], active: &[bool]) -> (Vec<usize>, Vec<usize>) {
    let n = adjacency.len();
    let mut size = vec![0usize; n];
    let mut canonical = vec![usize::MAX; n];
    let mut visited = vec![false; n];
    for seed in 0..n {
        if !active[seed] || visited[seed] {
            continue;
        }
        let mut component = vec![seed];
        visited[seed] = true;
        let mut head = 0;
        while head < component.len() {
            let vertex = component[head];
            head += 1;
            for &neighbor in &adjacency[vertex] {
                if active[neighbor] && !visited[neighbor] {
                    visited[neighbor] = true;
                    component.push(neighbor);
                }
            }
        }
        let low = *component.iter().min().expect("non-empty component");
        for &vertex in &component {
            size[vertex] = component.len();
            canonical[vertex] = low;
        }
    }
    (size, canonical)
}

fn subset_observables(
    adjacency: &[Vec<usize>],
    degrees: &[usize],
    bins: usize,
    active: &[bool],
) -> SubsetObservables {
    let n = adjacency.len();
    let (size, canonical) = bfs_components(adjacency, active);

    // Largest component, F4 tie policy: equal sizes resolve to the lowest
    // canonical identity.
    let mut largest_size = 0usize;
    let mut largest_id = usize::MAX;
    let mut component_sizes: HashMap<usize, usize> = HashMap::new();
    let mut raw_second_moment: u128 = 0;
    for vertex in 0..n {
        if !active[vertex] {
            continue;
        }
        let entry = component_sizes.entry(canonical[vertex]).or_insert(0usize);
        *entry += 1;
    }
    for (&id, &component_size) in &component_sizes {
        raw_second_moment += (component_size as u128).pow(2);
        if component_size > largest_size || (component_size == largest_size && id < largest_id) {
            largest_size = component_size;
            largest_id = id;
        }
    }

    let mut count_per_degree = vec![0usize; bins];
    for &degree in degrees {
        count_per_degree[degree] += 1;
    }
    let mut active_count = vec![0u64; bins];
    let mut giant_count = vec![0u64; bins];
    let mut first_moment = vec![0u64; bins];
    let mut second_moment = vec![0u128; bins];
    for vertex in 0..n {
        if !active[vertex] {
            continue;
        }
        let bin = degrees[vertex];
        active_count[bin] += 1;
        first_moment[bin] += size[vertex] as u64;
        second_moment[bin] += (size[vertex] as u128).pow(2);
        if canonical[vertex] == largest_id {
            giant_count[bin] += 1;
        }
    }
    let ratio = |counts: &Vec<u64>| -> Vec<f64> {
        (0..bins)
            .map(|bin| {
                if count_per_degree[bin] > 0 {
                    counts[bin] as f64 / count_per_degree[bin] as f64
                } else {
                    0.0
                }
            })
            .collect()
    };

    let active_total: usize = active.iter().filter(|&&a| a).count();
    let finite = active_total.saturating_sub(largest_size);
    let susceptibility = if finite > 0 {
        let largest_squared = (largest_size as u128).pow(2);
        Some((raw_second_moment - largest_squared) as f64 / finite as f64)
    } else {
        None
    };

    SubsetObservables {
        giant_fraction: largest_size as f64 / n as f64,
        susceptibility,
        active_fraction_by_degree: ratio(&active_count),
        giant_membership_by_degree: ratio(&giant_count),
        first_moment: ratio(&first_moment),
        second_moment: second_moment
            .iter()
            .enumerate()
            .map(|(bin, value)| {
                if count_per_degree[bin] > 0 {
                    *value as f64 / count_per_degree[bin] as f64
                } else {
                    0.0
                }
            })
            .collect(),
    }
}

/// Exact ensemble expectations by enumerating all `2^V` subsets.
struct ExactExpectations {
    giant_fraction: f64,
    /// Conditional mean susceptibility over defined samples.
    susceptibility: f64,
    active_fraction_by_degree: Vec<f64>,
    giant_membership_by_degree: Vec<f64>,
    first_moment: Vec<f64>,
    second_moment: Vec<f64>,
}

fn exact_expectations(
    adjacency: &[Vec<usize>],
    degrees: &[usize],
    bins: usize,
    keep: f64,
) -> ExactExpectations {
    let n = adjacency.len();
    let mut total_weight = 0.0;
    let mut defined_weight = 0.0;
    let mut giant = 0.0;
    let mut susceptibility = 0.0;
    let mut active_fraction = vec![0.0; bins];
    let mut giant_membership = vec![0.0; bins];
    let mut first = vec![0.0; bins];
    let mut second = vec![0.0; bins];
    for mask in 0..(1usize << n) {
        let active: Vec<bool> = (0..n).map(|vertex| mask & (1 << vertex) != 0).collect();
        let occupied = active.iter().filter(|&&a| a).count();
        let weight = keep.powi(occupied as i32) * (1.0 - keep).powi((n - occupied) as i32);
        let observables = subset_observables(adjacency, degrees, bins, &active);
        total_weight += weight;
        giant += weight * observables.giant_fraction;
        if let Some(chi) = observables.susceptibility {
            defined_weight += weight;
            susceptibility += weight * chi;
        }
        for bin in 0..bins {
            active_fraction[bin] += weight * observables.active_fraction_by_degree[bin];
            giant_membership[bin] += weight * observables.giant_membership_by_degree[bin];
            first[bin] += weight * observables.first_moment[bin];
            second[bin] += weight * observables.second_moment[bin];
        }
    }
    assert!(
        (total_weight - 1.0).abs() < 1e-12,
        "subset weights must sum to 1"
    );
    // K4 never admits a defined sample: every non-empty active subset is one
    // complete component, so V_active == S_max always. A zero defined weight
    // is legitimate; the Monte Carlo gate then skips via empty chi_means.
    if defined_weight > 0.0 {
        susceptibility /= defined_weight;
    }
    for array in [
        &mut active_fraction,
        &mut giant_membership,
        &mut first,
        &mut second,
    ] {
        for value in array.iter_mut() {
            *value /= total_weight;
        }
    }
    giant /= total_weight;
    ExactExpectations {
        giant_fraction: giant,
        susceptibility,
        active_fraction_by_degree: active_fraction,
        giant_membership_by_degree: giant_membership,
        first_moment: first,
        second_moment: second,
    }
}

fn adjacency_of(lattice: &cmc_rs::CsrLattice) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); lattice.n_sites];
    for bond in &lattice.edges {
        adjacency[bond.source].push(bond.target);
        adjacency[bond.target].push(bond.source);
    }
    adjacency
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: random failure reduces to the site law.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn random_failure_activity_is_bitwise_site_law() {
    let lattice = star_lattice(7);
    let plan = ObservablePlan::all(Vec::new());
    let mut facade = RandomFailureMC::new(lattice.clone(), 0.61, plan).expect("facade constructs");
    let mut adapter = StaticPercolationMC::new(
        lattice.clone(),
        StaticLaw::site(0.61).expect("probability"),
        ObservablePlan::new([StaticObservable::LargestSize], Vec::new()),
    )
    .expect("adapter constructs");

    let mut facade_ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 4);
    let mut adapter_ctx =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 4);
    for sweep in 0..32 {
        facade.sweep(&mut facade_ctx);
        adapter.sweep(&mut adapter_ctx);
        for index in 0..lattice.n_sites {
            let vertex = lattice.vertex_id(index).expect("dense index");
            assert_eq!(
                facade.configuration().vertex_active(vertex),
                adapter.configuration().vertex_active(vertex),
                "vertex {index} diverged at sweep {sweep}"
            );
        }
    }

    // Also match a raw F2 law sampling on an identically-seeded stream.
    let mut law_ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 4);
    let mut facade_ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 4);
    let mut configuration = StaticConfiguration::new(lattice.n_sites, lattice.n_edges());
    let mut fresh = RandomFailureMC::new(lattice, 0.61, ObservablePlan::all(Vec::new()))
        .expect("facade constructs");
    for sweep in 0..32 {
        SiteBernoulli::new(cmc_rs::ProbabilityField::Uniform(
            cmc_rs::Probability::new(0.61).expect("probability"),
        ))
        .sample(fresh.lattice(), &mut configuration, &mut law_ctx.rng)
        .expect("domains match");
        fresh.sweep(&mut facade_ctx);
        for index in 0..fresh.lattice().n_sites {
            let vertex = fresh.lattice().vertex_id(index).expect("dense index");
            assert_eq!(
                configuration.vertex_active(vertex),
                fresh.configuration().vertex_active(vertex),
                "raw-law divergence at sweep {sweep}, vertex {index}"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: hand graphs versus exhaustive enumeration.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn hand_graphs_match_exhaustive_enumeration() {
    struct Case {
        name: &'static str,
        lattice: cmc_rs::CsrLattice,
    }
    let cases = [
        Case {
            name: "path5",
            lattice: build_chain(5, false),
        },
        Case {
            name: "cycle6",
            lattice: cycle_lattice(6),
        },
        Case {
            name: "star6",
            lattice: star_lattice(5),
        },
        Case {
            name: "complete4",
            lattice: complete_lattice(4),
        },
    ];

    let n_seeds = zscore_seed_count(16);
    let sweeps = 600;
    let mut failures = Vec::new();

    for case in &cases {
        let adjacency = adjacency_of(&case.lattice);
        let profile = NetworkProfile::from_graph(&case.lattice).expect("profile");
        let bins = profile.bins().len();
        let degrees: Vec<usize> = profile
            .degrees()
            .iter()
            .map(|&degree| degree as usize)
            .collect();

        for &keep in &[0.3_f64, 0.5, 0.7] {
            let exact = exact_expectations(&adjacency, &degrees, bins, keep);
            let mut z_max = 0.0_f64;
            let mut z_sum = 0.0_f64;
            let mut count = 0usize;

            let mut check = |name: String, mean: f64, stderr: f64, expected: f64| {
                let stderr = stderr.max(1e-12);
                let z = (mean - expected) / stderr;
                z_max = z_max.max(z.abs());
                z_sum += z;
                count += 1;
                if z.abs() >= 4.0 {
                    failures.push(format!(
                        "{}/keep={keep} {name}: mean {mean:.5} vs exact {expected:.5}, |z|={:.2}",
                        case.name,
                        z.abs()
                    ));
                }
            };

            // Pool all seeds in one pass for exact per-observable gates.
            let mut giant_means = Vec::new();
            let mut chi_means = Vec::new();
            let mut active_means = vec![Vec::new(); bins];
            let mut giant_deg_means = vec![Vec::new(); bins];
            let mut first_means = vec![Vec::new(); bins];
            let mut second_means = vec![Vec::new(); bins];
            for seed in 0..n_seeds as u64 {
                let mut mc = RandomFailureMC::new(
                    case.lattice.clone(),
                    keep,
                    ObservablePlan::all(Vec::new()),
                )
                .expect("facade constructs");
                let (scalars, arrays) = run_failure(&mut mc, 0x4e32 + 7919 * seed, sweeps);
                giant_means.push((
                    scalars["GiantFraction"].mean,
                    scalars["GiantFraction"].stderr,
                ));
                // A seed where every sample is undefined records no
                // susceptibility at all (missing, never NaN); conditioning
                // on defined samples mirrors the enumeration reference.
                if let Some(chi) = scalars.get("FiniteClusterSusceptibility") {
                    chi_means.push((chi.mean, chi.stderr));
                }
                let read = |arrays: &HashMap<String, carlo_rs::ArrayEstimate>,
                            name: &str,
                            bin: usize|
                 -> (f64, f64) {
                    let element = arrays[name].element(bin).expect("bin");
                    (element.mean, element.stderr)
                };
                for bin in 0..bins {
                    active_means[bin].push(read(&arrays, "ActiveFractionByDegree", bin));
                    giant_deg_means[bin].push(read(&arrays, "GiantMembershipByDegree", bin));
                    first_means[bin].push(read(&arrays, "DegreeFirstMoment", bin));
                    second_means[bin].push(read(&arrays, "DegreeSecondMoment", bin));
                }
            }

            let pooled = |estimates: &[(f64, f64)]| -> (f64, f64) {
                let n = estimates.len() as f64;
                let mean = estimates.iter().map(|&(m, _)| m).sum::<f64>() / n;
                // Independent seeds: combine stderrs in quadrature.
                let variance = estimates.iter().map(|&(_, se)| se * se).sum::<f64>() / (n * n);
                (mean, variance.sqrt())
            };

            let (giant_mean, giant_stderr) = pooled(&giant_means);
            check(
                "GiantFraction".into(),
                giant_mean,
                giant_stderr,
                exact.giant_fraction,
            );
            // Skip the susceptibility gate only when the process almost
            // never produces finite clusters (e.g. K4 at keep = 0.3: an
            // all-or-one cluster dominates); the defined-indicator
            // convention is then the observable under test.
            if !chi_means.is_empty() {
                let (chi_mean, chi_stderr) = pooled(&chi_means);
                check(
                    "FiniteClusterSusceptibility".into(),
                    chi_mean,
                    chi_stderr,
                    exact.susceptibility,
                );
            } else {
                eprintln!(
                    "[N2 {}/keep={keep}] susceptibility undefined in every sample; gate skipped",
                    case.name
                );
            }
            for bin in 0..bins {
                let (mean, stderr) = pooled(&active_means[bin]);
                check(
                    format!("ActiveFractionByDegree[{bin}]"),
                    mean,
                    stderr,
                    exact.active_fraction_by_degree[bin],
                );
                let (mean, stderr) = pooled(&giant_deg_means[bin]);
                check(
                    format!("GiantMembershipByDegree[{bin}]"),
                    mean,
                    stderr,
                    exact.giant_membership_by_degree[bin],
                );
                let (mean, stderr) = pooled(&first_means[bin]);
                check(
                    format!("DegreeFirstMoment[{bin}]"),
                    mean,
                    stderr,
                    exact.first_moment[bin],
                );
                let (mean, stderr) = pooled(&second_means[bin]);
                check(
                    format!("DegreeSecondMoment[{bin}]"),
                    mean,
                    stderr,
                    exact.second_moment[bin],
                );
            }

            let mean_z = z_sum / count as f64;
            eprintln!(
                "[N2 {}/keep={keep}] max|z|={z_max:.2} mean z={mean_z:+.2} ({count} observables)",
                case.name
            );
            assert!(
                mean_z.abs() < 2.0,
                "{}/keep={keep}: pooled mean |z| = {mean_z:.2} exceeds 2",
                case.name
            );
        }
    }
    assert!(
        failures.is_empty(),
        "enumeration gates failed:\n{}",
        failures.join("\n")
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: relabeling invariance.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn degree_observables_are_relabeling_invariant() {
    let base = star_lattice(9);
    // Reverse labeling: vertex i becomes 9 - i.
    let permutation: Vec<usize> = (0..base.n_sites).rev().collect();
    let relabeled = permutation_lattice(&base, &permutation);
    assert_eq!(base.n_sites, relabeled.n_sites);
    assert_eq!(base.n_edges(), relabeled.n_edges());

    let sweeps = 2000;
    let mut identity =
        RandomFailureMC::new(base, 0.5, ObservablePlan::all(Vec::new())).expect("identity facade");
    let mut shuffled = RandomFailureMC::new(relabeled, 0.5, ObservablePlan::all(Vec::new()))
        .expect("relabeled facade");
    let (scalars_a, arrays_a) = run_failure(&mut identity, 11, sweeps);
    let (scalars_b, arrays_b) = run_failure(&mut shuffled, 12, sweeps);

    // Constant arrays must be identical: the degree multiset is a labeling
    // invariant, so bin shapes and occupancies agree exactly.
    for name in ["DegreeVertexCount", "DegreeOccupied"] {
        let left = arrays_a[name].elements();
        let right = arrays_b[name].elements();
        for (l, r) in left.iter().zip(right) {
            assert_eq!(l.mean, r.mean, "{name} differs under relabeling");
            assert_eq!(l.stderr, r.stderr, "{name} stderr differs");
        }
    }

    let compare = |name: &str,
                   left: &HashMap<String, carlo_rs::ArrayEstimate>,
                   right: &HashMap<String, carlo_rs::ArrayEstimate>| {
        for (l, r) in left[name].elements().iter().zip(right[name].elements()) {
            let sigma = (l.stderr * l.stderr + r.stderr * r.stderr)
                .sqrt()
                .max(1e-12);
            let z = (l.mean - r.mean) / sigma;
            assert!(
                z.abs() < 4.0,
                "{name} failed relabeling invariance: {z:.2} sigma"
            );
        }
    };
    for name in [
        "ActiveFractionByDegree",
        "GiantMembershipByDegree",
        "DegreeFirstMoment",
        "DegreeSecondMoment",
    ] {
        compare(name, &arrays_a, &arrays_b);
    }

    let giant_a = &scalars_a["GiantFraction"];
    let giant_b = &scalars_b["GiantFraction"];
    let sigma = (giant_a.stderr * giant_a.stderr + giant_b.stderr * giant_b.stderr).sqrt();
    assert!(
        ((giant_a.mean - giant_b.mean) / sigma.max(1e-12)).abs() < 4.0,
        "GiantFraction failed relabeling invariance"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: empty degree bins.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn empty_degree_bins_have_explicit_semantics() {
    // Star with 5 vertices: degrees {4, 1, 1, 1, 1}; bins 0, 2, 3 are empty.
    let lattice = star_lattice(4);
    let profile = NetworkProfile::from_graph(&lattice).expect("profile");
    assert_eq!(profile.bins(), &[0, 4, 0, 0, 1]);
    assert_eq!(profile.max_degree(), 4);
    assert_eq!(profile.vertices_with_degree(1), 4);
    assert_eq!(profile.vertices_with_degree(0), 0);

    let mut mc =
        RandomFailureMC::new(lattice, 0.5, ObservablePlan::all(Vec::new())).expect("facade");
    let (_, arrays) = run_failure(&mut mc, 3, 64);

    for name in [
        "DegreeVertexCount",
        "DegreeOccupied",
        "ActiveFractionByDegree",
        "GiantMembershipByDegree",
        "DegreeFirstMoment",
        "DegreeSecondMoment",
    ] {
        let array = &arrays[name];
        assert_eq!(array.shape(), &[5], "{name} must keep the full bin shape");
        for empty_bin in [0usize, 2, 3] {
            let element = array.element(empty_bin).expect("empty bin exists");
            assert_eq!(element.mean, 0.0, "{name}[{empty_bin}] must be exactly 0");
            assert_eq!(element.stderr, 0.0, "{name}[{empty_bin}] has no variance");
        }
    }
    assert_eq!(arrays["DegreeVertexCount"].element(1).unwrap().mean, 4.0);
    assert_eq!(arrays["DegreeVertexCount"].element(4).unwrap().mean, 1.0);
    assert_eq!(arrays["DegreeOccupied"].element(1).unwrap().mean, 1.0);
    assert_eq!(arrays["DegreeOccupied"].element(2).unwrap().mean, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5: literature anchors (Callaway et al. 2000; Cohen et al. 2000).
// ─────────────────────────────────────────────────────────────────────────

/// Sparse G(n, m) Erdos-Renyi graph with distinct non-loop edges.
fn erdos_renyi(n: usize, edge_count: usize, seed: u64) -> cmc_rs::CsrLattice {
    use rand::RngExt;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut pairs = Vec::with_capacity(edge_count);
    while pairs.len() < edge_count {
        let left = rng.random_range(0..n);
        let right = rng.random_range(0..n);
        if left != right {
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }
        }
    }
    let edges = pairs
        .into_iter()
        .map(|(left, right)| cmc_rs::Bond::new(left, right, cmc_rs::BondType::Generic, 1.0))
        .collect();
    cmc_rs::CsrLattice::from_edges(n, edges)
}

/// Callaway et al. Eq. (5) giant fraction for uniform keep `q` on a
/// configuration-model graph with the empirical degree distribution `bins`.
///
/// Solves `u = 1 - q + q G1(u)` by fixed-point iteration and returns
/// `S = q (1 - G0(u))`, plus the critical keep probability
/// `q_c = 1 / G1'(1)` implied by the same distribution.
fn callaway_giant_fraction(bins: &[usize], q: f64) -> (f64, f64) {
    let total: f64 = bins.iter().sum::<usize>() as f64;
    let g0 = |x: f64| -> f64 {
        bins.iter()
            .enumerate()
            .map(|(k, &count)| count as f64 / total * x.powi(k as i32))
            .sum()
    };
    // G1(x) = sum_k (k p_k / z) x^{k-1}; the k = 0 term drops out.
    let g1 = |x: f64| -> f64 {
        let z: f64 = bins
            .iter()
            .enumerate()
            .map(|(k, &count)| k as f64 * count as f64)
            .sum::<f64>()
            / total;
        bins.iter()
            .enumerate()
            .skip(1)
            .map(|(k, &count)| k as f64 * count as f64 / total / z * x.powi(k as i32 - 1))
            .sum()
    };
    // G1'(1) = (sum k(k-1) p_k) / z.
    let sum_kk: f64 = bins
        .iter()
        .enumerate()
        .map(|(k, &count)| k as f64 * (k as f64 - 1.0) * count as f64)
        .sum::<f64>()
        / total;
    let z: f64 = bins
        .iter()
        .enumerate()
        .map(|(k, &count)| k as f64 * count as f64)
        .sum::<f64>()
        / total;
    let q_c = z / sum_kk;

    let mut u = 1.0 - q;
    for _ in 0..200_000 {
        let next = 1.0 - q + q * g1(u);
        if (next - u).abs() < 1e-15 {
            u = next;
            break;
        }
        u = next;
    }
    (q * (1.0 - g0(u)), q_c)
}

#[test]
fn random_failure_matches_callaway_cohen_predictions() {
    // ER configuration graph, mean degree z = 2.5 -> ideal q_c = 0.4
    // (Cohen et al.: f_c = 1 - 1/(kappa - 1) with kappa = <k^2>/<k>).
    const N: usize = 100_000;
    const Z: f64 = 2.5;
    const SEEDS: u64 = 8;
    const SWEEPS: u64 = 250;

    let lattice = erdos_renyi(N, (Z * N as f64 / 2.0) as usize, 0x05ee_d0e2);
    let profile = NetworkProfile::from_graph(&lattice).expect("profile");
    let (_, kappa) = profile.degree_moments().expect("non-empty graph");
    // Cohen et al. critical removal fraction from the empirical moments.
    let q_cohen = 1.0 / (kappa - 1.0);
    let (mean_degree, _) = profile.degree_moments().expect("non-empty graph");
    eprintln!(
        "[N2 literature] ER N={N} <k>={mean_degree:.4} kappa={kappa:.4} \
         -> Cohen q_c = 1/(kappa-1) = {q_cohen:.4} (ideal Poisson: 0.4)"
    );
    assert!(
        (q_cohen - 0.4).abs() < 0.01,
        "empirical kappa must stay close to the Poisson ideal, got q_c = {q_cohen}"
    );

    let giant_fraction = |keep: f64| -> (f64, f64) {
        let mut means = Vec::new();
        let mut stderrs = Vec::new();
        for seed in 0..SEEDS {
            let mut mc = RandomFailureMC::new(
                lattice.clone(),
                keep,
                ObservablePlan::new([StaticObservable::GiantFraction], Vec::new()),
            )
            .expect("facade constructs");
            let (scalars, _) = run_failure(&mut mc, 0x1_a + 61 * seed, SWEEPS);
            let estimate = &scalars["GiantFraction"];
            means.push(estimate.mean);
            stderrs.push(estimate.stderr);
        }
        let n = means.len() as f64;
        let mean = means.iter().sum::<f64>() / n;
        let stderr = (stderrs.iter().map(|se| se * se).sum::<f64>() / (n * n)).sqrt();
        (mean, stderr)
    };

    // (a) Above threshold: measured S agrees with the Callaway et al.
    //     generating-function solution built from the graph's own degree
    //     distribution. Slack 0.006 absorbs finite-N corrections; the
    //     statistical gate is 4 combined sigma.
    for &keep in &[0.5_f64, 0.65, 0.8] {
        let (theory, _) = callaway_giant_fraction(profile.bins(), keep);
        let (mean, stderr) = giant_fraction(keep);
        let slack = 0.006;
        let bound = 4.0 * stderr + slack;
        eprintln!(
            "[N2 literature] keep={keep}: S_mc={mean:.4} ± {stderr:.4} vs \
             Callaway S={theory:.4} (bound {bound:.4})"
        );
        assert!(
            (mean - theory).abs() < bound,
            "keep={keep}: |S_mc - S_theory| = {:.4} exceeds {bound:.4}",
            (mean - theory).abs()
        );
    }

    // (b) Cohen et al. threshold bracket: below q_c the giant is
    //     statistically extinguished; above it a giant exists and follows
    //     the Callaway et al. curve. Both points sit outside the finite-size
    //     rounding window (~N^{-1/3} ~ 0.02 in q at N = 1e5): the below
    //     point at q_c - 0.03, the above point at q_c + 0.10 where mean-field
    //     theory applies cleanly. The prediction is a statistical bracket,
    //     never an exact threshold claim.
    let below = q_cohen - 0.03;
    let (below_mean, below_stderr) = giant_fraction(below);
    eprintln!("[N2 literature] keep=q_c-0.03={below:.4}: S={below_mean:.4} ± {below_stderr:.4}");
    assert!(
        below_mean < 4.0 * below_stderr + 0.004,
        "below the Cohen threshold the giant must vanish: S = {below_mean:.4} ± {below_stderr:.4}"
    );
    let above = q_cohen + 0.10;
    let (above_theory, _) = callaway_giant_fraction(profile.bins(), above);
    let (above_mean, above_stderr) = giant_fraction(above);
    eprintln!(
        "[N2 literature] keep=q_c+0.10={above:.4}: S={above_mean:.4} ± {above_stderr:.4} vs \
         Callaway S={above_theory:.4}"
    );
    assert!(
        above_mean > 0.1,
        "above the Cohen threshold a giant must exist: S = {above_mean:.4}"
    );
    assert!(
        (above_mean - above_theory).abs() < 4.0 * above_stderr + 0.006,
        "above-threshold giant fraction disagrees with Callaway: |{:.4}| exceeds bound",
        (above_mean - above_theory).abs()
    );
}

#[test]
fn degree_moment_gate_rejects_beyond_exact_f64_bound() {
    // Review follow-up (2026-10-22): V^3 must stay within the exact f64
    // integer range because degree second moments sum s(v)^2 over vertices
    // (sum_c s_c^3 <= V^3). 208063^3 = 9,007,091,372,906,047 <= 2^53 <
    // 208064^3 = 9,007,221,244,166,144, so exactly V = 208064 is the first
    // rejected size. Edgeless lattices make the construction cheap; the gate
    // fires before any workspace is prepared.
    let plan = ObservablePlan::all(Vec::new());

    let boundary = cmc_rs::CsrLattice::from_edges(208_063, Vec::new());
    assert!(RandomFailureMC::new(boundary, 0.5, plan.clone()).is_ok());

    let rejected = cmc_rs::CsrLattice::from_edges(208_064, Vec::new());
    match RandomFailureMC::new(rejected, 0.5, plan) {
        Err(cmc_rs::NetworkPercolationError::InexactDegreeMoment { total_vertices }) => {
            assert_eq!(total_vertices, 208_064)
        }
        _ => panic!("V = 208064 must be rejected by the exactness gate"),
    }
}
