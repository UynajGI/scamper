//! L2 chemical-distance validation (roadmap v2 §4, stage L2).
//!
//! Gates, verbatim from the stage spec (regime semantics corrected:
//! inside (1,2) the TYPICAL distance is the sharp main theorem; at
//! s <= 1 only the DIAMETER envelope is proved; s = 2 is explicitly not
//! understood and gets NO gate — see the module documentation of
//! `percolation::long_range_distances`):
//!
//! 1. **BFS correctness** — per-realization distance profiles vs a
//!    test-local independently written BFS on small boxes across seeds
//!    and s values, exact per-vertex equality including the INFINITY
//!    sentinels; the exact diameter vs a test-local all-pairs reference;
//!    occupancy parity with the L1 core's `occupied_pairs()`; estimator
//!    mean/SEM/reachable-fraction vs manual statistics (the reachable
//!    conditioning exercised at beta < 1);
//! 2. **Bitwise** — same seed, identical profiles and diameter bitwise;
//!    the draw stream equals the L1 adapter's (same occupied graph);
//! 3. **The (1,2) regime** — typical distance between pairs at
//!    |x - y| = r (center source, both directions pooled) at beta = 1:
//!    strictly monotone growth, the empirical D-vs-r log-log slope CI
//!    bounded below 0.5 at these box sizes (a finite-window statement,
//!    not a proof of sub-polynomial growth), the fitted polylog exponent
//!    ordered in s exactly as the paper's Delta, approaching it from
//!    below in the r-window (the honest finite-size direction — the
//!    asymptotic Delta itself is NOT reached at these box sizes; every
//!    fitted slope is reported against the paper value);
//! 4. **s > 2 linear** — at s = 2.5 the finite-window log-log slope sits
//!    at or below 1 (the O(1/r) deficit), so the gate is PAIRED across
//!    realizations: per-realization slopes and `D/r = A + c/r` fits,
//!    then 95% CIs requiring the slope to exclude the polylog regime,
//!    the deficit c to be positive, and the extrapolated time constant A
//!    to lie inside (0, 1), with a non-increasing measured ratio;
//! 5. **s <= 1 diameter envelopes** — s = 0.5: the box diameter is
//!    bracketed as a bounded constant against Eq. (1.8)'s
//!    `d/(d-s)` scale (no power growth); s = 1: Eq. (1.9)'s upper
//!    envelope only (diameter at most logarithmic — NO lower-bound
//!    claim);
//! 6. **s = 2 not bracketed** — no gate, no run (documented);
//! 7. **Throughput report** — realizations/second at two sizes.

use std::time::Instant;

use cmc_rs::{
    chemical_distance_run, long_range_run, ChemicalDistanceError, ChemicalDistanceEstimator,
    ChemicalDistanceSampler, ChemicalDistanceWorkspace, DiameterTracking, LongRangeBox,
    LongRangePercolationMC, LongRangeWorkspace, DIAMETER_VERTEX_CAP,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ─────────────────────────────────────────────────────────────────────────
// Test-local references and helpers
// ─────────────────────────────────────────────────────────────────────────

/// Draw one uniform field per the L1 contract: one `f64` per pair in
/// ascending pair order.
fn draw_uniforms(pairs: usize, rng: &mut Xoshiro256PlusPlus) -> Vec<f64> {
    (0..pairs).map(|_| rng.random::<f64>()).collect()
}

/// Test-local occupancy replay: the occupied endpoint pairs of a uniform
/// field under the law, decided pair-by-pair through the public geometry
/// API (independent of the production scan structure).
fn reference_occupied(box_: &LongRangeBox, beta: f64, s: f64, uniforms: &[f64]) -> Vec<usize> {
    let mut occupied = Vec::new();
    for (pair, &uniform) in uniforms.iter().enumerate() {
        let (first, second) = box_.pair_endpoints(pair).expect("dense pair");
        let probability = LongRangeBox::pair_probability(box_.distance(first, second), beta, s);
        if uniform < probability {
            occupied.push(pair);
        }
    }
    occupied
}

/// Test-local multi-source BFS over an explicit adjacency (usize::MAX =
/// unreachable).
fn reference_profile(vertices: usize, adjacency: &[Vec<usize>], sources: &[usize]) -> Vec<usize> {
    let mut distance = vec![usize::MAX; vertices];
    let mut queue = Vec::with_capacity(vertices);
    for &source in sources {
        if distance[source] != 0 {
            distance[source] = 0;
            queue.push(source);
        }
    }
    let mut head = 0;
    while head < queue.len() {
        let vertex = queue[head];
        head += 1;
        for &neighbor in &adjacency[vertex] {
            if distance[neighbor] == usize::MAX {
                distance[neighbor] = distance[vertex] + 1;
                queue.push(neighbor);
            }
        }
    }
    distance
}

/// Test-local exact diameter: eccentricity of every vertex, maximum
/// finite distance over all pairs including the diagonal (the production
/// definition, independently implemented).
fn reference_diameter(vertices: usize, adjacency: &[Vec<usize>]) -> usize {
    let mut diameter = 0;
    for source in 0..vertices {
        let distance = reference_profile(vertices, adjacency, &[source]);
        let eccentricity = distance
            .iter()
            .copied()
            .filter(|&d| d != usize::MAX)
            .max()
            .unwrap_or(0);
        diameter = diameter.max(eccentricity);
    }
    diameter
}

/// The occupied pairs of a production run, recovered from the workspace
/// CSR as sorted pair indices (the L1 consideration order).
fn workspace_occupied_pairs(
    box_: &LongRangeBox,
    workspace: &ChemicalDistanceWorkspace,
) -> Vec<usize> {
    let offsets = workspace.adjacency_offsets();
    let neighbors = workspace.adjacency_neighbors();
    let mut pairs: Vec<usize> = Vec::with_capacity(workspace.occupied_edge_count());
    for vertex in 0..offsets.len() - 1 {
        for &neighbor in &neighbors[offsets[vertex]..offsets[vertex + 1]] {
            if neighbor > vertex {
                pairs.push(
                    box_.pair_index(vertex, neighbor)
                        .expect("distinct vertices"),
                );
            }
        }
    }
    pairs.sort_unstable();
    pairs
}

/// Ordinary least squares of `y ~ x`: (slope, intercept, slope SE, R^2).
fn ols(points: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|point| point.0).sum();
    let sy: f64 = points.iter().map(|point| point.1).sum();
    let sxx: f64 = points.iter().map(|point| point.0 * point.0).sum();
    let sxy: f64 = points.iter().map(|point| point.0 * point.1).sum();
    let syy: f64 = points.iter().map(|point| point.1 * point.1).sum();
    let denominator = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denominator;
    let intercept = (sy - slope * sx) / n;
    let r_squared = ((n * sxy - sx * sy).powi(2)) / (denominator * (n * syy - sy * sy));
    let residuals: f64 = points
        .iter()
        .map(|point| {
            let error = point.1 - (intercept + slope * point.0);
            error * error
        })
        .sum();
    let slope_se = (residuals / (n - 2.0) / (sxx - sx * sx / n)).sqrt();
    (slope, intercept, slope_se, r_squared)
}

/// Biskup's Delta for d = 1: Delta^-1 = log2(2/s).
fn biskup_delta(s: f64) -> f64 {
    std::f64::consts::LN_2 / (2.0 / s).ln()
}

// ─────────────────────────────────────────────────────────────────────────
// Pre-run workspace state
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn pre_run_workspace_reports_empty_adjacency() {
    // Review follow-up (M2): `prepare` used to resize the neighbor buffer
    // to `2 * pairs`, so `occupied_edge_count()` reported the pair count
    // before any run. Reserve capacity only; the accessors must report
    // an empty adjacency and an all-INFINITY profile.
    let box_ = LongRangeBox::new(&[5]).expect("box");
    let mut workspace = ChemicalDistanceWorkspace::new();
    workspace.prepare(5, 10, 4).expect("prepare");
    assert_eq!(workspace.occupied_edge_count(), 0);
    assert_eq!(workspace.adjacency_offsets(), &[0; 6]);
    assert!(workspace.adjacency_neighbors().is_empty());
    assert_eq!(workspace.profile().len(), 5);
    assert!(workspace
        .profile()
        .iter()
        .all(|value| *value == f64::INFINITY));

    let sampler = ChemicalDistanceSampler::new(box_, 1.0, 1.5, &[0]).expect("sampler");
    assert_eq!(sampler.workspace().occupied_edge_count(), 0);
    assert!(sampler.workspace().adjacency_neighbors().is_empty());
    assert!(sampler
        .profile()
        .iter()
        .all(|value| *value == f64::INFINITY));
    eprintln!("L2 pre-run workspace: empty adjacency, zero occupied edges, INFINITY profile");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: BFS correctness vs independent references
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn profiles_and_diameter_match_independent_references() {
    struct Config {
        sides: Vec<usize>,
        beta: f64,
        s: f64,
        source_sets: Vec<Vec<usize>>,
    }
    // (0.3, 0.8) sits below saturation: disconnected configurations and
    // INFINITY sentinels occur; the others are beta >= 1 (always
    // connected in d = 1). The 5x5 box exercises the 2D distance.
    let configs = [
        Config {
            sides: vec![16],
            beta: 0.3,
            s: 0.8,
            source_sets: vec![vec![8], vec![0, 15]],
        },
        Config {
            sides: vec![24],
            beta: 1.0,
            s: 1.5,
            source_sets: vec![vec![12], vec![0], vec![0, 23]],
        },
        Config {
            sides: vec![5, 5],
            beta: 2.2,
            s: 0.6,
            source_sets: vec![vec![12], vec![0, 24]],
        },
    ];
    let mut profile_cases = 0_usize;
    let mut diameter_cases = 0_usize;
    for config in &configs {
        let box_ = LongRangeBox::new(&config.sides).expect("box");
        let vertices = box_.vertices();
        for seed in 0..5_u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004D_1101 + seed);
            let uniforms = draw_uniforms(box_.pairs(), &mut rng);

            // Occupancy parity with the L1 core through its public
            // occupied-pair exposure.
            let mut l1_workspace = LongRangeWorkspace::new();
            long_range_run(&box_, config.beta, config.s, &uniforms, &mut l1_workspace)
                .expect("L1 parity run");
            let expected_occupied = l1_workspace.occupied_pairs().to_vec();
            assert_eq!(
                expected_occupied,
                reference_occupied(&box_, config.beta, config.s, &uniforms),
                "the L1 core's occupied set must equal the test-local law replay \
                 (sides {:?}, seed {seed})",
                config.sides
            );

            // Independent adjacency from the pair-universe replay.
            let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
            for &pair in &expected_occupied {
                let (first, second) = box_.pair_endpoints(pair).expect("dense pair");
                adjacency[first].push(second);
                adjacency[second].push(first);
            }
            let expected_diameter = reference_diameter(vertices, &adjacency);

            for sources in &config.source_sets {
                let mut workspace = ChemicalDistanceWorkspace::new();
                let outcome = chemical_distance_run(
                    &box_,
                    config.beta,
                    config.s,
                    &uniforms,
                    sources,
                    DiameterTracking::AllVertexBfs,
                    &mut workspace,
                )
                .expect("production run");
                assert_eq!(
                    workspace_occupied_pairs(&box_, &workspace),
                    expected_occupied,
                    "occupied-set parity with the L1 core (sides {:?}, seed {seed})",
                    config.sides
                );
                assert_eq!(outcome.occupied_edges, expected_occupied.len());

                let reference = reference_profile(vertices, &adjacency, sources);
                for (vertex, &distance) in reference.iter().enumerate() {
                    let expected = if distance == usize::MAX {
                        f64::INFINITY
                    } else {
                        distance as f64
                    };
                    assert_eq!(
                        workspace.profile()[vertex].to_bits(),
                        expected.to_bits(),
                        "profile entry {vertex} (sides {:?}, beta {}, s {}, seed {seed}, \
                         sources {sources:?})",
                        config.sides,
                        config.beta,
                        config.s
                    );
                }
                let finite = reference.iter().filter(|&&d| d != usize::MAX).count();
                assert_eq!(outcome.source_component_size, finite);
                assert_eq!(outcome.unreachable_count, vertices - finite);
                assert_eq!(outcome.diameter, Some(expected_diameter));
                profile_cases += 1;
                diameter_cases += 1;
            }
        }
    }

    // Estimator statistics vs manual computation, including the beta < 1
    // reachable conditioning (INFINITY entries excluded, fraction < 1).
    let box_ = LongRangeBox::new(&[16]).expect("box");
    let mut sampler = ChemicalDistanceSampler::new(box_.clone(), 0.3, 0.8, &[8]).expect("sampler");
    let mut estimator = ChemicalDistanceEstimator::new(box_.vertices()).expect("estimator");
    let mut manual: Vec<Vec<f64>> = Vec::new();
    for seed in 0..40_u64 {
        let profile = sampler
            .realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004D_1102 + seed))
            .expect("realization");
        estimator.add_profile(profile).expect("add");
        manual.push(profile.to_vec());
    }
    let estimate = estimator.estimate().expect("estimate");
    for vertex in 0..box_.vertices() {
        let finite: Vec<f64> = manual
            .iter()
            .map(|profile| profile[vertex])
            .filter(|value| value.is_finite())
            .collect();
        let fraction = finite.len() as f64 / manual.len() as f64;
        assert_eq!(estimate.reachable_fraction[vertex], fraction);
        if finite.is_empty() {
            assert!(estimate.profile_mean[vertex].is_nan());
            assert_eq!(
                estimate.profile_standard_error[vertex], 0.0,
                "never-reached entries carry zero error, not NaN"
            );
        } else {
            let mean = finite.iter().sum::<f64>() / finite.len() as f64;
            assert!((estimate.profile_mean[vertex] - mean).abs() < 1e-12);
            if finite.len() >= 2 {
                let variance = finite
                    .iter()
                    .map(|value| (value - mean).powi(2))
                    .sum::<f64>()
                    / (finite.len() - 1) as f64;
                let standard_error = (variance / finite.len() as f64).sqrt();
                assert!((estimate.profile_standard_error[vertex] - standard_error).abs() < 1e-12);
            }
        }
    }
    // The far corner is unreachable in some beta = 0.3 realizations —
    // the conditioning is live, not vacuous.
    assert!(
        estimate
            .reachable_fraction
            .iter()
            .any(|&fraction| fraction < 1.0),
        "the beta = 0.3 spot check must produce disconnections"
    );
    assert!(
        estimate.reachable_fraction[8] == 1.0,
        "the source is always reachable"
    );
    eprintln!(
        "L2 references: {profile_cases} profiles + {diameter_cases} diameters vs independent \
         BFS on 3 boxes x 5 seeds; estimator parity on 40 beta<1 realizations (min reachable \
         fraction {:.2})",
        estimate
            .reachable_fraction
            .iter()
            .cloned()
            .fold(1.0_f64, f64::min)
    );

    // Never-reached entries: beta = 0 leaves the source isolated. The
    // mean is NaN, the error is 0.0 (not NaN), and the reachable
    // fraction is 0 — the documented conventions, exercised directly.
    let empty_box = LongRangeBox::new(&[4]).expect("box");
    let mut empty_sampler =
        ChemicalDistanceSampler::new(empty_box, 0.0, 1.0, &[0]).expect("sampler");
    let mut empty_estimator = ChemicalDistanceEstimator::new(4).expect("estimator");
    for seed in 0..3_u64 {
        empty_sampler
            .realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004D_1103 + seed))
            .expect("isolated realization");
        empty_estimator
            .add_profile(empty_sampler.profile())
            .expect("add");
    }
    let empty = empty_estimator.estimate().expect("estimate");
    assert_eq!(empty.profile_mean[0], 0.0);
    assert!(empty.profile_mean[1].is_nan());
    assert_eq!(empty.profile_standard_error[1], 0.0);
    assert_eq!(empty.reachable_fraction[0], 1.0);
    assert_eq!(empty.reachable_fraction[1], 0.0);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: bitwise draw contract
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn bitwise_contract_and_l1_stream_parity() {
    let box_ = LongRangeBox::new(&[32]).expect("box");
    let (beta, s) = (0.7_f64, 1.3_f64);
    let seed = 0x004D_2101_u64;

    let mut first_sampler =
        ChemicalDistanceSampler::new(box_.clone(), beta, s, &[16]).expect("sampler");
    let first_outcome = first_sampler
        .realize_with_diameter(&mut Xoshiro256PlusPlus::seed_from_u64(seed))
        .expect("realization");
    let first = (
        first_sampler.profile().to_vec(),
        first_outcome,
        workspace_occupied_pairs(&box_, first_sampler.workspace()),
    );
    let mut second_sampler =
        ChemicalDistanceSampler::new(box_.clone(), beta, s, &[16]).expect("sampler");
    let second_outcome = second_sampler
        .realize_with_diameter(&mut Xoshiro256PlusPlus::seed_from_u64(seed))
        .expect("realization");
    let second = (
        second_sampler.profile().to_vec(),
        second_outcome,
        workspace_occupied_pairs(&box_, second_sampler.workspace()),
    );
    assert_eq!(first.0.len(), box_.vertices());
    assert!(
        first
            .0
            .iter()
            .zip(&second.0)
            .all(|(left, right)| left.to_bits() == right.to_bits()),
        "bit-identical profiles at the same seed"
    );
    assert_eq!(
        first.1, second.1,
        "bit-identical outcome (diameter included)"
    );
    assert_eq!(first.2, second.2, "bit-identical adjacency");

    // The draw stream is the L1 adapter's: the same seed realizes the
    // same occupied graph through the L1 core.
    let mut l1_process = LongRangePercolationMC::new(box_.clone(), beta, s, 4).expect("L1 adapter");
    l1_process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    assert_eq!(
        first.2,
        l1_process.workspace().occupied_pairs(),
        "the L2 occupancy replay equals the L1 core's occupied set"
    );

    // Different seeds realize different graphs (deterministic for these
    // fixed seeds; a collision would fail loudly).
    first_sampler
        .realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed + 1))
        .expect("other seed");
    assert_ne!(
        first_sampler.profile().to_vec(),
        first.0,
        "different seeds differ"
    );
    eprintln!(
        "L2 bitwise: profile/diameter/adjacency bit-identical at one seed, L1 stream parity ok"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: the (1,2) regime — typical distance (the sharp Biskup theorem)
// ─────────────────────────────────────────────────────────────────────────

/// One (s, r-grid) row: pooled mean and conservative SEM of the typical
/// distance at each r (both box directions at |x - y| = r, Cauchy-Schwarz
/// bound (sem+ + sem-)/2 on the average), plus the reachable fraction.
fn typical_distance_row(
    s: f64,
    n: usize,
    realizations: u64,
    grid: &[usize],
    seed: u64,
) -> (Vec<f64>, Vec<f64>, f64) {
    let beta = 1.0_f64;
    let center = n / 2;
    let box_ = LongRangeBox::new(&[n]).expect("box");
    let mut sampler =
        ChemicalDistanceSampler::new(box_.clone(), beta, s, &[center]).expect("sampler");
    let mut estimator = ChemicalDistanceEstimator::new(n).expect("estimator");
    for index in 0..realizations {
        sampler
            .realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed + index))
            .expect("realization");
        assert_eq!(
            sampler.last_outcome().unreachable_count,
            0,
            "beta = 1 keeps the 1D box connected (s = {s}, realization {index})"
        );
        estimator.add_profile(sampler.profile()).expect("add");
    }
    let estimate = estimator.estimate().expect("estimate");
    let mut means = Vec::with_capacity(grid.len());
    let mut sems = Vec::with_capacity(grid.len());
    let mut min_fraction = 1.0_f64;
    for &r in grid {
        let (upper, lower) = (center + r, center - r);
        for vertex in [upper, lower] {
            assert_eq!(
                estimate.reachable_fraction[vertex], 1.0,
                "endpoint pairs stay connected at beta = 1"
            );
            min_fraction = min_fraction.min(estimate.reachable_fraction[vertex]);
        }
        means.push(0.5 * (estimate.profile_mean[upper] + estimate.profile_mean[lower]));
        sems.push(
            0.5 * (estimate.profile_standard_error[upper] + estimate.profile_standard_error[lower]),
        );
    }
    (means, sems, min_fraction)
}

/// Per-realization pooled distances at each grid radius (both box
/// directions averaged), for statistics that must be paired across
/// realizations.
fn typical_distance_samples(
    s: f64,
    n: usize,
    realizations: u64,
    grid: &[usize],
    seed: u64,
) -> Vec<Vec<f64>> {
    let center = n / 2;
    let box_ = LongRangeBox::new(&[n]).expect("box");
    let mut sampler = ChemicalDistanceSampler::new(box_, 1.0, s, &[center]).expect("sampler");
    let mut samples = Vec::with_capacity(realizations as usize);
    for index in 0..realizations {
        sampler
            .realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed + index))
            .expect("realization");
        assert_eq!(
            sampler.last_outcome().unreachable_count,
            0,
            "beta = 1 keeps the 1D box connected (s = {s}, realization {index})"
        );
        let profile = sampler.profile();
        samples.push(
            grid.iter()
                .map(|&r| 0.5 * (profile[center + r] + profile[center - r]))
                .collect(),
        );
    }
    samples
}

/// Mean and SEM of a non-empty sample.
fn mean_sem(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    (mean, (variance / n).sqrt())
}

#[test]
fn typical_distance_polylog_regime() {
    // beta = 1: nearest neighbors surely occupied, zero disconnections
    // (asserted per realization inside the row helper). n = 2048 keeps
    // the O(pairs) replay affordable in the debug test profile.
    const N: usize = 2048;
    const REALIZATIONS: u64 = 150;
    const GRID: [usize; 7] = [8, 16, 32, 64, 128, 256, 512];
    // t_{0.975} on df = 5 (7 points minus two parameters) and df = 2
    // (the 4-point windows).
    const T975_DF5: f64 = 2.5706;
    const T975_DF2: f64 = 4.3027;

    let mut full_slopes: Vec<(f64, f64)> = Vec::new(); // (slope, CI half-width)
    for s in [1.3_f64, 1.5, 1.8] {
        let (means, sems, min_fraction) =
            typical_distance_row(s, N, REALIZATIONS, &GRID, 0x004D_3100 + (s * 10.0) as u64);
        let delta = biskup_delta(s);

        // Strictly monotone growth, every step resolved against the
        // combined standard errors.
        for index in 0..GRID.len() - 1 {
            let margin = 2.0 * (sems[index].powi(2) + sems[index + 1].powi(2)).sqrt();
            assert!(
                means[index + 1] - means[index] > margin,
                "s = {s}: D must grow from r = {} to {} ({} vs {})",
                GRID[index],
                GRID[index + 1],
                means[index],
                means[index + 1]
            );
        }

        // Fits. Raw polylog exponent: ln D vs ln ln r, full grid and the
        // top-4 window (the least-biased window available in the box).
        let full_points: Vec<(f64, f64)> = GRID
            .iter()
            .zip(&means)
            .map(|(&r, &mean)| ((r as f64).ln().ln(), mean.ln()))
            .collect();
        let (full_slope, _, full_se, full_r2) = ols(&full_points);
        let top_points: Vec<(f64, f64)> = GRID[GRID.len() - 4..]
            .iter()
            .zip(&means[means.len() - 4..])
            .map(|(&r, &mean)| ((r as f64).ln().ln(), mean.ln()))
            .collect();
        let (top_slope, _, top_se, top_r2) = ols(&top_points);

        // Finite-window bound: a power law with exponent >= 0.5 would
        // not fit these D(r) values (the polylog prediction is far
        // flatter); this is not a proof of sub-polynomial behavior.
        let power_points: Vec<(f64, f64)> = GRID
            .iter()
            .zip(&means)
            .map(|(&r, &mean)| ((r as f64).ln(), mean.ln()))
            .collect();
        let (power_slope, _, power_se, power_r2) = ols(&power_points);
        assert!(
            power_slope + T975_DF5 * power_se < 0.5,
            "s = {s}: D vs r log-log slope {:.3} + {:.3} must stay below 0.5",
            power_slope,
            T975_DF5 * power_se
        );

        // The honest finite-size statement: the fitted polylog exponent
        // sits strictly BELOW the asymptotic Delta (the o(1) corrections
        // at ln r <= 7 are large) and rises toward it as the window
        // moves up.
        assert!(
            full_slope + T975_DF5 * full_se < delta,
            "s = {s}: fitted exponent {:.3} must sit below Delta = {delta:.3}",
            full_slope
        );
        assert!(
            top_slope > full_slope,
            "s = {s}: the top-window exponent {:.3} must exceed the full-window {:.3} \
             (approach from below)",
            top_slope,
            full_slope
        );

        eprintln!(
            "L2 typical s={s} beta=1 n={N} R={REALIZATIONS} min_reachable={min_fraction:.3}: \
             D(r)={}",
            means
                .iter()
                .map(|mean| format!("{mean:.3}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        eprintln!(
            "L2 fit s={s}: Delta_paper={delta:.3} | polylog exponent full={full_slope:.3}±{:.3} \
             (R2 {full_r2:.4}) top4(r>={})={top_slope:.3}±{:.3} (R2 {top_r2:.4}) | D-vs-r slope \
             {power_slope:.3}±{:.3} (R2 {power_r2:.4})",
            T975_DF5 * full_se,
            GRID[GRID.len() - 4],
            T975_DF2 * top_se,
            T975_DF5 * power_se
        );
        full_slopes.push((full_slope, T975_DF5 * full_se));
    }

    // The paper's ORDERING: Delta(s) strictly increasing in s, and so
    // are the fitted exponents, with disjoint CIs.
    for (left, right) in full_slopes.windows(2).map(|pair| (pair[0], pair[1])) {
        assert!(
            left.0 + left.1 < right.0 - right.1,
            "fitted polylog exponents must be CI-disjoint and increasing in s: \
             {:.3}±{:.3} then {:.3}±{:.3}",
            left.0,
            left.1,
            right.0,
            right.1
        );
    }
    eprintln!(
        "L2 typical ordering: Delta(1.3)={:.3} < Delta(1.5)={:.3} < Delta(1.8)={:.3}; fitted \
         {:.3} < {:.3} < {:.3} (all below their Delta — the documented finite-size deviation)",
        biskup_delta(1.3),
        biskup_delta(1.5),
        biskup_delta(1.8),
        full_slopes[0].0,
        full_slopes[1].0,
        full_slopes[2].0
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: s > 2 — linear typical distance with positive time constant
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn linear_regime_slope_at_s_2_5() {
    // n = 8192 is a profile-mode run (the 2^12 cap applies to the
    // all-vertex diameter only). At finite r the effective log-log slope
    // sits below 1 (an O(1/r) deficit), so the gate is paired across
    // realizations instead of asking one aggregate fit to straddle 1:
    // each realization contributes its own slope and its own
    // `D/r = A + c/r` fit, and the aggregate 95% CIs must exclude the
    // polylog regime, keep the slope consistent with linear growth,
    // require a positive deficit c, and place A inside (0, 1).
    const N: usize = 8192;
    const REALIZATIONS: u64 = 24;
    const GRID: [usize; 4] = [256, 512, 1024, 2048];
    // t_{0.975} on df = 23 (24 realizations minus one).
    const T975_DF23: f64 = 2.0687;
    let s = 2.5_f64;
    let samples = typical_distance_samples(s, N, REALIZATIONS, &GRID, 0x004D_4100);

    let mut columns = vec![Vec::with_capacity(REALIZATIONS as usize); GRID.len()];
    let mut slopes = Vec::with_capacity(REALIZATIONS as usize);
    let mut deficits = Vec::with_capacity(REALIZATIONS as usize);
    let mut time_constants = Vec::with_capacity(REALIZATIONS as usize);
    for sample in &samples {
        for (index, &value) in sample.iter().enumerate() {
            columns[index].push(value);
        }
        let log_points: Vec<(f64, f64)> = GRID
            .iter()
            .zip(sample)
            .map(|(&r, &mean)| ((r as f64).ln(), mean.ln()))
            .collect();
        let (slope, _, _, _) = ols(&log_points);
        slopes.push(slope);
        let ratio_points: Vec<(f64, f64)> = GRID
            .iter()
            .zip(sample)
            .map(|(&r, &mean)| (1.0 / r as f64, mean / r as f64))
            .collect();
        let (deficit, time_constant, _, _) = ols(&ratio_points);
        deficits.push(deficit);
        time_constants.push(time_constant);
    }
    let aggregate: Vec<(f64, f64)> = columns.iter().map(|values| mean_sem(values)).collect();
    let means: Vec<f64> = aggregate.iter().map(|&(mean, _)| mean).collect();
    let (slope_mean, slope_se) = mean_sem(&slopes);
    let (deficit_mean, deficit_se) = mean_sem(&deficits);
    let (time_constant_mean, time_constant_se) = mean_sem(&time_constants);
    let slope_half = T975_DF23 * slope_se;
    let deficit_half = T975_DF23 * deficit_se;
    let time_constant_half = T975_DF23 * time_constant_se;

    eprintln!(
        "L2 linear s={s} beta=1 n={N} R={REALIZATIONS}: D(r)={} ratios={}",
        means
            .iter()
            .map(|mean| format!("{mean:.1}"))
            .collect::<Vec<_>>()
            .join(", "),
        GRID.iter()
            .zip(&means)
            .map(|(&r, &mean)| format!("{:.3}", mean / r as f64))
            .collect::<Vec<_>>()
            .join(", ")
    );
    eprintln!(
        "L2 linear paired: slope={slope_mean:.3}±{slope_half:.3} \
         CI=({:.3}, {:.3}) | A={time_constant_mean:.4}±{time_constant_half:.4} \
         CI=({:.3}, {:.3}) | c={deficit_mean:.3}±{deficit_half:.3} \
         CI=({:.3}, {:.3})",
        slope_mean - slope_half,
        slope_mean + slope_half,
        time_constant_mean - time_constant_half,
        time_constant_mean + time_constant_half,
        deficit_mean - deficit_half,
        deficit_mean + deficit_half
    );

    assert!(
        slope_mean - slope_half > 0.8,
        "s = {s}: paired slope CI ({:.3}, {:.3}) lower bound must exclude the polylog regime \
         (measured slopes <= 0.35 there)",
        slope_mean - slope_half,
        slope_mean + slope_half
    );
    assert!(
        slope_mean + slope_half < 1.05,
        "s = {s}: paired slope CI ({:.3}, {:.3}) must stay consistent with linear growth \
         (the finite-window slope sits at or below 1)",
        slope_mean - slope_half,
        slope_mean + slope_half
    );
    assert!(
        deficit_mean - deficit_half > 0.0,
        "s = {s}: the paired 1/r deficit c={deficit_mean:.3}±{deficit_half:.3} (95% CI) \
         must be positive (D/r approaches A from above)"
    );
    assert!(
        time_constant_mean - time_constant_half > 0.15
            && time_constant_mean + time_constant_half < 1.0,
        "s = {s}: paired extrapolated time constant A={time_constant_mean:.3}±\
         {time_constant_half:.3} (95% CI) must be positive and inside (0, 1)"
    );
    // The measured ratio does not increase toward A beyond the combined
    // noise, and strictly decreases over the window.
    for index in 0..GRID.len() - 1 {
        let left = aggregate[index].0 / GRID[index] as f64;
        let right = aggregate[index + 1].0 / GRID[index + 1] as f64;
        let left_sem = aggregate[index].1 / GRID[index] as f64;
        let right_sem = aggregate[index + 1].1 / GRID[index + 1] as f64;
        let margin = 2.0 * (left_sem * left_sem + right_sem * right_sem).sqrt();
        assert!(
            right <= left + margin,
            "s = {s}: D(r)/r must not increase from r = {} to {} ({left:.4} -> {right:.4}, \
             margin {margin:.4})",
            GRID[index],
            GRID[index + 1]
        );
    }
    assert!(
        means[0] / GRID[0] as f64 > means[GRID.len() - 1] / GRID[GRID.len() - 1] as f64,
        "s = {s}: the time-constant ratio must decrease over the window"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5: s <= 1 diameter envelopes (never the typical distance)
// ─────────────────────────────────────────────────────────────────────────

/// Mean and SEM of the exact box diameter over `realizations` draws.
fn diameter_row(s: f64, n: usize, realizations: u64, seed: u64) -> (f64, f64) {
    let beta = 1.0_f64;
    let box_ = LongRangeBox::new(&[n]).expect("box");
    let mut sampler = ChemicalDistanceSampler::new(box_.clone(), beta, s, &[0]).expect("sampler");
    let mut estimator = ChemicalDistanceEstimator::new(n).expect("estimator");
    for index in 0..realizations {
        let outcome = sampler
            .realize_with_diameter(&mut Xoshiro256PlusPlus::seed_from_u64(seed + index))
            .expect("realization");
        assert_eq!(
            outcome.unreachable_count, 0,
            "beta = 1 keeps the box connected"
        );
        estimator.add_diameter(outcome.diameter.expect("tracking was requested"));
    }
    let estimate = estimator.estimate().expect("estimate");
    (estimate.diameter_mean, estimate.diameter_standard_error)
}

#[test]
fn diameter_envelope_below_s_one_is_bounded() {
    // Eq. (1.8): sup_{x,y in Z^d} D(x,y) = d/(d-s) a.s. — a CONSTANT
    // (2 at d = 1, s = 0.5). The box graph is a subgraph of Z^d with no
    // outside shortcut vertices, so finite-box diameters may sit above
    // that sup; Biskup's "similar asymptotic statement" for D_N is
    // quoted qualitatively, not as an equality. This is a BOUNDEDNESS
    // envelope, not a measurement of the constant 2: the lower edge is
    // nearly automatic for a connected non-complete box, and the upper
    // edge plus the absence of power growth is the actual gate.
    const SIZES: [usize; 4] = [256, 512, 1024, 2048];
    const REALIZATIONS: u64 = 8;
    let mut means = Vec::with_capacity(SIZES.len());
    for &n in &SIZES {
        let (mean, sem) = diameter_row(0.5, n, REALIZATIONS, 0x004D_5100 + n as u64);
        eprintln!(
            "L2 diameter s=0.5 beta=1 n={n} R={REALIZATIONS}: mean={mean:.3} sem={sem:.3} \
             (Eq 1.8 constant d/(d-s) = 2.00 reference scale)"
        );
        assert!(
            (2.0..=4.0).contains(&mean),
            "s = 0.5: diameter {mean:.3} must stay in the bounded envelope [2, 4] \
             (the finite-box value is not required to equal the a.s. sup 2)"
        );
        means.push(mean);
    }
    // At this scale every realization gives the same diameter (the
    // review confirmed 130 fresh diameters, all 3). A zero-residual OLS
    // has no defined R^2, so state the constant case explicitly rather
    // than reporting a degenerate fit.
    if means.windows(2).all(|pair| pair[0] == pair[1]) {
        eprintln!(
            "L2 diameter s=0.5: constant sequence (diameter {:.3} at every N; no power growth, \
             R2 undefined)",
            means[0]
        );
    } else {
        let points: Vec<(f64, f64)> = SIZES
            .iter()
            .zip(&means)
            .map(|(&n, &mean)| ((n as f64).ln(), mean.ln()))
            .collect();
        let (slope, _, slope_se, r_squared) = ols(&points);
        // t_{0.975} on df = 2 (four ladder points minus two parameters).
        const T975_DF2: f64 = 4.3027;
        let (low, high) = (slope - T975_DF2 * slope_se, slope + T975_DF2 * slope_se);
        eprintln!(
            "L2 diameter s=0.5 fit: slope={slope:.4}±{:.4} CI=({low:.4}, {high:.4}) \
             R2={r_squared:.4}",
            T975_DF2 * slope_se
        );
        assert!(
            low <= 0.0 && (0.0..0.5).contains(&high),
            "s = 0.5: diameter-vs-N slope CI ({low:.3}, {high:.3}) must contain 0 and stay \
             far below any positive power"
        );
    }
}

#[test]
fn diameter_envelope_at_s_one_is_at_most_logarithmic() {
    // Eq. (1.9): D_N = Theta(log N / log log N) for the DIAMETER; the
    // typical distance is at most logarithmic and the lower bound is
    // explicitly missing. Bracket the UPPER envelope only: growth
    // consistent with zero power slope and bounded by ln N; no lower
    // claim anywhere.
    const SIZES: [usize; 4] = [512, 1024, 2048, 4096];
    const REALIZATIONS: u64 = 3;
    let mut points = Vec::with_capacity(SIZES.len());
    for &n in &SIZES {
        let (mean, sem) = diameter_row(1.0, n, REALIZATIONS, 0x004D_5200 + n as u64);
        let envelope = (n as f64).ln();
        let scaled = mean / envelope.ln();
        eprintln!(
            "L2 diameter s=1.0 beta=1 n={n} R={REALIZATIONS}: mean={mean:.3} sem={sem:.3} \
             vs ln N={envelope:.2} (ratio {scaled:.2}; Eq 1.9 form ln N/ln ln N = {:.2})",
            envelope / envelope.ln()
        );
        assert!(
            mean <= envelope,
            "s = 1: diameter {mean:.3} exceeds the ln N = {envelope:.2} upper envelope"
        );
        points.push(((n as f64).ln(), mean.ln()));
    }
    let (slope, _, slope_se, r_squared) = ols(&points);
    const T975_DF2: f64 = 4.3027;
    let (low, high) = (slope - T975_DF2 * slope_se, slope + T975_DF2 * slope_se);
    eprintln!(
        "L2 diameter s=1.0 fit: slope={slope:.4}±{:.4} CI=({low:.4}, {high:.4}) R2={r_squared:.4}",
        T975_DF2 * slope_se
    );
    assert!(
        low <= 0.0 && (0.0..0.5).contains(&high),
        "s = 1: diameter-vs-N slope CI ({low:.3}, {high:.3}) must contain 0 and stay far \
         below 1 (at most logarithmic, upper envelope only)"
    );
}

// Gate 6 (s = 2 not bracketed) is a documentation claim, not a run: the
// module documentation of `percolation::long_range_distances` and the
// PERCOLATION.md section record why (Biskup ~line 255: the critical case
// s = 2d is "not very well understood"). No test may run s = 2 gates.

// ─────────────────────────────────────────────────────────────────────────
// Typed rejections
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn typed_rejections_leave_state_untouched() {
    let box_ = LongRangeBox::new(&[6]).expect("box");
    let mut workspace = ChemicalDistanceWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004D_6101);
    let uniforms = draw_uniforms(box_.pairs(), &mut rng);
    chemical_distance_run(
        &box_,
        0.8,
        0.5,
        &uniforms,
        &[2],
        DiameterTracking::AllVertexBfs,
        &mut workspace,
    )
    .expect("valid run");
    let before = (
        workspace.profile().to_vec(),
        workspace.occupied_edge_count(),
    );

    for beta in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            chemical_distance_run(&box_, beta, 0.5, &uniforms, &[2], DiameterTracking::Off, &mut workspace),
            Err(ChemicalDistanceError::BetaInvalid { value }) if value.to_bits() == beta.to_bits()
        ));
    }
    for s in [-0.5, f64::NAN] {
        assert!(matches!(
            chemical_distance_run(
                &box_,
                0.8,
                s,
                &uniforms,
                &[2],
                DiameterTracking::Off,
                &mut workspace
            ),
            Err(ChemicalDistanceError::ExponentInvalid { .. })
        ));
    }
    assert!(matches!(
        chemical_distance_run(
            &box_,
            0.8,
            0.5,
            &uniforms,
            &[],
            DiameterTracking::Off,
            &mut workspace
        ),
        Err(ChemicalDistanceError::NoSources)
    ));
    assert!(matches!(
        chemical_distance_run(
            &box_,
            0.8,
            0.5,
            &uniforms,
            &[2, 6],
            DiameterTracking::Off,
            &mut workspace
        ),
        Err(ChemicalDistanceError::SourceOutOfRange {
            index: 1,
            source: 6,
            vertices: 6
        })
    ));
    assert!(matches!(
        chemical_distance_run(
            &box_,
            0.8,
            0.5,
            &uniforms[..4],
            &[2],
            DiameterTracking::Off,
            &mut workspace
        ),
        Err(ChemicalDistanceError::UniformLengthMismatch {
            expected: 15,
            actual: 4
        })
    ));
    let mut bad = uniforms.clone();
    bad[3] = -0.25;
    assert!(matches!(
        chemical_distance_run(
            &box_,
            0.8,
            0.5,
            &bad,
            &[2],
            DiameterTracking::Off,
            &mut workspace
        ),
        Err(ChemicalDistanceError::UniformInvalid { index: 3, .. })
    ));

    // The all-vertex diameter cost cap: boxes above 2^12 vertices are a
    // typed rejection BEFORE any buffer work; profiles are not capped.
    let huge = LongRangeBox::new(&[DIAMETER_VERTEX_CAP + 1]).expect("box");
    assert!(matches!(
        chemical_distance_run(&huge, 1.0, 1.0, &[], &[0], DiameterTracking::AllVertexBfs, &mut workspace),
        Err(ChemicalDistanceError::DiameterCapExceeded { vertices, cap }) if vertices == DIAMETER_VERTEX_CAP + 1 && cap == DIAMETER_VERTEX_CAP
    ));

    // Sampler and estimator domains.
    assert!(matches!(
        ChemicalDistanceSampler::new(box_.clone(), 0.8, 0.5, &[]),
        Err(ChemicalDistanceError::NoSources)
    ));
    assert!(matches!(
        ChemicalDistanceSampler::new(box_.clone(), 0.8, 0.5, &[6]),
        Err(ChemicalDistanceError::SourceOutOfRange { source: 6, .. })
    ));
    assert!(matches!(
        ChemicalDistanceEstimator::new(0),
        Err(ChemicalDistanceError::NoVertices)
    ));
    let mut estimator = ChemicalDistanceEstimator::new(4).expect("estimator");
    assert!(matches!(
        estimator.add_profile(&[0.0; 3]),
        Err(ChemicalDistanceError::ProfileLength {
            expected: 4,
            actual: 3
        })
    ));
    assert_eq!(estimator.estimate(), None, "nothing accumulated yet");
    assert_eq!(estimator.profile_samples(), 0);

    // Public dimension/arithmetic overflow paths are typed rejections,
    // never panics, and are detected before any buffer mutation.
    assert!(matches!(
        workspace.prepare(usize::MAX, 0, 0),
        Err(ChemicalDistanceError::DimensionOverflow { field: "offsets" })
    ));
    assert!(matches!(
        workspace.prepare(1, 0, usize::MAX),
        Err(ChemicalDistanceError::DimensionOverflow {
            field: "probabilities"
        })
    ));
    assert!(matches!(
        workspace.prepare(1, usize::MAX, 0),
        Err(ChemicalDistanceError::DimensionOverflow { field: "neighbors" })
    ));
    let mut huge_estimator = ChemicalDistanceEstimator::new(4).expect("estimator");
    huge_estimator.add_diameter(usize::MAX);
    huge_estimator.add_diameter(usize::MAX);
    let huge_estimate = huge_estimator.estimate().expect("estimate");
    assert!(
        huge_estimate.diameter_mean.is_finite(),
        "the diameter square is taken in f64, so a usize-sized input cannot overflow"
    );

    // Every rejection left the last valid realization untouched.
    assert_eq!(workspace.profile(), &before.0[..]);
    assert_eq!(workspace.occupied_edge_count(), before.1);
    eprintln!("L2 typed rejections: 18 refusal cases, state untouched");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 8: throughput report (eprintln-timed, debug profile)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn throughput_report_two_sizes() {
    for length in [256_usize, 1024] {
        let box_ = LongRangeBox::new(&[length]).expect("1D box");
        let mut sampler = ChemicalDistanceSampler::new(box_.clone(), 1.0, 1.5, &[length / 2])
            .expect("throughput sampler");
        let mut estimator = ChemicalDistanceEstimator::new(length).expect("estimator");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004D_7101);
        for _ in 0..4 {
            estimator
                .add_profile(sampler.realize(&mut rng).expect("warmup"))
                .expect("warmup add");
        }
        let realizations = 48_usize;
        let start = Instant::now();
        for _ in 0..realizations {
            estimator
                .add_profile(sampler.realize(&mut rng).expect("realization"))
                .expect("add");
        }
        let elapsed = start.elapsed();
        let per_run = elapsed.as_secs_f64() / realizations as f64;
        eprintln!(
            "L2_THROUGHPUT box=1x{length} V={length} pairs={} profile_realizations={} \
             total={elapsed:?} per_realization={per_run_us:.0}us realizations_per_second={rps:.1} \
             pairs_per_second={pps:.3e}",
            box_.pairs(),
            realizations,
            per_run_us = per_run * 1.0e6,
            rps = 1.0 / per_run,
            pps = box_.pairs() as f64 / per_run,
        );

        // The all-vertex diameter cost at a small box (the O(V(V+m))
        // mode the cap is about).
        let small = LongRangeBox::new(&[256]).expect("small box");
        let mut diameter_sampler =
            ChemicalDistanceSampler::new(small.clone(), 1.0, 1.5, &[0]).expect("sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004D_7102);
        for _ in 0..2 {
            diameter_sampler
                .realize_with_diameter(&mut rng)
                .expect("warmup");
        }
        let start = Instant::now();
        let repetitions = 8_usize;
        for _ in 0..repetitions {
            diameter_sampler
                .realize_with_diameter(&mut rng)
                .expect("diameter realization");
        }
        let elapsed = start.elapsed();
        let per_run = elapsed.as_secs_f64() / repetitions as f64;
        eprintln!(
            "L2_THROUGHPUT box=1x256 diameter_mode repetitions={repetitions} total={elapsed:?} \
             per_realization={per_run_ms:.1}ms realizations_per_second={:.1}",
            1.0 / per_run,
            per_run_ms = per_run * 1.0e3,
        );
    }
}
