//! F2 first-passage observables validation (roadmap v2 §3 F2).
//!
//! Gates, verbatim from the stage brief:
//!
//! 1. **Exact enumeration** — corridors and segments at two-level dyadic
//!    weights {a, b} = {0.5, 1.0} (every path sum exact in f64): enumerate
//!    ALL 2^E box-edge assignments, compute the passage-time profile per
//!    assignment via a TEST-LOCAL Bellman-Ford over a test-built edge list
//!    (sharing no code with production), and verify the production profile
//!    BIT-EXACTLY per assignment. Includes tie census: assignments whose
//!    profile minimum is attained at several boundary rows are counted and
//!    asserted to exist (the profile is a vector of unique minima over
//!    path weights — path degeneracy never changes its content; the
//!    scalar minimum is a value, no row is singled out).
//! 2. **Bitwise draw-contract replay** — the documented contract (one
//!    `random::<f64>()` per box edge in dense edge order, then structural
//!    zeros) replayed from a seed reproduces the weight slice bitwise, for
//!    both laws and both geometries.
//! 3. **Statistical brackets** (bracket language ONLY, never exact
//!    claims): (a) Var(T_r) growth on aspect-1 corridors (w = r, so the
//!    transversal wandering r^{2/3} never meets the walls) for Exponential
//!    and Uniform laws — the fitted log-Var vs log-r exponent chi with CI
//!    must sit inside the rigorous open interval (0, 1/2) (Kesten's
//!    0 <= chi <= 1/2; the 1/3 prediction is never claimed); (b)
//!    subadditivity E[T_{2r}] <= 2 E[T_r] on reflection-symmetric segment
//!    boxes where it is EXACT in expectation (the reflection about the
//!    box's center line swaps origin and 2r-target and preserves the
//!    i.i.d. law); (c) the time-constant sequence E[T_r]/r on one segment
//!    box, decreasing along r multiples and stabilizing.
//! 4. **Zero-variance determinism** — same seed reproduces profiles
//!    bitwise; SEM shrinks like 1/sqrt(N).
//! 5. Throughput probe (allocations live in the
//!    `percolation_allocations` bench).

use cmc_rs::{
    EdgeWeightLaw, FirstPassageError, PassageEstimate, PassageEstimator, PassageGeometry,
    PassageObservablesError, PassageSampler,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Two-level enumeration weights: dyadic, so every path sum is exact in
/// f64 and production/reference compare BITWISE however the two
/// algorithms order their additions.
const ENUM_A: f64 = 0.5;
const ENUM_B: f64 = 1.0;

// ─────────────────────────────────────────────────────────────────────────
// Test-local references
// ─────────────────────────────────────────────────────────────────────────

/// Independent Bellman-Ford over a test-built edge list: repeated full
/// relaxation rounds to the fixpoint, both directions of every edge.
/// Shares NO code with the production Dijkstra.
fn bellman_ford(vertices: usize, edges: &[(usize, usize, f64)], source: usize) -> Vec<f64> {
    let mut dist = vec![f64::INFINITY; vertices];
    dist[source] = 0.0;
    for _ in 0..vertices {
        let mut changed = false;
        for &(u, v, w) in edges {
            if dist[u] + w < dist[v] {
                dist[v] = dist[u] + w;
                changed = true;
            }
            if dist[v] + w < dist[u] {
                dist[u] = dist[v] + w;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    dist
}

/// Test-built corridor topology (independent of the production builder):
/// all right and up bonds of the (r+1) x w box plus the super-source bonds
/// to the left column, in production's dense edge order (site-major X-then-Y
/// over the box, super-source bonds last) so an assignment mask indexes the
/// same edge in both. Returns (vertices, edges without weights).
fn corridor_topology(length: usize, width: usize) -> (usize, Vec<(usize, usize)>) {
    let stride = length + 1;
    let vertices = stride * width + 1;
    let super_source = stride * width;
    let mut edges = Vec::new();
    for y in 0..width {
        for x in 0..=length {
            if x < length {
                edges.push((x + y * stride, x + 1 + y * stride));
            }
            if y + 1 < width {
                edges.push((x + y * stride, x + (y + 1) * stride));
            }
        }
    }
    for y in 0..width {
        edges.push((super_source, y * stride));
    }
    (vertices, edges)
}

/// Test-built segment topology: the `box_width x rows` box plus ONE
/// super-source bond to the origin; targets s steps right of the origin.
fn segment_topology(
    steps: &[usize],
    margin: usize,
    half_height: usize,
) -> (usize, Vec<(usize, usize)>, usize, Vec<usize>) {
    let length = steps[steps.len() - 1];
    let box_width = length + 2 * margin + 1;
    let rows = 2 * half_height + 1;
    let vertices = box_width * rows + 1;
    let origin = margin + half_height * box_width;
    let super_source = box_width * rows;
    let mut edges = Vec::new();
    for y in 0..rows {
        for x in 0..box_width {
            if x + 1 < box_width {
                edges.push((x + y * box_width, x + 1 + y * box_width));
            }
            if y + 1 < rows {
                edges.push((x + y * box_width, x + (y + 1) * box_width));
            }
        }
    }
    edges.push((super_source, origin));
    let targets: Vec<usize> = steps.iter().map(|&step| origin + step).collect();
    (vertices, edges, super_source, targets)
}

/// Sorted endpoint-pair multiset of the production lattice (topology
/// parity pin: production must BE the test-built box plus super-source).
fn production_pairs(geometry: &PassageGeometry) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = geometry
        .lattice()
        .edges
        .iter()
        .map(|bond| {
            let (source, target) = (bond.source, bond.target);
            (source.min(target), source.max(target))
        })
        .collect();
    pairs.sort_unstable();
    pairs
}

fn test_pairs(edges: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = edges.iter().map(|&(u, v)| (u.min(v), u.max(v))).collect();
    pairs.sort_unstable();
    pairs
}

/// Ordinary least squares of y on x: (slope, intercept, slope standard
/// error, R^2). Test-local; the gates feed it (ln r, ln Var) points.
fn ols(points: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let n = points.len() as f64;
    let x_mean = points.iter().map(|&(x, _)| x).sum::<f64>() / n;
    let y_mean = points.iter().map(|&(_, y)| y).sum::<f64>() / n;
    let s_xx = points
        .iter()
        .map(|&(x, _)| (x - x_mean) * (x - x_mean))
        .sum::<f64>();
    let s_xy = points
        .iter()
        .map(|&(x, y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>();
    let slope = s_xy / s_xx;
    let intercept = y_mean - slope * x_mean;
    let rss = points
        .iter()
        .map(|&(x, y)| {
            let residual = y - (intercept + slope * x);
            residual * residual
        })
        .sum::<f64>();
    let tss = points
        .iter()
        .map(|&(_, y)| (y - y_mean) * (y - y_mean))
        .sum::<f64>();
    let slope_se = (rss / (n - 2.0) / s_xx).sqrt();
    (slope, intercept, slope_se, 1.0 - rss / tss)
}

/// Standard error of the per-target mean from the per-target variance.
fn target_sem(estimate: &PassageEstimate, target: usize) -> f64 {
    (estimate.profile_variance[target] / estimate.samples as f64).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: exact enumeration over all 2^E box-edge assignments
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn corridor_enumeration_matches_brute_force_per_assignment() {
    let mut total_assignments = 0_u64;
    let mut total_comparisons = 0_u64;
    let mut total_tied = 0_u64;
    for &(length, width) in &[(2_usize, 4_usize), (3, 3), (3, 2)] {
        let geometry = PassageGeometry::corridor(length, width);
        let (vertices, topology) = corridor_topology(length, width);
        assert_eq!(
            production_pairs(&geometry),
            test_pairs(&topology),
            "production topology must be the test-built box plus super-source"
        );
        let drawn = geometry.drawn_edge_count();
        let super_bonds = topology.len() - drawn;
        let mut sampler = PassageSampler::new(geometry.clone()).expect("enumeration sampler");
        let mut full_weights = vec![0.0; drawn + super_bonds];
        let mut assignments = 0_u64;
        let mut tied = 0_u64;
        for mask in 0_u64..1_u64 << drawn {
            for (edge, weight) in full_weights.iter_mut().take(drawn).enumerate() {
                *weight = if mask & (1_u64 << edge) != 0 {
                    ENUM_B
                } else {
                    ENUM_A
                };
            }
            // Structural super-source zeros (already 0.0; written explicitly
            // to mirror the production contract).
            for weight in &mut full_weights[drawn..] {
                *weight = 0.0;
            }
            let profile = sampler
                .realize_with_weights(&full_weights)
                .expect("enumeration run");
            let weighted: Vec<(usize, usize, f64)> = topology
                .iter()
                .zip(full_weights.iter())
                .map(|(&(u, v), &w)| (u, v, w))
                .collect();
            let reference = bellman_ford(vertices, &weighted, geometry.super_source());
            let mut minimum = f64::INFINITY;
            for (&target, &produced) in geometry.targets().iter().zip(profile) {
                let expected = reference[target];
                assert_eq!(
                    produced, expected,
                    "corridor {length}x{width} mask {mask:#x} target {target}: bit-exact profile"
                );
                minimum = minimum.min(expected);
            }
            // Tie census: how many boundary rows attain the profile minimum.
            let attaining = profile.iter().filter(|&&t| t == minimum).count();
            if attaining >= 2 {
                tied += 1;
            }
            assert_eq!(sampler.passage_time(), minimum, "scalar min over profile");
            assignments += 1;
        }
        eprintln!(
            "F2 enumeration corridor {length}x{width}: E_box={drawn} assignments={assignments} \
             profile entries compared={assignments}x{width} tied-min assignments={tied}"
        );
        assert!(
            tied > 0,
            "corridor {length}x{width}: tie assignments must occur at two-level weights"
        );
        total_assignments += assignments;
        total_comparisons += assignments * width as u64;
        total_tied += tied;
    }
    eprintln!(
        "F2 enumeration total: {total_assignments} assignments, {total_comparisons} profile-entry \
         comparisons, {total_tied} tied-minimum assignments, 0 mismatches"
    );
}

#[test]
fn segment_enumeration_matches_brute_force_per_assignment() {
    let cases: [(&[usize], usize, usize); 3] = [(&[2], 0, 1), (&[2, 3], 1, 0), (&[1, 2], 1, 1)];
    let mut total_assignments = 0_u64;
    for &(steps, margin, half_height) in &cases {
        let geometry = PassageGeometry::segment(steps, margin, half_height);
        let (vertices, topology, super_source, targets) =
            segment_topology(steps, margin, half_height);
        assert_eq!(super_source, geometry.super_source());
        assert_eq!(targets, geometry.targets());
        assert_eq!(
            production_pairs(&geometry),
            test_pairs(&topology),
            "production topology must be the test-built box plus super-source"
        );
        let drawn = geometry.drawn_edge_count();
        let mut sampler = PassageSampler::new(geometry.clone()).expect("enumeration sampler");
        let mut full_weights = vec![0.0; topology.len()];
        let mut assignments = 0_u64;
        for mask in 0_u64..1_u64 << drawn {
            for (edge, weight) in full_weights.iter_mut().take(drawn).enumerate() {
                *weight = if mask & (1_u64 << edge) != 0 {
                    ENUM_B
                } else {
                    ENUM_A
                };
            }
            for weight in &mut full_weights[drawn..] {
                *weight = 0.0;
            }
            let profile = sampler
                .realize_with_weights(&full_weights)
                .expect("segment enumeration run");
            let weighted: Vec<(usize, usize, f64)> = topology
                .iter()
                .zip(full_weights.iter())
                .map(|(&(u, v), &w)| (u, v, w))
                .collect();
            let reference = bellman_ford(vertices, &weighted, super_source);
            for (&target, &produced) in targets.iter().zip(profile) {
                assert_eq!(
                    produced, reference[target],
                    "segment {steps:?} margin {margin} hh {half_height} mask {mask:#x}"
                );
            }
            assignments += 1;
        }
        total_assignments += assignments;
        eprintln!(
            "F2 enumeration segment steps={steps:?} margin={margin} half_height={half_height}: \
             E_box={drawn} assignments={assignments}, 0 mismatches"
        );
    }
    eprintln!("F2 segment enumeration total: {total_assignments} assignments, 0 mismatches");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: bitwise draw-contract replay
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn weight_draw_contract_replays_bitwise() {
    let laws = [
        EdgeWeightLaw::exponential(1.0),
        EdgeWeightLaw::exponential(2.5),
        EdgeWeightLaw::Uniform,
    ];
    let mut comparisons = 0_usize;
    for (index, &law) in laws.iter().enumerate() {
        // Corridor geometry.
        let geometry = PassageGeometry::corridor(6, 5);
        let drawn = geometry.drawn_edge_count();
        let mut sampler = PassageSampler::new(geometry).expect("replay sampler");
        let seed = 0x00F2_0000 + index as u64;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        sampler.realize(law, &mut rng).expect("replay realization");
        let produced = sampler.weights().to_vec();
        let mut replay = Xoshiro256PlusPlus::seed_from_u64(seed);
        for (edge, &weight) in produced.iter().enumerate() {
            let expected = if edge < drawn {
                law.draw(replay.random::<f64>())
            } else {
                0.0
            };
            assert_eq!(
                weight.to_bits(),
                expected.to_bits(),
                "law {law:?} edge {edge}: bitwise draw-contract replay"
            );
            comparisons += 1;
        }
        // Segment geometry (one structural zero at the tail).
        let geometry = PassageGeometry::segment(&[3, 6], 2, 2);
        let drawn = geometry.drawn_edge_count();
        let mut sampler = PassageSampler::new(geometry).expect("replay sampler");
        let seed = 0x00F2_8000 + index as u64;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        sampler.realize(law, &mut rng).expect("replay realization");
        let produced = sampler.weights().to_vec();
        let mut replay = Xoshiro256PlusPlus::seed_from_u64(seed);
        for (edge, &weight) in produced.iter().enumerate() {
            let expected = if edge < drawn {
                law.draw(replay.random::<f64>())
            } else {
                0.0
            };
            assert_eq!(
                weight.to_bits(),
                expected.to_bits(),
                "segment law {law:?} edge {edge}"
            );
            comparisons += 1;
        }
    }
    eprintln!("F2 draw-contract replay: {comparisons} weight entries reproduced bitwise, 0 diffs");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: determinism and SEM scaling
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn same_seed_reproduces_profiles_bitwise() {
    for law in [EdgeWeightLaw::exponential(1.0), EdgeWeightLaw::Uniform] {
        let mut first = PassageSampler::new(PassageGeometry::corridor(16, 8)).expect("sampler a");
        let mut second = PassageSampler::new(PassageGeometry::corridor(16, 8)).expect("sampler b");
        let mut rng_a = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9001);
        let mut rng_b = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9001);
        let mut chain_one = Vec::new();
        for realization in 0..50 {
            let a = first.realize(law, &mut rng_a).expect("realization a");
            let b = second.realize(law, &mut rng_b).expect("realization b");
            if realization == 0 {
                chain_one = b.to_vec();
            }
            for (index, (x, y)) in a.iter().zip(b).enumerate() {
                assert_eq!(
                    x.to_bits(),
                    y.to_bits(),
                    "law {law:?} realization {realization} entry {index}"
                );
            }
            assert_eq!(
                first.passage_time().to_bits(),
                second.passage_time().to_bits()
            );
        }

        // A re-seeded stream restarts at chain one — and the history-laden
        // sampler must reproduce that first chain bitwise (no state leak).
        let mut fresh = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9001);
        let replay = first.realize(law, &mut fresh).expect("reseed realization");
        for (index, (x, y)) in replay.iter().zip(chain_one.iter()).enumerate() {
            assert_eq!(x.to_bits(), y.to_bits(), "reseeded entry {index}");
        }

        // Segment geometry determinism (profile over two targets).
        let mut seg_a =
            PassageSampler::new(PassageGeometry::segment(&[8, 16], 4, 4)).expect("seg a");
        let mut seg_b =
            PassageSampler::new(PassageGeometry::segment(&[8, 16], 4, 4)).expect("seg b");
        let mut rng_d = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9002);
        let mut rng_e = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9002);
        for _ in 0..20 {
            let a = seg_a.realize(law, &mut rng_d).expect("seg run a");
            let b = seg_b.realize(law, &mut rng_e).expect("seg run b");
            assert_eq!(a, b, "segment profile determinism");
        }
    }
    eprintln!(
        "F2 determinism: same-seed profiles bitwise identical (2 laws x 50 corridor runs on two \
         samplers, reseed-replays-chain-one included, 2 x 20 segment runs)"
    );
}

#[test]
fn standard_error_scales_as_one_over_sqrt_n() {
    let law = EdgeWeightLaw::exponential(1.0);
    let geometry = PassageGeometry::corridor(32, 16);
    let mut short = PassageSampler::new(geometry.clone()).expect("short sampler");
    let mut long = PassageSampler::new(geometry.clone()).expect("long sampler");
    let mut short_estimator = PassageEstimator::new(16).expect("short estimator");
    let mut long_estimator = PassageEstimator::new(16).expect("long estimator");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9100);
    for _ in 0..128 {
        short_estimator
            .add(short.realize(law, &mut rng).expect("short run"))
            .expect("short add");
    }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_9200);
    for _ in 0..512 {
        long_estimator
            .add(long.realize(law, &mut rng).expect("long run"))
            .expect("long add");
    }
    let short = short_estimator.estimate().expect("short estimate");
    let long = long_estimator.estimate().expect("long estimate");
    let ratio = long.standard_error / short.standard_error;
    eprintln!(
        "F2 SEM scaling: N=128 SEM={:.5}, N=512 SEM={:.5}, ratio={:.3} (1/sqrt(4) = 0.5)",
        short.standard_error, long.standard_error, ratio
    );
    assert!(
        long.standard_error < short.standard_error,
        "SEM must shrink"
    );
    assert!(
        ratio > 0.35 && ratio < 0.65,
        "SEM ratio {ratio} must bracket 1/sqrt(4) = 0.5 within noise"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3a: Var(T_r) growth on corridors (Exponential and Uniform)
// ─────────────────────────────────────────────────────────────────────────

/// One (r, Var, mean, SEM) row: aspect-1 corridor (w = r), R realizations.
fn sample_corridor(law: EdgeWeightLaw, r: usize, realizations: u64, seed: u64) -> PassageEstimate {
    let geometry = PassageGeometry::corridor(r, r);
    let mut sampler = PassageSampler::new(geometry).expect("corridor sampler");
    let mut estimator = PassageEstimator::new(r).expect("corridor estimator");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for _ in 0..realizations {
        estimator
            .add(sampler.realize(law, &mut rng).expect("corridor run"))
            .expect("corridor add");
    }
    estimator.estimate().expect("corridor estimate")
}

#[test]
fn variance_growth_exponent_bracket_on_corridors() {
    const LADDER: [usize; 5] = [8, 16, 32, 64, 128];
    const REALIZATIONS: u64 = 600;
    // t_{0.975} on df = 3 (five ladder points minus two parameters).
    const T975_DF3: f64 = 3.1824;
    for (law_index, &law) in [EdgeWeightLaw::exponential(1.0), EdgeWeightLaw::Uniform]
        .iter()
        .enumerate()
    {
        let mut points = Vec::with_capacity(LADDER.len());
        for &r in &LADDER {
            let estimate = sample_corridor(
                law,
                r,
                REALIZATIONS,
                0x00F2_A000 + law_index as u64 * 0x100 + r as u64,
            );
            eprintln!(
                "F2 chi ladder law={law_index} r={r:3} w={r:3} N={REALIZATIONS}: mean={:.4} \
                 var={:.5} sem={:.5}",
                estimate.mean, estimate.variance, estimate.standard_error
            );
            points.push(((r as f64).ln(), estimate.variance.ln()));
        }
        let (slope, _intercept, slope_se, r_squared) = ols(&points);
        let chi = slope / 2.0;
        let half_width = T975_DF3 * slope_se / 2.0;
        eprintln!(
            "F2 chi fit law={law_index}: fitted chi={chi:.4} CI=({:.4}, {:.4}) over r in \
             [8, 128], R={REALIZATIONS} per r, R^2={r_squared:.4} (slope 2*chi={slope:.4} ± \
             {:.4})",
            chi - half_width,
            chi + half_width,
            T975_DF3 * slope_se
        );
        assert!(
            chi - half_width > 0.0 && chi + half_width < 0.5,
            "law {law_index}: fitted chi CI ({:.4}, {:.4}) must sit inside the rigorous open \
             interval (0, 1/2)",
            chi - half_width,
            chi + half_width
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3b: subadditivity E[T_{2r}] <= 2 E[T_r] on symmetric segment boxes
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn subadditivity_bracket_on_segments() {
    let law = EdgeWeightLaw::exponential(1.0);
    for &r in &[8_usize, 16, 32, 64] {
        // Box symmetric under the reflection about its center line: margin
        // = r longitudinally, half_height = r transversally (wandering
        // r^{2/3} never meets the walls at these sizes).
        let geometry = PassageGeometry::segment(&[r, 2 * r], r, r);
        let realizations = if r <= 32 { 1200 } else { 500 };
        let mut sampler = PassageSampler::new(geometry).expect("segment sampler");
        let mut estimator = PassageEstimator::new(2).expect("segment estimator");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_B000 + r as u64);
        for _ in 0..realizations {
            estimator
                .add(sampler.realize(law, &mut rng).expect("segment run"))
                .expect("segment add");
        }
        let estimate = estimator.estimate().expect("segment estimate");
        let mu_r = estimate.profile_mean[0] / r as f64;
        let mu_2r = estimate.profile_mean[1] / (2 * r) as f64;
        let sem_r = target_sem(&estimate, 0) / r as f64;
        let sem_2r = target_sem(&estimate, 1) / (2 * r) as f64;
        let difference = mu_2r - mu_r;
        let z = difference / (sem_r * sem_r + sem_2r * sem_2r).sqrt();
        eprintln!(
            "F2 subadditivity r={r:3} N={realizations}: E[T_r]={:.4} E[T_2r]={:.4} \
             mu_r={:.6} mu_2r={:.6} diff={:+.6} ± {:.6} z={z:+.2}",
            estimate.profile_mean[0],
            estimate.profile_mean[1],
            mu_r,
            mu_2r,
            difference,
            (sem_r * sem_r + sem_2r * sem_2r).sqrt(),
        );
        assert!(
            z < 3.0,
            "r={r}: apparent violation z={z:.2} exceeds the consistency threshold"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3c: the time-constant sequence E[T_r]/r decreases and stabilizes
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn time_constant_sequence_decreases_and_stabilizes() {
    let law = EdgeWeightLaw::exponential(1.0);
    let steps = [8_usize, 16, 32, 64, 128];
    // One fixed generous box: margin = 64 (half the largest distance),
    // half_height = 64 — every target measured in the same geometry.
    let geometry = PassageGeometry::segment(&steps, 64, 64);
    let realizations = 800_u64;
    let mut sampler = PassageSampler::new(geometry).expect("ladder sampler");
    let mut estimator = PassageEstimator::new(steps.len()).expect("ladder estimator");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_C000);
    for _ in 0..realizations {
        estimator
            .add(sampler.realize(law, &mut rng).expect("ladder run"))
            .expect("ladder add");
    }
    let estimate = estimator.estimate().expect("ladder estimate");
    let mut previous_mu = f64::INFINITY;
    let mut previous_sem = f64::NAN;
    let mut sequence = String::new();
    for (index, &r) in steps.iter().enumerate() {
        let mu = estimate.profile_mean[index] / r as f64;
        let sem = target_sem(&estimate, index) / r as f64;
        sequence.push_str(&format!("mu({r})={mu:.6}±{sem:.6} "));
        if index > 0 {
            let difference = mu - previous_mu;
            let noise = (sem * sem + previous_sem * previous_sem).sqrt();
            let z = difference / noise;
            eprintln!(
                "F2 time constant r={r:3}: mu={mu:.6}±{sem:.6} diff vs r/2={difference:+.6} \
                 (z={z:+.2})"
            );
            assert!(
                difference < 3.0 * noise,
                "r={r}: mu increased by {difference:.6} (z={z:.2}) — must decrease within noise"
            );
        }
        previous_mu = mu;
        previous_sem = sem;
    }
    eprintln!(
        "F2 time constant sequence (segment box margin 64, hh 64, N={realizations}): {sequence}"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 6: throughput probe at two sizes
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn throughput_probe_two_corridor_sizes() {
    let law = EdgeWeightLaw::exponential(1.0);
    for (r, realizations) in [(64_usize, 300_u64), (128, 150)] {
        let geometry = PassageGeometry::corridor(r, r);
        let mut sampler = PassageSampler::new(geometry.clone()).expect("throughput sampler");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_D000 + r as u64);
        // Warmup (workspace and caches).
        for _ in 0..8 {
            std::hint::black_box(sampler.realize(law, &mut rng).expect("warmup"));
        }
        let start = std::time::Instant::now();
        for _ in 0..realizations {
            std::hint::black_box(sampler.realize(law, &mut rng).expect("measured run"));
        }
        let elapsed = start.elapsed();
        eprintln!(
            "F2 throughput corridor {r}x{r} (V={} E={}): {:.0} us/realization, {:.0} \
             realizations/s [debug test profile]",
            geometry.lattice().n_sites,
            geometry.total_edge_count(),
            elapsed.as_micros() as f64 / realizations as f64,
            realizations as f64 / elapsed.as_secs_f64()
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Typed rejections
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn typed_rejections_leave_last_profile_untouched() {
    let geometry = PassageGeometry::corridor(4, 3);
    let mut sampler = PassageSampler::new(geometry.clone()).expect("sampler");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F2_E000);
    let profile = sampler
        .realize(EdgeWeightLaw::exponential(1.0), &mut rng)
        .expect("valid run")
        .to_vec();
    let total = geometry.total_edge_count();

    // Wrong weight length.
    let short_slice: Vec<f64> = sampler.weights()[..total - 1].to_vec();
    let error = sampler.realize_with_weights(&short_slice).unwrap_err();
    assert!(matches!(
        error,
        PassageObservablesError::WeightLength {
            expected: 25,
            actual: 24
        }
    ));
    assert_eq!(
        sampler.profile(),
        profile,
        "rejection leaves profile intact"
    );

    // Non-finite / negative weights are the wrapped F1 rejections.
    let mut bad = sampler.weights().to_vec();
    bad[3] = -0.25;
    let error = sampler.realize_with_weights(&bad).unwrap_err();
    assert!(matches!(
        error,
        PassageObservablesError::Passage(FirstPassageError::NegativeWeight { edge: 3, .. })
    ));
    bad[3] = f64::NAN;
    assert!(matches!(
        sampler.realize_with_weights(&bad).unwrap_err(),
        PassageObservablesError::Passage(FirstPassageError::WeightNotFinite { edge: 3, .. })
    ));

    // Invalid exponential rates, named by value (bitwise: NaN payloads).
    for rate in [0.0_f64, -1.0, f64::NAN, f64::INFINITY] {
        let error = sampler
            .realize(EdgeWeightLaw::exponential(rate), &mut rng)
            .unwrap_err();
        assert!(
            matches!(error, PassageObservablesError::InvalidRate { rate: r } if r.to_bits() == rate.to_bits()),
            "rate {rate}"
        );
    }

    // Estimator boundary mismatch adds nothing.
    let mut estimator = PassageEstimator::new(3).expect("estimator");
    estimator.add(&profile).expect("matching add");
    let error = estimator.add(&profile[..2]).unwrap_err();
    assert!(matches!(
        error,
        PassageObservablesError::ProfileLength {
            expected: 3,
            actual: 2
        }
    ));
    assert_eq!(estimator.samples(), 1);
    assert!(PassageEstimator::new(3).unwrap().estimate().is_none());
}
