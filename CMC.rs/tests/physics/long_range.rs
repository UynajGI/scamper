//! L1 long-range edge-law validation (roadmap v2 §4, stage L1).
//!
//! Gates, verbatim from the stage spec:
//!
//! 1. **Edge independence + law correctness** — the per-pair probability
//!    verified through the public draw contract (threshold replay: for
//!    each pair, recompute `min(1, beta r^{-s})` and check the occupancy
//!    threshold from both sides, bit-adjacent to p), and empirical
//!    per-pair occupancy frequencies matching the law within CI on small
//!    boxes (per-pair z-scores + pooled Pearson chi-square, 600 >= 400
//!    realizations), plus an occupancy-covariance spot check between
//!    pairs (independence);
//! 2. **Degree moments vs the LRP measure** — the closed-form expectation
//!    `E[deg(x)] = sum_{y != x} min(1, beta |x-y|^{-s})` computed in the
//!    test and verified by the empirical mean across realizations
//!    (z-score gate) for (beta, s) covering s < 1, s = 1, s > 2 in d = 1;
//! 3. **Bitwise draw contract** — same seed, identical trajectory bitwise;
//!    the adapter's draw order pinned against the documented per-pair
//!    ascending order by manual replay;
//! 4. **Trajectory correctness** — per-step parity against a test-local
//!    independent union-find (the G1 pattern) and final component
//!    structure parity against a test-local flood-fill on the occupied
//!    edge set;
//! 5. **Saturation** — beta >= 1 makes every r = 1 pair surely occupied;
//!    mid radii saturate at `beta r^{-s} >= 1`; huge beta yields the
//!    complete box graph;
//! 6. **Typed rejections** — non-finite/negative beta and s, empty box,
//!    zero sides, wrong-length and out-of-domain uniform fields;
//! 7. **Throughput report** — eprintln-timed boxes/second at two sizes.

use std::time::Instant;

use cmc_rs::{
    long_range_run, LongRangeBox, LongRangeError, LongRangePercolationMC, LongRangeWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ─────────────────────────────────────────────────────────────────────────
// Test-local references and helpers
// ─────────────────────────────────────────────────────────────────────────

/// The next f64 strictly below a non-negative `value` (bit-adjacent).
fn next_below(value: f64) -> f64 {
    assert!(value > 0.0);
    f64::from_bits(value.to_bits() - 1)
}

/// The next f64 strictly above `value` (bit-adjacent).
fn next_above(value: f64) -> f64 {
    f64::from_bits(value.to_bits() + 1)
}

/// Draw one uniform field per the documented contract: one `f64` per pair
/// in ascending pair order.
fn draw_uniforms(pairs: usize, rng: &mut Xoshiro256PlusPlus) -> Vec<f64> {
    (0..pairs).map(|_| rng.random::<f64>()).collect()
}

/// Test-local union-find over all box vertices (active isolates from the
/// start), tracking the G1 online observables independently.
struct ReferenceUnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
    component_count: usize,
    largest: usize,
    second_moment: u128,
}

impl ReferenceUnionFind {
    fn new(vertices: usize) -> Self {
        Self {
            parent: (0..vertices).collect(),
            size: vec![1; vertices],
            component_count: vertices,
            largest: 1,
            second_moment: vertices as u128,
        }
    }

    fn find(&mut self, node: usize) -> usize {
        let mut root = node;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        self.parent[node] = root;
        root
    }

    fn union(&mut self, left: usize, right: usize) {
        let (left_root, right_root) = (self.find(left), self.find(right));
        if left_root == right_root {
            return;
        }
        let (big, small) = if self.size[left_root] >= self.size[right_root] {
            (left_root, right_root)
        } else {
            (right_root, left_root)
        };
        self.second_moment += 2 * (self.size[big] as u128) * (self.size[small] as u128);
        self.size[big] += self.size[small];
        self.parent[small] = big;
        self.component_count -= 1;
        self.largest = self.largest.max(self.size[big]);
    }
}

/// Test-local flood fill over an explicit occupied edge set: returns a
/// component label per vertex (labels are the minimum vertex index of the
/// component — the F3 canonical-identity convention).
fn flood_components(vertices: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    for &(left, right) in edges {
        adjacency[left].push(right);
        adjacency[right].push(left);
    }
    let mut labels = vec![usize::MAX; vertices];
    for start in 0..vertices {
        if labels[start] != usize::MAX {
            continue;
        }
        labels[start] = start;
        let mut stack = vec![start];
        while let Some(vertex) = stack.pop() {
            for &neighbor in &adjacency[vertex] {
                if labels[neighbor] == usize::MAX {
                    labels[neighbor] = start;
                    stack.push(neighbor);
                }
            }
        }
    }
    labels
}

/// One-sided chi-square critical value via the Wilson-Hilferty
/// approximation at z = 3.09 (upper tail alpha ~ 0.001).
fn chi2_critical_999(degrees: usize) -> f64 {
    let df = degrees as f64;
    let scale = 2.0 / (9.0 * df);
    df * (1.0 - scale + 3.09 * scale.sqrt()).powi(3)
}

/// Every pair index of a box with its endpoints (ascending pair order).
fn all_pairs(box_: &LongRangeBox) -> Vec<(usize, usize, usize)> {
    (0..box_.pairs())
        .map(|pair| {
            let (first, second) = box_.pair_endpoints(pair).expect("dense pair index");
            (pair, first, second)
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1a: threshold replay — the law and the pair order through the
// public draw contract
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn threshold_replay_pins_the_law_and_pair_order() {
    let box_ = LongRangeBox::new(&[8]).expect("1x8 box");
    let pairs = all_pairs(&box_);
    assert_eq!(pairs.len(), 28);
    // The dense order is the coordinate order in d=1; the documented pair
    // order is lexicographic (x, y) — pin the enumeration itself.
    for (pair, first, second) in &pairs {
        assert!(first < second);
        assert_eq!(
            box_.pair_index(*first, *second),
            Some(*pair),
            "pair_index/pair_endpoints must be inverse"
        );
        // ℓ1 distance == coordinate difference in d=1 (the norm-choice
        // promise for the stage-L1 target dimension).
        assert_eq!(box_.distance(*first, *second), second - first);
    }
    for (beta, s) in [(0.8, 0.5), (0.3, 1.0), (2.5, 1.7)] {
        let mut workspace = LongRangeWorkspace::new();
        for (pair, first, second) in &pairs {
            let r = box_.distance(*first, *second);
            let p = LongRangeBox::pair_probability(r, beta, s);
            // The contract expression agrees with the direct power-law
            // reading `beta / r^s` to rounding.
            let direct = (beta / (r as f64).powf(s)).min(1.0);
            assert!(
                ((p - direct) / direct).abs() < 1e-12,
                "p({r}) = {p} vs direct {direct}"
            );

            // Every other pair gets u = 1.0, which never occupies
            // (u < min(1, .) is false), so the occupied set isolates pair k.
            // Below the threshold (bit-adjacent): occupied.
            let mut uniforms = vec![1.0_f64; box_.pairs()];
            uniforms[*pair] = if p == 1.0 {
                next_below(1.0)
            } else {
                next_below(p)
            };
            let outcome = long_range_run(&box_, beta, s, &uniforms, &mut workspace)
                .expect("replay below threshold");
            assert_eq!(
                workspace.occupied_pairs(),
                &[*pair][..],
                "pair {pair} at r={r}: u just below p must occupy (and only it)"
            );
            assert_eq!(outcome.added_edges, 1);

            // At/above the threshold: not occupied. u = 1.0 exactly is the
            // saturation boundary (p = 1 still rejects it: strict <).
            let mut uniforms = vec![1.0_f64; box_.pairs()];
            uniforms[*pair] = 1.0;
            let outcome = long_range_run(&box_, beta, s, &uniforms, &mut workspace)
                .expect("replay at saturation boundary");
            assert!(
                workspace.occupied_pairs().is_empty(),
                "u = 1.0 must never occupy, even at p = 1"
            );
            assert_eq!(outcome.added_edges, 0);
            if p < 1.0 {
                let above = next_above(p);
                assert!(above <= 1.0);
                let mut uniforms = vec![1.0_f64; box_.pairs()];
                uniforms[*pair] = above;
                long_range_run(&box_, beta, s, &uniforms, &mut workspace)
                    .expect("replay above threshold");
                assert!(
                    workspace.occupied_pairs().is_empty(),
                    "pair {pair} at r={r}: u just above p must not occupy"
                );
            } else {
                // p == 1: verify saturation itself (beta r^{-s} >= 1).
                assert!(beta * (r as f64).powf(-s) >= 1.0);
            }
        }
    }
    eprintln!("long-range threshold replay: 28 pairs x 3 combos x both sides ok");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1b: per-pair occupancy frequencies vs the law (z + chi-square) and
// an independence spot check
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn per_pair_occupancy_frequencies_match_the_law() {
    struct Config {
        sides: Vec<usize>,
        beta: f64,
        s: f64,
    }
    let configs = [
        Config {
            sides: vec![8],
            beta: 0.8,
            s: 0.5,
        },
        Config {
            sides: vec![8],
            beta: 0.3,
            s: 1.0,
        },
        Config {
            sides: vec![8],
            beta: 2.5,
            s: 1.7,
        },
        // A 2D box exercises the multi-dimensional ℓ1 distance formula.
        Config {
            sides: vec![4, 4],
            beta: 0.5,
            s: 1.2,
        },
    ];
    let realizations = 600_usize;
    for config in &configs {
        let box_ = LongRangeBox::new(&config.sides).expect("small box");
        let pairs = all_pairs(&box_);
        let mut counts = vec![0_usize; box_.pairs()];
        let mut process =
            LongRangePercolationMC::new(box_.clone(), config.beta, config.s, 8).expect("adapter");
        // Pair-of-pairs occupancy products for the independence check.
        let probes: [(usize, usize); 3] = [
            (
                box_.pair_index(0, 1).unwrap(),
                box_.pair_index(2, 5).unwrap(),
            ),
            (
                box_.pair_index(0, box_.vertices() - 1).unwrap(),
                box_.pair_index(3, 4).unwrap(),
            ),
            (
                box_.pair_index(1, 2).unwrap(),
                box_.pair_index(5, 6).unwrap(),
            ),
        ];
        let mut joint = vec![0_usize; probes.len()];
        for seed in 0..realizations as u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_3131 + seed);
            process.realize(&mut rng);
            let occupied = process.workspace().occupied_pairs();
            for &pair in occupied {
                counts[pair] += 1;
            }
            for (position, &(left, right)) in probes.iter().enumerate() {
                if occupied.binary_search(&left).is_ok() && occupied.binary_search(&right).is_ok() {
                    joint[position] += 1;
                }
            }
        }

        let mut max_z = 0.0_f64;
        let mut chi2 = 0.0_f64;
        let mut degrees = 0_usize;
        let mut saturated = 0_usize;
        for (pair, first, second) in &pairs {
            let p = LongRangeBox::pair_probability(
                box_.distance(*first, *second),
                config.beta,
                config.s,
            );
            if p >= 1.0 {
                saturated += 1;
                assert_eq!(
                    counts[*pair], realizations,
                    "saturated pair {pair} must ALWAYS be occupied"
                );
                continue;
            }
            let expected = realizations as f64 * p;
            let variance = realizations as f64 * p * (1.0 - p);
            let z = (counts[*pair] as f64 - expected) / variance.sqrt();
            assert!(
                z.abs() < 5.0,
                "pair {pair} (r={}): frequency {} vs p={p:.4}, z={z:.2}",
                box_.distance(*first, *second),
                counts[*pair] as f64 / realizations as f64
            );
            max_z = max_z.max(z.abs());
            chi2 += (counts[*pair] as f64 - expected).powi(2) / expected
                + ((realizations - counts[*pair]) as f64 - realizations as f64 * (1.0 - p)).powi(2)
                    / (realizations as f64 * (1.0 - p));
            degrees += 1;
        }
        let critical = chi2_critical_999(degrees);
        assert!(
            chi2 < critical,
            "chi-square {chi2:.1} vs critical {critical:.1} at df={degrees}"
        );
        eprintln!(
            "long-range law frequencies sides={:?} beta={} s={}: pairs={} ({} saturated) \
             realizations={} max|z|={max_z:.2} chi2={chi2:.1} df={degrees} critical={critical:.1}",
            config.sides,
            config.beta,
            config.s,
            pairs.len(),
            saturated,
            realizations
        );

        // Independence spot check: the joint occupancy rate of two pairs
        // matches the product law within CI (z against the exact product
        // q = p_left * p_right, whose Bernoulli-product variance is
        // q (1 - q) / R).
        for (position, &(left, right)) in probes.iter().enumerate() {
            let (p_left, p_right) = (
                {
                    let (a, b) = box_.pair_endpoints(left).unwrap();
                    LongRangeBox::pair_probability(box_.distance(a, b), config.beta, config.s)
                },
                {
                    let (a, b) = box_.pair_endpoints(right).unwrap();
                    LongRangeBox::pair_probability(box_.distance(a, b), config.beta, config.s)
                },
            );
            let q = p_left * p_right;
            let joint_rate = joint[position] as f64 / realizations as f64;
            if q >= 1.0 {
                // Both probe pairs saturated: the joint occupancy is sure.
                assert_eq!(joint[position], realizations, "saturated probe pairs");
                continue;
            }
            let stderr = (q * (1.0 - q) / realizations as f64).sqrt();
            let z = (joint_rate - q) / stderr;
            assert!(
                z.abs() < 4.0,
                "pairs {left}/{right}: joint rate {joint_rate:.5} vs q={q:.5}, z={z:.2}"
            );
            eprintln!(
                "long-range independence sides={:?}: pairs {left}/{right} joint={joint_rate:.4} \
                 q={q:.4} z={z:.2}",
                config.sides
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: center degree moments vs the closed-form LRP expectation
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn center_degree_matches_closed_form_expectation() {
    let length = 16_usize;
    let box_ = LongRangeBox::new(&[length]).expect("1x16 box");
    let center = 8_usize;
    // The closed form per box: E[deg(x)] = sum over y != x of
    // min(1, beta |x-y|^{-s}) — summed explicitly, and cross-checked
    // against the hand formula for the d=1 center (r = 1..8 left,
    // r = 1..7 right).
    let realizations = 4000_usize;
    for (beta, s) in [(0.8, 0.5), (0.3, 1.0), (0.25, 2.6)] {
        let mut expected = 0.0_f64;
        for other in 0..length {
            if other != center {
                expected += LongRangeBox::pair_probability(box_.distance(center, other), beta, s);
            }
        }
        let hand: f64 = (1..=8)
            .map(|r| LongRangeBox::pair_probability(r, beta, s))
            .sum::<f64>()
            + (1..=7)
                .map(|r| LongRangeBox::pair_probability(r, beta, s))
                .sum::<f64>();
        assert!((expected - hand).abs() < 1e-12, "{expected} vs {hand}");

        let mut process = LongRangePercolationMC::new(box_.clone(), beta, s, 1).expect("adapter");
        let mut sum = 0.0_f64;
        let mut sum_squares = 0.0_f64;
        for seed in 0..realizations as u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_4450 + seed);
            process.realize(&mut rng);
            let mut degree = 0_usize;
            for &pair in process.workspace().occupied_pairs() {
                let (first, second) = box_.pair_endpoints(pair).expect("dense pair");
                if first == center || second == center {
                    degree += 1;
                }
            }
            sum += degree as f64;
            sum_squares += (degree * degree) as f64;
        }
        let mean = sum / realizations as f64;
        let variance = (sum_squares - sum * sum / realizations as f64) / (realizations - 1) as f64;
        let z = (mean - expected) / (variance / realizations as f64).sqrt();
        assert!(
            z.abs() < 4.0,
            "beta={beta} s={s}: mean degree {mean:.4} vs E={expected:.4}, z={z:.2}"
        );
        eprintln!(
            "long-range degree moment beta={beta} s={s}: E[deg]={expected:.4} mean={mean:.4} \
             var={variance:.3} R={realizations} z={z:.2}"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: bitwise draw contract — same seed identical, adapter order pinned
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn bitwise_draw_contract_and_adapter_order() {
    let box_ = LongRangeBox::new(&[8]).expect("1x8 box");
    let (beta, s) = (0.7_f64, 1.3_f64);
    let steps = 12_usize;

    // Same seed twice: bit-identical history, occupied pairs, and outcome.
    let first = {
        let mut process = LongRangePercolationMC::new(box_.clone(), beta, s, steps).unwrap();
        process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_4231));
        (
            process.last_history().to_vec(),
            process.workspace().occupied_pairs().to_vec(),
            process.last_outcome(),
            process
                .workspace()
                .substrate()
                .trajectory()
                .summaries()
                .to_vec(),
            process
                .workspace()
                .substrate()
                .trajectory()
                .steps()
                .to_vec(),
        )
    };
    let second = {
        let mut process = LongRangePercolationMC::new(box_.clone(), beta, s, steps).unwrap();
        process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_4231));
        (
            process.last_history().to_vec(),
            process.workspace().occupied_pairs().to_vec(),
            process.last_outcome(),
            process
                .workspace()
                .substrate()
                .trajectory()
                .summaries()
                .to_vec(),
            process
                .workspace()
                .substrate()
                .trajectory()
                .steps()
                .to_vec(),
        )
    };
    for (left, right) in first.3.iter().zip(&second.3) {
        assert_eq!(left, right);
    }
    assert_eq!(first.0, second.0, "bit-identical history");
    assert_eq!(first.1, second.1, "bit-identical occupied pairs");
    assert_eq!(first.2, second.2);
    assert_eq!(first.4, second.4);

    // The adapter's draw order is the documented one: manually draw one
    // f64 per pair in ascending order from an equally seeded stream and
    // run the RNG-free core — everything must match bitwise.
    let mut manual_rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_4231);
    let uniforms = draw_uniforms(box_.pairs(), &mut manual_rng);
    let mut workspace = LongRangeWorkspace::new();
    long_range_run(&box_, beta, s, &uniforms, &mut workspace).expect("manual replay");
    assert_eq!(workspace.occupied_pairs(), &first.1[..]);
    assert_eq!(workspace.substrate().trajectory().steps(), &first.4[..]);
    assert_eq!(workspace.substrate().trajectory().summaries(), &first.3[..]);

    // Different seeds realize different boxes (deterministic for these
    // fixed seeds; a collision would fail loudly).
    let mut process = LongRangePercolationMC::new(box_, beta, s, steps).unwrap();
    process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_4232));
    assert_ne!(process.workspace().occupied_pairs(), &first.1[..]);

    // C0 padding convention (G3 shape rule): entry i is the largest size
    // after the (i+1)-th added edge; the tail repeats the final value.
    let added = first.2.added_edges;
    let summaries = &first.3;
    let base = 8_usize;
    for (step, &value) in first.0.iter().enumerate() {
        let expected = if step < added {
            summaries[base + step].largest_component_size as f64
        } else {
            first.2.largest_component_size as f64
        };
        assert_eq!(value.to_bits(), expected.to_bits(), "history entry {step}");
    }
    eprintln!(
        "long-range bitwise contract: {added} added edges, history len {} (padded), ok",
        first.0.len()
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: per-step union-find parity + final flood-fill parity
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn trajectory_and_component_parity_against_references() {
    let configs: Vec<(Vec<usize>, f64, f64)> = vec![
        (vec![8], 0.8, 0.5),
        (vec![8], 0.3, 1.0),
        (vec![13], 2.5, 1.7),
        (vec![13], 1.5, 0.9),
        (vec![4, 4], 0.5, 1.2),
        (vec![3, 5], 1.1, 2.6),
    ];
    let mut parity_cases = 0_usize;
    let mut structure_cases = 0_usize;
    for (sides, beta, s) in &configs {
        let box_ = LongRangeBox::new(sides).expect("box");
        let vertices = box_.vertices();
        for seed in 0..5_u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_5041 + seed);
            let uniforms = draw_uniforms(box_.pairs(), &mut rng);
            let mut workspace = LongRangeWorkspace::new();
            let outcome =
                long_range_run(&box_, *beta, *s, &uniforms, &mut workspace).expect("parity run");

            // Independent replay through the public geometry API: decide
            // occupancy per pair in consideration order, drive a
            // test-local union-find, compare every trajectory entry.
            let trajectory = workspace.substrate().trajectory();
            let summaries = trajectory.summaries();
            let steps = trajectory.steps();
            let mut reference = ReferenceUnionFind::new(vertices);
            let mut expected_steps = vec![0_usize; vertices];
            let mut occupied: Vec<(usize, usize)> = Vec::new();
            let mut added = 0_usize;
            for (pair, &uniform) in uniforms.iter().enumerate() {
                let (first, second) = box_.pair_endpoints(pair).expect("dense pair");
                let p = LongRangeBox::pair_probability(box_.distance(first, second), *beta, *s);
                if uniform < p {
                    added += 1;
                    expected_steps.push(added);
                    reference.union(first, second);
                    occupied.push((first, second));
                    let summary = summaries[vertices + added - 1];
                    assert_eq!(summary.active_vertex_count, vertices);
                    assert_eq!(summary.active_edge_count, added);
                    assert_eq!(summary.component_count, reference.component_count);
                    assert_eq!(summary.largest_component_size, reference.largest);
                    assert_eq!(summary.raw_second_moment, reference.second_moment);
                    parity_cases += 1;
                }
            }
            assert_eq!(steps, expected_steps.as_slice(), "steps = [0;n]+[1..=m]");
            assert_eq!(
                workspace.occupied_pairs(),
                occupied
                    .iter()
                    .map(|&(first, second)| box_.pair_index(first, second).unwrap())
                    .collect::<Vec<_>>()
                    .as_slice()
            );
            // The vertex stage: j+1 isolates after activation j.
            for (index, summary) in summaries[..vertices].iter().enumerate() {
                assert_eq!(summary.active_vertex_count, index + 1);
                assert_eq!(summary.component_count, index + 1);
                assert_eq!(summary.largest_component_size, 1);
                assert_eq!(summary.raw_second_moment, (index + 1) as u128);
                assert_eq!(summary.active_edge_count, 0);
            }

            // Final structure parity vs a test-local flood fill on the
            // occupied edge set.
            let labels = flood_components(vertices, &occupied);
            let mut sizes = vec![0_usize; vertices];
            for &label in &labels {
                sizes[label] += 1;
            }
            let component_count = labels
                .iter()
                .collect::<std::collections::HashSet<_>>()
                .len();
            let largest = *sizes.iter().max().unwrap();
            let second_moment: u128 = sizes.iter().map(|&size| size as u128 * size as u128).sum();
            assert_eq!(outcome.added_edges, added);
            assert_eq!(outcome.component_count, component_count);
            assert_eq!(outcome.largest_component_size, largest);
            assert_eq!(outcome.raw_second_moment, second_moment);
            for vertex in 0..vertices {
                assert_eq!(
                    workspace.component_size(vertex),
                    Some(sizes[labels[vertex]]),
                    "component size of vertex {vertex}"
                );
            }
            // The canonical identity: lowest dense index of a largest
            // component (F3 tie rule).
            let identity = labels
                .iter()
                .enumerate()
                .filter(|&(_, &label)| sizes[label] == largest)
                .map(|(vertex, _)| vertex)
                .min()
                .unwrap();
            assert_eq!(
                summaries
                    .last()
                    .expect("vertices are active")
                    .largest_component,
                Some(identity)
            );
            structure_cases += 1;
        }
    }
    eprintln!(
        "long-range parity: {parity_cases} added-edge steps vs independent union-find; \
         {structure_cases} final structures vs flood fill"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5: saturation regimes
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn saturation_regimes_are_sure_occupied() {
    // beta >= 1 saturates every r = 1 pair: nearest neighbors are ALWAYS
    // present, across seeds and exponents.
    let box_ = LongRangeBox::new(&[12]).expect("1x12 box");
    let neighbor_pairs: Vec<usize> = (0..11)
        .map(|index| box_.pair_index(index, index + 1).unwrap())
        .collect();
    for (beta, s) in [(1.0, 2.6), (1.5, 0.9), (40.0, 3.0)] {
        let mut process = LongRangePercolationMC::new(box_.clone(), beta, s, 4).unwrap();
        for seed in 0..40_u64 {
            process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_5341 + seed));
            for &pair in &neighbor_pairs {
                assert!(
                    process
                        .workspace()
                        .occupied_pairs()
                        .binary_search(&pair)
                        .is_ok(),
                    "beta={beta} s={s} seed={seed}: neighbor pair {pair} missing"
                );
            }
        }
    }

    // Mid-radius saturation: beta = 4, s = 1 saturates every pair with
    // r <= 4 (beta r^{-s} >= 1) — always present; pairs beyond stay
    // probabilistic.
    let (beta, s) = (4.0_f64, 1.0_f64);
    let saturated: Vec<usize> = (0..box_.pairs())
        .filter(|&pair| {
            let (first, second) = box_.pair_endpoints(pair).unwrap();
            LongRangeBox::pair_probability(box_.distance(first, second), beta, s) >= 1.0
        })
        .collect();
    assert_eq!(saturated.len(), 11 + 10 + 9 + 8, "r <= 1..4 shells");
    let mut process = LongRangePercolationMC::new(box_.clone(), beta, s, 4).unwrap();
    for seed in 0..40_u64 {
        process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_5342 + seed));
        for &pair in &saturated {
            assert!(
                process
                    .workspace()
                    .occupied_pairs()
                    .binary_search(&pair)
                    .is_ok(),
                "seed={seed}: saturated pair {pair} missing"
            );
        }
    }

    // Huge beta: the complete box graph, every realization.
    let small = LongRangeBox::new(&[8]).expect("1x8 box");
    let mut process =
        LongRangePercolationMC::new(small.clone(), 1.0e6, 2.0, 4).expect("complete adapter");
    for seed in 0..8_u64 {
        process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_5343 + seed));
        let outcome = process.last_outcome();
        assert_eq!(outcome.added_edges, small.pairs());
        assert_eq!(outcome.component_count, 1);
        assert_eq!(outcome.largest_component_size, small.vertices());
        assert_eq!(outcome.raw_second_moment, 64);
    }

    // Degenerate laws stay degenerate: beta = 0 adds nothing; a
    // single-vertex box has an empty pair universe.
    let mut empty = LongRangePercolationMC::new(box_.clone(), 0.0, 1.0, 4).unwrap();
    empty.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_5344));
    assert_eq!(empty.last_outcome().added_edges, 0);
    assert_eq!(empty.last_outcome().component_count, 12);
    assert_eq!(empty.last_history(), &[1.0; 4]);
    let singleton = LongRangeBox::new(&[1]).expect("1x1 box");
    assert_eq!(singleton.pairs(), 0);
    let mut process = LongRangePercolationMC::new(singleton, 0.9, 0.5, 3).unwrap();
    process.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x004C_5345));
    assert_eq!(process.last_outcome().added_edges, 0);
    assert_eq!(process.last_outcome().largest_component_size, 1);
    eprintln!(
        "long-range saturation: r=1 always at beta>=1 (3 laws x 40 seeds); \
               38 saturated pairs at beta=4/s=1 (40 seeds); complete graph at beta=1e6"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 6: typed rejections
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn typed_rejections_leave_state_untouched() {
    assert!(matches!(
        LongRangeBox::new(&[]),
        Err(LongRangeError::EmptyBox)
    ));
    assert!(matches!(
        LongRangeBox::new(&[4, 0, 3]),
        Err(LongRangeError::ZeroSide { axis: 1 })
    ));

    let box_ = LongRangeBox::new(&[6]).expect("1x6 box");
    let mut workspace = LongRangeWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
    let uniforms = draw_uniforms(box_.pairs(), &mut rng);
    long_range_run(&box_, 0.8, 0.5, &uniforms, &mut workspace).expect("valid run");
    let before = (
        workspace.occupied_pairs().to_vec(),
        workspace.substrate().summary(),
    );

    for beta in [-1.0, f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert!(matches!(
            long_range_run(&box_, beta, 0.5, &uniforms, &mut workspace),
            Err(LongRangeError::BetaInvalid { value }) if value.to_bits() == beta.to_bits()
        ));
    }
    for s in [-0.5, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            long_range_run(&box_, 0.8, s, &uniforms, &mut workspace),
            Err(LongRangeError::ExponentInvalid { .. })
        ));
    }
    assert!(matches!(
        long_range_run(
            &box_,
            0.8,
            0.5,
            &uniforms[..box_.pairs() - 1],
            &mut workspace
        ),
        Err(LongRangeError::UniformLengthMismatch {
            expected: 15,
            actual: 14
        })
    ));
    let mut bad = uniforms.clone();
    bad[3] = 1.5;
    assert!(matches!(
        long_range_run(&box_, 0.8, 0.5, &bad, &mut workspace),
        Err(LongRangeError::UniformInvalid { index: 3, .. })
    ));
    bad[3] = -0.25;
    assert!(matches!(
        long_range_run(&box_, 0.8, 0.5, &bad, &mut workspace),
        Err(LongRangeError::UniformInvalid { index: 3, .. })
    ));
    bad[3] = f64::NAN;
    assert!(matches!(
        long_range_run(&box_, 0.8, 0.5, &bad, &mut workspace),
        Err(LongRangeError::UniformInvalid { index: 3, .. })
    ));
    assert!(matches!(
        LongRangePercolationMC::new(box_.clone(), 0.8, 0.5, 0),
        Err(LongRangeError::NoSteps)
    ));
    assert!(matches!(
        LongRangePercolationMC::new(box_.clone(), -0.8, 0.5, 8),
        Err(LongRangeError::BetaInvalid { .. })
    ));

    // Every rejection left the last valid realization untouched.
    assert_eq!(workspace.occupied_pairs(), &before.0[..]);
    assert_eq!(workspace.substrate().summary(), before.1);
    eprintln!("long-range typed rejections: 13 refusal cases, state untouched");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 7: throughput report (eprintln-timed, debug profile)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn throughput_report_two_sizes() {
    for length in [256_usize, 512] {
        let box_ = LongRangeBox::new(&[length]).expect("1D box");
        let mut process =
            LongRangePercolationMC::new(box_.clone(), 0.8, 1.0, 64).expect("throughput adapter");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004C_5431);
        for _ in 0..4 {
            process.realize(&mut rng);
        }
        let realizations = 64_usize;
        let start = Instant::now();
        for _ in 0..realizations {
            process.realize(&mut rng);
        }
        let elapsed = start.elapsed();
        let per_box = elapsed.as_secs_f64() / realizations as f64;
        eprintln!(
            "LONG_RANGE_THROUGHPUT box=1x{length} V={V} pairs={} realizations={realizations} \
             total={elapsed:?} per_box={per_box_us:.1}us boxes_per_second={bps:.0} \
             pairs_per_second={pps:.3e}",
            box_.pairs(),
            V = length,
            per_box_us = per_box * 1.0e6,
            bps = 1.0 / per_box,
            pps = box_.pairs() as f64 / per_box,
        );
    }
}

#[test]
fn coordinates_round_trip_across_dimensions_and_domains() {
    // Review follow-up (M1): coordinates() was axis-reversed for d >= 2
    // (the mixed-radix digits were peeled with the sides reversed). The
    // gate pins the inverse: every dense index decomposes and recomposes
    // to itself, every decomposition is inside the box, and the domain
    // edge rejects — across 1D/2D/3D boxes, symmetric and asymmetric.
    for sides in [&[16_usize][..], &[3, 4], &[4, 4], &[2, 3, 4], &[5, 1, 2]] {
        let box_ = LongRangeBox::new(sides).expect("valid sides");
        let mut hand_checked = 0;
        for index in 0..box_.vertices() {
            let coords = box_.coordinates(index).expect("in-domain index");
            assert_eq!(
                coords.len(),
                sides.len(),
                "dimension of coordinates({index}) on {sides:?}"
            );
            for (&coordinate, &side) in coords.iter().zip(sides.iter()) {
                assert!(
                    coordinate < side,
                    "coordinates({index}) = {coords:?} outside {sides:?}"
                );
            }
            assert_eq!(
                box_.vertex_index(&coords),
                Some(index),
                "round trip failed for index {index} on {sides:?}: {coords:?}"
            );
            hand_checked += 1;
        }
        assert!(
            box_.coordinates(box_.vertices()).is_none(),
            "domain edge on {sides:?}"
        );
        assert!(
            box_.vertex_index(&[sides[0], 0]).is_none(),
            "out-of-range or wrong-length coordinate must reject on {sides:?}"
        );
        // Spot pins from the review: [3,4] index 1 -> [0,1], index 3 -> [0,3].
        if sides == [3, 4] {
            assert_eq!(box_.coordinates(1).unwrap(), vec![0, 1]);
            assert_eq!(box_.coordinates(3).unwrap(), vec![0, 3]);
            assert_eq!(box_.coordinates(4).unwrap(), vec![1, 0]);
        }
        eprintln!("coordinates round trip: {sides:?} x {hand_checked} indices ok");
    }
}
