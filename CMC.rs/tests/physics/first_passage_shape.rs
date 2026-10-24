//! F3 shape-theorem and geodesic validation (roadmap v2 §3 F3).
//!
//! Anchor: Auffinger-Damron-Hanson, *50 years of first-passage
//! percolation*, AMS ULS 68 (2017) (archived). The survey's statements are
//! quoted in the module documentation of
//! `percolation::first_passage_shape`; this file gates them as finite-box
//! statistical brackets only, never as equalities.
//!
//! Gates, verbatim from the stage brief:
//!
//! 1. **Directional time constants** `mu(u)` on open Z^2 boxes for
//!    `(1,0)`, `(0,1)`, `(1,1)` and `(1,-1)` over a distance ladder:
//!    positive finite `mu`, `T(nu)/n` non-increasing within noise
//!    (subadditivity gives convergence from above), a finite-ladder
//!    stabilization band, rotational symmetry of the two axis
//!    directions, and the paired shape inequality
//!    `mu(u) + mu(v) >= mu(u+v)` within CI. NO closed-form shape claim.
//! 2. **Geodesic counting** through the tight-edge predecessor DAG:
//!    exact vs brute-force simple-path enumeration on small graphs at
//!    two-level dyadic weights, exact binomial counts for uniform
//!    weights on larger boxes (counts reported), and a reconstructed
//!    path whose weight sum is the F1 distance bitwise.
//! 3. **Exploratory aging statistic** (the survey has no crisp
//!    finite-horizon aging law; Question 30 is the closest anchor):
//!    correctness is gated (every selected path is a true simple
//!    geodesic, persistence in [0,1]) and convergence is reported with a
//!    deliberately generous exploratory band; limitations are explicit.
//! 4. **Exact small-graph shape checks**: uniform weights reproduce the
//!    l1 ball exactly; two-level weights match an independent Bellman-
//!    Ford / brute-force count.
//! 5. **Typed rejections with state isolation.**

use cmc_rs::{
    build_square, first_passage_run, geodesic_aging, geodesic_count_run, geodesic_path, Bond,
    BondType, CsrLattice, Direction, EdgeWeightLaw, FirstPassageError, FirstPassageWorkspace,
    GeodesicCountWorkspace, GeodesicTracking, ShapeError, ShapeEstimator, ShapeGeometry,
    ShapeSampler, SHAPE_SIDE_CAP,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ─────────────────────────────────────────────────────────────────────────
// Test-local statistics and references
// ─────────────────────────────────────────────────────────────────────────

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

/// Exact binomial coefficient in `u128` (C(n, k)).
fn binomial(n: usize, k: usize) -> u128 {
    let mut result = 1_u128;
    for step in 0..k {
        result = result * (n - step) as u128 / (step + 1) as u128;
    }
    result
}

/// Brute-force geodesic count: enumerate all simple paths from `source`
/// to `target`, pruning prefixes whose weight already exceeds the target
/// distance, and count the paths whose total equals it exactly. Dyadic
/// weights make the sums exact. Shares no code with the DAG method.
fn brute_force_geodesics(
    lattice: &CsrLattice,
    weights: &[f64],
    source: usize,
    target: usize,
    distance: f64,
) -> u128 {
    let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); lattice.n_sites];
    for (edge, bond) in lattice.edges.iter().enumerate() {
        adjacency[bond.source].push((bond.target, weights[edge]));
        adjacency[bond.target].push((bond.source, weights[edge]));
    }
    fn dfs(
        vertex: usize,
        target: usize,
        sum: f64,
        distance: f64,
        adjacency: &[Vec<(usize, f64)>],
        visited: &mut [bool],
    ) -> u128 {
        if vertex == target {
            return u128::from(sum == distance);
        }
        let mut count = 0_u128;
        for &(neighbor, weight) in &adjacency[vertex] {
            if visited[neighbor] {
                continue;
            }
            let next = sum + weight;
            if next > distance {
                continue;
            }
            visited[neighbor] = true;
            count += dfs(neighbor, target, next, distance, adjacency, visited);
            visited[neighbor] = false;
        }
        count
    }
    let mut visited = vec![false; lattice.n_sites];
    visited[source] = true;
    dfs(source, target, 0.0, distance, &adjacency, &mut visited)
}

/// Sum of a path's physical-edge weights in traversal order.
fn path_weight(edges: &[usize], weights: &[f64]) -> f64 {
    edges.iter().map(|&edge| weights[edge]).sum()
}

/// Deterministic two-level weight field: 0.5 or 1.0 per edge.
fn two_level_weights(edges: usize, seed: u64) -> Vec<f64> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..edges)
        .map(|_| if rng.random::<f64>() < 0.5 { 0.5 } else { 1.0 })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────
// Uniform shape: exact l1 ball and exact binomial geodesic counts
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn uniform_weights_reproduce_the_l1_ball_exactly() {
    // Open 6x6 box, all weights 1.0: T(0, (x, y)) = |x| + |y| exactly.
    let side = 6_usize;
    let lattice = build_square(side, side, false);
    let weights = vec![1.0_f64; lattice.n_edges()];
    let mut workspace = FirstPassageWorkspace::new();
    first_passage_run(&lattice, &weights, 0, GeodesicTracking::Off, &mut workspace)
        .expect("uniform pass");
    for vertex in 0..lattice.n_sites {
        let x = vertex % side;
        let y = vertex / side;
        assert_eq!(
            workspace.distances()[vertex],
            (x + y) as f64,
            "uniform weights: T(0, {vertex}) must equal the l1 norm"
        );
    }

    // The same check through the ShapeGeometry C0 profile, all four
    // directions at three distances.
    let geometry = ShapeGeometry::new(
        &[
            Direction::East,
            Direction::North,
            Direction::Northeast,
            Direction::Southeast,
        ],
        &[1, 2, 4],
        2,
    )
    .expect("geometry");
    let mut sampler = ShapeSampler::new(geometry.clone()).expect("sampler");
    let weights = vec![1.0_f64; geometry.edge_count()];
    let profile = sampler
        .realize_with_weights(&weights, GeodesicTracking::Off)
        .expect("uniform profile");
    for (index, target) in geometry.targets().iter().enumerate() {
        let (dx, dy) = target.direction.step();
        let expected = (dx.abs() + dy.abs()) as f64 * target.distance as f64;
        assert_eq!(
            profile[index], expected,
            "target {index} ({:?}, n = {}) must be the l1 distance",
            target.direction, target.distance
        );
    }
    eprintln!(
        "F3 uniform shape: 36 vertex distances equal l1; 12 directional targets exact on 4 \
         directions x 3 distances"
    );

    // Geodesic counts at uniform weights are binomial: corner to corner
    // on an n x n grid has C(2(n-1), n-1) monotone geodesics.
    let mut count_workspace = GeodesicCountWorkspace::new();
    for side in [6_usize, 10, 12] {
        let lattice = build_square(side, side, false);
        let weights = vec![1.0_f64; lattice.n_edges()];
        let source = 0;
        let target = side * side - 1;
        let outcome = geodesic_count_run(&lattice, &weights, source, target, &mut count_workspace)
            .expect("uniform count");
        let expected = binomial(2 * (side - 1), side - 1);
        assert_eq!(
            outcome.geodesic_count, expected,
            "uniform {side}x{side} corner-to-corner geodesic count"
        );
        assert!(!outcome.saturated);
        // An interior target: C(x + y, x).
        let (x, y) = (side - 1, side - 2);
        let interior = x + y * side;
        let interior_outcome =
            geodesic_count_run(&lattice, &weights, source, interior, &mut count_workspace)
                .expect("interior count");
        assert_eq!(interior_outcome.geodesic_count, binomial(x + y, x));
        eprintln!(
            "F3 uniform counts {side}x{side}: corner C({}, {}) = {}, interior C({}, {}) = {}",
            2 * (side - 1),
            side - 1,
            expected,
            x + y,
            x,
            binomial(x + y, x)
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Two-level weights: exact DAG counts vs brute-force enumeration
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn two_level_geodesic_counts_match_brute_force() {
    let cases = [(3_usize, 3_usize), (4, 3)];
    let mut total_graphs = 0_usize;
    let mut total_count = 0_u128;
    let mut tied_cases = 0_usize;
    for &(width, height) in &cases {
        let lattice = build_square(width, height, false);
        let source = 0;
        let target = width * height - 1;
        for seed in 0..8_u64 {
            let weights = two_level_weights(lattice.n_edges(), 0x00F3_0000 + seed);
            let mut passage = FirstPassageWorkspace::new();
            first_passage_run(
                &lattice,
                &weights,
                source,
                GeodesicTracking::Predecessors,
                &mut passage,
            )
            .expect("passage pass");
            let distance = passage.distances()[target];
            assert!(distance.is_finite());
            let expected = brute_force_geodesics(&lattice, &weights, source, target, distance);
            let mut workspace = GeodesicCountWorkspace::new();
            let outcome = geodesic_count_run(&lattice, &weights, source, target, &mut workspace)
                .expect("DAG count");
            assert_eq!(
                outcome.geodesic_count, expected,
                "{width}x{height} seed {seed}: DAG count must equal brute-force simple-path \
                 enumeration"
            );
            if expected > 1 {
                tied_cases += 1;
            }
            total_count += expected;
            total_graphs += 1;

            // The reconstructed predecessor path is a true geodesic: its
            // weight sum equals the F1 distance bitwise, it is simple, and
            // it uses exactly `vertices - 1` edges.
            let mut path = Vec::new();
            let mut edges = Vec::new();
            let path_outcome = geodesic_path(
                passage.predecessor_vertices(),
                passage.predecessor_edges(),
                source,
                target,
                &mut path,
                &mut edges,
            )
            .expect("path reconstruction");
            assert_eq!(path_outcome.vertices, path.len());
            assert_eq!(path_outcome.edges, path.len() - 1);
            assert_eq!(path.first().copied(), Some(source));
            assert_eq!(path.last().copied(), Some(target));
            let mut seen = vec![false; lattice.n_sites];
            for &vertex in &path {
                assert!(!seen[vertex], "path must be simple");
                seen[vertex] = true;
            }
            assert_eq!(
                path_weight(&edges, &weights).to_bits(),
                distance.to_bits(),
                "reconstructed path weight must equal the F1 distance bitwise"
            );
        }

        // All-0.5 weights: exactly the monotone corner paths are
        // geodesics, so the count is C(w + h - 2, w - 1) and ties are
        // guaranteed.
        let weights = vec![0.5_f64; lattice.n_edges()];
        let expected = binomial(width + height - 2, width - 1);
        let mut workspace = GeodesicCountWorkspace::new();
        let outcome = geodesic_count_run(&lattice, &weights, source, target, &mut workspace)
            .expect("all-half count");
        assert_eq!(
            outcome.geodesic_count, expected,
            "all-0.5 {width}x{height}: every monotone path is a geodesic"
        );
        assert!(expected > 1);
        eprintln!(
            "F3 all-0.5 {width}x{height}: geodesics = C({}, {}) = {expected}",
            width + height - 2,
            width - 1
        );
    }
    assert!(
        tied_cases > 0,
        "two-level fields must produce non-unique geodesics"
    );
    eprintln!(
        "F3 geodesic counting: {total_graphs} two-level graphs (3x3, 4x3) x 8 seeds; DAG == \
         brute force everywhere; {tied_cases} graphs with non-unique geodesics; total counted \
         geodesics {total_count}"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: directional time constants, convergence from above, shape
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn directional_time_constants_converge_from_above() {
    const DISTANCES: [usize; 3] = [8, 16, 32];
    const DIRECTIONS: [Direction; 4] = [
        Direction::East,
        Direction::North,
        Direction::Northeast,
        Direction::Southeast,
    ];
    const REALIZATIONS: usize = 150;
    // t_{0.975} on df = 149.
    const T975_DF149: f64 = 1.976;
    let geometry = ShapeGeometry::new(&DIRECTIONS, &DISTANCES, 16).expect("geometry");
    let level = |direction: usize, distance: usize| direction * DISTANCES.len() + distance;
    let mut sampler = ShapeSampler::new(geometry).expect("sampler");
    let mut estimator = ShapeEstimator::new(DIRECTIONS.len() * DISTANCES.len()).expect("estimator");
    let mut profiles = Vec::with_capacity(REALIZATIONS);
    for realization in 0..REALIZATIONS {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F3_1000 + realization as u64);
        let profile = sampler
            .realize(EdgeWeightLaw::exponential(1.0), &mut rng)
            .expect("realization")
            .to_vec();
        estimator.add(&profile).expect("add");
        profiles.push(profile);
    }
    let estimate = estimator.estimate().expect("estimate");

    let mut mu = vec![0.0_f64; DIRECTIONS.len() * DISTANCES.len()];
    let mut mu_sem = vec![0.0_f64; DIRECTIONS.len() * DISTANCES.len()];
    for (direction, &name) in DIRECTIONS.iter().enumerate() {
        let mut sequence = String::new();
        for (distance_index, &distance) in DISTANCES.iter().enumerate() {
            let index = level(direction, distance_index);
            let value = estimate.means[index] / distance as f64;
            let sem = estimate.standard_errors[index] / distance as f64;
            mu[index] = value;
            mu_sem[index] = sem;
            sequence.push_str(&format!("mu({distance})={value:.4}±{sem:.4} "));
            assert!(
                value - T975_DF149 * sem > 0.0,
                "{:?} n = {distance}: mu must be positive (CI lower bound)",
                name
            );
            assert!(
                value + T975_DF149 * sem < 1.0,
                "{:?} n = {distance}: subadditivity gives mu <= E tau_e = 1, strictly here \
                 (Theorem 2.4 for a non-trivial F)",
                name
            );
        }
        for distance_index in 1..DISTANCES.len() {
            let previous = level(direction, distance_index - 1);
            let current = level(direction, distance_index);
            let difference = mu[current] - mu[previous];
            let noise = (mu_sem[current].powi(2) + mu_sem[previous].powi(2)).sqrt();
            assert!(
                difference < 3.0 * noise,
                "{:?}: T(nu)/n increased by {difference:.5} (z = {:.2}) from n = {} to {} — \
                 convergence from above must hold within noise",
                name,
                difference / noise,
                DISTANCES[distance_index - 1],
                DISTANCES[distance_index]
            );
        }
        let last = level(direction, DISTANCES.len() - 1);
        let previous = level(direction, DISTANCES.len() - 2);
        let first = level(direction, 0);
        let recent_decrement = mu[previous] - mu[last];
        let earlier_decrement = mu[first] - mu[previous];
        let stabilization_noise = (mu_sem[last].powi(2) + mu_sem[previous].powi(2)).sqrt();
        eprintln!(
            "F3 direction {:?} (R = {REALIZATIONS}, exponential rate 1): {sequence}| \
             decrements 8->16 {earlier_decrement:.5}, 16->32 {recent_decrement:.5} \
             (noise {stabilization_noise:.5})",
            name
        );
        // Stabilization as flattening: the recent decrement must not grow
        // beyond the earlier one (within noise) and must stay inside a
        // finite finite-ladder band. The ladder is short, so this is
        // convergence evidence, not a limit claim.
        assert!(
            recent_decrement <= earlier_decrement + 3.0 * stabilization_noise,
            "{:?}: the 16->32 decrement {recent_decrement:.5} must not exceed the 8->16 \
             decrement {earlier_decrement:.5} beyond noise (flattening)",
            name
        );
        assert!(
            recent_decrement <= 0.06,
            "{:?}: the finite-ladder 16->32 decrement {recent_decrement:.5} must stay inside \
             the exploratory band 0.06",
            name
        );
    }

    // Rotational symmetry of the two axis directions (the square box and
    // centered origin are invariant under x <-> y; i.i.d. weights make the
    // distributions equal).
    let east = level(0, DISTANCES.len() - 1);
    let north = level(1, DISTANCES.len() - 1);
    let symmetry_noise = (mu_sem[east].powi(2) + mu_sem[north].powi(2)).sqrt();
    assert!(
        (mu[east] - mu[north]).abs() <= 3.0 * symmetry_noise,
        "axis symmetry: mu(east) = {:.4} and mu(north) = {:.4} must agree within noise",
        mu[east],
        mu[north]
    );

    // Shape consistency, paired per realization at the largest distance:
    // mu(east) + mu(north) >= mu(northeast) and >= mu(southeast).
    for (label, diagonal) in [("northeast", 2_usize), ("southeast", 3_usize)] {
        let values: Vec<f64> = profiles
            .iter()
            .map(|profile| {
                let n = DISTANCES[DISTANCES.len() - 1] as f64;
                profile[level(0, DISTANCES.len() - 1)] / n
                    + profile[level(1, DISTANCES.len() - 1)] / n
                    - profile[level(diagonal, DISTANCES.len() - 1)] / n
            })
            .collect();
        let (difference, sem) = mean_sem(&values);
        eprintln!(
            "F3 shape consistency mu(east)+mu(north)-mu({label}) = {difference:.5}±{sem:.5} \
             (R = {REALIZATIONS}, paired)"
        );
        assert!(
            difference - 3.0 * sem >= 0.0,
            "shape consistency: mu(east)+mu(north) must dominate mu({label}) within CI \
             (difference {difference:.5} ± {sem:.5})"
        );
    }
    eprintln!(
        "F3 time constants (exponential rate 1, box 97x97, margin 16, R = {REALIZATIONS}): \
         mu(east)=({:.4}, {:.4}, {:.4}) mu(north)=({:.4}, {:.4}, {:.4}) \
         mu(northeast)=({:.4}, {:.4}, {:.4}) mu(southeast)=({:.4}, {:.4}, {:.4})",
        mu[level(0, 0)],
        mu[level(0, 1)],
        mu[level(0, 2)],
        mu[level(1, 0)],
        mu[level(1, 1)],
        mu[level(1, 2)],
        mu[level(2, 0)],
        mu[level(2, 1)],
        mu[level(2, 2)],
        mu[level(3, 0)],
        mu[level(3, 1)],
        mu[level(3, 2)]
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: exploratory finite-horizon geodesic aging
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn geodesic_aging_is_exact_and_exploratory() {
    // Hand cases pin the statistic definition: shared leading prefix
    // edges over the short path's edge count.
    let aging = geodesic_aging(&[0, 1, 2], &[0, 1, 2, 3]).expect("prefix");
    assert_eq!(aging.shared_prefix_edges, 2);
    assert_eq!(aging.persistence, 1.0);
    let aging = geodesic_aging(&[0, 1, 2], &[0, 1, 3]).expect("shared leading edge");
    assert_eq!(aging.shared_prefix_edges, 1);
    assert_eq!(aging.persistence, 0.5);
    let aging = geodesic_aging(&[0, 1, 2], &[0, 3, 4]).expect("divergent");
    assert_eq!(aging.persistence, 0.0);
    assert!(geodesic_aging(&[0], &[0, 1]).is_none());
    assert!(geodesic_aging(&[0, 1], &[5, 1]).is_none());

    const RADII: [usize; 3] = [8, 16, 24];
    const REALIZATIONS: usize = 60;
    let mut means = Vec::with_capacity(RADII.len());
    let mut sems = Vec::with_capacity(RADII.len());
    for &radius in &RADII {
        let geometry = ShapeGeometry::new(&[Direction::East], &[radius, 2 * radius], radius)
            .expect("aging geometry");
        let mut sampler = ShapeSampler::new(geometry.clone()).expect("aging sampler");
        let mut persistence = Vec::with_capacity(REALIZATIONS);
        let mut path_checks = 0_usize;
        for realization in 0..REALIZATIONS {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                0x00F3_2000 + radius as u64 * 0x100 + realization as u64,
            );
            sampler
                .realize_tracked(EdgeWeightLaw::exponential(1.0), &mut rng)
                .expect("tracked realization");
            let workspace = sampler.workspace();
            let weights = sampler.weights();
            let mut short = Vec::new();
            let mut short_edges = Vec::new();
            let mut long = Vec::new();
            let mut long_edges = Vec::new();
            geodesic_path(
                workspace.predecessor_vertices(),
                workspace.predecessor_edges(),
                geometry.origin(),
                geometry.targets()[0].vertex,
                &mut short,
                &mut short_edges,
            )
            .expect("short path");
            geodesic_path(
                workspace.predecessor_vertices(),
                workspace.predecessor_edges(),
                geometry.origin(),
                geometry.targets()[1].vertex,
                &mut long,
                &mut long_edges,
            )
            .expect("long path");
            for (path, edges, target) in [
                (&short, &short_edges, geometry.targets()[0].vertex),
                (&long, &long_edges, geometry.targets()[1].vertex),
            ] {
                assert_eq!(path.first().copied(), Some(geometry.origin()));
                assert_eq!(path.last().copied(), Some(target));
                let mut seen = vec![false; geometry.lattice().n_sites];
                for &vertex in path {
                    assert!(!seen[vertex], "aging path must be simple");
                    seen[vertex] = true;
                }
                assert_eq!(
                    path_weight(edges, weights).to_bits(),
                    workspace.distances()[target].to_bits(),
                    "aging path must be a true geodesic (weight sum == F1 distance)"
                );
                path_checks += 1;
            }
            let aging = geodesic_aging(&short, &long).expect("same origin");
            assert!(
                (0.0..=1.0).contains(&aging.persistence),
                "persistence must lie in [0, 1]"
            );
            assert_eq!(aging.short_edges, short.len() - 1);
            assert_eq!(aging.long_edges, long.len() - 1);
            persistence.push(aging.persistence);
        }
        let (mean, sem) = mean_sem(&persistence);
        means.push(mean);
        sems.push(sem);
        eprintln!(
            "F3 aging east r={radius} R={REALIZATIONS}: persistence {mean:.4}±{sem:.4} \
             (EXPLORATORY; {path_checks} exact geodesic path checks)"
        );
    }
    // Exploratory convergence band only: the two largest radii must not
    // drift apart by more than 0.10 beyond 3 combined standard errors.
    // No asymptotic limit is claimed (the survey has no crisp
    // finite-horizon aging law; Question 30 is the closest anchor).
    let last = RADII.len() - 1;
    let drift = (means[last] - means[last - 1]).abs();
    let noise = (sems[last].powi(2) + sems[last - 1].powi(2)).sqrt();
    eprintln!(
        "F3 aging exploratory drift r={} vs r={}: {drift:.4} (noise {noise:.4}); sequence \
         {} — labelled exploratory, no asymptotic claim",
        RADII[last],
        RADII[last - 1],
        means
            .iter()
            .map(|mean| format!("{mean:.4}"))
            .collect::<Vec<_>>()
            .join(", ")
    );
    assert!(
        drift <= 3.0 * noise + 0.10,
        "aging: the exploratory drift {drift:.4} between the two largest radii exceeds the \
         generous band 0.10 + 3 noise {:.4}",
        3.0 * noise
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Draw contract and determinism
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn draw_contract_and_same_seed_determinism() {
    let geometry =
        ShapeGeometry::new(&[Direction::East, Direction::North], &[4, 8], 2).expect("geometry");
    let law = EdgeWeightLaw::exponential(1.0);
    let seed = 0x00F3_8000_u64;
    let mut first = ShapeSampler::new(geometry.clone()).expect("first");
    let mut second = ShapeSampler::new(geometry.clone()).expect("second");
    let mut first_rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut second_rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut reference_profile = Vec::new();
    for realization in 0..5 {
        let left = first
            .realize(law, &mut first_rng)
            .expect("first realization")
            .to_vec();
        let right = second
            .realize(law, &mut second_rng)
            .expect("second realization")
            .to_vec();
        for (index, (a, b)) in left.iter().zip(&right).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "realization {realization} entry {index}: same seed must be bit-identical"
            );
        }
        if realization == 0 {
            reference_profile = left;
        }
    }

    // The draw contract: one uniform per box edge in dense build_square
    // order, transformed by the law.
    let mut sampler = ShapeSampler::new(geometry.clone()).expect("replay sampler");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    sampler.realize(law, &mut rng).expect("replay run");
    let produced = sampler.weights().to_vec();
    let mut replay = Xoshiro256PlusPlus::seed_from_u64(seed);
    for (edge, &weight) in produced.iter().enumerate() {
        assert_eq!(
            weight.to_bits(),
            law.draw(replay.random::<f64>()).to_bits(),
            "edge {edge}: draw-contract replay"
        );
    }
    let mut tracked = ShapeSampler::new(geometry.clone()).expect("tracked sampler");
    let tracked_profile = tracked
        .realize_tracked(law, &mut Xoshiro256PlusPlus::seed_from_u64(seed))
        .expect("tracked run")
        .to_vec();
    let replayed_profile = tracked
        .realize_with_weights(&produced, GeodesicTracking::Off)
        .expect("weight replay")
        .to_vec();
    for (index, (a, b)) in tracked_profile.iter().zip(&replayed_profile).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "entry {index}: weight replay");
    }
    for (index, (a, b)) in reference_profile.iter().zip(&tracked_profile).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "entry {index}: tracking must not change the profile"
        );
    }
    eprintln!(
        "F3 draw contract: {} weight entries replayed bitwise; 5 realizations bit-identical \
         across samplers; tracked and weight-replayed profiles identical",
        produced.len()
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Typed rejections and state isolation
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn typed_rejections_leave_state_untouched() {
    // Geometry domains.
    assert!(matches!(
        ShapeGeometry::new(&[], &[1], 0),
        Err(ShapeError::NoDirections)
    ));
    assert!(matches!(
        ShapeGeometry::new(&[Direction::East], &[], 0),
        Err(ShapeError::NoDistances)
    ));
    assert!(matches!(
        ShapeGeometry::new(&[Direction::East], &[0, 1], 0),
        Err(ShapeError::ZeroDistance { index: 0 })
    ));
    assert!(matches!(
        ShapeGeometry::new(&[Direction::East], &[1], SHAPE_SIDE_CAP),
        Err(ShapeError::SideAboveCap { .. })
    ));
    assert!(matches!(
        ShapeGeometry::new(&[Direction::East], &[1], usize::MAX),
        Err(ShapeError::SideOverflow)
    ));

    // Sampler domains; every rejection leaves the last profile intact.
    let geometry = ShapeGeometry::new(&[Direction::East], &[4], 2).expect("geometry");
    let mut sampler = ShapeSampler::new(geometry.clone()).expect("sampler");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00F3_9000);
    let profile = sampler
        .realize(EdgeWeightLaw::exponential(1.0), &mut rng)
        .expect("valid")
        .to_vec();
    assert!(matches!(
        sampler.realize_with_weights(&[1.0; 2], GeodesicTracking::Off),
        Err(ShapeError::WeightLength { .. })
    ));
    assert_eq!(
        sampler.profile(),
        &profile[..],
        "length rejection keeps profile"
    );
    for rate in [0.0_f64, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            sampler.realize(EdgeWeightLaw::exponential(rate), &mut rng),
            Err(ShapeError::Law(_))
        ));
    }
    assert_eq!(
        sampler.profile(),
        &profile[..],
        "law rejection keeps profile"
    );
    let mut bad = sampler.weights().to_vec();
    bad[2] = -0.25;
    assert!(matches!(
        sampler.realize_with_weights(&bad, GeodesicTracking::Off),
        Err(ShapeError::Passage(FirstPassageError::NegativeWeight {
            edge: 2,
            ..
        }))
    ));
    bad[2] = f64::NAN;
    assert!(matches!(
        sampler.realize_with_weights(&bad, GeodesicTracking::Off),
        Err(ShapeError::Passage(FirstPassageError::WeightNotFinite {
            edge: 2,
            ..
        }))
    ));
    assert_eq!(
        sampler.profile(),
        &profile[..],
        "weight rejections keep profile"
    );

    // Estimator domains.
    assert!(matches!(ShapeEstimator::new(0), Err(ShapeError::NoTargets)));
    let mut estimator = ShapeEstimator::new(1).expect("estimator");
    estimator.add(&profile).expect("matching add");
    assert_eq!(estimator.samples(), 1);
    assert!(estimator.estimate().is_none());
    assert!(matches!(
        estimator.add(&[0.0; 2]),
        Err(ShapeError::ProfileLength {
            expected: 1,
            actual: 2
        })
    ));
    assert_eq!(estimator.samples(), 1, "rejection adds no sample");

    // Geodesic counting: zero weights are a typed rejection with the
    // previous counts and distances untouched.
    let lattice = build_square(3, 3, false);
    let weights = two_level_weights(lattice.n_edges(), 0x00F3_9100);
    let mut count_workspace = GeodesicCountWorkspace::new();
    let source = 0;
    let target = lattice.n_sites - 1;
    let valid = geodesic_count_run(&lattice, &weights, source, target, &mut count_workspace)
        .expect("valid count");
    let ways_before = count_workspace.ways().to_vec();
    let distances_before = count_workspace.distances().to_vec();
    let mut zero = weights.clone();
    zero[5] = 0.0;
    assert!(matches!(
        geodesic_count_run(&lattice, &zero, source, target, &mut count_workspace),
        Err(ShapeError::ZeroWeight { edge: 5 })
    ));
    assert_eq!(count_workspace.ways(), &ways_before[..]);
    assert_eq!(count_workspace.distances(), &distances_before[..]);
    assert!(matches!(
        geodesic_count_run(
            &lattice,
            &weights[..3],
            source,
            target,
            &mut count_workspace
        ),
        Err(ShapeError::WeightLength { .. })
    ));
    assert!(matches!(
        geodesic_count_run(
            &lattice,
            &weights,
            lattice.n_sites,
            target,
            &mut count_workspace
        ),
        Err(ShapeError::SourceDomain { .. })
    ));
    assert!(matches!(
        geodesic_count_run(
            &lattice,
            &weights,
            source,
            lattice.n_sites,
            &mut count_workspace
        ),
        Err(ShapeError::TargetDomain { .. })
    ));
    assert_eq!(count_workspace.ways(), &ways_before[..]);

    // Disconnected graph: target unreachable -> count 0 (not an error);
    // the path reconstruction refuses it.
    let disconnected = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
        ],
    );
    let pair_weights = vec![1.0_f64; 2];
    let outcome = geodesic_count_run(&disconnected, &pair_weights, 0, 3, &mut count_workspace)
        .expect("disconnected count");
    assert_eq!(outcome.geodesic_count, 0);
    assert_eq!(outcome.reached, 2);
    let mut passage = FirstPassageWorkspace::new();
    first_passage_run(
        &disconnected,
        &pair_weights,
        0,
        GeodesicTracking::Predecessors,
        &mut passage,
    )
    .expect("disconnected pass");
    let mut path = Vec::new();
    let mut path_edges = Vec::new();
    assert!(matches!(
        geodesic_path(
            passage.predecessor_vertices(),
            passage.predecessor_edges(),
            0,
            3,
            &mut path,
            &mut path_edges,
        ),
        Err(ShapeError::UnreachableTarget { target: 3 })
    ));
    assert!(matches!(
        geodesic_path(&[1], &[], 0, 0, &mut path, &mut path_edges),
        Err(ShapeError::PredecessorLength { .. })
    ));
    assert!(matches!(
        geodesic_path(
            passage.predecessor_vertices(),
            passage.predecessor_edges(),
            0,
            4,
            &mut path,
            &mut path_edges,
        ),
        Err(ShapeError::TargetDomain { .. })
    ));
    assert!(valid.geodesic_count >= 1);
    eprintln!(
        "F3 typed rejections: geometry (5), sampler (7), estimator (3), counting (6), path (3) \
         — state untouched on every counting rejection"
    );
}
