//! T4 continuum winding validation: the integer-cell displacement source
//! (P1 `PairGraph` + `ContinuumEmbedding`) feeding the unchanged W2
//! `analyze_winding`.
//!
//! Gates (roadmap v2 §2 T4):
//!
//! 1. **Exact negation**: for every physical edge, the embedding's
//!    displacement from `i` to `j` is the exact negation of `j` to `i` (the
//!    `PeriodicEmbedding` contract) and equals an independently recomputed
//!    `round((d - r) / L)` table, across sizes, seeds, and dimensions 2 and 3.
//! 2. **Hand-derived small cases**: integer tables for tiny configurations
//!    (a ring crossing once; a two-particle pair at wrap distance; an open
//!    chain containing a wrap edge) asserted exactly.
//! 3. **P2 independent reference (the physics gate)**: periodic Poisson disk
//!    configurations (Boolean model) at spanning and sub-spanning densities —
//!    for EVERY component, winding rank and wraps verdict from production
//!    must agree exactly with a test-local BFS unfolding.
//! 4. **Precondition boundary**: `r_max >= L/2` rejected with the axis named;
//!    a cutoff just under `L/2` works end to end.
//! 5. **Tolerance rejection**: hand-built pair lists whose implied cell
//!    offset is fractional or unrepresentable fail with the typed error.
//!
//! # What the physics gate compares (and why)
//!
//! The reference unfolds each component independently: BFS over the pair
//! graph accumulating float minimum-image displacements into unwrapped
//! coordinates, then reads the winding subgroup off the NON-tree edges — for
//! every cycle-closing edge `(u, v)`, the unwrapped delta minus the
//! minimum-image delta is `k * L` per axis, and the collected integer
//! vectors generate the winding subgroup. Rank and wraps from that
//! float/tree computation must equal the production integer-potential
//! verdicts component by component; this is exact (no statistical window).
//!
//! The unwrapped **extent** (max - min of an unwrapped axis) is logically
//! INDEPENDENT of winding, in both directions, and the census counts both
//! events instead of asserting either away. A tree-shaped tendril can
//! chain around the torus and reach extent `>= L_c` with no cycle at all
//! (spanning-without-wrapping, "crossing is not wrapping"); conversely a
//! compact cycle can wind with every representative inside an extent
//! strictly below `L_c` — the boundary-exact hand ring below is exactly
//! that (cycle sum `L`, unwrapped extent `8.5` in a box of `10`; review
//! 2026-10-23 demonstrated the same violation on 3/20 fresh small-box
//! Boolean configurations). Both are real physics, not analyzer
//! disagreements; at the committed far-from-threshold densities both
//! counters happen to read zero, but they are reported, never asserted.

use std::collections::VecDeque;

use cmc_rs::{
    analyze_winding, CandidatePairGraph, ContinuumEmbedding, ContinuumWindingError, EdgeActivity,
    EdgeId, GraphView, OrthorhombicCell, PairGraph, PeriodicEmbedding, UndirectedGraphView,
    VertexActivity, VertexId, WindingResult, WindingWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const SEEDS: u64 = 20;

#[derive(Clone, Debug)]
struct VertexFlags(Vec<bool>);

impl VertexActivity for VertexFlags {
    fn vertex_count(&self) -> usize {
        self.0.len()
    }
    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.0[vertex.index()]
    }
}

#[derive(Clone, Debug)]
struct EdgeFlags(Vec<bool>);

impl EdgeActivity for EdgeFlags {
    fn edge_count(&self) -> usize {
        self.0.len()
    }
    fn edge_active(&self, edge: EdgeId) -> bool {
        self.0[edge.index()]
    }
}

fn all_active(graph: &PairGraph) -> (VertexFlags, EdgeFlags) {
    (
        VertexFlags(vec![true; graph.vertex_count()]),
        EdgeFlags(vec![true; graph.edge_count()]),
    )
}

fn production<'a>(
    graph: &PairGraph,
    embedding: &ContinuumEmbedding,
    workspace: &'a mut WindingWorkspace,
) -> WindingResult<'a> {
    let (vertices, edges) = all_active(graph);
    analyze_winding(graph, embedding, &vertices, &edges, workspace)
        .unwrap_or_else(|error| panic!("production winding analysis failed: {error}"))
}

fn flatten<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
    positions.iter().flat_map(|p| p.iter().copied()).collect()
}

fn random_positions<const D: usize>(n: usize, lengths: [f64; D], seed: u64) -> Vec<[f64; D]> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut position = [0.0; D];
            for axis in 0..D {
                position[axis] = rng.random_range(0.0..lengths[axis]);
            }
            position
        })
        .collect()
}

/// Test-local minimum image (division form, unlike the production
/// reciprocal-multiply form) — the reference never calls the cell.
fn min_image(raw: f64, length: f64) -> f64 {
    raw - length * (raw / length).round()
}

// ── Gate 1: exact negation and independent recomputation ─────────────────

fn negation_sweep<const D: usize>(
    lengths: [f64; D],
    size: usize,
    seed: u64,
    cutoff_scale: f64,
    counters: &mut (usize, usize),
) {
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let shortest = lengths.iter().copied().fold(f64::INFINITY, f64::min);
    let cutoff = cutoff_scale * shortest / 2.0; // strictly under L/2
    let positions = random_positions(size, lengths, seed);
    let coordinates = flatten(&positions);
    let graph = CandidatePairGraph::new(&coordinates, &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();
    let embedding =
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, cutoff).unwrap();
    for (edge_index, &[left, right]) in graph.pairs().iter().enumerate() {
        let edge = graph.edge_id(edge_index).expect("dense edge index");
        let forward = embedding
            .edge_displacement(edge, graph.vertex_id(left).expect("dense index"))
            .expect("left endpoint");
        let backward = embedding
            .edge_displacement(edge, graph.vertex_id(right).expect("dense index"))
            .expect("right endpoint");
        for axis in 0..D {
            let forward_axis = forward.axis(axis).expect("axis in range");
            let backward_axis = backward.axis(axis).expect("axis in range");
            assert_eq!(
                forward_axis, -backward_axis,
                "size {size} seed {seed:#x} edge {edge_index} axis {axis}: \
                 i->j must negate j->i exactly"
            );
            // Independent recomputation of the stored table:
            // n = round((d - r)/L) with a division-form minimum image.
            let raw = positions[right][axis] - positions[left][axis];
            let expected = ((min_image(raw, lengths[axis]) - raw) / lengths[axis]).round();
            assert_eq!(
                forward_axis, expected as i32,
                "size {size} seed {seed:#x} edge {edge_index} axis {axis}: \
                 stored vector must match the recomputed cell offset"
            );
            counters.1 += 1;
        }
        // The (edge, endpoint) key rejects non-endpoints.
        if graph.vertex_count() > 2 {
            let outsider = (left + 1) % graph.vertex_count();
            if outsider != left && outsider != right {
                assert!(
                    embedding
                        .edge_displacement(edge, graph.vertex_id(outsider).expect("dense index"))
                        .is_none(),
                    "vertex {outsider} is not an endpoint of edge {edge_index}"
                );
            }
        }
        counters.0 += 1;
    }
}

#[test]
fn edge_displacements_are_exact_negations_across_sizes_and_dimensions() {
    let mut counters = (0usize, 0usize);
    for size in [40usize, 160, 500] {
        for seed in [0x51u64, 0x52, 0x53] {
            for cutoff_scale in [0.3, 0.95] {
                negation_sweep([12.0, 12.0], size, seed, cutoff_scale, &mut counters);
                negation_sweep([12.0, 12.0, 12.0], size, seed, cutoff_scale, &mut counters);
            }
        }
    }
    eprintln!(
        "CONTINUUM_WINDING_NEGATION dimensions=2,3 sizes=40,160,500 seeds=3 cutoff_scales=0.3,0.95 \
         checked_edges={} checked_axis_entries={} disagreements=0",
        counters.0, counters.1
    );
}

// ── Gate 2: hand-derived small configurations ────────────────────────────

#[test]
fn hand_ring_crosses_exactly_once_and_open_wrap_chains_do_not_wind() {
    // 2D ring of five particles along x (y pinned): pairs (0,1),(1,2),(2,3),
    // (3,4) in-box at distance exactly 2 and (0,4) through the wrap at
    // exactly the cutoff (raw 8.0, minimum image -2.0). Hand table, stored
    // i -> j: in-box edges 0, the (0,4) edge (-2.0 - 8.0)/10 = -1.
    let cell = OrthorhombicCell::new([10.0, 10.0]).unwrap();
    let ring = [[1.0, 5.0], [3.0, 5.0], [5.0, 5.0], [7.0, 5.0], [9.0, 5.0]];
    let coordinates = flatten(&ring);
    let graph = CandidatePairGraph::new(&coordinates, &cell, 2.0)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(
        graph.pairs(),
        &[[0usize, 1usize], [0, 4], [1, 2], [2, 3], [3, 4]],
        "hand ring pair set in P1's cell-major emission order"
    );
    let embedding =
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, 2.0).unwrap();
    for (edge_index, &[left, right]) in graph.pairs().iter().enumerate() {
        let edge = graph.edge_id(edge_index).expect("dense edge index");
        let vector = embedding
            .edge_displacement(edge, graph.vertex_id(left).expect("dense index"))
            .unwrap();
        let expected_x = -i32::from(left == 0 && right == 4);
        assert_eq!(
            (vector.axis(0), vector.axis(1)),
            (Some(expected_x), Some(0)),
            "edge {edge_index} ({left}, {right})"
        );
    }

    let mut workspace = WindingWorkspace::new();
    let result = production(&graph, &embedding, &mut workspace);
    assert_eq!(result.dimensions(), 2);
    assert_eq!(result.summary().component_count, 1);
    assert_eq!(result.summary().wrapping_component_count, 1);
    let info = result.components()[0];
    assert_eq!(
        (info.component().index(), info.size(), info.rank()),
        (0, 5, 1)
    );
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 0][..]],
        "single x wrap: generator is the +x unit vector"
    );

    // Independence pin (review follow-up 2026-10-23): this ring WINDS
    // (cycle sum 10 = L) while every unwrapped representative lies in
    // [1, 9] — extent 8 < L. Winding does not imply extent-spanning; the
    // converse is pinned by the wrap pair below (extent-spanning edge,
    // rank 0). Neither direction is assertable in the census.
    let unwrapped_extent = 9.0 - 1.0;
    assert!(
        info.rank() > 0 && unwrapped_extent < 10.0,
        "the hand ring must stay the compact-wrapping counterexample"
    );

    // Two particles at wrap distance in an 8-box (image -1.0): connected
    // through a crossing edge, but a tree has no cycle and cannot wind.
    let pair_cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let pair = [[0.5, 4.0], [7.5, 4.0]];
    let pair_coordinates = flatten(&pair);
    let pair_graph = CandidatePairGraph::new(&pair_coordinates, &pair_cell, 2.0)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(pair_graph.pairs(), &[[0usize, 1usize]]);
    let pair_embedding =
        ContinuumEmbedding::from_pairs(pair_graph.pairs(), &pair_coordinates, &pair_cell, 2.0)
            .unwrap();
    let vector = pair_embedding
        .edge_displacement(
            pair_graph.edge_id(0).expect("dense edge index"),
            pair_graph.vertex_id(0).expect("dense index"),
        )
        .unwrap();
    assert_eq!((vector.axis(0), vector.axis(1)), (Some(-1), Some(0)));
    let result = production(&pair_graph, &pair_embedding, &mut workspace);
    assert_eq!(result.summary().component_count, 1);
    assert_eq!(result.summary().wrapping_component_count, 0);
    assert_eq!(result.components()[0].rank(), 0);

    // Three-particle open chain whose second edge crosses the boundary:
    // x = 8.4 - 9.7 - 0.3 (wrap edge raw -9.4, image 0.6, n = +1), no
    // closing edge. Still rank zero.
    let chain = [[8.4, 5.0], [9.7, 5.0], [0.3, 5.0]];
    let chain_coordinates = flatten(&chain);
    let chain_graph = CandidatePairGraph::new(&chain_coordinates, &cell, 1.4)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(
        chain_graph.pairs(),
        &[[0usize, 1usize], [1, 2]],
        "only the consecutive chain pairs are within 1.4"
    );
    let chain_embedding =
        ContinuumEmbedding::from_pairs(chain_graph.pairs(), &chain_coordinates, &cell, 1.4)
            .unwrap();
    let vector = chain_embedding
        .edge_displacement(
            chain_graph.edge_id(1).expect("dense edge index"),
            chain_graph.vertex_id(1).expect("dense index"),
        )
        .unwrap();
    assert_eq!((vector.axis(0), vector.axis(1)), (Some(1), Some(0)));
    let result = production(&chain_graph, &chain_embedding, &mut workspace);
    assert_eq!(result.summary().wrapping_component_count, 0);
    assert_eq!(result.components()[0].rank(), 0);
}

// ── Gate 3: P2 independent reference (physics gate) ──────────────────────

/// One component's reference facts.
struct ReferenceComponent {
    /// Canonical identity (lowest vertex index).
    seed: usize,
    members: usize,
    /// Winding rank from the float BFS + non-tree cycle sums.
    rank: usize,
    /// Per axis: unwrapped extent `>= L - tol`.
    spans: Vec<bool>,
}

#[derive(Default)]
struct Census {
    configurations: usize,
    particles_total: usize,
    components: usize,
    wrapping_components: usize,
    spanning_by_extent: usize,
    /// Components with extent >= L on some axis but winding rank 0 (the
    /// tree-tendril event: spanning is not wrapping).
    spanning_without_wrapping: usize,
    /// Components winding around some axis whose unwrapped representatives
    /// all sit within an extent below L (the compact-cycle event: winding
    /// does not imply extent-spanning either).
    wrapping_without_spanning: usize,
}

/// Independent reference: BFS-unfold every component of a pair list,
/// accumulate unwrapped coordinates through float minimum-image steps, and
/// generate the winding subgroup from the non-tree (cycle-closing) edges.
/// Extent-based spanning is recorded per axis alongside.
fn unfold_reference<const D: usize>(
    pairs: &[[usize; 2]],
    positions: &[[f64; D]],
    lengths: &[f64; D],
) -> Vec<ReferenceComponent> {
    let n = positions.len();
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (edge_index, &[left, right]) in pairs.iter().enumerate() {
        adjacency[left].push((right, edge_index));
        adjacency[right].push((left, edge_index));
    }
    let mut unwrapped = vec![[f64::NAN; D]; n];
    let mut component_of = vec![usize::MAX; n];
    let mut tree_edge = vec![false; pairs.len()];
    let mut components = Vec::new();
    for seed in 0..n {
        if !unwrapped[seed][0].is_nan() {
            continue;
        }
        let component_index = components.len();
        unwrapped[seed] = positions[seed];
        component_of[seed] = component_index;
        let mut queue = VecDeque::new();
        queue.push_back(seed);
        let mut members = 1usize;
        let mut generators: Vec<[i64; D]> = Vec::new();
        let mut low = positions[seed];
        let mut high = positions[seed];
        while let Some(vertex) = queue.pop_front() {
            for &(neighbor, edge_index) in &adjacency[vertex] {
                if !unwrapped[neighbor][0].is_nan() {
                    continue;
                }
                unwrapped[neighbor] = unwrapped[vertex];
                for axis in 0..D {
                    let raw = positions[neighbor][axis] - positions[vertex][axis];
                    unwrapped[neighbor][axis] += min_image(raw, lengths[axis]);
                    low[axis] = low[axis].min(unwrapped[neighbor][axis]);
                    high[axis] = high[axis].max(unwrapped[neighbor][axis]);
                }
                component_of[neighbor] = component_index;
                tree_edge[edge_index] = true;
                members += 1;
                queue.push_back(neighbor);
            }
        }
        // Cycle-closing edges OF THIS COMPONENT: the unwrapped delta minus
        // the minimum-image delta is an integer number of box lengths per
        // axis — one generator observation each.
        for (edge_index, &[left, right]) in pairs.iter().enumerate() {
            if tree_edge[edge_index] || component_of[left] != component_index {
                continue;
            }
            let mut generator = [0i64; D];
            for axis in 0..D {
                let raw = positions[right][axis] - positions[left][axis];
                let quotient = (unwrapped[right][axis]
                    - unwrapped[left][axis]
                    - min_image(raw, lengths[axis]))
                    / lengths[axis];
                let rounded = quotient.round();
                assert!(
                    (quotient - rounded).abs() < 1e-6,
                    "reference unwrapping is not consistent on edge {edge_index} axis {axis}: \
                     left={left} right={right} quotient={quotient}"
                );
                generator[axis] = rounded as i64;
            }
            if generator.iter().any(|&component| component != 0) {
                generators.push(generator);
            }
        }
        components.push(ReferenceComponent {
            seed,
            members,
            rank: rational_rank(&generators),
            spans: (0..D)
                .map(|axis| high[axis] - low[axis] >= lengths[axis] - 1e-6 * lengths[axis])
                .collect(),
        });
    }
    components
}

/// Rank over the rationals of an integer vector set (equals the rank of the
/// subgroup they generate).
fn rational_rank<const D: usize>(vectors: &[[i64; D]]) -> usize {
    let mut rows: Vec<Vec<f64>> = vectors
        .iter()
        .map(|vector| vector.iter().map(|&component| component as f64).collect())
        .collect();
    let mut rank = 0;
    for column in 0..D {
        let pivot = (rank..rows.len())
            .filter(|&row| rows[row][column].abs() > 1e-9)
            .max_by(|&a, &b| {
                rows[a][column]
                    .abs()
                    .partial_cmp(&rows[b][column].abs())
                    .unwrap_or(core::cmp::Ordering::Equal)
            });
        let Some(pivot) = pivot else {
            continue;
        };
        rows.swap(rank, pivot);
        for row in (rank + 1)..rows.len() {
            let factor = rows[row][column] / rows[rank][column];
            let pivot_row = rows[rank].clone();
            for (value, pivot_value) in rows[row].iter_mut().zip(&pivot_row) {
                *value -= factor * pivot_value;
            }
        }
        rank += 1;
    }
    rank
}

/// Draw one Poisson count (Knuth multiplication method in 8 chunks so the
/// product threshold stays far above underflow; independent of the
/// production `PointProcess`).
fn poisson_count(rng: &mut Xoshiro256PlusPlus, mean: f64) -> usize {
    const CHUNKS: usize = 8;
    let part = mean / CHUNKS as f64;
    let threshold = (-part).exp();
    (0..CHUNKS)
        .map(|_| {
            let mut product = 1.0;
            let mut draws = 0;
            loop {
                draws += 1;
                product *= rng.random::<f64>();
                if product <= threshold {
                    return draws - 1;
                }
            }
        })
        .sum()
}

/// One Poisson Boolean-model configuration: draw the count, uniform wrapped
/// positions, P1 pairs at the contact cutoff `2r` (monodisperse overlap
/// criterion), embedding, production analysis, reference comparison.
fn run_boolean_configuration<const D: usize>(
    lengths: [f64; D],
    mean_particles: usize,
    filling_fraction: f64,
    seed: u64,
    census: &mut Census,
) -> usize {
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let unit_ball = match D {
        2 => std::f64::consts::PI,
        3 => 4.0 / 3.0 * std::f64::consts::PI,
        _ => unreachable!("gate covers 2D and 3D"),
    };
    let volume: f64 = lengths.iter().product();
    let radius =
        (filling_fraction * volume / (mean_particles as f64 * unit_ball)).powf(1.0 / D as f64);
    let cutoff = 2.0 * radius;
    assert!(
        cutoff < lengths.iter().copied().fold(f64::INFINITY, f64::min) / 2.0,
        "test geometry must satisfy the winding precondition"
    );

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let count = poisson_count(&mut rng, mean_particles as f64);
    let mut coordinates = Vec::with_capacity(count * D);
    for _ in 0..count {
        for &length in &lengths {
            coordinates.push(rng.random_range(0.0..length));
        }
    }
    let graph = CandidatePairGraph::new(&coordinates, &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();
    let embedding =
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, cutoff).unwrap();

    let positions: Vec<[f64; D]> = (0..count)
        .map(|particle| {
            let mut position = [0.0; D];
            position.copy_from_slice(&coordinates[particle * D..][..D]);
            position
        })
        .collect();
    let reference = unfold_reference(graph.pairs(), &positions, &lengths);

    let mut workspace = WindingWorkspace::new();
    let result = production(&graph, &embedding, &mut workspace);

    assert_eq!(
        result.components().len(),
        reference.len(),
        "component count, seed {seed:#x}"
    );
    census.configurations += 1;
    census.particles_total += count;
    census.components += reference.len();
    let mut wrapping_this_configuration = 0;
    for (info, reference_component) in result.components().iter().zip(reference.iter()) {
        assert_eq!(
            info.component().index(),
            reference_component.seed,
            "canonical identity, seed {seed:#x}"
        );
        assert_eq!(
            info.size(),
            reference_component.members,
            "size, seed {seed:#x}"
        );
        assert_eq!(
            info.rank(),
            reference_component.rank,
            "seed {seed:#x} component {}: winding rank (production integer potentials vs \
             reference float BFS cycle sums)",
            reference_component.seed
        );
        let index = result
            .component_index(info.component())
            .expect("component listed");
        assert_eq!(
            result.wraps(index),
            reference_component.rank > 0,
            "wraps verdict, seed {seed:#x} component {}",
            reference_component.seed
        );
        if info.rank() > 0 {
            census.wrapping_components += 1;
            wrapping_this_configuration += 1;
        }
        // Winding and per-representative extent are logically independent
        // (see the module doc): a compact cycle can wind with extent below
        // L, so the event is counted, never asserted away.
        let mut wraps_axis_without_extent = false;
        for generator in result.generators(index) {
            for (axis, &component) in generator.iter().enumerate() {
                if component != 0 && !reference_component.spans[axis] {
                    wraps_axis_without_extent = true;
                }
            }
        }
        if wraps_axis_without_extent {
            census.wrapping_without_spanning += 1;
        }
        if reference_component.spans.iter().any(|&spans| spans) {
            census.spanning_by_extent += 1;
            if info.rank() == 0 {
                census.spanning_without_wrapping += 1;
            }
        }
    }
    wrapping_this_configuration
}

fn boolean_census<const D: usize>(
    dimension: usize,
    lengths: [f64; D],
    mean_particles: usize,
    spanning_fraction: f64,
    sub_fraction: f64,
) {
    let mut spanning_census = Census::default();
    let mut sub_census = Census::default();
    let mut spanning_configurations_with_wrap = 0;
    for seed in 0..SEEDS {
        let wrapping = run_boolean_configuration(
            lengths,
            mean_particles,
            spanning_fraction,
            0x54_0000 + seed,
            &mut spanning_census,
        );
        spanning_configurations_with_wrap += usize::from(wrapping > 0);
        run_boolean_configuration(
            lengths,
            mean_particles,
            sub_fraction,
            0x55_0000 + seed,
            &mut sub_census,
        );
    }
    for (regime, eta, census) in [
        ("spanning", spanning_fraction, &spanning_census),
        ("sub", sub_fraction, &sub_census),
    ] {
        eprintln!(
            "CONTINUUM_WINDING_CENSUS regime={regime} dim={dimension} eta={eta} seeds={SEEDS} \
             mean_particles={mean_particles} configurations={count} particles={particles} \
             components={components} wrapping={wrapping} spanning_by_extent={spanning} \
             spanning_without_wrapping={snakes} wrapping_without_spanning={compacts} \
             disagreements=0",
            count = census.configurations,
            particles = census.particles_total,
            components = census.components,
            wrapping = census.wrapping_components,
            spanning = census.spanning_by_extent,
            snakes = census.spanning_without_wrapping,
            compacts = census.wrapping_without_spanning,
        );
    }
    // Density sanity: at the supercritical filling fraction every sampled
    // configuration contains a wrapping component (the gate must not be
    // vacuous); far below threshold no component wraps.
    assert_eq!(
        spanning_configurations_with_wrap, SEEDS as usize,
        "supercritical ensemble must wrap in every configuration"
    );
    assert_eq!(
        sub_census.wrapping_components, 0,
        "sub-critical ensemble wrapped (filling fraction not sub-spanning?)"
    );
}

#[test]
fn boolean_model_2d_matches_bfs_unfolding_exactly() {
    // eta_c(2D disks) ~ 1.128 (Mertens & Moore 2012, the P2 anchor): 1.5 is
    // safely spanning, 0.7 safely sub-spanning at this size.
    boolean_census(2, [40.0, 40.0], 2500, 1.5, 0.7);
}

#[test]
fn boolean_model_3d_matches_bfs_unfolding_exactly() {
    // eta_c(3D spheres) ~ 0.34: 0.6 spanning, 0.15 sub-spanning.
    boolean_census(3, [15.0, 15.0, 15.0], 4000, 0.6, 0.15);
}

// ── Gate 4: precondition boundary end to end ─────────────────────────────

#[test]
fn half_box_cutoff_is_rejected_with_the_axis_and_just_under_works() {
    let cell = OrthorhombicCell::new([4.0, 10.0]).unwrap();
    let positions = [[0.1, 5.0], [2.0, 5.0], [3.95, 5.0]];
    let coordinates = flatten(&positions);
    // P1's cell list ACCEPTS a cutoff of exactly half the shortest side
    // (minimum-image tie); the winding embedding must reject it, axis named.
    let graph = CandidatePairGraph::new(&coordinates, &cell, 2.0)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, 2.0).unwrap_err(),
        ContinuumWindingError::CutoffSpansHalfBox {
            axis: 0,
            cutoff: 2.0,
            half_length: 2.0
        }
    );
    // Reaching half the long side (5.0) is already past half the short side,
    // so axis 0 is still the named violator.
    assert_eq!(
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, 5.0).unwrap_err(),
        ContinuumWindingError::CutoffSpansHalfBox {
            axis: 0,
            cutoff: 5.0,
            half_length: 2.0
        }
    );

    // Just under L_x/2 the full path works: the three particles form a ring
    // (0-1 at 1.9, 1-2 at 1.95, 0-2 through the wrap at 0.15) that winds x.
    let boundary = 2.0 - 1e-9;
    let graph = CandidatePairGraph::new(&coordinates, &cell, boundary)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(
        graph.pairs(),
        &[[0usize, 1usize], [0, 2], [1, 2]],
        "hand pair set at the boundary cutoff (cell-major order)"
    );
    let embedding =
        ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, boundary).unwrap();
    let mut workspace = WindingWorkspace::new();
    let result = production(&graph, &embedding, &mut workspace);
    assert_eq!(result.summary().wrapping_component_count, 1);
    assert_eq!(result.components()[0].rank(), 1);
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 0][..]],
        "the boundary ring winds x exactly once"
    );
}

// ── Gate 5: tolerance rejections through the public constructor ─────────

#[test]
fn inconsistent_representatives_fail_loudly_never_round_silently() {
    // Hand-built pair list (the constructor accepts any pair table, so the
    // crafting bypasses the normal P1 path). Representatives 1e8 box lengths
    // apart: the minimum-image arithmetic loses the integer to rounding
    // (residual ~1.5e-8 > 1e-9) — typed error, no silent rounding.
    let cell = OrthorhombicCell::new([0.1]).unwrap();
    let error =
        ContinuumEmbedding::from_pairs(&[[0usize, 1usize]], &[0.0, 10_000_000.05], &cell, 0.04)
            .unwrap_err();
    assert!(
        matches!(
            error,
            ContinuumWindingError::NonIntegerCellOffset {
                edge: 0,
                axis: 0,
                ..
            }
        ),
        "{error}"
    );

    // An integer-valued but unrepresentable offset (1e300 boxes): typed
    // error, never saturation to i32::MAX.
    let cell = OrthorhombicCell::new([1.0]).unwrap();
    let error =
        ContinuumEmbedding::from_pairs(&[[0usize, 1usize]], &[0.0, 1e300], &cell, 0.4).unwrap_err();
    assert!(
        matches!(
            error,
            ContinuumWindingError::NonIntegerCellOffset {
                edge: 0,
                axis: 0,
                ..
            }
        ),
        "{error}"
    );
}
