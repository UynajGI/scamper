//! P1 production validation: continuum candidate-pair adapter.
//!
//! Gates (roadmap §14 P1/P2 candidate-pair prerequisite), each against an
//! independently coded reference:
//!
//! 1. **Cell-list vs all-pairs parity**: candidate pairs on random 2D and
//!    3D periodic configurations match an exhaustive `O(N^2)` scan exactly
//!    (pair SET equality; the search is independent, while the
//!    minimum-image metric formula is shared — the gate targets the
//!    stencil search, not floating-point noise).
//! 2. **Inclusive boundary**: pairs at exactly the cutoff (direct and
//!    through a periodic image) are kept (`<=`, documented in the module);
//!    just-short cutoffs drop them.
//! 3. **Periodic-image pairs**: pairs close only through the wrapped image
//!    are found.
//! 4. **Degenerate sizes**: a single particle yields an empty pair set and
//!    a one-vertex graph; zero particles is a typed rejection.
//! 5. **Determinism**: identical inputs reproduce the identical pair
//!    sequence, across both entry points (flat adapter and configuration).
//! 6. **Translation invariance**: an exact lattice-vector shift (bitwise
//!    stable under wrap) reproduces the identical sequence; a generic
//!    shift preserves the pair set.
//! 7. **F3 bridge**: `analyze`/`analyze_with_labels` on the ephemeral
//!    `PairGraph` agree with an independent adjacency-BFS component
//!    analysis of the same pair set, summary field by field and label by
//!    label.
//! 8. **Prepared cell list**: reuse of a caller-built `CellList` matches
//!    the fresh build and rejects stale lists.
//! 9. **Throughput record**: N = 10_000 in 2D — measured pairs vs the
//!    Poisson expectation and a pairs/s line (debug-build record; release
//!    numbers live in `benches/continuum_pairs.rs`).

use std::time::Instant;

use cmc_rs::{
    analyze, analyze_with_labels, candidate_pairs_from_configuration,
    candidate_pairs_with_cell_list, CandidatePairGraph, CellList, ComponentWorkspace,
    ContinuumError, EdgeActivity, EdgeId, GraphView, OrthorhombicCell, PairGraph,
    ParticleConfiguration, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Uniform random wrapped positions in an orthorhombic box.
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

/// Uniform random positions quantized to a dyadic grid, so uniform shifts
/// by lattice multiples wrap back to bit-identical coordinates.
fn quantized_positions<const D: usize>(n: usize, lengths: [f64; D], seed: u64) -> Vec<[f64; D]> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut position = [0.0; D];
            for axis in 0..D {
                position[axis] = (rng.random_range(0u64..512) as f64) * lengths[axis] / 512.0;
            }
            position
        })
        .collect()
}

fn flatten<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
    positions.iter().flat_map(|p| p.iter().copied()).collect()
}

/// Exhaustive all-pairs reference with an independently coded minimum
/// image (same rounding path as the production cell so only the SEARCH
/// differs, not the metric arithmetic).
fn bruteforce_pairs<const D: usize>(
    positions: &[[f64; D]],
    lengths: &[f64; D],
    cutoff: f64,
) -> Vec<[usize; 2]> {
    let cutoff_squared = cutoff * cutoff;
    let mut pairs = Vec::new();
    for first in 0..positions.len() {
        for second in (first + 1)..positions.len() {
            let mut distance_squared = 0.0;
            for axis in 0..D {
                let raw = positions[second][axis] - positions[first][axis];
                let image = raw - lengths[axis] * (raw * lengths[axis].recip()).round();
                distance_squared += image * image;
            }
            if distance_squared <= cutoff_squared {
                pairs.push([first, second]);
            }
        }
    }
    pairs
}

fn sorted_pairs(graph: &PairGraph) -> Vec<[usize; 2]> {
    let mut pairs = graph.pairs().to_vec();
    pairs.sort_unstable();
    pairs
}

#[test]
fn random_configurations_match_bruteforce_in_2d_and_3d() {
    // 2D: mean coordination ~ rho * pi * r_c^2 ~ 3.1; 3D ~ 3.4.
    for seed in [0x11u64, 0x2b, 0x5d] {
        let lengths = [6.0, 5.0];
        let positions = random_positions(150, lengths, seed);
        let cell = OrthorhombicCell::new(lengths).unwrap();
        let graph = CandidatePairGraph::new(&flatten(&positions), &cell, 1.3)
            .unwrap()
            .pairs()
            .unwrap();
        let expected = bruteforce_pairs(&positions, &lengths, 1.3);
        assert!(!expected.is_empty(), "2D fixture should contain pairs");
        assert_eq!(sorted_pairs(&graph), expected, "2D seed {seed:#x}");
        assert_eq!(graph.vertex_count(), 150);

        let lengths = [5.0, 4.5, 4.0];
        let positions = random_positions(200, lengths, seed);
        let cell = OrthorhombicCell::new(lengths).unwrap();
        let graph = CandidatePairGraph::new(&flatten(&positions), &cell, 1.1)
            .unwrap()
            .pairs()
            .unwrap();
        let expected = bruteforce_pairs(&positions, &lengths, 1.1);
        assert!(!expected.is_empty(), "3D fixture should contain pairs");
        assert_eq!(sorted_pairs(&graph), expected, "3D seed {seed:#x}");
    }
}

#[test]
fn exact_cutoff_boundary_is_inclusive() {
    // All-dyadic coordinates: squared distances are exact in f64.
    let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let positions = [
        [1.0, 1.0],  // 0: exactly cutoff 2.0 away from 1 (direct)
        [3.0, 1.0],  // 1
        [1.0, 4.0],  // 2: exactly cutoff 2.0 away from 3 (through image)
        [7.0, 4.0],  // 3
        [1.0, 6.0],  // 4: far from everything
        [4.25, 6.0], // 5
    ];
    let graph = CandidatePairGraph::new(&flatten(&positions), &cell, 2.0)
        .unwrap()
        .pairs()
        .unwrap();
    // (0,1) sits exactly at the cutoff directly; (2,3) through a periodic
    // image; (2,4) exactly again (dyadic arithmetic is exact).
    assert_eq!(sorted_pairs(&graph), vec![[0usize, 1usize], [2, 3], [2, 4]]);

    // A cutoff a hair below 2.0 drops both exactly-at-boundary pairs.
    let tighter = 2.0_f64 * (1.0 - 4.0 * f64::EPSILON);
    let graph = CandidatePairGraph::new(&flatten(&positions), &cell, tighter)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(graph.pair_count(), 0);
}

#[test]
fn periodic_image_only_pairs_are_found() {
    let cell = OrthorhombicCell::new([10.0, 10.0]).unwrap();
    let positions = [
        [0.3, 0.3], // direct difference to particle 1 is 9.5; image 0.5
        [9.8, 0.4], //
        [5.0, 5.0], // isolated
        [9.9, 9.9], // image-close to particle 0 (distance ~ sqrt(0.08))
    ];
    let graph = CandidatePairGraph::new(&flatten(&positions), &cell, 1.0)
        .unwrap()
        .pairs()
        .unwrap();
    // (1,3) is close only through BOTH coordinates' images.
    assert_eq!(sorted_pairs(&graph), vec![[0usize, 1usize], [0, 3], [1, 3]]);
    let expected = bruteforce_pairs(&positions, &[10.0, 10.0], 1.0);
    assert_eq!(sorted_pairs(&graph), expected);
}

#[test]
fn single_particle_is_valid_and_zero_particles_is_typed() {
    let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
    let graph = CandidatePairGraph::new(&[0.7, 0.9], &cell, 1.0)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(graph.vertex_count(), 1);
    assert_eq!(graph.pair_count(), 0);

    let error = CandidatePairGraph::<2>::new(&[], &cell, 1.0).unwrap_err();
    assert_eq!(error, ContinuumError::EmptyCoordinates { dimension: 2 });
}

#[test]
fn enumeration_is_deterministic_across_entry_points() {
    let lengths = [6.0, 5.0];
    let positions = random_positions(120, lengths, 0x91);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let flat = flatten(&positions);
    let adapter = CandidatePairGraph::new(&flat, &cell, 1.3).unwrap();
    let first = adapter.pairs().unwrap();
    let second = adapter.pairs().unwrap();
    assert_eq!(first.pairs(), second.pairs(), "same input, same order");

    let configuration =
        ParticleConfiguration::new(positions.clone(), vec![0; positions.len()], cell).unwrap();
    let from_configuration = candidate_pairs_from_configuration(&configuration, 1.3).unwrap();
    assert_eq!(
        first.pairs(),
        from_configuration.pairs(),
        "flat and configuration entry points agree"
    );
}

#[test]
fn translation_invariance_lattice_and_generic_shifts() {
    let lengths = [8.0, 8.0];
    let cell = OrthorhombicCell::new(lengths).unwrap();
    // 1.24 is deliberately NOT a multiple of the 1/64 quantization grid,
    // so no pair distance can coincide with the cutoff.
    let cutoff = 1.24;
    let positions = quantized_positions(200, lengths, 0x7f);
    let reference = CandidatePairGraph::new(&flatten(&positions), &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();

    // Exact lattice shift: wrapped coordinates are bit-identical, so the
    // canonical order must reproduce the identical sequence.
    let shifted: Vec<[f64; 2]> = positions
        .iter()
        .map(|p| [(p[0] + 8.0).rem_euclid(8.0), (p[1] - 16.0).rem_euclid(8.0)])
        .collect();
    assert_eq!(shifted, positions, "dyadic quantization wraps exactly");
    let lattice = CandidatePairGraph::new(&flatten(&shifted), &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(lattice.pairs(), reference.pairs());

    // Generic shift: relative geometry is unchanged, so the pair SET is
    // invariant even though particles cross cell boundaries. Guard against
    // boundary-flip sensitivity by requiring every pair distance to sit a
    // safe margin away from the cutoff.
    let shift = [0.37, 2.71];
    let shifted: Vec<[f64; 2]> = positions
        .iter()
        .map(|p| {
            let mut q = [0.0; 2];
            for axis in 0..2 {
                q[axis] = (p[axis] + shift[axis]).rem_euclid(lengths[axis]);
            }
            q
        })
        .collect();
    for first in 0..shifted.len() {
        for second in (first + 1)..shifted.len() {
            let mut distance_squared = 0.0;
            for axis in 0..2 {
                let raw = shifted[second][axis] - shifted[first][axis];
                let image = raw - lengths[axis] * (raw * lengths[axis].recip()).round();
                distance_squared += image * image;
            }
            let distance = distance_squared.sqrt();
            assert!(
                (distance - cutoff).abs() > 1e-9,
                "fixture pair sits on the cutoff; choose another seed"
            );
        }
    }
    let generic = CandidatePairGraph::new(&flatten(&shifted), &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();
    assert_eq!(sorted_pairs(&generic), sorted_pairs(&reference));
}

struct SiteMask {
    active: Vec<bool>,
}

impl VertexActivity for SiteMask {
    fn vertex_count(&self) -> usize {
        self.active.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.active[vertex.index()]
    }
}

struct AllEdges {
    count: usize,
}

impl EdgeActivity for AllEdges {
    fn edge_count(&self) -> usize {
        self.count
    }

    fn edge_active(&self, _edge: EdgeId) -> bool {
        true
    }
}

/// Independent adjacency-BFS component reference over the same pair set.
fn bfs_components(vertex_count: usize, pairs: &[[usize; 2]], active: &[bool]) -> Vec<Vec<usize>> {
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
    for &[left, right] in pairs {
        if active[left] && active[right] {
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
    }
    let mut visited = vec![false; vertex_count];
    let mut components = Vec::new();
    for seed in 0..vertex_count {
        if !active[seed] || visited[seed] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![seed];
        visited[seed] = true;
        while let Some(vertex) = stack.pop() {
            component.push(vertex);
            for &neighbor in &adjacency[vertex] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

#[test]
fn f3_analysis_on_pair_graph_matches_bfs_reference() {
    let lengths = [5.0, 4.5, 4.0];
    let positions = random_positions(150, lengths, 0xc3);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let graph = CandidatePairGraph::new(&flatten(&positions), &cell, 1.2)
        .unwrap()
        .pairs()
        .unwrap();
    assert!(graph.pair_count() > 100, "fixture should be well connected");

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x3a);
    let mask = SiteMask {
        active: (0..graph.vertex_count())
            .map(|_| rng.random::<f64>() < 0.6)
            .collect(),
    };
    let edges = AllEdges {
        count: graph.pair_count(),
    };

    let mut workspace = ComponentWorkspace::new();
    let result = analyze(&graph, &mask, &edges, &[], &mut workspace).unwrap();
    let summary = result.summary();

    let components = bfs_components(graph.vertex_count(), graph.pairs(), &mask.active);
    let active_count: usize = mask.active.iter().map(|&flag| usize::from(flag)).sum();
    let active_edges: usize = graph
        .pairs()
        .iter()
        .filter(|&&[left, right]| mask.active[left] && mask.active[right])
        .count();
    let second_moment: u128 = components
        .iter()
        .map(|component| (component.len() as u128).pow(2))
        .sum();
    assert_eq!(summary.active_vertex_count, active_count);
    assert_eq!(summary.active_edge_count, active_edges);
    assert_eq!(summary.component_count, components.len());
    assert_eq!(
        summary.largest_component_size,
        components.iter().map(Vec::len).max().unwrap_or(0)
    );
    assert_eq!(summary.raw_second_moment, second_moment);

    // Label parity: canonical labels are the lowest active index of each
    // component; the F4 tie policy picks the lowest identity among
    // equally-largest components.
    let labeled = analyze_with_labels(&graph, &mask, &edges, &[], &mut workspace).unwrap();
    let mut canonical_of: Vec<Option<usize>> = vec![None; graph.vertex_count()];
    for component in &components {
        let canonical = component[0];
        for &vertex in component {
            canonical_of[vertex] = Some(canonical);
        }
    }
    for (index, (canonical, &label)) in canonical_of
        .iter()
        .zip(labeled.labels().unwrap())
        .enumerate()
    {
        let expected = canonical.and_then(|canonical| graph.vertex_id(canonical));
        assert_eq!(label, expected, "vertex {index}");
    }
    let largest_identity = components
        .iter()
        .filter(|c| c.len() == summary.largest_component_size)
        .map(|c| c[0])
        .min()
        .and_then(|canonical| graph.vertex_id(canonical));
    assert_eq!(summary.largest_component, largest_identity);
}

#[test]
fn prepared_cell_list_matches_fresh_build_and_rejects_stale_lists() {
    let lengths = [6.0, 5.0];
    let positions = random_positions(100, lengths, 0x1f);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let cutoff = 1.3;
    let configuration =
        ParticleConfiguration::new(positions.clone(), vec![0; positions.len()], cell).unwrap();
    let list = CellList::new(&configuration, cutoff * cutoff).unwrap();

    let via_list = candidate_pairs_with_cell_list(&configuration, &list).unwrap();
    let fresh = candidate_pairs_from_configuration(&configuration, cutoff).unwrap();
    assert_eq!(via_list.pairs(), fresh.pairs());

    let moved = random_positions(100, lengths, 0x2e);
    let stale = ParticleConfiguration::new(moved, vec![0; 100], cell).unwrap();
    assert!(matches!(
        candidate_pairs_with_cell_list(&stale, &list),
        Err(ContinuumError::CellListInconsistent(_))
    ));
}

#[test]
fn invalid_domains_are_typed_rejections() {
    let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
    for cutoff in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        // NaN payloads compare through `matches!`: derived equality on the
        // f64 field would report NaN != NaN.
        assert!(
            matches!(
                CandidatePairGraph::new(&[1.0, 1.0], &cell, cutoff).unwrap_err(),
                ContinuumError::InvalidCutoff { .. }
            ),
            "cutoff {cutoff} must be rejected"
        );
    }
    assert_eq!(
        CandidatePairGraph::<2>::new(&[1.0, 2.0, 3.0], &cell, 1.0).unwrap_err(),
        ContinuumError::CoordinateCountMismatch {
            length: 3,
            dimension: 2
        }
    );
    assert_eq!(
        CandidatePairGraph::<2>::new(&[1.0, f64::NAN], &cell, 1.0).unwrap_err(),
        ContinuumError::NonFiniteCoordinate {
            particle: 0,
            axis: 1
        }
    );
    // Geometry limits surface at enumeration time (the reused CellList
    // owns that check): cutoff 2.5 exceeds half the shortest side 2.0.
    let error = CandidatePairGraph::new(&[1.0, 1.0], &cell, 2.5)
        .unwrap()
        .pairs()
        .unwrap_err();
    assert!(matches!(error, ContinuumError::CutoffTooLarge { .. }));
}

#[test]
fn throughput_record_10k_particles_2d() {
    let n = 10_000_usize;
    let lengths = [100.0, 100.0];
    let positions = random_positions(n, lengths, 0x5eed);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let cutoff = 1.0;
    let start = Instant::now();
    let graph = CandidatePairGraph::new(&flatten(&positions), &cell, cutoff)
        .unwrap()
        .pairs()
        .unwrap();
    let elapsed = start.elapsed();
    let measured = graph.pair_count();
    // Poisson expectation: N * (rho * pi * r_c^2) / 2 with rho = N / V = 1.
    let expected = (n as f64) * std::f64::consts::PI * cutoff * cutoff / 2.0;
    eprintln!(
        "CONTINUUM_PAIRS n={n} cutoff={cutoff} pairs={measured} \
         poisson_expectation={expected:.0} ratio={:.4} pairs_per_second={:.3e} \
         (debug-build record; release numbers in benches/continuum_pairs.rs)",
        measured as f64 / expected,
        measured as f64 / elapsed.as_secs_f64(),
    );
    assert!(
        (measured as f64 / expected - 1.0).abs() < 0.1,
        "measured pair count {measured} vs Poisson expectation {expected:.0}"
    );
}
