//! W2 winding/wrapping validation: hand-constructed torus cycles, an
//! independent finite-cover lift-BFS reference, gauge/reversal/labeling
//! invariances, typed overflow rejections, full bond enumeration on 2x2 and
//! 3x3 tori, and scheduler-level wrapping observables with 16-seed z-scores.

use std::sync::OnceLock;

use carlo_rs::{Context, FromParams, Params, Run, RunConfig, RunId, TaskId};
use cmc_rs::{
    analyze, analyze_winding, build_chain_with_embedding, build_hypercubic_with_embedding,
    build_square, build_square_with_embedding, Bond, BondType, BoundaryQuery, ComponentWorkspace,
    CsrLattice, EdgeActivity, EdgeId, GraphView, LatticeEmbedding, ObservablePlan,
    PeriodicEmbedding, StaticLaw, StaticObservable, StaticPercolationError, StaticPercolationMC,
    UndirectedGraphView, VertexActivity, VertexId, WindingError, WindingResult, WindingWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::zscore_seed_count;

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

fn all_active(lattice: &CsrLattice) -> (VertexFlags, EdgeFlags) {
    (
        VertexFlags(vec![true; lattice.vertex_count()]),
        EdgeFlags(vec![true; lattice.edge_count()]),
    )
}

fn production<'a>(
    lattice: &CsrLattice,
    embedding: &LatticeEmbedding,
    vertices: &VertexFlags,
    edges: &EdgeFlags,
    workspace: &'a mut WindingWorkspace,
) -> WindingResult<'a> {
    analyze_winding(lattice, embedding, vertices, edges, workspace)
        .unwrap_or_else(|error| panic!("production winding analysis failed: {error}"))
}

// ── Independent reference: finite-cover lift-BFS ─────────────────────────
//
// Replicates the torus into a 5^D cover of cells (radius 2), BFS's each
// component's lift from its seed pinned in the center cell, and reads the
// winding subgroup off the cells in which the seed vertex reappears. Rank is
// then computed by an independent method: plain f64 Gaussian elimination on
// the collected cell offsets (the real rank of integer vectors equals the
// rank of the subgroup they generate). This shares no algorithmic path with
// the production weighted union-find. Radius 2 provably suffices on the
// tested tori: a simple cycle winding k times around one axis of an L-torus
// needs at least k*L edges, and the 3x3 (9 sites) and 2x2 (4 sites) tori
// cannot host simple cycles with any winding component above 2.

const REFERENCE_RADIUS: i32 = 2;

/// (canonical vertex index, size, winding rank) per component, ascending
/// canonical order — the identity policy is shared bookkeeping, not the
/// algorithm under test.
struct LiftReference(Vec<(usize, usize, usize)>);

fn lift_reference(
    lattice: &CsrLattice,
    embedding: &LatticeEmbedding,
    vertices: &VertexFlags,
    edges: &EdgeFlags,
    dimensions: usize,
) -> LiftReference {
    let mut adjacency = vec![Vec::<(usize, Vec<i32>)>::new(); lattice.vertex_count()];
    for edge_index in 0..lattice.edge_count() {
        let edge = lattice.edge_id(edge_index).expect("dense edge index");
        if !edges.edge_active(edge) {
            continue;
        }
        let [source, target] = lattice.edge_endpoints(edge);
        if !vertices.vertex_active(source) || !vertices.vertex_active(target) {
            continue;
        }
        let forward: Vec<i32> = embedding
            .edge_displacement(edge, source)
            .expect("source endpoint")
            .iter()
            .collect();
        let backward: Vec<i32> = forward.iter().map(|&component| -component).collect();
        adjacency[source.index()].push((target.index(), forward));
        adjacency[target.index()].push((source.index(), backward));
    }

    let mut seen = vec![false; lattice.vertex_count()];
    let mut components = Vec::new();
    for seed in 0..lattice.vertex_count() {
        if seen[seed] || !vertices.0[seed] {
            continue;
        }
        seen[seed] = true;
        let mut members = 1usize;
        let mut stack = vec![seed];
        while let Some(vertex) = stack.pop() {
            for &(neighbor, _) in &adjacency[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    members += 1;
                    stack.push(neighbor);
                }
            }
        }

        // Lift BFS over the (2R+1)^D cover; steps leaving the box are
        // skipped (radius argument above makes that lossless for these tori).
        let side = (2 * REFERENCE_RADIUS + 1) as usize;
        let cell_index = |cell: &[i32]| -> usize {
            cell.iter()
                .enumerate()
                .map(|(axis, &offset)| (offset + REFERENCE_RADIUS) as usize * side.pow(axis as u32))
                .sum()
        };
        let total_cells = side.pow(dimensions as u32);
        let mut lifted_seen = vec![false; lattice.vertex_count() * total_cells];
        let mut offsets = Vec::new();
        lifted_seen[seed * total_cells + cell_index(&vec![0; dimensions])] = true;
        let mut queue = vec![(seed, vec![0i32; dimensions])];
        let mut head = 0;
        while head < queue.len() {
            let (vertex, cell) = queue[head].clone();
            head += 1;
            if vertex == seed && cell.iter().any(|&offset| offset != 0) {
                offsets.push(cell.clone());
            }
            for (neighbor, displacement) in &adjacency[vertex] {
                let mut next_cell = cell.clone();
                let mut inside = true;
                for (axis, &step) in displacement.iter().enumerate() {
                    next_cell[axis] += step;
                    inside &= next_cell[axis].abs() <= REFERENCE_RADIUS;
                }
                if !inside {
                    continue;
                }
                let slot = *neighbor * total_cells + cell_index(&next_cell);
                if !lifted_seen[slot] {
                    lifted_seen[slot] = true;
                    queue.push((*neighbor, next_cell));
                }
            }
        }
        components.push((seed, members, real_rank(&offsets, dimensions)));
    }
    LiftReference(components)
}

/// Rank over the reals of an integer vector set (equals the subgroup rank).
fn real_rank(vectors: &[Vec<i32>], dimensions: usize) -> usize {
    let mut rows: Vec<Vec<f64>> = vectors
        .iter()
        .map(|vector| vector.iter().map(|&component| component as f64).collect())
        .collect();
    let mut rank = 0;
    for column in 0..dimensions {
        let pivot = (rank..rows.len())
            .filter(|&row| rows[row][column].abs() > 1e-9)
            .max_by(|&a, &b| {
                rows[a][column]
                    .abs()
                    .partial_cmp(&rows[b][column].abs())
                    .unwrap_or(core::cmp::Ordering::Equal)
            });
        let Some(pivot) = pivot else { continue };
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

/// Test-side canonical row Hermite normal form (small vectors): reduces a
/// generating set to a unique echelon form so two generator sets can be
/// compared as lattices, independent of each set's basis choice.
fn canonical_lattice(vectors: &[Vec<i32>], dimensions: usize) -> Vec<Vec<i64>> {
    let mut rows: Vec<Vec<i64>> = vectors
        .iter()
        .map(|vector| vector.iter().map(|&component| component as i64).collect())
        .collect();
    if rows.is_empty() {
        return rows;
    }
    let mut rank = 0usize;
    for column in 0..dimensions {
        if rank == rows.len() || rows[rank..].iter().all(|row| row[column] == 0) {
            continue;
        }
        let pivot = (rank..rows.len())
            .find(|&row| rows[row][column] != 0)
            .unwrap();
        rows.swap(rank, pivot);
        for other in (rank + 1)..rows.len() {
            while rows[other][column] != 0 {
                let quotient = rows[rank][column] / rows[other][column];
                if quotient != 0 {
                    let other_row = rows[other].clone();
                    for (value, other_value) in rows[rank].iter_mut().zip(&other_row) {
                        *value -= quotient * other_value;
                    }
                }
                rows.swap(rank, other);
            }
        }
        if rows[rank][column] < 0 {
            for component in &mut rows[rank] {
                *component = -*component;
            }
        }
        for above in 0..rank {
            let quotient = rows[above][column] / rows[rank][column];
            let pivot_row = rows[rank].clone();
            for (value, pivot_value) in rows[above].iter_mut().zip(&pivot_row) {
                *value -= quotient * pivot_value;
            }
        }
        rank += 1;
    }
    rows.truncate(rank);
    rows
}

/// Per-component generator lists (components are in identical canonical
/// order across the compared runs, so zip aligns them).
fn component_generators(result: &WindingResult) -> Vec<Vec<Vec<i32>>> {
    (0..result.components().len())
        .map(|index| {
            result
                .generators(index)
                .map(|generator| generator.to_vec())
                .collect::<Vec<Vec<i32>>>()
        })
        .collect()
}

fn assert_same_lattice(left: &[Vec<i32>], right: &[Vec<i32>], dimensions: usize, label: &str) {
    assert_eq!(
        canonical_lattice(left, dimensions),
        canonical_lattice(right, dimensions),
        "{label}: generated winding lattices differ"
    );
}

fn assert_matches_reference(
    lattice: &CsrLattice,
    embedding: &LatticeEmbedding,
    vertices: &VertexFlags,
    edges: &EdgeFlags,
    dimensions: usize,
    label: &str,
) {
    let mut workspace = WindingWorkspace::new();
    let result = production(lattice, embedding, vertices, edges, &mut workspace);
    let reference = lift_reference(lattice, embedding, vertices, edges, dimensions);
    assert_eq!(
        result.components().len(),
        reference.0.len(),
        "{label}: component count"
    );
    for (info, (canonical, size, rank)) in result.components().iter().zip(reference.0.iter()) {
        assert_eq!(info.component().index(), *canonical, "{label}: identity");
        assert_eq!(info.size(), *size, "{label}: size");
        assert_eq!(info.rank(), *rank, "{label}: winding rank");
        assert_eq!(
            result.wraps(result.component_index(info.component()).unwrap()),
            *rank > 0,
            "{label}: wraps is rank > 0"
        );
    }
}

// ── Hand-constructed cycles ──────────────────────────────────────────────

#[test]
fn interior_plaquette_is_contractible() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    // Activate only the (0,0)-(1,0)-(1,1)-(0,1) plaquette sites; the loop
    // stays inside one fundamental cell and can never wind.
    let mut vertices = VertexFlags(vec![false; 9]);
    for site in [0, 1, 3, 4] {
        vertices.0[site] = true;
    }
    let edges = EdgeFlags(vec![true; lattice.edge_count()]);
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().component_count, 1);
    assert_eq!(
        result.summary().wrapping_component_count,
        0,
        "contractible loop"
    );
    assert_eq!(result.summary().largest_winding_rank, 0);
    assert_matches_reference(&lattice, &embedding, &vertices, &edges, 2, "plaquette");
}

#[test]
fn single_axis_ring_winds_with_exact_unit_generator() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    // Row-0 ring: sites 0, 1, 2 and their three horizontal bonds (the third
    // wraps).
    let mut vertices = VertexFlags(vec![false; 9]);
    for site in 0..3 {
        vertices.0[site] = true;
    }
    let mut edges = EdgeFlags(vec![false; lattice.edge_count()]);
    for (index, bond) in lattice.edges.iter().enumerate() {
        if bond.source < 3 && bond.target < 3 {
            edges.0[index] = true;
        }
    }
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().wrapping_component_count, 1);
    let info = result.components()[0];
    assert_eq!(
        (info.component().index(), info.size(), info.rank()),
        (0, 3, 1)
    );
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 0][..]],
        "single-axis wrap generator is exactly the x unit vector"
    );
    assert_matches_reference(&lattice, &embedding, &vertices, &edges, 2, "row ring");
}

#[test]
fn double_axis_wrap_reports_rank_two_unit_generators() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    // Row-0 ring (sites 0,1,2) plus column-0 ring (sites 0,3,6), sharing
    // site 0: one component winding both axes.
    let mut vertices = VertexFlags(vec![false; 9]);
    for site in [0, 1, 2, 3, 6] {
        vertices.0[site] = true;
    }
    let mut edges = EdgeFlags(vec![false; lattice.edge_count()]);
    for (index, bond) in lattice.edges.iter().enumerate() {
        let ring = matches!(
            (bond.source, bond.target),
            (0, 1) | (1, 2) | (2, 0) | (0, 3) | (3, 6) | (6, 0)
        );
        edges.0[index] = ring;
    }
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().wrapping_component_count, 1);
    let info = result.components()[0];
    assert_eq!((info.size(), info.rank()), (5, 2));
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 0][..], &[0, 1][..]],
        "two-axis wrap reports both unit generators"
    );
    assert_matches_reference(&lattice, &embedding, &vertices, &edges, 2, "cross rings");
}

#[test]
fn diagonal_wrap_is_rank_one_with_diagonal_generator() {
    // Hand-built 4-cycle whose lift accumulates (1, 1): winding rank 1 like a
    // single-axis wrap, but the generator distinguishes them.
    let lattice = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
            Bond::new(3, 0, BondType::Generic, 1.0),
        ],
    );
    let embedding =
        LatticeEmbedding::try_from_edge_displacements(&lattice, 2, &[0, 0, 0, 1, 0, 0, 1, 0])
            .unwrap();
    let (vertices, edges) = all_active(&lattice);
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().largest_winding_rank, 1);
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 1][..]],
        "diagonal wrap generator is (1, 1), not an axis unit vector"
    );
    assert_matches_reference(&lattice, &embedding, &vertices, &edges, 2, "diagonal");
}

#[test]
fn two_ring_parallel_bonds_wrap_only_through_the_cycle() {
    // L=2 ring: two parallel bonds, one in-cell (displacement 0) and one
    // wrapping. This pins the (edge, endpoint) key path for parallel edges.
    let (lattice, embedding) = build_chain_with_embedding(2, true).unwrap();
    let (vertices, edges) = all_active(&lattice);
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().component_count, 1);
    assert_eq!(result.summary().wrapping_component_count, 1);
    assert_eq!(result.summary().largest_winding_rank, 1);

    // With only the wrapping bond the component is an open path: connected,
    // but a path has no cycles and cannot wind.
    let mut wrap_only = edges.clone();
    wrap_only.0 = vec![false, true];
    let result = production(&lattice, &embedding, &vertices, &wrap_only, &mut workspace);
    assert_eq!(
        result.summary().component_count,
        1,
        "still connected through the wrap bond"
    );
    assert_eq!(result.summary().wrapping_component_count, 0);
    assert_matches_reference(
        &lattice,
        &embedding,
        &vertices,
        &wrap_only,
        1,
        "open 2-path",
    );
}

#[test]
fn cubic_torus_full_lattice_winds_all_three_axes() {
    let (lattice, embedding) = build_hypercubic_with_embedding(
        &[2, 3, 4],
        &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
        true,
    )
    .unwrap();
    let (vertices, edges) = all_active(&lattice);
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
    assert_eq!(result.summary().largest_winding_rank, 3);
    assert_eq!(
        result.generators(0).collect::<Vec<_>>(),
        [&[1, 0, 0][..], &[0, 1, 0][..], &[0, 0, 1][..]]
    );
}

// ── Invariances ──────────────────────────────────────────────────────────

#[test]
fn winding_is_invariant_under_gauge_shift_edge_reversal_and_labeling() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    let dimensions = 2;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5749_4e44);
    let stored: Vec<i32> = (0..lattice.edge_count())
        .flat_map(|index| {
            let edge = lattice.edge_id(index).unwrap();
            let [source, _target] = lattice.edge_endpoints(edge);
            embedding
                .edge_displacement(edge, source)
                .unwrap()
                .iter()
                .collect::<Vec<i32>>()
        })
        .collect();
    let endpoints: Vec<[usize; 2]> = (0..lattice.edge_count())
        .map(|index| {
            let edge = lattice.edge_id(index).unwrap();
            let [source, target] = lattice.edge_endpoints(edge);
            [source.index(), target.index()]
        })
        .collect();

    for sample in 0..8 {
        let vertices = VertexFlags((0..9).map(|_| rng.random::<bool>()).collect());
        let edges = EdgeFlags((0..18).map(|_| rng.random::<bool>()).collect());
        let mut workspace = WindingWorkspace::new();
        let (baseline_summary, baseline_generators) = {
            let baseline = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
            (baseline.summary(), component_generators(&baseline))
        };
        let _ = sample;

        // (a) Representative/origin shift: rebuild the table under a gauge
        // transform d'(s->t) = d(s->t) + phi(t) - phi(s). Cycle sums are
        // unchanged (phi telescopes), so the winding subgroup is identical.
        let phi: Vec<[i32; 2]> = (0..9usize)
            .map(|vertex| [((vertex * 7) % 5) as i32 - 2, ((vertex * 3) % 5) as i32 - 2])
            .collect();
        let shifted: Vec<i32> = endpoints
            .iter()
            .enumerate()
            .flat_map(|(index, &[source, target])| {
                (0..dimensions)
                    .map(|axis| {
                        stored[index * dimensions + axis] + phi[target][axis] - phi[source][axis]
                    })
                    .collect::<Vec<i32>>()
            })
            .collect();
        let gauge_embedding =
            LatticeEmbedding::try_from_edge_displacements(&lattice, dimensions, &shifted).unwrap();
        let (gauge_summary, gauge_generators) = {
            let gauge_result = production(
                &lattice,
                &gauge_embedding,
                &vertices,
                &edges,
                &mut workspace,
            );
            (gauge_result.summary(), component_generators(&gauge_result))
        };
        assert_eq!(
            gauge_summary, baseline_summary,
            "gauge/representative shift must not change winding facts"
        );
        for (baseline_gens, gauge_gens) in baseline_generators.iter().zip(&gauge_generators) {
            assert_same_lattice(baseline_gens, gauge_gens, dimensions, "gauge shift");
        }

        // (b) Edge direction reversal: negate every stored vector. The
        // generated subgroup is negation-closed, hence identical.
        let reversed: Vec<i32> = stored.iter().map(|&component| -component).collect();
        let reversed_embedding =
            LatticeEmbedding::try_from_edge_displacements(&lattice, dimensions, &reversed).unwrap();
        let (reversed_summary, reversed_generators) = {
            let reversed_result = production(
                &lattice,
                &reversed_embedding,
                &vertices,
                &edges,
                &mut workspace,
            );
            (
                reversed_result.summary(),
                component_generators(&reversed_result),
            )
        };
        assert_eq!(reversed_summary, baseline_summary);
        for (baseline_gens, reversed_gens) in baseline_generators.iter().zip(&reversed_generators) {
            assert_same_lattice(baseline_gens, reversed_gens, dimensions, "edge reversal");
        }
    }

    // (c) Origin/label translation: relabel the winding 4-ring cyclically;
    // the same physical ring must wind identically.
    for shift in 0..4 {
        let relabel = |site: usize| (site + shift) % 4;
        let lattice = CsrLattice::from_edges(
            4,
            vec![
                Bond::new(relabel(0), relabel(1), BondType::Generic, 1.0),
                Bond::new(relabel(1), relabel(2), BondType::Generic, 1.0),
                Bond::new(relabel(2), relabel(3), BondType::Generic, 1.0),
                Bond::new(relabel(3), relabel(0), BondType::Generic, 1.0),
            ],
        );
        let embedding =
            LatticeEmbedding::try_from_edge_displacements(&lattice, 1, &[0, 0, 0, 1]).unwrap();
        let (vertices, edges) = all_active(&lattice);
        let mut workspace = WindingWorkspace::new();
        let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
        assert_eq!(result.summary().largest_winding_rank, 1, "shift {shift}");
        assert_eq!(
            result.generators(0).collect::<Vec<_>>(),
            [&[1][..]],
            "shift {shift}: generator is the unit vector"
        );
    }
}

// ── Wrapping on open/non-embedded graphs is an explicit error ────────────

#[test]
fn wrapping_requests_without_an_embedding_are_explicit_errors() {
    let open = build_square(2, 2, false);
    let plan = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
    let error = StaticPercolationMC::new(open.clone(), StaticLaw::bond(0.5).unwrap(), plan.clone())
        .err()
        .expect("open boundary + winding plan must be rejected");
    assert!(matches!(
        error,
        StaticPercolationError::WindingWithoutEmbedding
    ));

    // The parameter path rejects the same request loudly.
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", 2);
    params.set("Ly", 2);
    params.set("p", 0.5);
    params.set("observables", "wrapping");
    let error =
        StaticPercolationMC::from_params(&params, &mut Xoshiro256PlusPlus::seed_from_u64(1))
            .err()
            .expect("default-open parameter topology + winding must be rejected");
    assert!(
        error.to_string().contains("embedded periodic topology"),
        "{error}"
    );

    // The W1 builders refuse open embeddings outright.
    assert!(build_chain_with_embedding(4, false).is_err());
}

// ── Scheduler-level physical sampling ────────────────────────────────────

/// Run one directly constructed adapter through `Run::from_parts`.
fn run_adapter(
    lattice: &CsrLattice,
    embedding: &LatticeEmbedding,
    law: StaticLaw,
    plan: ObservablePlan,
    seed: u64,
    sweeps: u64,
    binsize: usize,
) -> carlo_rs::Results {
    let adapter =
        StaticPercolationMC::with_embedding(lattice.clone(), embedding.clone(), law, plan).unwrap();
    let context = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, binsize);
    let mut run = Run::from_parts(
        context,
        adapter,
        TaskId::new(0),
        RunId::new(0),
        RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: sweeps,
            binsize,
            base_seed: seed,
            ..Default::default()
        },
    );
    run.run(sweeps);
    run.finalize(seed)
}

#[test]
fn bond_p_one_always_wraps_and_p_zero_never_does() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    let make_plan = |probability: f64| {
        (
            StaticLaw::bond(probability).unwrap(),
            ObservablePlan::new(
                [StaticObservable::Wrapping, StaticObservable::WindingRank],
                Vec::new(),
            ),
        )
    };
    let (law, plan) = make_plan(1.0);
    let saturated = run_adapter(&lattice, &embedding, law, plan, 3, 8, 8);
    assert_eq!(saturated.get("Wrapping").unwrap().mean, 1.0);
    assert_eq!(saturated.get("WindingRank").unwrap().mean, 2.0);

    let (law, plan) = make_plan(0.0);
    let empty = run_adapter(&lattice, &embedding, law, plan, 3, 8, 8);
    assert_eq!(empty.get("Wrapping").unwrap().mean, 0.0);
    assert_eq!(empty.get("WindingRank").unwrap().mean, 0.0);
}

#[test]
fn site_p_one_wraps_through_both_axes() {
    let (lattice, embedding) = build_square_with_embedding(2, 2, true).unwrap();
    let plan = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
    let results = run_adapter(
        &lattice,
        &embedding,
        StaticLaw::site(1.0).unwrap(),
        plan,
        5,
        4,
        4,
    );
    assert_eq!(results.get("Wrapping").unwrap().mean, 1.0);
}

#[test]
fn fixed_seed_reproduces_winding_observables_bitwise() {
    let (lattice, embedding) = build_square_with_embedding(3, 3, true).unwrap();
    let query = BoundaryQuery::new(&lattice, &[0], &[1]).unwrap();
    let plan = || {
        ObservablePlan::new(
            [StaticObservable::Wrapping, StaticObservable::WindingRank],
            vec![query.clone()],
        )
    };
    let run = |seed: u64| {
        run_adapter(
            &lattice,
            &embedding,
            StaticLaw::bond(0.5).unwrap(),
            plan(),
            seed,
            64,
            8,
        )
    };
    let first = run(11);
    let second = run(11);
    for name in ["Wrapping", "WindingRank", "BoundaryCrossing0"] {
        let left = first.get(name).expect(name);
        let right = second.get(name).expect(name);
        assert_eq!(
            (left.mean, left.n_bins),
            (right.mean, right.n_bins),
            "{name}"
        );
    }
    let shifted = run(12);
    assert_ne!(
        first.get("Wrapping").unwrap().stderr,
        shifted.get("Wrapping").unwrap().stderr
    );
}

// ── Full enumeration on 2x2 and 3x3 bond tori ───────────────────────────

/// Enumerates every bond configuration of an `size` x `size` torus, checks
/// each against the lift-BFS reference component by component, and returns
/// `(P(wrapping), mean largest-component winding rank)`.
fn enumerated_wrapping_probability(size: usize) -> (f64, f64) {
    let (lattice, embedding) = build_square_with_embedding(size, size, true).unwrap();
    let edge_count = lattice.edge_count();
    let vertices = VertexFlags(vec![true; lattice.vertex_count()]);
    let mut workspace = WindingWorkspace::new();
    let mut wrapping = 0u64;
    let mut rank_sum = 0u64;
    for mask in 0..1u64 << edge_count {
        let edges = EdgeFlags(
            (0..edge_count)
                .map(|edge| mask & (1 << edge) != 0)
                .collect(),
        );
        let result = production(&lattice, &embedding, &vertices, &edges, &mut workspace);
        let reference = lift_reference(&lattice, &embedding, &vertices, &edges, 2);
        assert_eq!(
            result.components().len(),
            reference.0.len(),
            "mask {mask:#x}: component count"
        );
        for (info, (canonical, component_size, rank)) in
            result.components().iter().zip(reference.0.iter())
        {
            assert_eq!(
                (info.component().index(), info.size(), info.rank()),
                (*canonical, *component_size, *rank),
                "mask {mask:#x}"
            );
        }
        wrapping += u64::from(result.summary().wrapping_component_count > 0);
        rank_sum += result.summary().largest_winding_rank as u64;
    }
    let total = 1u64 << edge_count;
    (
        wrapping as f64 / total as f64,
        rank_sum as f64 / total as f64,
    )
}

#[test]
fn two_by_two_torus_enumeration_matches_reference_and_hand_checks() {
    // 8 bonds -> 256 configurations, each parity-checked inside.
    let (probability, mean_rank) = enumerated_wrapping_probability(2);
    eprintln!("2x2 bond torus: P(wrapping) = {probability}, mean largest rank = {mean_rank}");
    // Exact rational: 256 configurations in total.
    assert_eq!(probability * 256.0, (probability * 256.0).round());
    assert!(probability > 0.0 && probability < 1.0);

    // Hand check: the four horizontal bonds (one in-cell + one wrapping per
    // row) form two disjoint 2-cycles, each winding once around x. The two
    // wrapping bonds alone form open pairs: connected, but a path has no
    // cycles and cannot wind.
    let (lattice, embedding) = build_square_with_embedding(2, 2, true).unwrap();
    let mut horizontal = EdgeFlags(vec![false; 8]);
    let mut wrap_only = EdgeFlags(vec![false; 8]);
    for (index, bond) in lattice.edges.iter().enumerate() {
        let pair = (bond.source.min(bond.target), bond.source.max(bond.target));
        let is_horizontal = matches!(pair, (0, 1) | (2, 3));
        let edge = lattice.edge_id(index).unwrap();
        let source = lattice.vertex_id(bond.source).unwrap();
        let displacement = embedding.edge_displacement(edge, source).unwrap();
        let wraps_x = displacement.axis(0) != Some(0) && displacement.axis(1) == Some(0);
        horizontal.0[index] = is_horizontal;
        wrap_only.0[index] = is_horizontal && wraps_x;
    }
    assert_eq!(
        wrap_only.0.iter().filter(|&&flag| flag).count(),
        2,
        "exactly the two horizontal wrap bonds are selected"
    );
    let vertices = VertexFlags(vec![true; 4]);
    let mut workspace = WindingWorkspace::new();
    let result = production(&lattice, &embedding, &vertices, &horizontal, &mut workspace);
    assert_eq!(result.summary().component_count, 2, "two disjoint 2-cycles");
    assert_eq!(result.summary().wrapping_component_count, 2);
    for index in 0..2 {
        assert_eq!(result.components()[index].rank(), 1);
        assert_eq!(result.components()[index].size(), 2);
        let generator = result.generators(index).next().unwrap();
        assert_eq!(generator[1], 0, "winding is purely horizontal");
    }
    let result = production(&lattice, &embedding, &vertices, &wrap_only, &mut workspace);
    assert_eq!(result.summary().component_count, 2, "two open pairs");
    assert_eq!(
        result.summary().wrapping_component_count,
        0,
        "paths do not wind"
    );
}

#[test]
fn three_by_three_torus_enumeration_matches_reference() {
    // 18 bonds -> 262144 configurations, each parity-checked against the
    // lift-BFS reference.
    let (probability, mean_rank) = enumerated_wrapping_probability(3);
    eprintln!("3x3 bond torus: P(wrapping) = {probability}, mean largest rank = {mean_rank}");
    assert!(probability > 0.0 && probability < 1.0);
}

// ── 16-seed z-scores against the enumerated 2x2 reference ────────────────

fn two_by_two_reference() -> &'static (f64, f64) {
    static REFERENCE: OnceLock<(f64, f64)> = OnceLock::new();
    REFERENCE.get_or_init(|| enumerated_wrapping_probability(2))
}

#[test]
fn scheduler_wrapping_probability_matches_enumeration_within_z_scores() {
    let (exact_probability, exact_rank) = *two_by_two_reference();
    let seeds = zscore_seed_count(16);
    let sweeps = 4000u64;
    let mut z_wrapping = Vec::new();
    let mut z_rank = Vec::new();
    for seed in 0..seeds as u64 {
        let (lattice, embedding) = build_square_with_embedding(2, 2, true).unwrap();
        let plan = ObservablePlan::new(
            [StaticObservable::Wrapping, StaticObservable::WindingRank],
            Vec::new(),
        );
        let results = run_adapter(
            &lattice,
            &embedding,
            StaticLaw::bond(0.5).unwrap(),
            plan,
            0x5A5A_0000 + seed,
            sweeps,
            100,
        );
        for (name, exact, z_target) in [
            ("Wrapping", exact_probability, &mut z_wrapping),
            ("WindingRank", exact_rank, &mut z_rank),
        ] {
            let estimate = results
                .get(name)
                .unwrap_or_else(|| panic!("missing {name}"));
            z_target.push((estimate.mean - exact) / estimate.stderr.max(1e-12));
        }
    }
    for (label, z_scores) in [("Wrapping", &z_wrapping), ("WindingRank", &z_rank)] {
        for (index, &z) in z_scores.iter().enumerate() {
            assert!(z.abs() < 4.0, "{label} seed {index}: z = {z}");
        }
        let mean_z = z_scores.iter().sum::<f64>() / seeds as f64;
        assert!(mean_z.abs() < 1.5, "{label} mean z = {mean_z}");
    }
    assert!(
        z_wrapping.iter().any(|&z| z > 0.0) && z_wrapping.iter().any(|&z| z < 0.0),
        "one-sided wrapping bias"
    );
}

// ── Cross-consistency with F3 ────────────────────────────────────────────

#[test]
fn winding_summary_agrees_with_f3_component_summary() {
    let (lattice, embedding) = build_square_with_embedding(4, 4, true).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x465f_3346);
    let mut winding_workspace = WindingWorkspace::new();
    let mut component_workspace = ComponentWorkspace::new();
    component_workspace
        .prepare(lattice.vertex_count(), 0, false)
        .unwrap();
    for _ in 0..32 {
        let vertices = VertexFlags((0..16).map(|_| rng.random::<bool>()).collect());
        let edges = EdgeFlags((0..32).map(|_| rng.random::<bool>()).collect());
        let winding = production(
            &lattice,
            &embedding,
            &vertices,
            &edges,
            &mut winding_workspace,
        );
        let components =
            analyze(&lattice, &vertices, &edges, &[], &mut component_workspace).unwrap();
        let summary = components.summary();
        assert_eq!(
            winding.summary().active_vertex_count,
            summary.active_vertex_count
        );
        assert_eq!(
            winding.summary().active_edge_count,
            summary.active_edge_count
        );
        assert_eq!(winding.summary().component_count, summary.component_count);
        assert_eq!(
            winding.summary().largest_component,
            summary.largest_component
        );
        assert_eq!(
            winding.summary().largest_component_size,
            summary.largest_component_size
        );
    }
}

// ── Adversarial magnitudes ───────────────────────────────────────────────

#[test]
fn adversarial_displacement_tables_return_typed_overflow_errors() {
    // 8-cycle of 6*10^8 displacements: cycle sums reach 4.8*10^9 > i32::MAX.
    // The analyzer must return a typed overflow error attributed to an edge,
    // never a silently wrapped value; the adapter's static gate rejects the
    // same table at construction instead of panicking inside measure().
    let length = 8usize;
    let lattice = CsrLattice::from_edges(
        length,
        (0..length - 1)
            .map(|site| Bond::new(site, site + 1, BondType::Generic, 1.0))
            .chain([Bond::new(length - 1, 0, BondType::Generic, 1.0)])
            .collect(),
    );
    let table: Vec<i32> = (0..length).flat_map(|_| [600_000_000i32]).collect();
    let embedding = LatticeEmbedding::try_from_edge_displacements(&lattice, 1, &table).unwrap();
    let (vertices, edges) = all_active(&lattice);
    let mut workspace = WindingWorkspace::new();
    let error =
        analyze_winding(&lattice, &embedding, &vertices, &edges, &mut workspace).unwrap_err();
    assert!(
        matches!(error, WindingError::DisplacementOverflow { .. }),
        "{error}"
    );

    let plan = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
    assert!(matches!(
        StaticPercolationMC::with_embedding(
            lattice,
            embedding,
            StaticLaw::bond(1.0).unwrap(),
            plan
        )
        .err()
        .expect("magnitude gate rejects the adversarial table"),
        StaticPercolationError::EmbeddingMagnitude { .. }
    ));
}
