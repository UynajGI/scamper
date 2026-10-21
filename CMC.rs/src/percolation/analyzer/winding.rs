//! Winding-capable component analysis over a periodic embedding (W2).
//!
//! One pass over physical edges runs a union-find with integer `Z^D`
//! potentials beside the ordinary connectivity bookkeeping: each vertex
//! stores its displacement toward its union-find parent, every union across
//! an active edge updates that potential with checked `i32` arithmetic, and
//! every same-root closure (a cycle) contributes the accumulated potential
//! mismatch as a winding generator of that component. Generators are folded
//! incrementally into a per-root echelon basis (at most `D` linearly
//! independent vectors), so the analysis reports, per component, the winding
//! rank and a generating set of the winding subgroup — never only a Boolean.
//!
//! The analysis is mode-free exactly like F3: it consumes the
//! [`UndirectedGraphView`] plus [`VertexActivity`]/[`EdgeActivity`]
//! capabilities and one [`PeriodicEmbedding`](crate::PeriodicEmbedding)
//! table. It is an independent implementation on its own
//! [`WindingWorkspace`] rather than a layer over the F3 analyzer: potentials
//! require their own path-compressed union-find, so composing on F3 results
//! would still re-scan every edge and materialize labels while changing no
//! work; the single fused scan here is the roadmap's shared-traversal design
//! orientation. F3's semantics and workspace are untouched.

use core::fmt;

use super::super::{EdgeActivity, VertexActivity, WindingWorkspace};
use crate::topology::{PeriodicEmbedding, VertexId};
use crate::UndirectedGraphView;

/// Aggregate winding facts for one active configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WindingSummary {
    /// Active vertices (identical to F3 `ComponentSummary::active_vertex_count`).
    pub active_vertex_count: usize,
    /// Active physical edges with both endpoints active (F3 parity).
    pub active_edge_count: usize,
    /// Active components including isolates.
    pub component_count: usize,
    /// Components with winding rank > 0 (wrapping components).
    pub wrapping_component_count: usize,
    /// Canonical identity (lowest active vertex ID) of the largest component,
    /// with the F3 tie policy (equal sizes resolved by lowest identity).
    pub largest_component: Option<VertexId>,
    /// Size of [`Self::largest_component`], or zero when nothing is active.
    pub largest_component_size: usize,
    /// Winding rank of the canonical largest component (`0..=D`).
    pub largest_winding_rank: usize,
}

impl WindingSummary {
    const EMPTY: Self = Self {
        active_vertex_count: 0,
        active_edge_count: 0,
        component_count: 0,
        wrapping_component_count: 0,
        largest_component: None,
        largest_component_size: 0,
        largest_winding_rank: 0,
    };
}

/// Winding facts of one component.
///
/// Identities use the F3 canonical policy: a component's identity is its
/// lowest active vertex ID. Components are listed in ascending canonical
/// identity order.
#[derive(Clone, Copy, Debug)]
pub struct ComponentWinding {
    component: VertexId,
    size: usize,
    rank: usize,
    root: usize,
}

impl ComponentWinding {
    /// Canonical component identity (lowest active vertex ID in it).
    #[inline]
    pub const fn component(self) -> VertexId {
        self.component
    }

    /// Number of active vertices in the component.
    #[inline]
    pub const fn size(self) -> usize {
        self.size
    }

    /// Rank of the winding subgroup: the number of independent periodic
    /// directions the component winds around. Zero means contractible
    /// (non-wrapping).
    #[inline]
    pub const fn rank(self) -> usize {
        self.rank
    }
}

/// Borrowed winding analysis result; the workspace stays mutably locked for
/// the lifetime of this result.
#[derive(Clone, Copy, Debug)]
pub struct WindingResult<'a> {
    summary: WindingSummary,
    components: &'a [ComponentWinding],
    basis: &'a [i32],
    dimensions: usize,
}

impl WindingResult<'_> {
    #[inline]
    pub const fn summary(&self) -> WindingSummary {
        self.summary
    }

    /// All active components in ascending canonical-identity order.
    #[inline]
    pub const fn components(&self) -> &[ComponentWinding] {
        self.components
    }

    /// Embedding dimension `D` of this analysis.
    #[inline]
    pub const fn dimensions(&self) -> usize {
        self.dimensions
    }

    /// Linearly independent winding generators of component `component` (an
    /// index into [`Self::components`]): exactly `rank` vectors of `D` signed
    /// cell components each, in ascending first-nonzero-axis order with
    /// positive leading components. Empty when the component does not wrap.
    ///
    /// The vectors generate the full winding subgroup lattice; they are an
    /// echelon generating set, not a canonical Hermite normal form.
    pub fn generators(&self, component: usize) -> impl Iterator<Item = &[i32]> + '_ {
        let info = &self.components[component];
        let block = self.dimensions * self.dimensions;
        self.basis[info.root * block..][..info.rank * self.dimensions].chunks_exact(self.dimensions)
    }

    /// Whether component `component` (an index into [`Self::components`])
    /// wraps: winding rank greater than zero.
    pub fn wraps(&self, component: usize) -> bool {
        self.components[component].rank > 0
    }

    /// Index into [`Self::components`] of the component with canonical
    /// identity `component`, if present.
    pub fn component_index(&self, component: VertexId) -> Option<usize> {
        self.components
            .iter()
            .position(|info| info.component == component)
    }
}

/// Analyze winding and connectivity of one active configuration in a single
/// physical-edge scan.
///
/// `graph`, `vertex_activity`, and `edge_activity` follow the F3 `analyze`
/// contract (dense IDs, activity domains matching the graph, active edge =
/// edge active with both endpoints active). `embedding` is the W1 periodic
/// displacement table for the *same* graph; displacements are looked up by
/// the directed incidence `(edge, endpoint)` and never inferred.
///
/// All `i32` potential arithmetic is checked: an accumulation that would
/// silently wrap instead returns [`WindingError::DisplacementOverflow`],
/// attributed to the physical edge being scanned. The check is conservative
/// — it can reject an adversarial displacement table whose intermediate sums
/// exceed `i32` even when the mathematical winding subgroup is small — and
/// never wraps silently.
///
/// The result assumes the W1 producer contract (displacements form a
/// consistent cocycle). For an arbitrary antisymmetric table the analysis is
/// still well-defined: it reports the subgroup generated by cycle sums of
/// the supplied assignment.
pub fn analyze_winding<'a, G, V, E, P>(
    graph: &G,
    embedding: &P,
    vertex_activity: &V,
    edge_activity: &E,
    workspace: &'a mut WindingWorkspace,
) -> Result<WindingResult<'a>, WindingError>
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
    P: PeriodicEmbedding,
{
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();
    validate_winding_domains(
        vertex_count,
        edge_count,
        vertex_activity,
        edge_activity,
        embedding,
    )?;
    let dimensions = embedding.periodic_dimensions();
    workspace
        .resize_for(vertex_count, dimensions)
        .map_err(WindingError::Capacity)?;

    let WindingWorkspace {
        parent,
        size,
        canonical,
        potential,
        basis,
        basis_counts,
        root_slot,
        scratch,
        chain,
        components,
        dimensions: scratch_dimensions,
    } = &mut *workspace;
    debug_assert_eq!(*scratch_dimensions, dimensions);
    let mut buffers = WindingBuffers {
        parent,
        size,
        canonical,
        potential,
        basis,
        basis_counts,
        root_slot,
        scratch,
        chain,
        components,
        dimensions,
    };

    let active_vertex_count = initialize_vertices(graph, vertex_activity, &mut buffers);
    let active_edge_count = scan_edges(graph, embedding, edge_activity, &mut buffers)?;
    let summary = collect_components(graph, &mut buffers, active_vertex_count, active_edge_count);

    let WindingBuffers {
        components,
        basis,
        dimensions,
        ..
    } = buffers;
    let components: &[ComponentWinding] = components;
    let basis: &[i32] = basis;
    Ok(WindingResult {
        summary,
        components,
        basis,
        dimensions,
    })
}

/// Reusable per-analysis buffers handed to each stage of the single scan.
struct WindingBuffers<'a> {
    parent: &'a mut [usize],
    size: &'a mut [usize],
    canonical: &'a mut [usize],
    potential: &'a mut [i32],
    basis: &'a mut [i32],
    basis_counts: &'a mut [u32],
    root_slot: &'a mut [u32],
    scratch: &'a mut [i32],
    chain: &'a mut Vec<usize>,
    components: &'a mut Vec<ComponentWinding>,
    dimensions: usize,
}

fn validate_winding_domains<V, E, P>(
    vertex_count: usize,
    edge_count: usize,
    vertex_activity: &V,
    edge_activity: &E,
    embedding: &P,
) -> Result<(), WindingError>
where
    V: VertexActivity,
    E: EdgeActivity,
    P: PeriodicEmbedding,
{
    if vertex_activity.vertex_count() != vertex_count {
        return Err(WindingError::VertexActivityLength {
            expected: vertex_count,
            actual: vertex_activity.vertex_count(),
        });
    }
    if edge_activity.edge_count() != edge_count {
        return Err(WindingError::EdgeActivityLength {
            expected: edge_count,
            actual: edge_activity.edge_count(),
        });
    }
    if embedding.periodic_dimensions() == 0 {
        return Err(WindingError::EmbeddingDimensionZero);
    }
    Ok(())
}

fn initialize_vertices<G, V>(
    graph: &G,
    vertex_activity: &V,
    buffers: &mut WindingBuffers<'_>,
) -> usize
where
    G: UndirectedGraphView,
    V: VertexActivity,
{
    let mut active = 0;
    for index in 0..graph.vertex_count() {
        let vertex = graph.vertex_id(index).expect("dense vertex index");
        // Basis counts are per-root state from the previous analysis on a
        // reused workspace; every entry resets here, roots included.
        buffers.basis_counts[index] = 0;
        if vertex_activity.vertex_active(vertex) {
            buffers.parent[index] = index;
            buffers.size[index] = 1;
            buffers.canonical[index] = index;
            active += 1;
        } else {
            buffers.parent[index] = usize::MAX;
            buffers.size[index] = 0;
            buffers.canonical[index] = usize::MAX;
        }
    }
    active
}

fn scan_edges<G, E, P>(
    graph: &G,
    embedding: &P,
    edge_activity: &E,
    buffers: &mut WindingBuffers<'_>,
) -> Result<usize, WindingError>
where
    G: UndirectedGraphView,
    E: EdgeActivity,
    P: PeriodicEmbedding,
{
    let mut active = 0;
    for edge_index in 0..graph.edge_count() {
        active += process_active_edge(graph, embedding, edge_activity, edge_index, buffers)?;
    }
    Ok(active)
}

fn process_active_edge<G, E, P>(
    graph: &G,
    embedding: &P,
    edge_activity: &E,
    edge_index: usize,
    buffers: &mut WindingBuffers<'_>,
) -> Result<usize, WindingError>
where
    G: UndirectedGraphView,
    E: EdgeActivity,
    P: PeriodicEmbedding,
{
    let edge = graph
        .edge_id(edge_index)
        .expect("dense physical edge index");
    if !edge_activity.edge_active(edge) {
        return Ok(0);
    }
    let [left, right] = graph.edge_endpoints(edge);
    let (left_index, right_index) = (left.index(), right.index());
    if buffers.parent[left_index] == usize::MAX || buffers.parent[right_index] == usize::MAX {
        return Ok(0);
    }
    let displacement = embedding
        .edge_displacement(edge, left)
        .ok_or(WindingError::MissingDisplacement { edge: edge_index })?;

    let WindingBuffers {
        parent,
        size,
        canonical,
        potential,
        basis,
        basis_counts,
        scratch,
        chain,
        dimensions,
        ..
    } = buffers;
    let dimensions = *dimensions;
    let (left_potential, scratch) = scratch.split_at_mut(dimensions);
    let (right_potential, scratch) = scratch.split_at_mut(dimensions);
    let (vector, work) = scratch.split_at_mut(dimensions);

    let root_left = find_root_with_potential(
        parent,
        potential,
        chain,
        dimensions,
        left_index,
        edge_index,
        left_potential,
    )?;
    let root_right = find_root_with_potential(
        parent,
        potential,
        chain,
        dimensions,
        right_index,
        edge_index,
        right_potential,
    )?;

    // `vector` = displacement(root_left -> root_right) through this edge:
    // down to `left`, across it (`displacement`), up from `right`.
    for axis in 0..dimensions {
        let component = displacement
            .axis(axis)
            .expect("embedding reports its own dimension");
        vector[axis] = component
            .checked_add(right_potential[axis])
            .ok_or(WindingError::DisplacementOverflow { edge: edge_index })?;
    }

    if root_left == root_right {
        // Same-root closure: the mismatch around the cycle is a winding
        // generator observation for this component.
        record_cycle(
            basis,
            basis_counts,
            vector,
            left_potential,
            work,
            root_left,
            edge_index,
        )?;
    } else if size[root_left] >= size[root_right] {
        // Attach the right tree under the left root; the new potential of
        // the old right root is -(root_left -> root_right).
        set_root_potential(
            potential,
            dimensions,
            root_right,
            left_potential,
            vector,
            true,
            edge_index,
        )?;
        parent[root_right] = root_left;
        size[root_left] += size[root_right];
        canonical[root_left] = canonical[root_left].min(canonical[root_right]);
        merge_basis_blocks(
            basis,
            basis_counts,
            dimensions,
            root_left,
            root_right,
            work,
            edge_index,
        )?;
    } else {
        set_root_potential(
            potential,
            dimensions,
            root_left,
            left_potential,
            vector,
            false,
            edge_index,
        )?;
        parent[root_left] = root_right;
        size[root_right] += size[root_left];
        canonical[root_right] = canonical[root_right].min(canonical[root_left]);
        merge_basis_blocks(
            basis,
            basis_counts,
            dimensions,
            root_right,
            root_left,
            work,
            edge_index,
        )?;
    }
    Ok(1)
}

fn record_cycle(
    basis: &mut [i32],
    basis_counts: &mut [u32],
    vector: &mut [i32],
    left_potential: &[i32],
    work: &mut [i32],
    root: usize,
    edge: usize,
) -> Result<(), WindingError> {
    let dimensions = vector.len();
    let mut nonzero = false;
    for axis in 0..dimensions {
        let mismatch = vector[axis]
            .checked_sub(left_potential[axis])
            .ok_or(WindingError::DisplacementOverflow { edge })?;
        nonzero |= mismatch != 0;
        vector[axis] = mismatch;
    }
    if !nonzero {
        return Ok(());
    }
    let block = dimensions * dimensions;
    let count = &mut basis_counts[root];
    basis_insert(
        &mut basis[root * block..][..block],
        count,
        vector,
        work,
        edge,
    )
}

fn collect_components<G>(
    graph: &G,
    buffers: &mut WindingBuffers<'_>,
    active_vertex_count: usize,
    active_edge_count: usize,
) -> WindingSummary
where
    G: UndirectedGraphView,
{
    let WindingBuffers {
        parent,
        size,
        canonical,
        basis_counts,
        root_slot,
        components,
        ..
    } = buffers;
    // Component listing in ascending canonical order: scanning vertices
    // ascending, the first vertex of each root is that component's minimum,
    // so no sort is needed. Strictly-greater size updates keep the F3
    // canonical-largest tie policy (lowest identity among equal sizes).
    components.clear();
    root_slot.fill(u32::MAX);
    let mut summary = WindingSummary {
        active_vertex_count,
        active_edge_count,
        ..WindingSummary::EMPTY
    };
    let mut largest: Option<(usize, usize, usize)> = None;
    for index in 0..graph.vertex_count() {
        if parent[index] == usize::MAX {
            continue;
        }
        let root = find_root(parent, index);
        if root_slot[root] != u32::MAX {
            continue;
        }
        root_slot[root] = components.len() as u32;
        let identity = canonical[root];
        let component_size = size[root];
        let rank = basis_counts[root] as usize;
        if rank > 0 {
            summary.wrapping_component_count += 1;
        }
        if largest.is_none_or(|(_, largest_size, _)| component_size > largest_size) {
            largest = Some((identity, component_size, rank));
        }
        components.push(ComponentWinding {
            component: graph.vertex_id(identity).expect("dense vertex index"),
            size: component_size,
            rank,
            root,
        });
    }
    summary.component_count = components.len();
    if let Some((identity, component_size, rank)) = largest {
        summary.largest_component = graph.vertex_id(identity);
        summary.largest_component_size = component_size;
        summary.largest_winding_rank = rank;
    }
    summary
}

/// Writes the potential of a root that is being attached under another root.
///
/// With `negate`, the stored potential is `left_potential - vector`
/// (attach right root under left: displacement right-root -> left-root);
/// otherwise `vector - left_potential` (attach left under right).
fn set_root_potential(
    potential: &mut [i32],
    dimensions: usize,
    root: usize,
    left_potential: &[i32],
    vector: &[i32],
    negate: bool,
    edge: usize,
) -> Result<(), WindingError> {
    for axis in 0..dimensions {
        let value = if negate {
            left_potential[axis].checked_sub(vector[axis])
        } else {
            vector[axis].checked_sub(left_potential[axis])
        }
        .ok_or(WindingError::DisplacementOverflow { edge })?;
        potential[root * dimensions + axis] = value;
    }
    Ok(())
}

/// Path-compressed find that also returns the node-to-root potential and
/// rewrites every visited potential relative to the root.
///
/// Potentials are recomputed as suffix sums walking from the root, so no
/// subtraction of partial sums is ever needed. Returns
/// [`WindingError::DisplacementOverflow`] (attributed to `edge`) when a
/// suffix sum leaves `i32`.
fn find_root_with_potential(
    parent: &mut [usize],
    potential: &mut [i32],
    chain: &mut Vec<usize>,
    dimensions: usize,
    node: usize,
    edge: usize,
    out: &mut [i32],
) -> Result<usize, WindingError> {
    chain.clear();
    let mut root = node;
    chain.push(node);
    while parent[root] != root {
        root = parent[root];
        chain.push(root);
    }
    out.iter_mut().for_each(|value| *value = 0);
    for index in (0..chain.len() - 1).rev() {
        let vertex = chain[index];
        for axis in 0..dimensions {
            let step = potential[vertex * dimensions + axis];
            out[axis] = step
                .checked_add(out[axis])
                .ok_or(WindingError::DisplacementOverflow { edge })?;
        }
        parent[vertex] = root;
        potential[vertex * dimensions..][..dimensions].copy_from_slice(out);
    }
    Ok(root)
}

/// Plain path-compressed find without potentials (final component pass).
fn find_root(parent: &mut [usize], node: usize) -> usize {
    let mut root = node;
    while parent[root] != root {
        root = parent[root];
    }
    let mut current = node;
    while parent[current] != root {
        let next = parent[current];
        parent[current] = root;
        current = next;
    }
    root
}

/// Writes `vector` into one per-root basis `block` (length `D * D`) keeping
/// the block a first-nonzero-pivot echelon generating set of the same
/// lattice: rows sorted by pivot axis, distinct pivots, positive leading
/// components, at most `D` rows.
///
/// Each outer iteration combines `vector` against the row sharing its
/// first-nonzero axis so that the vector's leading component becomes exactly
/// zero; its first-nonzero axis therefore strictly increases, terminating
/// within `D + 1` iterations with either an appended row (rank grew) or the
/// zero vector (already generated).
fn basis_insert(
    block: &mut [i32],
    count: &mut u32,
    vector: &[i32],
    work: &mut [i32],
    edge: usize,
) -> Result<(), WindingError> {
    let dimensions = vector.len();
    debug_assert_eq!(block.len(), dimensions * dimensions);
    debug_assert_eq!(work.len(), dimensions);
    work.copy_from_slice(vector);
    loop {
        let Some(pivot) = first_nonzero(work) else {
            return Ok(());
        };
        let mut target = None;
        let mut insertion = 0;
        for row in 0..*count as usize {
            let row_pivot = first_nonzero(&block[row * dimensions..][..dimensions])
                .expect("stored basis rows are nonzero");
            if row_pivot == pivot {
                target = Some(row);
                break;
            }
            if row_pivot < pivot {
                insertion += 1;
            }
        }
        let Some(row) = target else {
            if work[pivot] < 0 {
                for component in work.iter_mut() {
                    *component = 0i32
                        .checked_sub(*component)
                        .ok_or(WindingError::DisplacementOverflow { edge })?;
                }
            }
            block.copy_within(
                insertion * dimensions..*count as usize * dimensions,
                (insertion + 1) * dimensions,
            );
            block[insertion * dimensions..][..dimensions].copy_from_slice(work);
            *count += 1;
            return Ok(());
        };
        let (leading, incoming) = (block[row * dimensions + pivot] as i64, work[pivot] as i64);
        let (g, x, y) = extended_gcd(leading, incoming);
        for axis in 0..dimensions {
            let row_component = block[row * dimensions + axis] as i64;
            let work_component = work[axis] as i64;
            let new_row = x
                .checked_mul(row_component)
                .and_then(|value| value.checked_add(y.checked_mul(work_component)?))
                .ok_or(WindingError::DisplacementOverflow { edge })?;
            let new_work = (leading / g)
                .checked_mul(work_component)
                .and_then(|value| value.checked_sub((incoming / g).checked_mul(row_component)?))
                .ok_or(WindingError::DisplacementOverflow { edge })?;
            block[row * dimensions + axis] =
                i32::try_from(new_row).map_err(|_| WindingError::DisplacementOverflow { edge })?;
            work[axis] =
                i32::try_from(new_work).map_err(|_| WindingError::DisplacementOverflow { edge })?;
        }
    }
}

/// Folds the absorbed root's basis rows into the surviving root's block.
fn merge_basis_blocks(
    basis: &mut [i32],
    basis_counts: &mut [u32],
    dimensions: usize,
    keep: usize,
    drop: usize,
    work: &mut [i32],
    edge: usize,
) -> Result<(), WindingError> {
    debug_assert_ne!(keep, drop);
    let block = dimensions * dimensions;
    let dropped_rows = basis_counts[drop];
    basis_counts[drop] = 0;
    let (keep_block, drop_block) = if keep < drop {
        let (low, high) = basis.split_at_mut(drop * block);
        (&mut low[keep * block..][..block], &mut high[..block])
    } else {
        let (low, high) = basis.split_at_mut(keep * block);
        (&mut high[..block], &mut low[drop * block..][..block])
    };
    let count = &mut basis_counts[keep];
    for row in 0..dropped_rows as usize {
        basis_insert(
            keep_block,
            count,
            &drop_block[row * dimensions..][..dimensions],
            work,
            edge,
        )?;
    }
    Ok(())
}

fn first_nonzero(vector: &[i32]) -> Option<usize> {
    vector.iter().position(|&component| component != 0)
}

/// Extended Euclid on `i64`: returns `(g, x, y)` with `a*x + b*y = g` and
/// `g = gcd(|a|, |b|) > 0` for not-both-zero inputs.
fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1i64, 0i64);
    let (mut old_t, mut t) = (0i64, 1i64);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }
    if old_r < 0 {
        (-old_r, -old_s, -old_t)
    } else {
        (old_r, old_s, old_t)
    }
}

/// Validation failure of a winding analysis.
#[derive(Debug)]
pub enum WindingError {
    VertexActivityLength {
        expected: usize,
        actual: usize,
    },
    EdgeActivityLength {
        expected: usize,
        actual: usize,
    },
    /// The embedding reported zero periodic dimensions, violating the
    /// `PeriodicEmbedding` contract (`D >= 1`).
    EmbeddingDimensionZero,
    /// The embedding had no displacement for one endpoint of an active edge:
    /// the table does not describe this graph.
    MissingDisplacement {
        edge: usize,
    },
    /// A checked `i32` potential or basis accumulation overflowed while
    /// scanning `edge`. Conservative: fires on intermediate sums leaving
    /// `i32`, never wraps silently.
    DisplacementOverflow {
        edge: usize,
    },
    Capacity(std::collections::TryReserveError),
}

impl fmt::Display for WindingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::VertexActivityLength { expected, actual } => write!(
                formatter,
                "vertex activity length is {actual}, expected {expected}"
            ),
            Self::EdgeActivityLength { expected, actual } => write!(
                formatter,
                "edge activity length is {actual}, expected {expected}"
            ),
            Self::EmbeddingDimensionZero => {
                write!(
                    formatter,
                    "periodic embedding dimension must be at least one"
                )
            }
            Self::MissingDisplacement { edge } => write!(
                formatter,
                "embedding has no displacement for endpoint of active edge {edge}"
            ),
            Self::DisplacementOverflow { edge } => write!(
                formatter,
                "i32 displacement accumulation overflowed while scanning edge {edge}"
            ),
            Self::Capacity(_) => formatter.write_str("could not reserve winding workspace"),
        }
    }
}

impl std::error::Error for WindingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        build_chain_with_embedding, build_square_with_embedding, Bond, BondType, CsrLattice,
        EdgeId, GraphView, LatticeEmbedding, StaticConfiguration,
    };

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

    /// Embedding wrapper whose displacement lookup fails for one edge.
    struct EmbeddingMissingEdge<'a> {
        inner: &'a LatticeEmbedding,
        missing: usize,
    }

    impl PeriodicEmbedding for EmbeddingMissingEdge<'_> {
        fn periodic_dimensions(&self) -> usize {
            self.inner.periodic_dimensions()
        }
        fn edge_displacement(
            &self,
            edge: EdgeId,
            from: VertexId,
        ) -> Option<crate::DirectedDisplacement<'_>> {
            if edge.index() == self.missing {
                None
            } else {
                self.inner.edge_displacement(edge, from)
            }
        }
    }

    struct ZeroDimensionEmbedding;

    impl PeriodicEmbedding for ZeroDimensionEmbedding {
        fn periodic_dimensions(&self) -> usize {
            0
        }
        fn edge_displacement(
            &self,
            _edge: EdgeId,
            _from: VertexId,
        ) -> Option<crate::DirectedDisplacement<'_>> {
            None
        }
    }

    #[test]
    fn ring_winds_exactly_once_with_unit_generator() {
        let (lattice, embedding) = build_chain_with_embedding(4, true).expect("ring embedding");
        let (vertices, edges) = all_active(&lattice);
        let mut workspace = WindingWorkspace::new();
        let result =
            analyze_winding(&lattice, &embedding, &vertices, &edges, &mut workspace).unwrap();
        assert_eq!(result.dimensions(), 1);
        assert_eq!(result.summary().component_count, 1);
        assert_eq!(result.summary().wrapping_component_count, 1);
        let info = result.components()[0];
        assert_eq!(
            (info.component().index(), info.size(), info.rank()),
            (0, 4, 1)
        );
        assert_eq!(
            result.generators(0).collect::<Vec<_>>(),
            [&[1][..]],
            "single-axis wrap generator is the unit vector"
        );
    }

    #[test]
    fn activity_domain_mismatches_are_typed_rejections() {
        let (lattice, embedding) = build_chain_with_embedding(3, true).expect("ring embedding");
        let (vertices, edges) = all_active(&lattice);
        let mut workspace = WindingWorkspace::new();
        let short_vertices = VertexFlags(vec![true; 2]);
        assert!(matches!(
            analyze_winding(
                &lattice,
                &embedding,
                &short_vertices,
                &edges,
                &mut workspace
            )
            .unwrap_err(),
            WindingError::VertexActivityLength {
                expected: 3,
                actual: 2
            }
        ));
        let short_edges = EdgeFlags(vec![true; 1]);
        assert!(matches!(
            analyze_winding(
                &lattice,
                &embedding,
                &vertices,
                &short_edges,
                &mut workspace
            )
            .unwrap_err(),
            WindingError::EdgeActivityLength {
                expected: 3,
                actual: 1
            }
        ));
    }

    #[test]
    fn broken_embedding_capabilities_are_typed_rejections() {
        let (lattice, embedding) = build_chain_with_embedding(3, true).expect("ring embedding");
        let (vertices, edges) = all_active(&lattice);
        let mut workspace = WindingWorkspace::new();
        assert!(matches!(
            analyze_winding(
                &lattice,
                &ZeroDimensionEmbedding,
                &vertices,
                &edges,
                &mut workspace
            )
            .unwrap_err(),
            WindingError::EmbeddingDimensionZero
        ));
        let broken = EmbeddingMissingEdge {
            inner: &embedding,
            missing: 1,
        };
        assert!(matches!(
            analyze_winding(&lattice, &broken, &vertices, &edges, &mut workspace).unwrap_err(),
            WindingError::MissingDisplacement { edge: 1 }
        ));
    }

    #[test]
    fn adversarial_displacements_overflow_with_typed_error_and_small_magnitudes_do_not() {
        // 3-cycle with D=1. Cycle sums of 2 * 10^9 fit i32 exactly
        // (contractible); 3 * 10^9 do not and must be a typed overflow.
        let lattice = CsrLattice::from_edges(
            3,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
                Bond::new(2, 0, BondType::Generic, 1.0),
            ],
        );
        let (vertices, edges) = all_active(&lattice);
        let mut workspace = WindingWorkspace::new();

        let contractible = LatticeEmbedding::try_from_edge_displacements(
            &lattice,
            1,
            &[1_000_000_000, 1_000_000_000, -2_000_000_000],
        )
        .expect("valid table");
        let result =
            analyze_winding(&lattice, &contractible, &vertices, &edges, &mut workspace).unwrap();
        assert_eq!(
            result.summary().wrapping_component_count,
            0,
            "zero cycle sum"
        );

        let overflowing = LatticeEmbedding::try_from_edge_displacements(
            &lattice,
            1,
            &[1_000_000_000, 1_000_000_000, 1_000_000_000],
        )
        .expect("valid table");
        assert!(matches!(
            analyze_winding(&lattice, &overflowing, &vertices, &edges, &mut workspace).unwrap_err(),
            WindingError::DisplacementOverflow { edge: 2 }
        ));
    }

    #[test]
    fn workspace_reuse_and_dimension_switch_stay_consistent() {
        let (ring, ring_embedding) = build_chain_with_embedding(5, true).expect("ring embedding");
        let (square, square_embedding) =
            build_square_with_embedding(3, 3, true).expect("torus embedding");
        let mut workspace = WindingWorkspace::new();
        for _ in 0..3 {
            let (vertices, edges) = all_active(&ring);
            let result =
                analyze_winding(&ring, &ring_embedding, &vertices, &edges, &mut workspace).unwrap();
            assert_eq!(result.summary().largest_winding_rank, 1);
        }
        for _ in 0..3 {
            let (vertices, edges) = all_active(&square);
            let result = analyze_winding(
                &square,
                &square_embedding,
                &vertices,
                &edges,
                &mut workspace,
            )
            .unwrap();
            assert_eq!(result.summary().largest_winding_rank, 2);
            assert_eq!(
                result.generators(0).collect::<Vec<_>>(),
                [&[1, 0][..], &[0, 1][..]],
                "full torus winds along both unit axes"
            );
        }
    }

    #[test]
    fn empty_configuration_has_empty_winding_summary() {
        let (lattice, embedding) = build_square_with_embedding(2, 2, true).expect("torus");
        let vertices = VertexFlags(vec![false; 4]);
        let edges = EdgeFlags(vec![false; 8]);
        let mut workspace = WindingWorkspace::new();
        let result =
            analyze_winding(&lattice, &embedding, &vertices, &edges, &mut workspace).unwrap();
        assert_eq!(result.summary(), WindingSummary::EMPTY);
        assert!(result.components().is_empty());
        assert_eq!(result.component_index(lattice.vertex_id(0).unwrap()), None);
    }

    #[test]
    fn static_configuration_activity_composes_with_winding_analysis() {
        let (lattice, embedding) = build_square_with_embedding(3, 3, true).expect("torus");
        let configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
        let mut workspace = WindingWorkspace::new();
        let result = analyze_winding(
            &lattice,
            &embedding,
            &configuration,
            &configuration,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(result.summary(), WindingSummary::EMPTY);
    }

    #[test]
    fn extended_gcd_reference_values() {
        assert_eq!(extended_gcd(2, 3), (1, -1, 1));
        assert_eq!(extended_gcd(12, 8), (4, 1, -1));
        assert_eq!(extended_gcd(-4, 6), (2, 1, 1));
        assert_eq!(extended_gcd(5, 0), (5, 1, 0));
        let (g, x, y) = extended_gcd(-9, -6);
        assert_eq!(g, 3);
        assert_eq!(-9 * x + -6 * y, 3);
    }
}
