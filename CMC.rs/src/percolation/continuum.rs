//! Continuum candidate-pair adapter (P track, roadmap §14 P1).
//!
//! This is the geometry/candidate-pair infrastructure that Boolean
//! percolation (P2) will consume — it deliberately implements **no** overlap
//! or percolation observables of its own. The design rule is
//! borrow-not-duplicate: coordinates, the periodic geometry, and the
//! `CellList` stencil search are the particle module's
//! ([`crate::particle`]) and are reused verbatim; nothing here forks a
//! second geometry implementation, and the particle module is not modified.
//!
//! # Geometry judgment
//!
//! The adapter takes a concrete [`OrthorhombicCell<D>`], not a generic
//! `SimulationCell<D>`: the reused `CellList` bucketing
//! (`cell_of_position`/`rebuild`) is defined only over orthorhombic side
//! lengths, so accepting the wider trait would let callers supply a
//! triclinic geometry the stencil search cannot support. Minimum-image
//! geometry flows through the particle module's own metric
//! (`OrthorhombicCell::distance_squared` over the
//! [`crate::particle::SimulationCell`] displacement convention)
//! trait methods on that concrete type. Supporting other cell shapes means
//! generalizing the particle module's cell list first — out of P1/P2 scope
//! (the roadmap's initial continuum scope is orthorhombic cells).
//!
//! # Ephemeral `PairGraph`, not a persistent CSR
//!
//! The roadmap forbids materializing a *persistent* CSR for every particle
//! snapshot: an equilibrium MC simulation would rebuild it every sweep,
//! paying `O(E)` construction churn on top of each analysis.
//! [`PairGraph`] is the deliberate opposite — a lightweight, **transient**
//! snapshot built per analysis call from the cell list (a flat
//! `vertex_count` plus a `Vec` of pairs), handed to the F3 analyzer, and
//! dropped. It never survives the analysis that requested it, and no
//! incidence/adjacency index is built for it: `incidences` are derived by
//! scanning the pair list on demand (see the method docs).
//!
//! # Pair order contract
//!
//! Pairs `(i, j)` with `i < j` are emitted in a canonical deterministic
//! order:
//!
//! 1. **Cell-major:** the outer loop visits packed cells in ascending
//!    linear index (row-major over the `CellList` grid, axis 0 fastest —
//!    its strides make axis 0 the least significant; review follow-up
//!    2026-10-23).
//! 2. **Within-cell ascending:** for each cell, particles whose *own* cell
//!    it is are visited in ascending particle index `i` (the first
//!    endpoint), and each accepted partner `j > i` is emitted in ascending
//!    index order.
//!
//! The sequence is therefore sorted by `(cell(i), i, j)` — it is **not**
//! globally sorted by `i`, because particle indices are not cell-major.
//! Consequences: an exact lattice-vector translation (bit-identical wrapped
//! coordinates) reproduces the identical sequence; a generic translation
//! preserves the pair *set* but may reorder it when particles cross cell
//! boundaries. A pair is emitted exactly once, from the first endpoint's
//! own cell, with `j` guaranteed to lie in that cell's `3^D` stencil — the
//! stencil provably contains every partner within the cutoff because the
//! cell width is at least the cutoff.
//!
//! # Distances
//!
//! A pair is a candidate when its minimum-image squared distance is
//! **`<=` the squared cutoff** (boundary-inclusive). The cutoff must be
//! finite, positive, and at most half the shortest cell side (at exactly
//! half, two images of a pair can both be within range; the
//! minimum-image tie picks one via round-half-away-from-zero and the
//! unordered pair set stays well-defined — the pair qualifies via either
//! image)
//! (minimum-image limit, enforced by the reused `CellList`).
//!
//! Zero particles are rejected with a typed error rather than analyzed:
//! a candidate-pair analysis of an empty coordinate buffer is almost
//! certainly a caller bug, and the rejection makes that loud.

use crate::particle::{CellList, OrthorhombicCell, ParticleConfiguration, ParticleError};
use crate::topology::{EdgeId, GraphView, Incidence, UndirectedGraphView, VertexId};

/// Typed rejection of a continuum candidate-pair request.
#[derive(Debug, Clone, PartialEq)]
pub enum ContinuumError {
    /// The const dimension is zero; no geometry exists.
    ZeroDimension,
    /// The cutoff was non-finite or not strictly positive.
    InvalidCutoff { cutoff: f64 },
    /// The coordinate buffer holds no particles.
    EmptyCoordinates { dimension: usize },
    /// The flat coordinate length is not a multiple of the dimension.
    CoordinateCountMismatch { length: usize, dimension: usize },
    /// One coordinate was non-finite.
    NonFiniteCoordinate { particle: usize, axis: usize },
    /// The cutoff exceeds the minimum-image maximum of the cell (half the
    /// shortest side).
    CutoffTooLarge { cutoff: f64, maximum: f64 },
    /// The reused cell list could not be built or did not validate against
    /// the configuration it was claimed to describe.
    CellListInconsistent(String),
}

impl core::fmt::Display for ContinuumError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ZeroDimension => write!(formatter, "particle dimension must be positive"),
            Self::InvalidCutoff { cutoff } => write!(
                formatter,
                "candidate-pair cutoff must be finite and positive, got {cutoff}"
            ),
            Self::EmptyCoordinates { dimension } => write!(
                formatter,
                "coordinate buffer holds no {dimension}-dimensional particles"
            ),
            Self::CoordinateCountMismatch { length, dimension } => write!(
                formatter,
                "flat coordinate length {length} is not a multiple of dimension {dimension}"
            ),
            Self::NonFiniteCoordinate { particle, axis } => write!(
                formatter,
                "particle {particle} coordinate {axis} is non-finite"
            ),
            Self::CutoffTooLarge { cutoff, maximum } => write!(
                formatter,
                "cutoff {cutoff} exceeds minimum-image maximum {maximum}"
            ),
            Self::CellListInconsistent(reason) => {
                write!(
                    formatter,
                    "cell list inconsistent with the request: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for ContinuumError {}

/// Borrowed per-analysis adapter from flat particle coordinates to candidate
/// pairs (P1).
///
/// Holds only borrows — the flat `&[f64]` coordinates, the periodic cell,
/// and the cutoff — so constructing it is free. Each [`Self::pairs`] call
/// performs one transient analysis: it copies the coordinates into a
/// throwaway [`ParticleConfiguration`], builds a [`CellList`] over it,
/// enumerates candidate pairs in the documented canonical order, and returns
/// an ephemeral owned [`PairGraph`]. Nothing persists between calls.
///
/// Callers already holding the particle module's
/// [`ParticleConfiguration`] (for example an equilibrium simulation
/// measuring connectivity of its accepted coordinates) should prefer the
/// zero-copy [`candidate_pairs_from_configuration`] or the prepared-list
/// [`candidate_pairs_with_cell_list`] instead of re-flattening positions.
#[derive(Clone, Copy, Debug)]
pub struct CandidatePairGraph<'a, const D: usize> {
    coordinates: &'a [f64],
    cell: &'a OrthorhombicCell<D>,
    cutoff: f64,
}

impl<'a, const D: usize> CandidatePairGraph<'a, D> {
    /// Validate borrowed flat coordinates, an orthorhombic periodic cell,
    /// and a cutoff radius.
    ///
    /// `coordinates` is laid out particle-major (`[x_0, y_0, x_1, y_1, ...]`
    /// in 2D); its length must be a positive multiple of `D` and every
    /// entry finite. Positions outside the primary cell are wrapped by the
    /// reused particle machinery during enumeration, matching the minimum
    /// image convention.
    pub fn new(
        coordinates: &'a [f64],
        cell: &'a OrthorhombicCell<D>,
        cutoff: f64,
    ) -> Result<Self, ContinuumError> {
        if D == 0 {
            return Err(ContinuumError::ZeroDimension);
        }
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(ContinuumError::InvalidCutoff { cutoff });
        }
        if !coordinates.len().is_multiple_of(D) {
            return Err(ContinuumError::CoordinateCountMismatch {
                length: coordinates.len(),
                dimension: D,
            });
        }
        if coordinates.is_empty() {
            return Err(ContinuumError::EmptyCoordinates { dimension: D });
        }
        for (flat, &value) in coordinates.iter().enumerate() {
            if !value.is_finite() {
                return Err(ContinuumError::NonFiniteCoordinate {
                    particle: flat / D,
                    axis: flat % D,
                });
            }
        }
        Ok(Self {
            coordinates,
            cell,
            cutoff,
        })
    }

    /// Cutoff radius this adapter was validated for.
    #[inline]
    pub const fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Number of `D`-dimensional particles in the borrowed buffer.
    #[inline]
    pub const fn particle_count(&self) -> usize {
        self.coordinates.len() / D
    }

    /// Enumerate candidate pairs into an ephemeral [`PairGraph`].
    ///
    /// One transient `O(N)` coordinate copy plus the cell-list build and
    /// stencil walk; see the [module docs](self) for the exact pair order.
    /// Deterministic: the same borrowed inputs yield the identical pair
    /// sequence on every call.
    pub fn pairs(&self) -> Result<PairGraph, ContinuumError> {
        let configuration = configuration_from_flat(self.coordinates, self.cell)?;
        build_from_configuration(&configuration, self.cutoff)
    }
}

/// Enumerate candidate pairs directly from the particle module's accepted
/// configuration (zero-copy path).
///
/// Equivalent to [`CandidatePairGraph::pairs`] on the flattened
/// `configuration.positions()`, without the transient coordinate copy.
pub fn candidate_pairs_from_configuration<const D: usize>(
    configuration: &ParticleConfiguration<D>,
    cutoff: f64,
) -> Result<PairGraph, ContinuumError> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(ContinuumError::InvalidCutoff { cutoff });
    }
    build_from_configuration(configuration, cutoff)
}

/// Enumerate candidate pairs reusing a caller-prepared [`CellList`].
///
/// The list must describe exactly this configuration; that contract is
/// verified with the particle module's own [`CellList::validate`] before
/// any geometry is computed (a typed
/// [`ContinuumError::CellListInconsistent`] otherwise). This skips the
/// per-call cell-list rebuild when a caller already maintains one at the
/// same cutoff.
pub fn candidate_pairs_with_cell_list<const D: usize>(
    configuration: &ParticleConfiguration<D>,
    list: &CellList<D>,
) -> Result<PairGraph, ContinuumError> {
    if configuration.is_empty() {
        return Err(ContinuumError::EmptyCoordinates { dimension: D });
    }
    list.validate(configuration)
        .map_err(|error| ContinuumError::CellListInconsistent(error.to_string()))?;
    Ok(enumerate_pairs(
        configuration.positions(),
        configuration.cell(),
        list,
    ))
}

/// Ephemeral owned snapshot of one candidate-pair analysis.
///
/// A vertex per particle and one physical edge per candidate pair `(i, j)`,
/// `i < j`, in the module's canonical order. This is the deliberate
/// transient counterpart to a persistent per-snapshot CSR (see the module
/// docs): it exists for the duration of one F3 [`crate::percolation::analyze`]
/// pass and is then dropped. Endpoint order is normalized `left < right`;
/// self-loops and parallel edges cannot occur because each unordered pair
/// is emitted exactly once.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairGraph {
    vertex_count: usize,
    pairs: Vec<[usize; 2]>,
}

impl PairGraph {
    /// Number of candidate pairs (physical edges).
    #[inline]
    pub const fn pair_count(&self) -> usize {
        self.pairs.len()
    }

    /// Candidate pairs in the canonical order documented at the module.
    #[inline]
    pub fn pairs(&self) -> &[[usize; 2]] {
        &self.pairs
    }

    /// Consume the snapshot into its flat pair list.
    #[inline]
    pub fn into_pairs(self) -> Vec<[usize; 2]> {
        self.pairs
    }
}

impl GraphView for PairGraph {
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }
}

impl UndirectedGraphView for PairGraph {
    fn edge_count(&self) -> usize {
        self.pairs.len()
    }

    fn edge_endpoints(&self, edge: EdgeId) -> [VertexId; 2] {
        let [left, right] = self.pairs[edge.index()];
        [
            VertexId::from_valid_index(left),
            VertexId::from_valid_index(right),
        ]
    }

    /// Incidences are derived by scanning the pair list — no adjacency
    /// index is materialized. One call costs `O(pair_count)`; callers
    /// needing repeated adjacency should build their own index. The F3
    /// component analyzer never calls this (it scans the physical edge
    /// table), which is exactly why the transient snapshot stays cheap.
    fn incidences(&self, vertex: VertexId) -> impl Iterator<Item = Incidence> + '_ {
        let target = vertex.index();
        self.pairs
            .iter()
            .enumerate()
            .filter_map(move |(index, &[left, right])| {
                if left == target {
                    Some(Incidence {
                        neighbor: VertexId::from_valid_index(right),
                        edge: EdgeId::from_valid_index(index),
                    })
                } else if right == target {
                    Some(Incidence {
                        neighbor: VertexId::from_valid_index(left),
                        edge: EdgeId::from_valid_index(index),
                    })
                } else {
                    None
                }
            })
    }
}

fn build_from_configuration<const D: usize>(
    configuration: &ParticleConfiguration<D>,
    cutoff: f64,
) -> Result<PairGraph, ContinuumError> {
    if configuration.is_empty() {
        return Err(ContinuumError::EmptyCoordinates { dimension: D });
    }
    let list = CellList::new(configuration, cutoff * cutoff).map_err(map_cell_list_error)?;
    Ok(enumerate_pairs(
        configuration.positions(),
        configuration.cell(),
        &list,
    ))
}

/// Copy flat coordinates into a throwaway particle configuration so the
/// reused `CellList` machinery can bucket them.
fn configuration_from_flat<const D: usize>(
    coordinates: &[f64],
    cell: &OrthorhombicCell<D>,
) -> Result<ParticleConfiguration<D>, ContinuumError> {
    let mut positions = Vec::with_capacity(coordinates.len() / D);
    for particle in 0..coordinates.len() / D {
        let mut position = [0.0; D];
        position.copy_from_slice(&coordinates[particle * D..(particle + 1) * D]);
        positions.push(position);
    }
    ParticleConfiguration::new(positions, vec![0; coordinates.len() / D], *cell).map_err(|error| {
        match error {
            ParticleError::NonFinitePosition { particle, axis } => {
                ContinuumError::NonFiniteCoordinate { particle, axis }
            }
            other => ContinuumError::CellListInconsistent(other.to_string()),
        }
    })
}

/// The stencil walk: cell-major, own-cell first endpoint ascending, partner
/// ascending. See the module docs for the full order contract and the
/// stencil-completeness argument.
fn enumerate_pairs<const D: usize>(
    positions: &[[f64; D]],
    cell: &OrthorhombicCell<D>,
    list: &CellList<D>,
) -> PairGraph {
    let cutoff_squared = list.cutoff_squared();
    let total_cells: usize = list.cell_counts().iter().product();
    let mut candidates: Vec<usize> = Vec::new();
    let mut pairs: Vec<[usize; 2]> = Vec::new();
    for cell_index in 0..total_cells {
        list.fill_candidates(cell_index, &mut candidates);
        if candidates.len() < 2 {
            continue;
        }
        // Buckets are disjoint, so sorting the neighborhood once per cell
        // yields ascending own-cell particles and ascending partners.
        candidates.sort_unstable();
        for (slot, &first) in candidates.iter().enumerate() {
            if list.particle_cell(first) != cell_index {
                continue;
            }
            let origin = positions[first];
            for &second in &candidates[slot + 1..] {
                // Borrow the particle module's own metric (review follow-up
                // 2026-10-23): identical fold order, zero duplicated formula.
                let distance_squared = cell.distance_squared(&origin, &positions[second]);
                if distance_squared <= cutoff_squared {
                    pairs.push([first, second]);
                }
            }
        }
    }
    PairGraph {
        vertex_count: positions.len(),
        pairs,
    }
}

fn map_cell_list_error(error: ParticleError) -> ContinuumError {
    match error {
        ParticleError::CutoffTooLarge { cutoff, maximum } => {
            ContinuumError::CutoffTooLarge { cutoff, maximum }
        }
        other => ContinuumError::CellListInconsistent(other.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::OrthorhombicCell;

    fn flat<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
        positions.iter().flat_map(|p| p.iter().copied()).collect()
    }

    #[test]
    fn order_is_cell_major_then_ascending() {
        // Dyadic coordinates (box 2.5, cutoff 1.25 -> 2x2 cells of width
        // 1.25), so every expectation below is exact in floating point.
        // Particle 5 and particle 0 sit in cell 0 and cell 3 respectively,
        // so first endpoints are NOT globally ascending — cell-major order
        // is the contract.
        let cell = OrthorhombicCell::new([2.5, 2.5]).unwrap();
        let positions = [
            [1.75, 1.75],   // particle 0 -> cell (1,1) = 3
            [0.25, 2.375],  // particle 1 -> cell (0,1) = 2
            [0.25, 0.25],   // particle 2 -> cell (0,0) = 0
            [1.75, 1.625],  // particle 3 -> cell (1,1) = 3
            [2.375, 0.125], // particle 4 -> cell (1,0) = 1
            [0.375, 0.375], // particle 5 -> cell (0,0) = 0
        ];
        let graph = CandidatePairGraph::new(&flat(&positions), &cell, 1.25)
            .unwrap()
            .pairs()
            .unwrap();
        // Hand-enumerated with minimum-image wrapping (dyadic, exact):
        // cell 0 emits (2,4),(2,5); cell 1 emits (4,5); cell 2 emits
        // (1,2),(1,3),(1,4),(1,5); cell 3 emits (0,1),(0,3),(0,4),(3,4).
        // Pair (1,3) sits exactly ON the cutoff through a periodic image
        // (displacement (-1.0, -0.75), squared length 1.5625) — an
        // incidental check that the boundary is inclusive.
        assert_eq!(
            graph.pairs(),
            &[
                [2usize, 4usize],
                [2, 5],
                [4, 5],
                [1, 2],
                [1, 3],
                [1, 4],
                [1, 5],
                [0, 1],
                [0, 3],
                [0, 4],
                [3, 4],
            ]
        );
        let first_endpoints: Vec<usize> = graph.pairs().iter().map(|p| p[0]).collect();
        let cells: Vec<usize> = first_endpoints
            .iter()
            .map(|&p| graph_cell(&positions[p]))
            .collect();
        assert_eq!(cells, vec![0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 3]);
        assert!(cells.windows(2).all(|w| w[0] <= w[1]));
    }

    /// Linear cell of a wrapped position on the 2x2 grid above, matching
    /// the `CellList` stride convention (axis 0 least significant).
    fn graph_cell(position: &[f64; 2]) -> usize {
        (position[0] / 1.25) as usize + 2 * (position[1] / 1.25) as usize
    }

    #[test]
    fn derived_incidences_agree_with_endpoints() {
        let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
        let positions = [[0.5, 0.5], [1.0, 0.5], [3.9, 3.9]];
        let graph = CandidatePairGraph::new(&flat(&positions), &cell, 1.0)
            .unwrap()
            .pairs()
            .unwrap();
        assert_eq!(graph.pair_count(), 2);
        let mut total = 0;
        for vertex in 0..graph.vertex_count() {
            for incidence in graph.incidences(VertexId::from_valid_index(vertex)) {
                let [left, right] = graph.edge_endpoints(incidence.edge);
                let other = if left.index() == vertex {
                    right
                } else {
                    debug_assert_eq!(right.index(), vertex);
                    left
                };
                assert_eq!(other, incidence.neighbor);
                total += 1;
            }
        }
        assert_eq!(total, 2 * graph.edge_count());
    }

    #[test]
    fn rejects_empty_and_mismatched_buffers() {
        let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
        assert!(matches!(
            CandidatePairGraph::<2>::new(&[], &cell, 1.0),
            Err(ContinuumError::EmptyCoordinates { dimension: 2 })
        ));
        assert!(matches!(
            CandidatePairGraph::<2>::new(&[1.0, 2.0, 3.0], &cell, 1.0),
            Err(ContinuumError::CoordinateCountMismatch {
                length: 3,
                dimension: 2
            })
        ));
        assert!(matches!(
            CandidatePairGraph::<2>::new(&[1.0, f64::NAN], &cell, 1.0),
            Err(ContinuumError::NonFiniteCoordinate {
                particle: 0,
                axis: 1
            })
        ));
    }
}
