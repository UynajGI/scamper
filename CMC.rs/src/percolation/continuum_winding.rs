//! Continuum winding embedding (roadmap v2 §2 T4) — the W2 displacement
//! source for P1's ephemeral pair graphs.
//!
//! [`analyze_winding`](crate::analyze_winding) is generic over a
//! [`PeriodicEmbedding`]: it needs one integer `Z^D` cell vector per physical
//! edge and nothing else. Lattice builders record those vectors at bond
//! emission; a continuum configuration has no bonds, only positions, so the
//! vector must be *derived* — this module defines and validates that
//! derivation. It stores data and implements no analysis (W2's analyzer is
//! consumed unchanged, exactly as the v1 W-track note demanded: continuum
//! wrapping may reuse the W2 contract "only after their own geometry and
//! physical configurations are defined" — the definition below is that
//! geometry).
//!
//! # The integer cell offset of one pair edge
//!
//! For a pair `(i, j)` of particles at caller coordinates `x_i`, `x_j` in an
//! [`OrthorhombicCell`] with side lengths `L`, define per axis `c`
//!
//! - the raw in-buffer delta `r_c = x_j[c] - x_i[c]` (the caller's
//!   representatives, wrapped or not), and
//! - the minimum-image displacement `d_c` = `cell.displacement(x_i, x_j)[c]`
//!   (the particle module's `to - from` convention, `|d_c| <= L_c / 2`).
//!
//! Both describe the same pair, so they differ by an integer number of box
//! lengths: `n_c = (d_c - r_c) / L_c` must be (numerically) an integer — for
//! in-box representatives, the number of `c`-boundaries crossed by the
//! *minimum-image* connection from `x_i` to `x_j` (with general
//! representatives it shifts by the boundaries of the representative path
//! itself; see the gauge note below). The stored vector is
//! `n = (n_0, ..., n_{D-1})` in the `i -> j` direction.
//!
//! ## Rounding and tolerance convention
//!
//! Each `n_c` is `round((d_c - r_c) / L_c)` with a **tolerance**, never a
//! silent rounding: the dimensionless ratio must lie within `1e-9` of an
//! integer, else construction fails with
//! [`ContinuumWindingError::NonIntegerCellOffset`] naming the edge and axis.
//! For coordinates produced by any consistent wrapping convention the ratio
//! sits within a few ulps of its integer (the subtraction telescopes through
//! the same minimum-image arithmetic), so `1e-9` is generous headroom; it is
//! only exceeded by inconsistent representatives separated by ~`1e8` box
//! lengths or more, which is numerical garbage and fails loudly. Ratios that
//! are non-finite, or whose integer exceeds `2 * 10^9` in magnitude (far
//! beyond any physical crossing count; also the exact-negation bound for
//! `i32`), are rejected by the same typed error — no saturation, no wrap.
//!
//! ## Negation convention
//!
//! The vector is stored once, in the pair's `i -> j` direction (P1 emits
//! `i < j`); querying [`PeriodicEmbedding::edge_displacement`] from `j`
//! returns the exact negation through the lazy
//! [`DirectedDisplacement`](crate::DirectedDisplacement) view — the W1
//! antisymmetry contract, structural rather than re-validated.
//!
//! ## Representative (gauge) freedom
//!
//! The coordinates handed to [`ContinuumEmbedding::from_pairs`] must be the
//! *same buffer* used to enumerate the pairs, so every particle has exactly
//! one representative. Which representative each particle uses is free:
//! re-wrapping any particle shifts the `n` of its incident edges by a
//! telescoping amount, leaves every cycle sum (hence every winding fact)
//! unchanged, and is the continuum analogue of the W2 gauge invariance. No
//! cocycle consistency is re-validated here — it is a producer property, as
//! for [`LatticeEmbedding`](crate::LatticeEmbedding).
//!
//! # Load-bearing precondition: `r_max < L_c / 2` on every axis
//!
//! The constructor rejects a cutoff `>= L_c / 2` on any axis `c`, naming the
//! axis. With `r_max < L_c / 2`, a pair within cutoff has a *unique* closest
//! periodic image, so the minimum-image displacement is the displacement of
//! the one actual physical connection and `n` records exactly how many
//! boundaries that connection crosses; the winding subgroup computed from
//! these vectors is then the winding of the cluster as a geometric object in
//! the torus (in particular, winding around axis `c` implies a
//! BFS-unwrapped extent `>= L_c` on that axis — the direction pinned by the
//! T4 validation gates). Once `r_max >= L_c / 2` on some axis, two images of
//! the same pair can both lie within cutoff: the pair graph records ONE edge
//! whose two geometric crossings are distinct, and the minimum-image tie
//! silently selects one. Integer winding remains a well-defined cochain, but
//! it no longer describes a definite physical connection — sub-ulp coordinate
//! noise at the tie can flip verdicts — so the regime is rejected outright.
//! (P1's cell list accepts `r_max = L_min / 2` exactly, the minimum-image tie;
//! this contract is deliberately stricter and documents why.)
//!
//! # Allocation accounting
//!
//! Construction is a per-realization transient like the [`PairGraph`]
//! itself: two owned vectors (`[usize; 2]` endpoints plus one `i32` per axis
//! per edge — `24 + 4 D` bytes per edge) reserved once and dropped with the
//! analysis. Lookups are O(1) table reads; steady-state winding analysis on a
//! fixed pair graph and embedding allocates nothing (the T4 probe in
//! `benches/percolation_allocations.rs`).

use crate::particle::{OrthorhombicCell, SimulationCell};
use crate::topology::{DirectedDisplacement, EdgeId, PeriodicEmbedding, VertexId};

/// Tolerance on the dimensionless ratio `(d_c - r_c) / L_c` around its
/// nearest integer. Consistent coordinates land within a few ulps; anything
/// beyond this is inconsistent input and a typed rejection.
const OFFSET_TOLERANCE: f64 = 1e-9;

/// Largest accepted per-axis cell offset. Far above any physical crossing
/// count, and small enough that exact negation stays inside `i32`.
const MAX_CELL_OFFSET: f64 = 2.0e9;

/// Typed rejection of a continuum winding-embedding request.
#[derive(Debug, Clone, PartialEq)]
pub enum ContinuumWindingError {
    /// The const dimension is zero; no geometry exists.
    ZeroDimension,
    /// The cutoff was non-finite or not strictly positive.
    InvalidCutoff { cutoff: f64 },
    /// The flat coordinate length is not a multiple of the dimension.
    CoordinateCountMismatch { length: usize, dimension: usize },
    /// The coordinate buffer holds no particles.
    EmptyCoordinates { dimension: usize },
    /// One coordinate was non-finite.
    NonFiniteCoordinate { particle: usize, axis: usize },
    /// The cutoff reaches half a box side or beyond on `axis`: the
    /// `r_max < L_c / 2` winding precondition fails (see the module docs for
    /// why double-crossing edges make winding physically ambiguous).
    CutoffSpansHalfBox {
        axis: usize,
        cutoff: f64,
        half_length: f64,
    },
    /// A pair endpoint is not a particle of the coordinate buffer.
    EndpointOutOfRange {
        edge: usize,
        endpoint: usize,
        particle_count: usize,
    },
    /// A self-loop pair cannot carry a directed displacement keyed by
    /// endpoint (the W1 rejection, same reason as
    /// [`crate::EmbeddingError::SelfLoopEdge`]).
    SelfLoopEdge { edge: usize },
    /// The implied cell offset `(d_c - r_c) / L_c` of one (edge, axis) is not
    /// within [`OFFSET_TOLERANCE`] of an integer, is non-finite, or exceeds
    /// [`MAX_CELL_OFFSET`] in magnitude: numerical garbage in, loud failure
    /// out — never a silent rounding or saturation.
    NonIntegerCellOffset {
        edge: usize,
        axis: usize,
        ratio: f64,
    },
    /// The embedding table could not reserve its per-edge storage.
    Capacity(std::collections::TryReserveError),
}

impl core::fmt::Display for ContinuumWindingError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ZeroDimension => write!(formatter, "particle dimension must be positive"),
            Self::InvalidCutoff { cutoff } => write!(
                formatter,
                "winding-embedding cutoff must be finite and positive, got {cutoff}"
            ),
            Self::CoordinateCountMismatch { length, dimension } => write!(
                formatter,
                "flat coordinate length {length} is not a multiple of dimension {dimension}"
            ),
            Self::EmptyCoordinates { dimension } => write!(
                formatter,
                "coordinate buffer holds no {dimension}-dimensional particles"
            ),
            Self::NonFiniteCoordinate { particle, axis } => {
                write!(
                    formatter,
                    "particle {particle} coordinate {axis} is non-finite"
                )
            }
            Self::CutoffSpansHalfBox {
                axis,
                cutoff,
                half_length,
            } => write!(
                formatter,
                "cutoff {cutoff} reaches half the box side {half_length} on axis {axis}; \
                 winding requires cutoff < L/2 on every axis (double-crossing edges are \
                 physically ambiguous)"
            ),
            Self::EndpointOutOfRange {
                edge,
                endpoint,
                particle_count,
            } => write!(
                formatter,
                "pair {edge} endpoint {endpoint} is outside the {particle_count} particles \
                 of the coordinate buffer"
            ),
            Self::SelfLoopEdge { edge } => {
                write!(
                    formatter,
                    "pair {edge} is a self-loop and cannot carry a displacement"
                )
            }
            Self::NonIntegerCellOffset { edge, axis, ratio } => write!(
                formatter,
                "cell offset ratio {ratio} of pair {edge} axis {axis} is not an integer \
                 within tolerance {OFFSET_TOLERANCE} (or is non-finite/unrepresentable)"
            ),
            Self::Capacity(_) => formatter.write_str("could not reserve winding embedding"),
        }
    }
}

impl std::error::Error for ContinuumWindingError {}

/// Owned integer displacement table for one continuum pair graph (P1
/// [`PairGraph`](crate::PairGraph) or P2 `OverlapGraph` pair list).
///
/// Built from the same coordinate buffer the pairs were enumerated on, plus
/// the enumeration cutoff (needed to enforce the `r_max < L/2` precondition).
/// Implements [`PeriodicEmbedding`], so the unchanged
/// [`analyze_winding`](crate::analyze_winding) accepts it next to the pair
/// graph. Like [`LatticeEmbedding`](crate::LatticeEmbedding), the table is
/// detached data: graph and embedding must stay unchanged together.
#[derive(Clone, Debug)]
pub struct ContinuumEmbedding {
    dimensions: usize,
    endpoints: Vec<[usize; 2]>,
    displacements: Vec<i32>,
}

impl ContinuumEmbedding {
    /// Derive the per-edge integer cell vectors of a pair list.
    ///
    /// `pairs` is the physical edge table of a continuum graph in any order
    /// (`PairGraph::pairs` or `OverlapGraph::pairs`); `coordinates` must be
    /// the same flat buffer those pairs were enumerated from (one
    /// representative per particle — see the module's gauge note);
    /// `cutoff` is the enumeration cutoff, used only to enforce the
    /// `cutoff < L_c/2` winding precondition.
    ///
    /// Rejections, all before any state is kept: zero dimension, invalid
    /// cutoff, buffer length not a multiple of `D`, empty buffer, non-finite
    /// coordinates, the half-box precondition (axis named), out-of-range or
    /// self-loop endpoints, non-integer cell offsets (edge and axis named),
    /// unresolvable capacity.
    pub fn from_pairs<const D: usize>(
        pairs: &[[usize; 2]],
        coordinates: &[f64],
        cell: &OrthorhombicCell<D>,
        cutoff: f64,
    ) -> Result<Self, ContinuumWindingError> {
        if D == 0 {
            return Err(ContinuumWindingError::ZeroDimension);
        }
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(ContinuumWindingError::InvalidCutoff { cutoff });
        }
        if !coordinates.len().is_multiple_of(D) {
            return Err(ContinuumWindingError::CoordinateCountMismatch {
                length: coordinates.len(),
                dimension: D,
            });
        }
        if coordinates.is_empty() {
            return Err(ContinuumWindingError::EmptyCoordinates { dimension: D });
        }
        for (flat, &value) in coordinates.iter().enumerate() {
            if !value.is_finite() {
                return Err(ContinuumWindingError::NonFiniteCoordinate {
                    particle: flat / D,
                    axis: flat % D,
                });
            }
        }
        for (axis, &length) in cell.lengths().iter().enumerate() {
            if cutoff >= length / 2.0 {
                return Err(ContinuumWindingError::CutoffSpansHalfBox {
                    axis,
                    cutoff,
                    half_length: length / 2.0,
                });
            }
        }
        let particle_count = coordinates.len() / D;
        let mut endpoints = Vec::new();
        endpoints
            .try_reserve_exact(pairs.len())
            .map_err(ContinuumWindingError::Capacity)?;
        let mut table = Vec::new();
        table
            .try_reserve_exact(pairs.len() * D)
            .map_err(ContinuumWindingError::Capacity)?;
        for (edge, &[left, right]) in pairs.iter().enumerate() {
            if left == right {
                return Err(ContinuumWindingError::SelfLoopEdge { edge });
            }
            if left >= particle_count || right >= particle_count {
                return Err(ContinuumWindingError::EndpointOutOfRange {
                    edge,
                    endpoint: left.max(right),
                    particle_count,
                });
            }
            let mut origin = [0.0; D];
            let mut target = [0.0; D];
            origin.copy_from_slice(&coordinates[left * D..][..D]);
            target.copy_from_slice(&coordinates[right * D..][..D]);
            let displacement = cell.displacement(&origin, &target);
            for axis in 0..D {
                let raw = target[axis] - origin[axis];
                let ratio = (displacement[axis] - raw) / cell.lengths()[axis];
                table.push(integer_cell_offset(ratio, edge, axis)?);
            }
            endpoints.push([left, right]);
        }
        Ok(Self {
            dimensions: D,
            endpoints,
            displacements: table,
        })
    }

    /// Number of physical edges this table carries displacements for.
    ///
    /// Pair it with the graph before the infallible per-edge lookup (an
    /// out-of-range edge ID would panic there), the same consumer pattern as
    /// [`LatticeEmbedding::edge_count`](crate::LatticeEmbedding::edge_count).
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.endpoints.len()
    }
}

impl PeriodicEmbedding for ContinuumEmbedding {
    #[inline(always)]
    fn periodic_dimensions(&self) -> usize {
        self.dimensions
    }

    #[inline]
    fn edge_displacement(&self, edge: EdgeId, from: VertexId) -> Option<DirectedDisplacement<'_>> {
        let index = edge.index();
        let [source, target] = self.endpoints[index];
        let vector = &self.displacements[index * self.dimensions..][..self.dimensions];
        if from.index() == source {
            Some(DirectedDisplacement::new(vector, false))
        } else if from.index() == target {
            Some(DirectedDisplacement::new(vector, true))
        } else {
            None
        }
    }
}

/// Round one dimensionless cell-offset ratio to its integer, rejecting
/// anything not within [`OFFSET_TOLERANCE`] of an integer and any magnitude
/// beyond [`MAX_CELL_OFFSET`] (non-finite ratios fold into the same
/// rejection).
fn integer_cell_offset(ratio: f64, edge: usize, axis: usize) -> Result<i32, ContinuumWindingError> {
    if !ratio.is_finite() || ratio.abs() > MAX_CELL_OFFSET {
        return Err(ContinuumWindingError::NonIntegerCellOffset { edge, axis, ratio });
    }
    let rounded = ratio.round();
    if (ratio - rounded).abs() > OFFSET_TOLERANCE {
        return Err(ContinuumWindingError::NonIntegerCellOffset { edge, axis, ratio });
    }
    // |rounded| <= 2e9 < 2^31: the cast cannot saturate and the negated
    // direction of the same vector stays exactly representable.
    Ok(rounded as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
        positions.iter().flat_map(|p| p.iter().copied()).collect()
    }

    fn displacement_table(embedding: &ContinuumEmbedding) -> Vec<Vec<i32>> {
        let dimensions = embedding.periodic_dimensions();
        (0..embedding.edge_count())
            .map(|edge| embedding.displacements[edge * dimensions..][..dimensions].to_vec())
            .collect()
    }

    #[test]
    fn hand_derived_ring_in_1d_crosses_exactly_once() {
        // Box 10, cutoff 2: chain 1.0 - 3.0 - 5.0 - 7.0 - 9.5 - (wrap) 1.0.
        // Only the last pair crosses the boundary. It is stored in the
        // 0 -> 4 direction: raw 9.5 - 1.0 = 8.5, minimum image -1.5,
        // n = (-1.5 - 8.5)/10 = -1 (the wrap traversal 4 -> 0 is +1).
        let cell = OrthorhombicCell::new([10.0]).unwrap();
        let positions = [[1.0], [3.0], [5.0], [7.0], [9.5]];
        let pairs = [[0usize, 1usize], [1, 2], [2, 3], [3, 4], [0, 4]];
        let embedding =
            ContinuumEmbedding::from_pairs(&pairs, &flat(&positions), &cell, 2.0).unwrap();
        assert_eq!(
            displacement_table(&embedding),
            vec![vec![0], vec![0], vec![0], vec![0], vec![-1]],
            "hand table: pair [0, 4] is stored 0 -> 4, which crosses -1 \
             (the wrap traversal 4 -> 0 is its negation, +1)"
        );
    }

    #[test]
    fn hand_derived_two_particle_wrap_pair_and_2d_chain() {
        // 2D two-particle pair at wrap distance: raw delta (7.0, 0.0),
        // minimum image (-1.0, 0.0), n = (-1, 0).
        let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
        let pair_positions = [[0.5, 4.0], [7.5, 4.0]];
        let embedding =
            ContinuumEmbedding::from_pairs(&[[0usize, 1usize]], &flat(&pair_positions), &cell, 2.0)
                .unwrap();
        assert_eq!(displacement_table(&embedding), vec![vec![-1, 0]]);

        // A 4-chain along x that closes around the torus exactly once:
        // consecutive gaps 3, 3, 2.8, then 9.8 -> 0.2 wraps (image 0.4, n=+1),
        // and 0.2 -> 1.0 (0.8, n=0) closes the cycle.
        let chain = [[1.0, 5.0], [4.0, 5.0], [7.0, 5.0], [9.8, 5.0], [0.2, 5.0]];
        let pairs = [[0usize, 1usize], [1, 2], [2, 3], [3, 4], [0, 4]];
        let embedding = ContinuumEmbedding::from_pairs(&pairs, &flat(&chain), &cell, 3.1).unwrap();
        assert_eq!(
            displacement_table(&embedding),
            vec![vec![0, 0], vec![0, 0], vec![0, 0], vec![1, 0], vec![0, 0]],
            "hand table: exactly the 9.8 -> 0.2 edge crosses +1"
        );
    }

    #[test]
    fn displacement_is_exact_negation_and_non_endpoint_is_none() {
        let cell = OrthorhombicCell::new([6.0, 5.0]).unwrap();
        let positions = [[0.3, 0.2], [5.9, 4.8], [3.0, 2.5]];
        let pairs = [[0usize, 1usize], [0, 2]];
        let embedding =
            ContinuumEmbedding::from_pairs(&pairs, &flat(&positions), &cell, 2.0).unwrap();
        for (edge_index, &[left, right]) in pairs.iter().enumerate() {
            let edge = EdgeId::from_valid_index(edge_index);
            let forward = embedding
                .edge_displacement(edge, VertexId::from_valid_index(left))
                .expect("left endpoint");
            let backward = embedding
                .edge_displacement(edge, VertexId::from_valid_index(right))
                .expect("right endpoint");
            for axis in 0..2 {
                assert_eq!(
                    forward.axis(axis).unwrap(),
                    -backward.axis(axis).unwrap(),
                    "edge {edge_index} axis {axis}: exact negation"
                );
            }
        }
        assert!(
            embedding
                .edge_displacement(EdgeId::from_valid_index(0), VertexId::from_valid_index(2))
                .is_none(),
            "vertex 2 is not an endpoint of edge 0"
        );
    }

    #[test]
    fn cutoff_reaching_half_a_side_is_rejected_with_the_axis() {
        let cell = OrthorhombicCell::new([4.0, 10.0]).unwrap();
        let positions = [[0.0, 0.0], [1.0, 1.0]];
        let coordinates = flat(&positions);
        let pairs = [[0usize, 1usize]];
        // Exactly half the x side: rejected even though P1's cell list would
        // accept the tie.
        assert_eq!(
            ContinuumEmbedding::from_pairs(&pairs, &coordinates, &cell, 2.0).unwrap_err(),
            ContinuumWindingError::CutoffSpansHalfBox {
                axis: 0,
                cutoff: 2.0,
                half_length: 2.0
            }
        );
        // Way past half the x side: axis 0, not the longer axis 1.
        assert_eq!(
            ContinuumEmbedding::from_pairs(&pairs, &coordinates, &cell, 5.0).unwrap_err(),
            ContinuumWindingError::CutoffSpansHalfBox {
                axis: 0,
                cutoff: 5.0,
                half_length: 2.0
            }
        );
        // Comfortably inside x but reaching half the y side: axis 1 named.
        let tall_cell = OrthorhombicCell::new([10.0, 4.0]).unwrap();
        assert_eq!(
            ContinuumEmbedding::from_pairs(&pairs, &coordinates, &tall_cell, 2.0).unwrap_err(),
            ContinuumWindingError::CutoffSpansHalfBox {
                axis: 1,
                cutoff: 2.0,
                half_length: 2.0
            }
        );
        // Just under half the shortest side still constructs.
        ContinuumEmbedding::from_pairs(&pairs, &coordinates, &cell, 2.0 - 1e-9)
            .expect("cutoff just under L/2 is the boundary of the valid regime");
    }

    #[test]
    fn fractional_and_unrepresentable_offsets_are_typed_rejections() {
        // Public-construction reachability of the tolerance: representatives
        // 1e8 box lengths apart make the minimum-image arithmetic lose the
        // integer to rounding (residual ~1.5e-8 > 1e-9). Deterministic float
        // behavior, verified against the exact production formula.
        let cell = OrthorhombicCell::new([0.1]).unwrap();
        let positions = [[0.0], [10_000_000.05]];
        let error =
            ContinuumEmbedding::from_pairs(&[[0usize, 1usize]], &flat(&positions), &cell, 0.04)
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

        // Public-construction reachability of the magnitude bound: a
        // representative 1e300 boxes away yields an integer-valued but
        // unrepresentable ratio — rejected, never saturated to i32::MAX.
        let cell = OrthorhombicCell::new([1.0]).unwrap();
        let positions = [[0.0], [1e300]];
        let error =
            ContinuumEmbedding::from_pairs(&[[0usize, 1usize]], &flat(&positions), &cell, 0.4)
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
    }

    #[test]
    fn offset_helper_rounds_within_tolerance_and_rejects_garbage() {
        use ContinuumWindingError::NonIntegerCellOffset;
        assert_eq!(integer_cell_offset(0.0, 0, 0), Ok(0));
        assert_eq!(integer_cell_offset(-3.0, 0, 0), Ok(-3));
        assert_eq!(integer_cell_offset(2.0 + 5e-10, 0, 0), Ok(2));
        assert_eq!(integer_cell_offset(-1.0 - 5e-10, 0, 0), Ok(-1));
        // The task's crafted case: half-away-from-integer residual.
        assert_eq!(
            integer_cell_offset(3.4999, 7, 2),
            Err(NonIntegerCellOffset {
                edge: 7,
                axis: 2,
                ratio: 3.4999
            })
        );
        assert_eq!(
            integer_cell_offset(-2.5000001, 0, 0).unwrap_err(),
            NonIntegerCellOffset {
                edge: 0,
                axis: 0,
                ratio: -2.5000001
            }
        );
        // Non-finite and out-of-range magnitudes share the typed rejection.
        assert!(matches!(
            integer_cell_offset(f64::NAN, 0, 0),
            Err(NonIntegerCellOffset { edge: 0, axis: 0, ratio }) if ratio.is_nan()
        ));
        assert!(matches!(
            integer_cell_offset(f64::INFINITY, 0, 0),
            Err(NonIntegerCellOffset { .. })
        ));
        assert!(matches!(
            integer_cell_offset(3.0e9, 0, 0),
            Err(NonIntegerCellOffset { .. })
        ));
    }

    #[test]
    fn structurally_invalid_requests_are_rejected() {
        let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
        let coordinates = [0.0, 0.0, 1.0, 1.0];
        assert_eq!(
            ContinuumEmbedding::from_pairs::<2>(&[], &coordinates, &cell, 0.0).unwrap_err(),
            ContinuumWindingError::InvalidCutoff { cutoff: 0.0 }
        );
        assert_eq!(
            ContinuumEmbedding::from_pairs::<2>(&[], &[0.0, 0.0, 1.0], &cell, 1.0).unwrap_err(),
            ContinuumWindingError::CoordinateCountMismatch {
                length: 3,
                dimension: 2
            }
        );
        assert_eq!(
            ContinuumEmbedding::from_pairs::<2>(&[], &[], &cell, 1.0).unwrap_err(),
            ContinuumWindingError::EmptyCoordinates { dimension: 2 }
        );
        assert_eq!(
            ContinuumEmbedding::from_pairs::<2>(&[], &[0.0, f64::NAN], &cell, 1.0).unwrap_err(),
            ContinuumWindingError::NonFiniteCoordinate {
                particle: 0,
                axis: 1
            }
        );
        assert_eq!(
            ContinuumEmbedding::from_pairs(&[[0usize, 0usize]], &coordinates, &cell, 1.0)
                .unwrap_err(),
            ContinuumWindingError::SelfLoopEdge { edge: 0 }
        );
        assert_eq!(
            ContinuumEmbedding::from_pairs(&[[0usize, 2usize]], &coordinates, &cell, 1.0)
                .unwrap_err(),
            ContinuumWindingError::EndpointOutOfRange {
                edge: 0,
                endpoint: 2,
                particle_count: 2
            }
        );
    }

    #[test]
    fn representative_rewrapping_leaves_the_cycle_sums_unchanged() {
        // Gauge check at the data level: re-wrapping particle 1 by one box
        // length changes its incident edges' stored vectors by a telescoping
        // amount, so the cycle sum of the triangle is identical.
        let cell = OrthorhombicCell::new([10.0]).unwrap();
        let positions = [[1.0], [5.0], [9.0]];
        let pairs = [[0usize, 1usize], [1, 2], [0, 2]];
        let base = ContinuumEmbedding::from_pairs(&pairs, &flat(&positions), &cell, 4.1).unwrap();
        let rewrapped_positions = [[1.0], [15.0], [9.0]];
        let rewrapped =
            ContinuumEmbedding::from_pairs(&pairs, &flat(&rewrapped_positions), &cell, 4.1)
                .unwrap();
        let base_table = displacement_table(&base);
        let rewrapped_table = displacement_table(&rewrapped);
        // The triangle 0 -> 1 -> 2 -> 0 traverses the stored vectors 0->1,
        // 1->2 forward and 0->2 backward.
        let cycle = |table: &[Vec<i32>]| table[0][0] + table[1][0] - table[2][0];
        assert_eq!(cycle(&base_table), cycle(&rewrapped_table));
        // Sanity on the hand table: only 0 -> 2 wraps (raw 8.0, image -2.0),
        // so the cycle crosses the boundary exactly once.
        assert_eq!(base_table, vec![vec![0], vec![0], vec![-1]]);
        assert_eq!(cycle(&base_table), 1);
        // The rewrapped table itself differs (edges at particle 1 shifted):
        assert_eq!(rewrapped_table, vec![vec![-1], vec![1], vec![-1]]);
    }
}
