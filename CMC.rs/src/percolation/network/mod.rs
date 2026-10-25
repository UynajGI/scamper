//! N2: degree-conditioned observables and random network failure.

use core::fmt;

mod attack;
mod failure;
mod input;
#[cfg(feature = "io")]
mod io;
mod multiplex;
mod profile;

pub use attack::{AttackGrid, AttackPolicy, DegreeAttackMC};
pub use failure::{DegreeScratch, RandomFailureMC};
pub use input::{
    csr_from_edge_list, DuplicateEdgePolicy, EdgeListSpec, NetworkInputError, SelfLoopPolicy,
};
#[cfg(feature = "io")]
pub use io::{
    NetworkDocument, NetworkEdge, NetworkIoError, EDGE_LIST_HEADER_V1, NETWORK_JSON_FORMAT,
    NETWORK_JSON_VERSION,
};
pub use multiplex::{multiplex_from_edges, MultiplexError, MultiplexNetwork, MutualComponent};
pub use profile::{NetworkError, NetworkProfile};

/// Typed rejections for the N2 facade.
#[derive(Debug)]
pub enum NetworkPercolationError {
    /// Keep probability outside `[0, 1]`.
    Probability(&'static str, super::ProbabilityError),
    /// The lattice representation failed validation.
    Topology(String),
    /// The plan requested an unsupported observable (winding) or violated
    /// the exact-conversion moment gate.
    Observable(super::ObservablePlanError),
    /// The plan requested a winding observable, which needs a periodic
    /// embedding this facade does not carry; use [`crate::StaticPercolationMC::with_embedding`].
    WindingUnsupported,
    /// A boundary query was built for a different vertex count.
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// Degree-conditioned second moments sum component sizes squared per
    /// vertex (`sum_c s_c^3 <= V^3`), so `V^3` must stay within the exact
    /// `f64` integer range for every recorded value to stay exact.
    InexactDegreeMoment { total_vertices: usize },
    /// A scratch buffer could not be reserved.
    Capacity,
    /// The graph has no vertices to attack.
    EmptyGraph,
    /// An `AttackGrid::Every` step was zero.
    InvalidGrid(&'static str),
    /// An `AttackGrid::At` removal count lies outside `0..=V`.
    GridDomain { removal: usize, vertices: usize },
}

impl fmt::Display for NetworkPercolationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Probability(field, source) => {
                write!(formatter, "invalid {field}: {source}")
            }
            Self::Topology(reason) => write!(formatter, "invalid topology: {reason}"),
            Self::WindingUnsupported => write!(
                formatter,
                "winding observables need a periodic embedding; \
                 use StaticPercolationMC::with_embedding"
            ),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} was validated for {actual} vertices, expected {expected}"
            ),
            Self::InexactDegreeMoment { total_vertices } => write!(
                formatter,
                "degree second moments reach V^3 > 2^53 for V={total_vertices}; \
                 exact f64 conversion is not guaranteed"
            ),
            Self::Capacity => write!(formatter, "could not reserve network scratch storage"),
            Self::EmptyGraph => write!(formatter, "attack needs at least one vertex"),
            Self::InvalidGrid(reason) => write!(formatter, "invalid attack grid: {reason}"),
            Self::GridDomain { removal, vertices } => write!(
                formatter,
                "attack grid checkpoint at removal {removal} lies outside 0..={vertices}"
            ),
        }
    }
}

impl std::error::Error for NetworkPercolationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Probability(_, source) => Some(source),
            Self::Observable(source) => Some(source),
            _ => None,
        }
    }
}
