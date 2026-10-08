use core::fmt;

/// Validation failure for a borrowed undirected CSR representation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TopologyError {
    OffsetLength {
        expected: usize,
        actual: usize,
    },
    OffsetStart {
        actual: usize,
    },
    OffsetsNotMonotone {
        index: usize,
    },
    OffsetEnd {
        expected: usize,
        actual: usize,
    },
    IncidenceLength {
        neighbors: usize,
        edge_ids: usize,
    },
    NeighborOutOfRange {
        incidence: usize,
        neighbor: usize,
        vertex_count: usize,
    },
    EdgeIdOutOfRange {
        incidence: usize,
        edge: usize,
        edge_count: usize,
    },
    EndpointOutOfRange {
        edge: usize,
        endpoint: usize,
        vertex_count: usize,
    },
    IncidenceMismatch {
        incidence: usize,
        vertex: usize,
        neighbor: usize,
        edge: usize,
    },
    IncidenceMultiplicity {
        edge: usize,
        actual: usize,
    },
    EndpointIncidenceMultiplicity {
        edge: usize,
        source_seen: bool,
        target_seen: bool,
    },
    /// Validation scratch could not reserve one byte per physical edge.
    ///
    /// This covers capacity overflow and allocator errors reported by
    /// `try_reserve_exact`; an allocator configured to abort on OOM may still
    /// terminate the process before Rust can return an error.
    ValidationCapacity {
        edge_count: usize,
    },
}

impl fmt::Display for TopologyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::OffsetLength { expected, actual } => {
                write!(
                    formatter,
                    "offset length mismatch: expected {expected}, got {actual}"
                )
            }
            Self::OffsetStart { actual } => {
                write!(formatter, "offsets must start at zero, got {actual}")
            }
            Self::OffsetsNotMonotone { index } => {
                write!(formatter, "offsets are decreasing at index {index}")
            }
            Self::OffsetEnd { expected, actual } => write!(
                formatter,
                "last offset mismatch: expected {expected} incidences, got {actual}"
            ),
            Self::IncidenceLength {
                neighbors,
                edge_ids,
            } => write!(
                formatter,
                "incidence array length mismatch: {neighbors} neighbors, {edge_ids} edge IDs"
            ),
            Self::NeighborOutOfRange {
                incidence,
                neighbor,
                vertex_count,
            } => write!(
                formatter,
                "incidence {incidence} neighbor {neighbor} is outside [0, {vertex_count})"
            ),
            Self::EdgeIdOutOfRange {
                incidence,
                edge,
                edge_count,
            } => write!(
                formatter,
                "incidence {incidence} edge ID {edge} is outside [0, {edge_count})"
            ),
            Self::EndpointOutOfRange {
                edge,
                endpoint,
                vertex_count,
            } => write!(
                formatter,
                "edge {edge} endpoint {endpoint} is outside [0, {vertex_count})"
            ),
            Self::IncidenceMismatch {
                incidence,
                vertex,
                neighbor,
                edge,
            } => write!(
                formatter,
                "incidence {incidence} ({vertex}, {neighbor}) does not match edge {edge}"
            ),
            Self::IncidenceMultiplicity { edge, actual } => write!(
                formatter,
                "physical edge {edge} must have exactly two incidences, got {actual}"
            ),
            Self::EndpointIncidenceMultiplicity {
                edge,
                source_seen,
                target_seen,
            } => write!(
                formatter,
                "physical edge {edge} has invalid endpoint incidence state \
                 [source_seen={source_seen}, target_seen={target_seen}]"
            ),
            Self::ValidationCapacity { edge_count } => write!(
                formatter,
                "could not reserve validation state for {edge_count} physical edges"
            ),
        }
    }
}

impl std::error::Error for TopologyError {}
