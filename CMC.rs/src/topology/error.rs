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
    /// An out-incidence references an arc that does not start there:
    /// `arc_endpoints[arc][0] != vertex` or the neighbor is not the target.
    OutArcMismatch {
        incidence: usize,
        vertex: usize,
        neighbor: usize,
        arc: usize,
    },
    /// An in-incidence references an arc that does not end there:
    /// `arc_endpoints[arc][1] != vertex` or the neighbor is not the source.
    InArcMismatch {
        incidence: usize,
        vertex: usize,
        neighbor: usize,
        arc: usize,
    },
    /// An arc appears `actual != 1` times across the out incidences
    /// (`0` = missing, `>= 2` = duplicated).
    ArcOutMultiplicity {
        arc: usize,
        actual: usize,
    },
    /// An arc appears `actual != 1` times across the transpose incidences.
    ArcInMultiplicity {
        arc: usize,
        actual: usize,
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
            Self::OutArcMismatch {
                incidence,
                vertex,
                neighbor,
                arc,
            } => write!(
                formatter,
                "out-incidence {incidence} at vertex {vertex} references arc {arc} \
                 that does not run {vertex} -> {neighbor}"
            ),
            Self::InArcMismatch {
                incidence,
                vertex,
                neighbor,
                arc,
            } => write!(
                formatter,
                "in-incidence {incidence} at vertex {vertex} references arc {arc} \
                 that does not run {neighbor} -> {vertex}"
            ),
            Self::ArcOutMultiplicity { arc, actual } => write!(
                formatter,
                "arc {arc} appears {actual} times in the out incidences (expected exactly 1)"
            ),
            Self::ArcInMultiplicity { arc, actual } => write!(
                formatter,
                "arc {arc} appears {actual} times in the transpose incidences (expected exactly 1)"
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
