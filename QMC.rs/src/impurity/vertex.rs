//! Retarded four-leg vertices.

use super::error::ImpurityError;

/// Stable slot identifier for a vertex in the configuration.
///
/// Unlike a raw `Vec` index, a `VertexId` remains valid across insertions
/// and deletions because vertices occupy fixed slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(pub usize);

/// Identifies one of the two endpoints (A or B) of a retarded vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EndpointId {
    /// The vertex that owns this endpoint.
    pub vertex: VertexId,
    /// `0` for endpoint A, `1` for endpoint B.
    pub endpoint: u8,
}

/// Which side of an endpoint a leg belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LegSide {
    /// Incoming leg (toward earlier imaginary time).
    Incoming,
    /// Outgoing leg (toward later imaginary time).
    Outgoing,
}

/// Fully qualified identifier for one of the four legs of a retarded vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LegId {
    /// The endpoint that owns this leg.
    pub endpoint: EndpointId,
    /// Whether this is the incoming or outgoing leg.
    pub side: LegSide,
}

impl LegId {
    /// Construct a `LegId` from a vertex and a local leg index.
    ///
    /// Local leg numbering:
    /// - `0` = `A_IN`  (endpoint A, incoming)
    /// - `1` = `A_OUT` (endpoint A, outgoing)
    /// - `2` = `B_IN`  (endpoint B, incoming)
    /// - `3` = `B_OUT` (endpoint B, outgoing)
    pub fn from_local(vertex: VertexId, local_leg: usize) -> Self {
        match local_leg {
            A_IN => Self {
                endpoint: EndpointId {
                    vertex,
                    endpoint: 0,
                },
                side: LegSide::Incoming,
            },
            A_OUT => Self {
                endpoint: EndpointId {
                    vertex,
                    endpoint: 0,
                },
                side: LegSide::Outgoing,
            },
            B_IN => Self {
                endpoint: EndpointId {
                    vertex,
                    endpoint: 1,
                },
                side: LegSide::Incoming,
            },
            B_OUT => Self {
                endpoint: EndpointId {
                    vertex,
                    endpoint: 1,
                },
                side: LegSide::Outgoing,
            },
            _ => panic!("invalid local leg index: {local_leg}"),
        }
    }

    /// Convert back to the local leg index (0..4).
    pub fn local_leg(self) -> usize {
        match (self.endpoint.endpoint, self.side) {
            (0, LegSide::Incoming) => A_IN,
            (0, LegSide::Outgoing) => A_OUT,
            (1, LegSide::Incoming) => B_IN,
            (1, LegSide::Outgoing) => B_OUT,
            _ => unreachable!("invalid endpoint value: {}", self.endpoint.endpoint),
        }
    }
}

/// Spin state stored on a worldline leg (`-1` or `+1`, corresponding to
/// `sigma_z`).
pub type Spin = i8;

/// Local leg numbering of one retarded vertex.
pub const A_IN: usize = 0;
/// Outgoing leg at endpoint A.
pub const A_OUT: usize = 1;
/// Incoming leg at endpoint B.
pub const B_IN: usize = 2;
/// Outgoing leg at endpoint B.
pub const B_OUT: usize = 3;
/// Number of legs on a retarded vertex.
pub const LEGS_PER_VERTEX: usize = 4;

/// Immutable local vertex type supplied by a impurity model.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexKind {
    name: String,
    legs: [Spin; LEGS_PER_VERTEX],
    weight: f64,
    diagonal: bool,
}

impl VertexKind {
    /// Construct a positive local vertex type.
    pub fn new(
        name: impl Into<String>,
        legs: [Spin; LEGS_PER_VERTEX],
        weight: f64,
        diagonal: bool,
    ) -> Result<Self, ImpurityError> {
        if legs.iter().any(|spin| !matches!(spin, -1 | 1)) {
            return Err(ImpurityError::parameter(
                "vertex legs",
                "spin-1/2 legs must be encoded as -1 or +1",
            ));
        }
        if !weight.is_finite() || weight <= 0.0 {
            return Err(ImpurityError::parameter(
                "vertex weight",
                format!("must be finite and positive, got {weight}"),
            ));
        }
        let inferred_diagonal = legs[A_IN] == legs[A_OUT] && legs[B_IN] == legs[B_OUT];
        if inferred_diagonal != diagonal {
            return Err(ImpurityError::parameter(
                "diagonal",
                "diagonal flag does not match the leg pattern",
            ));
        }
        Ok(Self {
            name: name.into(),
            legs,
            weight,
            diagonal,
        })
    }

    /// Human-readable vertex name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Four local spin legs.
    pub fn legs(&self) -> &[Spin; LEGS_PER_VERTEX] {
        &self.legs
    }

    /// Positive local matrix-element weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// Whether both endpoint operators are diagonal.
    pub fn is_diagonal(&self) -> bool {
        self.diagonal
    }

    /// Spin on one local leg.
    pub fn spin(&self, leg: usize) -> Spin {
        self.legs[leg]
    }
}

/// One sampled retarded interaction vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    /// First endpoint time.
    pub tau_a: f64,
    /// Second endpoint time.
    pub tau_b: f64,
    /// Sampled bath frequency.
    pub omega: f64,
    /// Interaction-channel index.
    pub interaction: usize,
    /// Local kind index inside the interaction channel.
    pub kind: usize,
}

/// One time-ordered endpoint in the worldline index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Event {
    /// Imaginary time.
    pub time: f64,
    /// Vertex index.
    pub vertex: usize,
    /// Endpoint (`0` for A, `1` for B).
    pub endpoint: usize,
}

impl Event {
    /// Incoming global leg index.
    pub fn incoming_leg(self) -> usize {
        LEGS_PER_VERTEX * self.vertex + 2 * self.endpoint
    }

    /// Outgoing global leg index.
    pub fn outgoing_leg(self) -> usize {
        self.incoming_leg() + 1
    }
}
