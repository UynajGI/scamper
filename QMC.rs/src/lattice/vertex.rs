//! Sparse positive local operators and sampled continuous-time vertices.

use crate::local_space::BasisState;

use super::error::LatticeQmcError;

/// Positive local matrix element represented by its incoming/outgoing legs.
///
/// Legs are ordered `[site0_in, site0_out, site1_in, site1_out, ...]`.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexKind {
    name: String,
    legs: Box<[BasisState]>,
    weight: f64,
    diagonal: bool,
}

impl VertexKind {
    /// Construct a positive local matrix element.
    pub fn new(
        name: impl Into<String>,
        legs: Vec<BasisState>,
        weight: f64,
    ) -> Result<Self, LatticeQmcError> {
        if legs.is_empty() || !legs.len().is_multiple_of(2) {
            return Err(LatticeQmcError::InvalidModel(
                "a vertex needs two legs per local site".into(),
            ));
        }
        if !weight.is_finite() || weight <= 0.0 {
            return Err(LatticeQmcError::InvalidModel(format!(
                "vertex weight must be finite and positive, got {weight}"
            )));
        }
        let diagonal = legs.chunks_exact(2).all(|pair| pair[0] == pair[1]);
        Ok(Self {
            name: name.into(),
            legs: legs.into_boxed_slice(),
            weight,
            diagonal,
        })
    }

    /// Human-readable matrix-element name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Local leg states.
    pub fn legs(&self) -> &[BasisState] {
        &self.legs
    }

    /// Positive matrix-element weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }

    /// Whether every incoming state equals its outgoing state.
    pub fn is_diagonal(&self) -> bool {
        self.diagonal
    }

    /// State on one leg.
    pub fn state(&self, leg: usize) -> BasisState {
        self.legs[leg]
    }
}

/// One sampled operator insertion.
#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    /// Imaginary time in `[0,beta)`.
    pub tau: f64,
    /// Operator-term index in the model catalog.
    pub term: usize,
    /// Matrix-element kind inside the term.
    pub kind: usize,
}

/// One site endpoint in a time-ordered worldline index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Event {
    /// Imaginary time.
    pub tau: f64,
    /// Vertex index.
    pub vertex: usize,
    /// Position of the site inside the local operator term.
    pub local_site: usize,
    /// Incoming global leg.
    pub incoming_leg: usize,
    /// Outgoing global leg.
    pub outgoing_leg: usize,
}
