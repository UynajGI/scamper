//! Network input (N track): text edge lists to owned CSR graphs.
//!
//! The roadmap's network-input deliverable: an edge-list parser that builds
//! an owned [`CsrLattice`] with **explicit** parse errors — offending line
//! number, endpoint range violations, and the configured duplicate and
//! self-loop policies — plus the already-validated zero-copy path for
//! complete external CSR arrays ([`crate::BorrowedUndirectedCsr`]).
//! No binary or memory-mapped format is introduced: million-edge ingestion
//! benchmarks (see `benches/percolation_network_input.rs`) must justify one
//! first.
//!
//! Format: one edge per line as two whitespace-separated `usize` endpoints
//! (`left right`, undirected; order within a line is irrelevant). Empty
//! lines and lines whose first non-whitespace character is `#` are skipped.
//! Any other content (extra tokens, non-integers, trailing junk) is a typed
//! parse error carrying the 1-based line number. `vertex_count` must be at
//! least 1. CRLF input parses (line splitting and trimming both drop `\r`);
//! a UTF-8 BOM is not special-cased and yields a typed parse error.

use std::collections::HashSet;

use crate::{Bond, BondType, CsrLattice};

/// Policy for repeated edges in an edge list.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DuplicateEdgePolicy {
    /// Reject the input: the first repeated edge is a typed error.
    Reject,
    /// Keep every repetition as a distinct physical edge (parallel edges
    /// are a valid multigraph representation and the F1/F3 layers support
    /// them natively).
    AllowParallel,
}

/// Policy for `left == right` lines.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SelfLoopPolicy {
    /// Reject the input: the first self-loop is a typed error.
    Reject,
    /// Keep self-loops (they contribute two incidences at their vertex,
    /// matching the CSR convention; the analyzer treats them as union
    /// no-ops).
    Allow,
}

/// Explicit input contract for [`csr_from_edge_list`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct EdgeListSpec {
    /// Declared vertex count; every endpoint must lie in `0..vertex_count`.
    pub vertex_count: usize,
    pub duplicate_policy: DuplicateEdgePolicy,
    pub self_loop_policy: SelfLoopPolicy,
}

impl EdgeListSpec {
    /// Strict defaults: `vertex_count` vertices, no duplicates, no
    /// self-loops — the contract of a simple-graph edge list.
    pub const fn simple(vertex_count: usize) -> Self {
        Self {
            vertex_count,
            duplicate_policy: DuplicateEdgePolicy::Reject,
            self_loop_policy: SelfLoopPolicy::Reject,
        }
    }

    /// Multigraph contract: parallel edges and self-loops allowed.
    pub const fn multigraph(vertex_count: usize) -> Self {
        Self {
            vertex_count,
            duplicate_policy: DuplicateEdgePolicy::AllowParallel,
            self_loop_policy: SelfLoopPolicy::Allow,
        }
    }
}

/// Typed edge-list rejection, always carrying the offending line (1-based)
/// or a construction-level reason.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NetworkInputError {
    /// The spec declared zero vertices: no graph can be built. (Review
    /// follow-up 2026-10-22: this used to fall through to a construction
    /// panic.)
    InvalidVertexCount { vertex_count: usize },
    /// A line is not `left right` with two integer endpoints.
    Parse { line: usize, reason: &'static str },
    /// An endpoint lies outside `0..vertex_count`.
    Endpoint {
        line: usize,
        value: usize,
        vertex_count: usize,
    },
    /// A repeated edge under [`DuplicateEdgePolicy::Reject`].
    DuplicateEdge {
        line: usize,
        left: usize,
        right: usize,
    },
    /// A `left == right` line under [`SelfLoopPolicy::Reject`].
    SelfLoop { line: usize, vertex: usize },
}

impl core::fmt::Display for NetworkInputError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidVertexCount { vertex_count } => {
                write!(
                    formatter,
                    "vertex count must be at least 1, got {vertex_count}"
                )
            }
            Self::Parse { line, reason } => {
                write!(formatter, "line {line}: {reason}")
            }
            Self::Endpoint {
                line,
                value,
                vertex_count,
            } => write!(
                formatter,
                "line {line}: endpoint {value} outside 0..{vertex_count}"
            ),
            Self::DuplicateEdge { line, left, right } => {
                write!(formatter, "line {line}: repeated edge {left}-{right}")
            }
            Self::SelfLoop { line, vertex } => {
                write!(formatter, "line {line}: self-loop at vertex {vertex}")
            }
        }
    }
}

impl std::error::Error for NetworkInputError {}

/// Parse a text edge list under `spec` and build an owned CSR lattice.
///
/// Physical edges are emitted in file order (so CSR edge IDs follow input
/// order); each undirected pair is normalized to `(min, max)` internally
/// for duplicate detection only — the emitted bond keeps normalized
/// endpoints. Validation is line-local and typed: nothing is silently
/// dropped or merged.
pub fn csr_from_edge_list(spec: EdgeListSpec, text: &str) -> Result<CsrLattice, NetworkInputError> {
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    let mut edges: Vec<Bond> = Vec::new();
    for (index, raw_line) in text.lines().enumerate() {
        let line = index + 1;
        let trimmed = raw_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut tokens = trimmed.split_whitespace();
        let (Some(left), Some(right), None) = (tokens.next(), tokens.next(), tokens.next()) else {
            return Err(NetworkInputError::Parse {
                line,
                reason: "expected exactly two whitespace-separated integer endpoints",
            });
        };
        let parse = |token: &str| -> Result<usize, NetworkInputError> {
            token
                .parse::<usize>()
                .map_err(|_| NetworkInputError::Parse {
                    line,
                    reason: "endpoints must be non-negative decimal integers",
                })
        };
        let left = parse(left)?;
        let right = parse(right)?;
        for &endpoint in &[left, right] {
            if endpoint >= spec.vertex_count {
                return Err(NetworkInputError::Endpoint {
                    line,
                    value: endpoint,
                    vertex_count: spec.vertex_count,
                });
            }
        }
        if left == right && spec.self_loop_policy == SelfLoopPolicy::Reject {
            return Err(NetworkInputError::SelfLoop { line, vertex: left });
        }
        // Duplicate detection covers every edge, self-loops included: a
        // repeated edge is a repeated edge. (Review follow-up 2026-10-22:
        // self-loops used to be silently exempt.)
        if spec.duplicate_policy == DuplicateEdgePolicy::Reject {
            let pair = (left.min(right), left.max(right));
            if !seen.insert(pair) {
                return Err(NetworkInputError::DuplicateEdge {
                    line,
                    left: pair.0,
                    right: pair.1,
                });
            }
        }
        edges.push(Bond::new(
            left.min(right),
            left.max(right),
            BondType::Generic,
            1.0,
        ));
    }
    Ok(CsrLattice::from_edges(spec.vertex_count, edges))
}
