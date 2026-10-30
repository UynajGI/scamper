//! Versioned network I/O (X1): weighted edge lists and multiplex documents.
//!
//! Exactly two formats are supported, both explicit about their version:
//!
//! - **Edge-list v1** (text): the first line must be
//!   `# cmc-edge-list v1`. Required directives are `# nodes N` and
//!   `# layer NAME` (exactly once each); other `#` lines are comments.
//!   Each remaining line is `source target weight`, whitespace-separated.
//!   The graph is undirected; a line order within a pair is irrelevant for
//!   topology conversion but preserved for round-tripping.
//! - **Network JSON v1**: a document with `format: "cmc.network"`,
//!   `version: 1`, a node count, named layers, and weighted edges with an
//!   explicit, strictly increasing layer-index list. Unknown fields and
//!   unknown versions are typed rejections.
//!
//! Nothing is coerced silently: floats are never truncated to integers,
//! weights must be finite and non-negative, layer indices must be in
//! range and strictly increasing, and `to_csr_layer` requires the caller's
//! [`EdgeListSpec`] vertex count to equal the document node count. Duplicate
//! and self-loop policies are applied only when building a CSR topology, so
//! a document can faithfully represent a multigraph.

use std::collections::HashSet;

use serde::{Deserialize, Serialize};

use super::{DuplicateEdgePolicy, EdgeListSpec, NetworkInputError, SelfLoopPolicy};
use crate::{Bond, BondType, CsrLattice};

/// JSON `format` discriminator for [`NetworkDocument`].
pub const NETWORK_JSON_FORMAT: &str = "cmc.network";
/// JSON schema version understood by this crate.
pub const NETWORK_JSON_VERSION: u32 = 1;
/// Mandatory first line of the text edge-list format.
pub const EDGE_LIST_HEADER_V1: &str = "# cmc-edge-list v1";

/// One weighted, layer-tagged undirected edge.
///
/// `layers` is a strictly increasing list of indices into
/// [`NetworkDocument::layers`]; an edge belongs to every listed layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkEdge {
    pub source: usize,
    pub target: usize,
    pub weight: f64,
    pub layers: Vec<usize>,
}

impl NetworkEdge {
    /// Build an edge, validating weight and layer-list shape (not ranges).
    pub fn new(
        source: usize,
        target: usize,
        weight: f64,
        layers: impl Into<Vec<usize>>,
    ) -> Result<Self, NetworkIoError> {
        let edge = Self {
            source,
            target,
            weight,
            layers: layers.into(),
        };
        if !weight.is_finite() || weight < 0.0 {
            return Err(NetworkIoError::EdgeWeight {
                edge: 0,
                value: weight,
            });
        }
        if edge.layers.is_empty() {
            return Err(NetworkIoError::EmptyEdgeLayers { edge: 0 });
        }
        if !edge.layers.windows(2).all(|pair| pair[0] < pair[1]) {
            return Err(NetworkIoError::EdgeLayerOrder { edge: 0 });
        }
        Ok(edge)
    }
}

/// Validated multiplex network document (the in-memory form of both formats).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NetworkDocument {
    format: String,
    version: u32,
    nodes: usize,
    layers: Vec<String>,
    edges: Vec<NetworkEdge>,
}

impl NetworkDocument {
    /// A document with `nodes` dense vertex IDs and named layers.
    pub fn new(
        nodes: usize,
        layers: impl IntoIterator<Item = impl Into<String>>,
    ) -> Result<Self, NetworkIoError> {
        let document = Self {
            format: NETWORK_JSON_FORMAT.to_owned(),
            version: NETWORK_JSON_VERSION,
            nodes,
            layers: layers.into_iter().map(Into::into).collect(),
            edges: Vec::new(),
        };
        document.validate()?;
        Ok(document)
    }

    /// Append an edge after checking weight, layer-list shape, and all
    /// document-relative ranges.
    pub fn push_edge(&mut self, edge: NetworkEdge) -> Result<(), NetworkIoError> {
        let index = self.edges.len();
        self.check_edge(index, &edge)?;
        self.edges.push(edge);
        Ok(())
    }

    pub fn nodes(&self) -> usize {
        self.nodes
    }

    pub fn layers(&self) -> &[String] {
        &self.layers
    }

    pub fn edges(&self) -> &[NetworkEdge] {
        &self.edges
    }

    /// Re-check every invariant; [`Self::from_json`] and deserialization
    /// paths call this before returning a document.
    pub fn validate(&self) -> Result<(), NetworkIoError> {
        if self.format != NETWORK_JSON_FORMAT {
            return Err(NetworkIoError::Format {
                expected: NETWORK_JSON_FORMAT,
                actual: self.format.clone(),
            });
        }
        if self.version != NETWORK_JSON_VERSION {
            return Err(NetworkIoError::Version {
                supported: NETWORK_JSON_VERSION,
                actual: self.version,
            });
        }
        if self.nodes == 0 {
            return Err(NetworkIoError::NodeCount { nodes: 0 });
        }
        if self.layers.is_empty() {
            return Err(NetworkIoError::NoLayers);
        }
        let mut names = HashSet::with_capacity(self.layers.len());
        for (index, name) in self.layers.iter().enumerate() {
            if name.is_empty() {
                return Err(NetworkIoError::EmptyLayerName { layer: index });
            }
            if !names.insert(name.as_str()) {
                return Err(NetworkIoError::DuplicateLayerName { name: name.clone() });
            }
        }
        for (index, edge) in self.edges.iter().enumerate() {
            self.check_edge(index, edge)?;
        }
        Ok(())
    }

    fn check_edge(&self, index: usize, edge: &NetworkEdge) -> Result<(), NetworkIoError> {
        for &endpoint in &[edge.source, edge.target] {
            if endpoint >= self.nodes {
                return Err(NetworkIoError::EdgeEndpoint {
                    edge: index,
                    endpoint,
                    nodes: self.nodes,
                });
            }
        }
        if !edge.weight.is_finite() || edge.weight < 0.0 {
            return Err(NetworkIoError::EdgeWeight {
                edge: index,
                value: edge.weight,
            });
        }
        if edge.layers.is_empty() {
            return Err(NetworkIoError::EmptyEdgeLayers { edge: index });
        }
        if !edge.layers.windows(2).all(|pair| pair[0] < pair[1]) {
            return Err(NetworkIoError::EdgeLayerOrder { edge: index });
        }
        for &layer in &edge.layers {
            if layer >= self.layers.len() {
                return Err(NetworkIoError::EdgeLayer {
                    edge: index,
                    layer,
                    layers: self.layers.len(),
                });
            }
        }
        Ok(())
    }

    /// Canonical JSON for this document.
    pub fn to_json(&self) -> Result<String, NetworkIoError> {
        self.validate()?;
        serde_json::to_string_pretty(self)
            .map_err(|source| NetworkIoError::Json(source.to_string()))
    }

    /// Parse and fully validate a JSON document.
    pub fn from_json(text: &str) -> Result<Self, NetworkIoError> {
        let document: Self = serde_json::from_str(text)
            .map_err(|source| NetworkIoError::Json(source.to_string()))?;
        document.validate()?;
        Ok(document)
    }

    /// Serialize a single-layer document as edge-list v1 text.
    pub fn to_edge_list_v1(&self) -> Result<String, NetworkIoError> {
        self.validate()?;
        if self.layers.len() != 1 {
            return Err(NetworkIoError::SingleLayerOnly {
                layers: self.layers.len(),
            });
        }
        let mut text = String::new();
        text.push_str(EDGE_LIST_HEADER_V1);
        text.push('\n');
        text.push_str(&format!("# nodes {}\n", self.nodes));
        text.push_str(&format!("# layer {}\n", self.layers[0]));
        for edge in &self.edges {
            // Debug formatting is the shortest round-tripping `f64` form.
            text.push_str(&format!(
                "{} {} {:?}\n",
                edge.source, edge.target, edge.weight
            ));
        }
        Ok(text)
    }

    /// Parse edge-list v1 text into a single-layer document.
    pub fn from_edge_list_v1(text: &str) -> Result<Self, NetworkIoError> {
        let mut lines = text.lines();
        let header = lines.next().ok_or(NetworkIoError::EdgeListHeader {
            line: 1,
            reason: "missing `# cmc-edge-list v1` header",
        })?;
        if header.trim() != EDGE_LIST_HEADER_V1 {
            return Err(NetworkIoError::EdgeListHeader {
                line: 1,
                reason: "first line must be `# cmc-edge-list v1`",
            });
        }
        let mut nodes = None;
        let mut layer = None;
        let mut edges = Vec::new();
        for (index, raw_line) in lines.enumerate() {
            let line = index + 2;
            let trimmed = raw_line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(directive) = trimmed.strip_prefix('#') {
                let tokens: Vec<&str> = directive.split_whitespace().collect();
                match tokens.as_slice() {
                    ["nodes", ..] if tokens.len() != 2 => {
                        return Err(NetworkIoError::EdgeListLine {
                            line,
                            reason: "`# nodes` expects exactly one integer",
                        });
                    }
                    ["layer", ..] if tokens.len() != 2 => {
                        return Err(NetworkIoError::EdgeListLine {
                            line,
                            reason: "`# layer` expects exactly one name",
                        });
                    }
                    ["nodes", value] => {
                        if nodes.is_some() {
                            return Err(NetworkIoError::DuplicateDirective {
                                line,
                                name: "nodes",
                            });
                        }
                        nodes = Some(value.parse::<usize>().map_err(|_| {
                            NetworkIoError::EdgeListLine {
                                line,
                                reason: "`# nodes` expects a non-negative integer",
                            }
                        })?);
                    }
                    ["layer", name] => {
                        if layer.is_some() {
                            return Err(NetworkIoError::DuplicateDirective {
                                line,
                                name: "layer",
                            });
                        }
                        layer = Some((*name).to_owned());
                    }
                    _ => {
                        // Plain comment line.
                    }
                }
                continue;
            }
            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            let [source, target, weight] = tokens.as_slice() else {
                return Err(NetworkIoError::EdgeListLine {
                    line,
                    reason: "expected `source target weight`",
                });
            };
            let parse_index = |token: &str| {
                token
                    .parse::<usize>()
                    .map_err(|_| NetworkIoError::EdgeListLine {
                        line,
                        reason: "endpoints must be non-negative decimal integers",
                    })
            };
            let weight = weight
                .parse::<f64>()
                .map_err(|_| NetworkIoError::EdgeListLine {
                    line,
                    reason: "weight must be a decimal number",
                })?;
            edges.push(NetworkEdge {
                source: parse_index(source)?,
                target: parse_index(target)?,
                weight,
                layers: vec![0],
            });
        }
        let nodes = nodes.ok_or(NetworkIoError::MissingDirective { name: "nodes" })?;
        let layer = layer.ok_or(NetworkIoError::MissingDirective { name: "layer" })?;
        let document = Self {
            format: NETWORK_JSON_FORMAT.to_owned(),
            version: NETWORK_JSON_VERSION,
            nodes,
            layers: vec![layer],
            edges,
        };
        document.validate()?;
        Ok(document)
    }

    /// Project one layer into an owned CSR topology.
    ///
    /// The spec's vertex count must equal [`Self::nodes`]; duplicate and
    /// self-loop policies apply exactly as in [`super::csr_from_edge_list`].
    /// Edge weights are carried into the [`CsrLattice::edges`] bonds.
    pub fn to_csr_layer(
        &self,
        layer: usize,
        spec: EdgeListSpec,
    ) -> Result<CsrLattice, NetworkIoError> {
        self.validate()?;
        if layer >= self.layers.len() {
            return Err(NetworkIoError::LayerIndex {
                layer,
                layers: self.layers.len(),
            });
        }
        if spec.vertex_count != self.nodes {
            return Err(NetworkIoError::NodeCountMismatch {
                document: self.nodes,
                spec: spec.vertex_count,
            });
        }
        let mut seen = HashSet::new();
        let mut bonds = Vec::new();
        for (index, edge) in self.edges.iter().enumerate() {
            if !edge.layers.contains(&layer) {
                continue;
            }
            let line = index + 1;
            if edge.source == edge.target && spec.self_loop_policy == SelfLoopPolicy::Reject {
                return Err(NetworkIoError::Csr(NetworkInputError::SelfLoop {
                    line,
                    vertex: edge.source,
                }));
            }
            let pair = (edge.source.min(edge.target), edge.source.max(edge.target));
            if spec.duplicate_policy == DuplicateEdgePolicy::Reject && !seen.insert(pair) {
                return Err(NetworkIoError::Csr(NetworkInputError::DuplicateEdge {
                    line,
                    left: pair.0,
                    right: pair.1,
                }));
            }
            bonds.push(Bond::new(pair.0, pair.1, BondType::Generic, edge.weight));
        }
        Ok(CsrLattice::from_edges(spec.vertex_count, bonds))
    }
}

/// Typed network-I/O rejection; nothing is silently repaired.
#[derive(Clone, Debug, PartialEq)]
pub enum NetworkIoError {
    Format {
        expected: &'static str,
        actual: String,
    },
    Version {
        supported: u32,
        actual: u32,
    },
    NodeCount {
        nodes: usize,
    },
    NoLayers,
    EmptyLayerName {
        layer: usize,
    },
    DuplicateLayerName {
        name: String,
    },
    EdgeEndpoint {
        edge: usize,
        endpoint: usize,
        nodes: usize,
    },
    EdgeWeight {
        edge: usize,
        value: f64,
    },
    EmptyEdgeLayers {
        edge: usize,
    },
    EdgeLayerOrder {
        edge: usize,
    },
    EdgeLayer {
        edge: usize,
        layer: usize,
        layers: usize,
    },
    LayerIndex {
        layer: usize,
        layers: usize,
    },
    NodeCountMismatch {
        document: usize,
        spec: usize,
    },
    SingleLayerOnly {
        layers: usize,
    },
    MissingDirective {
        name: &'static str,
    },
    DuplicateDirective {
        line: usize,
        name: &'static str,
    },
    EdgeListHeader {
        line: usize,
        reason: &'static str,
    },
    EdgeListLine {
        line: usize,
        reason: &'static str,
    },
    Json(String),
    Csr(NetworkInputError),
}

impl core::fmt::Display for NetworkIoError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Format { expected, actual } => {
                write!(
                    formatter,
                    "unsupported network format {actual:?} (expected {expected:?})"
                )
            }
            Self::Version { supported, actual } => {
                write!(
                    formatter,
                    "unsupported network schema version {actual} (supported {supported})"
                )
            }
            Self::NodeCount { nodes } => {
                write!(formatter, "node count must be at least 1, got {nodes}")
            }
            Self::NoLayers => write!(formatter, "a network document needs at least one layer"),
            Self::EmptyLayerName { layer } => write!(formatter, "layer {layer} has an empty name"),
            Self::DuplicateLayerName { name } => write!(formatter, "duplicate layer name {name:?}"),
            Self::EdgeEndpoint {
                edge,
                endpoint,
                nodes,
            } => write!(
                formatter,
                "edge {edge}: endpoint {endpoint} outside 0..{nodes}"
            ),
            Self::EdgeWeight { edge, value } => {
                write!(
                    formatter,
                    "edge {edge}: weight {value} must be finite and non-negative"
                )
            }
            Self::EmptyEdgeLayers { edge } => write!(formatter, "edge {edge} belongs to no layer"),
            Self::EdgeLayerOrder { edge } => {
                write!(
                    formatter,
                    "edge {edge}: layer indices must be strictly increasing"
                )
            }
            Self::EdgeLayer {
                edge,
                layer,
                layers,
            } => write!(formatter, "edge {edge}: layer {layer} outside 0..{layers}"),
            Self::LayerIndex { layer, layers } => {
                write!(formatter, "layer index {layer} outside 0..{layers}")
            }
            Self::NodeCountMismatch { document, spec } => write!(
                formatter,
                "edge-list spec declares {spec} vertices but the document has {document}"
            ),
            Self::SingleLayerOnly { layers } => write!(
                formatter,
                "edge-list v1 supports one layer, document has {layers}"
            ),
            Self::MissingDirective { name } => {
                write!(
                    formatter,
                    "edge-list v1 is missing the `# {name}` directive"
                )
            }
            Self::DuplicateDirective { line, name } => {
                write!(formatter, "line {line}: duplicate `# {name}` directive")
            }
            Self::EdgeListHeader { line, reason } => write!(formatter, "line {line}: {reason}"),
            Self::EdgeListLine { line, reason } => write!(formatter, "line {line}: {reason}"),
            Self::Json(reason) => write!(formatter, "invalid network JSON: {reason}"),
            Self::Csr(source) => write!(formatter, "cannot build topology: {source}"),
        }
    }
}

impl std::error::Error for NetworkIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Csr(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GraphView;

    fn document() -> NetworkDocument {
        let mut document = NetworkDocument::new(4, ["main", "aux"]).expect("document");
        document
            .push_edge(NetworkEdge::new(0, 1, 1.0, [0]).expect("edge"))
            .expect("edge");
        document
            .push_edge(NetworkEdge::new(1, 2, 2.5, [0, 1]).expect("edge"))
            .expect("edge");
        document
            .push_edge(NetworkEdge::new(2, 3, 0.0, [1]).expect("edge"))
            .expect("edge");
        document
    }

    #[test]
    fn json_round_trip_preserves_the_document() {
        let document = document();
        let text = document.to_json().expect("serialize");
        assert_eq!(NetworkDocument::from_json(&text).expect("parse"), document);
    }

    #[test]
    fn edge_list_round_trip_preserves_single_layer_documents() {
        let mut document = NetworkDocument::new(3, ["main"]).expect("document");
        for (source, target, weight) in [(0, 1, 0.1), (1, 2, 2.5), (0, 2, 0.0)] {
            document
                .push_edge(NetworkEdge::new(source, target, weight, [0]).expect("edge"))
                .expect("edge");
        }
        let text = document.to_edge_list_v1().expect("serialize");
        assert_eq!(
            NetworkDocument::from_edge_list_v1(&text).expect("parse"),
            document
        );
    }

    #[test]
    fn json_rejections_are_typed_and_versioned() {
        let wrong_format = r#"{"format":"other","version":1,"nodes":2,"layers":["a"],"edges":[]}"#;
        assert!(matches!(
            NetworkDocument::from_json(wrong_format),
            Err(NetworkIoError::Format { .. })
        ));
        let wrong_version =
            r#"{"format":"cmc.network","version":2,"nodes":2,"layers":["a"],"edges":[]}"#;
        assert!(matches!(
            NetworkDocument::from_json(wrong_version),
            Err(NetworkIoError::Version { actual: 2, .. })
        ));
        let unknown_field =
            r#"{"format":"cmc.network","version":1,"nodes":2,"layers":["a"],"edges":[],"extra":1}"#;
        assert!(matches!(
            NetworkDocument::from_json(unknown_field),
            Err(NetworkIoError::Json(_))
        ));
        let out_of_range = r#"{"format":"cmc.network","version":1,"nodes":2,"layers":["a"],"edges":[{"source":0,"target":2,"weight":1.0,"layers":[0]}]}"#;
        assert!(matches!(
            NetworkDocument::from_json(out_of_range),
            Err(NetworkIoError::EdgeEndpoint { endpoint: 2, .. })
        ));
        let unsorted = r#"{"format":"cmc.network","version":1,"nodes":2,"layers":["a","b"],"edges":[{"source":0,"target":1,"weight":1.0,"layers":[1,0]}]}"#;
        assert!(matches!(
            NetworkDocument::from_json(unsorted),
            Err(NetworkIoError::EdgeLayerOrder { .. })
        ));
    }

    #[test]
    fn edge_list_rejections_are_typed_and_line_local() {
        assert!(matches!(
            NetworkDocument::from_edge_list_v1("# cmc-edge-list v2\n# nodes 2\n# layer a\n"),
            Err(NetworkIoError::EdgeListHeader { line: 1, .. })
        ));
        assert!(matches!(
            NetworkDocument::from_edge_list_v1("# cmc-edge-list v1\n# layer a\n"),
            Err(NetworkIoError::MissingDirective { name: "nodes" })
        ));
        assert!(matches!(
            NetworkDocument::from_edge_list_v1("# cmc-edge-list v1\n# nodes 2\n# layer a\n0 1\n"),
            Err(NetworkIoError::EdgeListLine { line: 4, .. })
        ));
        assert!(matches!(
            NetworkDocument::from_edge_list_v1(
                "# cmc-edge-list v1\n# nodes 2\n# layer a\n0 1 -1.0\n"
            ),
            Err(NetworkIoError::EdgeWeight { .. })
        ));
        assert!(matches!(
            NetworkDocument::from_edge_list_v1(
                "# cmc-edge-list v1\n# nodes 2\n# nodes 3\n# layer a\n"
            ),
            Err(NetworkIoError::DuplicateDirective {
                line: 3,
                name: "nodes"
            })
        ));
    }

    #[test]
    fn layer_projection_carries_weights_and_applies_policies() {
        let document = document();
        let lattice = document
            .to_csr_layer(0, EdgeListSpec::simple(4))
            .expect("layer main");
        assert_eq!(lattice.vertex_count(), 4);
        assert_eq!(lattice.edges.len(), 2);
        assert_eq!(lattice.edges[0].weight, 1.0);
        assert_eq!(lattice.edges[1].weight, 2.5);
        // Edge 2-3 lives only in the aux layer.
        assert!(lattice
            .edges
            .iter()
            .all(|bond| !(bond.source == 2 && bond.target == 3)));

        assert!(matches!(
            document.to_csr_layer(0, EdgeListSpec::simple(3)),
            Err(NetworkIoError::NodeCountMismatch {
                document: 4,
                spec: 3
            })
        ));
        assert!(matches!(
            document.to_csr_layer(2, EdgeListSpec::simple(4)),
            Err(NetworkIoError::LayerIndex {
                layer: 2,
                layers: 2
            })
        ));

        let mut duplicate = NetworkDocument::new(2, ["main"]).expect("document");
        duplicate
            .push_edge(NetworkEdge::new(0, 1, 1.0, [0]).expect("edge"))
            .expect("edge");
        duplicate
            .push_edge(NetworkEdge::new(1, 0, 1.0, [0]).expect("edge"))
            .expect("edge");
        assert!(matches!(
            duplicate.to_csr_layer(0, EdgeListSpec::simple(2)),
            Err(NetworkIoError::Csr(NetworkInputError::DuplicateEdge { .. }))
        ));
        assert!(duplicate
            .to_csr_layer(0, EdgeListSpec::multigraph(2))
            .is_ok());
    }

    #[test]
    fn layer_projection_matches_the_reference_edge_list_graph() {
        let reference =
            super::super::csr_from_edge_list(EdgeListSpec::simple(4), "0 1\n1 2\n2 3\n0 3\n")
                .expect("reference graph");
        let document = NetworkDocument::from_edge_list_v1(
            "# cmc-edge-list v1\n# nodes 4\n# layer main\n0 1 1.0\n1 2 1.0\n2 3 1.0\n0 3 1.0\n",
        )
        .expect("document");
        let projected = document
            .to_csr_layer(0, EdgeListSpec::simple(4))
            .expect("projection");
        assert_eq!(projected.edges, reference.edges);
    }

    #[test]
    fn multi_layer_documents_reject_edge_list_serialization() {
        assert!(matches!(
            document().to_edge_list_v1(),
            Err(NetworkIoError::SingleLayerOnly { layers: 2 })
        ));
    }
}
