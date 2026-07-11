//! Graph topology for lattice QMC.
//!
//! The quantum algorithms operate on an undirected weighted multigraph. A
//! square lattice, a hypercubic lattice in arbitrary dimension, a molecular
//! graph, and an irregular finite cluster are therefore the same runtime
//! object. The storage is CSR-like for cache-friendly site traversal while a
//! separate unique edge table supports bond-operator sampling.

use std::collections::HashSet;

use thiserror::Error;

/// Error raised while constructing or validating a graph.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum GraphError {
    /// A graph must contain at least one site.
    #[error("graph must contain at least one site")]
    Empty,
    /// An edge endpoint is outside the graph.
    #[error("edge {edge} endpoint {site} is outside 0..{n_sites}")]
    EndpointOutOfRange {
        /// Edge position in the input sequence.
        edge: usize,
        /// Invalid endpoint.
        site: usize,
        /// Number of sites.
        n_sites: usize,
    },
    /// Self loops are not supported by the current local-operator backend.
    #[error("self loop at site {site} is not supported")]
    SelfLoop {
        /// Site carrying the self loop.
        site: usize,
    },
    /// Exact duplicate edges are rejected to avoid accidental double counting.
    #[error("duplicate edge ({src}, {target}, kind={kind})")]
    DuplicateEdge {
        /// Canonical first endpoint.
        src: usize,
        /// Canonical second endpoint.
        target: usize,
        /// User-defined edge kind.
        kind: u16,
    },
    /// Edge weight is not finite.
    #[error("edge {edge} has non-finite weight {weight}")]
    InvalidWeight {
        /// Edge position in the input sequence.
        edge: usize,
        /// Invalid weight.
        weight: f64,
    },
    /// Hypercubic dimensions are invalid.
    #[error("all hypercubic dimensions must be positive")]
    InvalidDimensions,
    /// Raw CSR arrays do not describe valid adjacency rows.
    #[error("invalid CSR graph: {0}")]
    InvalidCsr(String),
    /// An adjacency list was not symmetric.
    #[error("adjacency is not symmetric: {src} lists {target}, but not vice versa")]
    AsymmetricAdjacency {
        /// Source site.
        src: usize,
        /// Missing reverse endpoint.
        target: usize,
    },
}

/// Construction-time edge specification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeSpec {
    /// First endpoint.
    pub source: usize,
    /// Second endpoint.
    pub target: usize,
    /// User-defined bond class. Algorithms do not assign geometry semantics.
    pub kind: u16,
    /// Multiplicative coupling scale attached to this edge.
    pub weight: f64,
}

impl EdgeSpec {
    /// Construct an untyped unit-weight edge.
    pub const fn new(source: usize, target: usize) -> Self {
        Self {
            source,
            target,
            kind: 0,
            weight: 1.0,
        }
    }

    /// Construct a typed weighted edge.
    pub const fn typed(source: usize, target: usize, kind: u16, weight: f64) -> Self {
        Self {
            source,
            target,
            kind,
            weight,
        }
    }
}

/// Canonically oriented unique undirected edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Edge {
    /// Smaller endpoint.
    pub source: usize,
    /// Larger endpoint.
    pub target: usize,
    /// User-defined bond class.
    pub kind: u16,
    /// Multiplicative coupling scale.
    pub weight: f64,
}

/// One CSR adjacency entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Neighbor {
    /// Neighboring site.
    pub site: usize,
    /// Index into [`CsrGraph::edges`].
    pub edge: usize,
}

/// Undirected graph stored as a CSR adjacency plus a unique edge table.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrGraph {
    offsets: Vec<usize>,
    neighbors: Vec<Neighbor>,
    edges: Vec<Edge>,
    n_sites: usize,
}

impl CsrGraph {
    /// Build a graph from unique undirected edge specifications.
    pub fn from_edges(
        n_sites: usize,
        specs: impl IntoIterator<Item = EdgeSpec>,
    ) -> Result<Self, GraphError> {
        if n_sites == 0 {
            return Err(GraphError::Empty);
        }
        let mut edges = Vec::new();
        let mut seen = HashSet::new();
        for (edge_id, spec) in specs.into_iter().enumerate() {
            for site in [spec.source, spec.target] {
                if site >= n_sites {
                    return Err(GraphError::EndpointOutOfRange {
                        edge: edge_id,
                        site,
                        n_sites,
                    });
                }
            }
            if spec.source == spec.target {
                return Err(GraphError::SelfLoop { site: spec.source });
            }
            if !spec.weight.is_finite() {
                return Err(GraphError::InvalidWeight {
                    edge: edge_id,
                    weight: spec.weight,
                });
            }
            let (source, target) = if spec.source < spec.target {
                (spec.source, spec.target)
            } else {
                (spec.target, spec.source)
            };
            if !seen.insert((source, target, spec.kind)) {
                return Err(GraphError::DuplicateEdge {
                    src: source,
                    target,
                    kind: spec.kind,
                });
            }
            edges.push(Edge {
                source,
                target,
                kind: spec.kind,
                weight: spec.weight,
            });
        }
        edges.sort_by_key(|edge| (edge.source, edge.target, edge.kind));

        let mut adjacency = vec![Vec::<Neighbor>::new(); n_sites];
        for (edge_id, edge) in edges.iter().enumerate() {
            adjacency[edge.source].push(Neighbor {
                site: edge.target,
                edge: edge_id,
            });
            adjacency[edge.target].push(Neighbor {
                site: edge.source,
                edge: edge_id,
            });
        }
        for row in &mut adjacency {
            row.sort_by_key(|neighbor| (neighbor.site, neighbor.edge));
        }

        let mut offsets = Vec::with_capacity(n_sites + 1);
        let mut neighbors = Vec::with_capacity(2 * edges.len());
        offsets.push(0);
        for row in adjacency {
            neighbors.extend(row);
            offsets.push(neighbors.len());
        }
        Ok(Self {
            offsets,
            neighbors,
            edges,
            n_sites,
        })
    }

    /// Build from a symmetric adjacency list.
    pub fn from_adjacency(adjacency: &[Vec<usize>]) -> Result<Self, GraphError> {
        if adjacency.is_empty() {
            return Err(GraphError::Empty);
        }
        for (source, row) in adjacency.iter().enumerate() {
            for &target in row {
                if target >= adjacency.len() {
                    return Err(GraphError::EndpointOutOfRange {
                        edge: source,
                        site: target,
                        n_sites: adjacency.len(),
                    });
                }
                if !adjacency[target].contains(&source) {
                    return Err(GraphError::AsymmetricAdjacency {
                        src: source,
                        target,
                    });
                }
            }
        }
        let specs = adjacency.iter().enumerate().flat_map(|(source, row)| {
            row.iter()
                .copied()
                .filter(move |&target| source < target)
                .map(move |target| EdgeSpec::new(source, target))
        });
        Self::from_edges(adjacency.len(), specs)
    }

    /// Build from raw CSR arrays compatible with CMC-style adjacency storage.
    ///
    /// `offsets` must start at zero, be nondecreasing, and end at
    /// `neighbors.len()`. The adjacency must be undirected and symmetric.
    pub fn from_csr(offsets: &[usize], neighbors: &[usize]) -> Result<Self, GraphError> {
        if offsets.len() < 2 {
            return Err(GraphError::Empty);
        }
        if offsets[0] != 0 {
            return Err(GraphError::InvalidCsr(
                "the first offset must be zero".into(),
            ));
        }
        if offsets.windows(2).any(|pair| pair[0] > pair[1]) {
            return Err(GraphError::InvalidCsr(
                "offsets must be nondecreasing".into(),
            ));
        }
        if offsets.last().copied() != Some(neighbors.len()) {
            return Err(GraphError::InvalidCsr(
                "the final offset must equal neighbors.len()".into(),
            ));
        }
        let adjacency: Vec<Vec<usize>> = offsets
            .windows(2)
            .map(|range| neighbors[range[0]..range[1]].to_vec())
            .collect();
        Self::from_adjacency(&adjacency)
    }

    /// Export plain symmetric adjacency rows.
    pub fn to_adjacency(&self) -> Vec<Vec<usize>> {
        (0..self.site_count())
            .map(|site| {
                self.neighbors(site)
                    .iter()
                    .map(|entry| entry.site)
                    .collect()
            })
            .collect()
    }

    /// One-dimensional chain or ring.
    pub fn chain(n_sites: usize, periodic: bool) -> Result<Self, GraphError> {
        if n_sites == 0 {
            return Err(GraphError::Empty);
        }
        let mut edges = Vec::new();
        for site in 0..n_sites.saturating_sub(1) {
            edges.push(EdgeSpec::typed(site, site + 1, 0, 1.0));
        }
        if periodic && n_sites == 2 {
            edges[0].weight = 2.0;
        } else if periodic && n_sites > 2 {
            edges.push(EdgeSpec::typed(n_sites - 1, 0, 0, 1.0));
        }
        Self::from_edges(n_sites, edges)
    }

    /// Hypercubic graph in any positive dimension.
    pub fn hypercubic(dimensions: &[usize], periodic: bool) -> Result<Self, GraphError> {
        if dimensions.is_empty() || dimensions.contains(&0) {
            return Err(GraphError::InvalidDimensions);
        }
        let n_sites = dimensions.iter().product();
        let mut strides = vec![1_usize; dimensions.len()];
        for axis in 1..dimensions.len() {
            strides[axis] = strides[axis - 1] * dimensions[axis - 1];
        }
        let mut specs = Vec::new();
        for site in 0..n_sites {
            for axis in 0..dimensions.len() {
                let coordinate = (site / strides[axis]) % dimensions[axis];
                if dimensions[axis] == 2 && periodic {
                    if coordinate == 0 {
                        specs.push(EdgeSpec::typed(
                            site,
                            site + strides[axis],
                            axis as u16,
                            2.0,
                        ));
                    }
                    continue;
                }
                let target_coordinate = if coordinate + 1 < dimensions[axis] {
                    Some(coordinate + 1)
                } else if periodic && dimensions[axis] > 2 {
                    Some(0)
                } else {
                    None
                };
                if let Some(target_coordinate) = target_coordinate {
                    let target = if target_coordinate > coordinate {
                        site + strides[axis]
                    } else {
                        site - coordinate * strides[axis]
                    };
                    specs.push(EdgeSpec::typed(site, target, axis as u16, 1.0));
                }
            }
        }
        Self::from_edges(n_sites, specs)
    }

    /// Two-dimensional square graph, a special case of [`Self::hypercubic`].
    pub fn square(width: usize, height: usize, periodic: bool) -> Result<Self, GraphError> {
        Self::hypercubic(&[width, height], periodic)
    }

    /// Number of sites.
    pub const fn site_count(&self) -> usize {
        self.n_sites
    }

    /// Number of unique undirected edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Unique edge table.
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Edge by index.
    pub fn edge(&self, edge: usize) -> Edge {
        self.edges[edge]
    }

    /// CSR neighbors of a site.
    pub fn neighbors(&self, site: usize) -> &[Neighbor] {
        &self.neighbors[self.offsets[site]..self.offsets[site + 1]]
    }

    /// Degree of a site.
    pub fn degree(&self, site: usize) -> usize {
        self.offsets[site + 1] - self.offsets[site]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arbitrary_adjacency_round_trips() {
        let graph = CsrGraph::from_adjacency(&[vec![1, 2], vec![0, 3], vec![0, 3], vec![1, 2]])
            .expect("graph");
        assert_eq!(graph.site_count(), 4);
        assert_eq!(graph.edge_count(), 4);
        assert_eq!(graph.degree(0), 2);
    }

    #[test]
    fn raw_csr_matches_adjacency_builder() {
        let graph = CsrGraph::from_csr(&[0, 2, 4, 6], &[1, 2, 0, 2, 0, 1]).expect("csr graph");
        assert_eq!(graph.edge_count(), 3);
        assert_eq!(
            graph.to_adjacency(),
            vec![vec![1, 2], vec![0, 2], vec![0, 1]]
        );
    }

    #[test]
    fn hypercubic_is_dimension_agnostic() {
        let graph = CsrGraph::hypercubic(&[3, 4, 2], false).expect("graph");
        assert_eq!(graph.site_count(), 24);
        assert_eq!(graph.edge_count(), 46);
    }

    #[test]
    fn rejects_asymmetric_adjacency() {
        let error = CsrGraph::from_adjacency(&[vec![1], Vec::new()]).expect_err("error");
        assert!(matches!(error, GraphError::AsymmetricAdjacency { .. }));
    }
}
