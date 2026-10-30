use core::fmt;
use std::collections::TryReserveError;

use crate::UndirectedGraphView;

/// Cached original-graph degree data (N2 `NetworkProfile`).
///
/// Degrees are counted from physical edge endpoints once at construction and
/// are never updated by later activity: every degree-conditioned observable
/// is conditioned on **original** degree, never residual degree after
/// removals. A self-loop contributes two incidences at its vertex, matching
/// the two directed incidences it contributes to the CSR representation.
///
/// `bins()[k]` is the number of vertices with original degree `k`;
/// `bins().len() == max_degree() + 1`, so `k = 0 ..= max_degree()` covers
/// every bin including degrees absent from the graph (`bins()[k] == 0`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NetworkProfile {
    degrees: Vec<u32>,
    bins: Vec<usize>,
}

impl NetworkProfile {
    /// Count original degrees from `graph`'s physical edge table.
    pub fn from_graph(graph: &(impl UndirectedGraphView + ?Sized)) -> Result<Self, NetworkError> {
        let vertex_count = graph.vertex_count();
        let mut degrees = vec![0u32; vertex_count];
        for edge_index in 0..graph.edge_count() {
            let edge = graph
                .edge_id(edge_index)
                .ok_or(NetworkError::DenseEdgeId { edge: edge_index })?;
            for endpoint in graph.edge_endpoints(edge) {
                degrees[endpoint.index()] += 1;
            }
        }
        let max_degree = degrees.iter().copied().max().unwrap_or(0) as usize;
        let mut bins = vec![0usize; max_degree + 1];
        for &degree in &degrees {
            bins[degree as usize] += 1;
        }
        Ok(Self { degrees, bins })
    }

    /// Original degree of the vertex at dense index `index`.
    ///
    /// # Panics
    /// Panics when `index >= vertex_count()`.
    pub fn degree(&self, index: usize) -> u32 {
        self.degrees[index]
    }

    /// Original degree of every vertex, indexed by dense vertex index.
    pub fn degrees(&self) -> &[u32] {
        &self.degrees
    }

    /// Vertices per original degree value; `bins()[k] == 0` marks a degree
    /// value absent from the graph (an empty degree bin).
    pub fn bins(&self) -> &[usize] {
        &self.bins
    }

    /// Largest original degree; `bins().len() - 1`.
    pub fn max_degree(&self) -> usize {
        self.bins.len() - 1
    }

    /// Number of vertices with original degree `k` (`0` beyond the bins).
    pub fn vertices_with_degree(&self, degree: usize) -> usize {
        self.bins.get(degree).copied().unwrap_or(0)
    }

    /// Mean and second-over-first moment ratio of the original degree
    /// distribution: `(<k>, kappa)` with `kappa = <k^2>/<k>`.
    ///
    /// `kappa` is the quantity in Cohen et al., PRL 85, 4626 (2000):
    /// random-failure removal breaks the giant component at
    /// `f_c = 1 - 1/(kappa - 1)`.
    pub fn degree_moments(&self) -> Option<(f64, f64)> {
        let total: usize = self.bins.iter().sum();
        if total == 0 {
            return None;
        }
        let (mut first, mut second) = (0.0_f64, 0.0_f64);
        for (degree, &count) in self.bins.iter().enumerate() {
            let degree = degree as f64;
            first += degree * count as f64;
            second += degree * degree * count as f64;
        }
        first /= total as f64;
        second /= total as f64;
        Some((first, second / first))
    }

    /// Total vertex count the profile was built for.
    pub fn vertex_count(&self) -> usize {
        self.degrees.len()
    }
}

#[derive(Debug)]
pub enum NetworkError {
    /// The graph returned no dense ID for a physical edge index.
    DenseEdgeId { edge: usize },
    /// Degree storage could not be reserved.
    Capacity(TryReserveError),
}

impl From<TryReserveError> for NetworkError {
    fn from(source: TryReserveError) -> Self {
        Self::Capacity(source)
    }
}

impl fmt::Display for NetworkError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DenseEdgeId { edge } => {
                write!(
                    formatter,
                    "graph returned no dense ID for edge index {edge}"
                )
            }
            Self::Capacity(_) => write!(formatter, "could not reserve degree storage"),
        }
    }
}

impl std::error::Error for NetworkError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            _ => None,
        }
    }
}
