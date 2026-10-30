//! X3: hypergraph percolation components.
//!
//! Definitional anchor (no literature value is claimed): a hypergraph has a
//! dense vertex set and a list of hyperedges (vertex sets). Under Bernoulli
//! site and hyperedge occupation, two active vertices are connected when
//! they share an active hyperedge; components are the connected components
//! of that active incidence subgraph. This is the natural hypergraph
//! generalization of the F3 component analysis and reduces to bond
//! percolation when every hyperedge has size two.
//!
//! Canonical identity rule: every component is identified by its smallest
//! dense vertex ID (the same lowest-ID convention as [`crate::ComponentSummary`]);
//! equal-size largest components therefore resolve to the one containing the
//! lowest vertex. Hyperedges are sets: members must be strictly increasing
//! at construction, so repeated members are a typed rejection rather than a
//! silent de-duplication. Repeated hyperedges are allowed (a multi-hypergraph)
//! and each carries its own occupation bit.

use core::fmt;

/// A hypergraph over a dense vertex set.
///
/// ```
/// use cmc_rs::Hypergraph;
///
/// let graph = Hypergraph::new(4, [vec![0, 1, 2]]).expect("hypergraph");
/// let summary = graph
///     .analyze(&[true, true, true, false], &[true])
///     .expect("analysis");
/// assert_eq!(summary.largest_component_size, 3);
/// assert_eq!(summary.component_count, 1);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct Hypergraph {
    vertices: usize,
    edges: Vec<Vec<usize>>,
}

impl Hypergraph {
    /// Build a hypergraph; every hyperedge must be a non-empty, strictly
    /// increasing list of vertex IDs below `vertices`.
    pub fn new(
        vertices: usize,
        edges: impl IntoIterator<Item = Vec<usize>>,
    ) -> Result<Self, HypergraphError> {
        if vertices == 0 {
            return Err(HypergraphError::NoVertices);
        }
        let edges: Vec<Vec<usize>> = edges.into_iter().collect();
        for (index, edge) in edges.iter().enumerate() {
            if edge.is_empty() {
                return Err(HypergraphError::EmptyHyperedge { edge: index });
            }
            for &member in edge {
                if member >= vertices {
                    return Err(HypergraphError::Member {
                        edge: index,
                        member,
                        vertices,
                    });
                }
            }
            if !edge.windows(2).all(|pair| pair[0] < pair[1]) {
                return Err(HypergraphError::MemberOrder { edge: index });
            }
        }
        Ok(Self { vertices, edges })
    }

    pub const fn vertices(&self) -> usize {
        self.vertices
    }

    pub fn edges(&self) -> &[Vec<usize>] {
        &self.edges
    }

    /// Analyze one occupation realization.
    ///
    /// `vertices` and `edges` are dense activity arrays; both must match the
    /// hypergraph domains. The summary reports active counts, components,
    /// the largest component, its canonical lowest-ID root, and per-vertex
    /// labels (`None` for inactive vertices).
    pub fn analyze(
        &self,
        vertices: &[bool],
        edges: &[bool],
    ) -> Result<HypergraphSummary, HypergraphError> {
        if vertices.len() != self.vertices {
            return Err(HypergraphError::Activity {
                entity: "vertex",
                expected: self.vertices,
                actual: vertices.len(),
            });
        }
        if edges.len() != self.edges.len() {
            return Err(HypergraphError::Activity {
                entity: "hyperedge",
                expected: self.edges.len(),
                actual: edges.len(),
            });
        }
        let mut union = UnionFind::new(self.vertices);
        let mut active_vertex_count = 0;
        for (vertex, &active) in vertices.iter().enumerate() {
            if active {
                active_vertex_count += 1;
                union.enable(vertex);
            }
        }
        let mut active_edge_count = 0;
        for (members, &active) in self.edges.iter().zip(edges) {
            if !active {
                continue;
            }
            active_edge_count += 1;
            let mut first_active = None;
            for &member in members {
                if !vertices[member] {
                    continue;
                }
                match first_active {
                    Some(first) => union.union(first, member),
                    None => first_active = Some(member),
                }
            }
        }
        union.finish(vertices, active_vertex_count, active_edge_count)
    }
}

/// Component summary of one hypergraph realization.
#[derive(Clone, Debug, PartialEq)]
pub struct HypergraphSummary {
    pub active_vertex_count: usize,
    pub active_edge_count: usize,
    pub component_count: usize,
    pub largest_component_size: usize,
    /// Canonical lowest vertex ID of the largest component (ties: lowest ID).
    pub largest_component: Option<usize>,
    pub raw_second_moment: u128,
    /// Per-vertex canonical component ID (`None` when inactive).
    pub labels: Vec<Option<usize>>,
}

struct UnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
    active: Vec<bool>,
}

impl UnionFind {
    fn new(vertices: usize) -> Self {
        Self {
            parent: (0..vertices).collect(),
            size: vec![0; vertices],
            active: vec![false; vertices],
        }
    }

    fn enable(&mut self, vertex: usize) {
        self.active[vertex] = true;
        self.size[vertex] = 1;
    }

    fn root(&mut self, mut vertex: usize) -> usize {
        while self.parent[vertex] != vertex {
            self.parent[vertex] = self.parent[self.parent[vertex]];
            vertex = self.parent[vertex];
        }
        vertex
    }

    fn union(&mut self, left: usize, right: usize) {
        let mut left = self.root(left);
        let mut right = self.root(right);
        if left == right {
            return;
        }
        // Keep the smallest vertex ID as the root: canonical components.
        if left > right {
            core::mem::swap(&mut left, &mut right);
        }
        self.parent[right] = left;
        self.size[left] += self.size[right];
        self.size[right] = 0;
    }

    fn finish(
        mut self,
        active_vertices: &[bool],
        active_vertex_count: usize,
        active_edge_count: usize,
    ) -> Result<HypergraphSummary, HypergraphError> {
        let vertices = self.parent.len();
        let mut component_count = 0;
        let mut largest_component_size = 0;
        let mut largest_component = None;
        let mut raw_second_moment = 0_u128;
        for vertex in 0..vertices {
            if !self.active[vertex] {
                continue;
            }
            let root = self.root(vertex);
            if root == vertex {
                let size = self.size[vertex];
                component_count += 1;
                raw_second_moment += (size as u128) * (size as u128);
                if size > largest_component_size
                    || (size == largest_component_size
                        && largest_component.is_none_or(|current| vertex < current))
                {
                    largest_component_size = size;
                    largest_component = Some(vertex);
                }
            }
        }
        debug_assert_eq!(
            self.active.iter().filter(|&&active| active).count(),
            active_vertex_count
        );
        let labels = (0..vertices)
            .map(|vertex| {
                if active_vertices[vertex] {
                    Some(self.root(vertex))
                } else {
                    None
                }
            })
            .collect();
        Ok(HypergraphSummary {
            active_vertex_count,
            active_edge_count,
            component_count,
            largest_component_size,
            largest_component,
            raw_second_moment,
            labels,
        })
    }
}

/// Typed hypergraph rejection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HypergraphError {
    NoVertices,
    EmptyHyperedge {
        edge: usize,
    },
    Member {
        edge: usize,
        member: usize,
        vertices: usize,
    },
    MemberOrder {
        edge: usize,
    },
    Activity {
        entity: &'static str,
        expected: usize,
        actual: usize,
    },
}

impl fmt::Display for HypergraphError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoVertices => write!(formatter, "hypergraph needs at least one vertex"),
            Self::EmptyHyperedge { edge } => write!(formatter, "hyperedge {edge} has no members"),
            Self::Member {
                edge,
                member,
                vertices,
            } => write!(
                formatter,
                "hyperedge {edge}: member {member} outside 0..{vertices}"
            ),
            Self::MemberOrder { edge } => write!(
                formatter,
                "hyperedge {edge}: members must be strictly increasing vertex IDs"
            ),
            Self::Activity {
                entity,
                expected,
                actual,
            } => write!(
                formatter,
                "{entity} activity has length {actual}, expected {expected}"
            ),
        }
    }
}

impl std::error::Error for HypergraphError {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;

    #[test]
    fn one_hyperedges_joins_its_active_members() {
        let graph = Hypergraph::new(4, [vec![0, 1, 2]]).expect("hypergraph");
        let summary = graph
            .analyze(&[true, true, true, false], &[true])
            .expect("analysis");
        assert_eq!(summary.active_vertex_count, 3);
        assert_eq!(summary.active_edge_count, 1);
        assert_eq!(summary.component_count, 1);
        assert_eq!(summary.largest_component_size, 3);
        assert_eq!(summary.largest_component, Some(0));
        assert_eq!(summary.raw_second_moment, 9);
        assert_eq!(summary.labels, vec![Some(0), Some(0), Some(0), None]);
    }

    #[test]
    fn inactive_members_are_excluded_from_the_hyperedge() {
        let graph = Hypergraph::new(4, [vec![0, 1, 2]]).expect("hypergraph");
        let summary = graph
            .analyze(&[true, false, true, true], &[true])
            .expect("analysis");
        // 0 and 2 are joined; 3 is separate (vertex 1 is inactive).
        assert_eq!(summary.component_count, 2);
        assert_eq!(summary.largest_component_size, 2);
        assert_eq!(summary.labels, vec![Some(0), None, Some(0), Some(3)]);
    }

    #[test]
    fn shared_members_merge_hyperedges_and_ties_use_the_lowest_id() {
        let graph = Hypergraph::new(6, [vec![0, 1], vec![1, 2], vec![4, 5]]).expect("hypergraph");
        let summary = graph
            .analyze(&[true; 6], &[true, true, false])
            .expect("analysis");
        // Components {0,1,2} and isolated 3,4,5. Largest is unique.
        assert_eq!(summary.component_count, 4);
        assert_eq!(summary.largest_component_size, 3);
        assert_eq!(summary.largest_component, Some(0));

        let tie = Hypergraph::new(4, [vec![2, 3], vec![0, 1]])
            .expect("hypergraph")
            .analyze(&[true; 4], &[true, true])
            .expect("analysis");
        assert_eq!(tie.component_count, 2);
        assert_eq!(tie.largest_component_size, 2);
        assert_eq!(tie.largest_component, Some(0));
    }

    #[test]
    fn exact_small_enumeration_matches_an_independent_incidence_bfs() {
        fn reference(
            vertices: usize,
            edges: &[Vec<usize>],
            active_vertices: &[bool],
            active_edges: &[bool],
        ) -> (usize, usize, u128, Vec<Option<usize>>) {
            let mut adjacency = vec![Vec::new(); vertices];
            for (members, &active) in edges.iter().zip(active_edges) {
                if !active {
                    continue;
                }
                for &left in members {
                    for &right in members {
                        if left != right && active_vertices[left] && active_vertices[right] {
                            adjacency[left].push(right);
                        }
                    }
                }
            }
            let mut labels = vec![None; vertices];
            let mut component_count = 0;
            let mut largest = 0;
            let mut raw_second_moment = 0_u128;
            for start in 0..vertices {
                if !active_vertices[start] || labels[start].is_some() {
                    continue;
                }
                let mut queue = VecDeque::from([start]);
                labels[start] = Some(start);
                let mut size = 0usize;
                while let Some(vertex) = queue.pop_front() {
                    size += 1;
                    for &neighbor in &adjacency[vertex] {
                        if labels[neighbor].is_none() {
                            labels[neighbor] = Some(start);
                            queue.push_back(neighbor);
                        }
                    }
                }
                component_count += 1;
                largest = largest.max(size);
                raw_second_moment += (size as u128) * (size as u128);
            }
            (component_count, largest, raw_second_moment, labels)
        }

        let edges = vec![vec![0, 1, 2], vec![2, 3], vec![4, 5], vec![0, 5]];
        let graph = Hypergraph::new(6, edges.clone()).expect("hypergraph");
        for vertex_bits in 0..(1 << 6) {
            for edge_bits in 0..(1 << edges.len()) {
                let active_vertices: Vec<bool> = (0..6)
                    .map(|vertex| vertex_bits >> vertex & 1 == 1)
                    .collect();
                let active_edges: Vec<bool> = (0..edges.len())
                    .map(|edge| edge_bits >> edge & 1 == 1)
                    .collect();
                let summary = graph
                    .analyze(&active_vertices, &active_edges)
                    .expect("analysis");
                let (component_count, largest, raw_second_moment, labels) =
                    reference(6, &edges, &active_vertices, &active_edges);
                assert_eq!(summary.component_count, component_count);
                assert_eq!(summary.largest_component_size, largest);
                assert_eq!(summary.raw_second_moment, raw_second_moment);
                assert_eq!(summary.labels, labels);
            }
        }
    }

    #[test]
    fn rejections_are_typed() {
        assert!(matches!(
            Hypergraph::new(0, Vec::<Vec<usize>>::new()),
            Err(HypergraphError::NoVertices)
        ));
        assert!(matches!(
            Hypergraph::new(2, [vec![]]),
            Err(HypergraphError::EmptyHyperedge { edge: 0 })
        ));
        assert!(matches!(
            Hypergraph::new(2, [vec![0, 2]]),
            Err(HypergraphError::Member {
                edge: 0,
                member: 2,
                vertices: 2
            })
        ));
        assert!(matches!(
            Hypergraph::new(3, [vec![1, 1]]),
            Err(HypergraphError::MemberOrder { edge: 0 })
        ));
        let graph = Hypergraph::new(2, [vec![0, 1]]).expect("hypergraph");
        assert!(matches!(
            graph.analyze(&[true], &[true]),
            Err(HypergraphError::Activity {
                entity: "vertex",
                expected: 2,
                actual: 1
            })
        ));
    }
}
