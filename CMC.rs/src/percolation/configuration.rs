use super::activity::ActivityMask;
use super::{EdgeActivity, VertexActivity};
use crate::{EdgeId, VertexId};

/// Reusable activity-only configuration for static percolation.
///
/// It owns no topology, law, probability, mode, analyzer, or workspace.
#[derive(Clone, Debug)]
pub struct StaticConfiguration {
    vertices: ActivityMask,
    edges: ActivityMask,
}

impl StaticConfiguration {
    /// Create an all-inactive configuration with the requested dense domains.
    pub const fn new(vertex_count: usize, edge_count: usize) -> Self {
        Self {
            vertices: ActivityMask::none(vertex_count),
            edges: ActivityMask::none(edge_count),
        }
    }

    #[inline]
    pub const fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    #[inline]
    pub const fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Resize both activity domains and clear all logical activity.
    pub fn resize(&mut self, vertex_count: usize, edge_count: usize) {
        self.vertices.reconfigure(vertex_count, false);
        self.edges.reconfigure(edge_count, false);
    }

    /// Bytes reserved by both private dense activity backends.
    #[cfg(test)]
    #[inline]
    pub(crate) fn owned_mask_bytes(&self) -> usize {
        self.vertices.owned_bytes() + self.edges.owned_bytes()
    }

    /// Set one vertex's activity.
    ///
    /// Public so a third-party [`crate::OccupationLaw`] can write its own
    /// realization without access to the private dense masks. The ID must
    /// come from the same topology the configuration was sized for.
    #[inline]
    pub fn set_vertex_active(&mut self, vertex: VertexId, active: bool) {
        assert!(
            self.vertices.set(vertex.index(), active),
            "vertex ID outside configuration domain"
        );
    }

    /// Set one physical edge's activity; see [`Self::set_vertex_active`].
    #[inline]
    pub fn set_edge_active(&mut self, edge: EdgeId, active: bool) {
        assert!(
            self.edges.set(edge.index(), active),
            "edge ID outside configuration domain"
        );
    }

    pub(crate) fn vertex_mask_mut(&mut self) -> &mut ActivityMask {
        &mut self.vertices
    }

    pub(crate) fn edge_mask_mut(&mut self) -> &mut ActivityMask {
        &mut self.edges
    }
}

impl VertexActivity for StaticConfiguration {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.vertex_count()
    }

    #[inline]
    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.vertices
            .get(vertex.index())
            .expect("vertex ID outside configuration domain")
    }
}

impl EdgeActivity for StaticConfiguration {
    #[inline]
    fn edge_count(&self) -> usize {
        self.edge_count()
    }

    #[inline]
    fn edge_active(&self, edge: EdgeId) -> bool {
        self.edges
            .get(edge.index())
            .expect("edge ID outside configuration domain")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GraphView, UndirectedGraphView};

    #[test]
    fn zero_size_and_resize_are_well_defined() {
        let mut configuration = StaticConfiguration::new(0, 0);
        assert_eq!(configuration.vertex_count(), 0);
        assert_eq!(configuration.edge_count(), 0);
        assert_eq!(configuration.owned_mask_bytes(), 0);

        configuration.resize(4, 3);
        assert_eq!(configuration.vertex_count(), 4);
        assert_eq!(configuration.edge_count(), 3);
        configuration.resize(0, 0);
        assert_eq!(configuration.vertex_count(), 0);
        assert_eq!(configuration.edge_count(), 0);
    }

    #[test]
    fn public_setters_drive_the_activity_traits() {
        let graph = crate::build_chain(4, false);
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        configuration.set_vertex_active(graph.vertex_id(1).expect("vertex"), true);
        configuration.set_edge_active(graph.edge_id(2).expect("edge"), true);
        assert!(configuration.vertex_active(graph.vertex_id(1).expect("vertex")));
        assert!(!configuration.vertex_active(graph.vertex_id(0).expect("vertex")));
        assert!(configuration.edge_active(graph.edge_id(2).expect("edge")));
        assert!(!configuration.edge_active(graph.edge_id(0).expect("edge")));
    }

    #[test]
    #[should_panic(expected = "vertex ID outside configuration domain")]
    fn public_setter_rejects_a_foreign_vertex_id() {
        let graph = crate::build_chain(4, false);
        let foreign = graph.vertex_id(3).expect("vertex");
        let mut configuration = StaticConfiguration::new(3, graph.edge_count());
        configuration.set_vertex_active(foreign, true);
    }
}
