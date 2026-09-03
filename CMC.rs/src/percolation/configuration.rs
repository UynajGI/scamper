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
}
