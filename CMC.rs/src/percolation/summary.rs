use crate::VertexId;

/// Raw component facts for one active undirected configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ComponentSummary {
    /// Vertices for which `VertexActivity::vertex_active` returned true.
    pub active_vertex_count: usize,
    /// Active physical edges whose two endpoints are active. Self-loops count once.
    pub active_edge_count: usize,
    /// Connected components among active vertices, including active isolates.
    pub component_count: usize,
    /// Canonical identity of the selected largest component.
    ///
    /// A component identity is its lowest vertex ID. Equal-size largest
    /// components are resolved by the lowest such identity.
    pub largest_component: Option<VertexId>,
    /// Number of active vertices in `largest_component`, or zero when empty.
    pub largest_component_size: usize,
    /// Raw `sum(s_i^2)` over all active components.
    pub raw_second_moment: u128,
}

impl ComponentSummary {
    pub(crate) const EMPTY: Self = Self {
        active_vertex_count: 0,
        active_edge_count: 0,
        component_count: 0,
        largest_component: None,
        largest_component_size: 0,
        raw_second_moment: 0,
    };
}
