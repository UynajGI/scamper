use core::fmt;
use std::collections::TryReserveError;

use super::super::{
    BoundaryQuery, ComponentLabel, ComponentSummary, ComponentWorkspace, EdgeActivity,
    VertexActivity,
};
use crate::UndirectedGraphView;

/// Borrowed outputs from one component analysis.
///
/// The borrow prevents another analysis from mutating the workspace while this
/// result or its optional labels are in use. Query outcomes preserve input order.
#[derive(Clone, Copy, Debug)]
pub struct AnalysisResult<'a> {
    summary: ComponentSummary,
    query_outcomes: &'a [bool],
    labels: Option<&'a [ComponentLabel]>,
}

impl AnalysisResult<'_> {
    #[inline]
    pub const fn summary(&self) -> ComponentSummary {
        self.summary
    }

    #[inline]
    pub const fn query_outcomes(&self) -> &[bool] {
        self.query_outcomes
    }

    #[inline]
    pub const fn labels(&self) -> Option<&[ComponentLabel]> {
        self.labels
    }
}

/// Analyze without materializing per-vertex labels.
pub fn analyze<'a, G, V, E>(
    graph: &G,
    vertex_activity: &V,
    edge_activity: &E,
    queries: &[BoundaryQuery],
    workspace: &'a mut ComponentWorkspace,
) -> Result<AnalysisResult<'a>, AnalysisError>
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
{
    analyze_inner(
        graph,
        vertex_activity,
        edge_activity,
        queries,
        workspace,
        false,
    )
}

/// Analyze and materialize canonical component labels for every vertex.
pub fn analyze_with_labels<'a, G, V, E>(
    graph: &G,
    vertex_activity: &V,
    edge_activity: &E,
    queries: &[BoundaryQuery],
    workspace: &'a mut ComponentWorkspace,
) -> Result<AnalysisResult<'a>, AnalysisError>
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
{
    analyze_inner(
        graph,
        vertex_activity,
        edge_activity,
        queries,
        workspace,
        true,
    )
}

fn analyze_inner<'a, G, V, E>(
    graph: &G,
    vertex_activity: &V,
    edge_activity: &E,
    queries: &[BoundaryQuery],
    workspace: &'a mut ComponentWorkspace,
    with_labels: bool,
) -> Result<AnalysisResult<'a>, AnalysisError>
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
{
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();
    validate_domains(
        vertex_count,
        edge_count,
        vertex_activity,
        edge_activity,
        queries,
    )?;
    workspace
        .resize_for(vertex_count, queries.len(), with_labels)
        .map_err(AnalysisError::Capacity)?;

    let mut summary = ComponentSummary::EMPTY;
    for index in 0..vertex_count {
        let vertex = graph.vertex_id(index).expect("dense vertex index");
        if vertex_activity.vertex_active(vertex) {
            workspace.parent[index] = index;
            workspace.size[index] = 1;
            workspace.canonical[index] = index;
            summary.active_vertex_count += 1;
        } else {
            workspace.parent[index] = usize::MAX;
            workspace.size[index] = 0;
            workspace.canonical[index] = usize::MAX;
        }
    }

    for index in 0..edge_count {
        let edge = graph.edge_id(index).expect("dense physical edge index");
        if !edge_activity.edge_active(edge) {
            continue;
        }
        let [left, right] = graph.edge_endpoints(edge);
        let left = left.index();
        let right = right.index();
        if workspace.parent[left] == usize::MAX || workspace.parent[right] == usize::MAX {
            continue;
        }
        summary.active_edge_count += 1;
        union(workspace, left, right);
    }

    for index in 0..vertex_count {
        if workspace.parent[index] != index {
            continue;
        }
        let size = workspace.size[index];
        let identity = workspace.canonical[index];
        summary.component_count += 1;
        summary.raw_second_moment += (size as u128) * (size as u128);
        if size > summary.largest_component_size
            || (size == summary.largest_component_size
                && summary
                    .largest_component
                    .is_none_or(|current| identity < current.index()))
        {
            summary.largest_component_size = size;
            summary.largest_component = graph.vertex_id(identity);
        }
    }

    evaluate_queries(queries, workspace);

    let labels = if with_labels {
        for index in 0..vertex_count {
            workspace.labels[index] = if workspace.parent[index] == usize::MAX {
                None
            } else {
                let root = find(&mut workspace.parent, index);
                graph.vertex_id(workspace.canonical[root])
            };
        }
        Some(workspace.labels.as_slice())
    } else {
        None
    };

    Ok(AnalysisResult {
        summary,
        query_outcomes: &workspace.query_outcomes,
        labels,
    })
}

fn validate_domains<V: VertexActivity, E: EdgeActivity>(
    vertex_count: usize,
    edge_count: usize,
    vertex_activity: &V,
    edge_activity: &E,
    queries: &[BoundaryQuery],
) -> Result<(), AnalysisError> {
    if vertex_activity.vertex_count() != vertex_count {
        return Err(AnalysisError::VertexActivityLength {
            expected: vertex_count,
            actual: vertex_activity.vertex_count(),
        });
    }
    if edge_activity.edge_count() != edge_count {
        return Err(AnalysisError::EdgeActivityLength {
            expected: edge_count,
            actual: edge_activity.edge_count(),
        });
    }
    if let Some((query, actual)) = queries.iter().enumerate().find_map(|(index, query)| {
        (query.vertex_count() != vertex_count).then_some((index, query.vertex_count()))
    }) {
        return Err(AnalysisError::QueryVertexCount {
            query,
            expected: vertex_count,
            actual,
        });
    }
    Ok(())
}

fn find(parent: &mut [usize], node: usize) -> usize {
    let mut root = node;
    while parent[root] != root {
        root = parent[root];
    }
    let mut current = node;
    while parent[current] != root {
        let next = parent[current];
        parent[current] = root;
        current = next;
    }
    root
}

fn union(workspace: &mut ComponentWorkspace, left: usize, right: usize) {
    let left = find(&mut workspace.parent, left);
    let right = find(&mut workspace.parent, right);
    if left == right {
        return;
    }
    let (big, small) = if workspace.size[left] >= workspace.size[right] {
        (left, right)
    } else {
        (right, left)
    };
    workspace.parent[small] = big;
    workspace.size[big] += workspace.size[small];
    workspace.canonical[big] = workspace.canonical[big].min(workspace.canonical[small]);
}

fn evaluate_queries(queries: &[BoundaryQuery], workspace: &mut ComponentWorkspace) {
    for (query_index, query) in queries.iter().enumerate() {
        let stamp = workspace.next_stamp();
        for &vertex in query.from() {
            let index = vertex.index();
            if workspace.parent[index] != usize::MAX {
                let root = find(&mut workspace.parent, index);
                workspace.query_root_stamps[root] = stamp;
            }
        }
        workspace.query_outcomes[query_index] = query.to().iter().any(|&vertex| {
            let index = vertex.index();
            if workspace.parent[index] == usize::MAX {
                return false;
            }
            let root = find(&mut workspace.parent, index);
            workspace.query_root_stamps[root] == stamp
        });
    }
}

#[derive(Debug)]
pub enum AnalysisError {
    VertexActivityLength {
        expected: usize,
        actual: usize,
    },
    EdgeActivityLength {
        expected: usize,
        actual: usize,
    },
    QueryVertexCount {
        query: usize,
        expected: usize,
        actual: usize,
    },
    Capacity(TryReserveError),
}

impl fmt::Display for AnalysisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::VertexActivityLength { expected, actual } => write!(
                formatter,
                "vertex activity length is {actual}, expected {expected}"
            ),
            Self::EdgeActivityLength { expected, actual } => write!(
                formatter,
                "edge activity length is {actual}, expected {expected}"
            ),
            Self::QueryVertexCount {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} was validated for {actual} vertices, expected {expected}"
            ),
            Self::Capacity(_) => write!(formatter, "could not reserve component workspace"),
        }
    }
}

impl std::error::Error for AnalysisError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
            _ => None,
        }
    }
}
