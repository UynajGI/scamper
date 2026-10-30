//! Shared plan recording for scheduler adapters.
//!
//! Every owned adapter runs the same tail in `measure()`: analyze the
//! current realization (materializing labels only when the plan asks),
//! borrow the result into a [`SampleView`], and record the plan. The
//! samplers differ; this tail does not. Keeping it in one place means the
//! label branch and the `SampleView` construction cannot drift between
//! adapters.

use carlo_rs::Context;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{
    analyze, analyze_bonds, analyze_bonds_with_labels, analyze_with_labels, ComponentSummary,
    ComponentWorkspace, EdgeActivity, ObservablePlan, SampleView, VertexActivity, WindingSummary,
};
use crate::UndirectedGraphView;

/// Analyze one vertex/edge activity pair and record the plan.
///
/// Returns the component summary so adapters that maintain incremental
/// invariants can cross-check the analyzer.
pub(crate) fn record_component_analysis<G, V, E>(
    graph: &G,
    vertices: &V,
    edges: &E,
    plan: &ObservablePlan,
    workspace: &mut ComponentWorkspace,
    winding: Option<WindingSummary>,
    ctx: &mut Context<Xoshiro256PlusPlus>,
) -> ComponentSummary
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
{
    let result = if plan.requests_labels() {
        analyze_with_labels(graph, vertices, edges, plan.boundary_queries(), workspace)
    } else {
        analyze(graph, vertices, edges, plan.boundary_queries(), workspace)
    }
    .expect("constructor prepared matching analysis domains");
    let summary = result.summary();
    let sample = SampleView::from_analysis(&result, graph.vertex_count(), winding);
    plan.record(&sample, ctx);
    summary
}

/// Analyze one bond mask (every vertex active) and record the plan.
pub(crate) fn record_bond_analysis<G>(
    graph: &G,
    open_bonds: &[bool],
    plan: &ObservablePlan,
    workspace: &mut ComponentWorkspace,
    ctx: &mut Context<Xoshiro256PlusPlus>,
) -> ComponentSummary
where
    G: UndirectedGraphView,
{
    let result = if plan.requests_labels() {
        analyze_bonds_with_labels(graph, open_bonds, plan.boundary_queries(), workspace)
    } else {
        analyze_bonds(graph, open_bonds, plan.boundary_queries(), workspace)
    }
    .expect("constructor prepared matching analysis domains");
    let summary = result.summary();
    let sample = SampleView::from_analysis(&result, graph.vertex_count(), None);
    plan.record(&sample, ctx);
    summary
}
