use core::fmt;
use std::collections::TryReserveError;
#[cfg(feature = "allocation-probe")]
use std::mem::size_of;

use rand::{Rng, RngExt};

use crate::UndirectedGraphView;

/// Occupation family of one Newman-Ziff scan.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScanMode {
    Site,
    Bond,
}

/// Component facts after a prefix of the activation order.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ScanPoint {
    /// Active vertices (site) or active physical edges (bond).
    pub occupation: usize,
    /// Active physical edges with both endpoints active.
    pub active_edge_count: usize,
    /// Connected components among tracked vertices, including isolates.
    pub component_count: usize,
    /// Largest component size; zero when nothing is active.
    pub largest_component_size: usize,
    /// Raw `sum(s_i^2)` over tracked components.
    pub raw_second_moment: u128,
}

impl ScanPoint {
    const EMPTY: Self = Self {
        occupation: 0,
        active_edge_count: 0,
        component_count: 0,
        largest_component_size: 0,
        raw_second_moment: 0,
    };
}

/// Borrowed occupation-count curve from one Newman-Ziff scan.
///
/// `points[k]` describes the configuration after the first `k` entries of the
/// activation order. The order contains dense vertex indices for site scans and
/// dense physical-edge indices for bond scans.
#[derive(Clone, Copy, Debug)]
pub struct ScanCurve<'a> {
    mode: ScanMode,
    trials: usize,
    total_vertices: usize,
    order: &'a [usize],
    points: &'a [ScanPoint],
}

impl<'a> ScanCurve<'a> {
    #[inline]
    pub const fn mode(&self) -> ScanMode {
        self.mode
    }

    /// Number of activation steps: `V` for site scans and `E` for bond scans.
    #[inline]
    pub const fn trials(&self) -> usize {
        self.trials
    }

    /// Total topology vertices, also the giant-fraction denominator.
    #[inline]
    pub const fn total_vertices(&self) -> usize {
        self.total_vertices
    }

    /// Activation order (dense vertex or edge indices, depending on the mode).
    #[inline]
    pub const fn activation_order(&self) -> &'a [usize] {
        self.order
    }

    /// Curve points; exactly `trials + 1` entries.
    #[inline]
    pub const fn points(&self) -> &'a [ScanPoint] {
        self.points
    }

    #[inline]
    pub fn point(&self, occupation: usize) -> Option<&'a ScanPoint> {
        self.points.get(occupation)
    }

    pub(crate) const fn from_parts(
        mode: ScanMode,
        trials: usize,
        total_vertices: usize,
        order: &'a [usize],
        points: &'a [ScanPoint],
    ) -> Self {
        Self {
            mode,
            trials,
            total_vertices,
            order,
            points,
        }
    }
}

/// Reusable storage for Newman-Ziff scans.
///
/// Call [`Self::prepare`] once for a fixed topology; repeated site and bond
/// scans then reuse the union-find, activation and curve buffers without
/// allocating.
#[derive(Clone, Debug, Default)]
pub struct NewmanZiffWorkspace {
    parent: Vec<usize>,
    size: Vec<usize>,
    vertex_active: Vec<u8>,
    edge_counted: Vec<u8>,
    order: Vec<usize>,
    points: Vec<ScanPoint>,
    prepared_vertices: usize,
    prepared_edges: usize,
}

impl NewmanZiffWorkspace {
    pub const fn new() -> Self {
        Self {
            parent: Vec::new(),
            size: Vec::new(),
            vertex_active: Vec::new(),
            edge_counted: Vec::new(),
            order: Vec::new(),
            points: Vec::new(),
            prepared_vertices: 0,
            prepared_edges: 0,
        }
    }

    /// Reserve storage for every scan on a topology with these sizes.
    pub fn prepare(
        &mut self,
        vertex_count: usize,
        edge_count: usize,
    ) -> Result<(), TryReserveError> {
        reserve_len(&mut self.parent, vertex_count)?;
        reserve_len(&mut self.size, vertex_count)?;
        reserve_len(&mut self.vertex_active, vertex_count)?;
        reserve_len(&mut self.edge_counted, edge_count)?;
        let trials = vertex_count.max(edge_count);
        reserve_len(&mut self.order, trials)?;
        reserve_len(&mut self.points, trials.saturating_add(1))?;
        self.prepared_vertices = vertex_count;
        self.prepared_edges = edge_count;
        Ok(())
    }

    /// Prepared topology sizes.
    #[inline]
    pub const fn prepared_sizes(&self) -> (usize, usize) {
        (self.prepared_vertices, self.prepared_edges)
    }

    #[inline]
    pub fn vertex_capacity(&self) -> usize {
        self.parent
            .capacity()
            .min(self.size.capacity())
            .min(self.vertex_active.capacity())
    }

    #[inline]
    pub fn edge_capacity(&self) -> usize {
        self.edge_counted.capacity()
    }

    #[inline]
    pub fn order_capacity(&self) -> usize {
        self.order.capacity()
    }

    #[inline]
    pub fn point_capacity(&self) -> usize {
        self.points.capacity()
    }

    /// Audit-only retained-capacity accounting for the allocation probe.
    #[cfg(feature = "allocation-probe")]
    #[doc(hidden)]
    pub fn capacity_audit(&self) -> WorkspaceScanCapacityAudit {
        WorkspaceScanCapacityAudit {
            parent: vector_capacity_bytes(&self.parent),
            size: vector_capacity_bytes(&self.size),
            vertex_active: vector_capacity_bytes(&self.vertex_active),
            edge_counted: vector_capacity_bytes(&self.edge_counted),
            order: vector_capacity_bytes(&self.order),
            points: vector_capacity_bytes(&self.points),
        }
    }

    fn resize_for(
        &mut self,
        mode: ScanMode,
        vertex_count: usize,
        edge_count: usize,
    ) -> Result<(), TryReserveError> {
        self.prepare(vertex_count, edge_count)?;
        let trials = match mode {
            ScanMode::Site => vertex_count,
            ScanMode::Bond => edge_count,
        };
        self.parent.resize(vertex_count, usize::MAX);
        self.size.resize(vertex_count, 0);
        self.vertex_active.resize(vertex_count, 0);
        self.edge_counted.resize(edge_count, 0);
        self.order.resize(trials, 0);
        self.points
            .resize(trials.saturating_add(1), ScanPoint::EMPTY);
        Ok(())
    }
}

/// Run a site-percolation Newman-Ziff scan.
///
/// Vertices are activated in a uniformly random order. An edge becomes active
/// when both endpoints are active, so `points[k].active_edge_count` is the
/// induced active-edge count after `k` vertex activations.
pub fn scan_site<'a, G, R>(
    graph: &G,
    workspace: &'a mut NewmanZiffWorkspace,
    rng: &mut R,
) -> Result<ScanCurve<'a>, ScanError>
where
    G: UndirectedGraphView,
    R: Rng + ?Sized,
{
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();
    workspace
        .resize_for(ScanMode::Site, vertex_count, edge_count)
        .map_err(ScanError::Capacity)?;

    let NewmanZiffWorkspace {
        parent,
        size,
        vertex_active,
        edge_counted,
        order,
        points,
        ..
    } = workspace;
    parent.fill(usize::MAX);
    size.fill(0);
    vertex_active.fill(0);
    edge_counted.fill(0);
    fill_and_shuffle(order, rng);

    points[0] = ScanPoint::EMPTY;
    let mut components = 0_usize;
    let mut second_moment = 0_u128;
    let mut largest = 0_usize;
    let mut active_edges = 0_usize;

    for (position, &vertex_index) in order.iter().enumerate() {
        let vertex = graph.vertex_id(vertex_index).expect("dense vertex index");
        vertex_active[vertex_index] = 1;
        parent[vertex_index] = vertex_index;
        size[vertex_index] = 1;
        components += 1;
        second_moment += 1;
        largest = largest.max(1);

        for incidence in graph.incidences(vertex) {
            let neighbor = incidence.neighbor.index();
            let edge = incidence.edge.index();
            if vertex_active[neighbor] == 0 || edge_counted[edge] != 0 {
                continue;
            }
            edge_counted[edge] = 1;
            active_edges += 1;
            if let Some((left_size, right_size)) = union(parent, size, vertex_index, neighbor) {
                components -= 1;
                second_moment += 2 * (left_size as u128) * (right_size as u128);
                largest = largest.max(left_size + right_size);
            }
        }

        points[position + 1] = ScanPoint {
            occupation: position + 1,
            active_edge_count: active_edges,
            component_count: components,
            largest_component_size: largest,
            raw_second_moment: second_moment,
        };
    }

    Ok(ScanCurve::from_parts(
        ScanMode::Site,
        vertex_count,
        vertex_count,
        &order[..vertex_count],
        &points[..vertex_count + 1],
    ))
}

/// Run a bond-percolation Newman-Ziff scan.
///
/// Physical edges are activated in a uniformly random order. Every vertex is
/// tracked from the start, so isolated vertices count as components exactly as
/// in the F3 bond analyzer.
pub fn scan_bond<'a, G, R>(
    graph: &G,
    workspace: &'a mut NewmanZiffWorkspace,
    rng: &mut R,
) -> Result<ScanCurve<'a>, ScanError>
where
    G: UndirectedGraphView,
    R: Rng + ?Sized,
{
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();
    workspace
        .resize_for(ScanMode::Bond, vertex_count, edge_count)
        .map_err(ScanError::Capacity)?;

    let NewmanZiffWorkspace {
        parent,
        size,
        order,
        points,
        ..
    } = workspace;
    for index in 0..vertex_count {
        parent[index] = index;
        size[index] = 1;
    }
    fill_and_shuffle(order, rng);

    let mut components = vertex_count;
    let mut second_moment = vertex_count as u128;
    let mut largest = usize::from(vertex_count > 0);
    points[0] = ScanPoint {
        occupation: 0,
        active_edge_count: 0,
        component_count: components,
        largest_component_size: largest,
        raw_second_moment: second_moment,
    };

    for (position, &edge_index) in order.iter().enumerate() {
        let edge = graph
            .edge_id(edge_index)
            .expect("dense physical edge index");
        let [left, right] = graph.edge_endpoints(edge);
        if let Some((left_size, right_size)) = union(parent, size, left.index(), right.index()) {
            components -= 1;
            second_moment += 2 * (left_size as u128) * (right_size as u128);
            largest = largest.max(left_size + right_size);
        }
        points[position + 1] = ScanPoint {
            occupation: position + 1,
            active_edge_count: position + 1,
            component_count: components,
            largest_component_size: largest,
            raw_second_moment: second_moment,
        };
    }

    Ok(ScanCurve::from_parts(
        ScanMode::Bond,
        edge_count,
        vertex_count,
        &order[..edge_count],
        &points[..edge_count + 1],
    ))
}

fn fill_and_shuffle<R: Rng + ?Sized>(order: &mut [usize], rng: &mut R) {
    for (index, slot) in order.iter_mut().enumerate() {
        *slot = index;
    }
    for index in (1..order.len()).rev() {
        let swap_with = rng.random_range(0..=index);
        order.swap(index, swap_with);
    }
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

/// Unions two roots and returns their sizes before the merge when distinct.
fn union(
    parent: &mut [usize],
    size: &mut [usize],
    left: usize,
    right: usize,
) -> Option<(usize, usize)> {
    let left_root = find(parent, left);
    let right_root = find(parent, right);
    if left_root == right_root {
        return None;
    }
    let (big, small) = if size[left_root] >= size[right_root] {
        (left_root, right_root)
    } else {
        (right_root, left_root)
    };
    let (big_size, small_size) = (size[big], size[small]);
    parent[small] = big;
    size[big] = big_size + small_size;
    Some((big_size, small_size))
}

fn reserve_len<T>(storage: &mut Vec<T>, len: usize) -> Result<(), TryReserveError> {
    if storage.capacity() < len {
        storage.try_reserve_exact(len - storage.len())?;
    }
    Ok(())
}

/// Capacity failure of a Newman-Ziff scan.
#[derive(Debug)]
pub enum ScanError {
    Capacity(TryReserveError),
}

impl fmt::Display for ScanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Capacity(_) => write!(formatter, "could not reserve Newman-Ziff scan workspace"),
        }
    }
}

impl std::error::Error for ScanError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Capacity(source) => Some(source),
        }
    }
}

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WorkspaceScanCapacityAudit {
    pub parent: usize,
    pub size: usize,
    pub vertex_active: usize,
    pub edge_counted: usize,
    pub order: usize,
    pub points: usize,
}

#[cfg(feature = "allocation-probe")]
impl WorkspaceScanCapacityAudit {
    pub const fn total_bytes(self) -> usize {
        self.parent
            .saturating_add(self.size)
            .saturating_add(self.vertex_active)
            .saturating_add(self.edge_counted)
            .saturating_add(self.order)
            .saturating_add(self.points)
    }
}

#[cfg(feature = "allocation-probe")]
fn vector_capacity_bytes<T>(storage: &Vec<T>) -> usize {
    storage.capacity().saturating_mul(size_of::<T>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, Bond, BondType, BorrowedUndirectedCsr, CsrLattice};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    fn rng(seed: u64) -> Xoshiro256PlusPlus {
        Xoshiro256PlusPlus::seed_from_u64(seed)
    }

    #[test]
    fn empty_graph_yields_one_empty_point() {
        let graph = BorrowedUndirectedCsr::new(&[0], &[], &[], &[]).unwrap();
        let mut workspace = NewmanZiffWorkspace::new();
        let mut rng = rng(7);

        let site = scan_site(&graph, &mut workspace, &mut rng).unwrap();
        assert_eq!(site.mode(), ScanMode::Site);
        assert_eq!(site.trials(), 0);
        assert_eq!(site.total_vertices(), 0);
        assert_eq!(site.points(), &[ScanPoint::EMPTY]);
        assert!(site.activation_order().is_empty());

        let bond = scan_bond(&graph, &mut workspace, &mut rng).unwrap();
        assert_eq!(bond.mode(), ScanMode::Bond);
        assert_eq!(bond.trials(), 0);
        assert_eq!(bond.points(), &[ScanPoint::EMPTY]);
    }

    #[test]
    fn self_loop_counts_once_in_site_scan() {
        let graph =
            CsrLattice::try_from_edges(2, vec![Bond::new(0, 0, BondType::Generic, 1.0)]).unwrap();
        let mut workspace = NewmanZiffWorkspace::new();
        let mut rng = rng(11);
        let curve = scan_site(&graph, &mut workspace, &mut rng).unwrap();

        let activated = curve
            .activation_order()
            .iter()
            .position(|&vertex| vertex == 0)
            .unwrap();
        for (occupation, point) in curve.points().iter().enumerate() {
            let expected = usize::from(occupation > activated);
            assert_eq!(
                point.active_edge_count, expected,
                "self-loop count at occupation {occupation}"
            );
        }
    }

    #[test]
    fn parallel_edges_count_separately() {
        let graph = CsrLattice::try_from_edges(
            2,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 1.0),
            ],
        )
        .unwrap();
        let mut workspace = NewmanZiffWorkspace::new();
        let mut rng = rng(13);

        let site = scan_site(&graph, &mut workspace, &mut rng).unwrap();
        assert_eq!(site.points()[2].active_edge_count, 2);
        assert_eq!(site.points()[2].component_count, 1);
        assert_eq!(site.points()[2].largest_component_size, 2);
        assert_eq!(site.points()[2].raw_second_moment, 4);

        let bond = scan_bond(&graph, &mut workspace, &mut rng).unwrap();
        assert_eq!(bond.points()[2].active_edge_count, 2);
        assert_eq!(bond.points()[2].component_count, 1);
        assert_eq!(bond.points()[2].largest_component_size, 2);
        assert_eq!(bond.points()[2].raw_second_moment, 4);
    }

    #[test]
    fn bond_scan_tracks_isolates_and_hand_checks_first_union() {
        let graph = build_chain(4, false);
        let mut workspace = NewmanZiffWorkspace::new();
        let mut rng = rng(17);
        let curve = scan_bond(&graph, &mut workspace, &mut rng).unwrap();

        assert_eq!(
            curve.points()[0],
            ScanPoint {
                occupation: 0,
                active_edge_count: 0,
                component_count: 4,
                largest_component_size: 1,
                raw_second_moment: 4,
            }
        );
        assert_eq!(
            curve.points()[1],
            ScanPoint {
                occupation: 1,
                active_edge_count: 1,
                component_count: 3,
                largest_component_size: 2,
                raw_second_moment: 6,
            }
        );
    }

    #[test]
    fn seeded_scans_are_reproducible_and_workspace_is_reusable() {
        let graph = build_chain(8, false);
        let mut workspace = NewmanZiffWorkspace::new();
        workspace.prepare(graph.n_sites, graph.n_edges()).unwrap();

        let mut first_rng = rng(21);
        let first = scan_site(&graph, &mut workspace, &mut first_rng).unwrap();
        let first_order = first.activation_order().to_vec();
        let first_points = first.points().to_vec();

        let mut second_rng = rng(21);
        let second = scan_site(&graph, &mut workspace, &mut second_rng).unwrap();
        assert_eq!(second.activation_order(), first_order.as_slice());
        assert_eq!(second.points(), first_points.as_slice());

        let mut third_rng = rng(22);
        let bond = scan_bond(&graph, &mut workspace, &mut third_rng).unwrap();
        assert_eq!(bond.mode(), ScanMode::Bond);
        assert_eq!(bond.trials(), graph.n_edges());
    }
}
