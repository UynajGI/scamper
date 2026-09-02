//! Union-find cluster analysis of occupied subgraphs.

use super::state::{OccupancyState, PercolationMode};
use crate::lattice::graph::CsrLattice;

/// Disjoint-set forest with union by size and path compression.
#[derive(Debug, Clone)]
pub struct UnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl UnionFind {
    /// `n` singleton clusters.
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    /// Root of `x`, compressing the traversed path.
    pub fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut node = x;
        while self.parent[node] != root {
            let next = self.parent[node];
            self.parent[node] = root;
            node = next;
        }
        root
    }

    /// Merge the clusters of `left` and `right` (no-op when already equal).
    pub fn union(&mut self, left: usize, right: usize) {
        let left = self.find(left);
        let right = self.find(right);
        if left == right {
            return;
        }
        let (big, small) = if self.size[left] >= self.size[right] {
            (left, right)
        } else {
            (right, left)
        };
        self.parent[small] = big;
        self.size[big] += self.size[small];
    }

    /// Cluster size at a root, as returned by [`UnionFind::find`].
    pub fn size_of(&self, root: usize) -> usize {
        self.size[root]
    }
}

/// Cluster observables of one occupancy sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClusterStats {
    /// Largest cluster size in sites.
    pub max_size: usize,
    /// Sum of squared cluster sizes `sum(s_i^2)`.
    pub second_moment: u64,
    /// Number of clusters among tracked sites.
    pub n_clusters: usize,
    /// Whether one cluster contains sites of both spanning sets.
    pub spanning: bool,
}

/// Cluster statistics of the occupied subgraph of `lattice`.
///
/// Site percolation connects an edge when both endpoints are open and counts
/// clusters over open sites only. Bond percolation connects open edges and
/// counts clusters over all sites, so isolated sites are singleton clusters.
/// `from`/`to` are spanning test sets: [`ClusterStats::spanning`] is set when
/// a single cluster contains at least one site of each. All set entries must
/// be in range; disjoint sets keep `spanning = false` at `p = 0` meaningful.
pub fn cluster_stats(
    lattice: &CsrLattice,
    occupancy: &OccupancyState,
    from: &[usize],
    to: &[usize],
) -> ClusterStats {
    let n_sites = lattice.n_sites;
    let mut uf = UnionFind::new(n_sites);
    for (edge_id, edge) in lattice.edges.iter().enumerate() {
        let connected = match occupancy.mode {
            PercolationMode::Site => {
                occupancy.site_open[edge.source] && occupancy.site_open[edge.target]
            }
            PercolationMode::Bond => occupancy.bond_open[edge_id],
        };
        if connected {
            uf.union(edge.source, edge.target);
        }
    }

    let mut in_from = vec![false; n_sites];
    let mut in_to = vec![false; n_sites];
    for &site in from {
        in_from[site] = true;
    }
    for &site in to {
        in_to[site] = true;
    }

    let mut counted = vec![false; n_sites];
    let mut root_touches_from = vec![false; n_sites];
    let mut root_touches_to = vec![false; n_sites];
    let mut stats = ClusterStats {
        max_size: 0,
        second_moment: 0,
        n_clusters: 0,
        spanning: false,
    };
    for site in 0..n_sites {
        let tracked = match occupancy.mode {
            PercolationMode::Site => occupancy.site_open[site],
            PercolationMode::Bond => true,
        };
        if !tracked {
            continue;
        }
        let root = uf.find(site);
        if !counted[root] {
            counted[root] = true;
            stats.n_clusters += 1;
            let size = uf.size_of(root);
            stats.max_size = stats.max_size.max(size);
            stats.second_moment += (size as u64) * (size as u64);
        }
        root_touches_from[root] |= in_from[site];
        root_touches_to[root] |= in_to[site];
        if root_touches_from[root] && root_touches_to[root] {
            stats.spanning = true;
        }
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::graph::{build_chain, build_square, Bond, CsrLattice};

    #[test]
    fn union_find_tracks_sizes() {
        let mut uf = UnionFind::new(6);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(3, 4);
        let root_left = uf.find(0);
        let root_mid = uf.find(3);
        let root_free = uf.find(5);
        assert_eq!(root_left, uf.find(2));
        assert_eq!(root_mid, uf.find(4));
        assert_ne!(root_left, root_mid);
        assert_ne!(root_left, root_free);
        assert_eq!(uf.size_of(root_left), 3);
        assert_eq!(uf.size_of(root_mid), 2);
        assert_eq!(uf.size_of(root_free), 1);
        // Redundant unions are no-ops.
        uf.union(2, 0);
        uf.union(4, 3);
        let root = uf.find(0);
        assert_eq!(uf.size_of(root), 3);
    }

    #[test]
    fn site_percolation_semantics_on_two_by_two() {
        //  0 1
        //  2 3   edges: (0,1) (0,2) (1,3) (2,3)
        let lattice = build_square(2, 2, false);
        let from = [0, 2];
        let to = [1, 3];

        let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Site);
        occupancy
            .site_open
            .copy_from_slice(&[true, false, false, true]);
        let stats = cluster_stats(&lattice, &occupancy, &from, &to);
        assert_eq!(
            stats,
            ClusterStats {
                max_size: 1,
                second_moment: 2,
                n_clusters: 2,
                spanning: false,
            }
        );

        occupancy
            .site_open
            .copy_from_slice(&[true, true, true, false]);
        let stats = cluster_stats(&lattice, &occupancy, &from, &to);
        // Cluster {0,1,2} touches column 0 (0, 2) and column 1 (1).
        assert_eq!(
            stats,
            ClusterStats {
                max_size: 3,
                second_moment: 9,
                n_clusters: 1,
                spanning: true,
            }
        );
    }

    #[test]
    fn bond_percolation_counts_isolated_sites() {
        // Chain 0 - 1 - 2 with only the first bond open.
        let lattice = build_chain(3, false);
        let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Bond);
        occupancy.bond_open.copy_from_slice(&[true, false]);
        let stats = cluster_stats(&lattice, &occupancy, &[0], &[2]);
        assert_eq!(
            stats,
            ClusterStats {
                max_size: 2,
                second_moment: 5,
                n_clusters: 2,
                spanning: false,
            }
        );

        occupancy.bond_open.copy_from_slice(&[true, true]);
        let stats = cluster_stats(&lattice, &occupancy, &[0], &[2]);
        assert!(stats.spanning);
        assert_eq!(stats.max_size, 3);
        assert_eq!(stats.n_clusters, 1);
    }

    #[test]
    fn extremes_match_percolation_physics() {
        let lattice = build_square(3, 3, false);
        let from = vec![0, 3, 6];
        let to = vec![2, 5, 8];
        for mode in [PercolationMode::Site, PercolationMode::Bond] {
            let occupancy = OccupancyState::new(&lattice, mode);
            let stats = cluster_stats(&lattice, &occupancy, &from, &to);
            assert!(!stats.spanning, "p = 0 cannot span");
            match mode {
                // Site: nothing occupied, so no clusters exist at all.
                PercolationMode::Site => {
                    assert_eq!(stats.max_size, 0);
                    assert_eq!(stats.second_moment, 0);
                    assert_eq!(stats.n_clusters, 0);
                }
                // Bond: every site is an isolated singleton cluster.
                PercolationMode::Bond => {
                    assert_eq!(stats.max_size, 1);
                    assert_eq!(stats.second_moment, 9);
                    assert_eq!(stats.n_clusters, 9);
                }
            }

            let mut all_open = OccupancyState::new(&lattice, mode);
            all_open.site_open.fill(true);
            all_open.bond_open.fill(true);
            let stats = cluster_stats(&lattice, &all_open, &from, &to);
            assert!(stats.spanning, "p = 1 must span");
            assert_eq!(stats.max_size, 9);
            assert_eq!(stats.second_moment, 81);
            assert_eq!(stats.n_clusters, 1, "p = 1 leaves one connected cluster");
        }
    }

    #[test]
    fn arbitrary_graphs_work_via_edges() {
        use crate::lattice::graph::BondType;
        // Star graph: center 0 with leaves 1..4; center plus leaves 1, 2 open.
        let lattice = CsrLattice::try_from_edges(
            5,
            vec![
                Bond::new(0, 1, BondType::ChainX, 1.0),
                Bond::new(0, 2, BondType::ChainX, 1.0),
                Bond::new(0, 3, BondType::ChainX, 1.0),
                Bond::new(0, 4, BondType::ChainX, 1.0),
            ],
        )
        .expect("star graph is valid");
        let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Site);
        occupancy
            .site_open
            .copy_from_slice(&[true, true, true, false, false]);
        let stats = cluster_stats(&lattice, &occupancy, &[1], &[2]);
        assert_eq!(
            stats,
            ClusterStats {
                max_size: 3,
                second_moment: 9,
                n_clusters: 1,
                spanning: true,
            }
        );
    }
}
