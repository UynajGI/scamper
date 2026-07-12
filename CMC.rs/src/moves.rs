//! Lattice-spin move and incremental-cache patch types.

use smallvec::SmallVec;

/// Small-vector spin storage used by the compatibility lattice backend.
pub type Spin = SmallVec<[f64; 3]>;

/// Replace one site's spin.
#[derive(Debug, Clone)]
pub struct SiteSpinMove {
    pub site: usize,
    pub spin: Spin,
}

impl SiteSpinMove {
    #[inline]
    pub fn new(site: usize, spin: Spin) -> Self {
        Self { site, spin }
    }
}

/// Replace an arbitrary set of distinct sites atomically.
#[derive(Debug, Clone, Default)]
pub struct BatchSpinMove {
    spin_dim: usize,
    sites: Vec<usize>,
    spins: Vec<f64>,
}

impl BatchSpinMove {
    pub fn new(spin_dim: usize) -> Self {
        assert!(spin_dim > 0, "spin dimension must be positive");
        Self {
            spin_dim,
            sites: Vec::new(),
            spins: Vec::new(),
        }
    }

    pub fn with_capacity(spin_dim: usize, sites: usize) -> Self {
        assert!(spin_dim > 0, "spin dimension must be positive");
        Self {
            spin_dim,
            sites: Vec::with_capacity(sites),
            spins: Vec::with_capacity(sites.saturating_mul(spin_dim)),
        }
    }

    pub fn reset(&mut self, spin_dim: usize) {
        assert!(spin_dim > 0, "spin dimension must be positive");
        self.spin_dim = spin_dim;
        self.sites.clear();
        self.spins.clear();
    }

    pub fn push(&mut self, site: usize, spin: &[f64]) {
        assert_eq!(spin.len(), self.spin_dim, "batch spin dimension mismatch");
        assert!(
            spin.iter().all(|component| component.is_finite()),
            "batch spin contains a non-finite component"
        );
        self.sites.push(site);
        self.spins.extend_from_slice(spin);
    }

    #[inline]
    pub const fn spin_dim(&self) -> usize {
        self.spin_dim
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    #[inline]
    pub fn sites(&self) -> &[usize] {
        &self.sites
    }

    #[inline]
    pub fn spin(&self, index: usize) -> &[f64] {
        let base = index * self.spin_dim;
        &self.spins[base..base + self.spin_dim]
    }
}

/// Cache update for a one-site move.
#[derive(Debug, Clone, Copy, Default)]
pub struct EnergyPatch {
    pub delta_energy: f64,
}

/// Reusable scratch for exact incremental energy of a batch move.
#[derive(Debug, Clone, Default)]
pub struct BatchEnergyWorkspace {
    generation: usize,
    site_stamp: Vec<usize>,
    site_change_index: Vec<usize>,
    edge_stamp: Vec<usize>,
    scratch_spins: Vec<f64>,
}

impl BatchEnergyWorkspace {
    pub const fn new() -> Self {
        Self {
            generation: 0,
            site_stamp: Vec::new(),
            site_change_index: Vec::new(),
            edge_stamp: Vec::new(),
            scratch_spins: Vec::new(),
        }
    }

    /// Start one batch evaluation and index all changed sites.
    ///
    /// Custom lattice Hamiltonians overriding `batch_delta_energy` may reuse
    /// this workspace instead of allocating maps or visited sets per move.
    pub fn prepare(&mut self, n_sites: usize, n_edges: usize, movement: &BatchSpinMove) {
        if self.site_stamp.len() != n_sites {
            self.site_stamp.resize(n_sites, 0);
            self.site_change_index.resize(n_sites, 0);
        }
        if self.edge_stamp.len() != n_edges {
            self.edge_stamp.resize(n_edges, 0);
        }

        if self.generation == usize::MAX {
            self.site_stamp.fill(0);
            self.edge_stamp.fill(0);
            self.generation = 1;
        } else {
            self.generation += 1;
            if self.generation == 0 {
                self.generation = 1;
            }
        }

        for (change_index, &site) in movement.sites().iter().enumerate() {
            assert!(site < n_sites, "batch move site out of range");
            assert_ne!(
                self.site_stamp[site], self.generation,
                "batch move contains duplicate site {site}"
            );
            self.site_stamp[site] = self.generation;
            self.site_change_index[site] = change_index;
        }
    }

    /// Return the movement index for a changed site in the active generation.
    #[inline]
    pub fn change_index(&self, site: usize) -> Option<usize> {
        (self.site_stamp[site] == self.generation).then_some(self.site_change_index[site])
    }

    /// Mark a physical edge for the active generation.
    ///
    /// Returns `true` only on the first visit, allowing both endpoints to be
    /// scanned without double-counting the edge.
    #[inline]
    pub fn mark_edge_once(&mut self, edge_id: usize) -> bool {
        if self.edge_stamp[edge_id] == self.generation {
            false
        } else {
            self.edge_stamp[edge_id] = self.generation;
            true
        }
    }

    /// Materialize a proposed configuration in reusable flat scratch storage.
    pub fn scratch_configuration<'a>(
        &'a mut self,
        spins: &[f64],
        movement: &BatchSpinMove,
    ) -> &'a [f64] {
        self.scratch_spins.clear();
        self.scratch_spins.extend_from_slice(spins);
        for (index, &site) in movement.sites().iter().enumerate() {
            let base = site * movement.spin_dim();
            self.scratch_spins[base..base + movement.spin_dim()]
                .copy_from_slice(movement.spin(index));
        }
        &self.scratch_spins
    }
}

/// Cache patch for an atomic batch move.
#[derive(Debug, Clone, Default)]
pub struct BatchEnergyPatch {
    pub delta_energy: f64,
    pub workspace: BatchEnergyWorkspace,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "duplicate site")]
    fn duplicate_batch_sites_are_rejected() {
        let mut movement = BatchSpinMove::new(1);
        movement.push(0, &[1.0]);
        movement.push(0, &[-1.0]);
        let mut workspace = BatchEnergyWorkspace::default();
        workspace.prepare(1, 0, &movement);
    }
}
