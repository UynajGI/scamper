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
