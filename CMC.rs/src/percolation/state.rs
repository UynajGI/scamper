//! Occupancy configurations for site, bond and mixed site-bond percolation.

use crate::lattice::graph::CsrLattice;
use rand::{Rng, RngExt};

/// Percolation mode: occupied sites, occupied bonds, or both independently.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PercolationMode {
    #[default]
    Site,
    Bond,
    /// Mixed site-bond: a bond connects only when the bond itself and both
    /// endpoint sites are open.
    SiteBond,
}

impl PercolationMode {
    /// Parse the `mode` parameter label.
    pub fn from_label(label: &str) -> Result<Self, String> {
        match label {
            "site" => Ok(Self::Site),
            "bond" => Ok(Self::Bond),
            "site-bond" => Ok(Self::SiteBond),
            other => Err(format!(
                "unknown percolation mode `{other}`, expected `site`, `bond` or `site-bond`"
            )),
        }
    }

    /// Parameter label of this mode.
    pub const fn as_label(self) -> &'static str {
        match self {
            Self::Site => "site",
            Self::Bond => "bond",
            Self::SiteBond => "site-bond",
        }
    }

    /// Whether this mode samples site occupation.
    pub const fn samples_sites(self) -> bool {
        !matches!(self, Self::Bond)
    }

    /// Whether this mode samples bond occupation.
    pub const fn samples_bonds(self) -> bool {
        !matches!(self, Self::Site)
    }
}

/// Independent occupancy configuration of one percolation sample.
///
/// Site percolation uses `site_open`; bond percolation uses `bond_open`
/// (indexed by physical edge id); mixed site-bond percolation uses both.
/// Arrays irrelevant to the mode stay all-`false` and are ignored by
/// [`super::cluster_stats`].
#[derive(Debug, Clone)]
pub struct OccupancyState {
    pub mode: PercolationMode,
    /// Site occupation, one flag per lattice site.
    pub site_open: Vec<bool>,
    /// Bond occupation, one flag per physical edge.
    pub bond_open: Vec<bool>,
}

impl OccupancyState {
    /// All-closed configuration for `lattice`.
    pub fn new(lattice: &CsrLattice, mode: PercolationMode) -> Self {
        Self {
            mode,
            site_open: vec![false; lattice.n_sites],
            bond_open: vec![false; lattice.n_edges()],
        }
    }

    /// Redraw every sampled element independently open: sites with
    /// probability `p_site`, bonds with probability `p_bond` (pure modes
    /// ignore the irrelevant probability).
    pub fn resample<R: Rng + ?Sized>(&mut self, p_site: f64, p_bond: f64, rng: &mut R) {
        if self.mode.samples_sites() {
            for open in &mut self.site_open {
                *open = rng.random::<f64>() < p_site;
            }
        }
        if self.mode.samples_bonds() {
            for open in &mut self.bond_open {
                *open = rng.random::<f64>() < p_bond;
            }
        }
    }

    /// Number of open sites.
    pub fn occupied_sites(&self) -> usize {
        self.site_open.iter().filter(|&&open| open).count()
    }

    /// Number of open bonds.
    pub fn occupied_bonds(&self) -> usize {
        self.bond_open.iter().filter(|&&open| open).count()
    }

    /// Mode-relevant occupied element count: open sites (site mode), open
    /// bonds (bond mode) or both summed (site-bond mode).
    pub fn occupied(&self) -> usize {
        match self.mode {
            PercolationMode::Site => self.occupied_sites(),
            PercolationMode::Bond => self.occupied_bonds(),
            PercolationMode::SiteBond => self.occupied_sites() + self.occupied_bonds(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::graph::build_square;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn mode_labels_round_trip() {
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            assert_eq!(PercolationMode::from_label(mode.as_label()), Ok(mode));
        }
        assert!(PercolationMode::from_label("edge").is_err());
    }

    #[test]
    fn new_configuration_is_all_closed() {
        let lattice = build_square(2, 2, false);
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            let occupancy = OccupancyState::new(&lattice, mode);
            assert_eq!(occupancy.occupied(), 0);
            assert_eq!(occupancy.occupied_sites(), 0);
            assert_eq!(occupancy.occupied_bonds(), 0);
            assert_eq!(occupancy.site_open.len(), 4);
            assert_eq!(occupancy.bond_open.len(), lattice.n_edges());
        }
    }

    #[test]
    fn resample_extremes_are_deterministic() {
        let lattice = build_square(2, 2, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            let mut occupancy = OccupancyState::new(&lattice, mode);
            occupancy.resample(0.0, 0.0, &mut rng);
            assert_eq!(occupancy.occupied(), 0, "p = 0 must occupy nothing");
            occupancy.resample(1.0, 1.0, &mut rng);
            let sampled = match mode {
                PercolationMode::Site => lattice.n_sites,
                PercolationMode::Bond => lattice.n_edges(),
                PercolationMode::SiteBond => lattice.n_sites + lattice.n_edges(),
            };
            assert_eq!(
                occupancy.occupied(),
                sampled,
                "p = 1 must occupy everything"
            );
        }
    }

    #[test]
    fn mixed_resample_ignores_unused_probability() {
        // Bond mode must sample bonds even with a degenerate p_site, and
        // vice versa; mixed mode samples both arrays.
        let lattice = build_square(2, 2, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(13);

        let mut bond = OccupancyState::new(&lattice, PercolationMode::Bond);
        bond.resample(0.0, 1.0, &mut rng);
        assert_eq!(bond.occupied_bonds(), lattice.n_edges());
        assert_eq!(bond.occupied_sites(), 0);

        let mut site = OccupancyState::new(&lattice, PercolationMode::Site);
        site.resample(1.0, 0.0, &mut rng);
        assert_eq!(site.occupied_sites(), lattice.n_sites);
        assert_eq!(site.occupied_bonds(), 0);

        let mut mixed = OccupancyState::new(&lattice, PercolationMode::SiteBond);
        mixed.resample(1.0, 1.0, &mut rng);
        assert_eq!(mixed.occupied_sites(), lattice.n_sites);
        assert_eq!(mixed.occupied_bonds(), lattice.n_edges());
    }

    #[test]
    fn resample_respects_independence_and_bounds() {
        let lattice = build_square(8, 8, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Site);
        let trials = 200;
        let mut total_open = 0usize;
        for _ in 0..trials {
            occupancy.resample(0.5, 0.0, &mut rng);
            let open = occupancy.occupied();
            assert!(open <= lattice.n_sites);
            total_open += open;
        }
        // 64 sites x 200 trials at p = 0.5: mean 6400, sigma ~ 39.6.
        let z = (total_open as f64 - 6400.0) / (64.0_f64 * 0.5 * 200.0_f64).sqrt();
        assert!(z.abs() < 4.0, "occupancy rate deviates: z = {z:.2}");
    }
}
