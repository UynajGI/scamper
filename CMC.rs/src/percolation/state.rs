//! Occupancy configurations for site and bond percolation.

use crate::lattice::graph::CsrLattice;
use rand::{Rng, RngExt};

/// Percolation mode: occupied sites or occupied bonds.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PercolationMode {
    #[default]
    Site,
    Bond,
}

impl PercolationMode {
    /// Parse the `mode` parameter label.
    pub fn from_label(label: &str) -> Result<Self, String> {
        match label {
            "site" => Ok(Self::Site),
            "bond" => Ok(Self::Bond),
            other => Err(format!(
                "unknown percolation mode `{other}`, expected `site` or `bond`"
            )),
        }
    }

    /// Parameter label of this mode.
    pub const fn as_label(self) -> &'static str {
        match self {
            Self::Site => "site",
            Self::Bond => "bond",
        }
    }
}

/// Independent occupancy configuration of one percolation sample.
///
/// Site percolation uses `site_open`; bond percolation uses `bond_open`
/// (indexed by physical edge id). The irrelevant array stays all-`false`
/// and is ignored by [`super::cluster_stats`].
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

    /// Redraw every relevant element independently open with probability `p`.
    pub fn resample<R: Rng + ?Sized>(&mut self, p: f64, rng: &mut R) {
        match self.mode {
            PercolationMode::Site => {
                for open in &mut self.site_open {
                    *open = rng.random::<f64>() < p;
                }
            }
            PercolationMode::Bond => {
                for open in &mut self.bond_open {
                    *open = rng.random::<f64>() < p;
                }
            }
        }
    }

    /// Number of occupied elements: open sites (site mode) or open bonds
    /// (bond mode).
    pub fn occupied(&self) -> usize {
        match self.mode {
            PercolationMode::Site => self.site_open.iter().filter(|&&open| open).count(),
            PercolationMode::Bond => self.bond_open.iter().filter(|&&open| open).count(),
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
        for mode in [PercolationMode::Site, PercolationMode::Bond] {
            assert_eq!(PercolationMode::from_label(mode.as_label()), Ok(mode));
        }
        assert!(PercolationMode::from_label("edge").is_err());
    }

    #[test]
    fn new_configuration_is_all_closed() {
        let lattice = build_square(2, 2, false);
        for mode in [PercolationMode::Site, PercolationMode::Bond] {
            let occupancy = OccupancyState::new(&lattice, mode);
            assert_eq!(occupancy.occupied(), 0);
            assert_eq!(occupancy.site_open.len(), 4);
            assert_eq!(occupancy.bond_open.len(), lattice.n_edges());
        }
    }

    #[test]
    fn resample_extremes_are_deterministic() {
        let lattice = build_square(2, 2, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        for mode in [PercolationMode::Site, PercolationMode::Bond] {
            let mut occupancy = OccupancyState::new(&lattice, mode);
            occupancy.resample(0.0, &mut rng);
            assert_eq!(occupancy.occupied(), 0, "p = 0 must occupy nothing");
            occupancy.resample(1.0, &mut rng);
            let total = match mode {
                PercolationMode::Site => lattice.n_sites,
                PercolationMode::Bond => lattice.n_edges(),
            };
            assert_eq!(occupancy.occupied(), total, "p = 1 must occupy everything");
        }
    }

    #[test]
    fn resample_respects_independence_and_bounds() {
        let lattice = build_square(8, 8, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        let mut occupancy = OccupancyState::new(&lattice, PercolationMode::Site);
        let trials = 200;
        let mut total_open = 0usize;
        for _ in 0..trials {
            occupancy.resample(0.5, &mut rng);
            let open = occupancy.occupied();
            assert!(open <= lattice.n_sites);
            total_open += open;
        }
        // 64 sites x 200 trials at p = 0.5: mean 6400, sigma ~ 39.6.
        let z = (total_open as f64 - 6400.0) / (64.0_f64 * 0.5 * 200.0_f64).sqrt();
        assert!(z.abs() < 4.0, "occupancy rate deviates: z = {z:.2}");
    }
}
