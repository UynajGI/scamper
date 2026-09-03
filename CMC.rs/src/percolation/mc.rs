//! Carlo.rs adapter for independent-sample percolation.

use super::state::{OccupancyState, PercolationMode};
use super::{cluster_stats, ClusterStats};
use crate::classical_mc::{build_lattice_from_params, parse_param};
use crate::lattice::graph::CsrLattice;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Scheduler-ready independent-sample percolation.
///
/// Every sweep draws a fresh occupancy configuration. Samples are i.i.d., so
/// there is no Markov-chain equilibration — set `thermalization_sweeps = 0`
/// and read every sweep as one independent sample.
///
/// Modes (`mode` parameter):
///
/// - `"site"` (default): sites open independently with probability `p`.
/// - `"bond"`: bonds open independently with probability `p`; clusters run
///   over all sites, so isolated sites are singleton clusters.
/// - `"site-bond"`: sites open with probability `p_site` and bonds with
///   probability `p_bond`; a bond connects only when it is open **and** both
///   endpoint sites are open.
///
/// Measured observables: `MaxCluster` (largest cluster in sites),
/// `SecondMoment` (`sum(s_i^2)`), `NClusters` and `Spanning` (0/1 indicator
/// whose mean is the crossing probability). Pure modes additionally measure
/// `Occupied` (open sites or open bonds); `"site-bond"` measures
/// `OccupiedSites` and `OccupiedBonds` instead.
///
/// Crossing is tested between the `spanning_from` and `spanning_to` site
/// sets. Defaults: `square` lattices use the left vs. right column
/// (row-major indexing, open or periodic), `chain` uses the two end sites;
/// every other lattice type requires explicit comma-separated
/// `spanning_from`/`spanning_to` site lists.
pub struct PercolationMC {
    lattice: CsrLattice,
    occupancy: OccupancyState,
    p_site: f64,
    p_bond: f64,
    from: Vec<usize>,
    to: Vec<usize>,
    last: ClusterStats,
}

impl PercolationMC {
    /// Build a percolation ensemble on `lattice`. Pure modes ignore the
    /// probability of the unsampled element kind.
    pub fn new(
        lattice: CsrLattice,
        mode: PercolationMode,
        p_site: f64,
        p_bond: f64,
        from: Vec<usize>,
        to: Vec<usize>,
    ) -> Result<Self, String> {
        for (name, probability) in [("p_site", p_site), ("p_bond", p_bond)] {
            if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
                return Err(format!(
                    "occupation probability `{name}` must be finite within [0, 1], got {probability}"
                ));
            }
        }
        if from.is_empty() || to.is_empty() {
            return Err("spanning site sets must not be empty".to_string());
        }
        for &site in from.iter().chain(to.iter()) {
            if site >= lattice.n_sites {
                return Err(format!(
                    "spanning site {site} out of range for {} sites",
                    lattice.n_sites
                ));
            }
        }
        let occupancy = OccupancyState::new(&lattice, mode);
        Ok(Self {
            lattice,
            occupancy,
            p_site,
            p_bond,
            from,
            to,
            last: ClusterStats {
                max_size: 0,
                second_moment: 0,
                n_clusters: 0,
                spanning: false,
            },
        })
    }

    /// Underlying lattice.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// Current occupancy configuration.
    pub const fn occupancy(&self) -> &OccupancyState {
        &self.occupancy
    }

    /// Site occupation probability (`p_site`).
    pub const fn site_probability(&self) -> f64 {
        self.p_site
    }

    /// Bond occupation probability (`p_bond`).
    pub const fn bond_probability(&self) -> f64 {
        self.p_bond
    }

    /// Spanning test site sets `(from, to)` used by `Spanning`.
    pub fn spanning_sets(&self) -> (&[usize], &[usize]) {
        (&self.from, &self.to)
    }

    /// Cluster statistics of the most recent measurement.
    pub const fn last_stats(&self) -> ClusterStats {
        self.last
    }
}

impl MonteCarlo for PercolationMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.occupancy
            .resample(self.p_site, self.p_bond, &mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let stats = cluster_stats(&self.lattice, &self.occupancy, &self.from, &self.to);
        match self.occupancy.mode {
            PercolationMode::SiteBond => {
                ctx.measure("OccupiedSites", self.occupancy.occupied_sites() as f64);
                ctx.measure("OccupiedBonds", self.occupancy.occupied_bonds() as f64);
            }
            _ => ctx.measure("Occupied", self.occupancy.occupied() as f64),
        }
        ctx.measure("MaxCluster", stats.max_size as f64);
        ctx.measure("SecondMoment", stats.second_moment as f64);
        ctx.measure("NClusters", stats.n_clusters as f64);
        ctx.measure("Spanning", u8::from(stats.spanning) as f64);
        self.last = stats;
    }

    fn name(&self) -> &'static str {
        "Percolation"
    }
}

/// Parse and range-check the occupation probabilities demanded by `mode`.
///
/// Pure modes take `p` and reject `p_site`/`p_bond`; the mixed mode takes
/// `p_site`/`p_bond` and rejects `p` as ambiguous.
fn parse_probabilities(params: &Params, mode: PercolationMode) -> Result<(f64, f64), CarloError> {
    let parse_probability = |name: &str| -> Result<f64, CarloError> {
        let value = parse_param::<f64>(params, name)?
            .ok_or_else(|| invalid(name, "missing required occupation probability"))?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(invalid(
                name,
                format!("must be finite within [0, 1], got {value}"),
            ));
        }
        Ok(value)
    };
    let reject_unused = |name: &str| -> Result<(), CarloError> {
        if params.contains(name) {
            return Err(invalid(
                name,
                format!(
                    "is only valid with `mode = \"site-bond\"`, not `{}`",
                    mode.as_label()
                ),
            ));
        }
        Ok(())
    };
    match mode {
        PercolationMode::Site | PercolationMode::Bond => {
            reject_unused("p_site")?;
            reject_unused("p_bond")?;
            let p = parse_probability("p")?;
            Ok((p, p))
        }
        PercolationMode::SiteBond => {
            if params.contains("p") {
                return Err(invalid(
                    "p",
                    "is ambiguous with `mode = \"site-bond\"`; use `p_site`/`p_bond`",
                ));
            }
            Ok((parse_probability("p_site")?, parse_probability("p_bond")?))
        }
    }
}

impl FromParams for PercolationMC {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let mode = parse_mode(params)?;
        parse_probabilities(params, mode)?;
        parse_index_list(params, "spanning_from")?;
        parse_index_list(params, "spanning_to")?;
        Ok(())
    }

    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let mode = parse_mode(params)?;
        let (p_site, p_bond) = parse_probabilities(params, mode)?;
        let pbc = parse_param::<bool>(params, "pbc")?.unwrap_or(false);
        let lattice = build_lattice_from_params(params, pbc)?;
        let (from, to) = resolve_spanning_sets(params, &lattice)?;
        Self::new(lattice, mode, p_site, p_bond, from, to)
            .map_err(|error| invalid("percolation", error))
    }
}

/// Parse the `mode` parameter (default `"site"`).
fn parse_mode(params: &Params) -> Result<PercolationMode, CarloError> {
    match parse_param::<String>(params, "mode")? {
        Some(label) => PercolationMode::from_label(&label).map_err(|error| invalid("mode", error)),
        None => Ok(PercolationMode::Site),
    }
}

/// Explicit `spanning_from`/`spanning_to` lists, or structural defaults.
///
/// Defaults: `square` (and a 2D `hypercubic` spec) uses the left vs. right
/// column under row-major site indexing, `chain` uses the two end sites.
/// Every other lattice type requires the explicit sets because a meaningful
/// crossing direction is not implied by the graph alone.
fn resolve_spanning_sets(
    params: &Params,
    lattice: &CsrLattice,
) -> Result<(Vec<usize>, Vec<usize>), CarloError> {
    let from = parse_index_list(params, "spanning_from")?;
    let to = parse_index_list(params, "spanning_to")?;
    match (from, to) {
        (Some(from), Some(to)) => return Ok((from, to)),
        (None, None) => {}
        (Some(_), None) | (None, Some(_)) => {
            return Err(invalid(
                "spanning_to",
                "`spanning_from` and `spanning_to` must be given together",
            ));
        }
    }

    let inferred = if params.contains("Lx") {
        "hypercubic"
    } else {
        "chain"
    };
    let lattice_type = params
        .get::<String>("lattice_type")
        .unwrap_or_else(|| inferred.to_string())
        .to_ascii_lowercase();
    match lattice_type.as_str() {
        "chain" => Ok((vec![0], vec![lattice.n_sites - 1])),
        "square" => square_columns(params, lattice.n_sites),
        "hypercubic" if params.contains("Lx") && !params.contains("Lz") => {
            square_columns(params, lattice.n_sites)
        }
        other => Err(invalid(
            "spanning_from",
            format!("lattice type `{other}` requires explicit `spanning_from`/`spanning_to`"),
        )),
    }
}

/// Left vs. right column of an `lx * ly` row-major square grid, mirroring
/// `build_lattice_from_params` dimension defaults.
fn square_columns(params: &Params, n_sites: usize) -> Result<(Vec<usize>, Vec<usize>), CarloError> {
    let lx = parse_param::<usize>(params, "Lx")?.unwrap_or(4);
    let ly = parse_param::<usize>(params, "Ly")?.unwrap_or(lx);
    if lx * ly != n_sites {
        return Err(invalid(
            "Lx",
            format!("{lx}x{ly} does not match the {n_sites} built sites"),
        ));
    }
    Ok((
        (0..ly).map(|row| row * lx).collect(),
        (0..ly).map(|row| row * lx + lx - 1).collect(),
    ))
}

/// Parse a comma-separated site-index list parameter.
fn parse_index_list(params: &Params, name: &str) -> Result<Option<Vec<usize>>, CarloError> {
    let Some(raw) = parse_param::<String>(params, name)? else {
        return Ok(None);
    };
    if raw.trim().is_empty() {
        return Err(invalid(name, "site list must not be empty"));
    }
    let mut sites = Vec::new();
    for token in raw.split(',') {
        let token = token.trim();
        let site = token
            .parse::<usize>()
            .map_err(|_| invalid(name, format!("cannot parse site index `{token}`")))?;
        sites.push(site);
    }
    Ok(Some(sites))
}

fn invalid(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn rng() -> Xoshiro256PlusPlus {
        Xoshiro256PlusPlus::from_seed([0; 32])
    }

    fn square_params() -> Params {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 2);
        params.set("Ly", 2);
        params.set("p", 0.5);
        params
    }

    #[test]
    fn from_params_builds_square_with_column_defaults() {
        let params = square_params();
        let mc = PercolationMC::from_params(&params, &mut rng()).expect("square params are valid");
        assert_eq!(mc.lattice().n_sites, 4);
        assert_eq!(mc.site_probability(), 0.5);
        assert_eq!(mc.bond_probability(), 0.5);
        assert_eq!(mc.occupancy().mode, PercolationMode::Site);
        // Column defaults for the 2x2 grid: {0, 2} vs. {1, 3}.
        assert_eq!(mc.spanning_sets(), (&[0, 2][..], &[1, 3][..]));
    }

    #[test]
    fn from_params_accepts_bond_mode_and_explicit_sets() {
        let mut params = square_params();
        params.set("mode", "bond");
        params.set("spanning_from", "0, 1");
        params.set("spanning_to", "2, 3");
        let mc = PercolationMC::from_params(&params, &mut rng()).expect("bond params are valid");
        assert_eq!(mc.occupancy().mode, PercolationMode::Bond);
        assert_eq!(mc.spanning_sets(), (&[0, 1][..], &[2, 3][..]));
    }

    #[test]
    fn from_params_builds_site_bond_mode() {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 2);
        params.set("Ly", 2);
        params.set("mode", "site-bond");
        params.set("p_site", 0.7);
        params.set("p_bond", 0.4);
        let mc = PercolationMC::from_params(&params, &mut rng()).expect("mixed params are valid");
        assert_eq!(mc.occupancy().mode, PercolationMode::SiteBond);
        assert_eq!(mc.site_probability(), 0.7);
        assert_eq!(mc.bond_probability(), 0.4);
    }

    #[test]
    fn from_params_rejects_invalid_probability_and_mode() {
        let mut params = square_params();
        params.set("p", 1.5);
        assert!(PercolationMC::validate_params(&params).is_err());
        params.set("p", 0.5);
        params.set("mode", "edge");
        assert!(PercolationMC::validate_params(&params).is_err());
        let missing_p = Params::new();
        assert!(PercolationMC::validate_params(&missing_p).is_err());
    }

    #[test]
    fn from_params_rejects_mismatched_probability_parameters() {
        // `p_site`/`p_bond` are only meaningful in mixed mode.
        let mut params = square_params();
        params.set("p_site", 0.5);
        assert!(PercolationMC::validate_params(&params).is_err());
        let mut params = square_params();
        params.set("p_bond", 0.5);
        assert!(PercolationMC::validate_params(&params).is_err());

        // `p` is ambiguous in mixed mode.
        let mut params = Params::new();
        params.set("mode", "site-bond");
        params.set("p", 0.5);
        params.set("p_site", 0.5);
        params.set("p_bond", 0.5);
        assert!(PercolationMC::validate_params(&params).is_err());

        // Mixed mode requires both probabilities.
        let mut params = Params::new();
        params.set("mode", "site-bond");
        params.set("p_site", 0.5);
        assert!(PercolationMC::validate_params(&params).is_err());
        params.set("p_bond", 1.5);
        assert!(PercolationMC::validate_params(&params).is_err());
    }

    #[test]
    fn from_params_rejects_one_sided_spanning_sets() {
        let mut params = square_params();
        params.set("spanning_from", "0");
        assert!(PercolationMC::from_params(&params, &mut rng()).is_err());
    }

    #[test]
    fn from_params_rejects_out_of_range_spanning_sites() {
        let mut params = square_params();
        params.set("spanning_from", "0");
        params.set("spanning_to", "4");
        assert!(PercolationMC::from_params(&params, &mut rng()).is_err());
    }

    #[test]
    fn chain_defaults_to_end_sites() {
        let mut params = Params::new();
        params.set("lattice_type", "chain");
        params.set("L", 5);
        params.set("p", 0.5);
        let mc = PercolationMC::from_params(&params, &mut rng()).expect("chain params are valid");
        assert_eq!(mc.spanning_sets(), (&[0][..], &[4][..]));
    }

    #[test]
    fn non_trivial_lattices_require_explicit_sets() {
        let mut params = Params::new();
        params.set("lattice_type", "triangular");
        params.set("Lx", 2);
        params.set("Ly", 2);
        params.set("pbc", true);
        params.set("p", 0.5);
        assert!(PercolationMC::from_params(&params, &mut rng()).is_err());
        params.set("spanning_from", "0");
        params.set("spanning_to", "3");
        let mc =
            PercolationMC::from_params(&params, &mut rng()).expect("explicit sets make it valid");
        assert_eq!(mc.lattice().n_sites, 4);
        assert_eq!(mc.spanning_sets(), (&[0][..], &[3][..]));
    }
}
