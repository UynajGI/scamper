//! Carlo.rs adapter for the explicit occupation worldline solver.

use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;

use crate::impurity::spin_boson::occupation::model::{
    CavityMode, OccupationModelKind, OccupationSpinBosonModel,
};
use crate::impurity::spin_boson::occupation::worldline::OccupationWorldlineSampler;
use crate::impurity::ImpurityError;

#[derive(Debug)]
pub struct OccupationWorldlineQmc {
    sampler: OccupationWorldlineSampler,
}

impl OccupationWorldlineQmc {
    pub const fn new(sampler: OccupationWorldlineSampler) -> Self {
        Self { sampler }
    }
    pub const fn sampler(&self) -> &OccupationWorldlineSampler {
        &self.sampler
    }
    pub fn sampler_mut(&mut self) -> &mut OccupationWorldlineSampler {
        &mut self.sampler
    }
}

impl MonteCarlo for OccupationWorldlineQmc {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sampler
            .sweep(&mut ctx.rng)
            .unwrap_or_else(|error| panic!("occupation worldline sweep failed: {error}"));
        ctx.record_attempts(self.sampler.slices() as u64);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let obs = self
            .sampler
            .measure()
            .unwrap_or_else(|error| panic!("occupation measurement failed: {error}"));
        ctx.measure("OccupationEnergy", obs.energy);
        ctx.measure("OccupationSigmaZ", obs.sigma_z);
        ctx.measure("OccupationSigmaX", obs.sigma_x);
        ctx.measure("OccupationBosonNumber", obs.total_boson_number);
        ctx.measure("OccupationParity", obs.parity);
        ctx.measure("OccupationReducedSpinPurity", obs.reduced_spin_purity);
        ctx.measure_array("OccupationModeNumber", &obs.mode_occupations);
        ctx.measure_array("OccupationModeNumberSquared", &obs.mode_number_squared);
        ctx.measure_array("OccupationModeFactorialMoment", &obs.mode_factorial_moments);
        ctx.measure_array("OccupationModeG2Zero", &obs.mode_g2_zero);
        ctx.measure_array(
            "OccupationModeNumberCorrelations",
            &obs.mode_cross_correlations,
        );
        ctx.measure_array(
            "OccupationSpinBosonCovarianceZn",
            &obs.spin_boson_covariance_z_n,
        );
        ctx.measure(
            "OccupationWorldlineChangeFraction",
            self.sampler.acceptance_fraction(),
        );
    }

    fn name(&self) -> &'static str {
        "OccupationWorldlineQmc"
    }
}

impl FromParams for OccupationWorldlineQmc {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let beta = required::<f64>(params, "beta")?;
        let slices = params.get::<usize>("slices").unwrap_or(8);
        let kind = model_kind_from_params(params)?;
        let spin_splitting = params
            .get::<f64>("spin_splitting")
            .or_else(|| params.get::<f64>("omega_q"))
            .or_else(|| params.get::<f64>("epsilon"))
            .unwrap_or(0.0);
        let modes = modes_from_params(params).map_err(to_carlo_error)?;
        let model =
            OccupationSpinBosonModel::new(kind, spin_splitting, modes).map_err(to_carlo_error)?;
        let dimension = model.basis().dimension();
        let initial_state = params.get::<usize>("initial_state").unwrap_or_else(|| {
            // Randomize the starting occupation state when none is requested,
            // matching the cluster backend's initial-spin convention.
            rng.random_range(0..dimension)
        });
        let sampler = OccupationWorldlineSampler::new(model, beta, slices, initial_state)
            .map_err(to_carlo_error)?;
        Ok(Self::new(sampler))
    }

    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let beta = required::<f64>(params, "beta")?;
        if !beta.is_finite() || beta <= 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "beta".into(),
                reason: format!("must be finite and positive, got {beta}"),
            });
        }
        // At least one mode must be supplied either through the single-mode
        // shortcut (`omega0`/`g`/`cutoff`) or through the indexed list
        // (`omega_0`/`g_0`/`cutoff_0`, ...). The detailed validation happens
        // in `modes_from_params`; here we only reject the case where neither
        // convention supplies any mode.
        if modes_from_params(params)
            .map_err(to_carlo_error)?
            .is_empty()
        {
            return Err(CarloError::InvalidConfig {
                field: "modes".into(),
                reason: "at least one cavity mode is required (set `omega0`/`g`/`cutoff` or \
                         indexed `omega_<i>`/`g_<i>`/`cutoff_<i>`)"
                    .into(),
            });
        }
        Ok(())
    }
}

fn model_kind_from_params(params: &Params) -> Result<OccupationModelKind, CarloError> {
    let kind = params
        .get::<String>("kind")
        .unwrap_or_else(|| "rabi".into())
        .to_ascii_lowercase();
    match kind.as_str() {
        "rabi" => Ok(OccupationModelKind::Rabi),
        "jaynes-cummings" | "jaynes_cummings" | "jc" => Ok(OccupationModelKind::JaynesCummings),
        other => Err(CarloError::InvalidConfig {
            field: "kind".into(),
            reason: format!("unsupported occupation model `{other}` (expected `rabi` or `jc`)"),
        }),
    }
}

/// Parse cavity modes from the flat `Params` map.
///
/// Two conventions are supported:
/// 1. Single-mode shortcut: `omega0`, `g`, `cutoff`.
/// 2. Indexed list: `omega_<i>`, `g_<i>`, `cutoff_<i>` for `i = 0, 1, 2, ...`.
///
/// If the single-mode shortcut is present it takes precedence and any indexed
/// entries are ignored, matching the cluster/wormhole single-mode convention.
fn modes_from_params(params: &Params) -> Result<Vec<CavityMode>, ImpurityError> {
    if params.contains("omega0") || params.contains("cutoff") {
        let omega = params.get::<f64>("omega0").unwrap_or(1.0);
        let coupling = params
            .get::<f64>("g")
            .or_else(|| params.get::<f64>("g_sigma"))
            .unwrap_or(0.0);
        let cutoff = params.get::<usize>("cutoff").unwrap_or(4);
        return Ok(vec![CavityMode::new(omega, coupling, cutoff)?]);
    }

    let mut modes = Vec::new();
    let mut index = 0;
    loop {
        let omega_key = format!("omega_{index}");
        if !params.contains(&omega_key) {
            break;
        }
        let cutoff_key = format!("cutoff_{index}");
        let g_key = format!("g_{index}");
        let omega = params.get::<f64>(&omega_key).unwrap_or(1.0);
        let coupling = params.get::<f64>(&g_key).unwrap_or(0.0);
        let cutoff = params.get::<usize>(&cutoff_key).unwrap_or(4);
        modes.push(CavityMode::new(omega, coupling, cutoff)?);
        index += 1;
    }
    Ok(modes)
}

fn required<T>(params: &Params, key: &str) -> Result<T, CarloError>
where
    T: std::str::FromStr,
{
    params
        .get::<T>(key)
        .ok_or_else(|| CarloError::InvalidConfig {
            field: key.into(),
            reason: format!("parameter `{key}` is required"),
        })
}

fn to_carlo_error(error: ImpurityError) -> CarloError {
    CarloError::InvalidConfig {
        field: "impurity.occupation".into(),
        reason: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;

    use super::*;

    #[test]
    fn from_params_builds_single_mode_rabi_solver() {
        let mut params = Params::new();
        params.set("beta", 2.0);
        params.set("kind", "rabi");
        params.set("omega_q", 0.8);
        params.set("omega0", 1.2);
        params.set("g", 0.3);
        params.set("cutoff", 5);
        params.set("slices", 6);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        let simulation =
            OccupationWorldlineQmc::from_params(&params, &mut rng).expect("simulation");
        assert_eq!(
            simulation.sampler().model().kind(),
            OccupationModelKind::Rabi
        );
        assert!((simulation.sampler().model().spin_splitting() - 0.8).abs() < f64::EPSILON);
        assert_eq!(simulation.sampler().model().modes().len(), 1);
        assert_eq!(simulation.sampler().slices(), 6);
    }

    #[test]
    fn from_params_builds_jaynes_cummings_solver_via_jc_alias() {
        let mut params = Params::new();
        params.set("beta", 1.5);
        params.set("kind", "jc");
        params.set("spin_splitting", 1.0);
        params.set("omega0", 1.0);
        params.set("g", 0.25);
        params.set("cutoff", 6);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2);
        let simulation =
            OccupationWorldlineQmc::from_params(&params, &mut rng).expect("simulation");
        assert_eq!(
            simulation.sampler().model().kind(),
            OccupationModelKind::JaynesCummings
        );
    }

    #[test]
    fn from_params_builds_multimode_solver_from_indexed_keys() {
        let mut params = Params::new();
        params.set("beta", 1.2);
        params.set("spin_splitting", 0.7);
        params.set("omega_0", 1.0);
        params.set("g_0", 0.15);
        params.set("cutoff_0", 4);
        params.set("omega_1", 1.4);
        params.set("g_1", 0.08);
        params.set("cutoff_1", 3);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        let simulation =
            OccupationWorldlineQmc::from_params(&params, &mut rng).expect("simulation");
        let modes = simulation.sampler().model().modes();
        assert_eq!(modes.len(), 2);
        assert!((modes[0].omega - 1.0).abs() < f64::EPSILON);
        assert!((modes[1].omega - 1.4).abs() < f64::EPSILON);
    }

    #[test]
    fn from_params_rejects_missing_beta() {
        let mut params = Params::new();
        params.set("omega0", 1.0);
        params.set("cutoff", 4);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        let error = OccupationWorldlineQmc::from_params(&params, &mut rng)
            .expect_err("missing beta must fail");
        assert!(matches!(error, CarloError::InvalidConfig { field, .. } if field == "beta"));
    }

    #[test]
    fn from_params_rejects_no_modes() {
        let mut params = Params::new();
        params.set("beta", 1.0);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        let error =
            OccupationWorldlineQmc::from_params(&params, &mut rng).expect_err("no modes must fail");
        assert!(matches!(error, CarloError::InvalidConfig { field, .. } if field == "modes"));
    }

    #[test]
    fn from_params_rejects_unknown_kind() {
        let mut params = Params::new();
        params.set("beta", 1.0);
        params.set("kind", "bogus");
        params.set("omega0", 1.0);
        params.set("cutoff", 4);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(6);
        let error = OccupationWorldlineQmc::from_params(&params, &mut rng)
            .expect_err("unknown kind must fail");
        assert!(matches!(error, CarloError::InvalidConfig { field, .. } if field == "kind"));
    }

    #[test]
    fn validate_params_accepts_single_mode_shortcut() {
        let mut params = Params::new();
        params.set("beta", 1.0);
        params.set("omega0", 1.0);
        params.set("cutoff", 4);
        assert!(OccupationWorldlineQmc::validate_params(&params).is_ok());
    }

    #[test]
    fn validate_params_rejects_non_positive_beta() {
        let mut params = Params::new();
        params.set("beta", 0.0);
        params.set("omega0", 1.0);
        params.set("cutoff", 4);
        assert!(OccupationWorldlineQmc::validate_params(&params).is_err());
    }
}
