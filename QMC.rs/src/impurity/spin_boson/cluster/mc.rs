//! Carlo.rs adapter for the longitudinal spin-boson continuous-time cluster
//! solver.

use carlo_rs::{CarloError, Context, Evaluator, FromParams, MonteCarlo, Params};
use rand::RngExt;

use crate::impurity::core::estimator::register_connected_susceptibility;
use crate::impurity::spin_boson::bath::{Bath, PowerLawBath, SingleModeBath, TabulatedBath};
use crate::impurity::spin_boson::cluster::cluster_builder::ContinuousTimeClusterEngine;
use crate::impurity::spin_boson::cluster::retarded_bonds::LongitudinalSpinBosonModel;
use crate::impurity::spin_boson::cluster::segments::LongitudinalWorldline;
use crate::impurity::ImpurityError;

/// Scalar measurements of a longitudinal continuous-time worldline.
#[derive(Debug, Clone, PartialEq)]
pub struct LongitudinalClusterObservables {
    pub magnetization_sigma_z: f64,
    pub magnetization_s_z: f64,
    pub magnetization_sigma_z_squared: f64,
    pub magnetization_sigma_z_fourth: f64,
    pub correlation_sigma_z: Vec<f64>,
    pub correlation_s_z: Vec<f64>,
    pub kink_count: f64,
    /// Exact expansion-order estimator for `-(Delta/2) sigma_x`.
    pub transverse_field_energy: f64,
    /// Configuration estimator of `(epsilon/2) sigma_z`.
    pub longitudinal_field_energy: f64,
}

/// Measure exact diagonal observables of a piecewise-constant worldline.
pub fn measure_cluster_observables(
    worldline: &LongitudinalWorldline,
    model: &LongitudinalSpinBosonModel,
    correlation_bins: usize,
) -> LongitudinalClusterObservables {
    let bins = correlation_bins.max(1);
    let magnetization_sigma_z = worldline.integrated_sigma_z();
    let correlation_sigma_z = (0..bins)
        .map(|index| {
            let delta = worldline.beta() * index as f64 / bins as f64;
            worldline.correlation_sigma_z(delta)
        })
        .collect::<Vec<_>>();
    let correlation_s_z = correlation_sigma_z
        .iter()
        .map(|value| 0.25 * value)
        .collect::<Vec<_>>();
    let kink_count = worldline.kink_count() as f64;
    LongitudinalClusterObservables {
        magnetization_sigma_z,
        magnetization_s_z: 0.5 * magnetization_sigma_z,
        magnetization_sigma_z_squared: magnetization_sigma_z * magnetization_sigma_z,
        magnetization_sigma_z_fourth: magnetization_sigma_z.powi(4),
        correlation_sigma_z,
        correlation_s_z,
        kink_count,
        transverse_field_energy: -kink_count / worldline.beta(),
        longitudinal_field_energy: 0.5 * model.bias() * magnetization_sigma_z,
    }
}

/// Register jackknife-safe connected susceptibilities for direct and
/// cluster-improved longitudinal estimators.
pub fn register_cluster_evaluables(evaluator: &mut Evaluator, beta: f64) -> Result<(), CarloError> {
    register_connected_susceptibility(
        evaluator,
        "ClusterChiSigmaZConnected",
        "ClusterMagnetizationSigmaZ",
        "ClusterM2SigmaZ",
        beta,
    )?;
    register_connected_susceptibility(
        evaluator,
        "ClusterChiSzConnected",
        "ClusterMagnetizationSz",
        "ClusterM2Sz",
        beta,
    )?;
    register_connected_susceptibility(
        evaluator,
        "ClusterImprovedChiSigmaZConnected",
        "ClusterImprovedMagnetizationSigmaZ",
        "ClusterImprovedM2SigmaZ",
        beta,
    )?;
    register_connected_susceptibility(
        evaluator,
        "ClusterImprovedChiSzConnected",
        "ClusterImprovedMagnetizationSz",
        "ClusterImprovedM2Sz",
        beta,
    )
}

/// Runnable longitudinal spin-boson continuous-time cluster simulation.
pub struct LongitudinalSpinBosonClusterQmc {
    worldline: LongitudinalWorldline,
    engine: ContinuousTimeClusterEngine,
    correlation_bins: usize,
}

impl LongitudinalSpinBosonClusterQmc {
    /// Construct from an already-built longitudinal model.
    pub fn new(
        model: LongitudinalSpinBosonModel,
        beta: f64,
        initial_spin: i8,
        correlation_bins: usize,
    ) -> Result<Self, ImpurityError> {
        Ok(Self {
            worldline: LongitudinalWorldline::new(beta, initial_spin)?,
            engine: ContinuousTimeClusterEngine::new(model),
            correlation_bins: correlation_bins.max(1),
        })
    }

    /// Current worldline.
    #[inline]
    pub const fn worldline(&self) -> &LongitudinalWorldline {
        &self.worldline
    }

    /// Cluster update engine.
    #[inline]
    pub const fn engine(&self) -> &ContinuousTimeClusterEngine {
        &self.engine
    }

    /// Mutable cluster update engine for diagnostics and validation settings.
    #[inline]
    pub fn engine_mut(&mut self) -> &mut ContinuousTimeClusterEngine {
        &mut self.engine
    }
}

impl MonteCarlo for LongitudinalSpinBosonClusterQmc {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.engine
            .update(&mut self.worldline, &mut ctx.rng)
            .unwrap_or_else(|error| {
                panic!("longitudinal spin-boson cluster sweep failed: {error}")
            });
        ctx.record_attempts(1);
        ctx.record_accepted_moves(1);
        // Per-sweep reports are stored and measured in `measure`; recording
        // them here would mix thermalization into production observables.
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let observables = measure_cluster_observables(
            &self.worldline,
            self.engine.model(),
            self.correlation_bins,
        );
        let report = self.engine.last_report();
        let beta = self.worldline.beta();

        ctx.measure(
            "ClusterMagnetizationSigmaZ",
            observables.magnetization_sigma_z,
        );
        ctx.measure("ClusterMagnetizationSz", observables.magnetization_s_z);
        ctx.measure("ClusterM2SigmaZ", observables.magnetization_sigma_z_squared);
        ctx.measure(
            "ClusterM2Sz",
            observables.magnetization_s_z * observables.magnetization_s_z,
        );
        ctx.measure("ClusterM4SigmaZ", observables.magnetization_sigma_z_fourth);
        ctx.measure(
            "ClusterChiSigmaZRaw",
            beta * observables.magnetization_sigma_z_squared,
        );
        ctx.measure(
            "ClusterChiSzRaw",
            beta * observables.magnetization_s_z * observables.magnetization_s_z,
        );
        ctx.measure_array("ClusterCorrelationSigmaZ", &observables.correlation_sigma_z);
        ctx.measure_array("ClusterCorrelationSz", &observables.correlation_s_z);
        ctx.measure("ClusterKinkCount", observables.kink_count);
        ctx.measure(
            "ClusterTransverseFieldEnergy",
            observables.transverse_field_energy,
        );
        ctx.measure(
            "ClusterLongitudinalFieldEnergy",
            observables.longitudinal_field_energy,
        );

        ctx.measure("ClusterInsertedCuts", report.inserted_cuts as f64);
        ctx.measure("ClusterSegments", report.segments as f64);
        ctx.measure("ClusterSameSpinPairs", report.same_spin_pairs as f64);
        ctx.measure("ClusterRetardedBonds", report.retarded_bonds as f64);
        ctx.measure("ClusterCount", report.clusters as f64);
        ctx.measure("ClusterFinalKinks", report.final_kinks as f64);
        ctx.measure(
            "ClusterImprovedMagnetizationSigmaZ",
            report.improved_magnetization_sigma_z,
        );
        ctx.measure("ClusterImprovedM2SigmaZ", report.improved_m2_sigma_z);
        ctx.measure(
            "ClusterImprovedMagnetizationSz",
            0.5 * report.improved_magnetization_sigma_z,
        );
        ctx.measure("ClusterImprovedM2Sz", 0.25 * report.improved_m2_sigma_z);
        ctx.measure(
            "ClusterImprovedChiSigmaZRaw",
            beta * report.improved_m2_sigma_z,
        );
        ctx.measure(
            "ClusterImprovedChiSzRaw",
            0.25 * beta * report.improved_m2_sigma_z,
        );
    }

    fn name(&self) -> &'static str {
        "LongitudinalSpinBosonClusterQmc"
    }
}

impl FromParams for LongitudinalSpinBosonClusterQmc {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let beta = required::<f64>(params, "beta")?;
        let bath = bath_from_params(params).map_err(to_carlo_error)?;
        let lambda = effective_longitudinal_coupling(params, &bath)?;
        let tunnelling = params
            .get::<f64>("tunnelling")
            .or_else(|| params.get::<f64>("delta"))
            .unwrap_or(0.0);
        let bias = params
            .get::<f64>("epsilon")
            .or_else(|| params.get::<f64>("bias"))
            .unwrap_or(0.0);
        let model = LongitudinalSpinBosonModel::new(
            bath,
            lambda,
            tunnelling,
            bias,
            params.get::<usize>("kernel_quadrature").unwrap_or(64),
        )
        .map_err(to_carlo_error)?;
        let initial_spin = if rng.random::<bool>() { 1 } else { -1 };
        let mut simulation = Self::new(
            model,
            beta,
            initial_spin,
            params.get::<usize>("correlation_bins").unwrap_or(64),
        )
        .map_err(to_carlo_error)?;
        simulation
            .engine
            .set_validate_each_sweep(params.get::<bool>("validate_each_sweep").unwrap_or(false));
        simulation.engine.set_max_auxiliary_cuts(
            params
                .get::<usize>("max_auxiliary_cuts")
                .unwrap_or(1_000_000),
        );
        Ok(simulation)
    }

    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let beta = required::<f64>(params, "beta")?;
        if !beta.is_finite() || beta <= 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "beta".into(),
                reason: format!("must be finite and positive, got {beta}"),
            });
        }
        Ok(())
    }
}

fn bath_from_params(params: &Params) -> Result<Bath, ImpurityError> {
    let bath_name = params
        .get::<String>("bath")
        .unwrap_or_else(|| "powerlaw".into())
        .to_ascii_lowercase();
    match bath_name.as_str() {
        "single" | "single_mode" | "single-mode" => {
            let omega = params.get::<f64>("omega0").unwrap_or(1.0);
            Ok(Bath::SingleMode(SingleModeBath::new(omega)?))
        }
        "powerlaw" | "power_law" | "power-law" => {
            let exponent = params.get::<f64>("s").unwrap_or(0.8);
            let cutoff = params.get::<f64>("omega_c").unwrap_or(1.0);
            Ok(Bath::PowerLaw(PowerLawBath::new(exponent, cutoff)?))
        }
        "tabulated" | "discrete" => {
            let frequencies = parse_csv(params, "bath_omegas")?;
            let weights = parse_csv(params, "bath_weights")?;
            Ok(Bath::Tabulated(TabulatedBath::new(frequencies, weights)?))
        }
        other => Err(ImpurityError::parameter(
            "bath",
            format!("unsupported bath `{other}`"),
        )),
    }
}

fn effective_longitudinal_coupling(params: &Params, bath: &Bath) -> Result<f64, CarloError> {
    if let Some(lambda) = params.get::<f64>("lambda") {
        return Ok(lambda);
    }
    match bath {
        Bath::SingleMode(mode) => {
            if let Some(coupling_s) = params.get::<f64>("g_s") {
                // `g_s S_z (b+b^dagger)` equals `(g_s/2) sigma_z (b+b^dagger)`.
                return Ok(0.25 * coupling_s * coupling_s / mode.omega());
            }
            let coupling_sigma = params
                .get::<f64>("g_sigma")
                .or_else(|| params.get::<f64>("g"))
                .unwrap_or(0.0);
            // `g` and `g_sigma` multiply sigma_z directly in the documented
            // longitudinal Hamiltonian.
            Ok(coupling_sigma * coupling_sigma / mode.omega())
        }
        Bath::PowerLaw(power) => {
            let alpha = params.get::<f64>("alpha").unwrap_or(0.0);
            Ok(2.0 * alpha * power.cutoff() / power.exponent())
        }
        Bath::Tabulated(_) => Err(CarloError::InvalidConfig {
            field: "lambda".into(),
            reason: "a tabulated normalized bath requires explicit `lambda`".into(),
        }),
    }
}

fn parse_csv(params: &Params, key: &str) -> Result<Vec<f64>, ImpurityError> {
    let raw = params
        .get::<String>(key)
        .ok_or_else(|| ImpurityError::parameter(key, "comma-separated values are required"))?;
    raw.split(',')
        .map(|part| {
            part.trim().parse::<f64>().map_err(|_| {
                ImpurityError::parameter(key, format!("cannot parse `{}` as f64", part.trim()))
            })
        })
        .collect()
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
        field: "impurity.cluster".into(),
        reason: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;

    use crate::impurity::spin_boson::bath::{Bath, SingleModeBath};

    use super::*;

    #[test]
    fn constant_worldline_observables_are_exact() {
        let model = LongitudinalSpinBosonModel::with_default_quadrature(
            Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
            0.0,
            0.0,
            0.4,
        )
        .expect("model");
        let worldline = LongitudinalWorldline::new(5.0, -1).expect("worldline");
        let observables = measure_cluster_observables(&worldline, &model, 8);
        assert_eq!(observables.magnetization_sigma_z, -1.0);
        assert!(observables
            .correlation_sigma_z
            .iter()
            .all(|value| (*value - 1.0).abs() < f64::EPSILON));
        assert!((observables.longitudinal_field_energy + 0.2).abs() < 1e-14);
    }

    #[test]
    fn params_construct_single_mode_cluster_solver() {
        let mut params = Params::new();
        params.set("beta", 4.0);
        params.set("bath", "single");
        params.set("omega0", 1.2);
        params.set("g", 0.3);
        params.set("tunnelling", 0.8);
        params.set("epsilon", 0.1);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        let simulation =
            LongitudinalSpinBosonClusterQmc::from_params(&params, &mut rng).expect("simulation");
        assert!((simulation.worldline().beta() - 4.0).abs() < f64::EPSILON);
        assert!((simulation.engine().model().tunnelling() - 0.8).abs() < f64::EPSILON);
        let expected_lambda = 0.3 * 0.3 / 1.2;
        assert!((simulation.engine().model().lambda() - expected_lambda).abs() < 1e-14);
    }
}
