//! Carlo.rs adapter for the generic impurity wormhole engine.

use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;

use crate::algorithm::{QmcKernel, UpdateSchedule};

use super::bath::{Bath, PowerLawBath, SingleModeBath, TabulatedBath};
use super::configuration::WormholeConfiguration;
use super::error::ImpurityError;
use super::model::{CouplingNormalization, ImpurityModel};
use super::observables::measure_observables;
use super::updates::{LoopStartPolicy, WormholeEngine};

/// Runnable quantum-impurity QMC simulation.
///
/// Carlo.rs calls [`MonteCarlo::sweep`] and [`MonteCarlo::measure`]; this type
/// delegates physics updates to [`WormholeEngine`] and records observables in
/// the Carlo.rs measurement context.
pub struct ImpurityQmc {
    configuration: WormholeConfiguration,
    engine: WormholeEngine,
    correlation_samples: usize,
    adaptive_schedule: bool,
    adaptation_interval: usize,
    adaptation_samples: usize,
    adaptation_order_sum: usize,
}

impl ImpurityQmc {
    /// Construct a runnable simulation from an already-built model.
    pub fn new(
        model: ImpurityModel,
        beta: f64,
        empty_spin: i8,
        schedule: UpdateSchedule,
        correlation_samples: usize,
    ) -> Result<Self, ImpurityError> {
        Ok(Self {
            configuration: WormholeConfiguration::new(beta, empty_spin)?,
            engine: WormholeEngine::new(model, schedule),
            correlation_samples: correlation_samples.max(1),
            adaptive_schedule: false,
            adaptation_interval: 100,
            adaptation_samples: 0,
            adaptation_order_sum: 0,
        })
    }

    /// Current sampled configuration.
    pub fn configuration(&self) -> &WormholeConfiguration {
        &self.configuration
    }

    /// Wormhole engine and model catalog.
    pub fn engine(&self) -> &WormholeEngine {
        &self.engine
    }

    /// Enable warmup-only work-count adaptation.
    pub fn set_adaptive_schedule(&mut self, enabled: bool, interval: usize) {
        self.adaptive_schedule = enabled;
        self.adaptation_interval = interval.max(1);
    }

    fn adapt_during_warmup(&mut self, thermalized: bool) {
        if !self.adaptive_schedule || thermalized {
            return;
        }
        self.adaptation_samples += 1;
        self.adaptation_order_sum += self.configuration.expansion_order();
        if self.adaptation_samples < self.adaptation_interval {
            return;
        }
        let mean_order = self.adaptation_order_sum as f64 / self.adaptation_samples as f64;
        let diagonal_proposals = mean_order.round().max(1.0) as usize;
        let directed_loops = mean_order.sqrt().ceil().max(1.0) as usize;
        let previous = self.engine.schedule();
        self.engine.set_schedule(UpdateSchedule::new(
            diagonal_proposals,
            directed_loops,
            previous.max_loop_steps_factor,
        ));
        self.adaptation_samples = 0;
        self.adaptation_order_sum = 0;
    }
}

impl MonteCarlo for ImpurityQmc {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.engine
            .sweep(&mut self.configuration, &mut ctx.rng)
            .unwrap_or_else(|error| panic!("impurity wormhole sweep failed: {error}"));
        self.adapt_during_warmup(ctx.is_thermalized());
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let observables = measure_observables(
            &self.configuration,
            self.engine.model(),
            self.correlation_samples,
            &mut ctx.rng,
        )
        .unwrap_or_else(|error| panic!("impurity measurement failed: {error}"));

        ctx.measure("MagnetizationSigmaZ", observables.magnetization_sigma_z);
        ctx.measure("MagnetizationSz", observables.magnetization_s_z);
        ctx.measure("M2SigmaZ", observables.magnetization_sigma_z_squared);
        ctx.measure("M4SigmaZ", observables.magnetization_sigma_z_fourth);
        ctx.measure("ChiZ", observables.susceptibility_z);
        ctx.measure(
            "CorrelationSigmaZHalf",
            observables.correlation_sigma_z_half,
        );
        ctx.measure("CorrelationSzHalf", observables.correlation_s_z_half);
        ctx.measure("ExpansionOrder", observables.expansion_order);
        ctx.measure("DiagonalOrder", observables.diagonal_order);
        ctx.measure("OffDiagonalOrder", observables.offdiagonal_order);
        ctx.measure(
            "ShiftedInteractionEnergy",
            observables.shifted_interaction_energy,
        );

        let stats = self.engine.stats();
        ctx.measure("DiagonalAcceptance", stats.diagonal_acceptance());
        ctx.measure("MeanLoopSteps", stats.mean_loop_steps());
        ctx.measure("BounceFraction", stats.bounce_fraction());
        ctx.measure("WormholeFraction", stats.wormhole_fraction());
        ctx.measure("LoopAbortFraction", stats.loop_abort_fraction());
    }

    fn name(&self) -> &'static str {
        "ImpurityQmc"
    }
}

impl FromParams for ImpurityQmc {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let beta = required::<f64>(params, "beta")?;
        let bath = bath_from_params(params).map_err(to_carlo_error)?;
        let model_name = params
            .get::<String>("model")
            .unwrap_or_else(|| "xxz".into())
            .to_ascii_lowercase();
        let constant = params.get::<f64>("C");
        let h_z = params.get::<f64>("h_z").unwrap_or(0.0);

        let model = match model_name.as_str() {
            "jc" | "jaynes_cummings" | "jaynes-cummings" => {
                let lambda = effective_coupling(params, &bath, "lambda", "g", "alpha")?;
                ImpurityModel::jaynes_cummings(bath, lambda, h_z, constant)
            }
            "rw_crw" | "rw-crw" | "weber" => {
                let vertex_scale = effective_rw_crw_scale(params, &bath)?;
                let crw_ratio = params
                    .get::<f64>("crw_ratio")
                    .or_else(|| params.get::<f64>("delta"))
                    .unwrap_or(0.0);
                let tunnelling = params
                    .get::<f64>("tunnelling")
                    .or_else(|| params.get::<f64>("h_x"))
                    .unwrap_or(0.0);
                let normalization = coupling_normalization_from_params(params)?;
                let rw_crw_constant = rw_crw_constant_from_params(params, tunnelling)?;
                ImpurityModel::rw_crw(
                    bath,
                    vertex_scale,
                    crw_ratio,
                    tunnelling,
                    normalization,
                    rw_crw_constant,
                )
            }
            "xxz" => {
                let lambda_xy = effective_coupling(params, &bath, "lambda_xy", "g_xy", "alpha_xy")?;
                let lambda_z = effective_coupling(params, &bath, "lambda_z", "g_z", "alpha_z")?;
                ImpurityModel::xxz(bath, lambda_xy, lambda_z, h_z, constant)
            }
            "xyz" => {
                let lambda_x = effective_coupling(params, &bath, "lambda_x", "g_x", "alpha_x")?;
                let lambda_y = effective_coupling(params, &bath, "lambda_y", "g_y", "alpha_y")?;
                let lambda_z = effective_coupling(params, &bath, "lambda_z", "g_z", "alpha_z")?;
                ImpurityModel::xyz(bath, lambda_x, lambda_y, lambda_z, h_z, constant)
            }
            "impurity" | "rabi" | "rotated_impurity" => {
                let lambda = effective_rabi_coupling(params, &bath)?;
                let tunnelling = params
                    .get::<f64>("tunnelling")
                    .or_else(|| params.get::<f64>("h_x"))
                    .unwrap_or(h_z);
                ImpurityModel::rotated_impurity(bath, lambda, tunnelling, constant)
            }
            other => Err(ImpurityError::parameter(
                "model",
                format!("unsupported model `{other}`"),
            )),
        }
        .map_err(to_carlo_error)?;

        let schedule = UpdateSchedule::new(
            params.get::<usize>("diagonal_proposals").unwrap_or(1),
            params.get::<usize>("directed_loops").unwrap_or(1),
            params.get::<usize>("max_loop_steps_factor").unwrap_or(16),
        );
        let empty_spin = if rng.random::<bool>() { 1 } else { -1 };
        let mut simulation = Self::new(
            model,
            beta,
            empty_spin,
            schedule,
            params.get::<usize>("correlation_samples").unwrap_or(16),
        )
        .map_err(to_carlo_error)?;
        simulation
            .engine
            .set_validate_each_sweep(params.get::<bool>("validate_each_sweep").unwrap_or(false));
        simulation
            .engine
            .set_loop_start_policy(loop_start_policy_from_params(params)?);
        simulation.set_adaptive_schedule(
            params.get::<bool>("adaptive_schedule").unwrap_or(true),
            params.get::<usize>("adaptation_interval").unwrap_or(100),
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

fn coupling_normalization_from_params(
    params: &Params,
) -> Result<CouplingNormalization, CarloError> {
    let name = params
        .get::<String>("coupling_normalization")
        .unwrap_or_else(|| "fixed-rw".into())
        .to_ascii_lowercase();
    match name.as_str() {
        "fixed-rw" | "fixed_rw" | "rw" => Ok(CouplingNormalization::FixedRw),
        "fixed-total" | "fixed_total" | "total" => Ok(CouplingNormalization::FixedTotal),
        "fixed-quadratic" | "fixed_quadratic" | "quadratic" => {
            Ok(CouplingNormalization::FixedQuadratic)
        }
        other => Err(CarloError::InvalidConfig {
            field: "coupling_normalization".into(),
            reason: format!("unsupported normalization `{other}`"),
        }),
    }
}

fn rw_crw_constant_from_params(
    params: &Params,
    tunnelling: f64,
) -> Result<Option<f64>, CarloError> {
    let explicit_constant = params.get::<f64>("C");
    let diagonal_shift = params.get::<f64>("diagonal_shift");
    if explicit_constant.is_some() && diagonal_shift.is_some() {
        return Err(CarloError::InvalidConfig {
            field: "C".into(),
            reason: "use either `C` or `diagonal_shift`, not both".into(),
        });
    }
    if let Some(shift) = diagonal_shift {
        if !shift.is_finite() || shift < 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "diagonal_shift".into(),
                reason: format!("must be finite and non-negative, got {shift}"),
            });
        }
        return Ok(Some(0.5 * tunnelling.abs() + shift + 16.0 * f64::EPSILON));
    }
    Ok(explicit_constant)
}

fn loop_start_policy_from_params(params: &Params) -> Result<LoopStartPolicy, CarloError> {
    let name = params
        .get::<String>("loop_start")
        .unwrap_or_else(|| "random-time".into())
        .to_ascii_lowercase();
    match name.as_str() {
        "random-time" | "random_time" | "time" => Ok(LoopStartPolicy::RandomTime),
        "random-leg" | "random_leg" | "leg" => Ok(LoopStartPolicy::RandomLeg),
        other => Err(CarloError::InvalidConfig {
            field: "loop_start".into(),
            reason: format!("unsupported loop-start policy `{other}`"),
        }),
    }
}

fn effective_rw_crw_scale(params: &Params, bath: &Bath) -> Result<f64, CarloError> {
    if let Some(vertex_scale) = params.get::<f64>("vertex_scale") {
        return Ok(vertex_scale);
    }
    if let Some(lambda) = params.get::<f64>("lambda") {
        return Ok(lambda);
    }
    effective_coupling(params, bath, "vertex_scale", "g", "alpha")
}

fn effective_rabi_coupling(params: &Params, bath: &Bath) -> Result<f64, CarloError> {
    if let Some(lambda) = params.get::<f64>("lambda") {
        return Ok(lambda);
    }
    match bath {
        Bath::SingleMode(mode) => {
            // Two explicit conventions are supported:
            //   g       multiplies S_z (or S_x after rotation),
            //   g_sigma multiplies sigma_z = 2 S_z.
            // Therefore g_sigma produces four times the retarded weight.
            if let Some(g_sigma) = params.get::<f64>("g_sigma") {
                return Ok(4.0 * g_sigma * g_sigma / mode.omega());
            }
            let coupling_s = params.get::<f64>("g").unwrap_or(0.0);
            Ok(coupling_s * coupling_s / mode.omega())
        }
        Bath::PowerLaw(power) => {
            let alpha = params.get::<f64>("alpha").unwrap_or(0.0);
            Ok(2.0 * alpha * power.cutoff() / power.exponent())
        }
        Bath::Tabulated(_) => Err(CarloError::InvalidConfig {
            field: "lambda".into(),
            reason: "a tabulated normalized Rabi bath requires explicit `lambda`".into(),
        }),
    }
}

fn effective_coupling(
    params: &Params,
    bath: &Bath,
    lambda_key: &str,
    single_key: &str,
    alpha_key: &str,
) -> Result<f64, CarloError> {
    if let Some(lambda) = params.get::<f64>(lambda_key) {
        return Ok(lambda);
    }
    match bath {
        Bath::SingleMode(mode) => {
            let coupling = params.get::<f64>(single_key).unwrap_or(0.0);
            Ok(coupling * coupling / mode.omega())
        }
        Bath::PowerLaw(power) => {
            let alpha = params.get::<f64>(alpha_key).unwrap_or(0.0);
            Ok(2.0 * alpha * power.cutoff() / power.exponent())
        }
        Bath::Tabulated(_) => {
            params
                .get::<f64>(lambda_key)
                .ok_or_else(|| CarloError::InvalidConfig {
                    field: lambda_key.into(),
                    reason: "a tabulated normalized bath requires an explicit effective coupling"
                        .into(),
                })
        }
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
        field: "impurity".into(),
        reason: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use carlo_rs::{RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

    use super::*;

    #[test]
    fn from_params_builds_each_model() {
        for model in ["jc", "rw_crw", "xxz", "xyz", "rabi"] {
            let mut params = Params::new();
            params.set("beta", 4.0);
            params.set("model", model);
            params.set("bath", "single");
            params.set("omega0", 1.0);
            params.set("g", 0.3);
            params.set("g_xy", 0.3);
            params.set("g_x", 0.3);
            params.set("crw_ratio", 0.2);
            params.set("tunnelling", 0.1);
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
            let simulation = ImpurityQmc::from_params(&params, &mut rng);
            assert!(simulation.is_ok(), "failed to build {model}");
        }
    }

    #[test]
    fn rw_crw_diagonal_shift_uses_the_weber_positivity_bound() {
        let mut params = Params::new();
        params.set("diagonal_shift", 0.2);
        let constant = rw_crw_constant_from_params(&params, -0.6)
            .expect("constant")
            .expect("explicit shifted constant");
        assert!((constant - 0.5).abs() < 1.0e-14);
    }

    #[test]
    fn standard_rabi_sigma_convention_has_factor_four() {
        let mut params = Params::new();
        params.set("g_sigma", 0.3);
        let bath = Bath::SingleMode(SingleModeBath::new(2.0).expect("mode"));
        let lambda = effective_rabi_coupling(&params, &bath).expect("coupling");
        assert!((lambda - 0.18).abs() < 1.0e-14);
    }

    #[test]
    fn scheduler_runs_end_to_end() {
        let mut params = Params::new();
        params.set("beta", 2.0);
        params.set("model", "jc");
        params.set("bath", "single");
        params.set("omega0", 1.0);
        params.set("g", 0.25);
        params.set("h_z", 0.2);
        params.set("validate_each_sweep", true);
        let config = RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 400,
            binsize: 20,
            base_seed: 123,
            ..Default::default()
        };
        let results = Scheduler::new(RayonBackend::new(1), config).run_one::<ImpurityQmc>(&params);
        assert!(results.get("ExpansionOrder").is_some());
        assert!(results.get("MagnetizationSz").is_some());
    }
}
