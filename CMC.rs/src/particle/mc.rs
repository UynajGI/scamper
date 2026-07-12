//! Carlo.rs adapter for Lennard-Jones NVT particle simulations.

use crate::algorithms::SimulationPhase;
use crate::particle::{
    CutoffTreatment, LennardJones, OrthorhombicCell, ParticleAlgorithm, ParticleConfiguration,
    ParticleError, ParticleMetropolisCore, ParticleSystem, SimulationCell, TranslateParticle,
};
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};
use std::{fmt::Display, str::FromStr};

/// Generic continuous-particle Monte Carlo chain.
pub struct ParticleMC<const D: usize, P, A>
where
    P: crate::particle::PairPotential,
    A: ParticleAlgorithm<D, P>,
{
    /// Accepted particle state.
    pub system: ParticleSystem<D>,
    /// Pair-potential model.
    pub potential: P,
    /// Transition kernel.
    pub algorithm: A,
}

impl<const D: usize, P, A> ParticleMC<D, P, A>
where
    P: crate::particle::PairPotential,
    A: ParticleAlgorithm<D, P>,
{
    /// Compose an accepted state, potential and particle kernel.
    pub const fn new(system: ParticleSystem<D>, potential: P, algorithm: A) -> Self {
        Self {
            system,
            potential,
            algorithm,
        }
    }
}

impl<const D: usize, P, A> MonteCarlo for ParticleMC<D, P, A>
where
    P: crate::particle::PairPotential,
    A: ParticleAlgorithm<D, P>,
{
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let phase = SimulationPhase::from_run_phase(context.phase());
        self.algorithm
            .sweep_with_phase(&mut self.system, &self.potential, &mut context.rng, phase);
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        let particles = self.system.len() as f64;
        context.measure("Energy", self.system.energy);
        context.measure(
            "EnergyPerParticle",
            if particles == 0.0 {
                0.0
            } else {
                self.system.energy / particles
            },
        );
        context.measure(
            "Density",
            particles / self.system.configuration().cell().volume(),
        );
    }

    fn name(&self) -> &'static str {
        self.algorithm.name()
    }
}

impl<const D: usize, P, A> ParallelTemperingCompatible for ParticleMC<D, P, A>
where
    P: crate::particle::PairPotential,
    A: ParticleAlgorithm<D, P>,
{
    fn log_weight_ratio(&self, parameter: &str, new_value: f64) -> f64 {
        match parameter {
            "beta" => (self.system.beta - new_value) * self.system.energy,
            _ => panic!("unsupported particle PT parameter: {parameter}"),
        }
    }

    fn change_parameter(&mut self, parameter: &str, new_value: f64) {
        match parameter {
            "beta" => self
                .system
                .set_beta(new_value)
                .expect("parallel-tempering beta must be finite and non-negative"),
            _ => panic!("unsupported particle PT parameter: {parameter}"),
        }
    }
}

/// Ready-to-schedule monatomic Lennard-Jones NVT simulation.
pub type LennardJonesNvt<const D: usize> = ParticleMC<D, LennardJones, ParticleMetropolisCore<D>>;

impl<const D: usize> FromParams for LennardJonesNvt<D> {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        parse_lj_params::<D>(params).map(|_| ())
    }

    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let parsed = parse_lj_params::<D>(params)?;
        let cell = OrthorhombicCell::new(parsed.lengths).map_err(particle_config_error)?;
        let positions = regular_grid_positions(parsed.n_particles, &cell);
        let configuration =
            ParticleConfiguration::new(positions, vec![0; parsed.n_particles], cell)
                .map_err(particle_config_error)?;
        let potential = LennardJones::with_treatment(
            parsed.sigma,
            parsed.epsilon,
            parsed.cutoff,
            parsed.treatment,
        )
        .map_err(particle_config_error)?;
        let system = ParticleSystem::new(configuration, &potential, parsed.beta)
            .map_err(particle_config_error)?;
        let translation = TranslateParticle::new(parsed.max_displacement).with_adaptation(
            parsed.target_acceptance,
            parsed.adaptation_interval,
            parsed.adaptation_gain,
            parsed.minimum_displacement,
            parsed.maximum_displacement,
        );
        let algorithm = ParticleMetropolisCore::new(parsed.max_displacement)
            .with_translation(translation)
            .with_energy_check_interval(parsed.energy_check_interval);
        Ok(Self::new(system, potential, algorithm))
    }
}

#[derive(Debug, Clone, Copy)]
struct ParsedLennardJones<const D: usize> {
    n_particles: usize,
    lengths: [f64; D],
    beta: f64,
    sigma: f64,
    epsilon: f64,
    cutoff: f64,
    treatment: CutoffTreatment,
    max_displacement: f64,
    target_acceptance: f64,
    adaptation_interval: u64,
    adaptation_gain: f64,
    minimum_displacement: f64,
    maximum_displacement: f64,
    energy_check_interval: u64,
}

fn parse_lj_params<const D: usize>(params: &Params) -> Result<ParsedLennardJones<D>, CarloError> {
    if D == 0 {
        return Err(invalid("dimension", "must be positive"));
    }
    let default_particles = if D >= 3 { 108usize } else { 64usize };
    let n_particles = required_positive(params, "n_particles", default_particles)?;
    let beta = finite_non_negative(params, "beta", 1.0)?;
    let sigma = finite_positive(params, "sigma", 1.0)?;
    let epsilon = finite_non_negative(params, "epsilon", 1.0)?;
    let density = finite_positive(params, "density", 0.7)?;
    let inferred_length = (n_particles as f64 / density).powf(1.0 / D as f64);
    let mut lengths = [0.0; D];
    for (axis, length) in lengths.iter_mut().enumerate() {
        let key = axis_key(axis);
        *length = if params.contains(key) {
            finite_positive(params, key, inferred_length)?
        } else if params.contains("box_length") {
            finite_positive(params, "box_length", inferred_length)?
        } else {
            inferred_length
        };
    }
    let cutoff = finite_positive(params, "cutoff", 2.5 * sigma)?;
    let treatment = match params
        .get::<String>("cutoff_treatment")
        .unwrap_or_else(|| "shifted_potential".to_string())
        .to_ascii_lowercase()
        .as_str()
    {
        "truncated" => CutoffTreatment::Truncated,
        "shifted" | "shifted_potential" => CutoffTreatment::ShiftedPotential,
        "shifted_force" => CutoffTreatment::ShiftedForce,
        _ => {
            return Err(invalid(
                "cutoff_treatment",
                "expected truncated, shifted_potential, or shifted_force",
            ))
        }
    };
    let max_displacement = finite_positive(params, "max_displacement", 0.1 * sigma)?;
    let target_acceptance = finite_open_unit(params, "target_acceptance", 0.5)?;
    let adaptation_interval = required_positive(params, "adaptation_interval", 20u64)?;
    let adaptation_gain = finite_positive(params, "adaptation_gain", 0.5)?;
    let minimum_displacement =
        finite_positive(params, "minimum_displacement", max_displacement * 1e-6)?;
    let maximum_displacement =
        finite_positive(params, "maximum_displacement", max_displacement * 1e6)?;
    if maximum_displacement < minimum_displacement {
        return Err(invalid(
            "maximum_displacement",
            "must be greater than or equal to minimum_displacement",
        ));
    }
    let energy_check_interval = parse_param(params, "energy_check_interval")?.unwrap_or(0u64);

    let cell = OrthorhombicCell::new(lengths).map_err(particle_config_error)?;
    let potential = LennardJones::with_treatment(sigma, epsilon, cutoff, treatment)
        .map_err(particle_config_error)?;
    crate::particle::CellList::new(
        &ParticleConfiguration::new(Vec::new(), Vec::new(), cell).map_err(particle_config_error)?,
        <LennardJones as crate::particle::PairPotential>::cutoff_squared(&potential),
    )
    .map_err(particle_config_error)?;

    Ok(ParsedLennardJones {
        n_particles,
        lengths,
        beta,
        sigma,
        epsilon,
        cutoff,
        treatment,
        max_displacement,
        target_acceptance,
        adaptation_interval,
        adaptation_gain,
        minimum_displacement,
        maximum_displacement,
        energy_check_interval,
    })
}

fn regular_grid_positions<const D: usize>(
    n_particles: usize,
    cell: &OrthorhombicCell<D>,
) -> Vec<[f64; D]> {
    let points_per_axis = (n_particles as f64).powf(1.0 / D as f64).ceil() as usize;
    let points_per_axis = points_per_axis.max(1);
    let mut positions = Vec::with_capacity(n_particles);
    for linear in 0..n_particles {
        let mut index = linear;
        let mut position = [0.0; D];
        for (axis, value) in position.iter_mut().enumerate() {
            let coordinate = index % points_per_axis;
            index /= points_per_axis;
            *value = (coordinate as f64 + 0.5) * cell.lengths()[axis] / points_per_axis as f64;
        }
        positions.push(position);
    }
    positions
}

fn axis_key(axis: usize) -> &'static str {
    const KEYS: [&str; 6] = ["Lx", "Ly", "Lz", "Lw", "Lv", "Lu"];
    KEYS.get(axis).copied().unwrap_or("box_length")
}

fn parse_param<T>(params: &Params, key: &str) -> Result<Option<T>, CarloError>
where
    T: FromStr,
    T::Err: Display,
{
    if !params.contains(key) {
        return Ok(None);
    }
    params
        .get::<T>(key)
        .map(Some)
        .ok_or_else(|| invalid(key, "could not parse value"))
}

fn required_positive<T>(params: &Params, key: &str, default: T) -> Result<T, CarloError>
where
    T: FromStr + Copy + PartialOrd + From<u8>,
    T::Err: Display,
{
    let value = parse_param(params, key)?.unwrap_or(default);
    if value > T::from(0) {
        Ok(value)
    } else {
        Err(invalid(key, "must be positive"))
    }
}

fn finite_positive(params: &Params, key: &str, default: f64) -> Result<f64, CarloError> {
    let value = parse_param(params, key)?.unwrap_or(default);
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(invalid(key, "must be finite and positive"))
    }
}

fn finite_non_negative(params: &Params, key: &str, default: f64) -> Result<f64, CarloError> {
    let value = parse_param(params, key)?.unwrap_or(default);
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(invalid(key, "must be finite and non-negative"))
    }
}

fn finite_open_unit(params: &Params, key: &str, default: f64) -> Result<f64, CarloError> {
    let value = parse_param(params, key)?.unwrap_or(default);
    if value.is_finite() && (0.0..1.0).contains(&value) {
        Ok(value)
    } else {
        Err(invalid(key, "must lie strictly between zero and one"))
    }
}

fn invalid(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}

fn particle_config_error(error: ParticleError) -> CarloError {
    invalid("particle_system", error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use carlo_rs::{RayonBackend, RunConfig, Scheduler};

    #[test]
    fn scheduler_runs_two_dimensional_lj_nvt() {
        let mut params = Params::new();
        params.set("n_particles", 25usize);
        params.set("density", 0.5);
        params.set("cutoff", 2.0);
        params.set("beta", 0.8);
        params.set("max_displacement", 0.15);
        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 20,
                measurement_sweeps: 40,
                binsize: 10,
                base_seed: 71,
                ..Default::default()
            },
        )
        .run_one::<LennardJonesNvt<2>>(&params);
        assert!(results.get("Energy").is_some());
        assert!(results.get("EnergyPerParticle").is_some());
        assert!(results.get("Density").is_some());
    }

    #[test]
    fn invalid_minimum_image_cutoff_is_rejected() {
        let mut params = Params::new();
        params.set("n_particles", 8usize);
        params.set("box_length", 4.0);
        params.set("cutoff", 2.1);
        assert!(LennardJonesNvt::<3>::validate_params(&params).is_err());
    }
}
