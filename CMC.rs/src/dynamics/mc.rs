//! Carlo.rs adapters for Stage 6 dynamics kernels.

use super::{
    BklIsingKernel, DynamicsError, HardSphereEventChain, KawasakiCore, KineticIsingModel,
    KineticRateLaw,
};
use crate::algorithms::{Algorithm, SimulationPhase};
use crate::classical_mc::{build_lattice_from_params, parse_param};
use crate::lattice::models::IsingModel;
use crate::lattice::state::System;
use crate::particle::{OrthorhombicCell, ParticleConfiguration, SimulationCell};
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Scheduler-ready conserved-order-parameter Ising dynamics.
pub struct KawasakiIsingMC {
    system: System,
    model: IsingModel,
    kernel: KawasakiCore,
}

impl KawasakiIsingMC {
    pub fn new(
        system: System,
        model: IsingModel,
        kernel: KawasakiCore,
    ) -> Result<Self, DynamicsError> {
        if model.spin_dim() != 1 {
            return Err(DynamicsError::new(
                "Kawasaki adapter requires scalar Ising spins",
            ));
        }
        let mut result = Self {
            system,
            model,
            kernel,
        };
        result.system.recompute_energy(&result.model);
        result
            .system
            .validate(&result.model)
            .map_err(DynamicsError::new)?;
        Ok(result)
    }

    #[inline]
    pub const fn system(&self) -> &System {
        &self.system
    }

    #[inline]
    pub const fn kernel(&self) -> &KawasakiCore {
        &self.kernel
    }
}

impl MonteCarlo for KawasakiIsingMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let phase = SimulationPhase::from_run_phase(context.phase());
        self.kernel
            .sweep_with_phase(&mut self.system, &self.model, &mut context.rng, phase);
        context.record_attempts(self.kernel.last_attempts());
        context.record_accepted_moves(self.kernel.last_accepts());
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        context.measure("Energy", self.system.energy);
        context.measure("Magnetization", signed_magnetization(&self.system));
        context.measure(
            "KawasakiAcceptance",
            ratio(self.kernel.last_accepts(), self.kernel.last_attempts()),
        );
        context.measure("AttemptClock", context.attempted_updates() as f64);
        context.measure("AcceptedMoveClock", context.accepted_moves() as f64);
    }

    fn name(&self) -> &'static str {
        "KawasakiIsing"
    }
}

impl FromParams for KawasakiIsingMC {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
        let lattice = build_lattice_from_params(params, parse_bool(params, "pbc", true)?)?;
        let mut system = System::new(lattice, 1, 1.0, beta);
        initialize_ising_spins(&mut system, params, rng)?;
        let attempts = parse_param::<usize>(params, "kawasaki_attempts_per_sweep")?.unwrap_or(0);
        let audit = parse_param::<u64>(params, "cache_audit_interval")?.unwrap_or(0);
        Self::new(
            system,
            IsingModel::new(coupling),
            KawasakiCore::new(attempts).with_cache_audit_interval(audit),
        )
        .map_err(|error| invalid("kawasaki", error.to_string()))
    }
}

/// Scheduler-ready rejection-free BKL Ising dynamics sampled at fixed event-time intervals.
pub struct KineticIsingBklMC {
    kernel: BklIsingKernel,
    event_time_per_sweep: f64,
    last_events: u64,
    last_delta_time: f64,
}

impl KineticIsingBklMC {
    pub fn new(kernel: BklIsingKernel, event_time_per_sweep: f64) -> Result<Self, DynamicsError> {
        if !event_time_per_sweep.is_finite() || event_time_per_sweep <= 0.0 {
            return Err(DynamicsError::new(
                "event_time_per_sweep must be finite and positive",
            ));
        }
        Ok(Self {
            kernel,
            event_time_per_sweep,
            last_events: 0,
            last_delta_time: 0.0,
        })
    }

    #[inline]
    pub const fn kernel(&self) -> &BklIsingKernel {
        &self.kernel
    }
}

impl MonteCarlo for KineticIsingBklMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let before_time = self.kernel.event_time();
        self.last_events = self
            .kernel
            .advance_by(self.event_time_per_sweep, &mut context.rng)
            .expect("BKL Ising event-time advance failed");
        self.last_delta_time = self.kernel.event_time() - before_time;
        context.record_attempts(self.last_events);
        context.record_accepted_moves(self.last_events);
        context.advance_event_time(self.last_delta_time);
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        context.measure("Energy", self.kernel.state().energy);
        context.measure("Magnetization", signed_magnetization(self.kernel.state()));
        context.measure("EventTime", context.event_time());
        context.measure("EventsPerWindow", self.last_events as f64);
        context.measure("TotalRate", self.kernel.total_rate());
    }

    fn name(&self) -> &'static str {
        "BklIsingKineticMC"
    }
}

impl FromParams for KineticIsingBklMC {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
        let lattice = build_lattice_from_params(params, parse_bool(params, "pbc", true)?)?;
        let mut system = System::new(lattice, 1, 1.0, beta);
        initialize_ising_spins(&mut system, params, rng)?;
        let frequency = parse_param::<f64>(params, "attempt_frequency")?.unwrap_or(1.0);
        let rate_name = params
            .get::<String>("kinetic_rate")
            .unwrap_or_else(|| "glauber".to_string());
        let rate_law = match rate_name.as_str() {
            "glauber" => KineticRateLaw::glauber(frequency),
            "metropolis" => KineticRateLaw::metropolis(frequency),
            other => Err(DynamicsError::new(format!(
                "unknown kinetic_rate `{other}`"
            ))),
        }
        .map_err(|error| invalid("kinetic_rate", error.to_string()))?;
        let model = KineticIsingModel::new(coupling, rate_law)
            .map_err(|error| invalid("kinetic_model", error.to_string()))?;
        let audit = parse_param::<u64>(params, "cache_audit_interval")?.unwrap_or(0);
        let kernel = BklIsingKernel::new(model, system, audit)
            .map_err(|error| invalid("kinetic_model", error.to_string()))?;
        let window = parse_param::<f64>(params, "event_time_per_sweep")?.unwrap_or(1.0);
        Self::new(kernel, window)
            .map_err(|error| invalid("event_time_per_sweep", error.to_string()))
    }
}

/// Scheduler-ready identical-hard-sphere event-chain simulation.
pub struct HardSphereEventChainMC<const D: usize> {
    kernel: HardSphereEventChain<D>,
    chains_per_sweep: usize,
    last_collisions: u64,
}

impl<const D: usize> HardSphereEventChainMC<D> {
    pub fn new(
        kernel: HardSphereEventChain<D>,
        chains_per_sweep: usize,
    ) -> Result<Self, DynamicsError> {
        if chains_per_sweep == 0 {
            return Err(DynamicsError::new("chains_per_sweep must be positive"));
        }
        Ok(Self {
            kernel,
            chains_per_sweep,
            last_collisions: 0,
        })
    }

    #[inline]
    pub const fn kernel(&self) -> &HardSphereEventChain<D> {
        &self.kernel
    }
}

impl<const D: usize> MonteCarlo for HardSphereEventChainMC<D> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let before = self.kernel.collisions();
        for _ in 0..self.chains_per_sweep {
            self.kernel
                .step(&mut context.rng)
                .expect("hard-sphere event chain failed");
        }
        self.last_collisions = self.kernel.collisions().saturating_sub(before);
        context.record_attempts(self.chains_per_sweep as u64);
        context.record_accepted_moves(self.chains_per_sweep as u64);
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        let configuration = self.kernel.configuration();
        let particle_volume = unit_ball_volume(D) * (0.5 * self.kernel.diameter()).powi(D as i32);
        let packing_fraction =
            configuration.len() as f64 * particle_volume / configuration.cell().volume();
        context.measure("PackingFraction", packing_fraction);
        context.measure("EventChainCollisions", self.last_collisions as f64);
        context.measure("LiftedDistance", self.kernel.lifted_distance());
    }

    fn name(&self) -> &'static str {
        "HardSphereEventChain"
    }
}

impl<const D: usize> FromParams for HardSphereEventChainMC<D> {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let particles = parse_param::<usize>(params, "n_particles")?.unwrap_or(16);
        let box_length = parse_param::<f64>(params, "box_length")?.unwrap_or(8.0);
        let diameter = parse_param::<f64>(params, "diameter")?.unwrap_or(1.0);
        let chain_length = parse_param::<f64>(params, "chain_length")?.unwrap_or(box_length);
        let chains_per_sweep =
            parse_param::<usize>(params, "chains_per_sweep")?.unwrap_or(particles.max(1));
        let audit = parse_param::<u64>(params, "cache_audit_interval")?.unwrap_or(0);
        let configuration = grid_configuration::<D>(particles, box_length, diameter)
            .map_err(|error| invalid("hard_sphere_configuration", error.to_string()))?;
        let kernel = HardSphereEventChain::new(configuration, diameter, chain_length, audit)
            .map_err(|error| invalid("event_chain", error.to_string()))?;
        Self::new(kernel, chains_per_sweep)
            .map_err(|error| invalid("event_chain", error.to_string()))
    }
}

fn initialize_ising_spins(
    system: &mut System,
    params: &Params,
    rng: &mut Xoshiro256PlusPlus,
) -> Result<(), CarloError> {
    let up_fraction = parse_param::<f64>(params, "up_fraction")?.unwrap_or(0.5);
    if !up_fraction.is_finite() || !(0.0..=1.0).contains(&up_fraction) {
        return Err(invalid("up_fraction", "must lie in [0, 1]"));
    }
    for spin in &mut system.spins {
        *spin = if rng.random::<f64>() < up_fraction {
            1.0
        } else {
            -1.0
        };
    }
    Ok(())
}

fn grid_configuration<const D: usize>(
    particles: usize,
    box_length: f64,
    diameter: f64,
) -> Result<ParticleConfiguration<D>, DynamicsError> {
    if D == 0 || particles == 0 {
        return Err(DynamicsError::new(
            "hard-sphere grid requires positive dimension and particle count",
        ));
    }
    if !box_length.is_finite() || box_length <= 0.0 {
        return Err(DynamicsError::new("box_length must be finite and positive"));
    }
    let mut cells_per_axis = 1usize;
    while cells_per_axis.saturating_pow(D as u32) < particles {
        cells_per_axis = cells_per_axis.saturating_add(1);
    }
    let spacing = box_length / cells_per_axis as f64;
    if diameter >= spacing {
        return Err(DynamicsError::new(
            "grid spacing must exceed the hard-sphere diameter",
        ));
    }
    let mut positions = Vec::with_capacity(particles);
    for index in 0..particles {
        let mut value = index;
        let mut position = [0.0; D];
        for coordinate in &mut position {
            let cell_index = value % cells_per_axis;
            value /= cells_per_axis;
            *coordinate = (cell_index as f64 + 0.5) * spacing;
        }
        positions.push(position);
    }
    let cell = OrthorhombicCell::new([box_length; D])
        .map_err(|error| DynamicsError::new(error.to_string()))?;
    ParticleConfiguration::new(positions, vec![0; particles], cell)
        .map_err(|error| DynamicsError::new(error.to_string()))
}

fn unit_ball_volume(dimension: usize) -> f64 {
    match dimension {
        0 => 1.0,
        1 => 2.0,
        _ => 2.0 * std::f64::consts::PI / dimension as f64 * unit_ball_volume(dimension - 2),
    }
}

fn signed_magnetization(system: &System) -> f64 {
    system.spins.iter().sum::<f64>() / system.n_sites().max(1) as f64
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn parse_bool(params: &Params, name: &str, default: bool) -> Result<bool, CarloError> {
    if !params.contains(name) {
        return Ok(default);
    }
    params
        .get::<bool>(name)
        .ok_or_else(|| invalid(name, "cannot parse value as bool"))
}

fn invalid(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}
