//! Isotropic volume proposals and NPT Metropolis-Hastings updates.

use crate::algorithms::SimulationPhase;
use crate::core::acceptance::MetropolisHastingsAcceptance;
use crate::core::ensemble::{IsothermalIsobaric, ThermodynamicDelta};
use crate::core::trial::{metropolis_hastings_step, ProposedMove, TrialEvaluator};
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::particle::{
    compute_total_energy, CellList, OrthorhombicCell, PairPotential, ParticleAlgorithm,
    ParticleEnergyPatch, ParticleSystem, SimulationCell, TranslateParticle,
};
use rand::{Rng, RngExt};

/// Isotropic change expressed as `ln(V_new / V_old)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsotropicVolumeChange {
    pub log_volume_ratio: f64,
}

impl IsotropicVolumeChange {
    pub fn new(log_volume_ratio: f64) -> Self {
        assert!(
            log_volume_ratio.is_finite(),
            "log-volume change must be finite"
        );
        Self { log_volume_ratio }
    }
}

/// Reusable full-state patch for an accepted volume change.
#[derive(Debug, Clone)]
pub struct VolumeChangePatch<const D: usize> {
    /// Physical potential-energy change.
    pub delta_energy: f64,
    /// Physical volume change.
    pub delta_volume: f64,
    new_energy: f64,
    new_cell: Option<OrthorhombicCell<D>>,
    new_positions: Vec<[f64; D]>,
    valid: bool,
}

impl<const D: usize> Default for VolumeChangePatch<D> {
    fn default() -> Self {
        Self {
            delta_energy: 0.0,
            delta_volume: 0.0,
            new_energy: 0.0,
            new_cell: None,
            new_positions: Vec::new(),
            valid: false,
        }
    }
}

impl<const D: usize> ParticleSystem<D> {
    pub(crate) fn evaluate_volume_change<P: PairPotential>(
        &self,
        potential: &P,
        movement: &IsotropicVolumeChange,
        patch: &mut VolumeChangePatch<D>,
    ) -> ThermodynamicDelta {
        patch.valid = false;
        patch.new_cell = None;
        patch.new_positions.clear();

        let scale = (movement.log_volume_ratio / D as f64).exp();
        if !scale.is_finite() || scale <= 0.0 {
            patch.delta_energy = f64::INFINITY;
            patch.delta_volume = 0.0;
            return ThermodynamicDelta::volume(f64::INFINITY, 0.0, 0.0);
        }
        let new_cell = match self.configuration.cell().scaled(scale) {
            Ok(cell) => cell,
            Err(_) => {
                patch.delta_energy = f64::INFINITY;
                patch.delta_volume = 0.0;
                return ThermodynamicDelta::volume(f64::INFINITY, 0.0, 0.0);
            }
        };
        let cutoff = potential.cutoff_squared().sqrt();
        if cutoff > 0.5 * new_cell.minimum_length() * (1.0 + 16.0 * f64::EPSILON) {
            patch.delta_energy = f64::INFINITY;
            patch.delta_volume = 0.0;
            return ThermodynamicDelta::volume(f64::INFINITY, 0.0, 0.0);
        }

        patch.new_positions.reserve(self.len());
        for position in self.configuration.positions() {
            let mut scaled = *position;
            for coordinate in &mut scaled {
                *coordinate *= scale;
            }
            patch.new_positions.push(scaled);
        }
        let trial_configuration = crate::particle::ParticleConfiguration::new(
            patch.new_positions.clone(),
            self.configuration.species().to_vec(),
            new_cell,
        )
        .expect("finite isotropic scaling must produce a valid configuration");
        let new_energy = compute_total_energy(&trial_configuration, potential);
        if !new_energy.is_finite() {
            patch.delta_energy = f64::INFINITY;
            patch.delta_volume = 0.0;
            return ThermodynamicDelta::volume(f64::INFINITY, 0.0, 0.0);
        }

        let old_volume = self.configuration.cell().volume();
        let new_volume = new_cell.volume();
        patch.delta_energy = new_energy - self.energy;
        patch.delta_volume = new_volume - old_volume;
        patch.new_energy = new_energy;
        patch.new_cell = Some(new_cell);
        patch.valid = true;
        ThermodynamicDelta::volume(
            patch.delta_energy,
            patch.delta_volume,
            self.len() as f64 * movement.log_volume_ratio,
        )
    }

    pub(crate) fn commit_volume_change(
        &mut self,
        _movement: &IsotropicVolumeChange,
        patch: &VolumeChangePatch<D>,
    ) {
        assert!(patch.valid, "invalid volume trial must never be committed");
        let new_cell = patch.new_cell.expect("valid volume patch has a cell");
        self.configuration
            .set_positions_and_cell(patch.new_positions.clone(), new_cell);
        self.cell_list = CellList::new(&self.configuration, self.cell_list.cutoff_squared())
            .expect("validated volume trial must rebuild its cell list");
        self.energy = patch.new_energy;
    }
}

impl<const D: usize, P: PairPotential> TrialEvaluator<P, IsotropicVolumeChange>
    for ParticleSystem<D>
{
    type Delta = ThermodynamicDelta;
    type Patch = VolumeChangePatch<D>;

    fn evaluate_trial(
        &self,
        model: &P,
        movement: &IsotropicVolumeChange,
        patch: &mut Self::Patch,
    ) -> Self::Delta {
        self.evaluate_volume_change(model, movement, patch)
    }

    fn commit_trial(&mut self, movement: &IsotropicVolumeChange, patch: &Self::Patch) {
        self.commit_volume_change(movement, patch);
    }
}

/// Symmetric random walk in logarithmic volume with warmup adaptation.
#[derive(Debug, Clone)]
pub struct LogVolumeScale {
    max_log_volume_change: f64,
    minimum: f64,
    maximum: f64,
    target_acceptance: f64,
    gain: f64,
    interval_sweeps: u64,
    window_sweeps: u64,
    window_attempts: u64,
    window_accepted: u64,
    total_attempts: u64,
    total_accepted: u64,
}

impl LogVolumeScale {
    pub fn new(max_log_volume_change: f64) -> Self {
        assert!(
            max_log_volume_change.is_finite() && max_log_volume_change > 0.0,
            "maximum log-volume change must be finite and positive"
        );
        Self {
            max_log_volume_change,
            minimum: max_log_volume_change * 1e-6,
            maximum: (max_log_volume_change * 1e6).min(100.0),
            target_acceptance: 0.3,
            gain: 0.5,
            interval_sweeps: 20,
            window_sweeps: 0,
            window_attempts: 0,
            window_accepted: 0,
            total_attempts: 0,
            total_accepted: 0,
        }
    }

    pub fn with_adaptation(
        mut self,
        target_acceptance: f64,
        interval_sweeps: u64,
        gain: f64,
        minimum: f64,
        maximum: f64,
    ) -> Self {
        assert!(target_acceptance.is_finite() && (0.0..1.0).contains(&target_acceptance));
        assert!(interval_sweeps > 0);
        assert!(gain.is_finite() && gain > 0.0);
        assert!(minimum.is_finite() && maximum.is_finite() && minimum > 0.0 && maximum >= minimum);
        self.target_acceptance = target_acceptance;
        self.interval_sweeps = interval_sweeps;
        self.gain = gain;
        self.minimum = minimum;
        self.maximum = maximum;
        self.max_log_volume_change = self.max_log_volume_change.clamp(minimum, maximum);
        self
    }

    #[inline]
    pub const fn max_log_volume_change(&self) -> f64 {
        self.max_log_volume_change
    }

    #[inline]
    pub fn acceptance_rate(&self) -> f64 {
        if self.total_attempts == 0 {
            0.0
        } else {
            self.total_accepted as f64 / self.total_attempts as f64
        }
    }

    pub fn propose(&self, rng: &mut impl Rng) -> ProposedMove<IsotropicVolumeChange> {
        let log_volume_ratio =
            rng.random_range(-self.max_log_volume_change..self.max_log_volume_change);
        // A symmetric random walk in ln V has q(V|V')/q(V'|V)=V'/V.
        ProposedMove::new(
            IsotropicVolumeChange::new(log_volume_ratio),
            log_volume_ratio,
        )
    }

    pub fn record_result(&mut self, accepted: bool) {
        self.window_attempts += 1;
        self.total_attempts += 1;
        if accepted {
            self.window_accepted += 1;
            self.total_accepted += 1;
        }
    }

    pub fn finish_sweep(&mut self, allows_adaptation: bool) {
        if !allows_adaptation {
            self.clear_window();
            return;
        }
        self.window_sweeps += 1;
        if self.window_sweeps < self.interval_sweeps || self.window_attempts == 0 {
            return;
        }
        let acceptance = self.window_accepted as f64 / self.window_attempts as f64;
        let factor = (self.gain * (acceptance - self.target_acceptance)).exp();
        self.max_log_volume_change =
            (self.max_log_volume_change * factor).clamp(self.minimum, self.maximum);
        self.clear_window();
    }

    fn clear_window(&mut self) {
        self.window_sweeps = 0;
        self.window_attempts = 0;
        self.window_accepted = 0;
    }
}

/// NPT kernel: one translation attempt per particle plus volume attempts.
#[derive(Debug, Clone)]
pub struct ParticleNptMetropolisCore<const D: usize> {
    translation: TranslateParticle<D>,
    volume: LogVolumeScale,
    pressure: f64,
    volume_attempts_per_sweep: u64,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    translation_patch: ParticleEnergyPatch,
    volume_patch: VolumeChangePatch<D>,
    energy_check_interval: u64,
    sweeps: u64,
}

impl<const D: usize> ParticleNptMetropolisCore<D> {
    pub fn new(max_displacement: f64, max_log_volume_change: f64, pressure: f64) -> Self {
        assert!(pressure.is_finite(), "pressure must be finite");
        Self {
            translation: TranslateParticle::new(max_displacement),
            volume: LogVolumeScale::new(max_log_volume_change),
            pressure,
            volume_attempts_per_sweep: 1,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            translation_patch: ParticleEnergyPatch::default(),
            volume_patch: VolumeChangePatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
        }
    }

    pub fn with_translation(mut self, translation: TranslateParticle<D>) -> Self {
        self.translation = translation;
        self
    }

    pub fn with_volume_proposal(mut self, volume: LogVolumeScale) -> Self {
        self.volume = volume;
        self
    }

    pub fn with_volume_attempts(mut self, attempts: u64) -> Self {
        self.volume_attempts_per_sweep = attempts;
        self
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }

    #[inline]
    pub const fn pressure(&self) -> f64 {
        self.pressure
    }

    #[inline]
    pub const fn translation(&self) -> &TranslateParticle<D> {
        &self.translation
    }

    #[inline]
    pub const fn volume(&self) -> &LogVolumeScale {
        &self.volume
    }
}

impl<const D: usize, P: PairPotential> ParticleAlgorithm<D, P> for ParticleNptMetropolisCore<D> {
    fn sweep_with_phase(
        &mut self,
        system: &mut ParticleSystem<D>,
        potential: &P,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        let target = IsothermalIsobaric::new(system.beta, self.pressure);
        let acceptance = MetropolisHastingsAcceptance;
        let particles = self.order.prepare(system.len(), self.visit_schedule, rng);
        for &particle in particles {
            let proposal = self
                .translation
                .propose(system.configuration(), particle, rng);
            let outcome = metropolis_hastings_step(
                system,
                potential,
                &proposal,
                &target,
                &acceptance,
                &mut self.translation_patch,
                rng,
            );
            self.translation.record_result(outcome.accepted);
        }
        for _ in 0..self.volume_attempts_per_sweep {
            let proposal = self.volume.propose(rng);
            let outcome = metropolis_hastings_step(
                system,
                potential,
                &proposal,
                &target,
                &acceptance,
                &mut self.volume_patch,
                rng,
            );
            self.volume.record_result(outcome.accepted);
        }
        self.translation.finish_sweep(phase.allows_adaptation());
        self.volume.finish_sweep(phase.allows_adaptation());

        self.sweeps = self.sweeps.wrapping_add(1);
        if self.energy_check_interval > 0 && self.sweeps.is_multiple_of(self.energy_check_interval)
        {
            system
                .validate(potential)
                .expect("NPT particle cache audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Particle NPT Metropolis-Hastings"
    }
}
