//! Single-particle translation moves and adaptive proposal scale.

use crate::core::trial::{ProposedMove, TrialEvaluator};
use crate::particle::{
    PairPotential, ParticleConfiguration, ParticleEnergyPatch, ParticleSystem, SimulationCell,
};
use rand::{Rng, RngExt};

/// Replace one particle position atomically.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParticleTranslation<const D: usize> {
    /// Particle index.
    pub particle: usize,
    /// Trial position; accepted commits wrap it into the primary cell.
    pub position: [f64; D],
}

impl<const D: usize> ParticleTranslation<D> {
    /// Create a translation move.
    #[inline]
    pub const fn new(particle: usize, position: [f64; D]) -> Self {
        Self { particle, position }
    }
}

/// Uniform-box single-particle translation proposal with warmup adaptation.
#[derive(Debug, Clone)]
pub struct TranslateParticle<const D: usize> {
    max_displacement: f64,
    minimum_displacement: f64,
    maximum_displacement: f64,
    target_acceptance: f64,
    adaptation_gain: f64,
    adaptation_interval_sweeps: u64,
    window_sweeps: u64,
    window_attempts: u64,
    window_accepted: u64,
    total_attempts: u64,
    total_accepted: u64,
}

impl<const D: usize> TranslateParticle<D> {
    /// Build a proposal with a uniform displacement in each Cartesian component.
    pub fn new(max_displacement: f64) -> Self {
        assert!(D > 0, "particle dimension must be positive");
        assert!(
            max_displacement.is_finite() && max_displacement > 0.0,
            "maximum displacement must be finite and positive"
        );
        let scaled_minimum = max_displacement * 1e-6;
        let minimum_displacement = if scaled_minimum > 0.0 {
            scaled_minimum
        } else {
            max_displacement
        };
        let scaled_maximum = max_displacement * 1e6;
        let maximum_displacement = if scaled_maximum.is_finite() {
            scaled_maximum
        } else {
            f64::MAX
        };
        Self {
            max_displacement,
            minimum_displacement,
            maximum_displacement,
            target_acceptance: 0.5,
            adaptation_gain: 0.5,
            adaptation_interval_sweeps: 20,
            window_sweeps: 0,
            window_attempts: 0,
            window_accepted: 0,
            total_attempts: 0,
            total_accepted: 0,
        }
    }

    /// Configure bounded multiplicative warmup adaptation.
    pub fn with_adaptation(
        mut self,
        target_acceptance: f64,
        interval_sweeps: u64,
        gain: f64,
        minimum_displacement: f64,
        maximum_displacement: f64,
    ) -> Self {
        assert!(
            target_acceptance.is_finite() && (0.0..1.0).contains(&target_acceptance),
            "target acceptance must lie strictly between zero and one"
        );
        assert!(interval_sweeps > 0, "adaptation interval must be positive");
        assert!(
            gain.is_finite() && gain > 0.0,
            "adaptation gain must be positive"
        );
        assert!(
            minimum_displacement.is_finite()
                && maximum_displacement.is_finite()
                && minimum_displacement > 0.0
                && maximum_displacement >= minimum_displacement,
            "invalid displacement adaptation bounds"
        );
        self.target_acceptance = target_acceptance;
        self.adaptation_interval_sweeps = interval_sweeps;
        self.adaptation_gain = gain;
        self.minimum_displacement = minimum_displacement;
        self.maximum_displacement = maximum_displacement;
        self.max_displacement = self
            .max_displacement
            .clamp(self.minimum_displacement, self.maximum_displacement);
        self
    }

    /// Current half-width of the uniform displacement box.
    #[inline]
    pub const fn max_displacement(&self) -> f64 {
        self.max_displacement
    }

    /// Acceptance fraction since statistics were last reset.
    #[inline]
    pub fn acceptance_rate(&self) -> f64 {
        if self.total_attempts == 0 {
            0.0
        } else {
            self.total_accepted as f64 / self.total_attempts as f64
        }
    }

    /// Attempt count since statistics were last reset.
    #[inline]
    pub const fn attempted(&self) -> u64 {
        self.total_attempts
    }

    /// Accepted count since statistics were last reset.
    #[inline]
    pub const fn accepted(&self) -> u64 {
        self.total_accepted
    }

    /// Reset reporting statistics without changing the proposal scale.
    pub fn reset_statistics(&mut self) {
        self.total_attempts = 0;
        self.total_accepted = 0;
    }

    /// Propose a symmetric wrapped translation for a selected particle.
    pub fn propose(
        &self,
        configuration: &ParticleConfiguration<D>,
        particle: usize,
        rng: &mut impl Rng,
    ) -> ProposedMove<ParticleTranslation<D>> {
        let mut position = *configuration.position(particle);
        for coordinate in &mut position {
            *coordinate += rng.random_range(-self.max_displacement..self.max_displacement);
        }
        configuration.cell().wrap(&mut position);
        ProposedMove::symmetric(ParticleTranslation::new(particle, position))
    }

    /// Record one trial result.
    pub fn record_result(&mut self, accepted: bool) {
        self.window_attempts += 1;
        self.total_attempts += 1;
        if accepted {
            self.window_accepted += 1;
            self.total_accepted += 1;
        }
    }

    /// Complete a sweep, adapting only when the run phase permits it.
    pub fn finish_sweep(&mut self, allows_adaptation: bool) {
        if !allows_adaptation {
            self.clear_adaptation_window();
            return;
        }
        self.window_sweeps += 1;
        if self.window_sweeps < self.adaptation_interval_sweeps || self.window_attempts == 0 {
            return;
        }
        let acceptance = self.window_accepted as f64 / self.window_attempts as f64;
        let factor = (self.adaptation_gain * (acceptance - self.target_acceptance)).exp();
        self.max_displacement = (self.max_displacement * factor)
            .clamp(self.minimum_displacement, self.maximum_displacement);
        self.clear_adaptation_window();
    }

    fn clear_adaptation_window(&mut self) {
        self.window_sweeps = 0;
        self.window_attempts = 0;
        self.window_accepted = 0;
    }
}

impl<const D: usize, P: PairPotential> TrialEvaluator<P, ParticleTranslation<D>>
    for ParticleSystem<D>
{
    type Delta = crate::core::ensemble::ThermodynamicDelta;
    type Patch = ParticleEnergyPatch;

    fn evaluate_trial(
        &self,
        model: &P,
        movement: &ParticleTranslation<D>,
        patch: &mut Self::Patch,
    ) -> Self::Delta {
        self.evaluate_translation(model, movement, patch)
    }

    fn commit_trial(&mut self, movement: &ParticleTranslation<D>, patch: &Self::Patch) {
        self.commit_translation(movement, patch);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warmup_adaptation_moves_scale_and_production_freezes_it() {
        let mut proposal = TranslateParticle::<3>::new(0.2).with_adaptation(0.5, 2, 1.0, 0.01, 1.0);
        for _ in 0..8 {
            proposal.record_result(true);
        }
        proposal.finish_sweep(true);
        proposal.finish_sweep(true);
        let adapted = proposal.max_displacement();
        assert!(adapted > 0.2);

        for _ in 0..8 {
            proposal.record_result(false);
        }
        proposal.finish_sweep(false);
        proposal.finish_sweep(false);
        assert_eq!(proposal.max_displacement(), adapted);
    }
}
