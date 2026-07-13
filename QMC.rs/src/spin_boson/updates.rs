//! Diagonal insertion/removal and wormhole directed-loop updates.

use carlo_rs::accept_log_probability;
use rand::Rng;
use rand::RngExt;

use crate::algorithm::{QmcKernel, UpdateSchedule};

use super::configuration::WormholeConfiguration;
use super::error::SpinBosonError;
use super::model::SpinBosonModel;
use super::vertex::{LegId, Vertex, VertexId};

/// Accumulated update diagnostics.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WormholeUpdateStats {
    /// Diagonal add/remove proposals.
    pub diagonal_proposals: u64,
    /// Accepted insertions.
    pub diagonal_add_accepts: u64,
    /// Accepted removals.
    pub diagonal_remove_accepts: u64,
    /// Directed loops started and closed.
    pub loops: u64,
    /// Total local scattering steps.
    pub loop_steps: u64,
    /// Bounce steps.
    pub bounces: u64,
    /// Exits through the other endpoint of a retarded vertex.
    pub wormholes: u64,
    /// Non-bounce exits through the same endpoint.
    pub same_endpoint_exits: u64,
}

impl WormholeUpdateStats {
    /// Accepted diagonal moves divided by proposals.
    pub fn diagonal_acceptance(&self) -> f64 {
        if self.diagonal_proposals == 0 {
            return 0.0;
        }
        (self.diagonal_add_accepts + self.diagonal_remove_accepts) as f64
            / self.diagonal_proposals as f64
    }

    /// Mean number of local steps per loop.
    pub fn mean_loop_steps(&self) -> f64 {
        if self.loops == 0 {
            return 0.0;
        }
        self.loop_steps as f64 / self.loops as f64
    }

    /// Fraction of loop steps that bounce.
    pub fn bounce_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            return 0.0;
        }
        self.bounces as f64 / self.loop_steps as f64
    }

    /// Fraction of loop steps that traverse a retarded wormhole.
    pub fn wormhole_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            return 0.0;
        }
        self.wormholes as f64 / self.loop_steps as f64
    }
}

/// Generic continuous-time spin-boson update engine.
#[derive(Debug, Clone, PartialEq)]
pub struct WormholeEngine {
    model: SpinBosonModel,
    schedule: UpdateSchedule,
    stats: WormholeUpdateStats,
    validate_each_sweep: bool,
}

impl WormholeEngine {
    /// Construct an engine for one model catalog.
    pub fn new(model: SpinBosonModel, schedule: UpdateSchedule) -> Self {
        Self {
            model,
            schedule,
            stats: WormholeUpdateStats::default(),
            validate_each_sweep: false,
        }
    }

    /// Model catalog.
    pub fn model(&self) -> &SpinBosonModel {
        &self.model
    }

    /// Fixed work schedule.
    pub fn schedule(&self) -> UpdateSchedule {
        self.schedule
    }

    /// Replace the fixed work schedule.
    pub fn set_schedule(&mut self, schedule: UpdateSchedule) {
        self.schedule = schedule;
    }

    /// Enable expensive invariant checking after every sweep.
    pub fn set_validate_each_sweep(&mut self, enabled: bool) {
        self.validate_each_sweep = enabled;
    }

    /// Update statistics.
    pub fn stats(&self) -> WormholeUpdateStats {
        self.stats
    }

    fn diagonal_update_block<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        rng: &mut R,
    ) -> Result<(), SpinBosonError> {
        for _ in 0..self.schedule.diagonal_proposals {
            self.stats.diagonal_proposals += 1;
            let diagonal_order = configuration.diagonal_order();
            if rng.random::<bool>() {
                self.try_add_diagonal(configuration, diagonal_order, rng)?;
            } else if diagonal_order > 0 {
                self.try_remove_diagonal(configuration, diagonal_order, rng)?;
            }
        }
        Ok(())
    }

    fn try_add_diagonal<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        diagonal_order: usize,
        rng: &mut R,
    ) -> Result<(), SpinBosonError> {
        let interaction_id = rng.random_range(0..self.model.interaction_count());
        let interaction = self.model.interaction(interaction_id);
        let tau_a = rng.random::<f64>() * configuration.beta();
        let sample =
            interaction
                .bath()
                .sample(configuration.beta(), interaction.direction(), rng)?;
        let tau_b = (tau_a - sample.delta_tau).rem_euclid(configuration.beta());
        let spin_a = configuration.spin_before(&self.model, tau_a)?;
        let spin_b = configuration.spin_before(&self.model, tau_b)?;
        let kind = interaction.diagonal_kind(spin_a, spin_b);
        let weight = interaction.kind(kind).weight();
        let interaction_probability = 1.0 / self.model.interaction_count() as f64;
        let log_acceptance = configuration.beta().ln() + weight.ln()
            - ((diagonal_order + 1) as f64).ln()
            - interaction_probability.ln();
        if accept_log_probability(log_acceptance, rng) {
            configuration.insert_vertex(
                Vertex {
                    tau_a,
                    tau_b,
                    omega: sample.omega,
                    interaction: interaction_id,
                    kind,
                },
                &self.model,
            )?;
            self.stats.diagonal_add_accepts += 1;
        }
        Ok(())
    }

    fn try_remove_diagonal<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        diagonal_order: usize,
        rng: &mut R,
    ) -> Result<(), SpinBosonError> {
        let selected = configuration.random_diagonal_vertex(rng)?;
        let vertex = configuration.vertex(selected)?;
        let interaction = self.model.interaction(vertex.interaction);
        let weight = interaction.kind(vertex.kind).weight();
        let interaction_probability = 1.0 / self.model.interaction_count() as f64;
        let log_acceptance = (diagonal_order as f64).ln() + interaction_probability.ln()
            - configuration.beta().ln()
            - weight.ln();
        if accept_log_probability(log_acceptance, rng) {
            configuration.remove_vertex(selected)?;
            self.stats.diagonal_remove_accepts += 1;
        }
        Ok(())
    }

    fn directed_loop_block<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        rng: &mut R,
    ) -> Result<(), SpinBosonError> {
        if configuration.expansion_order() == 0 {
            if rng.random::<bool>() {
                configuration.set_empty_spin(-configuration.empty_spin());
            }
            return Ok(());
        }

        for _ in 0..self.schedule.directed_loops {
            let start = configuration.random_leg(rng)?;
            let mut current = start;
            let leg_count = 4 * configuration.expansion_order();
            let limit = (self.schedule.max_loop_steps_factor * leg_count).max(32);
            let mut steps = 0_usize;

            // Rollback journal: records (vertex_id, old_kind) for each kind change.
            let mut journal: Vec<(VertexId, usize)> = Vec::new();

            loop {
                let vertex_id = current.endpoint.vertex;
                let entrance = current.local_leg();

                let vertex = configuration.vertex(vertex_id)?;
                let interaction_id = vertex.interaction;
                let old_kind = vertex.kind;

                let choice = self
                    .model
                    .interaction(interaction_id)
                    .scattering()
                    .sample(old_kind, entrance, rng);

                // Record before modifying.
                journal.push((vertex_id, old_kind));
                configuration.set_kind(vertex_id, choice.new_kind, &self.model)?;

                let exit = LegId::from_local(vertex_id, choice.exit_leg);

                self.stats.loop_steps += 1;
                steps += 1;
                if choice.exit_leg == entrance {
                    self.stats.bounces += 1;
                } else if choice.exit_leg / 2 != entrance / 2 {
                    self.stats.wormholes += 1;
                } else {
                    self.stats.same_endpoint_exits += 1;
                }

                let next = configuration.linked_leg(exit)?;
                if next == start {
                    break;
                }
                current = next;

                if steps > limit {
                    // Rollback all kind changes.
                    for (vid, old_kind) in journal.into_iter().rev() {
                        configuration.set_kind(vid, old_kind, &self.model)?;
                    }
                    return Err(SpinBosonError::LoopDidNotClose { steps, limit });
                }
            }
            self.stats.loops += 1;
        }
        Ok(())
    }
}

impl<R: Rng + ?Sized> QmcKernel<WormholeConfiguration, R> for WormholeEngine {
    type Error = SpinBosonError;
    type Diagnostics = WormholeUpdateStats;

    fn sweep(
        &mut self,
        configuration: &mut WormholeConfiguration,
        rng: &mut R,
    ) -> Result<(), Self::Error> {
        self.diagonal_update_block(configuration, rng)?;
        self.directed_loop_block(configuration, rng)?;
        if self.validate_each_sweep {
            <Self as QmcKernel<WormholeConfiguration, R>>::validate(self, configuration)?;
        }
        Ok(())
    }

    fn validate(&self, configuration: &WormholeConfiguration) -> Result<(), Self::Error> {
        configuration.validate(&self.model)
    }

    fn diagnostics(&self) -> &Self::Diagnostics {
        &self.stats
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use crate::spin_boson::bath::{Bath, SingleModeBath};
    use crate::spin_boson::model::SpinBosonModel;

    use super::*;

    #[test]
    fn mixed_updates_preserve_worldline() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = SpinBosonModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut engine = WormholeEngine::new(model, UpdateSchedule::new(4, 2, 32));
        engine.set_validate_each_sweep(true);
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(19);
        for _ in 0..500 {
            engine.sweep(&mut configuration, &mut rng).expect("sweep");
        }
        <WormholeEngine as QmcKernel<WormholeConfiguration, Xoshiro256PlusPlus>>::validate(
            &engine,
            &configuration,
        )
        .expect("valid configuration");
        assert!(engine.stats().loops > 0);
    }
}
