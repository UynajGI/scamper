//! Diagonal insertion/removal and wormhole directed-loop updates.

use carlo_rs::accept_log_probability;
use rand::Rng;
use rand::RngExt;

use crate::algorithm::{QmcKernel, UpdateSchedule};

use super::configuration::WormholeConfiguration;
use super::error::ImpurityError;
use super::model::ImpurityModel;
use super::vertex::{LegId, Vertex, VertexId};

/// Proposal used to choose the first directed-loop leg.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LoopStartPolicy {
    /// Sample imaginary time and propagation direction uniformly.
    #[default]
    RandomTime,
    /// Sample one of the existing vertex legs uniformly.
    RandomLeg,
}

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
    /// Loops rolled back after reaching the safety limit.
    pub loop_aborts: u64,
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

    /// Mean number of local steps per attempted loop.
    pub fn mean_loop_steps(&self) -> f64 {
        let attempts = self.loops + self.loop_aborts;
        if attempts == 0 {
            return 0.0;
        }
        self.loop_steps as f64 / attempts as f64
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

    /// Fraction of attempted loops that were rolled back.
    pub fn loop_abort_fraction(&self) -> f64 {
        let attempts = self.loops + self.loop_aborts;
        if attempts == 0 {
            return 0.0;
        }
        self.loop_aborts as f64 / attempts as f64
    }
}

/// Generic continuous-time impurity update engine.
#[derive(Debug, Clone, PartialEq)]
pub struct WormholeEngine {
    model: ImpurityModel,
    schedule: UpdateSchedule,
    stats: WormholeUpdateStats,
    validate_each_sweep: bool,
    loop_start_policy: LoopStartPolicy,
}

impl WormholeEngine {
    /// Construct an engine for one model catalog.
    pub fn new(model: ImpurityModel, schedule: UpdateSchedule) -> Self {
        Self {
            model,
            schedule,
            stats: WormholeUpdateStats::default(),
            validate_each_sweep: false,
            loop_start_policy: LoopStartPolicy::RandomTime,
        }
    }

    /// Model catalog.
    pub fn model(&self) -> &ImpurityModel {
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

    /// Select the directed-loop start proposal.
    pub fn set_loop_start_policy(&mut self, policy: LoopStartPolicy) {
        self.loop_start_policy = policy;
    }

    /// Current directed-loop start proposal.
    pub fn loop_start_policy(&self) -> LoopStartPolicy {
        self.loop_start_policy
    }

    /// Update statistics.
    pub fn stats(&self) -> WormholeUpdateStats {
        self.stats
    }

    fn diagonal_update_block<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        rng: &mut R,
    ) -> Result<(), ImpurityError> {
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
    ) -> Result<(), ImpurityError> {
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
    ) -> Result<(), ImpurityError> {
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
    ) -> Result<(), ImpurityError> {
        if configuration.expansion_order() == 0 {
            if rng.random::<bool>() {
                configuration.set_empty_spin(-configuration.empty_spin());
            }
            return Ok(());
        }

        for _ in 0..self.schedule.directed_loops {
            let leg_count = 4 * configuration.expansion_order();
            let limit = (self.schedule.max_loop_steps_factor * leg_count).max(32);
            self.one_directed_loop(configuration, rng, limit)?;
        }
        Ok(())
    }

    fn one_directed_loop<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut WormholeConfiguration,
        rng: &mut R,
        limit: usize,
    ) -> Result<bool, ImpurityError> {
        let start = match self.loop_start_policy {
            LoopStartPolicy::RandomTime => {
                let tau = rng.random::<f64>() * configuration.beta();
                configuration.start_leg_at_time(tau, rng.random::<bool>())?
            }
            LoopStartPolicy::RandomLeg => configuration.random_leg(rng)?,
        };
        let original_empty_spin = configuration.empty_spin();
        let mut current = start;
        let mut first = true;
        let mut steps = 0_usize;
        let mut journal: Vec<(VertexId, usize)> = Vec::new();

        loop {
            if steps >= limit {
                rollback_kinds(configuration, &self.model, journal)?;
                self.stats.loop_aborts += 1;
                return Ok(false);
            }

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
            let bounce = choice.exit_leg == entrance;

            if !bounce {
                journal.push((vertex_id, old_kind));
                configuration.set_kind(vertex_id, choice.new_kind, &self.model)?;
            }

            self.stats.loop_steps += 1;
            steps += 1;
            if bounce {
                self.stats.bounces += 1;
            } else if choice.exit_leg / 2 != entrance / 2 {
                self.stats.wormholes += 1;
            } else {
                self.stats.same_endpoint_exits += 1;
            }

            if first && bounce {
                self.stats.loops += 1;
                return Ok(true);
            }

            let exit = LegId::from_local(vertex_id, choice.exit_leg);
            let next = match configuration.linked_leg(exit) {
                Ok(next) => next,
                Err(error) => {
                    rollback_kinds(configuration, &self.model, journal)?;
                    return Err(error);
                }
            };
            if next == start {
                if let Err(error) = configuration.sync_empty_spin_from_worldline(&self.model) {
                    rollback_kinds(configuration, &self.model, journal)?;
                    configuration.set_empty_spin(original_empty_spin);
                    return Err(error);
                }
                #[cfg(debug_assertions)]
                if configuration.validate(&self.model).is_err() {
                    rollback_kinds(configuration, &self.model, journal)?;
                    configuration.set_empty_spin(original_empty_spin);
                    self.stats.loop_aborts += 1;
                    return Ok(false);
                }
                self.stats.loops += 1;
                return Ok(true);
            }
            current = next;
            first = false;
        }
    }
}

fn rollback_kinds(
    configuration: &mut WormholeConfiguration,
    model: &ImpurityModel,
    journal: Vec<(VertexId, usize)>,
) -> Result<(), ImpurityError> {
    for (vertex, old_kind) in journal.into_iter().rev() {
        configuration.set_kind(vertex, old_kind, model)?;
    }
    Ok(())
}

impl<R: Rng + ?Sized> QmcKernel<WormholeConfiguration, R> for WormholeEngine {
    type Error = ImpurityError;
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

    use crate::impurity::bath::{Bath, SingleModeBath};
    use crate::impurity::model::ImpurityModel;
    use crate::impurity::vertex::Vertex;

    use super::*;

    #[test]
    fn loop_limit_abort_restores_the_configuration() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xyz(bath, 0.8, 0.2, 0.0, 0.0, Some(0.4)).expect("model");
        let diagonal_kind = model.interaction(0).diagonal_kind(1, 1);
        let mut baseline = WormholeConfiguration::new(2.0, 1).expect("configuration");
        for (tau_a, tau_b) in [(0.2, 0.7), (1.1, 1.6)] {
            baseline
                .insert_vertex(
                    Vertex {
                        tau_a,
                        tau_b,
                        omega: 1.0,
                        interaction: 0,
                        kind: diagonal_kind,
                    },
                    &model,
                )
                .expect("insert vertex");
        }
        baseline.validate(&model).expect("valid baseline");

        let mut observed_abort = false;
        for seed in 0..512 {
            let mut configuration = baseline.clone();
            let mut engine = WormholeEngine::new(model.clone(), UpdateSchedule::new(0, 1, 1));
            engine.set_loop_start_policy(LoopStartPolicy::RandomLeg);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let closed = engine
                .one_directed_loop(&mut configuration, &mut rng, 1)
                .expect("loop attempt");
            if !closed {
                assert_eq!(configuration, baseline);
                configuration
                    .validate(&model)
                    .expect("rollback preserves worldline");
                assert_eq!(engine.stats().loop_aborts, 1);
                observed_abort = true;
                break;
            }
        }
        assert!(
            observed_abort,
            "test seeds did not exercise the rollback path"
        );
    }

    #[test]
    fn mixed_updates_preserve_worldline() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
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
