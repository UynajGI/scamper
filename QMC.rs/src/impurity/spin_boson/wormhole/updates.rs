//! Diagonal insertion/removal and wormhole directed-loop updates.

use carlo_rs::accept_log_probability;
use rand::Rng;
use rand::RngExt;

use crate::algorithm::{QmcKernel, UpdateSchedule};
use crate::impurity::core::estimator::{
    LoopSegment, SpinFlipOperator, TransverseCorrelationSample, TransverseLoopAccumulator,
};
use crate::impurity::core::imaginary_time::PropagationDirection;
use crate::impurity::spin_boson::model::ImpurityModel;
use crate::impurity::spin_boson::wormhole::configuration::{
    LegId, LegSide, Vertex, VertexId, WormholeConfiguration,
};
use crate::impurity::ImpurityError;

/// Proposal used to choose the first directed-loop leg.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LoopStartPolicy {
    /// Sample imaginary time and propagation direction uniformly. This policy
    /// is required for the on-the-fly transverse improved estimator.
    #[default]
    RandomTime,
    /// Sample one of the existing vertex legs uniformly.
    RandomLeg,
}

/// Accumulated update diagnostics.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WormholeUpdateStats {
    pub diagonal_proposals: u64,
    pub diagonal_add_accepts: u64,
    pub diagonal_remove_accepts: u64,
    pub loops: u64,
    pub loop_steps: u64,
    pub bounces: u64,
    pub wormholes: u64,
    pub same_endpoint_exits: u64,
    pub loop_aborts: u64,
}

impl WormholeUpdateStats {
    pub fn diagonal_acceptance(&self) -> f64 {
        if self.diagonal_proposals == 0 {
            return 0.0;
        }
        (self.diagonal_add_accepts + self.diagonal_remove_accepts) as f64
            / self.diagonal_proposals as f64
    }
    pub fn mean_loop_steps(&self) -> f64 {
        let attempts = self.loops + self.loop_aborts;
        if attempts == 0 {
            0.0
        } else {
            self.loop_steps as f64 / attempts as f64
        }
    }
    pub fn bounce_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            0.0
        } else {
            self.bounces as f64 / self.loop_steps as f64
        }
    }
    pub fn wormhole_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            0.0
        } else {
            self.wormholes as f64 / self.loop_steps as f64
        }
    }
    pub fn loop_abort_fraction(&self) -> f64 {
        let attempts = self.loops + self.loop_aborts;
        if attempts == 0 {
            0.0
        } else {
            self.loop_aborts as f64 / attempts as f64
        }
    }
}

/// Generic continuous-time spin-boson update engine.
#[derive(Debug, Clone, PartialEq)]
pub struct WormholeEngine {
    model: ImpurityModel,
    schedule: UpdateSchedule,
    stats: WormholeUpdateStats,
    validate_each_sweep: bool,
    loop_start_policy: LoopStartPolicy,
    transverse: TransverseLoopAccumulator,
}

impl WormholeEngine {
    pub fn new(model: ImpurityModel, schedule: UpdateSchedule) -> Self {
        Self {
            model,
            schedule,
            stats: WormholeUpdateStats::default(),
            validate_each_sweep: false,
            loop_start_policy: LoopStartPolicy::RandomTime,
            transverse: TransverseLoopAccumulator::default(),
        }
    }

    pub fn model(&self) -> &ImpurityModel {
        &self.model
    }
    pub fn schedule(&self) -> UpdateSchedule {
        self.schedule
    }
    pub fn set_schedule(&mut self, schedule: UpdateSchedule) {
        self.schedule = schedule;
    }
    pub fn set_validate_each_sweep(&mut self, enabled: bool) {
        self.validate_each_sweep = enabled;
    }
    pub fn set_loop_start_policy(&mut self, policy: LoopStartPolicy) {
        self.loop_start_policy = policy;
    }
    pub fn loop_start_policy(&self) -> LoopStartPolicy {
        self.loop_start_policy
    }
    pub fn set_transverse_bins(&mut self, bins: usize) {
        self.transverse = TransverseLoopAccumulator::new(bins);
    }
    pub fn transverse_bins(&self) -> usize {
        self.transverse.bins()
    }
    pub fn take_transverse_sample(&mut self, beta: f64) -> TransverseCorrelationSample {
        self.transverse.take_sample(beta)
    }
    pub fn clear_transverse_estimator(&mut self) {
        self.transverse.clear();
    }
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
            if self.loop_start_policy == LoopStartPolicy::RandomTime {
                for _ in 0..self.schedule.directed_loops {
                    self.transverse.commit_free_loop(configuration.beta());
                }
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
        let beta = configuration.beta();
        let (start, tail_tau, mut direction, tail_operator, mut head_operator, measure_path) =
            match self.loop_start_policy {
                LoopStartPolicy::RandomTime => {
                    let tau = rng.random::<f64>() * beta;
                    let forward = rng.random::<bool>();
                    let direction = if forward {
                        PropagationDirection::Forward
                    } else {
                        PropagationDirection::Backward
                    };
                    let start = configuration.start_leg_at_time(tau, forward)?;
                    let spin = configuration.spin_at(&self.model, tau)?;
                    let tail = if forward {
                        SpinFlipOperator::from_transition(spin, -spin)
                    } else {
                        SpinFlipOperator::from_transition(-spin, spin)
                    }
                    .expect("spin-1/2 discontinuity");
                    (start, tau, direction, tail, tail.opposite(), true)
                }
                LoopStartPolicy::RandomLeg => {
                    let start = configuration.random_leg(rng)?;
                    let direction = direction_for_leg(start);
                    (
                        start,
                        configuration.endpoint_time(start.endpoint)?,
                        direction,
                        SpinFlipOperator::Raise,
                        SpinFlipOperator::Lower,
                        false,
                    )
                }
            };

        let original_empty_spin = configuration.empty_spin();
        let mut current = start;
        let mut first = true;
        let mut steps = 0_usize;
        let mut kind_journal: Vec<(VertexId, usize)> = Vec::new();
        let mut measurement_journal: Vec<LoopSegment> = Vec::new();
        if measure_path {
            measurement_journal.push(LoopSegment {
                tail_tau,
                from_tau: tail_tau,
                to_tau: configuration.endpoint_time(start.endpoint)?,
                direction,
                normal: tail_operator != head_operator,
            });
        }

        loop {
            if steps >= limit {
                rollback_kinds(configuration, &self.model, kind_journal)?;
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
            let old_exit_spin = self
                .model
                .interaction(interaction_id)
                .kind(old_kind)
                .spin(choice.exit_leg);

            if !bounce {
                kind_journal.push((vertex_id, old_kind));
                configuration.set_kind(vertex_id, choice.new_kind, &self.model)?;
                let new_exit_spin = self
                    .model
                    .interaction(interaction_id)
                    .kind(choice.new_kind)
                    .spin(choice.exit_leg);
                let exit_side = LegId::from_local(vertex_id, choice.exit_leg).side;
                head_operator = head_after_exit(exit_side, old_exit_spin, new_exit_spin)?;
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

            let exit = LegId::from_local(vertex_id, choice.exit_leg);
            direction = direction_for_leg(exit);
            let exit_tau = configuration.endpoint_time(exit.endpoint)?;

            if first && bounce {
                if measure_path {
                    measurement_journal.push(LoopSegment {
                        tail_tau,
                        from_tau: exit_tau,
                        to_tau: tail_tau,
                        direction,
                        normal: tail_operator != head_operator,
                    });
                    self.transverse.commit_loop(beta, &measurement_journal);
                }
                self.stats.loops += 1;
                return Ok(true);
            }

            let next = match configuration.linked_leg(exit) {
                Ok(next) => next,
                Err(error) => {
                    rollback_kinds(configuration, &self.model, kind_journal)?;
                    return Err(error);
                }
            };
            if next == start {
                if measure_path {
                    measurement_journal.push(LoopSegment {
                        tail_tau,
                        from_tau: exit_tau,
                        to_tau: tail_tau,
                        direction,
                        normal: tail_operator != head_operator,
                    });
                }
                if let Err(error) = configuration.sync_empty_spin_from_worldline(&self.model) {
                    rollback_kinds(configuration, &self.model, kind_journal)?;
                    configuration.set_empty_spin(original_empty_spin);
                    return Err(error);
                }
                #[cfg(debug_assertions)]
                if configuration.validate(&self.model).is_err() {
                    rollback_kinds(configuration, &self.model, kind_journal)?;
                    configuration.set_empty_spin(original_empty_spin);
                    self.stats.loop_aborts += 1;
                    return Ok(false);
                }
                if measure_path {
                    self.transverse.commit_loop(beta, &measurement_journal);
                }
                self.stats.loops += 1;
                return Ok(true);
            }
            if measure_path {
                measurement_journal.push(LoopSegment {
                    tail_tau,
                    from_tau: exit_tau,
                    to_tau: configuration.endpoint_time(next.endpoint)?,
                    direction,
                    normal: tail_operator != head_operator,
                });
            }
            current = next;
            first = false;
        }
    }
}

fn direction_for_leg(leg: LegId) -> PropagationDirection {
    match leg.side {
        LegSide::Outgoing => PropagationDirection::Forward,
        LegSide::Incoming => PropagationDirection::Backward,
    }
}

fn head_after_exit(
    side: LegSide,
    old_spin: i8,
    new_spin: i8,
) -> Result<SpinFlipOperator, ImpurityError> {
    let operator = match side {
        LegSide::Outgoing => SpinFlipOperator::from_transition(new_spin, old_spin),
        LegSide::Incoming => SpinFlipOperator::from_transition(old_spin, new_spin),
    };
    operator.ok_or_else(|| {
        ImpurityError::InvalidConfiguration(
            "non-bounce scattering did not flip the selected exit leg".into(),
        )
    })
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

    use crate::impurity::spin_boson::bath::{Bath, SingleModeBath};

    use super::*;

    #[test]
    fn loop_limit_abort_restores_configuration_and_discards_estimator() {
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
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let closed = engine
                .one_directed_loop(&mut configuration, &mut rng, 1)
                .expect("loop attempt");
            if !closed {
                assert_eq!(configuration, baseline);
                let sample = engine.take_transverse_sample(2.0);
                assert_eq!(sample.completed_loops, 0);
                observed_abort = true;
                break;
            }
        }
        assert!(observed_abort, "test seeds did not exercise rollback");
    }

    #[test]
    fn empty_sector_has_exact_transverse_estimator() {
        let model = ImpurityModel::xxz(
            Bath::SingleMode(SingleModeBath::new(1.0).unwrap()),
            0.0,
            0.0,
            0.0,
            Some(1.0),
        )
        .unwrap();
        let mut engine = WormholeEngine::new(model, UpdateSchedule::new(0, 2, 16));
        let mut configuration = WormholeConfiguration::new(4.0, 1).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        engine.sweep(&mut configuration, &mut rng).unwrap();
        let sample = engine.take_transverse_sample(4.0);
        assert_eq!(sample.completed_loops, 2);
        assert!(sample
            .sampled_x
            .iter()
            .all(|value| (*value - 0.25).abs() < 1e-14));
    }
}
