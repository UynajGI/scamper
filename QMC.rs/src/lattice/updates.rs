//! Continuous-time diagonal and directed-loop updates.

use rand::Rng;
use rand::RngExt;

use crate::algorithm::{QmcKernel, UpdateSchedule};
use crate::local_space::{BasisState, LocalHilbertSpace};

use super::configuration::{LatticeConfiguration, WorldlineIndex};
use super::error::LatticeQmcError;
use super::model::{PositiveOperatorModel, SpinLatticeModel};
use super::vertex::Vertex;

/// Accumulated update diagnostics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LatticeUpdateStats {
    /// Diagonal insertion proposals.
    pub diagonal_add_proposals: u64,
    /// Accepted insertions.
    pub diagonal_add_accepts: u64,
    /// Diagonal removal proposals.
    pub diagonal_remove_proposals: u64,
    /// Accepted removals.
    pub diagonal_remove_accepts: u64,
    /// Closed directed loops.
    pub loops: u64,
    /// Local scattering events.
    pub loop_steps: u64,
    /// Bounce events.
    pub bounces: u64,
    /// Exits to another site of the same bond operator.
    pub spatial_exits: u64,
    /// Whole empty worldlines changed without a vertex.
    pub free_worldline_moves: u64,
    /// Proposed loops rolled back because they hit the tail with incompatible
    /// charge or exceeded the configured safety limit.
    pub aborted_loops: u64,
}

impl LatticeUpdateStats {
    /// Combined diagonal acceptance rate.
    pub fn diagonal_acceptance(&self) -> f64 {
        let proposals = self.diagonal_add_proposals + self.diagonal_remove_proposals;
        let accepts = self.diagonal_add_accepts + self.diagonal_remove_accepts;
        if proposals == 0 {
            0.0
        } else {
            accepts as f64 / proposals as f64
        }
    }

    /// Mean local scatterings per closed loop.
    pub fn mean_loop_steps(&self) -> f64 {
        if self.loops == 0 {
            0.0
        } else {
            self.loop_steps as f64 / self.loops as f64
        }
    }

    /// Fraction of local scatterings that bounce.
    pub fn bounce_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            0.0
        } else {
            self.bounces as f64 / self.loop_steps as f64
        }
    }

    /// Fraction of local scatterings crossing a spatial bond.
    pub fn spatial_exit_fraction(&self) -> f64 {
        if self.loop_steps == 0 {
            0.0
        } else {
            self.spatial_exits as f64 / self.loop_steps as f64
        }
    }
}

/// Generic continuous-time directed-loop engine for positive sparse operators.
#[derive(Debug, Clone)]
pub struct ContinuousLatticeEngine<M = SpinLatticeModel> {
    model: M,
    schedule: UpdateSchedule,
    stats: LatticeUpdateStats,
    validate_each_sweep: bool,
    strict_loop_limits: bool,
    diagonal_scratch: Vec<usize>,
}

impl<M: PositiveOperatorModel> ContinuousLatticeEngine<M> {
    /// Construct an engine.
    pub fn new(model: M, schedule: UpdateSchedule) -> Self {
        Self {
            model,
            schedule,
            stats: LatticeUpdateStats::default(),
            validate_each_sweep: false,
            strict_loop_limits: false,
            diagonal_scratch: Vec::new(),
        }
    }

    /// Compiled model.
    pub fn model(&self) -> &M {
        &self.model
    }

    /// Current fixed sweep schedule.
    pub fn schedule(&self) -> UpdateSchedule {
        self.schedule
    }

    /// Replace the sweep schedule. Adapt only during thermalization.
    pub fn set_schedule(&mut self, schedule: UpdateSchedule) {
        self.schedule = schedule;
    }

    /// Enable expensive invariant checks after every sweep.
    pub fn set_validate_each_sweep(&mut self, enabled: bool) {
        self.validate_each_sweep = enabled;
    }

    /// Convert a loop safety-limit rollback into a hard runtime error.
    pub fn set_strict_loop_limits(&mut self, enabled: bool) {
        self.strict_loop_limits = enabled;
    }

    /// Update diagnostics.
    pub fn stats(&self) -> &LatticeUpdateStats {
        &self.stats
    }

    fn diagonal_update_block<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        rng: &mut R,
    ) -> Result<(), LatticeQmcError> {
        if self.model.term_count() == 0 {
            return Ok(());
        }
        let mut index = WorldlineIndex::build(configuration, &self.model)?;
        for _ in 0..self.schedule.diagonal_proposals {
            let diagonal_order = configuration.diagonal_order(&self.model);
            // Keep the add/remove selector probability fixed at one half even
            // at expansion order zero. Turning the zero-order removal branch
            // into an automatic insertion changes the proposal ratio at the
            // boundary and biases the sampled expansion order.
            if rng.random::<bool>() {
                if diagonal_order > 0
                    && self.try_remove_diagonal(configuration, diagonal_order, rng)?
                {
                    index = WorldlineIndex::build(configuration, &self.model)?;
                }
            } else {
                self.try_add_diagonal(configuration, &index, diagonal_order, rng)?;
            }
        }
        Ok(())
    }

    fn try_add_diagonal<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        index: &WorldlineIndex,
        diagonal_order: usize,
        rng: &mut R,
    ) -> Result<(), LatticeQmcError> {
        self.stats.diagonal_add_proposals += 1;
        let Some(term_id) = self.model.sample_term(rng) else {
            return Ok(());
        };
        let term = self.model.term(term_id);
        let tau = rng.random::<f64>() * configuration.beta();
        let states: Vec<BasisState> = term
            .sites()
            .iter()
            .map(|&site| index.state_before(configuration, &self.model, site, tau))
            .collect();
        let kind_id =
            term.diagonal_kind(&states)
                .ok_or_else(|| LatticeQmcError::MissingDiagonalVertex {
                    term: term_id,
                    states: states.clone(),
                })?;
        let weight = term.kind(kind_id).weight();
        let proposal_probability = self.model.term_probability(term_id);
        let ratio =
            configuration.beta() * weight / ((diagonal_order + 1) as f64 * proposal_probability);
        if rng.random::<f64>() < ratio.min(1.0) {
            configuration.vertices_mut().push(Vertex {
                tau,
                term: term_id,
                kind: kind_id,
            });
            self.stats.diagonal_add_accepts += 1;
        }
        Ok(())
    }

    fn try_remove_diagonal<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        diagonal_order: usize,
        rng: &mut R,
    ) -> Result<bool, LatticeQmcError> {
        self.stats.diagonal_remove_proposals += 1;
        self.diagonal_scratch.clear();
        self.diagonal_scratch
            .extend(configuration.vertices().iter().enumerate().filter_map(
                |(vertex_id, vertex)| {
                    self.model
                        .term(vertex.term)
                        .kind(vertex.kind)
                        .is_diagonal()
                        .then_some(vertex_id)
                },
            ));
        let selected = self.diagonal_scratch[rng.random_range(0..self.diagonal_scratch.len())];
        let vertex = &configuration.vertices()[selected];
        let term = self.model.term(vertex.term);
        let weight = term.kind(vertex.kind).weight();
        let proposal_probability = self.model.term_probability(vertex.term);
        let ratio = diagonal_order as f64 * proposal_probability / (configuration.beta() * weight);
        if rng.random::<f64>() < ratio.min(1.0) {
            configuration.vertices_mut().swap_remove(selected);
            self.stats.diagonal_remove_accepts += 1;
            return Ok(true);
        }
        Ok(false)
    }

    fn directed_loop_block<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        rng: &mut R,
    ) -> Result<(), LatticeQmcError> {
        let index = WorldlineIndex::build(configuration, &self.model)?;
        self.update_eventless_worldlines(configuration, &index, rng);
        if index.leg_count() == 0 {
            return Ok(());
        }
        for _ in 0..self.schedule.directed_loops {
            self.directed_loop(configuration, &index, rng)?;
        }
        Ok(())
    }

    fn update_eventless_worldlines<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        index: &WorldlineIndex,
        rng: &mut R,
    ) {
        for _ in 0..self.schedule.directed_loops {
            let site = rng.random_range(0..self.model.graph().site_count());
            if !index.events(site).is_empty() {
                continue;
            }
            let dimension = self.model.space().dimension(site);
            configuration.initial_states_mut()[site] = rng.random_range(0..dimension) as BasisState;
            self.stats.free_worldline_moves += 1;
        }
    }

    fn directed_loop<R: Rng + ?Sized>(
        &mut self,
        configuration: &mut LatticeConfiguration,
        index: &WorldlineIndex,
        rng: &mut R,
    ) -> Result<(), LatticeQmcError> {
        let (start_leg, start_delta) = self.choose_loop_start(configuration, index, rng)?;
        let mut current_leg = start_leg;
        let mut delta = start_delta;
        let limit = (self.schedule.max_loop_steps_factor * index.leg_count()).max(64);
        let mut journal = Vec::<(usize, usize)>::new();
        // Self-intersecting worms are legal in specialized formulations only
        // when every discontinuity on an already modified segment is tracked
        // explicitly. This generic sparse-operator backend instead samples
        // simple loops and rejects a proposal before it revisits a leg. The
        // rejection is reversal symmetric and avoids spin-1/2 assumptions in
        // arbitrary Spin-S local spaces.
        let mut visited_legs = vec![false; index.leg_count()];
        visited_legs[start_leg] = true;

        for steps in 1..=limit {
            let vertex_id = index.vertex_of_leg(current_leg);
            let entrance = index.local_leg(current_leg);
            let vertex = &configuration.vertices()[vertex_id];
            let term_id = vertex.term;
            let old_kind = vertex.kind;
            let Some(choice) = self
                .model
                .term(term_id)
                .scattering()
                .sample(old_kind, entrance, delta, rng)
            else {
                rollback(configuration, &journal);
                return Err(LatticeQmcError::InvalidConfiguration(format!(
                    "missing scattering row for term {term_id}, kind {old_kind}, leg {entrance}, delta {delta}"
                )));
            };
            let exit_global = index.vertex_offset(vertex_id) + choice.exit_leg;
            if exit_global != current_leg && visited_legs[exit_global] {
                rollback(configuration, &journal);
                self.stats.aborted_loops += 1;
                return Ok(());
            }

            journal.push((vertex_id, old_kind));
            configuration.vertices_mut()[vertex_id].kind = choice.new_kind;
            visited_legs[exit_global] = true;
            let next_leg = index.linked_leg(exit_global);

            self.stats.loop_steps += 1;
            let is_bounce = choice.new_kind == old_kind
                && choice.exit_leg == entrance
                && choice.next_delta == delta;
            if is_bounce {
                self.stats.bounces += 1;
            }
            if choice.exit_leg / 2 != entrance / 2 {
                self.stats.spatial_exits += 1;
            }

            if next_leg == start_leg {
                if choice.next_delta == start_delta {
                    index.canonicalize_initial_states(configuration, &self.model);
                    self.stats.loops += 1;
                } else {
                    rollback(configuration, &journal);
                    self.stats.aborted_loops += 1;
                }
                return Ok(());
            }
            if visited_legs[next_leg] {
                rollback(configuration, &journal);
                self.stats.aborted_loops += 1;
                return Ok(());
            }
            visited_legs[next_leg] = true;
            current_leg = next_leg;
            delta = choice.next_delta;

            if steps == limit {
                rollback(configuration, &journal);
                self.stats.aborted_loops += 1;
                if self.strict_loop_limits {
                    return Err(LatticeQmcError::LoopDidNotClose { steps, limit });
                }
                return Ok(());
            }
        }
        unreachable!("loop range always returns on its last step")
    }

    fn choose_loop_start<R: Rng + ?Sized>(
        &self,
        configuration: &LatticeConfiguration,
        index: &WorldlineIndex,
        rng: &mut R,
    ) -> Result<(usize, i8), LatticeQmcError> {
        let attempts = (4 * index.leg_count()).max(32);
        for _ in 0..attempts {
            let leg = rng.random_range(0..index.leg_count());
            let vertex = &configuration.vertices()[index.vertex_of_leg(leg)];
            let term = self.model.term(vertex.term);
            let entrance = index.local_leg(leg);
            let first_delta = if rng.random::<bool>() { 1 } else { -1 };
            for delta in [first_delta, -first_delta] {
                if term.scattering().has_row(vertex.kind, entrance, delta) {
                    return Ok((leg, delta));
                }
            }
        }
        Err(LatticeQmcError::NoLoopStart)
    }
}

fn rollback(configuration: &mut LatticeConfiguration, journal: &[(usize, usize)]) {
    for &(vertex, old_kind) in journal.iter().rev() {
        configuration.vertices_mut()[vertex].kind = old_kind;
    }
}

impl<M, R> QmcKernel<LatticeConfiguration, R> for ContinuousLatticeEngine<M>
where
    M: PositiveOperatorModel,
    R: Rng + ?Sized,
{
    type Error = LatticeQmcError;
    type Diagnostics = LatticeUpdateStats;

    fn sweep(
        &mut self,
        configuration: &mut LatticeConfiguration,
        rng: &mut R,
    ) -> Result<(), Self::Error> {
        self.diagonal_update_block(configuration, rng)?;
        self.directed_loop_block(configuration, rng)?;
        if self.validate_each_sweep {
            configuration.validate(&self.model)?;
        }
        Ok(())
    }

    fn validate(&self, configuration: &LatticeConfiguration) -> Result<(), Self::Error> {
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

    use crate::graph::CsrGraph;
    use crate::lattice::model::SpinLatticeModel;

    use super::*;

    #[test]
    fn arbitrary_spin_updates_preserve_worldlines() {
        let graph = CsrGraph::chain(6, true).expect("graph");
        let model = SpinLatticeModel::xxz(graph, 4, -0.7, 0.3).expect("model");
        let mut configuration =
            LatticeConfiguration::new(6.0, vec![2; 6], &model).expect("configuration");
        let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
        engine.set_validate_each_sweep(true);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(2026);
        for _ in 0..200 {
            engine.sweep(&mut configuration, &mut rng).expect("sweep");
        }
        configuration.validate(engine.model()).expect("valid");
        assert!(engine.stats().loops > 0);
    }
}
