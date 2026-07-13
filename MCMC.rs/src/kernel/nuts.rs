use rand::{Rng, RngExt};
use serde::{Deserialize, Serialize};

use crate::adaptation::{HmcWarmup, MetricAdaptation, MetricUpdate, WarmupWindowConfig};
use crate::integrator::{LeapfrogIntegrator, PhasePoint};
use crate::metric::{DenseMetric, DiagonalMetric, Metric, MetricKind, UnitMetric};
use crate::{
    DifferentiableLogDensity, EuclideanState, McmcError, SamplingPhase, TransitionKernel,
    TransitionReport,
};

const MAX_SUPPORTED_TREE_DEPTH: u16 = 20;

/// Multinomial No-U-Turn Sampler with Euclidean metric adaptation.
///
/// A transition doubles a binary Hamiltonian trajectory until it detects a
/// U-turn, encounters a divergent numerical path, or exhausts
/// `max_tree_depth`. Candidate states are sampled in the log domain with
/// weights proportional to `exp(-H)`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Nuts<M> {
    metric: M,
    step_size: f64,
    max_tree_depth: u16,
    max_energy_error: f64,
    warmup: Option<HmcWarmup>,
    #[serde(skip, default)]
    initial: PhasePoint,
    #[serde(skip, default)]
    current_gradient: Vec<f64>,
    #[serde(skip, default)]
    integrator: LeapfrogIntegrator,
}

impl<M> Nuts<M>
where
    M: Metric,
{
    pub fn new(metric: M, step_size: f64, max_tree_depth: u16) -> Result<Self, McmcError> {
        validate_nuts_config(metric.dimension(), step_size, max_tree_depth, 1_000.0)?;
        let dimension = metric.dimension();
        Ok(Self {
            metric,
            step_size,
            max_tree_depth,
            max_energy_error: 1_000.0,
            warmup: None,
            initial: PhasePoint::with_dimension(dimension),
            current_gradient: vec![0.0; dimension],
            integrator: LeapfrogIntegrator::with_dimension(dimension),
        })
    }

    pub fn with_max_energy_error(mut self, max_energy_error: f64) -> Result<Self, McmcError> {
        if !max_energy_error.is_finite() || max_energy_error <= 0.0 {
            return Err(McmcError::InvalidConfig(
                "maximum NUTS energy error must be finite and positive".to_string(),
            ));
        }
        self.max_energy_error = max_energy_error;
        Ok(self)
    }

    pub fn with_dual_averaging(
        self,
        total_warmup: u64,
        target_acceptance: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::None,
            windows,
        )
    }

    pub fn with_warmup_adaptation(
        mut self,
        total_warmup: u64,
        target_acceptance: f64,
        metric_adaptation: MetricAdaptation,
        windows: WarmupWindowConfig,
    ) -> Result<Self, McmcError> {
        validate_metric_adaptation(self.metric.kind(), metric_adaptation)?;
        self.warmup = Some(HmcWarmup::new(
            self.metric.dimension(),
            total_warmup,
            self.step_size,
            target_acceptance,
            metric_adaptation,
            windows,
        )?);
        Ok(self)
    }

    pub const fn metric(&self) -> &M {
        &self.metric
    }

    pub fn metric_mut(&mut self) -> &mut M {
        &mut self.metric
    }

    pub const fn step_size(&self) -> f64 {
        self.step_size
    }

    pub const fn max_tree_depth(&self) -> u16 {
        self.max_tree_depth
    }

    pub const fn max_energy_error(&self) -> f64 {
        self.max_energy_error
    }

    pub fn adaptation_is_frozen(&self) -> bool {
        self.warmup.as_ref().is_none_or(HmcWarmup::is_frozen)
    }

    fn ensure_workspace(&mut self, dimension: usize) -> Result<(), McmcError> {
        if self.metric.dimension() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.metric.dimension(),
                actual: dimension,
            });
        }
        if self.initial.position.is_empty()
            && self.initial.momentum.is_empty()
            && self.initial.gradient.is_empty()
        {
            self.initial = PhasePoint::with_dimension(dimension);
        }
        self.initial.validate_dimension(dimension)?;
        ensure_vector(&mut self.current_gradient, dimension)
    }

    fn prepare_sampling(&mut self) -> Result<(), McmcError> {
        if let Some(warmup) = &mut self.warmup {
            self.step_size = warmup.finish()?;
        }
        Ok(())
    }

    fn apply_metric_update(&mut self, update: MetricUpdate) -> Result<(), McmcError> {
        match update {
            MetricUpdate::Diagonal(diagonal) => self.metric.set_diagonal_inverse_mass(&diagonal),
            MetricUpdate::Dense {
                dimension,
                covariance,
                jitter,
            } => self
                .metric
                .set_dense_inverse_mass(dimension, &covariance, jitter),
        }
    }

    fn trajectory_continues(
        &self,
        left: &PhasePoint,
        right: &PhasePoint,
    ) -> Result<bool, McmcError> {
        no_u_turn(&self.metric, left, right)
    }
}

impl Nuts<UnitMetric> {
    pub fn unit(dimension: usize, step_size: f64, max_tree_depth: u16) -> Result<Self, McmcError> {
        Self::new(UnitMetric::new(dimension)?, step_size, max_tree_depth)
    }
}

impl Nuts<DiagonalMetric> {
    pub fn diagonal(
        inverse_mass: Vec<f64>,
        step_size: f64,
        max_tree_depth: u16,
    ) -> Result<Self, McmcError> {
        Self::new(
            DiagonalMetric::new(inverse_mass)?,
            step_size,
            max_tree_depth,
        )
    }

    pub fn with_diagonal_adaptation(
        self,
        total_warmup: u64,
        target_acceptance: f64,
        regularization: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::Diagonal { regularization },
            windows,
        )
    }
}

impl Nuts<DenseMetric> {
    pub fn dense(
        dimension: usize,
        inverse_mass: &[f64],
        jitter: f64,
        step_size: f64,
        max_tree_depth: u16,
    ) -> Result<Self, McmcError> {
        Self::new(
            DenseMetric::from_inverse_mass(dimension, inverse_mass, jitter)?,
            step_size,
            max_tree_depth,
        )
    }

    pub fn with_dense_adaptation(
        self,
        total_warmup: u64,
        target_acceptance: f64,
        regularization: f64,
    ) -> Result<Self, McmcError> {
        let windows = WarmupWindowConfig::default_for(total_warmup)?;
        self.with_warmup_adaptation(
            total_warmup,
            target_acceptance,
            MetricAdaptation::Dense { regularization },
            windows,
        )
    }
}

impl<T, M> TransitionKernel<T> for Nuts<M>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
{
    fn transition<R>(
        &mut self,
        target: &mut T,
        state: &mut EuclideanState,
        rng: &mut R,
        phase: SamplingPhase,
    ) -> Result<TransitionReport, McmcError>
    where
        R: Rng + ?Sized,
    {
        if phase == SamplingPhase::Sampling && !self.adaptation_is_frozen() {
            self.prepare_sampling()?;
        }
        state.validate()?;
        let dimension = state.dimension();
        self.ensure_workspace(dimension)?;
        let used_step_size = self.step_size;

        let mut target_evaluations = 0_u32;
        let mut gradient_evaluations = 0_u32;
        if let Some(gradient) = state.cache().gradient() {
            self.current_gradient.copy_from_slice(gradient);
        } else {
            let log_density =
                target.log_density_and_gradient(state.position(), &mut self.current_gradient);
            target_evaluations = target_evaluations.saturating_add(1);
            gradient_evaluations = gradient_evaluations.saturating_add(1);
            if !log_density.is_finite()
                || self.current_gradient.iter().any(|value| !value.is_finite())
            {
                return Err(McmcError::InvalidConfig(
                    "differentiable target returned a non-finite accepted-state gradient"
                        .to_string(),
                ));
            }
            state.synchronize_gradient(log_density, &self.current_gradient)?;
        }

        self.initial.position.copy_from_slice(state.position());
        self.initial
            .gradient
            .copy_from_slice(&self.current_gradient);
        self.initial.log_density = state.log_density();
        self.metric
            .sample_momentum(&mut self.initial.momentum, rng)?;
        let initial_energy =
            -state.log_density() + self.metric.kinetic_energy(&self.initial.momentum)?;
        if !initial_energy.is_finite() {
            return Err(McmcError::InvalidConfig(
                "NUTS initial Hamiltonian is non-finite".to_string(),
            ));
        }

        let mut left = self.initial.clone();
        let mut right = self.initial.clone();
        let mut proposal = self.initial.clone();
        let mut log_weight = 0.0;
        let mut acceptance_sum = 0.0;
        let mut acceptance_count = 0_u32;
        let mut leapfrog_steps = 0_u32;
        let mut divergent = false;
        let mut energy_error = None;
        let mut depth_reached = 0_u16;
        let mut stopped_early = false;

        for depth in 0..self.max_tree_depth {
            let direction = if rng.random::<bool>() { 1.0 } else { -1.0 };
            let start = if direction > 0.0 { &right } else { &left };
            let subtree = build_tree(
                target,
                &self.metric,
                &mut self.integrator,
                start,
                depth,
                direction,
                used_step_size,
                initial_energy,
                self.max_energy_error,
                rng,
            )?;

            target_evaluations = target_evaluations.saturating_add(subtree.target_evaluations);
            gradient_evaluations =
                gradient_evaluations.saturating_add(subtree.gradient_evaluations);
            leapfrog_steps = leapfrog_steps.saturating_add(subtree.leapfrog_steps);
            acceptance_sum += subtree.acceptance_sum;
            acceptance_count = acceptance_count.saturating_add(subtree.acceptance_count);
            divergent |= subtree.divergent;
            energy_error = worst_energy_error(energy_error, subtree.energy_error);

            // A subtree that already terminated internally contributes
            // diagnostics, but not a candidate or weight to the valid tree.
            if !subtree.continue_tree {
                stopped_early = true;
                break;
            }

            depth_reached = depth.saturating_add(1);
            if let Some(candidate) = subtree.proposal.as_ref() {
                let combined_weight = log_add_exp(log_weight, subtree.log_weight);
                if select_right(log_weight, subtree.log_weight, combined_weight, rng) {
                    proposal = candidate.clone();
                }
                log_weight = combined_weight;
            }

            if direction > 0.0 {
                right = subtree.right;
            } else {
                left = subtree.left;
            }

            if !self.trajectory_continues(&left, &right)? {
                stopped_early = true;
                break;
            }
        }

        let moved = proposal.position.as_slice() != state.position();
        if moved {
            state.commit_hamiltonian_proposal(
                &mut proposal.position,
                proposal.log_density,
                &proposal.gradient,
            )?;
        } else {
            state.mark_rejected_transition();
        }

        let acceptance_statistic = if acceptance_count == 0 {
            0.0
        } else {
            (acceptance_sum / f64::from(acceptance_count)).clamp(0.0, 1.0)
        };
        if phase == SamplingPhase::Warmup {
            if let Some(warmup) = &mut self.warmup {
                let observation = warmup.observe(acceptance_statistic, state.position())?;
                self.step_size = observation.step_size;
                if let Some(update) = observation.metric_update {
                    self.apply_metric_update(update)?;
                }
            }
        }

        Ok(TransitionReport {
            accepted: None,
            log_acceptance_probability: None,
            acceptance_statistic: Some(acceptance_statistic),
            proposals: 0,
            acceptances: 0,
            target_evaluations,
            gradient_evaluations,
            divergent,
            energy: Some(initial_energy),
            energy_error,
            leapfrog_steps,
            tree_depth: Some(depth_reached),
            max_tree_depth_reached: !stopped_early && depth_reached == self.max_tree_depth,
            proposal_scale: Some(used_step_size),
            subtransitions: 1,
        })
    }

    fn on_phase_start(
        &mut self,
        _target: &mut T,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Sampling {
            self.prepare_sampling()?;
        }
        Ok(())
    }

    fn on_phase_end(
        &mut self,
        _target: &mut T,
        phase: SamplingPhase,
        _state: &EuclideanState,
    ) -> Result<(), McmcError> {
        if phase == SamplingPhase::Warmup {
            self.prepare_sampling()?;
        }
        Ok(())
    }

    fn name(&self, _target: &T) -> &'static str {
        "Nuts"
    }
}

#[derive(Debug, Clone)]
struct Tree {
    left: PhasePoint,
    right: PhasePoint,
    proposal: Option<PhasePoint>,
    log_weight: f64,
    acceptance_sum: f64,
    acceptance_count: u32,
    target_evaluations: u32,
    gradient_evaluations: u32,
    leapfrog_steps: u32,
    divergent: bool,
    energy_error: Option<f64>,
    continue_tree: bool,
}

#[allow(clippy::too_many_arguments)]
fn build_tree<T, M, R>(
    target: &mut T,
    metric: &M,
    integrator: &mut LeapfrogIntegrator,
    start: &PhasePoint,
    depth: u16,
    direction: f64,
    step_size: f64,
    initial_energy: f64,
    max_energy_error: f64,
    rng: &mut R,
) -> Result<Tree, McmcError>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
    R: Rng + ?Sized,
{
    if depth == 0 {
        return build_leaf(
            target,
            metric,
            integrator,
            start,
            direction * step_size,
            initial_energy,
            max_energy_error,
        );
    }

    let first = build_tree(
        target,
        metric,
        integrator,
        start,
        depth - 1,
        direction,
        step_size,
        initial_energy,
        max_energy_error,
        rng,
    )?;
    if !first.continue_tree {
        return Ok(first);
    }

    let second_start = if direction > 0.0 {
        &first.right
    } else {
        &first.left
    };
    let second = build_tree(
        target,
        metric,
        integrator,
        second_start,
        depth - 1,
        direction,
        step_size,
        initial_energy,
        max_energy_error,
        rng,
    )?;

    let mut combined = if direction > 0.0 {
        combine_trees(first, second, rng)
    } else {
        combine_trees(second, first, rng)
    };
    if combined.continue_tree {
        combined.continue_tree = no_u_turn(metric, &combined.left, &combined.right)?;
    }
    Ok(combined)
}

fn build_leaf<T, M>(
    target: &mut T,
    metric: &M,
    integrator: &mut LeapfrogIntegrator,
    start: &PhasePoint,
    step_size: f64,
    initial_energy: f64,
    max_energy_error: f64,
) -> Result<Tree, McmcError>
where
    T: DifferentiableLogDensity + ?Sized,
    M: Metric,
{
    let mut point = start.clone();
    let integration = integrator.integrate(target, metric, &mut point, step_size, 1)?;
    if integration.invalid_trajectory {
        return Ok(Tree {
            left: point.clone(),
            right: point,
            proposal: None,
            log_weight: f64::NEG_INFINITY,
            acceptance_sum: 0.0,
            acceptance_count: 1,
            target_evaluations: integration.target_evaluations,
            gradient_evaluations: integration.gradient_evaluations,
            leapfrog_steps: 1,
            divergent: true,
            energy_error: None,
            continue_tree: false,
        });
    }

    let final_energy = -point.log_density + metric.kinetic_energy(&point.momentum)?;
    let difference = final_energy - initial_energy;
    let finite = difference.is_finite();
    let divergent = !finite || difference.abs() > max_energy_error;
    let acceptance_probability = if divergent {
        0.0
    } else {
        (-difference).min(0.0).exp()
    };
    let proposal = (!divergent).then(|| point.clone());
    Ok(Tree {
        left: point.clone(),
        right: point,
        proposal,
        log_weight: if divergent {
            f64::NEG_INFINITY
        } else {
            -difference
        },
        acceptance_sum: acceptance_probability,
        acceptance_count: 1,
        target_evaluations: integration.target_evaluations,
        gradient_evaluations: integration.gradient_evaluations,
        leapfrog_steps: 1,
        divergent,
        energy_error: finite.then_some(difference),
        continue_tree: !divergent,
    })
}

fn combine_trees<R>(left: Tree, right: Tree, rng: &mut R) -> Tree
where
    R: Rng + ?Sized,
{
    let log_weight = log_add_exp(left.log_weight, right.log_weight);
    let proposal = match (left.proposal, right.proposal) {
        (Some(left_proposal), Some(right_proposal)) => {
            if select_right(left.log_weight, right.log_weight, log_weight, rng) {
                Some(right_proposal)
            } else {
                Some(left_proposal)
            }
        }
        (Some(proposal), None) | (None, Some(proposal)) => Some(proposal),
        (None, None) => None,
    };
    Tree {
        left: left.left,
        right: right.right,
        proposal,
        log_weight,
        acceptance_sum: left.acceptance_sum + right.acceptance_sum,
        acceptance_count: left.acceptance_count.saturating_add(right.acceptance_count),
        target_evaluations: left
            .target_evaluations
            .saturating_add(right.target_evaluations),
        gradient_evaluations: left
            .gradient_evaluations
            .saturating_add(right.gradient_evaluations),
        leapfrog_steps: left.leapfrog_steps.saturating_add(right.leapfrog_steps),
        divergent: left.divergent || right.divergent,
        energy_error: worst_energy_error(left.energy_error, right.energy_error),
        continue_tree: left.continue_tree && right.continue_tree,
    }
}

fn no_u_turn<M>(metric: &M, left: &PhasePoint, right: &PhasePoint) -> Result<bool, McmcError>
where
    M: Metric,
{
    let left_dot =
        metric.displacement_dot_velocity(&left.position, &right.position, &left.momentum)?;
    let right_dot =
        metric.displacement_dot_velocity(&left.position, &right.position, &right.momentum)?;
    Ok(left_dot.is_finite() && right_dot.is_finite() && left_dot >= 0.0 && right_dot >= 0.0)
}

fn select_right<R>(left: f64, right: f64, total: f64, rng: &mut R) -> bool
where
    R: Rng + ?Sized,
{
    if !right.is_finite() {
        return false;
    }
    if !left.is_finite() {
        return true;
    }
    rng.random::<f64>() < (right - total).exp()
}

fn log_add_exp(left: f64, right: f64) -> f64 {
    if left == f64::NEG_INFINITY {
        return right;
    }
    if right == f64::NEG_INFINITY {
        return left;
    }
    let maximum = left.max(right);
    maximum + ((left - maximum).exp() + (right - maximum).exp()).ln()
}

fn worst_energy_error(left: Option<f64>, right: Option<f64>) -> Option<f64> {
    match (left, right) {
        (Some(left), Some(right)) if right.abs() > left.abs() => Some(right),
        (Some(left), Some(_)) => Some(left),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn ensure_vector(vector: &mut Vec<f64>, dimension: usize) -> Result<(), McmcError> {
    if vector.is_empty() {
        vector.resize(dimension, 0.0);
        Ok(())
    } else if vector.len() == dimension {
        Ok(())
    } else {
        Err(McmcError::DimensionMismatch {
            expected: dimension,
            actual: vector.len(),
        })
    }
}

fn validate_nuts_config(
    dimension: usize,
    step_size: f64,
    max_tree_depth: u16,
    max_energy_error: f64,
) -> Result<(), McmcError> {
    if dimension == 0 {
        return Err(McmcError::InvalidConfig(
            "NUTS dimension must be positive".to_string(),
        ));
    }
    if !step_size.is_finite() || step_size <= 0.0 {
        return Err(McmcError::InvalidConfig(
            "NUTS step size must be finite and positive".to_string(),
        ));
    }
    if max_tree_depth == 0 || max_tree_depth > MAX_SUPPORTED_TREE_DEPTH {
        return Err(McmcError::InvalidConfig(format!(
            "NUTS tree depth must lie between 1 and {MAX_SUPPORTED_TREE_DEPTH}"
        )));
    }
    if !max_energy_error.is_finite() || max_energy_error <= 0.0 {
        return Err(McmcError::InvalidConfig(
            "maximum NUTS energy error must be finite and positive".to_string(),
        ));
    }
    Ok(())
}

fn validate_metric_adaptation(
    metric: MetricKind,
    adaptation: MetricAdaptation,
) -> Result<(), McmcError> {
    let compatible = matches!(adaptation, MetricAdaptation::None)
        || matches!(
            (metric, adaptation),
            (MetricKind::Diagonal, MetricAdaptation::Diagonal { .. })
                | (MetricKind::Dense, MetricAdaptation::Diagonal { .. })
                | (MetricKind::Dense, MetricAdaptation::Dense { .. })
        );
    if compatible {
        Ok(())
    } else {
        Err(McmcError::InvalidConfig(
            "requested metric adaptation is incompatible with the NUTS metric".to_string(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u_turn_criterion_uses_endpoint_velocities() {
        let metric = UnitMetric::new(1).unwrap();
        let mut left = PhasePoint::with_dimension(1);
        let mut right = PhasePoint::with_dimension(1);
        left.position[0] = -1.0;
        right.position[0] = 1.0;
        left.momentum[0] = 1.0;
        right.momentum[0] = 1.0;
        assert!(no_u_turn(&metric, &left, &right).unwrap());
        right.momentum[0] = -1.0;
        assert!(!no_u_turn(&metric, &left, &right).unwrap());
    }

    #[test]
    fn log_weight_combination_is_stable() {
        let combined = log_add_exp(1_000.0, 999.0);
        assert!((combined - (1_000.0 + (-1.0_f64).exp().ln_1p())).abs() < 1.0e-12);
        assert_eq!(log_add_exp(f64::NEG_INFINITY, -3.0), -3.0);
    }
}
