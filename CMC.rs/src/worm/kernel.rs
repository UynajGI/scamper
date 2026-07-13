//! Detailed-balance persistent worm driver.

use super::{WormError, WormModel, WormSector, WormState, WormStepDelta, WormStepProposal};
use crate::audit::should_audit_cache;
use carlo_rs::accept_log_probability;
use rand::{Rng, RngExt};

/// Validated transition parameters for a persistent worm chain.
#[derive(Debug, Clone, PartialEq)]
pub struct WormConfig {
    /// Number of local extended-space transitions per Monte Carlo sweep.
    pub local_updates_per_sweep: usize,
    /// Probability of proposing closure when head and tail coincide.
    pub close_probability: f64,
    /// Logarithm of the relative worm-sector fugacity `η`.
    pub log_worm_fugacity: f64,
    /// Record a dense ordered endpoint-pair histogram when supported.
    pub track_endpoint_pairs: bool,
    /// Explicit cache-audit interval. Zero selects the build-mode policy.
    pub cache_audit_interval: u64,
}

impl Default for WormConfig {
    fn default() -> Self {
        Self {
            local_updates_per_sweep: 1,
            close_probability: 0.25,
            log_worm_fugacity: 0.0,
            track_endpoint_pairs: false,
            cache_audit_interval: 0,
        }
    }
}

impl WormConfig {
    pub fn validate(&self) -> Result<(), WormError> {
        if self.local_updates_per_sweep == 0 {
            return Err(WormError::new(
                "worm local_updates_per_sweep must be positive",
            ));
        }
        if !self.close_probability.is_finite()
            || self.close_probability <= 0.0
            || self.close_probability >= 1.0
        {
            return Err(WormError::new(
                "worm close_probability must lie strictly between zero and one",
            ));
        }
        if !self.log_worm_fugacity.is_finite() {
            return Err(WormError::new("worm log_worm_fugacity must be finite"));
        }
        Ok(())
    }
}

/// Transition category returned by one local update.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WormTransition {
    Open { accepted: bool },
    Close { accepted: bool },
    Step { accepted: bool },
    Bounce,
}

/// Cumulative and last-sweep transition diagnostics.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WormTransitionStatistics {
    pub open_attempts: u64,
    pub open_accepts: u64,
    pub close_attempts: u64,
    pub close_accepts: u64,
    pub step_attempts: u64,
    pub step_accepts: u64,
    pub bounces: u64,
    pub physical_visits: u64,
    pub worm_visits: u64,
    pub completed_worms: u64,
    pub total_completed_worm_steps: u64,
    pub last_completed_worm_steps: u64,
}

impl WormTransitionStatistics {
    pub fn step_acceptance_fraction(self) -> f64 {
        ratio(self.step_accepts, self.step_attempts)
    }

    pub fn open_acceptance_fraction(self) -> f64 {
        ratio(self.open_accepts, self.open_attempts)
    }

    pub fn close_acceptance_fraction(self) -> f64 {
        ratio(self.close_accepts, self.close_attempts)
    }

    pub fn mean_completed_worm_steps(self) -> f64 {
        ratio(self.total_completed_worm_steps, self.completed_worms)
    }
}

#[inline]
fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

/// Dense ordered endpoint-pair visit histogram.
///
/// With a constant worm fugacity and equal endpoint proposal measure,
/// `count(tail, head) / count(tail, tail)` estimates the corresponding
/// two-point correlation of the represented model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EndpointPairHistogram {
    bins: usize,
    counts: Vec<u64>,
    samples: u64,
}

impl EndpointPairHistogram {
    pub fn new(bins: usize) -> Result<Self, WormError> {
        let len = bins
            .checked_mul(bins)
            .ok_or_else(|| WormError::new("endpoint histogram size overflow"))?;
        Ok(Self {
            bins,
            counts: vec![0; len],
            samples: 0,
        })
    }

    #[inline]
    pub const fn bins(&self) -> usize {
        self.bins
    }

    #[inline]
    pub const fn samples(&self) -> u64 {
        self.samples
    }

    #[inline]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    pub fn observe(&mut self, tail: usize, head: usize) -> Result<(), WormError> {
        if tail >= self.bins || head >= self.bins {
            return Err(WormError::new("worm endpoint histogram index out of range"));
        }
        let index = tail * self.bins + head;
        self.counts[index] = self.counts[index].saturating_add(1);
        self.samples = self.samples.saturating_add(1);
        Ok(())
    }

    pub fn count(&self, tail: usize, head: usize) -> Option<u64> {
        if tail >= self.bins || head >= self.bins {
            None
        } else {
            Some(self.counts[tail * self.bins + head])
        }
    }

    pub fn correlation_ratio(&self, tail: usize, head: usize) -> Option<f64> {
        let diagonal = self.count(tail, tail)?;
        if diagonal == 0 {
            None
        } else {
            Some(self.count(tail, head)? as f64 / diagonal as f64)
        }
    }

    pub(crate) fn from_counts(
        bins: usize,
        counts: Vec<u64>,
        samples: u64,
    ) -> Result<Self, WormError> {
        if counts.len() != bins.saturating_mul(bins) {
            return Err(WormError::new(
                "endpoint histogram checkpoint length mismatch",
            ));
        }
        let counted_samples = counts.iter().try_fold(0u64, |total, count| {
            total
                .checked_add(*count)
                .ok_or_else(|| WormError::new("endpoint histogram checkpoint count overflow"))
        })?;
        if counted_samples != samples {
            return Err(WormError::new(
                "endpoint histogram checkpoint sample count mismatch",
            ));
        }
        Ok(Self {
            bins,
            counts,
            samples,
        })
    }
}

/// Generic persistent physical/worm-sector Markov kernel.
pub struct WormKernel<M>
where
    M: WormModel,
{
    model: M,
    state: WormState<M::Configuration, M::Defect>,
    config: WormConfig,
    patch: M::Patch,
    statistics: WormTransitionStatistics,
    last_sweep_statistics: WormTransitionStatistics,
    endpoint_pairs: Option<EndpointPairHistogram>,
    sweeps: u64,
    current_worm_steps: u64,
}

impl<M> WormKernel<M>
where
    M: WormModel,
{
    pub fn new(
        model: M,
        configuration: M::Configuration,
        config: WormConfig,
    ) -> Result<Self, WormError> {
        config.validate()?;
        let state = WormState::new(configuration);
        model.validate_state(&state)?;
        let endpoint_pairs = if config.track_endpoint_pairs {
            let bins = model.endpoint_bin_count();
            if bins == 0 {
                return Err(WormError::new(
                    "endpoint-pair tracking requested for a model without endpoint bins",
                ));
            }
            Some(EndpointPairHistogram::new(bins)?)
        } else {
            None
        };
        Ok(Self {
            model,
            state,
            config,
            patch: M::Patch::default(),
            statistics: WormTransitionStatistics::default(),
            last_sweep_statistics: WormTransitionStatistics::default(),
            endpoint_pairs,
            sweeps: 0,
            current_worm_steps: 0,
        })
    }

    #[inline]
    pub const fn model(&self) -> &M {
        &self.model
    }

    #[inline]
    pub const fn state(&self) -> &WormState<M::Configuration, M::Defect> {
        &self.state
    }

    #[inline]
    pub const fn config(&self) -> &WormConfig {
        &self.config
    }

    #[inline]
    pub const fn statistics(&self) -> &WormTransitionStatistics {
        &self.statistics
    }

    #[inline]
    pub const fn last_sweep_statistics(&self) -> &WormTransitionStatistics {
        &self.last_sweep_statistics
    }

    #[inline]
    pub const fn endpoint_pairs(&self) -> Option<&EndpointPairHistogram> {
        self.endpoint_pairs.as_ref()
    }

    #[inline]
    pub const fn sweeps(&self) -> u64 {
        self.sweeps
    }

    #[inline]
    pub const fn current_worm_steps(&self) -> u64 {
        self.current_worm_steps
    }

    pub fn validate(&self) -> Result<(), WormError> {
        self.config.validate()?;
        self.model.validate_state(&self.state)?;
        if self.model.open_defect_count(self.state.configuration()) == 0 {
            return Err(WormError::new(
                "worm model exposes no valid opening defects",
            ));
        }
        validate_statistics(
            self.statistics,
            self.state.sector(),
            self.current_worm_steps,
        )?;
        if let Some(histogram) = &self.endpoint_pairs {
            if histogram.bins() != self.model.endpoint_bin_count() {
                return Err(WormError::new(
                    "endpoint histogram bin count disagrees with the model",
                ));
            }
            if histogram.samples() != self.statistics.worm_visits {
                return Err(WormError::new(
                    "endpoint histogram samples disagree with worm-sector visits",
                ));
            }
        }
        Ok(())
    }

    pub fn local_update(&mut self, rng: &mut impl Rng) -> Result<WormTransition, WormError> {
        let before = self.statistics;
        let transition = match self.state.sector() {
            WormSector::Physical => self.try_open(rng)?,
            WormSector::Worm => {
                let coincident = self.state.head() == self.state.tail();
                if coincident && rng.random::<f64>() < self.config.close_probability {
                    self.try_close(rng)?
                } else {
                    self.try_step(rng)?
                }
            }
        };
        self.record_sector_visit()?;
        self.last_sweep_statistics = difference(self.statistics, before);
        Ok(transition)
    }

    pub fn sweep(&mut self, rng: &mut impl Rng) -> Result<(), WormError> {
        let before = self.statistics;
        for _ in 0..self.config.local_updates_per_sweep {
            self.local_update(rng)?;
        }
        self.last_sweep_statistics = difference(self.statistics, before);
        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.config.cache_audit_interval) {
            self.validate()?;
        }
        Ok(())
    }

    fn try_open(&mut self, rng: &mut impl Rng) -> Result<WormTransition, WormError> {
        self.statistics.open_attempts = self.statistics.open_attempts.saturating_add(1);
        let count = self.model.open_defect_count(self.state.configuration());
        if count == 0 {
            return Err(WormError::new(
                "worm model exposes no valid opening defects",
            ));
        }
        let defect = self
            .model
            .open_defect(self.state.configuration(), rng.random_range(0..count))?;
        let log_acceptance = self.config.log_worm_fugacity
            + self.config.close_probability.ln()
            + (count as f64).ln();
        let accepted = accept_log_probability(log_acceptance, rng);
        if accepted {
            self.state.open(defect)?;
            self.current_worm_steps = 0;
            self.statistics.open_accepts = self.statistics.open_accepts.saturating_add(1);
        }
        Ok(WormTransition::Open { accepted })
    }

    fn try_close(&mut self, rng: &mut impl Rng) -> Result<WormTransition, WormError> {
        self.statistics.close_attempts = self.statistics.close_attempts.saturating_add(1);
        let count = self.model.open_defect_count(self.state.configuration());
        if count == 0 {
            return Err(WormError::new(
                "worm model exposes no valid reverse opening defects",
            ));
        }
        let log_acceptance = -self.config.log_worm_fugacity
            - self.config.close_probability.ln()
            - (count as f64).ln();
        let accepted = accept_log_probability(log_acceptance, rng);
        if accepted {
            self.state.close()?;
            self.statistics.close_accepts = self.statistics.close_accepts.saturating_add(1);
            self.statistics.completed_worms = self.statistics.completed_worms.saturating_add(1);
            self.statistics.last_completed_worm_steps = self.current_worm_steps;
            self.statistics.total_completed_worm_steps = self
                .statistics
                .total_completed_worm_steps
                .saturating_add(self.current_worm_steps);
            self.current_worm_steps = 0;
        }
        Ok(WormTransition::Close { accepted })
    }

    fn try_step(&mut self, rng: &mut impl Rng) -> Result<WormTransition, WormError> {
        let Some(WormStepProposal {
            step,
            log_reverse_over_forward,
        }) = self.model.propose_step(&self.state, rng)?
        else {
            self.statistics.bounces = self.statistics.bounces.saturating_add(1);
            return Ok(WormTransition::Bounce);
        };

        self.statistics.step_attempts = self.statistics.step_attempts.saturating_add(1);
        let old_coincident = self.state.head() == self.state.tail();
        let WormStepDelta {
            new_head,
            log_weight_ratio,
        } = self
            .model
            .evaluate_step(&self.state, &step, &mut self.patch)?;
        let new_coincident = self.state.tail().is_some_and(|tail| tail == &new_head);
        let old_step_probability = if old_coincident {
            1.0 - self.config.close_probability
        } else {
            1.0
        };
        let new_step_probability = if new_coincident {
            1.0 - self.config.close_probability
        } else {
            1.0
        };
        let log_acceptance =
            log_weight_ratio + log_reverse_over_forward + new_step_probability.ln()
                - old_step_probability.ln();
        let accepted = accept_log_probability(log_acceptance, rng);
        if accepted {
            self.model.commit_step(&mut self.state, &step, &self.patch);
            self.state.move_head(new_head)?;
            self.current_worm_steps = self.current_worm_steps.saturating_add(1);
            self.statistics.step_accepts = self.statistics.step_accepts.saturating_add(1);
        }
        Ok(WormTransition::Step { accepted })
    }

    fn record_sector_visit(&mut self) -> Result<(), WormError> {
        match self.state.sector() {
            WormSector::Physical => {
                self.statistics.physical_visits = self.statistics.physical_visits.saturating_add(1);
            }
            WormSector::Worm => {
                self.statistics.worm_visits = self.statistics.worm_visits.saturating_add(1);
                if let Some(histogram) = &mut self.endpoint_pairs {
                    let head = self
                        .state
                        .head()
                        .and_then(|defect| self.model.endpoint_bin(defect))
                        .ok_or_else(|| WormError::new("worm head has no endpoint bin"))?;
                    let tail = self
                        .state
                        .tail()
                        .and_then(|defect| self.model.endpoint_bin(defect))
                        .ok_or_else(|| WormError::new("worm tail has no endpoint bin"))?;
                    histogram.observe(tail, head)?;
                }
            }
        }
        Ok(())
    }

    pub(crate) fn restore_runtime(
        &mut self,
        state: WormState<M::Configuration, M::Defect>,
        statistics: WormTransitionStatistics,
        endpoint_pairs: Option<EndpointPairHistogram>,
        sweeps: u64,
        current_worm_steps: u64,
    ) -> Result<(), WormError> {
        self.model.validate_state(&state)?;
        match (&self.endpoint_pairs, &endpoint_pairs) {
            (Some(expected), Some(restored)) if expected.bins() == restored.bins() => {}
            (None, None) => {}
            _ => {
                return Err(WormError::new(
                    "endpoint histogram checkpoint configuration mismatch",
                ));
            }
        }
        self.state = state;
        self.statistics = statistics;
        self.last_sweep_statistics = WormTransitionStatistics::default();
        self.endpoint_pairs = endpoint_pairs;
        self.sweeps = sweeps;
        self.current_worm_steps = current_worm_steps;
        self.validate()
    }
}

fn validate_statistics(
    statistics: WormTransitionStatistics,
    sector: WormSector,
    current_worm_steps: u64,
) -> Result<(), WormError> {
    if statistics.open_accepts > statistics.open_attempts
        || statistics.close_accepts > statistics.close_attempts
        || statistics.step_accepts > statistics.step_attempts
    {
        return Err(WormError::new(
            "worm accepted-transition count exceeds attempted count",
        ));
    }
    let transition_count = [
        statistics.open_attempts,
        statistics.close_attempts,
        statistics.step_attempts,
        statistics.bounces,
    ]
    .into_iter()
    .try_fold(0u64, |total, count| total.checked_add(count))
    .ok_or_else(|| WormError::new("worm transition counter overflow"))?;
    let visit_count = statistics
        .physical_visits
        .checked_add(statistics.worm_visits)
        .ok_or_else(|| WormError::new("worm sector-visit counter overflow"))?;
    if transition_count != visit_count {
        return Err(WormError::new(
            "worm transition counters disagree with sector visits",
        ));
    }
    if statistics.completed_worms != statistics.close_accepts {
        return Err(WormError::new(
            "worm completion count disagrees with accepted closures",
        ));
    }
    if statistics.total_completed_worm_steps < statistics.last_completed_worm_steps {
        return Err(WormError::new(
            "worm completed-length counters are inconsistent",
        ));
    }
    let open_worms = statistics
        .open_accepts
        .checked_sub(statistics.close_accepts)
        .ok_or_else(|| WormError::new("worm closures exceed accepted openings"))?;
    let expected_open_worms = u64::from(matches!(sector, WormSector::Worm));
    if open_worms != expected_open_worms {
        return Err(WormError::new(
            "worm sector disagrees with accepted open/close counts",
        ));
    }
    if matches!(sector, WormSector::Physical) && current_worm_steps != 0 {
        return Err(WormError::new(
            "physical sector cannot retain an in-progress worm length",
        ));
    }
    Ok(())
}

fn difference(
    after: WormTransitionStatistics,
    before: WormTransitionStatistics,
) -> WormTransitionStatistics {
    WormTransitionStatistics {
        open_attempts: after.open_attempts.saturating_sub(before.open_attempts),
        open_accepts: after.open_accepts.saturating_sub(before.open_accepts),
        close_attempts: after.close_attempts.saturating_sub(before.close_attempts),
        close_accepts: after.close_accepts.saturating_sub(before.close_accepts),
        step_attempts: after.step_attempts.saturating_sub(before.step_attempts),
        step_accepts: after.step_accepts.saturating_sub(before.step_accepts),
        bounces: after.bounces.saturating_sub(before.bounces),
        physical_visits: after.physical_visits.saturating_sub(before.physical_visits),
        worm_visits: after.worm_visits.saturating_sub(before.worm_visits),
        completed_worms: after.completed_worms.saturating_sub(before.completed_worms),
        total_completed_worm_steps: after
            .total_completed_worm_steps
            .saturating_sub(before.total_completed_worm_steps),
        last_completed_worm_steps: after.last_completed_worm_steps,
    }
}
