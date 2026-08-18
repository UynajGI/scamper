//! Adaptive Wang–Landau density-of-states estimation.
//!
//! The estimator owns the adaptive lifecycle and checkpoint state, while the
//! update kernel owns local proposals and transactional state changes.  Once
//! the estimator freezes, the same walk becomes a fixed multicanonical
//! production run with weight `1 / g(E)`.

use crate::algorithms::{Algorithm, SimulationPhase};
use crate::audit::{audit_lattice_cache, audit_macrostate_bin, should_audit_cache};
use crate::classical_mc::{
    build_lattice_from_params, parse_bool, parse_param, ClassicalMC, FromHamiltonianParams,
};
use crate::core::cache::EnergyPatch;
use crate::core::r#move::SiteSpinMove;
use crate::core::trial::TrialEvaluator;
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::generalized::{GeneralizedError, Histogram, LogDensityOfStates, MacrostateAxis};
use crate::lattice::interaction::{Hamiltonian, Initializable, Proposable};
use crate::lattice::models::IsingModel;
use crate::lattice::proposal::{ProposalStrategy, StandardStrategy};
use crate::lattice::state::System;
use crate::observables::{DefaultObservableSet, ObservableSet};
use carlo_rs::{
    accept_log_probability, AdaptiveRunControl, CarloError, Context, FromParams, MonteCarlo,
    Params, RunDecision, RunPhase,
};
use rand::Rng;
use serde_json::{json, Value as Json};

/// Lifecycle of one Wang–Landau estimate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WangLandauPhase {
    /// Initial range exploration. The DOS is updated, but flatness is not tested.
    Discovery,
    /// Adaptive flat-histogram or `1/t` refinement.
    Adaptation,
    /// Frozen `1/g(E)` production sampling.
    FrozenProduction,
    /// Terminal state after convergence or a maximum-sweep guard.
    Finished,
}

impl WangLandauPhase {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Discovery => "discovery",
            Self::Adaptation => "adaptation",
            Self::FrozenProduction => "frozen_production",
            Self::Finished => "finished",
        }
    }

    fn parse(value: &str) -> Result<Self, GeneralizedError> {
        match value {
            "discovery" => Ok(Self::Discovery),
            "adaptation" => Ok(Self::Adaptation),
            "frozen_production" => Ok(Self::FrozenProduction),
            "finished" => Ok(Self::Finished),
            _ => Err(GeneralizedError::new(format!(
                "unknown Wang-Landau phase `{value}`"
            ))),
        }
    }
}

/// Active modification-factor schedule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WangLandauRefinement {
    FlatHistogram,
    OneOverT,
}

impl WangLandauRefinement {
    const fn as_str(self) -> &'static str {
        match self {
            Self::FlatHistogram => "flat_histogram",
            Self::OneOverT => "one_over_t",
        }
    }

    fn parse(value: &str) -> Result<Self, GeneralizedError> {
        match value {
            "flat_histogram" => Ok(Self::FlatHistogram),
            "one_over_t" => Ok(Self::OneOverT),
            _ => Err(GeneralizedError::new(format!(
                "unknown Wang-Landau refinement `{value}`"
            ))),
        }
    }
}

/// Why an adaptive estimate stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WangLandauTermination {
    Converged,
    MaximumSweeps,
    FrozenByDriver,
    /// The configured [`WangLandauConfig::minimum_visited_fraction`] demands
    /// more visited bins than the walk can ever reach: bin discovery has
    /// plateaued below the required count, so the flatness gate can never
    /// pass. Terminates loudly instead of running to the sweep guard with a
    /// silently unconverged density of states.
    UnreachableBins,
}

impl WangLandauTermination {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Converged => "converged",
            Self::MaximumSweeps => "maximum_sweeps",
            Self::FrozenByDriver => "frozen_by_driver",
            Self::UnreachableBins => "unreachable_bins",
        }
    }

    fn parse(value: &str) -> Result<Self, GeneralizedError> {
        match value {
            "converged" => Ok(Self::Converged),
            "maximum_sweeps" => Ok(Self::MaximumSweeps),
            "frozen_by_driver" => Ok(Self::FrozenByDriver),
            "unreachable_bins" => Ok(Self::UnreachableBins),
            _ => Err(GeneralizedError::new(format!(
                "unknown Wang-Landau termination `{value}`"
            ))),
        }
    }
}

/// Validated Wang–Landau refinement parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct WangLandauConfig {
    /// Initial additive modification factor `ln f`.
    pub initial_log_f: f64,
    /// Freeze once the current modification factor is at or below this value.
    pub final_log_f: f64,
    /// Minimum count divided by the mean visited-bin count.
    pub flatness: f64,
    /// Number of completed adaptive sweeps between flatness checks.
    pub flatness_check_interval: u64,
    /// Number of initial sweeps before flatness checks begin.
    pub discovery_sweeps: u64,
    /// Switch from geometric flat-histogram refinement at this `ln f`.
    /// Set to zero to disable `1/t` refinement.
    pub one_over_t_threshold: f64,
    /// Maximum discovery/adaptation sweeps. Zero disables the guard.
    pub max_adaptation_sweeps: u64,
    /// Required fraction of represented bins with at least one visit.
    ///
    /// The flatness gate cannot pass until the walk has discovered at least
    /// `ceil(fraction · bins)` bins. Bins that are physically unreachable on
    /// the axis (for example, energy ranges with zero density of states on a
    /// user-supplied [`BinnedAxis`](crate::BinnedAxis)) are never discovered,
    /// so a fraction above the reachable fraction makes convergence
    /// impossible. The estimator auto-derives the reachable set as the
    /// discovery plateau of the walk: if the required fraction stays
    /// unattainable over many consecutive flatness checks, the estimate
    /// terminates loudly with
    /// [`WangLandauTermination::UnreachableBins`] instead of silently
    /// running to the maximum-sweep guard with an unconverged density of
    /// states. Exact-enumeration axes (every bin occupied) keep the strict
    /// default of `1.0`.
    pub minimum_visited_fraction: f64,
}

impl Default for WangLandauConfig {
    fn default() -> Self {
        Self {
            initial_log_f: 1.0,
            final_log_f: 1e-8,
            flatness: 0.8,
            flatness_check_interval: 100,
            discovery_sweeps: 0,
            one_over_t_threshold: 1e-4,
            max_adaptation_sweeps: 10_000_000,
            minimum_visited_fraction: 1.0,
        }
    }
}

impl WangLandauConfig {
    pub fn validate(&self) -> Result<(), GeneralizedError> {
        if !self.initial_log_f.is_finite() || self.initial_log_f <= 0.0 {
            return Err(GeneralizedError::new(
                "initial Wang-Landau log_f must be finite and positive",
            ));
        }
        if !self.final_log_f.is_finite()
            || self.final_log_f <= 0.0
            || self.final_log_f >= self.initial_log_f
        {
            return Err(GeneralizedError::new(concat!(
                "final Wang-Landau log_f must be finite, positive and smaller than ",
                "the initial value",
            )));
        }
        if !self.flatness.is_finite() || !(0.0..=1.0).contains(&self.flatness) {
            return Err(GeneralizedError::new(
                "Wang-Landau flatness must lie in [0, 1]",
            ));
        }
        if self.flatness_check_interval == 0 {
            return Err(GeneralizedError::new(
                "Wang-Landau flatness check interval must be positive",
            ));
        }
        if !self.one_over_t_threshold.is_finite()
            || self.one_over_t_threshold < 0.0
            || self.one_over_t_threshold >= self.initial_log_f
        {
            return Err(GeneralizedError::new(
                "Wang-Landau 1/t threshold must be finite, non-negative and below initial log_f",
            ));
        }
        if !self.minimum_visited_fraction.is_finite()
            || !(0.0..=1.0).contains(&self.minimum_visited_fraction)
        {
            return Err(GeneralizedError::new(
                "minimum visited fraction must lie in [0, 1]",
            ));
        }
        Ok(())
    }
}

/// Number of consecutive flatness checks without a newly discovered bin that
/// establishes the discovery plateau. Once the plateau is established below
/// the required visited fraction, no amount of further sweeping can satisfy
/// the flatness gate, so the estimate terminates loudly. With the default
/// check interval of 100 sweeps this waits 50 000 adaptation sweeps before
/// declaring the fraction unattainable — far beyond the discovery time of
/// any reachable bin on axes of practical size, while still orders of
/// magnitude below the default 10-million-sweep guard.
const DISCOVERY_STALL_CHECK_LIMIT: u32 = 500;

/// Complete adaptive estimator and its production histogram.
#[derive(Debug, Clone, PartialEq)]
pub struct WangLandauState {
    config: WangLandauConfig,
    log_density: LogDensityOfStates,
    adaptation_histogram: Histogram,
    production_histogram: Histogram,
    phase: WangLandauPhase,
    refinement: WangLandauRefinement,
    termination: Option<WangLandauTermination>,
    log_f: f64,
    adaptation_sweeps: u64,
    production_sweeps: u64,
    refinement_visits: u64,
    flatness_checks: u64,
    flatness_passes: u64,
    /// Visited-bin count at the previous flatness check (discovery tracking).
    last_check_visited_bins: usize,
    /// Consecutive flatness checks without a newly discovered bin.
    discovery_stall_checks: u32,
}

impl WangLandauState {
    pub fn new(bins: usize, config: WangLandauConfig) -> Result<Self, GeneralizedError> {
        config.validate()?;
        Ok(Self {
            log_density: LogDensityOfStates::new(bins)?,
            adaptation_histogram: Histogram::new(bins)?,
            production_histogram: Histogram::new(bins)?,
            phase: if config.discovery_sweeps == 0 {
                WangLandauPhase::Adaptation
            } else {
                WangLandauPhase::Discovery
            },
            refinement: WangLandauRefinement::FlatHistogram,
            termination: None,
            log_f: config.initial_log_f,
            adaptation_sweeps: 0,
            production_sweeps: 0,
            refinement_visits: 0,
            flatness_checks: 0,
            flatness_passes: 0,
            last_check_visited_bins: 0,
            discovery_stall_checks: 0,
            config,
        })
    }

    #[inline]
    pub const fn config(&self) -> &WangLandauConfig {
        &self.config
    }

    #[inline]
    pub const fn log_density(&self) -> &LogDensityOfStates {
        &self.log_density
    }

    #[inline]
    pub const fn adaptation_histogram(&self) -> &Histogram {
        &self.adaptation_histogram
    }

    #[inline]
    pub const fn production_histogram(&self) -> &Histogram {
        &self.production_histogram
    }

    #[inline]
    pub const fn phase(&self) -> WangLandauPhase {
        self.phase
    }

    #[inline]
    pub const fn refinement(&self) -> WangLandauRefinement {
        self.refinement
    }

    #[inline]
    pub const fn termination(&self) -> Option<WangLandauTermination> {
        self.termination
    }

    #[inline]
    pub const fn log_f(&self) -> f64 {
        self.log_f
    }

    #[inline]
    pub const fn adaptation_sweeps(&self) -> u64 {
        self.adaptation_sweeps
    }

    #[inline]
    pub const fn production_sweeps(&self) -> u64 {
        self.production_sweeps
    }

    #[inline]
    pub const fn refinement_visits(&self) -> u64 {
        self.refinement_visits
    }

    #[inline]
    pub const fn flatness_checks(&self) -> u64 {
        self.flatness_checks
    }

    #[inline]
    pub const fn flatness_passes(&self) -> u64 {
        self.flatness_passes
    }

    #[inline]
    pub const fn is_adaptive(&self) -> bool {
        matches!(
            self.phase,
            WangLandauPhase::Discovery | WangLandauPhase::Adaptation
        )
    }

    #[inline]
    pub const fn is_frozen(&self) -> bool {
        matches!(self.phase, WangLandauPhase::FrozenProduction)
    }

    /// `ln[g(E_old)/g(E_new)]`, the log acceptance contribution of the walk.
    pub fn log_weight_ratio(&self, old_bin: usize, new_bin: usize) -> f64 {
        self.log_density.value(old_bin) - self.log_density.value(new_bin)
    }

    /// Record the state occupied after one attempted transition.
    pub fn record_visit(&mut self, bin: usize) {
        match self.phase {
            WangLandauPhase::Discovery | WangLandauPhase::Adaptation => {
                self.adaptation_histogram.record(bin);
                self.refinement_visits = self.refinement_visits.saturating_add(1);
                let increment = match self.refinement {
                    WangLandauRefinement::FlatHistogram => self.log_f,
                    WangLandauRefinement::OneOverT => self.one_over_t_increment().min(self.log_f),
                };
                self.log_f = increment;
                self.log_density.add_visit_weight(bin, increment);
            }
            WangLandauPhase::FrozenProduction => self.production_histogram.record(bin),
            WangLandauPhase::Finished => {}
        }
    }

    /// Finish one sweep and advance the adaptive lifecycle if needed.
    pub fn finish_sweep(&mut self) {
        match self.phase {
            WangLandauPhase::Discovery | WangLandauPhase::Adaptation => {
                self.adaptation_sweeps = self.adaptation_sweeps.saturating_add(1);
                if self.phase == WangLandauPhase::Discovery
                    && self.adaptation_sweeps >= self.config.discovery_sweeps
                {
                    self.phase = WangLandauPhase::Adaptation;
                    self.adaptation_histogram.clear();
                    self.refinement_visits = 0;
                }
                if self.phase == WangLandauPhase::Adaptation {
                    match self.refinement {
                        WangLandauRefinement::FlatHistogram => self.finish_flat_histogram_sweep(),
                        WangLandauRefinement::OneOverT => self.finish_one_over_t_sweep(),
                    }
                }
                if self.config.max_adaptation_sweeps > 0
                    && self.adaptation_sweeps >= self.config.max_adaptation_sweeps
                    && self.is_adaptive()
                {
                    self.log_density.normalize_max_zero();
                    self.phase = WangLandauPhase::Finished;
                    self.termination = Some(WangLandauTermination::MaximumSweeps);
                }
            }
            WangLandauPhase::FrozenProduction => {
                self.production_sweeps = self.production_sweeps.saturating_add(1);
            }
            WangLandauPhase::Finished => {}
        }
    }

    /// Freeze the current estimate and begin fixed-weight production.
    pub fn freeze_for_production(&mut self) {
        if self.is_adaptive() {
            self.log_density.normalize_max_zero();
            self.phase = WangLandauPhase::FrozenProduction;
            if self.termination.is_none() {
                self.termination = Some(WangLandauTermination::FrozenByDriver);
            }
            self.production_histogram.clear();
        }
    }

    /// Mark a completed production run terminal without changing its DOS.
    pub fn finish_production(&mut self) {
        if self.is_frozen() {
            self.phase = WangLandauPhase::Finished;
        }
    }

    /// Versioned JSON checkpoint of every adaptive field.
    pub fn save_snapshot(&self) -> Json {
        json!({
            "format": "cmc-rs-wang-landau-v1",
            "config": {
                "initial_log_f": self.config.initial_log_f,
                "final_log_f": self.config.final_log_f,
                "flatness": self.config.flatness,
                "flatness_check_interval": self.config.flatness_check_interval,
                "discovery_sweeps": self.config.discovery_sweeps,
                "one_over_t_threshold": self.config.one_over_t_threshold,
                "max_adaptation_sweeps": self.config.max_adaptation_sweeps,
                "minimum_visited_fraction": self.config.minimum_visited_fraction,
            },
            "log_density": self.log_density.values(),
            "visited": self.log_density.visited(),
            "adaptation_histogram": self.adaptation_histogram.counts(),
            "production_histogram": self.production_histogram.counts(),
            "phase": self.phase.as_str(),
            "refinement": self.refinement.as_str(),
            "termination": self.termination.map(WangLandauTermination::as_str),
            "log_f": self.log_f,
            "adaptation_sweeps": self.adaptation_sweeps,
            "production_sweeps": self.production_sweeps,
            "refinement_visits": self.refinement_visits,
            "flatness_checks": self.flatness_checks,
            "flatness_passes": self.flatness_passes,
            "last_check_visited_bins": self.last_check_visited_bins,
            "discovery_stall_checks": self.discovery_stall_checks,
        })
    }

    /// Restore and validate a version-1 JSON checkpoint.
    pub fn load_snapshot(snapshot: &Json) -> Result<Self, GeneralizedError> {
        if snapshot["format"].as_str() != Some("cmc-rs-wang-landau-v1") {
            return Err(GeneralizedError::new(
                "unknown or missing Wang-Landau checkpoint format",
            ));
        }
        let config_json = &snapshot["config"];
        let config = WangLandauConfig {
            initial_log_f: required_f64(config_json, "initial_log_f")?,
            final_log_f: required_f64(config_json, "final_log_f")?,
            flatness: required_f64(config_json, "flatness")?,
            flatness_check_interval: required_u64(config_json, "flatness_check_interval")?,
            discovery_sweeps: required_u64(config_json, "discovery_sweeps")?,
            one_over_t_threshold: required_f64(config_json, "one_over_t_threshold")?,
            max_adaptation_sweeps: required_u64(config_json, "max_adaptation_sweeps")?,
            minimum_visited_fraction: required_f64(config_json, "minimum_visited_fraction")?,
        };
        config.validate()?;

        let values = required_f64_array(snapshot, "log_density")?;
        let visited = required_bool_array(snapshot, "visited")?;
        let adaptation_counts = required_u64_array(snapshot, "adaptation_histogram")?;
        let production_counts = required_u64_array(snapshot, "production_histogram")?;
        let bins = values.len();
        if bins == 0
            || visited.len() != bins
            || adaptation_counts.len() != bins
            || production_counts.len() != bins
        {
            return Err(GeneralizedError::new(
                "Wang-Landau checkpoint buffers have inconsistent bin counts",
            ));
        }

        let phase = WangLandauPhase::parse(required_str(snapshot, "phase")?)?;
        let refinement = WangLandauRefinement::parse(required_str(snapshot, "refinement")?)?;
        let termination_value = snapshot
            .get("termination")
            .ok_or_else(|| GeneralizedError::new("missing Wang-Landau checkpoint termination"))?;
        let termination = match termination_value.as_str() {
            Some(value) => Some(WangLandauTermination::parse(value)?),
            None if termination_value.is_null() => None,
            None => {
                return Err(GeneralizedError::new(
                    "invalid Wang-Landau checkpoint termination",
                ));
            }
        };
        let log_f = required_f64(snapshot, "log_f")?;
        if log_f <= 0.0 {
            return Err(GeneralizedError::new(
                "checkpoint Wang-Landau log_f must be positive",
            ));
        }

        let state = Self {
            config,
            log_density: LogDensityOfStates::from_values(values, visited)?,
            adaptation_histogram: Histogram::from_counts(adaptation_counts)?,
            production_histogram: Histogram::from_counts(production_counts)?,
            phase,
            refinement,
            termination,
            log_f,
            adaptation_sweeps: required_u64(snapshot, "adaptation_sweeps")?,
            production_sweeps: required_u64(snapshot, "production_sweeps")?,
            refinement_visits: required_u64(snapshot, "refinement_visits")?,
            flatness_checks: required_u64(snapshot, "flatness_checks")?,
            flatness_passes: required_u64(snapshot, "flatness_passes")?,
            // Version-1 checkpoints predate discovery-stall tracking; absent
            // fields restart the plateau detection from zero.
            last_check_visited_bins: snapshot["last_check_visited_bins"]
                .as_u64()
                .map(|value| value as usize)
                .unwrap_or(0),
            discovery_stall_checks: snapshot["discovery_stall_checks"]
                .as_u64()
                .map(|value| value as u32)
                .unwrap_or(0),
        };
        state.validate_checkpoint_consistency()?;
        Ok(state)
    }

    pub fn validate_bin_count(&self, expected: usize) -> Result<(), GeneralizedError> {
        if self.log_density.bins() != expected
            || self.adaptation_histogram.bins() != expected
            || self.production_histogram.bins() != expected
        {
            return Err(GeneralizedError::new(
                "Wang-Landau DOS/histogram buffers disagree with the macrostate axis",
            ));
        }
        if self.refinement_visits != self.adaptation_histogram.total() {
            return Err(GeneralizedError::new(
                "Wang-Landau refinement visits disagree with the adaptive histogram total",
            ));
        }
        Ok(())
    }

    fn validate_checkpoint_consistency(&self) -> Result<(), GeneralizedError> {
        self.validate_bin_count(self.log_density.bins())?;
        if self.refinement_visits != self.adaptation_histogram.total() {
            return Err(GeneralizedError::new(
                "checkpoint refinement visits do not match the adaptive histogram total",
            ));
        }
        if self.flatness_passes > self.flatness_checks {
            return Err(GeneralizedError::new(
                "checkpoint flatness passes exceed flatness checks",
            ));
        }
        if self.phase == WangLandauPhase::Discovery
            && self.adaptation_sweeps >= self.config.discovery_sweeps
        {
            return Err(GeneralizedError::new(
                "discovery checkpoint has already reached its configured discovery length",
            ));
        }
        if self.phase == WangLandauPhase::Adaptation
            && self.adaptation_sweeps < self.config.discovery_sweeps
        {
            return Err(GeneralizedError::new(
                "adaptation checkpoint has not completed discovery",
            ));
        }
        if self.phase == WangLandauPhase::Discovery
            && self.refinement == WangLandauRefinement::OneOverT
        {
            return Err(GeneralizedError::new(
                "discovery checkpoint cannot use 1/t refinement",
            ));
        }
        if self.log_f > self.config.initial_log_f {
            return Err(GeneralizedError::new(
                "checkpoint Wang-Landau log_f exceeds its configured initial value",
            ));
        }
        if self.refinement == WangLandauRefinement::OneOverT
            && self.config.one_over_t_threshold == 0.0
        {
            return Err(GeneralizedError::new(
                "checkpoint enables 1/t refinement although the schedule disables it",
            ));
        }
        match self.phase {
            WangLandauPhase::Discovery | WangLandauPhase::Adaptation => {
                if self.termination.is_some() {
                    return Err(GeneralizedError::new(
                        "adaptive checkpoint phase cannot already have a termination reason",
                    ));
                }
                if self.production_sweeps > 0 || self.production_histogram.total() > 0 {
                    return Err(GeneralizedError::new(
                        "adaptive checkpoint cannot contain production progress",
                    ));
                }
            }
            WangLandauPhase::FrozenProduction => {
                if !matches!(
                    self.termination,
                    Some(WangLandauTermination::Converged | WangLandauTermination::FrozenByDriver)
                ) {
                    return Err(GeneralizedError::new(
                        "frozen-production checkpoint lacks a compatible termination reason",
                    ));
                }
            }
            WangLandauPhase::Finished => {
                if self.termination.is_none() {
                    return Err(GeneralizedError::new(
                        "finished checkpoint lacks a termination reason",
                    ));
                }
            }
        }
        if self.termination == Some(WangLandauTermination::MaximumSweeps) {
            if self.phase != WangLandauPhase::Finished {
                return Err(GeneralizedError::new(
                    "maximum-sweeps termination requires the finished phase",
                ));
            }
            if self.config.max_adaptation_sweeps == 0
                || self.adaptation_sweeps < self.config.max_adaptation_sweeps
            {
                return Err(GeneralizedError::new(
                    "maximum-sweeps termination is inconsistent with the configured guard",
                ));
            }
            if self.production_sweeps > 0 || self.production_histogram.total() > 0 {
                return Err(GeneralizedError::new(
                    "maximum-sweeps termination cannot contain production progress",
                ));
            }
        }
        if self.termination == Some(WangLandauTermination::UnreachableBins) {
            if self.phase != WangLandauPhase::Finished {
                return Err(GeneralizedError::new(
                    "unreachable-bins termination requires the finished phase",
                ));
            }
            if self.discovery_stall_checks < DISCOVERY_STALL_CHECK_LIMIT {
                return Err(GeneralizedError::new(
                    "unreachable-bins termination lacks the discovery-stall evidence",
                ));
            }
            let required = ((self.config.minimum_visited_fraction
                * self.adaptation_histogram.bins() as f64)
                .ceil() as usize)
                .max(1);
            if self.last_check_visited_bins >= required {
                return Err(GeneralizedError::new(
                    "unreachable-bins termination saw enough visited bins for the gate",
                ));
            }
            if self.production_sweeps > 0 || self.production_histogram.total() > 0 {
                return Err(GeneralizedError::new(
                    "unreachable-bins termination cannot contain production progress",
                ));
            }
        }
        if self.is_adaptive()
            && self.config.max_adaptation_sweeps > 0
            && self.adaptation_sweeps >= self.config.max_adaptation_sweeps
        {
            return Err(GeneralizedError::new(
                "adaptive checkpoint has already reached the maximum-sweep guard",
            ));
        }
        if self.termination == Some(WangLandauTermination::Converged)
            && self.log_f > self.config.final_log_f
        {
            return Err(GeneralizedError::new(
                "converged checkpoint has not reached the configured final log_f",
            ));
        }
        Ok(())
    }

    fn finish_flat_histogram_sweep(&mut self) {
        let refinement_sweeps = self
            .adaptation_sweeps
            .saturating_sub(self.config.discovery_sweeps);
        if refinement_sweeps == 0
            || !refinement_sweeps.is_multiple_of(self.config.flatness_check_interval)
        {
            return;
        }
        self.flatness_checks = self.flatness_checks.saturating_add(1);

        // Discovery-plateau tracking: once the walk stops finding new bins,
        // the visited set is the reachable set. If the configured minimum
        // visited fraction demands more than that, the flatness gate is
        // unattainable and the estimate would silently burn sweeps until the
        // maximum-sweep guard — terminate loudly instead.
        let visited = self.adaptation_histogram.visited_bins();
        let required = ((self.config.minimum_visited_fraction
            * self.adaptation_histogram.bins() as f64)
            .ceil() as usize)
            .max(1);
        if visited == self.last_check_visited_bins {
            self.discovery_stall_checks = self.discovery_stall_checks.saturating_add(1);
        } else {
            self.discovery_stall_checks = 0;
        }
        self.last_check_visited_bins = visited;
        if visited < required && self.discovery_stall_checks >= DISCOVERY_STALL_CHECK_LIMIT {
            eprintln!(
                "Wang-Landau minimum_visited_fraction {} requires {required} visited bins, but \
                 bin discovery has plateaued at {visited} of {} bins over {} consecutive \
                 flatness checks; terminating with UnreachableBins — lower the fraction or \
                 narrow the macrostate axis to the reachable range",
                self.config.minimum_visited_fraction,
                self.adaptation_histogram.bins(),
                self.discovery_stall_checks
            );
            self.log_density.normalize_max_zero();
            self.phase = WangLandauPhase::Finished;
            self.termination = Some(WangLandauTermination::UnreachableBins);
            return;
        }

        if !self
            .adaptation_histogram
            .is_flat(self.config.flatness, self.config.minimum_visited_fraction)
        {
            return;
        }

        self.flatness_passes = self.flatness_passes.saturating_add(1);
        self.log_f *= 0.5;
        self.log_density.normalize_max_zero();
        self.adaptation_histogram.clear();
        self.refinement_visits = 0;
        // The cleared histogram restarts discovery; force the next check to
        // count as progress (usize::MAX is never a valid visited count).
        self.last_check_visited_bins = usize::MAX;
        self.discovery_stall_checks = 0;
        if self.log_f <= self.config.final_log_f {
            self.termination = Some(WangLandauTermination::Converged);
            self.freeze_for_production();
            return;
        }
        if self.config.one_over_t_threshold > 0.0 && self.log_f <= self.config.one_over_t_threshold
        {
            self.refinement = WangLandauRefinement::OneOverT;
            self.refinement_visits = 0;
        }
    }

    fn finish_one_over_t_sweep(&mut self) {
        self.log_f = self.one_over_t_increment().min(self.log_f);
        if self.log_f <= self.config.final_log_f {
            self.termination = Some(WangLandauTermination::Converged);
            self.log_density.normalize_max_zero();
            self.freeze_for_production();
        }
    }

    fn one_over_t_increment(&self) -> f64 {
        let represented = self.log_density.visited_bins().max(1) as f64;
        represented / self.refinement_visits.max(1) as f64
    }
}

/// Local lattice Wang–Landau kernel over a scalar energy axis.
#[derive(Debug, Clone)]
pub struct WangLandauCore<A, S = StandardStrategy> {
    axis: A,
    estimator: WangLandauState,
    strategy: S,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: EnergyPatch,
    energy_check_interval: u64,
    sweeps: u64,
    out_of_range_proposals: u64,
    last_visited_bin: Option<usize>,
}

impl<A> WangLandauCore<A, StandardStrategy>
where
    A: MacrostateAxis,
{
    pub fn new(axis: A, config: WangLandauConfig) -> Result<Self, GeneralizedError> {
        Self::with_strategy(axis, config, StandardStrategy::new())
    }
}

impl<A, S> WangLandauCore<A, S>
where
    A: MacrostateAxis,
{
    pub fn with_strategy(
        axis: A,
        config: WangLandauConfig,
        strategy: S,
    ) -> Result<Self, GeneralizedError> {
        let estimator = WangLandauState::new(axis.bins(), config)?;
        Ok(Self {
            axis,
            estimator,
            strategy,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: EnergyPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
            out_of_range_proposals: 0,
            last_visited_bin: None,
        })
    }

    fn from_state(
        axis: A,
        estimator: WangLandauState,
        strategy: S,
    ) -> Result<Self, GeneralizedError> {
        if axis.bins() != estimator.log_density().bins() {
            return Err(GeneralizedError::new(
                "Wang-Landau axis and checkpoint have different bin counts",
            ));
        }
        Ok(Self {
            axis,
            estimator,
            strategy,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: EnergyPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
            out_of_range_proposals: 0,
            last_visited_bin: None,
        })
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
    pub const fn axis(&self) -> &A {
        &self.axis
    }

    #[inline]
    pub const fn estimator(&self) -> &WangLandauState {
        &self.estimator
    }

    #[inline]
    pub const fn estimator_mut(&mut self) -> &mut WangLandauState {
        &mut self.estimator
    }

    /// Number of trial energies rejected because they were outside the fixed axis.
    #[inline]
    pub const fn out_of_range_proposals(&self) -> u64 {
        self.out_of_range_proposals
    }

    /// Number of completed kernel sweeps, including frozen production.
    #[inline]
    pub const fn sweeps(&self) -> u64 {
        self.sweeps
    }

    /// Checkpoint the adaptive state together with the physical axis values.
    pub fn save_snapshot(&self) -> Json {
        json!({
            "format": "cmc-rs-wang-landau-kernel-v1",
            "axis_centers": self.axis.centers(),
            "visit_schedule": visit_schedule_name(self.visit_schedule),
            "energy_check_interval": self.energy_check_interval,
            "sweeps": self.sweeps,
            "out_of_range_proposals": self.out_of_range_proposals,
            "last_visited_bin": self.last_visited_bin,
            "estimator": self.estimator.save_snapshot(),
        })
    }

    /// Restore a kernel checkpoint and reject an axis with matching length but
    /// different physical bin values.
    pub fn from_snapshot(axis: A, strategy: S, snapshot: &Json) -> Result<Self, GeneralizedError> {
        if snapshot["format"].as_str() != Some("cmc-rs-wang-landau-kernel-v1") {
            return Err(GeneralizedError::new(
                "unknown or missing Wang-Landau kernel checkpoint format",
            ));
        }
        let stored_centers = required_f64_array(snapshot, "axis_centers")?;
        let centers = axis.centers();
        if stored_centers.len() != centers.len()
            || stored_centers
                .iter()
                .zip(&centers)
                .any(|(&stored, &current)| {
                    (stored - current).abs() > 1e-12 * (1.0 + stored.abs().max(current.abs()))
                })
        {
            return Err(GeneralizedError::new(
                "Wang-Landau checkpoint axis does not match the requested axis",
            ));
        }
        let estimator = WangLandauState::load_snapshot(&snapshot["estimator"])?;
        let mut kernel = Self::from_state(axis, estimator, strategy)?;
        kernel.visit_schedule = parse_visit_schedule(required_str(snapshot, "visit_schedule")?)?;
        kernel.energy_check_interval = required_u64(snapshot, "energy_check_interval")?;
        kernel.sweeps = required_u64(snapshot, "sweeps")?;
        kernel.out_of_range_proposals = required_u64(snapshot, "out_of_range_proposals")?;
        kernel.last_visited_bin = match snapshot.get("last_visited_bin") {
            None => None,
            Some(value) if value.is_null() => None,
            Some(value) => {
                let bin = value.as_u64().ok_or_else(|| {
                    GeneralizedError::new("invalid Wang-Landau cached macrostate bin")
                })? as usize;
                if bin >= kernel.axis.bins() {
                    return Err(GeneralizedError::new(
                        "Wang-Landau cached macrostate bin is outside the axis",
                    ));
                }
                Some(bin)
            }
        };
        Ok(kernel)
    }
}

impl<H, A, S> Algorithm<H> for WangLandauCore<A, S>
where
    H: Hamiltonian + Proposable,
    A: MacrostateAxis,
    S: ProposalStrategy<H>,
{
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        if phase == SimulationPhase::Measurement && self.estimator.is_adaptive() {
            self.estimator.freeze_for_production();
        }
        if self.estimator.phase() == WangLandauPhase::Finished {
            return;
        }

        let sites = self
            .order
            .prepare(system.n_sites(), self.visit_schedule, rng);
        for &site in sites {
            let old_bin = self
                .axis
                .bin(system.energy)
                .expect("accepted energy lies outside the Wang-Landau axis");
            let proposed_spin = self.strategy.propose(model, system, site, rng);
            let movement = SiteSpinMove::new(site, proposed_spin.spin);
            let delta = system.evaluate_trial(model, &movement, &mut self.patch);
            let new_bin = self.axis.bin(system.energy + delta.energy);
            if new_bin.is_none() {
                self.out_of_range_proposals = self.out_of_range_proposals.saturating_add(1);
            }
            let log_acceptance = new_bin.map_or(f64::NEG_INFINITY, |bin| {
                self.estimator.log_weight_ratio(old_bin, bin)
                    + proposed_spin.log_reverse_over_forward
            });
            let accepted = accept_log_probability(log_acceptance, rng);
            if accepted {
                <System as TrialEvaluator<H, SiteSpinMove>>::commit_trial(
                    system,
                    &movement,
                    &self.patch,
                );
            }
            self.strategy.record_result(accepted);
            let visited_bin = if accepted {
                new_bin.expect("accepted Wang-Landau trial has a bin")
            } else {
                old_bin
            };
            self.estimator.record_visit(visited_bin);
            self.last_visited_bin = Some(visited_bin);
        }
        self.strategy
            .finish_sweep(phase.allows_adaptation() && self.estimator.is_adaptive());
        self.estimator.finish_sweep();

        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_lattice_cache(system, model).expect("Wang-Landau cache audit failed");
            self.estimator
                .validate_bin_count(self.axis.bins())
                .expect("Wang-Landau histogram/DOS cache audit failed");
            if let Some(bin) = self.last_visited_bin {
                audit_macrostate_bin(&self.axis, system.energy, bin)
                    .expect("Wang-Landau macrostate cache audit failed");
            }
        }
    }

    fn finish_run(&mut self) {
        self.estimator.finish_production();
    }

    fn name(&self) -> &'static str {
        "Wang-Landau density-of-states Metropolis-Hastings"
    }
}

/// Scheduler-ready exact-axis Wang–Landau simulation for small Ising graphs.
///
/// The exact axis is deliberately limited to 24 sites. It is the reference
/// implementation used to validate the generalized-ensemble stack; larger or
/// continuous systems should construct [`ClassicalMC`] with a user-supplied
/// [`BinnedAxis`](crate::BinnedAxis) or [`DiscreteAxis`](crate::DiscreteAxis).
pub struct IsingWangLandau {
    pub chain: ClassicalMC<
        IsingModel,
        WangLandauCore<crate::generalized::DiscreteAxis>,
        DefaultObservableSet<IsingModel>,
    >,
}

impl IsingWangLandau {
    #[inline]
    pub const fn estimator(&self) -> &WangLandauState {
        self.chain.algorithm.estimator()
    }

    #[inline]
    pub const fn estimator_mut(&mut self) -> &mut WangLandauState {
        self.chain.algorithm.estimator_mut()
    }
}

impl MonteCarlo for IsingWangLandau {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.chain.sweep(ctx);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        self.chain.measure(ctx);
    }

    fn on_phase_start(&mut self, phase: RunPhase, ctx: &mut Context<Self::Rng>) {
        self.chain.on_phase_start(phase, ctx);
    }

    fn on_phase_end(&mut self, phase: RunPhase, ctx: &mut Context<Self::Rng>) {
        self.chain.on_phase_end(phase, ctx);
    }

    fn name(&self) -> &'static str {
        self.chain.name()
    }
}

impl FromParams for IsingWangLandau {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let pbc = parse_bool(params, "pbc", true)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        IsingModel::from_hamiltonian_params(params)?;
        if lattice.n_sites > 24 {
            return Err(invalid_config(
                "lattice",
                "exact-axis Ising Wang-Landau is limited to 24 sites",
            ));
        }
        wang_landau_config_from_params(params)?;
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(0.0);
        if !beta.is_finite() || beta < 0.0 {
            return Err(invalid_config("beta", "must be finite and non-negative"));
        }
        Ok(())
    }

    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let pbc = parse_bool(params, "pbc", true)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        let model = IsingModel::from_hamiltonian_params(params)?;
        let axis = enumerate_reference_axis(&lattice, &model)?;
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(0.0);
        let mut system = System::new(lattice, 1, 0.0, beta);
        let initial_state = params
            .get::<String>("initial_state")
            .unwrap_or_else(|| "hot".to_string())
            .to_ascii_lowercase();
        for site in 0..system.n_sites() {
            let spin = match initial_state.as_str() {
                "hot" | "random" => model.random_spin(rng),
                "cold" | "ordered" => model.ordered_spin(),
                _ => {
                    return Err(invalid_config(
                        "initial_state",
                        "expected hot/random or cold/ordered",
                    ));
                }
            };
            system.spin_at_mut(site, 1).copy_from_slice(&spin);
        }
        system.recompute_energy(&model);
        if axis.bin(system.energy).is_none() {
            return Err(invalid_config(
                "initial_state",
                "initial energy is not represented by the Wang-Landau axis",
            ));
        }
        let algorithm = WangLandauCore::new(axis, wang_landau_config_from_params(params)?)
            .map_err(generalized_config_error)?;
        Ok(Self {
            chain: ClassicalMC::new(system, model, algorithm),
        })
    }
}

/// Carlo.rs controller that runs adaptation to convergence, then a fixed
/// number of frozen multicanonical production sweeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WangLandauRunControl {
    production_sweeps: u64,
}

impl WangLandauRunControl {
    pub const fn new(production_sweeps: u64) -> Self {
        Self { production_sweeps }
    }

    #[inline]
    pub const fn production_sweeps(self) -> u64 {
        self.production_sweeps
    }
}

impl<H, A, S, O> AdaptiveRunControl<ClassicalMC<H, WangLandauCore<A, S>, O>>
    for WangLandauRunControl
where
    H: Hamiltonian + Proposable,
    A: MacrostateAxis,
    S: ProposalStrategy<H>,
    O: ObservableSet<H>,
{
    fn initial_phase(&self) -> RunPhase {
        RunPhase::Thermalization
    }

    fn after_sweep(
        &mut self,
        mc: &ClassicalMC<H, WangLandauCore<A, S>, O>,
        ctx: &Context<<ClassicalMC<H, WangLandauCore<A, S>, O> as carlo_rs::MonteCarlo>::Rng>,
    ) -> RunDecision {
        decision_for_state(
            mc.algorithm.estimator(),
            ctx.phase(),
            self.production_sweeps,
        )
    }
}

impl AdaptiveRunControl<IsingWangLandau> for WangLandauRunControl {
    fn initial_phase(&self) -> RunPhase {
        RunPhase::Thermalization
    }

    fn after_sweep(
        &mut self,
        mc: &IsingWangLandau,
        ctx: &Context<<IsingWangLandau as MonteCarlo>::Rng>,
    ) -> RunDecision {
        decision_for_state(mc.estimator(), ctx.phase(), self.production_sweeps)
    }
}

fn decision_for_state(
    state: &WangLandauState,
    phase: RunPhase,
    production_sweeps: u64,
) -> RunDecision {
    match phase {
        RunPhase::Initialization | RunPhase::Thermalization => match state.phase() {
            WangLandauPhase::Discovery | WangLandauPhase::Adaptation => {
                RunDecision::ContinueAdaptation
            }
            WangLandauPhase::FrozenProduction => {
                if production_sweeps == 0 {
                    RunDecision::Stop
                } else {
                    RunDecision::BeginProduction
                }
            }
            WangLandauPhase::Finished => RunDecision::Stop,
        },
        RunPhase::Measurement => {
            if state.production_sweeps() >= production_sweeps {
                RunDecision::Stop
            } else {
                RunDecision::ContinueProduction
            }
        }
        RunPhase::Finished => RunDecision::Stop,
    }
}

fn enumerate_reference_axis(
    lattice: &crate::lattice::graph::CsrLattice,
    model: &IsingModel,
) -> Result<crate::generalized::DiscreteAxis, CarloError> {
    crate::generalized::enumerate_ising_density_of_states(lattice, model)
        .and_then(|exact| exact.axis())
        .map_err(generalized_config_error)
}

fn wang_landau_config_from_params(params: &Params) -> Result<WangLandauConfig, CarloError> {
    let defaults = WangLandauConfig::default();
    let config = WangLandauConfig {
        initial_log_f: parse_param(params, "wl_initial_log_f")?.unwrap_or(defaults.initial_log_f),
        final_log_f: parse_param(params, "wl_final_log_f")?.unwrap_or(defaults.final_log_f),
        flatness: parse_param(params, "wl_flatness")?.unwrap_or(defaults.flatness),
        flatness_check_interval: parse_param(params, "wl_flatness_check_interval")?
            .unwrap_or(defaults.flatness_check_interval),
        discovery_sweeps: parse_param(params, "wl_discovery_sweeps")?
            .unwrap_or(defaults.discovery_sweeps),
        one_over_t_threshold: parse_param(params, "wl_one_over_t_threshold")?
            .unwrap_or(defaults.one_over_t_threshold),
        max_adaptation_sweeps: parse_param(params, "wl_max_adaptation_sweeps")?
            .unwrap_or(defaults.max_adaptation_sweeps),
        minimum_visited_fraction: parse_param(params, "wl_minimum_visited_fraction")?
            .unwrap_or(defaults.minimum_visited_fraction),
    };
    config.validate().map_err(generalized_config_error)?;
    Ok(config)
}

fn generalized_config_error(error: GeneralizedError) -> CarloError {
    invalid_config("wang_landau", error.to_string())
}

fn invalid_config(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}

fn visit_schedule_name(schedule: VisitSchedule) -> &'static str {
    match schedule {
        VisitSchedule::Sequential => "sequential",
        VisitSchedule::RandomPermutation => "random_permutation",
    }
}

fn parse_visit_schedule(value: &str) -> Result<VisitSchedule, GeneralizedError> {
    match value {
        "sequential" => Ok(VisitSchedule::Sequential),
        "random_permutation" => Ok(VisitSchedule::RandomPermutation),
        _ => Err(GeneralizedError::new(format!(
            "unknown Wang-Landau visit schedule `{value}`"
        ))),
    }
}

fn required_str<'a>(value: &'a Json, name: &str) -> Result<&'a str, GeneralizedError> {
    value[name].as_str().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })
}

fn required_f64(value: &Json, name: &str) -> Result<f64, GeneralizedError> {
    let number = value[name].as_f64().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })?;
    if !number.is_finite() {
        return Err(GeneralizedError::new(format!(
            "checkpoint field `{name}` is non-finite"
        )));
    }
    Ok(number)
}

fn required_u64(value: &Json, name: &str) -> Result<u64, GeneralizedError> {
    value[name].as_u64().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })
}

fn required_f64_array(value: &Json, name: &str) -> Result<Vec<f64>, GeneralizedError> {
    let array = value[name].as_array().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })?;
    array
        .iter()
        .enumerate()
        .map(|(index, entry)| {
            let number = entry.as_f64().ok_or_else(|| {
                GeneralizedError::new(format!("checkpoint field `{name}[{index}]` is invalid"))
            })?;
            if number.is_finite() {
                Ok(number)
            } else {
                Err(GeneralizedError::new(format!(
                    "checkpoint field `{name}[{index}]` is non-finite"
                )))
            }
        })
        .collect()
}

fn required_u64_array(value: &Json, name: &str) -> Result<Vec<u64>, GeneralizedError> {
    let array = value[name].as_array().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })?;
    array
        .iter()
        .enumerate()
        .map(|(index, entry)| {
            entry.as_u64().ok_or_else(|| {
                GeneralizedError::new(format!("checkpoint field `{name}[{index}]` is invalid"))
            })
        })
        .collect()
}

fn required_bool_array(value: &Json, name: &str) -> Result<Vec<bool>, GeneralizedError> {
    let array = value[name].as_array().ok_or_else(|| {
        GeneralizedError::new(format!("checkpoint field `{name}` is missing or invalid"))
    })?;
    array
        .iter()
        .enumerate()
        .map(|(index, entry)| {
            entry.as_bool().ok_or_else(|| {
                GeneralizedError::new(format!("checkpoint field `{name}[{index}]` is invalid"))
            })
        })
        .collect()
}
