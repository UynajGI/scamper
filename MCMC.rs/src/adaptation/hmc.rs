use serde::{Deserialize, Serialize};

use crate::McmcError;

/// Nesterov dual averaging for a positive HMC step size.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DualAveraging {
    target_acceptance: f64,
    gamma: f64,
    t0: f64,
    kappa: f64,
    mu: f64,
    log_step_size: f64,
    log_averaged_step_size: f64,
    error_sum: f64,
    iteration: u64,
    frozen: bool,
}

impl DualAveraging {
    pub fn new(initial_step_size: f64, target_acceptance: f64) -> Result<Self, McmcError> {
        validate_step_size(initial_step_size)?;
        if !target_acceptance.is_finite() || !(0.0..1.0).contains(&target_acceptance) {
            return Err(McmcError::InvalidConfig(
                "HMC target acceptance must lie strictly between zero and one".to_string(),
            ));
        }
        let log_step_size = initial_step_size.ln();
        Ok(Self {
            target_acceptance,
            gamma: 0.05,
            t0: 10.0,
            kappa: 0.75,
            mu: (10.0 * initial_step_size).ln(),
            log_step_size,
            log_averaged_step_size: log_step_size,
            error_sum: 0.0,
            iteration: 0,
            frozen: false,
        })
    }

    pub fn observe(&mut self, acceptance_probability: f64) -> Result<f64, McmcError> {
        if self.frozen {
            return Err(McmcError::AdaptationFrozen);
        }
        if !acceptance_probability.is_finite() || !(0.0..=1.0).contains(&acceptance_probability) {
            return Err(McmcError::InvalidConfig(
                "HMC acceptance statistic must lie between zero and one".to_string(),
            ));
        }
        self.iteration = self.iteration.saturating_add(1);
        let iteration = self.iteration as f64;
        let weight = 1.0 / (iteration + self.t0);
        self.error_sum = (1.0 - weight) * self.error_sum
            + weight * (self.target_acceptance - acceptance_probability);
        self.log_step_size = self.mu - iteration.sqrt() / self.gamma * self.error_sum;
        self.log_step_size = self.log_step_size.clamp(-20.0, 20.0);
        let averaging_weight = iteration.powf(-self.kappa);
        self.log_averaged_step_size = averaging_weight * self.log_step_size
            + (1.0 - averaging_weight) * self.log_averaged_step_size;
        Ok(self.current_step_size())
    }

    pub fn restart(&mut self, initial_step_size: f64) -> Result<(), McmcError> {
        validate_step_size(initial_step_size)?;
        let log_step_size = initial_step_size.ln();
        self.mu = (10.0 * initial_step_size).ln();
        self.log_step_size = log_step_size;
        self.log_averaged_step_size = log_step_size;
        self.error_sum = 0.0;
        self.iteration = 0;
        self.frozen = false;
        Ok(())
    }

    pub fn freeze(&mut self) -> f64 {
        self.frozen = true;
        self.log_step_size = self.log_averaged_step_size;
        self.current_step_size()
    }

    pub fn current_step_size(&self) -> f64 {
        self.log_step_size.exp()
    }

    pub fn averaged_step_size(&self) -> f64 {
        self.log_averaged_step_size.exp()
    }

    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }
}

/// Fast/slow/fast warmup partition used for mass-matrix windows.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct WarmupWindowConfig {
    pub initial_buffer: u64,
    pub terminal_buffer: u64,
    pub initial_window: u64,
}

impl WarmupWindowConfig {
    pub fn default_for(total_warmup: u64) -> Result<Self, McmcError> {
        if total_warmup == 0 {
            return Err(McmcError::InvalidConfig(
                "HMC warmup length must be positive".to_string(),
            ));
        }
        if total_warmup < 20 {
            return Ok(Self {
                initial_buffer: total_warmup,
                terminal_buffer: 0,
                initial_window: 1,
            });
        }
        let initial_buffer = (total_warmup * 15 / 100).max(5);
        let terminal_buffer = (total_warmup * 10 / 100).max(5);
        let slow = total_warmup
            .saturating_sub(initial_buffer)
            .saturating_sub(terminal_buffer);
        Ok(Self {
            initial_buffer,
            terminal_buffer,
            initial_window: slow.clamp(1, 25),
        })
    }

    pub fn validate(&self, total_warmup: u64, metric_adaptation: bool) -> Result<(), McmcError> {
        if total_warmup == 0 || self.initial_window == 0 {
            return Err(McmcError::InvalidConfig(
                "warmup length and initial metric window must be positive".to_string(),
            ));
        }
        let buffers = self.initial_buffer.saturating_add(self.terminal_buffer);
        if buffers > total_warmup {
            return Err(McmcError::InvalidConfig(
                "HMC warmup buffers exceed total warmup length".to_string(),
            ));
        }
        if metric_adaptation && (buffers >= total_warmup || self.terminal_buffer == 0) {
            return Err(McmcError::InvalidConfig(
                "metric adaptation requires a non-empty slow window and terminal buffer"
                    .to_string(),
            ));
        }
        Ok(())
    }
}

/// Position-covariance geometry adapted during slow warmup windows.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum MetricAdaptation {
    None,
    Diagonal { regularization: f64 },
    Dense { regularization: f64 },
}

impl MetricAdaptation {
    fn validate(self) -> Result<(), McmcError> {
        match self {
            Self::None => Ok(()),
            Self::Diagonal { regularization } | Self::Dense { regularization }
                if regularization.is_finite() && regularization > 0.0 =>
            {
                Ok(())
            }
            Self::Diagonal { .. } | Self::Dense { .. } => Err(McmcError::InvalidConfig(
                "metric regularization must be finite and positive".to_string(),
            )),
        }
    }

    const fn enabled(self) -> bool {
        !matches!(self, Self::None)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
enum CovarianceAccumulator {
    None,
    Diagonal {
        mean: Vec<f64>,
        m2: Vec<f64>,
        count: u64,
        regularization: f64,
    },
    Dense {
        dimension: usize,
        mean: Vec<f64>,
        m2: Vec<f64>,
        delta: Vec<f64>,
        count: u64,
        regularization: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum MetricUpdate {
    Diagonal(Vec<f64>),
    Dense {
        dimension: usize,
        covariance: Vec<f64>,
        jitter: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct WarmupObservation {
    pub step_size: f64,
    pub metric_update: Option<MetricUpdate>,
}

/// Serializable HMC warmup controller combining dual averaging and windows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HmcWarmup {
    dimension: usize,
    total_warmup: u64,
    iteration: u64,
    windows: WarmupWindowConfig,
    window_ends: Vec<u64>,
    next_window: usize,
    dual_averaging: DualAveraging,
    metric_adaptation: MetricAdaptation,
    covariance: CovarianceAccumulator,
    frozen: bool,
}

impl HmcWarmup {
    pub fn new(
        dimension: usize,
        total_warmup: u64,
        initial_step_size: f64,
        target_acceptance: f64,
        metric_adaptation: MetricAdaptation,
        windows: WarmupWindowConfig,
    ) -> Result<Self, McmcError> {
        if dimension == 0 {
            return Err(McmcError::InvalidConfig(
                "HMC adaptation dimension must be positive".to_string(),
            ));
        }
        metric_adaptation.validate()?;
        windows.validate(total_warmup, metric_adaptation.enabled())?;
        let window_ends = build_window_ends(total_warmup, windows);
        let covariance = CovarianceAccumulator::new(dimension, metric_adaptation);
        Ok(Self {
            dimension,
            total_warmup,
            iteration: 0,
            windows,
            window_ends,
            next_window: 0,
            dual_averaging: DualAveraging::new(initial_step_size, target_acceptance)?,
            metric_adaptation,
            covariance,
            frozen: false,
        })
    }

    pub(crate) fn observe(
        &mut self,
        acceptance_probability: f64,
        position: &[f64],
    ) -> Result<WarmupObservation, McmcError> {
        if self.frozen || self.iteration >= self.total_warmup {
            return Err(McmcError::AdaptationFrozen);
        }
        if position.len() != self.dimension {
            return Err(McmcError::DimensionMismatch {
                expected: self.dimension,
                actual: position.len(),
            });
        }
        let mut step_size = self.dual_averaging.observe(acceptance_probability)?;
        self.iteration = self.iteration.saturating_add(1);

        if self.in_slow_window() {
            self.covariance.observe(position)?;
        }

        let mut metric_update = None;
        if self
            .window_ends
            .get(self.next_window)
            .is_some_and(|end| *end == self.iteration)
        {
            metric_update = self.covariance.estimate();
            self.covariance
                .reset(self.dimension, self.metric_adaptation);
            self.next_window = self.next_window.saturating_add(1);
            if metric_update.is_some() {
                self.dual_averaging.restart(step_size)?;
                step_size = self.dual_averaging.current_step_size();
            }
        }

        if self.iteration == self.total_warmup {
            step_size = self.dual_averaging.freeze();
            self.frozen = true;
        }
        Ok(WarmupObservation {
            step_size,
            metric_update,
        })
    }

    pub fn finish(&mut self) -> Result<f64, McmcError> {
        if self.iteration != self.total_warmup {
            return Err(McmcError::InvalidConfig(format!(
                "HMC warmup expected {} transitions but observed {}",
                self.total_warmup, self.iteration
            )));
        }
        if !self.frozen {
            self.frozen = true;
            return Ok(self.dual_averaging.freeze());
        }
        Ok(self.dual_averaging.current_step_size())
    }

    pub const fn iteration(&self) -> u64 {
        self.iteration
    }

    pub const fn total_warmup(&self) -> u64 {
        self.total_warmup
    }

    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }

    pub fn step_size(&self) -> f64 {
        self.dual_averaging.current_step_size()
    }

    fn in_slow_window(&self) -> bool {
        let slow_start = self.windows.initial_buffer;
        let slow_end = self
            .total_warmup
            .saturating_sub(self.windows.terminal_buffer);
        self.iteration > slow_start && self.iteration <= slow_end
    }
}

impl CovarianceAccumulator {
    fn new(dimension: usize, adaptation: MetricAdaptation) -> Self {
        match adaptation {
            MetricAdaptation::None => Self::None,
            MetricAdaptation::Diagonal { regularization } => Self::Diagonal {
                mean: vec![0.0; dimension],
                m2: vec![0.0; dimension],
                count: 0,
                regularization,
            },
            MetricAdaptation::Dense { regularization } => Self::Dense {
                dimension,
                mean: vec![0.0; dimension],
                m2: vec![0.0; dimension.saturating_mul(dimension)],
                delta: vec![0.0; dimension],
                count: 0,
                regularization,
            },
        }
    }

    fn observe(&mut self, position: &[f64]) -> Result<(), McmcError> {
        if position.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "metric adaptation observation must be finite".to_string(),
            ));
        }
        match self {
            Self::None => Ok(()),
            Self::Diagonal {
                mean, m2, count, ..
            } => {
                if position.len() != mean.len() {
                    return Err(McmcError::DimensionMismatch {
                        expected: mean.len(),
                        actual: position.len(),
                    });
                }
                *count = (*count).saturating_add(1);
                let denominator = *count as f64;
                for ((mean, m2), value) in mean
                    .iter_mut()
                    .zip(m2.iter_mut())
                    .zip(position.iter().copied())
                {
                    let delta = value - *mean;
                    *mean += delta / denominator;
                    *m2 += delta * (value - *mean);
                }
                Ok(())
            }
            Self::Dense {
                dimension,
                mean,
                m2,
                delta,
                count,
                ..
            } => {
                if position.len() != *dimension {
                    return Err(McmcError::DimensionMismatch {
                        expected: *dimension,
                        actual: position.len(),
                    });
                }
                *count = (*count).saturating_add(1);
                let denominator = *count as f64;
                for (index, value) in position.iter().copied().enumerate() {
                    delta[index] = value - mean[index];
                    mean[index] += delta[index] / denominator;
                }
                for (row, delta_row) in delta.iter().copied().enumerate() {
                    for (column, value) in position.iter().copied().enumerate() {
                        m2[row * *dimension + column] += delta_row * (value - mean[column]);
                    }
                }
                Ok(())
            }
        }
    }

    fn estimate(&self) -> Option<MetricUpdate> {
        match self {
            Self::None => None,
            Self::Diagonal {
                m2,
                count,
                regularization,
                ..
            } if *count >= 2 => {
                let denominator = (*count - 1) as f64;
                Some(MetricUpdate::Diagonal(
                    m2.iter()
                        .map(|value| value / denominator + *regularization)
                        .collect(),
                ))
            }
            Self::Dense {
                dimension,
                m2,
                count,
                regularization,
                ..
            } if *count >= 2 => {
                let denominator = (*count - 1) as f64;
                let mut covariance = vec![0.0; (*dimension).saturating_mul(*dimension)];
                for row in 0..*dimension {
                    for column in 0..*dimension {
                        covariance[row * *dimension + column] = 0.5
                            * (m2[row * *dimension + column] + m2[column * *dimension + row])
                            / denominator;
                    }
                    covariance[row * *dimension + row] += *regularization;
                }
                Some(MetricUpdate::Dense {
                    dimension: *dimension,
                    covariance,
                    jitter: *regularization,
                })
            }
            Self::Diagonal { .. } | Self::Dense { .. } => None,
        }
    }

    fn reset(&mut self, dimension: usize, adaptation: MetricAdaptation) {
        *self = Self::new(dimension, adaptation);
    }
}

fn build_window_ends(total_warmup: u64, windows: WarmupWindowConfig) -> Vec<u64> {
    let slow_end = total_warmup.saturating_sub(windows.terminal_buffer);
    let mut start = windows.initial_buffer;
    let mut size = windows.initial_window;
    let mut ends = Vec::new();
    while start < slow_end {
        let remaining = slow_end - start;
        let current = if size.saturating_mul(2) > remaining {
            remaining
        } else {
            size.min(remaining)
        };
        if current == 0 {
            break;
        }
        start += current;
        ends.push(start);
        size = size.saturating_mul(2).max(1);
    }
    ends
}

fn validate_step_size(step_size: f64) -> Result<(), McmcError> {
    if step_size.is_finite() && step_size > 0.0 {
        Ok(())
    } else {
        Err(McmcError::InvalidConfig(
            "HMC step size must be finite and positive".to_string(),
        ))
    }
}
