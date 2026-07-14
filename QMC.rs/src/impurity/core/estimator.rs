//! Derived and directed-loop improved estimators shared by impurity solvers.

use carlo_rs::{CarloError, Evaluator};

use crate::impurity::core::imaginary_time::{directed_distance, PropagationDirection};

/// Connected static susceptibility from ensemble means of a time-averaged spin.
pub fn connected_susceptibility(beta: f64, mean_m: f64, mean_m_squared: f64) -> f64 {
    beta * (mean_m_squared - mean_m * mean_m)
}

/// Register a jackknife-safe connected susceptibility in Carlo.rs.
pub fn register_connected_susceptibility(
    evaluator: &mut Evaluator,
    name: &str,
    magnetization: &str,
    magnetization_squared: &str,
    beta: f64,
) -> Result<(), CarloError> {
    evaluator.evaluate(name, &[magnetization, magnetization_squared], move |args| {
        (&args[1] - &(&args[0] * &args[0])) * beta
    })
}

/// Whether a discontinuity carries `S+` or `S-`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinFlipOperator {
    Raise,
    Lower,
}

impl SpinFlipOperator {
    /// Operator that changes `incoming` into `outgoing`.
    pub fn from_transition(incoming: i8, outgoing: i8) -> Option<Self> {
        match (incoming, outgoing) {
            (-1, 1) => Some(Self::Raise),
            (1, -1) => Some(Self::Lower),
            _ => None,
        }
    }

    pub const fn opposite(self) -> Self {
        match self {
            Self::Raise => Self::Lower,
            Self::Lower => Self::Raise,
        }
    }
}

/// One segment swept out by a directed-loop head.  Its contribution is to the
/// normal correlator `G+- + G-+` when tail/head operators are opposite and to
/// the anomalous correlator `G++ + G--` when they are equal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoopSegment {
    pub tail_tau: f64,
    pub from_tau: f64,
    pub to_tau: f64,
    pub direction: PropagationDirection,
    pub normal: bool,
}

/// Histogram accumulator for on-the-fly transverse worm estimators.
#[derive(Debug, Clone, PartialEq)]
pub struct TransverseLoopAccumulator {
    bins: usize,
    normal: Vec<f64>,
    anomalous: Vec<f64>,
    completed_loops: u64,
    path_length: f64,
}

impl Default for TransverseLoopAccumulator {
    fn default() -> Self {
        Self::new(64)
    }
}

impl TransverseLoopAccumulator {
    pub fn new(bins: usize) -> Self {
        let bins = bins.max(1);
        Self {
            bins,
            normal: vec![0.0; bins],
            anomalous: vec![0.0; bins],
            completed_loops: 0,
            path_length: 0.0,
        }
    }

    pub fn bins(&self) -> usize {
        self.bins
    }

    pub fn clear(&mut self) {
        self.normal.fill(0.0);
        self.anomalous.fill(0.0);
        self.completed_loops = 0;
        self.path_length = 0.0;
    }

    /// Commit a successfully closed loop. Aborted-loop journals must not call
    /// this method, because they do not belong to the sampled extended ensemble.
    pub fn commit_loop(&mut self, beta: f64, segments: &[LoopSegment]) {
        if !beta.is_finite() || beta <= 0.0 {
            return;
        }
        for segment in segments {
            self.add_segment(beta, *segment);
        }
        self.completed_loops += 1;
    }

    /// Exact zero-vertex estimator: a free spin has
    /// `<Sx(tau)Sx(0)> = <Sy(tau)Sy(0)> = 1/4`.
    pub fn commit_free_loop(&mut self, beta: f64) {
        if !beta.is_finite() || beta <= 0.0 {
            return;
        }
        self.normal.iter_mut().for_each(|value| *value += 1.0);
        self.path_length += beta;
        self.completed_loops += 1;
    }

    fn add_segment(&mut self, beta: f64, segment: LoopSegment) {
        let length = directed_distance(segment.from_tau, segment.to_tau, beta, segment.direction);
        if length <= 0.0 {
            return;
        }
        self.path_length += length;
        let bin_width = beta / self.bins as f64;
        let start_delta = match segment.direction {
            PropagationDirection::Forward => (segment.from_tau - segment.tail_tau).rem_euclid(beta),
            PropagationDirection::Backward => {
                (segment.tail_tau - segment.from_tau).rem_euclid(beta)
            }
        };
        let target = if segment.normal {
            &mut self.normal
        } else {
            &mut self.anomalous
        };

        // Deposit exact overlap length in each periodic separation bin.  The
        // resulting histogram is independent of vertex segmentation.
        let mut remaining = length;
        let mut cursor = start_delta;
        while remaining > 0.0 {
            let wrapped = cursor.rem_euclid(beta);
            let bin = ((wrapped / bin_width).floor() as usize).min(self.bins - 1);
            let edge = (bin + 1) as f64 * bin_width;
            let available = (edge - wrapped).max(f64::EPSILON);
            let piece = remaining.min(available);
            target[bin] += piece / bin_width;
            cursor += piece;
            remaining -= piece;
            if piece <= f64::EPSILON {
                break;
            }
        }
    }

    /// Return and reset the measurements accumulated since the previous call.
    pub fn take_sample(&mut self, beta: f64) -> TransverseCorrelationSample {
        let completed = self.completed_loops;
        let normalization = completed.max(1) as f64;
        let normal: Vec<f64> = self.normal.iter().map(|v| v / normalization).collect();
        let anomalous: Vec<f64> = self.anomalous.iter().map(|v| v / normalization).collect();
        let sampled_x: Vec<f64> = normal
            .iter()
            .zip(&anomalous)
            .map(|(n, a)| 0.25 * (n + a))
            .collect();
        let sampled_y: Vec<f64> = normal
            .iter()
            .zip(&anomalous)
            .map(|(n, a)| 0.25 * (n - a))
            .collect();
        let bin_width = beta / self.bins as f64;
        let susceptibility_x = sampled_x.iter().sum::<f64>() * bin_width;
        let susceptibility_y = sampled_y.iter().sum::<f64>() * bin_width;
        let half = self.bins / 2;
        let sample = TransverseCorrelationSample {
            normal,
            anomalous,
            sampled_x_half: sampled_x[half],
            sampled_y_half: sampled_y[half],
            sampled_x,
            sampled_y,
            susceptibility_x,
            susceptibility_y,
            completed_loops: completed,
            mean_path_length_beta: if completed == 0 || beta <= 0.0 {
                0.0
            } else {
                self.path_length / (completed as f64 * beta)
            },
        };
        self.clear();
        sample
    }
}

/// One measurement block of transverse directed-loop estimators.
#[derive(Debug, Clone, PartialEq)]
pub struct TransverseCorrelationSample {
    pub normal: Vec<f64>,
    pub anomalous: Vec<f64>,
    pub sampled_x: Vec<f64>,
    pub sampled_y: Vec<f64>,
    pub sampled_x_half: f64,
    pub sampled_y_half: f64,
    pub susceptibility_x: f64,
    pub susceptibility_y: f64,
    pub completed_loops: u64,
    pub mean_path_length_beta: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn connected_subtracts_the_one_point_piece() {
        assert!((connected_susceptibility(4.0, 0.25, 0.125) - 0.25).abs() < 1e-14);
    }

    #[test]
    fn free_loop_gives_exact_transverse_spin_correlations() {
        let mut accumulator = TransverseLoopAccumulator::new(16);
        accumulator.commit_free_loop(8.0);
        let sample = accumulator.take_sample(8.0);
        assert!(sample.sampled_x.iter().all(|v| (*v - 0.25).abs() < 1e-14));
        assert!(sample.sampled_y.iter().all(|v| (*v - 0.25).abs() < 1e-14));
        assert!((sample.susceptibility_x - 2.0).abs() < 1e-14);
    }

    #[test]
    fn aborted_journal_does_not_change_accumulator() {
        let mut accumulator = TransverseLoopAccumulator::new(8);
        let _journal = [LoopSegment {
            tail_tau: 0.0,
            from_tau: 0.0,
            to_tau: 1.0,
            direction: PropagationDirection::Forward,
            normal: true,
        }];
        let sample = accumulator.take_sample(2.0);
        assert_eq!(sample.completed_loops, 0);
        assert!(sample.sampled_x.iter().all(|v| v.abs() < 1e-14));
    }
}
