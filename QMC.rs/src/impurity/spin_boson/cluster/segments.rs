//! Continuous-time piecewise-constant worldlines for the longitudinal
//! spin-boson cluster solver.

use crate::impurity::ImpurityError;

const TIME_EPSILON_FACTOR: f64 = 64.0;

/// A non-wrapping half-open imaginary-time interval `[start, end)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeInterval {
    start: f64,
    end: f64,
}

impl TimeInterval {
    /// Construct a non-empty interval inside `[0, beta]`.
    pub fn new(start: f64, end: f64, beta: f64) -> Result<Self, ImpurityError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        if !start.is_finite() || !end.is_finite() || start < 0.0 || end > beta || start >= end {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "invalid imaginary-time interval [{start}, {end}) for beta={beta}"
            )));
        }
        Ok(Self { start, end })
    }

    /// Left endpoint.
    #[inline]
    pub const fn start(self) -> f64 {
        self.start
    }

    /// Right endpoint.
    #[inline]
    pub const fn end(self) -> f64 {
        self.end
    }

    /// Interval length.
    #[inline]
    pub fn length(self) -> f64 {
        self.end - self.start
    }
}

/// One circular segment between consecutive auxiliary cut points.
///
/// Every segment has one interval except the segment crossing `tau=0`, which
/// has two non-wrapping pieces.
#[derive(Debug, Clone, PartialEq)]
pub struct WorldlineSegment {
    spin: i8,
    pieces: [TimeInterval; 2],
    piece_count: u8,
    length: f64,
}

impl WorldlineSegment {
    fn from_one(spin: i8, piece: TimeInterval) -> Result<Self, ImpurityError> {
        validate_spin(spin)?;
        Ok(Self {
            spin,
            pieces: [piece, piece],
            piece_count: 1,
            length: piece.length(),
        })
    }

    fn from_two(
        spin: i8,
        first: TimeInterval,
        second: TimeInterval,
    ) -> Result<Self, ImpurityError> {
        validate_spin(spin)?;
        let length = first.length() + second.length();
        if !length.is_finite() || length <= 0.0 {
            return Err(ImpurityError::InvalidConfiguration(
                "a circular segment must have positive finite length".into(),
            ));
        }
        Ok(Self {
            spin,
            pieces: [first, second],
            piece_count: 2,
            length,
        })
    }

    /// Segment spin in the sampled `sigma_z` basis.
    #[inline]
    pub const fn spin(&self) -> i8 {
        self.spin
    }

    /// Non-wrapping interval pieces.
    #[inline]
    pub fn pieces(&self) -> &[TimeInterval] {
        &self.pieces[..usize::from(self.piece_count)]
    }

    /// Total circular length.
    #[inline]
    pub const fn length(&self) -> f64 {
        self.length
    }
}

/// A periodic spin-1/2 worldline represented by its real spin-flip times.
#[derive(Debug, Clone, PartialEq)]
pub struct LongitudinalWorldline {
    beta: f64,
    spin_at_zero: i8,
    kinks: Vec<f64>,
}

impl LongitudinalWorldline {
    /// Construct a constant worldline.
    pub fn new(beta: f64, spin_at_zero: i8) -> Result<Self, ImpurityError> {
        let worldline = Self {
            beta,
            spin_at_zero,
            kinks: Vec::new(),
        };
        worldline.validate()?;
        Ok(worldline)
    }

    /// Construct from a sorted list of real spin-flip times.
    pub fn from_kinks(beta: f64, spin_at_zero: i8, kinks: Vec<f64>) -> Result<Self, ImpurityError> {
        let worldline = Self {
            beta,
            spin_at_zero,
            kinks,
        };
        worldline.validate()?;
        Ok(worldline)
    }

    /// Inverse temperature.
    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// Spin immediately to the right of `tau=0`.
    #[inline]
    pub const fn spin_at_zero(&self) -> i8 {
        self.spin_at_zero
    }

    /// Sorted real spin-flip times in `(0, beta)`.
    #[inline]
    pub fn kinks(&self) -> &[f64] {
        &self.kinks
    }

    /// Number of real spin flips.
    #[inline]
    pub fn kink_count(&self) -> usize {
        self.kinks.len()
    }

    /// Spin at periodic imaginary time `tau`.
    pub fn spin_at(&self, tau: f64) -> i8 {
        let wrapped = tau.rem_euclid(self.beta);
        let flips = self.kinks.partition_point(|time| *time <= wrapped);
        if flips % 2 == 0 {
            self.spin_at_zero
        } else {
            -self.spin_at_zero
        }
    }

    /// Exact imaginary-time average of sampled `sigma_z`.
    pub fn integrated_sigma_z(&self) -> f64 {
        let mut previous = 0.0;
        let mut spin = self.spin_at_zero;
        let mut integral = 0.0;
        for &time in &self.kinks {
            integral += f64::from(spin) * (time - previous);
            spin = -spin;
            previous = time;
        }
        integral += f64::from(spin) * (self.beta - previous);
        integral / self.beta
    }

    /// Exact periodic correlation average
    /// `beta^-1 int d tau sigma_z(tau) sigma_z(tau + delta_tau)`.
    pub fn correlation_sigma_z(&self, delta_tau: f64) -> f64 {
        if self.kinks.is_empty() {
            return 1.0;
        }
        let delta = delta_tau.rem_euclid(self.beta);
        let mut boundaries = Vec::with_capacity(2 * self.kinks.len() + 2);
        boundaries.push(0.0);
        boundaries.push(self.beta);
        boundaries.extend(self.kinks.iter().copied());
        boundaries.extend(
            self.kinks
                .iter()
                .map(|time| (time - delta).rem_euclid(self.beta))
                .filter(|time| *time > 0.0 && *time < self.beta),
        );
        sort_and_deduplicate_times(&mut boundaries, self.beta);

        boundaries
            .windows(2)
            .map(|window| {
                let midpoint = 0.5 * (window[0] + window[1]);
                let product = self.spin_at(midpoint) * self.spin_at(midpoint + delta);
                f64::from(product) * (window[1] - window[0])
            })
            .sum::<f64>()
            / self.beta
    }

    /// Validate periodic worldline invariants.
    pub fn validate(&self) -> Result<(), ImpurityError> {
        if !self.beta.is_finite() || self.beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                format!("must be finite and positive, got {}", self.beta),
            ));
        }
        validate_spin(self.spin_at_zero)?;
        if !self.kinks.len().is_multiple_of(2) {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "periodic spin worldline requires an even kink count, got {}",
                self.kinks.len()
            )));
        }
        let mut previous = 0.0;
        for (index, &time) in self.kinks.iter().enumerate() {
            if !time.is_finite() || time <= 0.0 || time >= self.beta {
                return Err(ImpurityError::InvalidConfiguration(format!(
                    "kink {index} at {time} lies outside (0, beta)"
                )));
            }
            if index > 0 && time <= previous {
                return Err(ImpurityError::InvalidConfiguration(
                    "kink times must be strictly increasing".into(),
                ));
            }
            previous = time;
        }
        Ok(())
    }

    pub(crate) fn replace_from_segments(
        &mut self,
        cut_times: &[f64],
        segment_spins: &[i8],
    ) -> Result<(), ImpurityError> {
        if cut_times.is_empty() {
            if segment_spins.len() != 1 {
                return Err(ImpurityError::InvalidConfiguration(
                    "cut-free reconstruction requires one segment spin".into(),
                ));
            }
            validate_spin(segment_spins[0])?;
            self.spin_at_zero = segment_spins[0];
            self.kinks.clear();
            return Ok(());
        }
        if cut_times.len() != segment_spins.len() {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "{} cut times but {} segment spins",
                cut_times.len(),
                segment_spins.len()
            )));
        }

        for &spin in segment_spins {
            validate_spin(spin)?;
        }
        self.spin_at_zero = *segment_spins.last().ok_or_else(|| {
            ImpurityError::InvalidConfiguration("missing wrap-segment spin".into())
        })?;
        self.kinks.clear();
        for index in 0..segment_spins.len() {
            let next = (index + 1) % segment_spins.len();
            if segment_spins[index] != segment_spins[next] {
                self.kinks.push(cut_times[next]);
            }
        }
        self.kinks.retain(|time| *time > 0.0 && *time < self.beta);
        sort_and_deduplicate_times(&mut self.kinks, self.beta);
        self.validate()
    }
}

/// Build circular segments from all real and auxiliary cuts.
pub fn build_segments(
    worldline: &LongitudinalWorldline,
    cut_times: &[f64],
) -> Result<Vec<WorldlineSegment>, ImpurityError> {
    let mut cuts = cut_times.to_vec();
    sort_and_deduplicate_times(&mut cuts, worldline.beta());
    let mut segments = Vec::with_capacity(cuts.len().max(1));
    build_segments_from_sorted_into(worldline, &cuts, &mut segments)?;
    Ok(segments)
}

pub(crate) fn build_segments_from_sorted_into(
    worldline: &LongitudinalWorldline,
    cuts: &[f64],
    segments: &mut Vec<WorldlineSegment>,
) -> Result<(), ImpurityError> {
    segments.clear();
    if cuts.is_empty() {
        segments.push(WorldlineSegment::from_one(
            worldline.spin_at_zero(),
            TimeInterval::new(0.0, worldline.beta(), worldline.beta())?,
        )?);
        return Ok(());
    }

    if cuts
        .iter()
        .any(|time| *time <= 0.0 || *time >= worldline.beta())
    {
        return Err(ImpurityError::InvalidConfiguration(
            "auxiliary cut times must lie in (0, beta)".into(),
        ));
    }
    if cuts.windows(2).any(|window| window[0] >= window[1]) {
        return Err(ImpurityError::InvalidConfiguration(
            "auxiliary cut times must be strictly increasing".into(),
        ));
    }

    segments.reserve(cuts.len());
    for window in cuts.windows(2) {
        let interval = TimeInterval::new(window[0], window[1], worldline.beta())?;
        let midpoint = 0.5 * (window[0] + window[1]);
        segments.push(WorldlineSegment::from_one(
            worldline.spin_at(midpoint),
            interval,
        )?);
    }

    let first = cuts[0];
    let last = *cuts
        .last()
        .ok_or_else(|| ImpurityError::InvalidConfiguration("missing final auxiliary cut".into()))?;
    let upper = TimeInterval::new(last, worldline.beta(), worldline.beta())?;
    let lower = TimeInterval::new(0.0, first, worldline.beta())?;
    segments.push(WorldlineSegment::from_two(
        worldline.spin_at_zero(),
        upper,
        lower,
    )?);
    Ok(())
}

pub(crate) fn sort_and_deduplicate_times(times: &mut Vec<f64>, beta: f64) {
    times.sort_by(f64::total_cmp);
    let tolerance = TIME_EPSILON_FACTOR * f64::EPSILON * beta.max(1.0);
    times.dedup_by(|left, right| (*left - *right).abs() <= tolerance);
}

fn validate_spin(spin: i8) -> Result<(), ImpurityError> {
    if matches!(spin, -1 | 1) {
        Ok(())
    } else {
        Err(ImpurityError::InvalidConfiguration(format!(
            "spin must be -1 or +1, got {spin}"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn periodic_worldline_requires_even_kinks() {
        let error = LongitudinalWorldline::from_kinks(4.0, 1, vec![1.0])
            .expect_err("odd kink count must fail");
        assert!(matches!(error, ImpurityError::InvalidConfiguration(_)));
    }

    #[test]
    fn exact_integrated_magnetization_and_correlation() {
        let worldline =
            LongitudinalWorldline::from_kinks(4.0, 1, vec![1.0, 3.0]).expect("worldline");
        assert!(worldline.integrated_sigma_z().abs() < 1e-14);
        assert!((worldline.correlation_sigma_z(2.0) + 1.0).abs() < 1e-14);
        assert!((worldline.correlation_sigma_z(4.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn wrap_segment_contains_both_ends_of_time_circle() {
        let worldline =
            LongitudinalWorldline::from_kinks(4.0, 1, vec![1.0, 3.0]).expect("worldline");
        let segments = build_segments(&worldline, &[0.5, 1.0, 2.0, 3.0]).expect("segments");
        let wrap = segments.last().expect("wrap segment");
        assert_eq!(wrap.pieces().len(), 2);
        assert!((wrap.length() - 1.5).abs() < 1e-14);
        assert_eq!(wrap.spin(), 1);
    }

    #[test]
    fn reconstruction_removes_redundant_auxiliary_cuts() {
        let mut worldline = LongitudinalWorldline::new(4.0, 1).expect("worldline");
        worldline
            .replace_from_segments(&[1.0, 2.0, 3.0], &[1, -1, 1])
            .expect("reconstruct");
        assert_eq!(worldline.kinks(), &[2.0, 3.0]);
        assert_eq!(worldline.spin_at_zero(), 1);
    }
}
