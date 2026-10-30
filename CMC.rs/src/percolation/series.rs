//! Fixed-shape C0 series buffer with explicit padding conventions.
//!
//! Carlo.rs records vector observables as fixed-length arrays: every sample
//! of an observable name must have the same length, so a realization whose
//! natural trajectory is shorter (or longer) than the recorded shape is
//! compressed into that shape. The process adapters (invasion, bootstrap,
//! k-core, long-range) all did this by hand and each documented the same two
//! conventions:
//!
//! * **repeat-last**: entries past the realized trajectory repeat its final
//!   value (an empty trajectory uses a caller-supplied fallback);
//! * **fixed tail**: entries past the realized trajectory carry a constant
//!   final summary (for example the final core size).
//!
//! Realizations longer than the buffer are truncated; the exact scalar
//! summaries remain available alongside the compressed series.
//!
//! [`FixedSeries`] owns the buffer once and rewrites it in place, so the
//! per-sweep path does not allocate.

use std::collections::TryReserveError;

/// How [`FixedSeries::fill`] pads entries past the realized trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SeriesPad {
    /// Repeat the last realized entry; use `fallback` when nothing was
    /// realized.
    RepeatLast { fallback: f64 },
    /// Pad every remaining entry with this constant.
    Value(f64),
}

/// A fixed-length `f64` buffer written from a variable-length realization.
///
/// ```
/// use cmc_rs::{FixedSeries, SeriesPad};
///
/// let mut series = FixedSeries::new(4);
/// series.fill([1.0, 2.0], SeriesPad::RepeatLast { fallback: 0.0 });
/// assert_eq!(series.as_slice(), &[1.0, 2.0, 2.0, 2.0]);
///
/// // Longer realizations truncate to the recorded shape.
/// series.fill([5.0, 6.0, 7.0, 8.0, 9.0], SeriesPad::Value(0.0));
/// assert_eq!(series.as_slice(), &[5.0, 6.0, 7.0, 8.0]);
/// ```
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FixedSeries {
    values: Vec<f64>,
}

impl FixedSeries {
    /// A zero-filled series of exactly `len` entries.
    pub fn new(len: usize) -> Self {
        Self {
            values: vec![0.0; len],
        }
    }

    /// A zero-filled series, surfacing reservation failure instead of
    /// aborting on allocation failure.
    pub fn try_new(len: usize) -> Result<Self, TryReserveError> {
        let mut values = Vec::new();
        values.try_reserve_exact(len)?;
        values.resize(len, 0.0);
        Ok(Self { values })
    }

    /// Recorded shape.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the recorded shape is zero-length.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// The series as recorded.
    pub fn as_slice(&self) -> &[f64] {
        &self.values
    }

    /// Rewrite the series from one realization, applying `pad` after the
    /// realized entries and truncating entries beyond the recorded shape.
    ///
    /// Only the first [`Self::len`] realized values are consumed, so a long
    /// iterator is not walked past the end.
    pub fn fill<I>(&mut self, realized: I, pad: SeriesPad)
    where
        I: IntoIterator<Item = f64>,
    {
        let mut realized = realized.into_iter();
        let mut written = 0;
        let mut last = match pad {
            SeriesPad::RepeatLast { fallback } => fallback,
            SeriesPad::Value(value) => value,
        };
        while written < self.values.len() {
            match realized.next() {
                Some(value) => {
                    self.values[written] = value;
                    last = value;
                    written += 1;
                }
                None => break,
            }
        }
        let tail = match pad {
            SeriesPad::RepeatLast { .. } => last,
            SeriesPad::Value(value) => value,
        };
        for slot in &mut self.values[written..] {
            *slot = tail;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeat_last_pads_with_the_final_realized_value() {
        let mut series = FixedSeries::new(5);
        series.fill([1.0, 2.0, 3.0], SeriesPad::RepeatLast { fallback: 9.0 });
        assert_eq!(series.len(), 5);
        assert_eq!(series.as_slice(), &[1.0, 2.0, 3.0, 3.0, 3.0]);
    }

    #[test]
    fn repeat_last_uses_the_fallback_for_an_empty_realization() {
        let mut series = FixedSeries::new(3);
        series.fill([], SeriesPad::RepeatLast { fallback: 0.5 });
        assert_eq!(series.as_slice(), &[0.5, 0.5, 0.5]);
    }

    #[test]
    fn fixed_tail_pads_with_its_constant_even_after_values() {
        let mut series = FixedSeries::new(4);
        series.fill([1.0, 2.0], SeriesPad::Value(7.0));
        assert_eq!(series.as_slice(), &[1.0, 2.0, 7.0, 7.0]);
    }

    #[test]
    fn long_realizations_truncate_and_short_writes_restore_the_shape() {
        let mut series = FixedSeries::new(2);
        series.fill([1.0, 2.0, 3.0], SeriesPad::Value(-1.0));
        assert_eq!(series.as_slice(), &[1.0, 2.0]);
        series.fill([4.0], SeriesPad::Value(-1.0));
        assert_eq!(series.as_slice(), &[4.0, -1.0]);
    }

    #[test]
    fn zero_length_series_accepts_every_pad() {
        let mut series = FixedSeries::new(0);
        assert!(series.is_empty());
        series.fill([1.0, 2.0], SeriesPad::RepeatLast { fallback: 0.0 });
        series.fill([], SeriesPad::Value(3.0));
        assert!(series.as_slice().is_empty());
    }
}
