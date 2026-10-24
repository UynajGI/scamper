//! L3: box-diameter transitions across Biskup's long-range regimes
//! (roadmap v2 §4 L3), over the L2 chemical-distance layer.
//!
//! The observable is the exact finite-box graph diameter `D_N` — the
//! maximum finite chemical distance over all vertex pairs of one L1
//! realization (isolates contribute 0; see
//! [`ChemicalDistanceOutcome::diameter`]). This module is a thin,
//! allocation-controlled analysis layer: every realization is produced
//! by [`ChemicalDistanceSampler::realize_with_diameter`], so the L1 draw
//! contract, the CSR build, the BFS machinery, and the
//! [`DIAMETER_VERTEX_CAP`] guard-rail are exactly L2's — no second
//! diameter implementation exists.
//!
//! # Regime semantics (anchor, regime-correct)
//!
//! Anchor: Biskup, *On the scaling of the chemical distance in
//! long-range percolation*, Ann. Probab. 32, 2938 (2004), archived at
//! `docs/ref/md/Biskup_2003_long_range_chemical_distance.md`. In `d = 1`
//! the boundaries fall at `s = 1` and `s = 2`; the statements that concern
//! the **diameter** are:
//!
//! - `s < 1`: Eq. (1.8) gives the a.s. SUP statement
//!   `sup_{x,y in Z^d} D(x,y) = d/(d-s)` — a constant; the box diameter
//!   is bounded (finite boxes may sit above the infinite-lattice value
//!   because restricting to the box removes shortcut vertices).
//! - `s = 1`: Eq. (1.9) is the diameter statement
//!   `D_N = Theta(log N / log log N)` (Biskup attributes it to
//!   Coppersmith-Gamarnik-Sviridenko); the typical distance is only
//!   "at most logarithmic" with no lower bound.
//! - `1 < s < 2`: the sharp statement is for the TYPICAL distance; for
//!   the diameter Biskup explicitly says the available bounds are not
//!   sharp enough for a definitive asymptotic statement. This layer
//!   therefore brackets only sub-polynomial growth, never a sharp
//!   exponent.
//! - `s > 2`: Biskup expects `D_N` linear in `N` (proved in `d = 1` by
//!   Benjamini-Berger); gate the approach of `D_N/N` to a positive
//!   constant from above.
//! - `s = 2`: explicitly not well understood — [`DiameterRegime::Critical`]
//!   and `supports_diameter_gate() == false`; no gated claim.
//!
//! These are finite-box statistical brackets of infinite-lattice
//! asymptotics. The validation gates live in
//! `tests/physics/long_range_diameter.rs`; see PERCOLATION.md for the
//! recorded numbers and finite-size caveats.
//!
//! # Draw contract and cost
//!
//! A ladder is measured row by row (`s` outer, `N` inner). Row `r`
//! draws its realization `i` from `Xoshiro256PlusPlus::seed_from_u64(seed
//! + r * realizations + i)` — disjoint per-row streams, deterministic and
//! documented so a caller can replay one row exactly. Each realization costs the L2 replay (`O(N^2)` pair decisions)
//! plus `N` BFS passes (`O(N (N + m))`).
//!
//! The sampler reserves its buffers at construction; the ladder itself
//! only allocates the recorded results.

use core::fmt;

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::long_range::{LongRangeBox, LongRangeError};
use super::long_range_distances::{
    ChemicalDistanceError, ChemicalDistanceSampler, DIAMETER_VERTEX_CAP,
};

/// Biskup regime of the box diameter in `d = 1`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DiameterRegime {
    /// `s < 1`: bounded diameter (Eq. 1.8 sup statement).
    Bounded,
    /// `s = 1`: `D_N = Theta(log N / log log N)` (Eq. 1.9).
    Logarithmic,
    /// `1 < s < 2`: no sharp diameter law is proved; only
    /// sub-polynomial growth is bracketed.
    Polylogarithmic,
    /// `s = 2`: the critical case, explicitly not well understood.
    Critical,
    /// `s > 2`: linear diameter with a positive constant.
    Linear,
}

impl DiameterRegime {
    /// Classify an exponent for `d = 1`; `None` for non-finite or
    /// negative `s`.
    pub fn classify(s: f64) -> Option<Self> {
        if !s.is_finite() || s < 0.0 {
            return None;
        }
        if s < 1.0 {
            Some(Self::Bounded)
        } else if s == 1.0 {
            Some(Self::Logarithmic)
        } else if s < 2.0 {
            Some(Self::Polylogarithmic)
        } else if s == 2.0 {
            Some(Self::Critical)
        } else {
            Some(Self::Linear)
        }
    }

    /// Whether stage L3 makes a gated diameter claim in this regime.
    /// The critical case `s = 2` deliberately does not.
    pub const fn supports_diameter_gate(self) -> bool {
        !matches!(self, Self::Critical)
    }

    /// Short stable label for reports.
    pub const fn label(self) -> &'static str {
        match self {
            Self::Bounded => "bounded",
            Self::Logarithmic => "logarithmic",
            Self::Polylogarithmic => "polylogarithmic",
            Self::Critical => "critical-no-gate",
            Self::Linear => "linear",
        }
    }
}

/// One ladder row: the exact diameters of `samples` independent
/// realizations of one `(s, N)` box, with their sample mean and SEM.
#[derive(Clone, Debug, PartialEq)]
pub struct DiameterRow {
    exponent: f64,
    size: usize,
    diameters: Vec<usize>,
    mean: f64,
    standard_error: f64,
}

impl DiameterRow {
    /// The decay exponent `s`.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The 1D box side length `N`.
    #[inline]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// The regime this row belongs to (`None` for an invalid exponent,
    /// impossible after construction).
    pub fn regime(&self) -> Option<DiameterRegime> {
        DiameterRegime::classify(self.exponent)
    }

    /// The exact diameters, one per realization in draw order.
    #[inline]
    pub fn diameters(&self) -> &[usize] {
        &self.diameters
    }

    /// Number of realizations.
    #[inline]
    pub fn samples(&self) -> usize {
        self.diameters.len()
    }

    /// Sample mean of the diameters.
    #[inline]
    pub const fn mean(&self) -> f64 {
        self.mean
    }

    /// Standard error of [`Self::mean`].
    #[inline]
    pub const fn standard_error(&self) -> f64 {
        self.standard_error
    }

    /// `mean / N` — the linear-regime ratio.
    #[inline]
    pub fn ratio(&self) -> f64 {
        self.mean / self.size as f64
    }

    /// Standard error of [`Self::ratio`].
    #[inline]
    pub fn ratio_standard_error(&self) -> f64 {
        self.standard_error / self.size as f64
    }
}

/// A complete `(exponent, size)` ladder of exact box-diameter samples.
#[derive(Clone, Debug, PartialEq)]
pub struct DiameterLadder {
    beta: f64,
    rows: Vec<DiameterRow>,
}

impl DiameterLadder {
    /// The intensity every row was measured at.
    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// Rows in measurement order (`s` outer, `N` inner).
    #[inline]
    pub fn rows(&self) -> &[DiameterRow] {
        &self.rows
    }

    /// The row for `(exponent, size)`, if measured.
    pub fn row(&self, exponent: f64, size: usize) -> Option<&DiameterRow> {
        self.rows
            .iter()
            .find(|row| row.exponent == exponent && row.size == size)
    }

    /// All rows for one exponent, in ladder order.
    pub fn rows_for(&self, exponent: f64) -> impl Iterator<Item = &DiameterRow> {
        self.rows.iter().filter(move |row| row.exponent == exponent)
    }
}

/// Measure exact box diameters over a size/exponent ladder on 1D boxes.
///
/// Every realization reuses [`ChemicalDistanceSampler`] (one source at
/// `N / 2`, exact all-vertex diameter). The seed contract is documented
/// on the module: row `r` uses `seed + r + i` for its `i`-th realization.
/// `realizations >= 2` is required so the recorded SEM is defined.
pub fn measure_diameter_ladder(
    sizes: &[usize],
    exponents: &[f64],
    beta: f64,
    realizations: u64,
    seed: u64,
) -> Result<DiameterLadder, DiameterLadderError> {
    if sizes.is_empty() {
        return Err(DiameterLadderError::EmptySizes);
    }
    if exponents.is_empty() {
        return Err(DiameterLadderError::EmptyExponents);
    }
    if realizations < 2 {
        return Err(DiameterLadderError::TooFewRealizations { realizations });
    }
    for &size in sizes {
        if size == 0 {
            return Err(DiameterLadderError::ZeroSize);
        }
        if size > DIAMETER_VERTEX_CAP {
            return Err(DiameterLadderError::SizeAboveCap {
                size,
                cap: DIAMETER_VERTEX_CAP,
            });
        }
    }
    for &exponent in exponents {
        if DiameterRegime::classify(exponent).is_none() {
            return Err(DiameterLadderError::ExponentInvalid { value: exponent });
        }
    }

    let realization_count = usize::try_from(realizations)
        .map_err(|_| DiameterLadderError::TooManyRealizations { realizations })?;
    let row_count = sizes
        .len()
        .checked_mul(exponents.len())
        .ok_or(DiameterLadderError::TooManyRows)?;
    let mut rows = Vec::new();
    rows.try_reserve_exact(row_count)
        .map_err(|_| DiameterLadderError::Capacity)?;
    let mut row_index = 0_u64;
    for &exponent in exponents {
        for &size in sizes {
            let box_ = LongRangeBox::new(&[size])?;
            let mut sampler = ChemicalDistanceSampler::new(box_, beta, exponent, &[size / 2])?;
            let mut diameters = Vec::new();
            diameters
                .try_reserve_exact(realization_count)
                .map_err(|_| DiameterLadderError::Capacity)?;
            // Disjoint per-row streams: row r starts realization_count
            // draws after row r-1, so no two rows share a seed.
            let base = seed.wrapping_add(row_index.wrapping_mul(realizations));
            for index in 0..realizations {
                let outcome = sampler.realize_with_diameter(
                    &mut Xoshiro256PlusPlus::seed_from_u64(base.wrapping_add(index)),
                )?;
                diameters.push(
                    outcome
                        .diameter
                        .expect("diameter tracking was requested by realize_with_diameter"),
                );
            }
            let (mean, standard_error) = sample_mean_sem(&diameters)?;
            rows.push(DiameterRow {
                exponent,
                size,
                diameters,
                mean,
                standard_error,
            });
            row_index += 1;
        }
    }
    Ok(DiameterLadder { beta, rows })
}

/// Sample mean and SEM of exact diameters (kept private; the diagnostic
/// is the mean/SEM pair on [`DiameterRow`]).
fn sample_mean_sem(diameters: &[usize]) -> Result<(f64, f64), DiameterLadderError> {
    if diameters.len() < 2 {
        return Err(DiameterLadderError::TooFewRealizations {
            realizations: diameters.len() as u64,
        });
    }
    let count = diameters.len() as f64;
    let mean = diameters.iter().map(|&value| value as f64).sum::<f64>() / count;
    let variance = diameters
        .iter()
        .map(|&value| {
            let delta = value as f64 - mean;
            delta * delta
        })
        .sum::<f64>()
        / (count - 1.0);
    Ok((mean, (variance / count).sqrt()))
}

/// Typed rejections of the L3 ladder.
#[derive(Debug)]
pub enum DiameterLadderError {
    /// Empty size ladder.
    EmptySizes,
    /// Empty exponent ladder.
    EmptyExponents,
    /// A zero side length.
    ZeroSize,
    /// A side length above L2's all-vertex diameter cap.
    SizeAboveCap { size: usize, cap: usize },
    /// Fewer than two realizations (SEM undefined).
    TooFewRealizations { realizations: u64 },
    /// The realization count does not fit `usize`.
    TooManyRealizations { realizations: u64 },
    /// The requested size x exponent row count overflows.
    TooManyRows,
    /// Result storage could not be reserved.
    Capacity,
    /// A non-finite or negative exponent.
    ExponentInvalid { value: f64 },
    /// A box-geometry rejection.
    Geometry(LongRangeError),
    /// An L2 chemical-distance rejection.
    Distance(ChemicalDistanceError),
}

impl From<LongRangeError> for DiameterLadderError {
    fn from(error: LongRangeError) -> Self {
        Self::Geometry(error)
    }
}

impl From<ChemicalDistanceError> for DiameterLadderError {
    fn from(error: ChemicalDistanceError) -> Self {
        Self::Distance(error)
    }
}

impl fmt::Display for DiameterLadderError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySizes => write!(formatter, "the size ladder must not be empty"),
            Self::EmptyExponents => write!(formatter, "the exponent ladder must not be empty"),
            Self::ZeroSize => write!(formatter, "box side lengths must be positive"),
            Self::SizeAboveCap { size, cap } => write!(
                formatter,
                "box side {size} exceeds the all-vertex diameter cap {cap}"
            ),
            Self::TooFewRealizations { realizations } => write!(
                formatter,
                "at least two realizations are required for a diameter SEM, got {realizations}"
            ),
            Self::TooManyRealizations { realizations } => write!(
                formatter,
                "realization count {realizations} does not fit the platform usize"
            ),
            Self::TooManyRows => write!(
                formatter,
                "the size x exponent ladder row count overflows usize"
            ),
            Self::Capacity => write!(formatter, "result storage could not be reserved"),
            Self::ExponentInvalid { value } => write!(
                formatter,
                "exponent {value} is negative or not finite (zero is legal but degenerate)"
            ),
            Self::Geometry(error) => write!(formatter, "box geometry rejected: {error}"),
            Self::Distance(error) => write!(formatter, "chemical-distance layer rejected: {error}"),
        }
    }
}

impl std::error::Error for DiameterLadderError {}
