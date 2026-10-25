//! C3: Weinrib-Halperin criterion map and the `nu(a)` flip (roadmap v2 §5,
//! stage C3).
//!
//! ## Anchors and the exact statement under test
//!
//! Pruessner et al., PRE 97, 062141 (2018), archived
//! `docs/ref/md/Pruessner_2018_correlated_thresholds.md`, Eq. (2)
//! (md lines ~106-112):
//!
//! ```text
//! nu_a = 2/a                 for a < a_H = 2/nu
//! nu_a = nu                  for a >= a_H = 2/nu,
//! ```
//!
//! i.e. the extended Harris criterion `a nu_iid = 2` separates a
//! long-range-correlated branch `nu = 2/a` from the uncorrelated value.
//! Their Table I (md lines ~596-624) lists, for the 2D square lattice,
//! `nu_a = 4/3` down to `a = 1.5`, then `8/5` at `a = 1.25`, `2` at
//! `a = 1`, `8/3` at `a = 0.75`, and `4` at `a = 0.5` — the flip sits
//! between `a = 1.5` and `a = 1.25`, consistent with
//! `a_H = 2/(4/3) = 1.5`. The criterion is "still slightly controversial"
//! in their words (md lines ~113-115), and the `2/a` value is the
//! Weinrib-Halperin conjecture (A. Weinrib, B. I. Halperin, *Critical
//! phenomena in systems with long-range-correlated quenched disorder*,
//! Phys. Rev. B 27, 413 (1983); formula confirmed by later work, e.g.
//! arXiv:2008.03169, whose abstract states "the conjecture of Weinrib and
//! Halperin: nu = 2/a").
//!
//! Schrenk et al., PRE 88, 052102 (2013), archived
//! `docs/ref/md/Schrenk_2013_correlated_disorder_percolation.md`
//! (md lines ~185-205) states the same criterion in Hurst language and
//! reports numerical support for `nu_H = -1/H`, which is `nu = 2/a` under
//! their `a = -2H` mapping (md line ~1422).
//!
//! ## Estimator (stated before use)
//!
//! The C3 validation uses [`ThresholdLadder::nu_from_pc_shift`] on a
//! [`ThresholdLadder::pool`] of independent ladders:
//!
//! ```text
//! p_c(inf) - p_c(L) = A L^{-1/nu},    p_c(inf) fixed
//! ```
//!
//! with `p_c(L)` the **order-parameter** susceptibility peak
//! `N Var(S_max/N)`, refined by a local log-quadratic fit of `ln chi_op`
//! over [`PEAK_FIT_RADIUS`](super::correlated_thresholds::PEAK_FIT_RADIUS).
//! The fixed `p_c(inf)` is the table value under
//! test (`0.592746` for the exact iid control, `0.54299/0.53895/0.53452`
//! for `a = 1.5/1.25/1.0`). Pooling eight independent ladders of 512
//! replicas each (4096 replicas per size, `L = 128/256/512/1024`) gives a
//! between-ladder peak error that the within-ladder bootstrap misses.
//!
//! **Control, 2026-10-25.** At `a = +infinity` the generator is exactly
//! iid site percolation and `p_c(inf)` is exact. The refined peaks are
//! `0.582975/0.586747/0.589137/0.590577` and the estimator returns
//! `nu = 1.379`, 95% CI `(1.325, 1.433)`, `chi2/dof = 0.11`: the interval
//! contains `4/3 = 1.333` with a `+3.4%` bias. The earlier
//! [`ThresholdLadder::nu_from_log_derivative`] returns `1.215`
//! `(0.739, 1.690)` on the same pooled ladder and `3.04` on a coarse
//! default window, so it is not a usable estimator; its maximum is
//! selected on a noisy flank with a one-grid-spacing difference. That
//! failure mode, not a physical branch shift, explains the earlier
//! `0.934/1.020/0.804` values.
//!
//! **Correlated rows, 2026-10-25.** The same pooled configuration gives
//!
//! | `a` | `p_c(inf)` table | refined peaks `L=128..1024` | `nu` (95% CI) | predicted |
//! |---|---|---|---|---|
//! | 1.5 | 0.54299 | 0.51481, 0.52226, 0.53057, 0.53499 | 1.592 (1.413, 1.771) | 4/3 |
//! | 1.25 | 0.53895 | 0.50223, 0.51571, 0.52093, 0.52754 | 1.819 (1.542, 2.095) | 8/5 |
//! | 1.0 | 0.53452 | 0.50230, 0.50922, 0.51448, 0.52117 | 2.529 (1.741, 3.318) | 2 |
//!
//! `chi2/dof = 1.42/1.47/0.58`. The ordering is correct and every value is
//! biased high, increasingly so as `a` decreases. The fit is
//! **prior-dominated**: `d nu / d p_c(inf) ~ 0.12` per `0.001`, so the
//! effective `p_c(inf)` would have to be only `0.002-0.003` below the
//! table (comparable to the finite-`L` zero-cutoff/reference-scale spread
//! already recorded by C2) for all three rows to agree with the
//! predictions. The free three-parameter fit `(p_c(inf), A, nu)` is
//! ill-conditioned at these sizes: on the iid control it returns
//! `p_c(inf) = 0.5932` and `nu = 1.43` instead of the exact
//! `0.592746/1.333`, and on the correlated rows it runs to the scan
//! boundary. The gap `nu(1.25) - nu(1.5) = 0.226 +/- 0.208` (`1.1 sigma`)
//! does **not** separate the branches.
//!
//! **Status: not certified.** C3 keeps the analytic map and the calibrated
//! control, but the finite-size `nu` flip is an explicit negative result:
//! `L <= 1024` with a fixed literature `p_c(inf)` cannot distinguish the
//! branches at 95%. What would be needed is either `L >= 2048-8192`
//! (Pruessner used `8192` with `10^5` replicas) to suppress the
//! zero-cutoff/finite-size corrections, or an independent `p_c(inf)` for
//! this finite-`L` FFM implementation at the `0.001` level (a
//! Binder-cumulant or very-large-`L` crossing), or both.
//!
//! [`HarrisBracket`] carries the fit's reduced chi-square so poor scaling
//! fits are visible. A bracket consistent with `4/3` above the criterion
//! and with `2/a` below it, plus a positive gap, would be evidence for the
//! flip at the ladder's resolution; it would not be a proof of the
//! conjecture.
//!
//! ## Harris criterion product
//!
//! [`WeinribHalperinMap::criterion_product`] returns `a nu_iid`; the
//! correlated branch is selected exactly when it is below `2`. For the 2D
//! square lattice `nu_iid = 4/3`, hence `a_H = 3/2`.

use core::fmt;

use super::correlated_thresholds::{
    CorrelatedThresholdError, IntervalEstimate, NuEstimate, ThresholdLadder,
};

/// Failure to build or evaluate a Weinrib-Halperin bracket.
#[derive(Clone, Debug, PartialEq)]
pub enum WeinribHalperinError {
    /// Measured exponents come from a ladder or fit that failed.
    Fit(CorrelatedThresholdError),
    /// The correlation exponent must be finite and strictly positive.
    Exponent { exponent: f64 },
    /// The uncorrelated correlation-length exponent must be finite and
    /// positive.
    NuIid { nu_iid: f64 },
    /// A bracket needs two ladders with different exponents.
    SameExponent,
}

impl fmt::Display for WeinribHalperinError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(formatter),
            Self::Exponent { exponent } => write!(
                formatter,
                "Weinrib-Halperin map needs a finite positive correlation exponent, got {exponent}"
            ),
            Self::NuIid { nu_iid } => write!(
                formatter,
                "uncorrelated correlation-length exponent must be finite and positive, got {nu_iid}"
            ),
            Self::SameExponent => write!(
                formatter,
                "a branch flip bracket needs ladders with different exponents"
            ),
        }
    }
}

impl std::error::Error for WeinribHalperinError {}

/// The extended Harris criterion `a nu_iid = 2` and its two branches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WeinribHalperinMap {
    exponent: f64,
    nu_iid: f64,
}

impl WeinribHalperinMap {
    /// The 2D square-lattice map with the exact uncorrelated value
    /// `nu_iid = 4/3`.
    pub fn square(exponent: f64) -> Result<Self, WeinribHalperinError> {
        Self::new(exponent, 4.0 / 3.0)
    }

    /// Map with a caller-supplied uncorrelated exponent (e.g. `0.875` in
    /// `d = 3`); no other dimension is assumed by this module.
    pub fn new(exponent: f64, nu_iid: f64) -> Result<Self, WeinribHalperinError> {
        if !(exponent.is_finite() && exponent > 0.0) {
            return Err(WeinribHalperinError::Exponent { exponent });
        }
        if !(nu_iid.is_finite() && nu_iid > 0.0) {
            return Err(WeinribHalperinError::NuIid { nu_iid });
        }
        Ok(Self { exponent, nu_iid })
    }

    /// The correlation exponent `a`.
    #[inline]
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }

    /// The uncorrelated correlation-length exponent.
    #[inline]
    pub const fn nu_iid(&self) -> f64 {
        self.nu_iid
    }

    /// The criterion product `a nu_iid`; the correlated branch applies
    /// strictly below `2`.
    #[inline]
    pub fn criterion_product(&self) -> f64 {
        self.exponent * self.nu_iid
    }

    /// The Harris threshold `a_H = 2/nu_iid`.
    #[inline]
    pub fn threshold(&self) -> f64 {
        2.0 / self.nu_iid
    }

    /// Whether the correlated branch `nu = 2/a` is predicted.
    #[inline]
    pub fn is_correlated_branch(&self) -> bool {
        self.criterion_product() < 2.0
    }

    /// The predicted `nu(a)`: `2/a` below the criterion, `nu_iid` above it.
    #[inline]
    pub fn predicted_nu(&self) -> f64 {
        if self.is_correlated_branch() {
            2.0 / self.exponent
        } else {
            self.nu_iid
        }
    }

    /// The predicted branch name (`"2/a"` or `"nu_iid"`), for reports.
    #[inline]
    pub fn branch(&self) -> &'static str {
        if self.is_correlated_branch() {
            "2/a"
        } else {
            "nu_iid"
        }
    }
}

/// Measured `nu(a)` against one predicted branch.
#[derive(Clone, Debug, PartialEq)]
pub struct HarrisBracket {
    /// The map whose prediction is compared.
    pub map: WeinribHalperinMap,
    /// The measured exponent with its bootstrap interval.
    pub measured: NuEstimate,
    /// The predicted value.
    pub predicted: f64,
    /// `(measured - predicted) / measured.standard_error` (the systematic
    /// fit-window error is not included; inspect `reduced_chi_square`).
    pub z_score: f64,
    /// Whether `|z_score| <= 2` (a bracket, never a proof).
    pub consistent: bool,
}

impl HarrisBracket {
    /// The measured minus predicted gap with a 95% interval.
    pub fn gap(&self) -> IntervalEstimate {
        IntervalEstimate {
            value: self.measured.nu.value - self.predicted,
            lower: self.measured.nu.lower - self.predicted,
            upper: self.measured.nu.upper - self.predicted,
            standard_error: self.measured.nu.standard_error,
        }
    }
}

/// Estimate `nu(a)` from a C2 ladder and bracket it against the
/// Weinrib-Halperin prediction.
pub fn bracket(
    ladder: &ThresholdLadder,
    minimum_side: usize,
    nu_iid: f64,
) -> Result<HarrisBracket, WeinribHalperinError> {
    let map = WeinribHalperinMap::new(ladder.exponent(), nu_iid)?;
    let measured = ladder
        .nu_from_log_derivative(minimum_side)
        .map_err(WeinribHalperinError::Fit)?;
    let predicted = map.predicted_nu();
    let error = measured.nu.standard_error.max(1e-12);
    let z_score = (measured.nu.value - predicted) / error;
    Ok(HarrisBracket {
        map,
        consistent: z_score.abs() <= 2.0,
        predicted,
        z_score,
        measured,
    })
}

/// The flip between the uncorrelated branch (`a >= a_H`) and the `2/a`
/// branch (`a < a_H`), from two independent ladders.
#[derive(Clone, Debug, PartialEq)]
pub struct FlipBracket {
    /// The ladder above the Harris threshold (predicted `nu_iid`).
    pub above: HarrisBracket,
    /// The ladder below the Harris threshold (predicted `2/a`).
    pub below: HarrisBracket,
    /// Measured `nu(below) - nu(above)` with a 95% interval.
    pub gap: IntervalEstimate,
    /// Whether the gap's lower bound is positive (the branches are ordered
    /// as predicted at this ladder's resolution).
    pub separated: bool,
}

/// Compare two ladders across the Harris threshold.
///
/// The gap's interval assumes the two ladders use independent RNG streams
/// (they do when built from independent specs); it is the sum of the two
/// bootstrap variances. No correlation between ladders is assumed.
pub fn flip_bracket(
    above: &ThresholdLadder,
    below: &ThresholdLadder,
    minimum_side: usize,
    nu_iid: f64,
) -> Result<FlipBracket, WeinribHalperinError> {
    if above.exponent() == below.exponent() {
        return Err(WeinribHalperinError::SameExponent);
    }
    let above_bracket = bracket(above, minimum_side, nu_iid)?;
    let below_bracket = bracket(below, minimum_side, nu_iid)?;
    let value = below_bracket.measured.nu.value - above_bracket.measured.nu.value;
    let variance = below_bracket.measured.nu.standard_error.powi(2)
        + above_bracket.measured.nu.standard_error.powi(2);
    let error = variance.sqrt();
    let gap = IntervalEstimate {
        value,
        lower: value - 1.96 * error,
        upper: value + 1.96 * error,
        standard_error: error,
    };
    Ok(FlipBracket {
        separated: gap.lower > 0.0,
        above: above_bracket,
        below: below_bracket,
        gap,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_map_reproduces_the_table_branches() {
        let cases = [
            (3.0, 4.0 / 3.0),
            (2.5, 4.0 / 3.0),
            (2.0, 4.0 / 3.0),
            (1.75, 4.0 / 3.0),
            (1.5, 4.0 / 3.0),
            (1.25, 8.0 / 5.0),
            (1.0, 2.0),
            (0.75, 8.0 / 3.0),
            (0.5, 4.0),
        ];
        for (exponent, expected) in cases {
            let map = WeinribHalperinMap::square(exponent).unwrap();
            assert!(
                (map.predicted_nu() - expected).abs() < 1e-12,
                "a={exponent}: {} vs {expected}",
                map.predicted_nu()
            );
        }
        assert!((WeinribHalperinMap::square(1.5).unwrap().threshold() - 1.5).abs() < 1e-12);
        assert_eq!(
            WeinribHalperinMap::square(1.5).unwrap().criterion_product(),
            2.0
        );
        assert!(!WeinribHalperinMap::square(1.5)
            .unwrap()
            .is_correlated_branch());
        assert!(WeinribHalperinMap::square(1.499)
            .unwrap()
            .is_correlated_branch());
        assert_eq!(WeinribHalperinMap::square(0.5).unwrap().branch(), "2/a");
    }

    #[test]
    fn invalid_inputs_are_typed() {
        assert!(matches!(
            WeinribHalperinMap::new(0.0, 4.0 / 3.0),
            Err(WeinribHalperinError::Exponent { .. })
        ));
        assert!(matches!(
            WeinribHalperinMap::new(f64::INFINITY, 4.0 / 3.0),
            Err(WeinribHalperinError::Exponent { .. })
        ));
        assert!(matches!(
            WeinribHalperinMap::new(1.0, 0.0),
            Err(WeinribHalperinError::NuIid { .. })
        ));
    }
}
