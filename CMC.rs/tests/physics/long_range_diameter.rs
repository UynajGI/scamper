//! L3 validation: box-diameter transitions across Biskup's regimes.
//!
//! Gates:
//! 1. regime classification (including the ungated critical case s = 2);
//! 2. ladder replay parity against direct L2 diameter draws at the same
//!    documented seeds;
//! 3. s > 2: D_N/N approaches a positive constant from above (aggregate
//!    95% CI, log-log slope consistent with 1, extrapolated intercept in
//!    (0, 1));
//! 4. 1 < s < 2: the diameter grows sub-polynomially at these sizes
//!    (log-log slope CI bounded below 0.5) — never a sharp exponent;
//! 5. s = 1: D_N is bracketed by the Eq. (1.9) log N / log log N scale
//!    (upper envelope asserted, ratio recorded);
//! 6. s < 1: bounded diameter, no positive growth;
//! 7. s = 2: no run, no gate (documented).
//!
//! Every number is a finite-box statistical bracket of the archived
//! Biskup asymptotics; see the module documentation.

use cmc_rs::{
    measure_diameter_ladder, ChemicalDistanceSampler, DiameterLadderError, DiameterRegime,
    LongRangeBox,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Ordinary least squares of `y ~ x`: (slope, intercept, slope SE, R^2).
fn ols(points: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|point| point.0).sum();
    let sy: f64 = points.iter().map(|point| point.1).sum();
    let sxx: f64 = points.iter().map(|point| point.0 * point.0).sum();
    let sxy: f64 = points.iter().map(|point| point.0 * point.1).sum();
    let syy: f64 = points.iter().map(|point| point.1 * point.1).sum();
    let denominator = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denominator;
    let intercept = (sy - slope * sx) / n;
    let r_squared = (n * sxy - sx * sy).powi(2) / (denominator * (n * syy - sy * sy));
    let residuals: f64 = points
        .iter()
        .map(|point| {
            let error = point.1 - (intercept + slope * point.0);
            error * error
        })
        .sum();
    let slope_se = (residuals / (n - 2.0) / (sxx - sx * sx / n)).sqrt();
    (slope, intercept, slope_se, r_squared)
}

#[test]
fn regime_classification_and_critical_case_is_ungated() {
    assert_eq!(DiameterRegime::classify(0.5), Some(DiameterRegime::Bounded));
    assert_eq!(
        DiameterRegime::classify(1.0),
        Some(DiameterRegime::Logarithmic)
    );
    assert_eq!(
        DiameterRegime::classify(1.5),
        Some(DiameterRegime::Polylogarithmic)
    );
    assert_eq!(
        DiameterRegime::classify(2.0),
        Some(DiameterRegime::Critical)
    );
    assert_eq!(DiameterRegime::classify(2.5), Some(DiameterRegime::Linear));
    assert_eq!(DiameterRegime::classify(f64::NAN), None);
    assert_eq!(DiameterRegime::classify(-1.0), None);
    assert!(!DiameterRegime::Critical.supports_diameter_gate());
    for regime in [
        DiameterRegime::Bounded,
        DiameterRegime::Logarithmic,
        DiameterRegime::Polylogarithmic,
        DiameterRegime::Linear,
    ] {
        assert!(regime.supports_diameter_gate(), "{}", regime.label());
    }
    eprintln!("L3 regimes: s<1 bounded, s=1 logarithmic, (1,2) polylog, s>2 linear; s=2 no gate");
}

#[test]
fn ladder_replays_direct_l2_draws_bitwise() {
    let sizes = [12_usize, 16];
    let exponents = [1.5_f64];
    let (beta, realizations, seed) = (1.0_f64, 5_u64, 0x4C33_1001_u64);
    let ladder =
        measure_diameter_ladder(&sizes, &exponents, beta, realizations, seed).expect("ladder");
    assert_eq!(ladder.rows().len(), 2);
    for (row_index, &size) in sizes.iter().enumerate() {
        let row = ladder.row(1.5, size).expect("row");
        assert_eq!(row.samples(), realizations as usize);
        assert_eq!(row.regime(), Some(DiameterRegime::Polylogarithmic));
        let box_ = LongRangeBox::new(&[size]).expect("box");
        let mut sampler =
            ChemicalDistanceSampler::new(box_, beta, 1.5, &[size / 2]).expect("sampler");
        for (index, &recorded) in row.diameters().iter().enumerate() {
            let outcome = sampler
                .realize_with_diameter(&mut Xoshiro256PlusPlus::seed_from_u64(
                    seed.wrapping_add(row_index as u64)
                        .wrapping_add(index as u64),
                ))
                .expect("direct draw");
            assert_eq!(
                outcome.diameter,
                Some(recorded),
                "row {row_index} realization {index}"
            );
        }
    }
    eprintln!("L3 ladder replay: 10 diameters match direct L2 draws bitwise");
}

#[test]
fn linear_regime_ratio_approaches_a_positive_constant() {
    const SIZES: [usize; 5] = [128, 256, 512, 1024, 2048];
    const REALIZATIONS: u64 = 24;
    // t_{0.975} on df = 3 (five points minus two parameters).
    const T975_DF3: f64 = 3.1824;
    let ladder =
        measure_diameter_ladder(&SIZES, &[2.5], 1.0, REALIZATIONS, 0x4C33_2001).expect("ladder");
    let rows: Vec<_> = ladder.rows_for(2.5).collect();
    assert_eq!(rows.len(), SIZES.len());
    for row in &rows {
        assert_eq!(row.regime(), Some(DiameterRegime::Linear));
        assert!(
            row.diameters()
                .iter()
                .all(|&diameter| diameter < row.size()),
            "diameters stay below N"
        );
    }
    let ratios: Vec<f64> = rows.iter().map(|row| row.ratio()).collect();
    let ratio_sems: Vec<f64> = rows.iter().map(|row| row.ratio_standard_error()).collect();
    eprintln!(
        "L3 linear s=2.5: means={:?} ratios={:?}",
        rows.iter().map(|row| row.mean()).collect::<Vec<_>>(),
        ratios
    );
    for index in 0..ratios.len() - 1 {
        let margin = 2.0 * (ratio_sems[index].powi(2) + ratio_sems[index + 1].powi(2)).sqrt();
        assert!(
            ratios[index + 1] <= ratios[index] + margin,
            "D_N/N must not increase: {:.4} -> {:.4} (margin {margin:.4})",
            ratios[index],
            ratios[index + 1]
        );
    }
    assert!(
        ratios[0] > ratios[ratios.len() - 1],
        "the linear-regime ratio must decrease toward the time constant"
    );

    let log_points: Vec<(f64, f64)> = rows
        .iter()
        .map(|row| ((row.size() as f64).ln(), row.mean().ln()))
        .collect();
    let (slope, _, slope_se, r_squared) = ols(&log_points);
    let (low, high) = (slope - T975_DF3 * slope_se, slope + T975_DF3 * slope_se);
    assert!(
        low > 0.8,
        "log-log slope CI ({low:.3}, {high:.3}) excludes the polylog regime"
    );
    assert!(
        high < 1.05,
        "log-log slope CI ({low:.3}, {high:.3}) stays consistent with linear growth"
    );
    eprintln!(
        "L3 linear fit: slope={slope:.3}±{slope_se:.3} CI=({low:.3}, {high:.3}) \
         R2={r_squared:.4}"
    );

    let ratio_points: Vec<(f64, f64)> = rows
        .iter()
        .map(|row| (1.0 / row.size() as f64, row.ratio()))
        .collect();
    let (deficit_slope, time_constant, deficit_se, _) = ols(&ratio_points);
    let n = ratio_points.len() as f64;
    let mean_square = ratio_points
        .iter()
        .map(|point| point.0 * point.0)
        .sum::<f64>()
        / n;
    let time_constant_se = deficit_se * mean_square.sqrt();
    let (a_low, a_high) = (
        time_constant - T975_DF3 * time_constant_se,
        time_constant + T975_DF3 * time_constant_se,
    );
    eprintln!(
        "L3 linear extrapolation D/N = A + c/N: A={time_constant:.4}±\
         {time_constant_se:.4} CI=({a_low:.3}, {a_high:.3}) c={deficit_slope:.3}"
    );
    assert!(
        a_low > 0.15 && a_high < 1.0,
        "extrapolated time constant A={time_constant:.3} CI=({a_low:.3}, {a_high:.3}) \
         must be positive and inside (0, 1)"
    );
}

#[test]
fn intermediate_regime_grows_subpolynomially() {
    const SIZES: [usize; 4] = [128, 256, 512, 1024];
    const REALIZATIONS: u64 = 12;
    let ladder =
        measure_diameter_ladder(&SIZES, &[1.5], 1.0, REALIZATIONS, 0x4C33_3001).expect("ladder");
    let rows: Vec<_> = ladder.rows_for(1.5).collect();
    let means: Vec<f64> = rows.iter().map(|row| row.mean()).collect();
    let sems: Vec<f64> = rows.iter().map(|row| row.standard_error()).collect();
    eprintln!("L3 (1,2) s=1.5: means={means:?} sems={sems:?}");
    for index in 0..means.len() - 1 {
        assert!(
            means[index + 1] > means[index],
            "diameter must grow with N: {} -> {}",
            means[index],
            means[index + 1]
        );
    }
    let points: Vec<(f64, f64)> = rows
        .iter()
        .map(|row| ((row.size() as f64).ln(), row.mean().ln()))
        .collect();
    let (slope, _, slope_se, r_squared) = ols(&points);
    // t_{0.975} on df = 2 (four points minus two parameters).
    const T975_DF2: f64 = 4.3027;
    let high = slope + T975_DF2 * slope_se;
    eprintln!(
        "L3 (1,2) fit: slope={slope:.3}±{slope_se:.3} upper={high:.3} \
         R2={r_squared:.4} (no sharp D_N exponent is claimed)"
    );
    assert!(
        high < 0.5,
        "sub-polynomial bracket: log-log slope upper CI {high:.3} must stay below 0.5"
    );
}

#[test]
fn logarithmic_boundary_brackets_ln_n_over_ln_ln_n() {
    const SIZES: [usize; 3] = [512, 1024, 2048];
    const REALIZATIONS: u64 = 4;
    let ladder =
        measure_diameter_ladder(&SIZES, &[1.0], 1.0, REALIZATIONS, 0x4C33_4001).expect("ladder");
    let mut previous = 0.0;
    for row in ladder.rows_for(1.0) {
        let n = row.size() as f64;
        let log = n.ln();
        let scale = log / log.ln();
        let ratio = row.mean() / scale;
        eprintln!(
            "L3 s=1 N={} mean={:.3} sem={:.3} lnN={log:.3} scale={scale:.3} ratio={ratio:.3}",
            row.size(),
            row.mean(),
            row.standard_error()
        );
        assert!(
            row.mean() <= log,
            "s=1 diameter {:.3} must stay below ln N = {log:.3}",
            row.mean()
        );
        assert!(
            (0.3..=5.0).contains(&ratio),
            "s=1 diameter sits at {ratio:.3} times ln N/ln ln N, outside the bracket"
        );
        assert!(
            row.mean() >= previous - 1.0,
            "s=1 diameter must be non-decreasing within one integer step: \
             {previous:.3} -> {:.3}",
            row.mean()
        );
        previous = row.mean();
    }
}

#[test]
fn bounded_regime_diameter_does_not_grow() {
    const SIZES: [usize; 3] = [256, 512, 1024];
    const REALIZATIONS: u64 = 4;
    let ladder =
        measure_diameter_ladder(&SIZES, &[0.5], 1.0, REALIZATIONS, 0x4C33_5001).expect("ladder");
    let rows: Vec<_> = ladder.rows_for(0.5).collect();
    let first = rows.first().expect("row").mean();
    for row in &rows {
        assert_eq!(row.regime(), Some(DiameterRegime::Bounded));
        assert!(
            (2.0..=4.0).contains(&row.mean()),
            "s<1 diameter {:.3} must stay in the bounded envelope [2, 4]",
            row.mean()
        );
        let margin = 2.0 * (row.standard_error().powi(2) + rows[0].standard_error().powi(2)).sqrt();
        assert!(
            row.mean() <= first + margin,
            "s<1 diameter must not grow with N: {first:.3} -> {:.3} (margin {margin:.3})",
            row.mean()
        );
    }
    eprintln!(
        "L3 bounded s=0.5: means={:?} (constant sequence)",
        rows.iter().map(|row| row.mean()).collect::<Vec<_>>()
    );
}

#[test]
fn ladder_rejections_are_typed_and_stateless() {
    assert!(matches!(
        measure_diameter_ladder(&[], &[1.5], 1.0, 4, 0),
        Err(DiameterLadderError::EmptySizes)
    ));
    assert!(matches!(
        measure_diameter_ladder(&[8], &[], 1.0, 4, 0),
        Err(DiameterLadderError::EmptyExponents)
    ));
    assert!(matches!(
        measure_diameter_ladder(&[0], &[1.5], 1.0, 4, 0),
        Err(DiameterLadderError::ZeroSize)
    ));
    assert!(matches!(
        measure_diameter_ladder(&[4097], &[1.5], 1.0, 4, 0),
        Err(DiameterLadderError::SizeAboveCap { size: 4097, .. })
    ));
    assert!(matches!(
        measure_diameter_ladder(&[8], &[1.5], 1.0, 1, 0),
        Err(DiameterLadderError::TooFewRealizations { realizations: 1 })
    ));
    assert!(matches!(
        measure_diameter_ladder(&[8], &[f64::NAN], 1.0, 4, 0),
        Err(DiameterLadderError::ExponentInvalid { .. })
    ));
    assert!(matches!(
        measure_diameter_ladder(&[8], &[1.5], -1.0, 4, 0),
        Err(DiameterLadderError::Distance(_))
    ));
    eprintln!("L3 ladder rejections: 7 typed refusal cases");
}
