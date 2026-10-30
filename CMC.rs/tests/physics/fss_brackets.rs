//! A5: generic FSS bracket reproduction.
//!
//! The v2 C2 calibration is the iid (white-noise) site-percolation threshold
//! on the square lattice. These tests rebuild that bracket with the generic
//! toolkit only: independent `SiteBernoulli` realizations feed
//! [`BinderEstimator`] and [`Moments`], the per-size curves are [`FssCurve`]s,
//! and [`crossing`] extracts the finite-size estimate. No stage-specific
//! estimator is involved.

use cmc_rs::{
    analyze, build_square, BinderEstimator, ComponentWorkspace, FssCurve, FssError, FssPoint,
    GraphView, Moments, Probability, SiteBernoulli, StaticConfiguration, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const IID_SITE_PC: f64 = 0.592_746;

struct Ladder {
    sides: Vec<usize>,
    largest: Vec<FssCurve>,
    binder: Vec<FssCurve>,
}

fn side_seed(seed: u64, side: usize, probability: f64) -> u64 {
    seed ^ (side as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ probability.to_bits()
}

fn sample_ladder(sides: &[usize], probabilities: &[f64], replicas: usize, seed: u64) -> Ladder {
    let mut largest_curves = Vec::with_capacity(sides.len());
    let mut binder_curves = Vec::with_capacity(sides.len());
    for &side in sides {
        let lattice = build_square(side, side, true);
        let vertices = lattice.vertex_count();
        let edges = lattice.edge_count();
        let mut largest_points = Vec::with_capacity(probabilities.len());
        let mut binder_points = Vec::with_capacity(probabilities.len());
        for &probability in probabilities {
            let law = SiteBernoulli::new(
                Probability::new(probability)
                    .expect("test probability inside [0, 1]")
                    .into(),
            );
            let mut configuration = StaticConfiguration::new(vertices, edges);
            let mut workspace = ComponentWorkspace::new();
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(side_seed(seed, side, probability));
            let mut moments = Moments::new();
            let mut binder = BinderEstimator::new();
            for _ in 0..replicas {
                law.sample(&lattice, &mut configuration, &mut rng)
                    .expect("law and configuration domains match");
                let result = analyze(
                    &lattice,
                    &configuration,
                    &configuration,
                    &[],
                    &mut workspace,
                )
                .expect("analysis domains match");
                let largest = result.summary().largest_component_size;
                moments.add(largest as f64 / vertices as f64);
                binder.add(largest, vertices).expect("largest <= V");
            }
            largest_points.push(
                FssPoint::new(
                    side,
                    probability,
                    moments.mean().expect("replicas > 0"),
                    moments.standard_error().unwrap_or(0.0),
                )
                .expect("valid largest-fraction point"),
            );
            let cumulant = binder.estimate().expect("two or more replicas");
            binder_points.push(
                FssPoint::new(side, probability, cumulant.mean, cumulant.standard_error)
                    .expect("valid Binder point"),
            );
        }
        largest_curves.push(FssCurve::new(largest_points).expect("one size, ascending grid"));
        binder_curves.push(FssCurve::new(binder_points).expect("one size, ascending grid"));
    }
    Ladder {
        sides: sides.to_vec(),
        largest: largest_curves,
        binder: binder_curves,
    }
}

fn adjacent_crossings(curves: &[FssCurve]) -> Vec<cmc_rs::CrossingEstimate> {
    curves
        .windows(2)
        .map(|pair| {
            cmc_rs::crossing(&pair[0], &pair[1]).unwrap_or_else(|error| {
                panic!(
                    "sizes {}/{} do not cross: {error:?}",
                    pair[0].size(),
                    pair[1].size()
                )
            })
        })
        .collect()
}

#[test]
fn generic_fss_brackets_the_iid_site_threshold() {
    let sides = [32, 48, 64];
    let mut probabilities = Vec::new();
    let mut probability = 0.54;
    while probability <= 0.65 + 1e-12 {
        probabilities.push(probability);
        probability += 0.005;
    }
    let ladder = sample_ladder(&sides, &probabilities, 128, 0x4135_0001);

    let largest = adjacent_crossings(&ladder.largest);
    let binder = adjacent_crossings(&ladder.binder);
    for (crossing, sizes) in binder.iter().zip(sides.windows(2)) {
        eprintln!(
            "A5 Binder crossing L={}/{}: {:.4} +/- {:.4} (exact {IID_SITE_PC:.6})",
            sizes[0], sizes[1], crossing.probability, crossing.standard_error
        );
    }
    for (crossing, sizes) in largest.iter().zip(sides.windows(2)) {
        eprintln!(
            "A5 largest-fraction crossing L={}/{}: {:.4} +/- {:.4} (exact {IID_SITE_PC:.6})",
            sizes[0], sizes[1], crossing.probability, crossing.standard_error
        );
        assert!(
            (crossing.probability - IID_SITE_PC).abs() < 0.06,
            "crossing {:.4} is outside the 0.06 finite-size bracket around {IID_SITE_PC:.6}",
            crossing.probability
        );
    }
    for crossing in &binder {
        assert!(
            (crossing.probability - IID_SITE_PC).abs() < 0.05,
            "Binder crossing {:.4} is outside the 0.05 bracket around {IID_SITE_PC:.6}",
            crossing.probability
        );
    }
    // Averaging the adjacent crossings cancels part of the first-sign-change
    // noise at small L; the mean must sit inside a tighter bracket.
    let mean_largest =
        largest.iter().map(|value| value.probability).sum::<f64>() / largest.len() as f64;
    eprintln!("A5 mean largest-fraction crossing: {mean_largest:.4}");
    assert!(
        (mean_largest - IID_SITE_PC).abs() < 0.03,
        "mean largest-fraction crossing {mean_largest:.4} is outside the 0.03 bracket"
    );
    // Finite-size crossings must land inside the sampled window, and the
    // largest-fraction estimator is the better-conditioned one here.
    assert!(largest
        .iter()
        .all(|crossing| crossing.probability > 0.53 && crossing.probability < 0.65));
    assert_eq!(ladder.sides, sides);
    assert_eq!(largest.len(), sides.len() - 1);
    assert_eq!(binder.len(), sides.len() - 1);
}

#[test]
#[ignore = "long finite-size FSS calibration (L=64/128, 32 replicas); run in nightly release"]
fn generic_fss_reproduces_the_v2_white_noise_calibration() {
    // v2's C2 crossing calibration (ThresholdLadder, a=+inf) landed at
    // 0.5987 against the exact 0.592746. This test uses an independent
    // sampler and the generic toolkit; it asserts the same 0.02 bracket, not
    // the v2 point value.
    let sides = [64, 128];
    let mut probabilities = Vec::new();
    let mut probability = 0.56;
    while probability <= 0.64 + 1e-12 {
        probabilities.push(probability);
        probability += 0.002;
    }
    let ladder = sample_ladder(&sides, &probabilities, 32, 0x4135_0002);
    let crossings = adjacent_crossings(&ladder.largest);
    let crossing = crossings[0];
    eprintln!(
        "A5 white-noise crossing L=64/128: {:.4} +/- {:.4} (exact {IID_SITE_PC:.6}; v2 0.5987)",
        crossing.probability, crossing.standard_error
    );
    assert!(
        (crossing.probability - IID_SITE_PC).abs() < 0.02,
        "white-noise crossing {:.4} is outside the v2 0.02 calibration bracket",
        crossing.probability
    );
    // The reported interval is a first-order bracket; it must be finite and
    // contain the point estimate.
    assert!(crossing.lower <= crossing.probability && crossing.probability <= crossing.upper);
    assert!(crossing.standard_error.is_finite());
    assert!(matches!(
        cmc_rs::crossing(&ladder.largest[0], &ladder.largest[0]),
        Err(FssError::SizeOrder { .. })
    ));
}
