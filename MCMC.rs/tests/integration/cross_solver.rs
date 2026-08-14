//! Cross-solver posterior agreement (tracker P1.5).
//!
//! Six samplers — adaptive random-walk Metropolis, componentwise Metropolis,
//! slice sampling, StaticHMC, NUTS and an exact-conditional Gibbs sweep —
//! target the same correlated two-dimensional Gaussian. For every solver the
//! posterior mean vector and the three distinct covariance elements are
//! estimated from independent fixed-seed chains. Each solver must reproduce
//! the analytic moments, and every pair of solvers must agree on each moment
//! within four combined Monte Carlo standard errors.
//!
//! Monte Carlo errors are estimated by independent-seed replication (the
//! multi-seed estimator is also valid for the covariance elements, which the
//! public ESS API does not cover): the pooled estimate is the mean over seeds
//! and its error is the seed-level standard error.

use mcmc_rs::proposal::standard_normal;
use mcmc_rs::{
    ComponentWiseMetropolis, DiagonalMetric, DifferentiableLogDensity, EuclideanState, GibbsKernel,
    GibbsUpdate, GibbsUpdateResult, LogDensity, McmcError, Nuts, RandomWalkMetropolis,
    SamplingPhase, SliceSampler, StaticHmc, Then, TransitionKernel,
};
use rand::Rng;
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Shared target mean: mu = (1, -2).
const MEAN: [f64; 2] = [1.0, -2.0];
/// Shared target covariance, row major: [[2.0, 0.6], [0.6, 1.0]].
const COVARIANCE: [f64; 4] = [2.0, 0.6, 0.6, 1.0];
const DETERMINANT: f64 = 2.0 * 1.0 - 0.6 * 0.6;
/// Row-major precision (inverse covariance) of the shared target.
const PRECISION: [f64; 4] = [
    1.0 / DETERMINANT,
    -0.6 / DETERMINANT,
    -0.6 / DETERMINANT,
    2.0 / DETERMINANT,
];
/// Distinct covariance elements in the order used by the moment estimates.
const ANALYTIC_COVARIANCE: [f64; 3] = [COVARIANCE[0], COVARIANCE[1], COVARIANCE[3]];
/// Every chain starts away from the mean so warmup has to find the mode.
const INITIAL: [f64; 2] = [3.0, 0.0];
/// Repo-wide z-score tolerance for Monte Carlo agreement assertions.
const Z_LIMIT: f64 = 4.0;

#[derive(Clone, Copy)]
struct CorrelatedGaussian;

impl LogDensity<[f64]> for CorrelatedGaussian {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        let delta0 = position[0] - MEAN[0];
        let delta1 = position[1] - MEAN[1];
        let quadratic = PRECISION[0] * delta0 * delta0
            + (PRECISION[1] + PRECISION[2]) * delta0 * delta1
            + PRECISION[3] * delta1 * delta1;
        -0.5 * quadratic
    }
}

impl DifferentiableLogDensity for CorrelatedGaussian {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        let delta0 = position[0] - MEAN[0];
        let delta1 = position[1] - MEAN[1];
        gradient[0] = -(PRECISION[0] * delta0 + PRECISION[1] * delta1);
        gradient[1] = -(PRECISION[2] * delta0 + PRECISION[3] * delta1);
        self.log_density(position)
    }
}

/// Exact draw from the Gaussian full conditional of one coordinate.
#[derive(Clone, Copy)]
struct ConditionalDraw {
    coordinate: usize,
}

impl GibbsUpdate<CorrelatedGaussian> for ConditionalDraw {
    fn update<R>(
        &mut self,
        _target: &mut CorrelatedGaussian,
        current: &EuclideanState,
        proposed_position: &mut [f64],
        rng: &mut R,
        _phase: SamplingPhase,
    ) -> Result<GibbsUpdateResult, McmcError>
    where
        R: Rng + ?Sized,
    {
        let other = 1 - self.coordinate;
        let row = 2 * self.coordinate;
        let conditional_mean = MEAN[self.coordinate]
            - PRECISION[row + other] / PRECISION[row + self.coordinate]
                * (current.position()[other] - MEAN[other]);
        let conditional_sd = (1.0 / PRECISION[row + self.coordinate]).sqrt();
        proposed_position[self.coordinate] =
            conditional_mean + conditional_sd * standard_normal(rng);
        Ok(GibbsUpdateResult::requiring_target_evaluation())
    }

    fn name(&self, _target: &CorrelatedGaussian) -> &'static str {
        "GaussianConditionalDraw"
    }
}

/// Per-chain posterior moments of the two coordinates.
struct ChainMoments {
    mean: [f64; 2],
    /// Distinct covariance elements: (c00, c01, c11), unbiased.
    covariance: [f64; 3],
}

fn chain_moments(draws: &[[f64; 2]]) -> ChainMoments {
    let count = draws.len() as f64;
    let mean = [
        draws.iter().map(|draw| draw[0]).sum::<f64>() / count,
        draws.iter().map(|draw| draw[1]).sum::<f64>() / count,
    ];
    let mut covariance = [0.0; 3];
    for draw in draws {
        let delta0 = draw[0] - mean[0];
        let delta1 = draw[1] - mean[1];
        covariance[0] += delta0 * delta0;
        covariance[1] += delta0 * delta1;
        covariance[2] += delta1 * delta1;
    }
    for element in &mut covariance {
        *element /= count - 1.0;
    }
    ChainMoments { mean, covariance }
}

/// Pooled estimate over independent seeds and its Monte Carlo standard error.
fn pooled(values: &[f64]) -> (f64, f64) {
    let count = values.len() as f64;
    let mean = values.iter().sum::<f64>() / count;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / (count - 1.0);
    (mean, (variance / count).sqrt())
}

/// Pooled moments of one solver with seed-level Monte Carlo errors.
struct SolverMoments {
    name: &'static str,
    mean: [f64; 2],
    mean_error: [f64; 2],
    covariance: [f64; 3],
    covariance_error: [f64; 3],
}

fn run_solver<K>(
    name: &'static str,
    make_kernel: impl Fn(usize) -> Result<K, McmcError>,
    warmup: usize,
    draws: usize,
    seeds: &[u64],
) -> SolverMoments
where
    K: TransitionKernel<CorrelatedGaussian>,
{
    let mut chains = Vec::with_capacity(seeds.len());
    for seed in seeds.iter().copied() {
        let mut target = CorrelatedGaussian;
        let mut kernel = make_kernel(warmup).expect("valid kernel configuration");
        let mut state = EuclideanState::initialize(&mut target, INITIAL.to_vec()).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        kernel
            .on_phase_start(&mut target, SamplingPhase::Warmup, &state)
            .unwrap();
        for _ in 0..warmup {
            kernel
                .transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)
                .unwrap();
        }
        kernel
            .on_phase_end(&mut target, SamplingPhase::Warmup, &state)
            .unwrap();
        kernel
            .on_phase_start(&mut target, SamplingPhase::Sampling, &state)
            .unwrap();
        let mut samples = Vec::with_capacity(draws);
        for _ in 0..draws {
            kernel
                .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
                .unwrap();
            samples.push([state.position()[0], state.position()[1]]);
        }
        chains.push(chain_moments(&samples));
    }

    let mut mean = [0.0; 2];
    let mut mean_error = [0.0; 2];
    let mut covariance = [0.0; 3];
    let mut covariance_error = [0.0; 3];
    for coordinate in 0..2 {
        let (estimate, error) = pooled(
            &chains
                .iter()
                .map(|chain| chain.mean[coordinate])
                .collect::<Vec<_>>(),
        );
        mean[coordinate] = estimate;
        mean_error[coordinate] = error;
    }
    for element in 0..3 {
        let (estimate, error) = pooled(
            &chains
                .iter()
                .map(|chain| chain.covariance[element])
                .collect::<Vec<_>>(),
        );
        covariance[element] = estimate;
        covariance_error[element] = error;
    }
    SolverMoments {
        name,
        mean,
        mean_error,
        covariance,
        covariance_error,
    }
}

fn assert_z(context: &str, estimate: f64, standard_error: f64, reference: f64) {
    assert!(
        standard_error.is_finite() && standard_error > 0.0,
        "{context}: Monte Carlo standard error must be positive, got {standard_error}"
    );
    let deviation = (estimate - reference).abs();
    assert!(
        deviation < Z_LIMIT * standard_error,
        "{context}: estimate {estimate:.5} vs reference {reference:.5} \
         (error {standard_error:.5}, z {:.2})",
        deviation / standard_error,
    );
}

/// Mixing-power guard: the pooled mean error must beat analytic_sd / 20,
/// i.e. the pooled chains carry at least ~400 effective draws.
fn assert_mixing_power(solver: &SolverMoments) {
    for (coordinate, (error, reference_variance)) in
        solver.mean_error.iter().zip(COVARIANCE.iter()).enumerate()
    {
        assert!(
            *error < reference_variance.sqrt() / 20.0,
            "{}: mean[{coordinate}] error {error} exceeds analytic sd / 20 — \
             insufficient effective sample size",
            solver.name
        );
    }
}

fn assert_matches_analytic(solver: &SolverMoments) {
    for (coordinate, ((estimate, error), reference)) in solver
        .mean
        .iter()
        .zip(&solver.mean_error)
        .zip(MEAN.iter())
        .enumerate()
    {
        assert_z(
            &format!("{} analytic mean[{coordinate}]", solver.name),
            *estimate,
            *error,
            *reference,
        );
    }
    for (element, ((estimate, error), reference)) in solver
        .covariance
        .iter()
        .zip(&solver.covariance_error)
        .zip(ANALYTIC_COVARIANCE.iter())
        .enumerate()
    {
        assert_z(
            &format!("{} analytic covariance[{element}]", solver.name),
            *estimate,
            *error,
            *reference,
        );
    }
}

fn assert_pairwise_agreement(left: &SolverMoments, right: &SolverMoments) {
    let pair = format!("{} vs {}", left.name, right.name);
    for (coordinate, ((left_estimate, left_error), (right_estimate, right_error))) in left
        .mean
        .iter()
        .zip(&left.mean_error)
        .zip(right.mean.iter().zip(&right.mean_error))
        .enumerate()
    {
        let error = (left_error.powi(2) + right_error.powi(2)).sqrt();
        assert_z(
            &format!("{pair} mean[{coordinate}]"),
            *left_estimate,
            error,
            *right_estimate,
        );
    }
    for (element, ((left_estimate, left_error), (right_estimate, right_error))) in left
        .covariance
        .iter()
        .zip(&left.covariance_error)
        .zip(right.covariance.iter().zip(&right.covariance_error))
        .enumerate()
    {
        let error = (left_error.powi(2) + right_error.powi(2)).sqrt();
        assert_z(
            &format!("{pair} covariance[{element}]"),
            *left_estimate,
            error,
            *right_estimate,
        );
    }
}

fn evaluate_all_solvers(seed_base: u64, seed_count: usize, long_run: bool) -> Vec<SolverMoments> {
    let seeds: Vec<u64> = (0..seed_count)
        .map(|index| seed_base + index as u64)
        .collect();
    let chain_scale = if long_run { 3 } else { 1 };
    vec![
        run_solver(
            "random-walk-metropolis",
            |_| {
                RandomWalkMetropolis::isotropic(2, 1.0)
                    .unwrap()
                    .with_scale_adaptation(0.234)
                    .unwrap()
                    .with_dense_covariance_adaptation(1.0e-4)
            },
            3_000 * chain_scale,
            20_000 * chain_scale,
            &seeds,
        ),
        run_solver(
            "componentwise-metropolis",
            |_| {
                ComponentWiseMetropolis::new(vec![1.4, 1.0])
                    .unwrap()
                    .with_scale_adaptation(0.44)
            },
            3_000 * chain_scale,
            20_000 * chain_scale,
            &seeds,
        ),
        run_solver(
            "slice-sampler",
            |_| SliceSampler::new(vec![2.5, 2.0]),
            1_000 * chain_scale,
            10_000 * chain_scale,
            &seeds,
        ),
        run_solver(
            "static-hmc",
            |warmup| {
                StaticHmc::new(DiagonalMetric::unit(2).unwrap(), 0.2, 8)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            400 * chain_scale,
            4_000 * chain_scale,
            &seeds,
        ),
        run_solver(
            "nuts",
            |warmup| {
                Nuts::new(DiagonalMetric::unit(2).unwrap(), 0.3, 8)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            400 * chain_scale,
            4_000 * chain_scale,
            &seeds,
        ),
        run_solver(
            "gibbs-sweep",
            |_| {
                Ok(Then::new(
                    GibbsKernel::new(ConditionalDraw { coordinate: 0 }),
                    GibbsKernel::new(ConditionalDraw { coordinate: 1 }),
                ))
            },
            100 * chain_scale,
            5_000 * chain_scale,
            &seeds,
        ),
    ]
}

fn assert_cross_solver_agreement(solvers: &[SolverMoments]) {
    for solver in solvers {
        assert_mixing_power(solver);
        assert_matches_analytic(solver);
    }
    for left in 0..solvers.len() {
        for right in (left + 1)..solvers.len() {
            assert_pairwise_agreement(&solvers[left], &solvers[right]);
        }
    }
}

#[test]
fn six_solvers_agree_on_correlated_gaussian_moments() {
    let solvers = evaluate_all_solvers(0x5EED_0000, 8, false);
    assert_cross_solver_agreement(&solvers);
}

#[test]
#[ignore = "high-power replication with more seeds and longer chains (nightly)"]
fn six_solvers_agree_on_correlated_gaussian_moments_long_run() {
    let solvers = evaluate_all_solvers(0x5EED_1000, 12, true);
    assert_cross_solver_agreement(&solvers);
}
