//! Non-Gaussian moment recovery (tracker P1.8).
//!
//! Two targets with analytic moments but non-Gaussian shape:
//!
//! * a symmetric bimodal Gaussian mixture `0.5 N(-1.5, 0.36) + 0.5 N(1.5, 0.36)`
//!   in one dimension (mean 0, variance 2.61), and
//! * a two-dimensional Student-t with `nu = 5` (mean and covariance
//!   `nu / (nu - 2) * scale`).
//!
//! Each target is sampled by at least three solvers with fixed seeds. Moment
//! recovery is asserted with z-scores against the analytic values, where the
//! Monte Carlo error comes from independent-seed replication (multi-seed
//! estimator).
//!
//! The mixture modes are separated by a 23x density valley: deep enough to be
//! clearly bimodal, shallow enough that every solver demonstrably crosses
//! between the modes within a chain (a deeper valley makes random-walk and
//! slice chains mode-sticky and puts HMC/NUTS in a rare-tunneling regime).
//! Chains are initialized alternating between the two mode centers, and the
//! lower-mode occupancy must be symmetric (0.5 within four Monte Carlo
//! standard errors and inside [0.35, 0.65] in absolute terms), which fails if
//! a solver systematically sticks to one mode.

use mcmc_rs::{
    DiagonalMetric, DifferentiableLogDensity, EuclideanState, LogDensity, McmcError, Nuts,
    RandomWalkMetropolis, SamplingPhase, SliceSampler, StaticHmc, TransitionKernel,
};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Repo-wide z-score tolerance for Monte Carlo agreement assertions.
const Z_LIMIT: f64 = 4.0;

/// Pooled estimate over independent replications with Monte Carlo error.
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

// ── Bimodal Gaussian mixture ────────────────────────────────────────────

const MODE_OFFSET: f64 = 1.5;
const MODE_STANDARD_DEVIATION: f64 = 0.6;
/// Var(x) = mode variance + squared mode separation: 0.36 + 2.25.
const MIXTURE_VARIANCE: f64 =
    MODE_OFFSET * MODE_OFFSET + MODE_STANDARD_DEVIATION * MODE_STANDARD_DEVIATION;

#[derive(Clone, Copy)]
struct BimodalMixture;

fn mixture_component_log_density(offset: f64, position: f64) -> f64 {
    -0.5 * ((position - offset) / MODE_STANDARD_DEVIATION).powi(2)
}

impl LogDensity<[f64]> for BimodalMixture {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        let left = mixture_component_log_density(-MODE_OFFSET, position[0]);
        let right = mixture_component_log_density(MODE_OFFSET, position[0]);
        let (high, low) = if left > right {
            (left, right)
        } else {
            (right, left)
        };
        high + (-(high - low)).exp().ln_1p() - std::f64::consts::LN_2
    }
}

impl DifferentiableLogDensity for BimodalMixture {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        let left = mixture_component_log_density(-MODE_OFFSET, position[0]);
        let right = mixture_component_log_density(MODE_OFFSET, position[0]);
        let shift = left.max(right);
        let left_weight = (left - shift).exp();
        let right_weight = (right - shift).exp();
        let total = left_weight + right_weight;
        let denominator = MODE_STANDARD_DEVIATION * MODE_STANDARD_DEVIATION;
        let left_gradient = -(position[0] + MODE_OFFSET) / denominator;
        let right_gradient = -(position[0] - MODE_OFFSET) / denominator;
        gradient[0] = (left_weight * left_gradient + right_weight * right_gradient) / total;
        self.log_density(position)
    }
}

/// Per-chain mixture statistics pooled over seeds below.
struct MixtureChain {
    mean: f64,
    variance: f64,
    /// Fraction of draws in the lower (x < 0) mode.
    lower_fraction: f64,
}

struct MixtureSummary {
    name: &'static str,
    mean: f64,
    mean_error: f64,
    variance: f64,
    variance_error: f64,
    lower_fraction: f64,
    lower_fraction_error: f64,
}

fn run_mixture_solver<K>(
    name: &'static str,
    make_kernel: impl Fn(usize) -> Result<K, McmcError>,
    warmup: usize,
    draws: usize,
    seeds: &[u64],
) -> MixtureSummary
where
    K: TransitionKernel<BimodalMixture>,
{
    let mut chains = Vec::with_capacity(seeds.len());
    for (index, seed) in seeds.iter().enumerate() {
        // Even chains start in the lower mode, odd chains in the upper mode.
        let initial = if index % 2 == 0 {
            -MODE_OFFSET
        } else {
            MODE_OFFSET
        };
        let mut target = BimodalMixture;
        let mut kernel = make_kernel(warmup).expect("valid kernel configuration");
        let mut state = EuclideanState::initialize(&mut target, vec![initial]).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(*seed);
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
        let mut sum = 0.0;
        let mut lower = 0_usize;
        let mut samples = Vec::with_capacity(draws);
        for _ in 0..draws {
            kernel
                .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
                .unwrap();
            let draw = state.position()[0];
            sum += draw;
            lower += usize::from(draw < 0.0);
            samples.push(draw);
        }
        let count = draws as f64;
        let mean = sum / count;
        let variance = samples
            .iter()
            .map(|draw| {
                let delta = draw - mean;
                delta * delta
            })
            .sum::<f64>()
            / (count - 1.0);
        chains.push(MixtureChain {
            mean,
            variance,
            lower_fraction: lower as f64 / count,
        });
    }

    let (mean, mean_error) = pooled(&chains.iter().map(|chain| chain.mean).collect::<Vec<_>>());
    let (variance, variance_error) = pooled(
        &chains
            .iter()
            .map(|chain| chain.variance)
            .collect::<Vec<_>>(),
    );
    let (lower_fraction, lower_fraction_error) = pooled(
        &chains
            .iter()
            .map(|chain| chain.lower_fraction)
            .collect::<Vec<_>>(),
    );
    MixtureSummary {
        name,
        mean,
        mean_error,
        variance,
        variance_error,
        lower_fraction,
        lower_fraction_error,
    }
}

fn assert_mixture_recovery(summary: &MixtureSummary) {
    assert_z(
        &format!("{} mixture mean", summary.name),
        summary.mean,
        summary.mean_error,
        0.0,
    );
    assert_z(
        &format!("{} mixture variance", summary.name),
        summary.variance,
        summary.variance_error,
        MIXTURE_VARIANCE,
    );
    // Symmetric mode occupancy: a solver that sticks to its initial mode
    // drags the pooled fraction away from 0.5 well beyond four errors.
    assert_z(
        &format!("{} mixture lower-mode fraction", summary.name),
        summary.lower_fraction,
        summary.lower_fraction_error,
        0.5,
    );
    // Both modes must be populated at all, independently of the error model.
    assert!(
        (0.35..=0.65).contains(&summary.lower_fraction),
        "{}: lower-mode fraction {} outside [0.35, 0.65]",
        summary.name,
        summary.lower_fraction
    );
}

fn evaluate_mixture(seed_base: u64, seed_count: usize, long_run: bool) -> Vec<MixtureSummary> {
    let seeds: Vec<u64> = (0..seed_count)
        .map(|index| seed_base + index as u64)
        .collect();
    let chain_scale = if long_run { 3 } else { 1 };
    vec![
        run_mixture_solver(
            "nuts",
            |warmup| {
                Nuts::new(DiagonalMetric::unit(1).unwrap(), 0.4, 6)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            400 * chain_scale,
            8_000 * chain_scale,
            &seeds,
        ),
        run_mixture_solver(
            "static-hmc",
            |warmup| {
                StaticHmc::new(DiagonalMetric::unit(1).unwrap(), 0.4, 8)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            400 * chain_scale,
            8_000 * chain_scale,
            &seeds,
        ),
        run_mixture_solver(
            "random-walk-metropolis",
            |_| {
                RandomWalkMetropolis::isotropic(1, 1.0)
                    .unwrap()
                    .with_scale_adaptation(0.44)
            },
            2_000 * chain_scale,
            15_000 * chain_scale,
            &seeds,
        ),
        run_mixture_solver(
            "slice-sampler",
            |_| SliceSampler::new(vec![1.5]),
            500 * chain_scale,
            15_000 * chain_scale,
            &seeds,
        ),
    ]
}

#[test]
fn mixture_solvers_recover_bimodal_moments() {
    let summaries = evaluate_mixture(0xB1A5_0000, 16, false);
    for summary in &summaries {
        assert_mixture_recovery(summary);
    }
}

#[test]
#[ignore = "high-power mixture replication with more seeds and longer chains (nightly)"]
fn mixture_solvers_recover_bimodal_moments_long_run() {
    let summaries = evaluate_mixture(0xB1A5_1000, 24, true);
    for summary in &summaries {
        assert_mixture_recovery(summary);
    }
}

// ── Two-dimensional Student-t ───────────────────────────────────────────

const T_MEAN: [f64; 2] = [0.5, -1.0];
const T_DEGREES_OF_FREEDOM: f64 = 5.0;
/// Row-major scale matrix Sigma: [[1.0, 0.3], [0.3, 0.5]].
const T_SCALE: [f64; 4] = [1.0, 0.3, 0.3, 0.5];
const T_SCALE_DETERMINANT: f64 = 1.0 * 0.5 - 0.3 * 0.3;
/// Row-major inverse scale matrix.
const T_SCALE_PRECISION: [f64; 4] = [
    0.5 / T_SCALE_DETERMINANT,
    -0.3 / T_SCALE_DETERMINANT,
    -0.3 / T_SCALE_DETERMINANT,
    1.0 / T_SCALE_DETERMINANT,
];
/// Exponent of the two-dimensional density: -(nu + d) / 2 = -3.5.
const T_EXPONENT: f64 = (T_DEGREES_OF_FREEDOM + 2.0) / 2.0;
/// Analytic covariance is nu / (nu - 2) * Sigma with nu = 5.
const T_COVARIANCE_FACTOR: f64 = T_DEGREES_OF_FREEDOM / (T_DEGREES_OF_FREEDOM - 2.0);
const ANALYTIC_T_COVARIANCE: [f64; 3] = [
    T_COVARIANCE_FACTOR * T_SCALE[0],
    T_COVARIANCE_FACTOR * T_SCALE[1],
    T_COVARIANCE_FACTOR * T_SCALE[3],
];
const ANALYTIC_T_SD: [f64; 2] = [1.290_994_448_735_805_6, 0.912_870_929_175_276_9];
/// Start away from the mean so warmup has to reach the bulk of the mass.
const T_INITIAL: [f64; 2] = [3.0, -3.0];

#[derive(Clone, Copy)]
struct StudentT;

impl StudentT {
    fn quadratic_form(position: &[f64]) -> f64 {
        let delta0 = position[0] - T_MEAN[0];
        let delta1 = position[1] - T_MEAN[1];
        T_SCALE_PRECISION[0] * delta0 * delta0
            + (T_SCALE_PRECISION[1] + T_SCALE_PRECISION[2]) * delta0 * delta1
            + T_SCALE_PRECISION[3] * delta1 * delta1
    }
}

impl LogDensity<[f64]> for StudentT {
    fn log_density(&mut self, position: &[f64]) -> f64 {
        -T_EXPONENT * (1.0 + Self::quadratic_form(position) / T_DEGREES_OF_FREEDOM).ln()
    }
}

impl DifferentiableLogDensity for StudentT {
    fn log_density_and_gradient(&mut self, position: &[f64], gradient: &mut [f64]) -> f64 {
        let delta0 = position[0] - T_MEAN[0];
        let delta1 = position[1] - T_MEAN[1];
        let scale = -T_EXPONENT / (T_DEGREES_OF_FREEDOM + Self::quadratic_form(position));
        gradient[0] = scale * 2.0 * (T_SCALE_PRECISION[0] * delta0 + T_SCALE_PRECISION[1] * delta1);
        gradient[1] = scale * 2.0 * (T_SCALE_PRECISION[2] * delta0 + T_SCALE_PRECISION[3] * delta1);
        self.log_density(position)
    }
}

struct StudentMoments {
    name: &'static str,
    mean: [f64; 2],
    mean_error: [f64; 2],
    covariance: [f64; 3],
    covariance_error: [f64; 3],
}

fn run_student_solver<K>(
    name: &'static str,
    make_kernel: impl Fn(usize) -> Result<K, McmcError>,
    warmup: usize,
    draws: usize,
    seeds: &[u64],
) -> StudentMoments
where
    K: TransitionKernel<StudentT>,
{
    let mut means = vec![[0.0; 2]; seeds.len()];
    let mut covariances = vec![[0.0; 3]; seeds.len()];
    for (index, seed) in seeds.iter().enumerate() {
        let mut target = StudentT;
        let mut kernel = make_kernel(warmup).expect("valid kernel configuration");
        let mut state = EuclideanState::initialize(&mut target, T_INITIAL.to_vec()).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(*seed);
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
        let mut sum = [0.0; 2];
        let mut samples = Vec::with_capacity(draws);
        for _ in 0..draws {
            kernel
                .transition(&mut target, &mut state, &mut rng, SamplingPhase::Sampling)
                .unwrap();
            let draw = [state.position()[0], state.position()[1]];
            sum[0] += draw[0];
            sum[1] += draw[1];
            samples.push(draw);
        }
        let count = draws as f64;
        let mean = [sum[0] / count, sum[1] / count];
        let mut covariance = [0.0; 3];
        for draw in &samples {
            let delta0 = draw[0] - mean[0];
            let delta1 = draw[1] - mean[1];
            covariance[0] += delta0 * delta0;
            covariance[1] += delta0 * delta1;
            covariance[2] += delta1 * delta1;
        }
        for element in &mut covariance {
            *element /= count - 1.0;
        }
        means[index] = mean;
        covariances[index] = covariance;
    }

    let mut pooled_mean = [0.0; 2];
    let mut pooled_mean_error = [0.0; 2];
    let mut pooled_covariance = [0.0; 3];
    let mut pooled_covariance_error = [0.0; 3];
    for coordinate in 0..2 {
        let (estimate, error) = pooled(
            &means
                .iter()
                .map(|mean| mean[coordinate])
                .collect::<Vec<_>>(),
        );
        pooled_mean[coordinate] = estimate;
        pooled_mean_error[coordinate] = error;
    }
    for element in 0..3 {
        let (estimate, error) = pooled(
            &covariances
                .iter()
                .map(|covariance| covariance[element])
                .collect::<Vec<_>>(),
        );
        pooled_covariance[element] = estimate;
        pooled_covariance_error[element] = error;
    }
    StudentMoments {
        name,
        mean: pooled_mean,
        mean_error: pooled_mean_error,
        covariance: pooled_covariance,
        covariance_error: pooled_covariance_error,
    }
}

fn assert_student_recovery(moments: &StudentMoments) {
    for (coordinate, ((estimate, error), reference)) in moments
        .mean
        .iter()
        .zip(&moments.mean_error)
        .zip(T_MEAN.iter())
        .enumerate()
    {
        // Mixing-power guard: at least ~400 effective draws per solver.
        assert!(
            *error < ANALYTIC_T_SD[coordinate] / 20.0,
            "{}: mean[{coordinate}] error {error} exceeds analytic sd / 20 — \
             insufficient effective sample size",
            moments.name
        );
        assert_z(
            &format!("{} student-t mean[{coordinate}]", moments.name),
            *estimate,
            *error,
            *reference,
        );
    }
    for (element, ((estimate, error), reference)) in moments
        .covariance
        .iter()
        .zip(&moments.covariance_error)
        .zip(ANALYTIC_T_COVARIANCE.iter())
        .enumerate()
    {
        assert_z(
            &format!("{} student-t covariance[{element}]", moments.name),
            *estimate,
            *error,
            *reference,
        );
    }
}

fn evaluate_student_t(seed_base: u64, seed_count: usize, long_run: bool) -> Vec<StudentMoments> {
    let seeds: Vec<u64> = (0..seed_count)
        .map(|index| seed_base + index as u64)
        .collect();
    let chain_scale = if long_run { 3 } else { 1 };
    vec![
        run_student_solver(
            "nuts",
            |warmup| {
                Nuts::new(DiagonalMetric::unit(2).unwrap(), 0.3, 8)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            500 * chain_scale,
            5_000 * chain_scale,
            &seeds,
        ),
        run_student_solver(
            "static-hmc",
            |warmup| {
                StaticHmc::new(DiagonalMetric::unit(2).unwrap(), 0.25, 8)
                    .unwrap()
                    .with_diagonal_adaptation(warmup as u64, 0.8, 1.0e-3)
            },
            500 * chain_scale,
            5_000 * chain_scale,
            &seeds,
        ),
        run_student_solver(
            "slice-sampler",
            |_| SliceSampler::new(vec![3.0, 2.0]),
            1_000 * chain_scale,
            10_000 * chain_scale,
            &seeds,
        ),
        run_student_solver(
            "random-walk-metropolis",
            |_| {
                RandomWalkMetropolis::isotropic(2, 1.0)
                    .unwrap()
                    .with_scale_adaptation(0.234)
                    .unwrap()
                    .with_dense_covariance_adaptation(1.0e-4)
            },
            4_000 * chain_scale,
            30_000 * chain_scale,
            &seeds,
        ),
    ]
}

#[test]
fn student_t_solvers_recover_heavy_tailed_moments() {
    let moments = evaluate_student_t(0x57D5_0000, 8, false);
    for summary in &moments {
        assert_student_recovery(summary);
    }
}

#[test]
#[ignore = "high-power Student-t replication with more seeds and longer chains (nightly)"]
fn student_t_solvers_recover_heavy_tailed_moments_long_run() {
    let moments = evaluate_student_t(0x57D5_1000, 16, true);
    for summary in &moments {
        assert_student_recovery(summary);
    }
}
