//! Statistical validation of the L0 VMC sampler on He-4-like droplets.
//!
//! Every asserted number traces to a derivation stated in a comment; no
//! literature energy windows are used. The defensible anchors are:
//!
//! 1. A rigorous operator lower bound: `T = <psi| -1/2 sum lap^2 |psi>
//!    = 1/2 int |grad psi|^2 >= 0`, and every Lennard-Jones pair
//!    `4 eps[(sigma/r)^12 - (sigma/r)^6] >= -eps` (minimum at
//!    `r = 2^(1/6) sigma`). Hence for ANY state
//!    `E = <T> + <V> >= -eps * N(N-1)/2`.
//! 2. The exact closed-form VMC energy of a Gaussian trial state in a
//!    harmonic trap at arbitrary width (derived below), which tests the
//!    Metropolis sampling itself, not just the estimator.
//! 3. Multi-seed self-consistency z-scores, stationarity (batch means),
//!    acceptance-ratio sanity, and a significantly negative sampled mean
//!    Lennard-Jones pair energy (the droplet is correlated through the
//!    attractive well).
//!
//! The He-4-like state is the confined McMillan droplet
//! `psi = GaussianTrap(alpha) x McMillanJastrow(b)` under
//! `H = trap(omega) + LennardJones(eps, sigma)`: the open L0 box has no
//! periodic boundaries, so the one-body Gaussian provides the confinement
//! that a bulk liquid would get from PBC (documented L0 debt). Reducing
//! the interaction entirely (pure Gaussian) is the dilute reference.

use qmc_rs::{
    ContinuumHamiltonian, GaussianTrap, HarmonicTrap, McMillanJastrow, PairPotential, Positions,
    Product, VmcKernel, WaveFunctionParams, DIM,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// The He-4-like droplet trial state and Hamiltonian used by the tests
/// below: N = 8, Lennard-Jones in reduced units (eps = sigma = 1),
/// McMillan b = 1.0, one-body Gaussian width alpha = 0.15 under a weak
/// trap omega = 0.05. Choices are stated, not asserted.
fn droplet_hamiltonian() -> ContinuumHamiltonian {
    ContinuumHamiltonian::new(
        Some(HarmonicTrap::new(0.05, [0.0; DIM]).unwrap()),
        Some(PairPotential::LennardJones {
            epsilon: 1.0,
            sigma: 1.0,
        }),
    )
    .unwrap()
}

fn droplet_wave_function() -> Product<GaussianTrap, McMillanJastrow> {
    Product::new(
        GaussianTrap::new(0.15, [0.0; DIM]).unwrap(),
        McMillanJastrow::new(1.0).unwrap(),
    )
}

/// Batch-means standard error of the mean (same estimator as the lattice
/// z-score tests).
fn batch_means_stderr(samples: &[f64], binsize: usize) -> f64 {
    let n_bins = samples.len() / binsize;
    assert!(n_bins >= 4, "need at least 4 bins for stderr estimate");
    let bin_means: Vec<f64> = (0..n_bins)
        .map(|b| {
            let start = b * binsize;
            samples[start..start + binsize].iter().sum::<f64>() / binsize as f64
        })
        .collect();
    let grand_mean = bin_means.iter().sum::<f64>() / n_bins as f64;
    let variance = bin_means
        .iter()
        .map(|&m| (m - grand_mean).powi(2))
        .sum::<f64>()
        / (n_bins - 1) as f64;
    (variance / n_bins as f64).sqrt()
}

/// One VMC run collecting per-sweep population-mean local energies.
/// Returns `(samples, acceptance_ratio)`.
#[allow(clippy::too_many_arguments)]
fn run_vmc_samples<W>(
    seed: u64,
    wave_function: W,
    hamiltonian: ContinuumHamiltonian,
    n_walkers: usize,
    n_particles: usize,
    initial_spread: f64,
    proposal_width: f64,
    thermalization: usize,
    measurement: usize,
) -> (Vec<f64>, f64)
where
    W: WaveFunctionParams<Config = Positions>,
{
    let mut rng = Rng::seed_from_u64(seed);
    let mut kernel = VmcKernel::new(
        wave_function,
        hamiltonian,
        n_walkers,
        n_particles,
        initial_spread,
        proposal_width,
        &mut rng,
    )
    .expect("valid kernel inputs");
    let mut samples = Vec::with_capacity(measurement);
    for sweep in 0..(thermalization + measurement) {
        let phase = if sweep < thermalization {
            carlo_rs::RngPhase::Thermalization
        } else {
            carlo_rs::RngPhase::Measurement
        };
        kernel.sweep_with_phase(&mut rng, phase);
        if sweep >= thermalization {
            samples.push(kernel.population_mean_local_energy().value);
        }
    }
    (samples, kernel.stats().acceptance_ratio())
}

fn mean_of(samples: &[f64]) -> f64 {
    samples.iter().sum::<f64>() / samples.len() as f64
}

/// Droplet run: collects per-sweep population-mean local energies AND
/// per-sweep population-mean Lennard-Jones pair energies (the pair part
/// isolated through a pair-only Hamiltonian).
fn run_vmc_droplet(
    seed: u64,
    n_walkers: usize,
    thermalization: usize,
    measurement: usize,
) -> (Vec<f64>, f64, Vec<f64>) {
    let pair_only = ContinuumHamiltonian::pair_only(PairPotential::LennardJones {
        epsilon: 1.0,
        sigma: 1.0,
    })
    .unwrap();
    let mut rng = Rng::seed_from_u64(seed);
    let mut kernel = VmcKernel::new(
        droplet_wave_function(),
        droplet_hamiltonian(),
        n_walkers,
        8,
        1.8,
        0.7,
        &mut rng,
    )
    .expect("valid kernel inputs");
    let mut samples = Vec::with_capacity(measurement);
    let mut pair_samples = Vec::with_capacity(measurement);
    for sweep in 0..(thermalization + measurement) {
        let phase = if sweep < thermalization {
            carlo_rs::RngPhase::Thermalization
        } else {
            carlo_rs::RngPhase::Measurement
        };
        kernel.sweep_with_phase(&mut rng, phase);
        if sweep >= thermalization {
            samples.push(kernel.population_mean_local_energy().value);
            let pair_mean = kernel
                .walkers()
                .iter()
                .map(|walker| pair_only.potential_energy(walker.configuration()))
                .sum::<f64>()
                / n_walkers as f64;
            pair_samples.push(pair_mean);
        }
    }
    (samples, kernel.stats().acceptance_ratio(), pair_samples)
}

#[test]
fn mcmillan_he4_droplet_respects_rigorous_bound_and_mixes() {
    // N = 8 confined Lennard-Jones bosons with the McMillan factor.
    let n_particles = 8;
    let (samples, acceptance, pair_samples) = run_vmc_droplet(0x0E4, 12, 600, 3000);

    let mean = mean_of(&samples);
    let stderr = batch_means_stderr(&samples, 100);

    // Rigorous bound (module docs): E >= -eps * N(N-1)/2 = -28 exactly,
    // for ANY state. (Loose at N = 8 — but it is a theorem, and the
    // z-score/stationarity checks below carry the statistical weight.)
    let lower_bound = -(n_particles as f64) * (n_particles - 1) as f64 / 2.0;
    assert!(
        mean - lower_bound >= -5.0 * stderr,
        "E_VMC = {mean} +/- {stderr} violates the rigorous bound {lower_bound}"
    );

    // Sanity band for a working Metropolis kernel: neither frozen
    // (acceptance ~ 0) nor free-diffusing (acceptance ~ 1).
    assert!(
        (0.15..0.95).contains(&acceptance),
        "implausible acceptance ratio {acceptance}"
    );

    // The sampled mean Lennard-Jones pair energy must be significantly
    // negative: the pair distances concentrate on the attractive side
    // (r > 2^(1/6) sigma), while the McMillan core suppresses the positive
    // r^-12 region. A positive mean would mean the sampler never explores
    // pair distances of order the well.
    let pair_mean = mean_of(&pair_samples);
    let pair_stderr = batch_means_stderr(&pair_samples, 100);
    assert!(
        pair_mean + 5.0 * pair_stderr < 0.0,
        "mean LJ pair energy {pair_mean} +/- {pair_stderr} is not negative"
    );

    // Stationarity: first vs second half of the measurement series must
    // agree within 4 combined sigmas (autocorrelation sanity).
    let half = samples.len() / 2;
    let (m1, m2) = (mean_of(&samples[..half]), mean_of(&samples[half..]));
    let (s1, s2) = (
        batch_means_stderr(&samples[..half], 100),
        batch_means_stderr(&samples[half..], 100),
    );
    let z = (m1 - m2).abs() / (s1 * s1 + s2 * s2).sqrt();
    assert!(z < 4.0, "non-stationary energy series: halves z = {z}");
}

#[test]
fn mcmillan_energy_multi_seed_self_consistency() {
    // Four independent seeds must agree within combined batch-means errors.
    let seeds = [0x1111_u64, 0x2222, 0x3333, 0x4444];
    let runs: Vec<(f64, f64)> = seeds
        .iter()
        .map(|&seed| {
            let (samples, _, _) = run_vmc_droplet(seed, 8, 300, 1500);
            (mean_of(&samples), batch_means_stderr(&samples, 100))
        })
        .collect();
    for i in 0..runs.len() {
        for j in (i + 1)..runs.len() {
            let (mi, si) = runs[i];
            let (mj, sj) = runs[j];
            let z = (mi - mj).abs() / (si * si + sj * sj).sqrt();
            assert!(z < 4.0, "seeds {i}/{j} inconsistent: z = {z}");
        }
    }
}

#[test]
fn gaussian_trial_state_samples_its_exact_closed_form_energy() {
    // Derivation (Gaussian trap, arbitrary alpha): under
    // |psi_alpha|^2 ∝ exp(-2 alpha |r - r0|^2) each Cartesian coordinate
    // has variance 1/(4 alpha), so <sum |r - r0|^2> = 3N/(4 alpha) and
    //   E(alpha) = 3N alpha + (1/2 w^2 - 2 alpha^2) * 3N/(4 alpha)
    //            = N [ (3/2) alpha + 3 w^2 / (8 alpha) ].
    // By AM-GM this is >= N * 3w/2 = E_0 with equality iff alpha = w/2 —
    // the variational principle for this family, in closed form. Sampling
    // alpha = w/4 must therefore reproduce N * 15w/8 and stay above E_0.
    // This checks the METROPOLIS SAMPLING of |psi|^2 (the Gaussian width
    // moments), not just the estimator.
    let omega = 1.2_f64;
    let alpha = omega / 4.0;
    let n_particles = 6;
    let closed_form = n_particles as f64 * (1.5 * alpha + 3.0 * omega * omega / (8.0 * alpha));
    let e0 = 1.5 * n_particles as f64 * omega;
    assert!((closed_form - 15.0 * n_particles as f64 * omega / 8.0).abs() < 1e-12);

    let (samples, _acceptance) = run_vmc_samples(
        0x60,
        GaussianTrap::new(alpha, [0.0; DIM]).unwrap(),
        ContinuumHamiltonian::trap_only(HarmonicTrap::new(omega, [0.0; DIM]).unwrap()).unwrap(),
        20,
        n_particles,
        1.5,
        0.8,
        200,
        4000,
    );
    let mean = mean_of(&samples);
    let stderr = batch_means_stderr(&samples, 50);

    let z = (mean - closed_form).abs() / stderr;
    assert!(
        z < 4.0,
        "E_VMC = {mean} +/- {stderr} vs closed form {closed_form} (z = {z})"
    );
    // Variational principle: above the exact ground-state energy.
    assert!(
        mean - e0 >= -5.0 * stderr,
        "E_VMC = {mean} below E_0 = {e0} (variational violation)"
    );
}
