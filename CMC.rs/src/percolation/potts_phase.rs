//! Q3: exploratory Potts/random-cluster phase scan across `q` at the exact
//! two-dimensional Fortuin-Kasteleyn self-dual point (roadmap v2 §6 Q3).
//!
//! # Anchor
//!
//! On the square lattice the FK random-cluster model with cluster weight
//! `q` is self-dual at
//!
//! ```text
//! p_c(q) = sqrt(q) / (1 + sqrt(q)),
//! ```
//!
//! the critical curve of the `q`-state Potts model (Fortuin-Kasteleyn 1972;
//! Wu, *The Potts model*, Rev. Mod. Phys. 54, 235 (1982); the same curve is
//! quoted in the archived random-cluster literature). `q = 4` is the
//! marginal case (continuous), `q > 4` is first order, `1 <= q < 4` is
//! continuous. This module uses the EXACT curve as an anchor: it never fits
//! a critical point.
//!
//! # What is measured
//!
//! [`SweenyMC`] runs at `(p_c(q), q)` on a periodic `L x L` square. The
//! recorded order parameter proxy is the open-bond density `|A| / E`
//! (`e = -|A|` for the Potts energy, so a latent-heat jump is a gap in this
//! density). [`phase_scan`] reports its mean, variance and the
//! **bimodality coefficient**
//!
//! ```text
//! b = (skewness^2 + 1) / kurtosis,
//! ```
//!
//! with `b > 5/9` the usual descriptive indication of a two-mode
//! distribution. This is an EXPLORATORY finite-size diagnostic, not a
//! latent-heat estimator and not a finite-size-scaling proof: at the `L`
//! reachable in this stage the two peaks may overlap. The validation gate
//! checks the exact critical curve and the machinery; the `q`-ordering of
//! `b` is reported, and no asymptotic claim is made.
//!
//! # Cost and determinism
//!
//! One scan constructs one [`SweenyMC`] and runs `thermalization + samples`
//! full sweeps (one proposal per physical edge, ascending order, one
//! uniform per proposal). The seed is explicit; the same seed and sizes
//! reproduce the scan bitwise. Construction allocates once; sweeps are
//! steady-state allocation-free because the Sweeny workspace is reserved.

use core::fmt;

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::sweeny::{SweenyAcceptance, SweenyError, SweenyMC};
use super::{ObservablePlan, Probability};
use crate::{build_square, UndirectedGraphView};

/// The exact 2D FK/Potts self-dual bond probability `sqrt(q)/(1+sqrt(q))`.
pub fn critical_probability(q: f64) -> Result<f64, PottsPhaseError> {
    if !q.is_finite() || q <= 0.0 {
        return Err(PottsPhaseError::ClusterWeight { q });
    }
    let root = q.sqrt();
    Ok(root / (1.0 + root))
}

/// One exploratory phase-scan configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PottsPhaseSpec {
    /// Cluster weight `q > 0`.
    pub q: f64,
    /// Periodic square side `L`; the scan uses `L x L` vertices.
    pub side: usize,
    /// Sweeps discarded before measurement.
    pub thermalization: u64,
    /// Measured sweeps.
    pub samples: u64,
    /// RNG seed.
    pub seed: u64,
}

impl PottsPhaseSpec {
    /// Validate a scan specification.
    pub fn new(
        q: f64,
        side: usize,
        thermalization: u64,
        samples: u64,
        seed: u64,
    ) -> Result<Self, PottsPhaseError> {
        critical_probability(q)?;
        if !(2..=256).contains(&side) {
            return Err(PottsPhaseError::Side { side });
        }
        if samples < 2 {
            return Err(PottsPhaseError::Samples { samples });
        }
        Ok(Self {
            q,
            side,
            thermalization,
            samples,
            seed,
        })
    }
}

/// One scan's descriptive statistics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PottsPhaseEstimate {
    /// Cluster weight.
    pub q: f64,
    /// Exact self-dual probability used.
    pub critical_probability: f64,
    /// Periodic side.
    pub side: usize,
    /// Measured sweeps.
    pub samples: u64,
    /// Mean open-bond density `|A|/E`.
    pub density_mean: f64,
    /// Unbiased variance of the density.
    pub density_variance: f64,
    /// Skewness of the density.
    pub skewness: f64,
    /// Excess kurtosis is NOT used; this is the raw fourth standardized
    /// moment `E[(x-mu)^4]/sigma^4` (the coefficient's denominator).
    pub kurtosis: f64,
    /// Exploratory bimodality coefficient `(skewness^2 + 1)/kurtosis`.
    /// Values above `5/9` are the usual descriptive bimodality indication.
    pub bimodality: f64,
}

/// Run one exploratory phase scan at the exact `p_c(q)`.
pub fn phase_scan(spec: PottsPhaseSpec) -> Result<PottsPhaseEstimate, PottsPhaseError> {
    Ok(phase_scan_samples(spec)?.estimate)
}

/// Raw bond-density samples of one scan plus the descriptive estimate.
#[derive(Clone, Debug, PartialEq)]
pub struct PottsPhaseSamples {
    pub estimate: PottsPhaseEstimate,
    /// Open-bond density `|A|/E` per measured sweep, in draw order.
    pub densities: Vec<f64>,
}

/// Run one scan and return the raw density samples.
pub fn phase_scan_samples(spec: PottsPhaseSpec) -> Result<PottsPhaseSamples, PottsPhaseError> {
    let probability = critical_probability(spec.q)?;
    let probability = Probability::new(probability).map_err(PottsPhaseError::Probability)?;
    let lattice = build_square(spec.side, spec.side, true);
    let mut chain = SweenyMC::new(
        lattice.clone(),
        spec.q,
        probability.get(),
        SweenyAcceptance::HeatBath,
        ObservablePlan::all(Vec::new()),
    )
    .map_err(PottsPhaseError::Sweeny)?;
    let edges = lattice.edge_count() as f64;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(spec.seed);
    for _ in 0..spec.thermalization {
        sweep(&mut chain, &mut rng);
    }
    let count = spec.samples as usize;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| PottsPhaseError::Capacity)?;
    for _ in 0..spec.samples {
        sweep(&mut chain, &mut rng);
        values.push(chain.open_edge_count() as f64 / edges);
    }
    let mean = values.iter().sum::<f64>() / count as f64;
    let central: Vec<f64> = values.iter().map(|value| value - mean).collect();
    let variance = central.iter().map(|value| value * value).sum::<f64>() / (count - 1) as f64;
    if variance <= 0.0 {
        return Err(PottsPhaseError::DegenerateVariance);
    }
    let sigma = variance.sqrt();
    let skewness = central
        .iter()
        .map(|value| (value / sigma).powi(3))
        .sum::<f64>()
        / count as f64;
    let kurtosis = central
        .iter()
        .map(|value| (value / sigma).powi(4))
        .sum::<f64>()
        / count as f64;
    let bimodality = (skewness * skewness + 1.0) / kurtosis;
    Ok(PottsPhaseSamples {
        estimate: PottsPhaseEstimate {
            q: spec.q,
            critical_probability: probability.get(),
            side: spec.side,
            samples: spec.samples,
            density_mean: mean,
            density_variance: variance,
            skewness,
            kurtosis,
            bimodality,
        },
        densities: values,
    })
}

/// Two-mode (latent-heat proxy) decomposition of density samples.
///
/// A deterministic one-dimensional two-means clustering splits the samples
/// into a lower and an upper mode. `gap` is the distance between the mode
/// centers and `separation` divides it by the pooled within-mode standard
/// deviation, so overlapping peaks score low. This is a definitional
/// estimator, not a thermodynamic latent-heat measurement: the finite-size
/// bracket is recorded, and a large `separation` alone does not certify a
/// first-order transition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoModeEstimate {
    /// Lower mode center.
    pub lower: f64,
    /// Upper mode center.
    pub upper: f64,
    /// `upper - lower`.
    pub gap: f64,
    /// `gap / sqrt(pooled within-mode variance)` (`+inf` for exactly
    /// separated constant modes).
    pub separation: f64,
    /// Fraction of samples in the lower mode.
    pub lower_weight: f64,
    /// Midpoint between the two mode centers.
    pub threshold: f64,
}

/// Split density samples into lower/upper modes with deterministic 2-means.
pub fn two_mode_estimate(samples: &[f64]) -> Result<TwoModeEstimate, PottsPhaseError> {
    if samples.len() < 4 {
        return Err(PottsPhaseError::TooFewSamples {
            samples: samples.len(),
        });
    }
    if let Some(&value) = samples.iter().find(|value| !value.is_finite()) {
        return Err(PottsPhaseError::NonFiniteSample { value });
    }
    let mut sorted = samples.to_vec();
    sorted.sort_by(f64::total_cmp);
    let minimum = sorted[0];
    let maximum = sorted[sorted.len() - 1];
    if maximum <= minimum {
        return Err(PottsPhaseError::DegenerateVariance);
    }
    let mut lower = sorted[sorted.len() / 4];
    let mut upper = sorted[3 * sorted.len() / 4];
    if lower == upper {
        lower = minimum;
        upper = maximum;
    }
    for _ in 0..100 {
        let mut lower_sum = 0.0;
        let mut upper_sum = 0.0;
        let mut lower_count = 0usize;
        let mut upper_count = 0usize;
        let threshold = 0.5 * (lower + upper);
        for &value in &sorted {
            if value <= threshold {
                lower_sum += value;
                lower_count += 1;
            } else {
                upper_sum += value;
                upper_count += 1;
            }
        }
        if lower_count == 0 || upper_count == 0 {
            return Err(PottsPhaseError::DegenerateTwoMode);
        }
        let next_lower = lower_sum / lower_count as f64;
        let next_upper = upper_sum / upper_count as f64;
        let moved = (next_lower - lower).abs() + (next_upper - upper).abs();
        lower = next_lower;
        upper = next_upper;
        if moved < 1e-12 {
            break;
        }
    }
    let threshold = 0.5 * (lower + upper);
    let lower_count = sorted.iter().filter(|&&value| value <= threshold).count();
    if lower_count == 0 || lower_count == sorted.len() {
        return Err(PottsPhaseError::DegenerateTwoMode);
    }
    let lower_weight = lower_count as f64 / sorted.len() as f64;
    let mut within = 0.0;
    for &value in &sorted {
        let center = if value <= threshold { lower } else { upper };
        within += (value - center) * (value - center);
    }
    within /= sorted.len() as f64;
    let gap = upper - lower;
    let separation = if within > 0.0 {
        gap / within.sqrt()
    } else {
        f64::INFINITY
    };
    Ok(TwoModeEstimate {
        lower,
        upper,
        gap,
        separation,
        lower_weight,
        threshold,
    })
}

fn sweep(chain: &mut SweenyMC, rng: &mut Xoshiro256PlusPlus) {
    for edge in 0..chain.lattice().edge_count() {
        chain
            .propose(edge, rng)
            .expect("edge index is in the topology domain");
    }
}

/// Typed rejections of the phase scan.
#[derive(Debug)]
pub enum PottsPhaseError {
    /// Cluster weight outside `(0, infinity)`.
    ClusterWeight { q: f64 },
    /// Side outside `[2, 256]`.
    Side { side: usize },
    /// Fewer than two measured sweeps.
    Samples { samples: u64 },
    /// The exact critical probability failed the probability domain.
    Probability(super::ProbabilityError),
    /// The Sweeny constructor rejected the configuration.
    Sweeny(SweenyError),
    /// The measured density has zero variance.
    DegenerateVariance,
    /// Fewer than four samples cannot support a two-mode split.
    TooFewSamples { samples: usize },
    /// A sampled density is not finite.
    NonFiniteSample { value: f64 },
    /// The two-means split emptied one mode.
    DegenerateTwoMode,
    /// Result storage could not be reserved.
    Capacity,
}

impl fmt::Display for PottsPhaseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ClusterWeight { q } => {
                write!(
                    formatter,
                    "cluster weight q={q} must be finite and positive"
                )
            }
            Self::Side { side } => {
                write!(formatter, "side {side} must lie in [2, 256]")
            }
            Self::Samples { samples } => {
                write!(
                    formatter,
                    "at least two measured sweeps are required, got {samples}"
                )
            }
            Self::Probability(error) => write!(formatter, "critical probability rejected: {error}"),
            Self::Sweeny(error) => write!(formatter, "Sweeny chain rejected: {error}"),
            Self::DegenerateVariance => write!(formatter, "the measured density has zero variance"),
            Self::TooFewSamples { samples } => write!(
                formatter,
                "two-mode estimation needs at least four samples, got {samples}"
            ),
            Self::NonFiniteSample { value } => {
                write!(formatter, "density sample {value} is not finite")
            }
            Self::DegenerateTwoMode => {
                write!(formatter, "the two-means split emptied one mode")
            }
            Self::Capacity => write!(formatter, "phase-scan storage could not be reserved"),
        }
    }
}

impl std::error::Error for PottsPhaseError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_mode_recovers_synthetic_peaks() {
        let mut samples = Vec::new();
        for index in 0..300 {
            samples.push(0.2 + 0.001 * (index % 5) as f64);
            samples.push(0.8 + 0.001 * (index % 5) as f64);
        }
        let estimate = two_mode_estimate(&samples).expect("two modes");
        assert!((estimate.lower - 0.202).abs() < 0.01);
        assert!((estimate.upper - 0.802).abs() < 0.01);
        assert!(estimate.gap > 0.55);
        assert!(estimate.separation > 10.0);
        assert!((estimate.lower_weight - 0.5).abs() < 0.01);
        assert!(estimate.threshold > estimate.lower && estimate.threshold < estimate.upper);
    }

    #[test]
    fn uniform_samples_are_not_strongly_separated() {
        let samples: Vec<f64> = (0..=1000).map(|index| index as f64 / 1000.0).collect();
        let estimate = two_mode_estimate(&samples).expect("split");
        assert!(
            estimate.separation < 4.0,
            "uniform separation {:.3} is not the bimodal regime",
            estimate.separation
        );
    }

    #[test]
    fn two_mode_rejections_are_typed() {
        assert!(matches!(
            two_mode_estimate(&[0.1, 0.2, 0.3]),
            Err(PottsPhaseError::TooFewSamples { samples: 3 })
        ));
        assert!(matches!(
            two_mode_estimate(&[0.1, 0.1, f64::NAN, 0.2]),
            Err(PottsPhaseError::NonFiniteSample { .. })
        ));
        assert!(matches!(
            two_mode_estimate(&[0.5, 0.5, 0.5, 0.5]),
            Err(PottsPhaseError::DegenerateVariance)
        ));
    }
}
