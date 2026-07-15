//! Continuous-time Swendsen-Wang cluster construction for a longitudinal
//! spin-boson worldline.

use rand::Rng;
use rand::RngExt;

use crate::algorithm::QmcKernel;
use crate::impurity::spin_boson::cluster::retarded_bonds::LongitudinalSpinBosonModel;
use crate::impurity::spin_boson::cluster::segments::{
    build_segments_from_sorted_into, sort_and_deduplicate_times, LongitudinalWorldline,
    WorldlineSegment,
};
use crate::impurity::ImpurityError;

const DEFAULT_MAX_AUXILIARY_CUTS: usize = 1_000_000;

/// Diagnostics accumulated across continuous-time cluster sweeps.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClusterDiagnostics {
    sweeps: u64,
    inserted_cuts: u64,
    segment_count: u64,
    same_spin_pairs: u64,
    retarded_bonds: u64,
    cluster_count: u64,
    final_kinks: u64,
}

impl ClusterDiagnostics {
    /// Completed cluster sweeps.
    #[inline]
    pub const fn sweeps(&self) -> u64 {
        self.sweeps
    }

    /// Mean number of newly generated Poisson cuts.
    pub fn mean_inserted_cuts(&self) -> f64 {
        mean(self.inserted_cuts, self.sweeps)
    }

    /// Mean number of auxiliary worldline segments.
    pub fn mean_segments(&self) -> f64 {
        mean(self.segment_count, self.sweeps)
    }

    /// Fraction of considered same-spin pairs that were bonded.
    pub fn bond_fraction(&self) -> f64 {
        if self.same_spin_pairs == 0 {
            0.0
        } else {
            self.retarded_bonds as f64 / self.same_spin_pairs as f64
        }
    }

    /// Mean number of connected components before cluster orientation.
    pub fn mean_clusters(&self) -> f64 {
        mean(self.cluster_count, self.sweeps)
    }

    /// Mean real kink count after redundant auxiliary cuts are removed.
    pub fn mean_final_kinks(&self) -> f64 {
        mean(self.final_kinks, self.sweeps)
    }
}

/// Per-sweep cluster report and improved longitudinal estimators.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ClusterUpdateReport {
    pub inserted_cuts: usize,
    pub segments: usize,
    pub same_spin_pairs: usize,
    pub retarded_bonds: usize,
    pub clusters: usize,
    pub final_kinks: usize,
    /// Conditional cluster estimator of the time-averaged sampled `sigma_z`.
    pub improved_magnetization_sigma_z: f64,
    /// Conditional cluster estimator of its square.
    pub improved_m2_sigma_z: f64,
}

/// Reusable continuous-time cluster update engine.
#[derive(Debug, Clone)]
pub struct ContinuousTimeClusterEngine {
    model: LongitudinalSpinBosonModel,
    diagnostics: ClusterDiagnostics,
    last_report: ClusterUpdateReport,
    max_auxiliary_cuts: usize,
    validate_each_sweep: bool,
    cuts: Vec<f64>,
    segments: Vec<WorldlineSegment>,
    segment_spins: Vec<i8>,
    roots: Vec<usize>,
    cluster_lengths: Vec<f64>,
    cluster_orientations: Vec<i8>,
    union_find: UnionFind,
}

impl ContinuousTimeClusterEngine {
    /// Construct an engine for a longitudinal spin-boson model.
    pub fn new(model: LongitudinalSpinBosonModel) -> Self {
        Self {
            model,
            diagnostics: ClusterDiagnostics::default(),
            last_report: ClusterUpdateReport::default(),
            max_auxiliary_cuts: DEFAULT_MAX_AUXILIARY_CUTS,
            validate_each_sweep: false,
            cuts: Vec::new(),
            segments: Vec::new(),
            segment_spins: Vec::new(),
            roots: Vec::new(),
            cluster_lengths: Vec::new(),
            cluster_orientations: Vec::new(),
            union_find: UnionFind::default(),
        }
    }

    /// Physical model.
    #[inline]
    pub const fn model(&self) -> &LongitudinalSpinBosonModel {
        &self.model
    }

    /// Last completed sweep report.
    #[inline]
    pub const fn last_report(&self) -> ClusterUpdateReport {
        self.last_report
    }

    /// Cumulative update diagnostics.
    #[inline]
    pub const fn diagnostics(&self) -> &ClusterDiagnostics {
        &self.diagnostics
    }

    /// Enable expensive invariant validation after every sweep.
    pub fn set_validate_each_sweep(&mut self, enabled: bool) {
        self.validate_each_sweep = enabled;
    }

    /// Set a safety limit on the total auxiliary cut count.
    pub fn set_max_auxiliary_cuts(&mut self, maximum: usize) {
        self.max_auxiliary_cuts = maximum.max(1);
    }

    /// Execute one complete continuous-time cluster update.
    pub fn update<R: Rng + ?Sized>(
        &mut self,
        worldline: &mut LongitudinalWorldline,
        rng: &mut R,
    ) -> Result<ClusterUpdateReport, ImpurityError> {
        worldline.validate()?;
        let inserted_cuts = self.prepare_cuts(worldline, rng)?;
        build_segments_from_sorted_into(worldline, &self.cuts, &mut self.segments)?;
        self.union_find.reset(self.segments.len());

        let mut same_spin_pairs = 0usize;
        let mut retarded_bonds = 0usize;
        for (left_index, left_segment) in self.segments.iter().enumerate() {
            for (right_index, right_segment) in
                self.segments.iter().enumerate().skip(left_index + 1)
            {
                if left_segment.spin() != right_segment.spin() {
                    continue;
                }
                same_spin_pairs += 1;
                let integrated = self.model.kernel().integrated_segments(
                    worldline.beta(),
                    left_segment,
                    right_segment,
                )?;
                if !integrated.is_finite() || integrated < 0.0 {
                    return Err(ImpurityError::InvalidConfiguration(format!(
                        "retarded segment integral must be finite and non-negative, got {integrated}"
                    )));
                }
                if integrated <= 0.0 {
                    continue;
                }
                let probability = -(-2.0 * integrated).exp_m1();
                if rng.random::<f64>() < probability.min(1.0) {
                    self.union_find.union(left_index, right_index);
                    retarded_bonds += 1;
                }
            }
        }

        let (clusters, improved_magnetization, improved_m2) =
            self.assign_cluster_orientations(worldline.beta(), rng);
        self.segment_spins.clear();
        self.segment_spins.reserve(self.segments.len());
        for &root in &self.roots {
            self.segment_spins.push(self.cluster_orientations[root]);
        }
        worldline.replace_from_segments(&self.cuts, &self.segment_spins)?;
        if self.validate_each_sweep {
            worldline.validate()?;
        }

        let report = ClusterUpdateReport {
            inserted_cuts,
            segments: self.segments.len(),
            same_spin_pairs,
            retarded_bonds,
            clusters,
            final_kinks: worldline.kink_count(),
            improved_magnetization_sigma_z: improved_magnetization,
            improved_m2_sigma_z: improved_m2,
        };
        self.record(report);
        Ok(report)
    }

    fn prepare_cuts<R: Rng + ?Sized>(
        &mut self,
        worldline: &LongitudinalWorldline,
        rng: &mut R,
    ) -> Result<usize, ImpurityError> {
        self.cuts.clear();
        self.cuts.extend_from_slice(worldline.kinks());
        let rate = self.model.transverse_rate();
        let mut inserted = 0usize;
        if rate > 0.0 {
            let mut time = 0.0;
            loop {
                let uniform = rng.random::<f64>().max(f64::MIN_POSITIVE);
                time += -uniform.ln() / rate;
                if time >= worldline.beta() {
                    break;
                }
                if self.cuts.len() >= self.max_auxiliary_cuts {
                    return Err(ImpurityError::InvalidConfiguration(format!(
                        "continuous-time cluster update exceeded {} auxiliary cuts",
                        self.max_auxiliary_cuts
                    )));
                }
                self.cuts.push(time);
                inserted = inserted.saturating_add(1);
            }
        }
        sort_and_deduplicate_times(&mut self.cuts, worldline.beta());
        Ok(inserted)
    }

    fn assign_cluster_orientations<R: Rng + ?Sized>(
        &mut self,
        beta: f64,
        rng: &mut R,
    ) -> (usize, f64, f64) {
        let count = self.segments.len();
        self.roots.resize(count, 0);
        self.cluster_lengths.clear();
        self.cluster_lengths.resize(count, 0.0);
        self.cluster_orientations.clear();
        self.cluster_orientations.resize(count, 1);

        for (index, segment) in self.segments.iter().enumerate() {
            let root = self.union_find.find(index);
            self.roots[index] = root;
            self.cluster_lengths[root] += segment.length();
        }

        let mut clusters = 0usize;
        let mut conditional_integral = 0.0;
        let mut conditional_variance = 0.0;
        for (root, &length) in self.cluster_lengths.iter().enumerate() {
            if length <= 0.0 {
                continue;
            }
            clusters += 1;
            let field_argument = self.model.bias() * length;
            let probability_plus = probability_plus_spin(field_argument);
            self.cluster_orientations[root] = if rng.random::<f64>() < probability_plus {
                1
            } else {
                -1
            };

            let mean_spin = -(0.5 * field_argument).tanh();
            conditional_integral += length * mean_spin;
            conditional_variance += length * length * (1.0 - mean_spin * mean_spin);
        }

        let improved_magnetization = conditional_integral / beta;
        let improved_m2 =
            (conditional_integral * conditional_integral + conditional_variance) / (beta * beta);
        (clusters, improved_magnetization, improved_m2)
    }

    fn record(&mut self, report: ClusterUpdateReport) {
        self.last_report = report;
        self.diagnostics.sweeps = self.diagnostics.sweeps.saturating_add(1);
        self.diagnostics.inserted_cuts = self
            .diagnostics
            .inserted_cuts
            .saturating_add(report.inserted_cuts as u64);
        self.diagnostics.segment_count = self
            .diagnostics
            .segment_count
            .saturating_add(report.segments as u64);
        self.diagnostics.same_spin_pairs = self
            .diagnostics
            .same_spin_pairs
            .saturating_add(report.same_spin_pairs as u64);
        self.diagnostics.retarded_bonds = self
            .diagnostics
            .retarded_bonds
            .saturating_add(report.retarded_bonds as u64);
        self.diagnostics.cluster_count = self
            .diagnostics
            .cluster_count
            .saturating_add(report.clusters as u64);
        self.diagnostics.final_kinks = self
            .diagnostics
            .final_kinks
            .saturating_add(report.final_kinks as u64);
    }
}

impl<R> QmcKernel<LongitudinalWorldline, R> for ContinuousTimeClusterEngine
where
    R: Rng + ?Sized,
{
    type Error = ImpurityError;
    type Diagnostics = ClusterDiagnostics;

    fn sweep(
        &mut self,
        configuration: &mut LongitudinalWorldline,
        rng: &mut R,
    ) -> Result<(), Self::Error> {
        self.update(configuration, rng).map(|_| ())
    }

    fn validate(&self, configuration: &LongitudinalWorldline) -> Result<(), Self::Error> {
        configuration.validate()
    }

    fn diagnostics(&self) -> &Self::Diagnostics {
        &self.diagnostics
    }
}

#[derive(Debug, Clone, Default)]
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn reset(&mut self, count: usize) {
        self.parent.clear();
        self.parent.extend(0..count);
        self.rank.clear();
        self.rank.resize(count, 0);
    }

    fn find(&mut self, index: usize) -> usize {
        let mut root = index;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut current = index;
        while self.parent[current] != current {
            let next = self.parent[current];
            self.parent[current] = root;
            current = next;
        }
        root
    }

    fn union(&mut self, left: usize, right: usize) {
        let mut left_root = self.find(left);
        let mut right_root = self.find(right);
        if left_root == right_root {
            return;
        }
        if self.rank[left_root] < self.rank[right_root] {
            std::mem::swap(&mut left_root, &mut right_root);
        }
        self.parent[right_root] = left_root;
        if self.rank[left_root] == self.rank[right_root] {
            self.rank[left_root] = self.rank[left_root].saturating_add(1);
        }
    }
}

fn probability_plus_spin(field_argument: f64) -> f64 {
    // p(+) = 1 / (1 + exp(epsilon * cluster_length)), evaluated
    // without overflow.
    if field_argument >= 0.0 {
        let small = (-field_argument).exp();
        small / (1.0 + small)
    } else {
        1.0 / (1.0 + field_argument.exp())
    }
}

fn mean(total: u64, count: u64) -> f64 {
    if count == 0 {
        0.0
    } else {
        total as f64 / count as f64
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use crate::impurity::spin_boson::bath::{Bath, SingleModeBath};

    use super::*;

    fn model(lambda: f64, tunnelling: f64, bias: f64) -> LongitudinalSpinBosonModel {
        LongitudinalSpinBosonModel::with_default_quadrature(
            Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
            lambda,
            tunnelling,
            bias,
        )
        .expect("model")
    }

    #[test]
    fn zero_tunnelling_keeps_a_constant_worldline() {
        let mut engine = ContinuousTimeClusterEngine::new(model(0.0, 0.0, 0.0));
        let mut worldline = LongitudinalWorldline::new(5.0, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        for _ in 0..100 {
            engine.update(&mut worldline, &mut rng).expect("update");
            assert_eq!(worldline.kink_count(), 0);
        }
    }

    #[test]
    fn positive_bias_favors_negative_sigma_z() {
        let mut engine = ContinuousTimeClusterEngine::new(model(0.0, 0.0, 1.5));
        let mut worldline = LongitudinalWorldline::new(4.0, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let mut mean = 0.0;
        for _ in 0..20_000 {
            let report = engine.update(&mut worldline, &mut rng).expect("update");
            mean += worldline.integrated_sigma_z();
            assert!((report.improved_magnetization_sigma_z + (3.0f64).tanh()).abs() < 1e-14);
        }
        mean /= 20_000.0;
        assert!((mean + (3.0f64).tanh()).abs() < 0.02);
    }

    #[test]
    fn zero_field_improved_second_moment_is_cluster_length_sum() {
        let mut engine = ContinuousTimeClusterEngine::new(model(0.0, 1.0, 0.0));
        let mut worldline = LongitudinalWorldline::new(3.0, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(8);
        let report = engine.update(&mut worldline, &mut rng).expect("update");
        assert!(report.improved_magnetization_sigma_z.abs() < 1e-14);
        assert!((0.0..=1.0).contains(&report.improved_m2_sigma_z));
    }

    #[test]
    fn free_two_level_kink_count_matches_exact_mean() {
        let beta = 3.0;
        let tunnelling = 1.2;
        let gamma = 0.5 * tunnelling;
        let mut engine = ContinuousTimeClusterEngine::new(model(0.0, tunnelling, 0.0));
        let mut worldline = LongitudinalWorldline::new(beta, 1).expect("worldline");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(17);
        let warmup = 2_000;
        let samples = 40_000;
        for _ in 0..warmup {
            engine.update(&mut worldline, &mut rng).expect("update");
        }
        let mean = (0..samples)
            .map(|_| {
                engine.update(&mut worldline, &mut rng).expect("update");
                worldline.kink_count() as f64
            })
            .sum::<f64>()
            / samples as f64;
        let exact = beta * gamma * (beta * gamma).tanh();
        assert!((mean - exact).abs() < 0.04, "mean={mean}, exact={exact}");
    }
}
