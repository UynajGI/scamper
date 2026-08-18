//! Equilibrium-distribution validation for the rigid-molecule solver
//! (`MolecularMetropolisCore`).
//!
//! VALIDATION.md lists the rigid-molecule kernel as experimental with
//! "NOT validated: Equilibrium distribution" — only rigid-geometry
//! preservation was covered. These tests sample the real solver
//! (translation + rotation moves, Metropolis acceptance, cell-list caches)
//! and compare equilibrium observables against exact Boltzmann references.
//!
//! ## Why quadrature references instead of a Langevin curve
//!
//! `PairPotential::energy(species_i, species_j, distance_squared)` is the
//! only energy channel in the particle API: there is no external-field or
//! harmonic-trap handle, so a rigid dipole in a uniform field
//! (⟨cosθ⟩ = coth x − 1/x) cannot be expressed. The field is therefore
//! replaced by a second molecule coupled through a test-local truncated
//! Gaussian well. For tiny systems the full Boltzmann integral then
//! reduces — by translation and rotation invariance of the periodic cell —
//! to low-dimensional quadrature evaluated in-process:
//!
//! 1. **Pair** (two single-atom molecules): translation sector only. The
//!    pair partition function is exactly 1D (radial Simpson quadrature);
//!    checks ⟨U⟩ and the bound-state probability.
//! 2. **Probe** (one rigid dumbbell + one atom): joint COM/orientation
//!    equilibrium of a rotor relative to a point ligand. The 5D integral
//!    reduces to a 2D midpoint quadrature in the rotation gauge θ = 0;
//!    checks the nematic probe-axis alignment ⟨cos 2α⟩ and ⟨U⟩.
//! 3. **Rotor pair** (two rigid dumbbells): mutual orientational locking.
//!    ⟨cos 2Δθ⟩ is invariant under global rotation, so translation moves
//!    cannot change it — it directly validates the *rotation* move's
//!    equilibrium across a coupling sweep ε (the analog of the Langevin
//!    field sweep x, from the linear-response regime ε → 0 up to a locked
//!    state). The 8D integral reduces to a 3D quadrature (relative
//!    displacement s and Δθ); the reference is validated separately by
//!    grid refinement.
//!
//! The well parameters keep every equilibrium feature inside the
//! minimum-image-safe region (support radius < L/2), so the quadrature
//! needs no wrap handling, and keep escape barriers ≈ ε so chains mix
//! without tunneling problems at β = 1.
//!
//! Statistical standard: `zscore_seed_count(8)` independent chains per
//! coupling; pooled z = |mean − exact| / (σ_seed/√n) and per-seed
//! z_i = |mean_i − exact| / σ_seed must both stay under 4, with σ_seed
//! from the seed spread. `SCUTTLE_ZSCORE_SEEDS` scales the count (nightly).

use super::common::zscore_seed_count;
use cmc_rs::{
    MolecularMetropolisCore, MoleculeTopology, OrthorhombicCell, PairPotential, ParticleAlgorithm,
    ParticleConfiguration, ParticleSystem, SimulationCell, SimulationPhase,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const BETA: f64 = 1.0;

// ── Test-local pair potential ─────────────────────────────────────────────

/// Truncated attractive-Gaussian well acting only between unlike species.
///
/// u(r) = -ε·exp(-(r - r_well)² / (2·width²)) for r < cutoff, else 0.
/// Like-species pairs (the intramolecular bond of a dumbbell) see zero, so
/// rigid constraints carry no residual intramolecular energy and the exact
/// references below only integrate cross terms.
#[derive(Debug, Clone, Copy)]
struct CrossWell {
    epsilon: f64,
    r_well: f64,
    cutoff_squared: f64,
    inverse_two_variance: f64,
}

impl CrossWell {
    fn new(epsilon: f64, r_well: f64, width: f64, cutoff: f64) -> Self {
        Self {
            epsilon,
            r_well,
            cutoff_squared: cutoff * cutoff,
            inverse_two_variance: 0.5 / (width * width),
        }
    }

    #[inline]
    fn value_at_distance(&self, distance: f64) -> f64 {
        let offset = distance - self.r_well;
        -self.epsilon * (-offset * offset * self.inverse_two_variance).exp()
    }

    #[inline]
    fn value_at_distance_squared(&self, distance_squared: f64) -> f64 {
        if distance_squared >= self.cutoff_squared {
            0.0
        } else {
            self.value_at_distance(distance_squared.sqrt())
        }
    }
}

impl PairPotential for CrossWell {
    fn cutoff_squared(&self) -> f64 {
        self.cutoff_squared
    }

    fn energy(&self, species_i: u16, species_j: u16, distance_squared: f64) -> f64 {
        if species_i == species_j {
            0.0
        } else {
            self.value_at_distance_squared(distance_squared)
        }
    }
}

// ── Chain driver and statistics ───────────────────────────────────────────

/// Fixed sweep plan for one Markov chain.
struct ChainPlan {
    seed: u64,
    thermalization: u64,
    measurement: u64,
    max_displacement: f64,
    max_angle: f64,
}

/// Thermalize, then measure: call `observe` after every measurement sweep
/// of the real rigid-molecule kernel, and audit the state at the end.
fn run_chain(
    plan: &ChainPlan,
    mut system: ParticleSystem<2>,
    topology: MoleculeTopology,
    potential: &CrossWell,
    mut observe: impl FnMut(&ParticleSystem<2>),
) {
    let mut kernel = MolecularMetropolisCore::new(topology, plan.max_displacement, plan.max_angle)
        .expect("rigid-molecule kernel parameters are valid");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(plan.seed);
    for _ in 0..plan.thermalization {
        kernel.sweep_with_phase(
            &mut system,
            potential,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }
    for _ in 0..plan.measurement {
        kernel.sweep_with_phase(
            &mut system,
            potential,
            &mut rng,
            SimulationPhase::Measurement,
        );
        observe(&system);
    }
    system
        .validate(potential)
        .expect("rigid-molecule chain must keep energy and cell-list caches consistent");
}

/// Multi-seed z-score assertion against an exact reference.
///
/// Requires pooled |mean − exact| ≤ 4·σ_seed/√n (or the absolute `floor`,
/// which only rescues degenerate seed-count overrides like
/// SCUTTLE_ZSCORE_SEEDS=1 where the spread estimate is meaningless) and
/// per-seed |z_i| < 4 against the seed spread.
fn assert_seeds_match_exact(per_seed: &[f64], exact: f64, floor: f64, label: &str) {
    let count = per_seed.len();
    let mean = per_seed.iter().sum::<f64>() / count as f64;
    let deviation = (mean - exact).abs();
    if count < 3 {
        assert!(
            deviation <= floor,
            "{label}: single-chain deviation {deviation:.5} exceeds floor {floor:.5}"
        );
        return;
    }
    let variance = per_seed
        .iter()
        .map(|&value| (value - mean).powi(2))
        .sum::<f64>()
        / (count - 1) as f64;
    let spread = variance.sqrt();
    let stderr = spread / (count as f64).sqrt();
    eprintln!(
        "[molecule-equilibrium] {label}: exact {exact:+.5}, sampled {mean:+.5} ± {stderr:.5} \
         (σ_seed {spread:.5}, n {count}, pooled z {:+.2})",
        (mean - exact) / stderr.max(1e-15)
    );
    assert!(
        deviation <= 4.0 * stderr || deviation <= floor,
        "{label}: pooled mean {mean:.5} vs exact {exact:.5}: z = {:.2} \
         (σ_seed {spread:.5}, n {count})",
        (mean - exact) / stderr.max(1e-15)
    );
    for (index, &value) in per_seed.iter().enumerate() {
        let z = (value - exact) / spread.max(1e-15);
        assert!(
            z.abs() < 4.0,
            "{label}: seed {index} z = {z:.2} \
             (mean {value:.5}, exact {exact:.5}, σ_seed {spread:.5})"
        );
    }
}

fn simpson_integral(a: f64, b: f64, intervals: usize, integrand: impl Fn(f64) -> f64) -> f64 {
    assert!(
        intervals.is_multiple_of(2),
        "Simpson rule needs an even interval count"
    );
    let step = (b - a) / intervals as f64;
    let mut total = integrand(a) + integrand(b);
    for index in 1..intervals {
        let weight = if index.is_multiple_of(2) { 2.0 } else { 4.0 };
        total += weight * integrand(a + index as f64 * step);
    }
    total * step / 3.0
}

// ── Case 1: two single-atom molecules (translation sector) ────────────────

const PAIR_LENGTH: f64 = 6.0;
const PAIR_R_WELL: f64 = 0.9;
const PAIR_WIDTH: f64 = 0.15;
const PAIR_CUTOFF: f64 = 1.5;
const PAIR_SEED_BASE: u64 = 0x0C0F_FFE0;
const PAIR_THERMALIZATION: u64 = 3_000;
const PAIR_MEASUREMENT: u64 = 12_000;
const PAIR_ENERGY_FLOOR: f64 = 0.12;
const PAIR_BOUND_FLOOR: f64 = 0.08;

fn pair_system(epsilon: f64) -> (ParticleSystem<2>, CrossWell, MoleculeTopology) {
    let cell = OrthorhombicCell::new([PAIR_LENGTH; 2]).expect("valid cell lengths");
    let configuration = ParticleConfiguration::new(vec![[1.5, 3.0], [4.5, 3.0]], vec![0, 1], cell)
        .expect("valid initial configuration");
    let potential = CrossWell::new(epsilon, PAIR_R_WELL, PAIR_WIDTH, PAIR_CUTOFF);
    let topology = MoleculeTopology::new(2, vec![vec![0], vec![1]])
        .expect("one atom per molecule is a valid topology");
    (
        ParticleSystem::new(configuration, &potential, BETA).expect("valid initial pair system"),
        potential,
        topology,
    )
}

/// Exact ⟨U⟩ and bound fraction P(r < cutoff) by radial Simpson quadrature.
///
/// Translational invariance reduces the two-particle partition function to
/// the minimum-image relative coordinate; because cutoff < L/2 every point
/// with |s| < cutoff has minimum-image distance |s|, so
/// Z_rel = L² + 2π∫₀^{rc}(e^{-βu} − 1) r dr exactly.
fn pair_reference(epsilon: f64, intervals: usize) -> (f64, f64) {
    let well = CrossWell::new(epsilon, PAIR_R_WELL, PAIR_WIDTH, PAIR_CUTOFF);
    let delta_integral = simpson_integral(0.0, PAIR_CUTOFF, intervals, |r| {
        ((-BETA * well.value_at_distance(r)).exp() - 1.0) * r
    });
    let energy_integral = simpson_integral(0.0, PAIR_CUTOFF, intervals, |r| {
        well.value_at_distance(r) * (-BETA * well.value_at_distance(r)).exp() * r
    });
    let partition = PAIR_LENGTH * PAIR_LENGTH + 2.0 * std::f64::consts::PI * delta_integral;
    let mean_energy = 2.0 * std::f64::consts::PI * energy_integral / partition;
    let bound = (std::f64::consts::PI * PAIR_CUTOFF * PAIR_CUTOFF
        + 2.0 * std::f64::consts::PI * delta_integral)
        / partition;
    (mean_energy, bound)
}

#[test]
fn pair_translation_equilibrium_matches_radial_quadrature() {
    let epsilon = 2.0;
    let (exact_energy, exact_bound) = pair_reference(epsilon, 4096);
    let seeds = zscore_seed_count(8);
    let mut energy_means = Vec::with_capacity(seeds);
    let mut bound_means = Vec::with_capacity(seeds);
    for seed_index in 0..seeds {
        let (system, potential, topology) = pair_system(epsilon);
        let plan = ChainPlan {
            seed: PAIR_SEED_BASE + seed_index as u64,
            thermalization: PAIR_THERMALIZATION,
            measurement: PAIR_MEASUREMENT,
            max_displacement: 0.35,
            max_angle: 0.5,
        };
        let mut energy_sum = 0.0;
        let mut bound_sum = 0.0;
        let mut sweeps = 0u64;
        run_chain(&plan, system, topology, &potential, |chain| {
            let cell = *chain.configuration().cell();
            let positions = chain.configuration().positions();
            let distance_squared = cell.distance_squared(&positions[0], &positions[1]);
            energy_sum += chain.energy;
            bound_sum += if distance_squared < potential.cutoff_squared {
                1.0
            } else {
                0.0
            };
            sweeps += 1;
        });
        energy_means.push(energy_sum / sweeps as f64);
        bound_means.push(bound_sum / sweeps as f64);
    }
    assert_seeds_match_exact(
        &energy_means,
        exact_energy,
        PAIR_ENERGY_FLOOR,
        "pair ⟨U⟩ (ε=2.0)",
    );
    assert_seeds_match_exact(
        &bound_means,
        exact_bound,
        PAIR_BOUND_FLOOR,
        "pair bound fraction (ε=2.0)",
    );
}

// ── Case 2: rigid dumbbell + probe atom (joint equilibrium) ───────────────

const PROBE_LENGTH: f64 = 6.0;
const PROBE_BOND: f64 = 1.0;
const PROBE_R_WELL: f64 = 1.0;
const PROBE_WIDTH: f64 = 0.15;
const PROBE_CUTOFF: f64 = 1.8;
const PROBE_SEED_BASE: u64 = 0x0C0F_FFB1;
const PROBE_THERMALIZATION: u64 = 3_000;
const PROBE_MEASUREMENT: u64 = 12_000;
const PROBE_ORDER_FLOOR: f64 = 0.12;
const PROBE_ENERGY_FLOOR: f64 = 0.5;

fn probe_system(epsilon: f64) -> (ParticleSystem<2>, CrossWell, MoleculeTopology) {
    let half = 0.5 * PROBE_BOND;
    let axis = [0.7_f64.cos(), 0.7_f64.sin()];
    let dumbbell_center = [3.0, 3.0];
    let positions = vec![
        [
            dumbbell_center[0] - half * axis[0],
            dumbbell_center[1] - half * axis[1],
        ],
        [
            dumbbell_center[0] + half * axis[0],
            dumbbell_center[1] + half * axis[1],
        ],
        [1.2, 4.8],
    ];
    let cell = OrthorhombicCell::new([PROBE_LENGTH; 2]).expect("valid cell lengths");
    let configuration =
        ParticleConfiguration::new(positions, vec![0, 0, 1], cell).expect("valid configuration");
    let potential = CrossWell::new(epsilon, PROBE_R_WELL, PROBE_WIDTH, PROBE_CUTOFF);
    let topology = MoleculeTopology::new(3, vec![vec![0, 1], vec![2]])
        .expect("dumbbell + monomer is a valid topology");
    (
        ParticleSystem::new(configuration, &potential, BETA).expect("valid initial probe system"),
        potential,
        topology,
    )
}

/// Exact ⟨cos 2α⟩ and ⟨U⟩ (α = angle between the bond axis and the
/// probe-to-COM direction) by 2D midpoint quadrature in the gauge θ = 0.
///
/// Global-rotation invariance fixes the dumbbell along x; the integral
/// then runs over the relative probe coordinate s alone. The grid square
/// covers the whole interaction support (|s| ≤ bond/2 + cutoff < L/2, so
/// no minimum-image handling is needed); outside it U = 0 exactly.
/// Because grid square and box are both centered squares, the remaining
/// domain is invariant under 90° rotations, giving ∫ m₂ = area/2 there.
fn probe_reference(epsilon: f64, step: f64) -> (f64, f64) {
    let well = CrossWell::new(epsilon, PROBE_R_WELL, PROBE_WIDTH, PROBE_CUTOFF);
    let half = 0.5 * PROBE_BOND;
    let support = half + PROBE_CUTOFF;
    let points = ((2.0 * support) / step).ceil() as usize;
    let spacing = 2.0 * support / points as f64;
    let cell_area = spacing * spacing;
    let mut weighted = 0.0;
    let mut energy_weighted = 0.0;
    let mut alignment_weighted = 0.0;
    let mut alignment_plain = 0.0;
    for index_x in 0..points {
        let x = -support + (index_x as f64 + 0.5) * spacing;
        for index_y in 0..points {
            let y = -support + (index_y as f64 + 0.5) * spacing;
            let energy = well.value_at_distance((x - half).hypot(y))
                + well.value_at_distance((x + half).hypot(y));
            let boltzmann = (-BETA * energy).exp();
            let alignment = x * x / (x * x + y * y);
            weighted += boltzmann * cell_area;
            energy_weighted += energy * boltzmann * cell_area;
            alignment_weighted += alignment * boltzmann * cell_area;
            alignment_plain += alignment * cell_area;
        }
    }
    let grid_area = (points as f64 * spacing).powi(2);
    let partition = weighted + PROBE_LENGTH * PROBE_LENGTH - grid_area;
    let mean_alignment =
        (alignment_weighted + PROBE_LENGTH * PROBE_LENGTH / 2.0 - alignment_plain) / partition;
    (2.0 * mean_alignment - 1.0, energy_weighted / partition)
}

#[test]
fn probe_dumbbell_joint_equilibrium_matches_2d_quadrature() {
    for (coupling_index, epsilon) in [2.0_f64, 4.0].iter().enumerate() {
        let (exact_order, exact_energy) = probe_reference(*epsilon, 0.02);
        let seeds = zscore_seed_count(8);
        let mut order_means = Vec::with_capacity(seeds);
        let mut energy_means = Vec::with_capacity(seeds);
        for seed_index in 0..seeds {
            let (system, potential, topology) = probe_system(*epsilon);
            let plan = ChainPlan {
                seed: PROBE_SEED_BASE + 100_000 * coupling_index as u64 + seed_index as u64,
                thermalization: PROBE_THERMALIZATION,
                measurement: PROBE_MEASUREMENT,
                max_displacement: 0.3,
                max_angle: 0.35,
            };
            let mut order_sum = 0.0;
            let mut energy_sum = 0.0;
            let mut sweeps = 0u64;
            run_chain(&plan, system, topology, &potential, |chain| {
                let cell = *chain.configuration().cell();
                let positions = chain.configuration().positions();
                let bond = cell.displacement(&positions[0], &positions[1]);
                let bond_norm = (bond[0] * bond[0] + bond[1] * bond[1]).sqrt();
                let axis = [bond[0] / bond_norm, bond[1] / bond_norm];
                let center = [
                    positions[0][0] + 0.5 * bond[0],
                    positions[0][1] + 0.5 * bond[1],
                ];
                let relative = cell.displacement(&center, &positions[2]);
                let norm_squared = relative[0] * relative[0] + relative[1] * relative[1];
                let alignment = if norm_squared < 1e-24 {
                    0.5
                } else {
                    let projection = (relative[0] * axis[0] + relative[1] * axis[1]).powi(2);
                    projection / norm_squared
                };
                order_sum += 2.0 * alignment - 1.0;
                energy_sum += chain.energy;
                sweeps += 1;
            });
            order_means.push(order_sum / sweeps as f64);
            energy_means.push(energy_sum / sweeps as f64);
        }
        assert_seeds_match_exact(
            &order_means,
            exact_order,
            PROBE_ORDER_FLOOR,
            &format!("probe ⟨cos2α⟩ (ε={epsilon})"),
        );
        assert_seeds_match_exact(
            &energy_means,
            exact_energy,
            PROBE_ENERGY_FLOOR,
            &format!("probe ⟨U⟩ (ε={epsilon})"),
        );
    }
}

// ── Case 3: two rigid dumbbells (orientational locking) ───────────────────

const ROTOR_LENGTH: f64 = 6.2;
const ROTOR_BOND: f64 = 1.6;
const ROTOR_R_WELL: f64 = 0.7;
const ROTOR_WIDTH: f64 = 0.15;
const ROTOR_CUTOFF: f64 = 1.3;
const ROTOR_SEED_BASE: u64 = 0x0C0F_FF72;
const ROTOR_THERMALIZATION: u64 = 5_000;
const ROTOR_MEASUREMENT: u64 = 12_000;
const ROTOR_ORDER_FLOOR: f64 = 0.12;
const ROTOR_ENERGY_FLOOR: f64 = 0.8;

fn rotor_system(epsilon: f64) -> (ParticleSystem<2>, CrossWell, MoleculeTopology) {
    let half = 0.5 * ROTOR_BOND;
    let axis_zero = [0.3_f64.cos(), 0.3_f64.sin()];
    let axis_one = [1.5_f64.cos(), 1.5_f64.sin()];
    let center_zero = [3.1, 3.1];
    let center_one = [1.6, 4.9];
    let positions = vec![
        [
            center_zero[0] - half * axis_zero[0],
            center_zero[1] - half * axis_zero[1],
        ],
        [
            center_zero[0] + half * axis_zero[0],
            center_zero[1] + half * axis_zero[1],
        ],
        [
            center_one[0] - half * axis_one[0],
            center_one[1] - half * axis_one[1],
        ],
        [
            center_one[0] + half * axis_one[0],
            center_one[1] + half * axis_one[1],
        ],
    ];
    let cell = OrthorhombicCell::new([ROTOR_LENGTH; 2]).expect("valid cell lengths");
    let configuration =
        ParticleConfiguration::new(positions, vec![0, 0, 1, 1], cell).expect("valid configuration");
    let potential = CrossWell::new(epsilon, ROTOR_R_WELL, ROTOR_WIDTH, ROTOR_CUTOFF);
    let topology = MoleculeTopology::new(4, vec![vec![0, 1], vec![2, 3]])
        .expect("two dumbbells form a valid topology");
    (
        ParticleSystem::new(configuration, &potential, BETA).expect("valid initial rotor system"),
        potential,
        topology,
    )
}

/// Exact ⟨cos 2Δθ⟩ and ⟨U⟩ for two dumbbells by 3D quadrature: a 2D
/// midpoint grid over the relative displacement s (support radius
/// bond + cutoff < L/2) times Simpson over Δθ ∈ [0, π/2].
///
/// The integrand is π-periodic (swapping a dumbbell's atoms is the same
/// rigid state) and even in Δθ, so the [−π, π] integral is 4× the
/// quarter-interval. Global-rotation invariance sets one axis along x.
/// Ã(Δθ) = ∫(e^{-βU} − 1) d²s vanishes outside the support, so the free
/// part of the partition function (2πL²) is exact, and cos 2Δθ averages
/// to zero against it — the numerator only receives support weight.
fn rotor_reference(epsilon: f64, step: f64) -> (f64, f64) {
    let well = CrossWell::new(epsilon, ROTOR_R_WELL, ROTOR_WIDTH, ROTOR_CUTOFF);
    let half = 0.5 * ROTOR_BOND;
    let support = ROTOR_BOND + ROTOR_CUTOFF;
    let points = ((2.0 * support) / step).ceil() as usize;
    let spacing = 2.0 * support / points as f64;
    let cell_area = spacing * spacing;
    let intervals = 32;
    let step_theta = std::f64::consts::FRAC_PI_2 / intervals as f64;
    let mut integral_a = 0.0;
    let mut integral_b = 0.0;
    let mut integral_numerator = 0.0;
    for index in 0..=intervals {
        let theta = index as f64 * step_theta;
        let axis_one = [half * theta.cos(), half * theta.sin()];
        let mut integrand_a = 0.0;
        let mut integrand_b = 0.0;
        for index_x in 0..points {
            let x = -support + (index_x as f64 + 0.5) * spacing;
            for index_y in 0..points {
                let y = -support + (index_y as f64 + 0.5) * spacing;
                let mut energy = 0.0;
                let mut interacting = false;
                for sign_axis in [-1.0, 1.0] {
                    for sign_other in [-1.0, 1.0] {
                        // Pair offsets: tau*a0 - sigma*a1 with a0 = (half, 0).
                        let offset_x = sign_axis * half - sign_other * axis_one[0];
                        let offset_y = -sign_other * axis_one[1];
                        let distance_squared =
                            (x - offset_x) * (x - offset_x) + (y - offset_y) * (y - offset_y);
                        if distance_squared < well.cutoff_squared {
                            energy += well.value_at_distance_squared(distance_squared);
                            interacting = true;
                        }
                    }
                }
                if !interacting {
                    continue;
                }
                let boltzmann = (-BETA * energy).exp();
                integrand_a += (boltzmann - 1.0) * cell_area;
                integrand_b += energy * boltzmann * cell_area;
            }
        }
        let weight = if index == 0 || index == intervals {
            1.0
        } else if index % 2 == 1 {
            4.0
        } else {
            2.0
        };
        integral_a += weight * integrand_a;
        integral_b += weight * integrand_b;
        integral_numerator += weight * (2.0 * theta).cos() * integrand_a;
    }
    integral_a *= step_theta / 3.0;
    integral_b *= step_theta / 3.0;
    integral_numerator *= step_theta / 3.0;
    let partition = 2.0 * std::f64::consts::PI * ROTOR_LENGTH * ROTOR_LENGTH + 4.0 * integral_a;
    (
        4.0 * integral_numerator / partition,
        4.0 * integral_b / partition,
    )
}

fn assert_rotor_equilibrium(
    couplings: &[f64],
    seed_base: u64,
    seeds: usize,
    thermalization: u64,
    measurement: u64,
) {
    for (coupling_index, epsilon) in couplings.iter().enumerate() {
        let (exact_order, exact_energy) = rotor_reference(*epsilon, 0.05);
        let mut order_means = Vec::with_capacity(seeds);
        let mut energy_means = Vec::with_capacity(seeds);
        for seed_index in 0..seeds {
            let (system, potential, topology) = rotor_system(*epsilon);
            let plan = ChainPlan {
                seed: seed_base + 100_000 * coupling_index as u64 + seed_index as u64,
                thermalization,
                measurement,
                max_displacement: 0.25,
                max_angle: 0.3,
            };
            let mut order_sum = 0.0;
            let mut energy_sum = 0.0;
            let mut sweeps = 0u64;
            run_chain(&plan, system, topology, &potential, |chain| {
                let cell = *chain.configuration().cell();
                let positions = chain.configuration().positions();
                let bond_zero = cell.displacement(&positions[0], &positions[1]);
                let bond_one = cell.displacement(&positions[2], &positions[3]);
                let norm_zero = (bond_zero[0] * bond_zero[0] + bond_zero[1] * bond_zero[1]).sqrt();
                let norm_one = (bond_one[0] * bond_one[0] + bond_one[1] * bond_one[1]).sqrt();
                let dot = (bond_zero[0] * bond_one[0] + bond_zero[1] * bond_one[1])
                    / (norm_zero * norm_one);
                order_sum += 2.0 * dot * dot - 1.0;
                energy_sum += chain.energy;
                sweeps += 1;
            });
            order_means.push(order_sum / sweeps as f64);
            energy_means.push(energy_sum / sweeps as f64);
        }
        assert_seeds_match_exact(
            &order_means,
            exact_order,
            ROTOR_ORDER_FLOOR,
            &format!("rotor ⟨cos2Δθ⟩ (ε={epsilon})"),
        );
        assert_seeds_match_exact(
            &energy_means,
            exact_energy,
            ROTOR_ENERGY_FLOOR,
            &format!("rotor ⟨U⟩ (ε={epsilon})"),
        );
    }
}

/// Orientational locking across three couplings (ε = 1, 2, 3). Δθ is
/// invisible to translation moves, so this curve validates the equilibrium
/// of the rotation move itself. The full seven-point sweep — including
/// the linear-response point ε = 0.5 and the saturated regime up to ε = 6 —
/// runs in `rotor_pair_locking_full_coupling_sweep_long`.
#[test]
fn rotor_pair_orientational_locking_matches_3d_quadrature() {
    assert_rotor_equilibrium(
        &[1.0, 2.0, 3.0],
        ROTOR_SEED_BASE,
        zscore_seed_count(8),
        ROTOR_THERMALIZATION,
        ROTOR_MEASUREMENT,
    );
}

/// Full coupling sweep including the linear-response point (ε = 0.5) and
/// the saturated regime (ε = 6), with longer chains and more seeds.
#[test]
#[ignore = "long: 7-coupling rotor locking sweep, 16 seeds (~2 min; nightly runs --ignored)"]
fn rotor_pair_locking_full_coupling_sweep_long() {
    assert_rotor_equilibrium(
        &[0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ROTOR_SEED_BASE + 7,
        zscore_seed_count(16),
        24_000,
        60_000,
    );
}

// ── Reference self-validation ─────────────────────────────────────────────

/// The quadrature references must be converged: Gaussian integrands under
/// midpoint/Simpson rules reach machine-level accuracy well before the
/// grid sizes used above (measured: step 0.05 → 0.025 changes ⟨cos2Δθ⟩
/// by < 1e-6), so MC deviations cannot be blamed on the reference.
#[test]
fn equilibrium_references_converge_under_refinement() {
    let (pair_energy_coarse, pair_bound_coarse) = pair_reference(2.0, 1024);
    let (pair_energy_fine, pair_bound_fine) = pair_reference(2.0, 4096);
    assert!((pair_energy_coarse - pair_energy_fine).abs() < 1e-12);
    assert!((pair_bound_coarse - pair_bound_fine).abs() < 1e-12);

    let (probe_coarse_order, probe_coarse_energy) = probe_reference(3.0, 0.05);
    let (probe_fine_order, probe_fine_energy) = probe_reference(3.0, 0.0125);
    assert!((probe_coarse_order - probe_fine_order).abs() < 1e-6);
    assert!((probe_coarse_energy - probe_fine_energy).abs() < 1e-6);

    let (rotor_coarse_order, rotor_coarse_energy) = rotor_reference(3.0, 0.05);
    let (rotor_fine_order, rotor_fine_energy) = rotor_reference(3.0, 0.025);
    assert!(
        (rotor_coarse_order - rotor_fine_order).abs() < 1e-6,
        "rotor order reference moved by {:.3e} under refinement",
        (rotor_coarse_order - rotor_fine_order).abs()
    );
    assert!(
        (rotor_coarse_energy - rotor_fine_energy).abs() < 1e-6,
        "rotor energy reference moved by {:.3e} under refinement",
        (rotor_coarse_energy - rotor_fine_energy).abs()
    );
}
