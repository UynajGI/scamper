//! P2 production validation: Boolean/continuum percolation.
//!
//! Gates (roadmap §14 P2), each against an independently coded reference or
//! a paper anchor used strictly as a statistical bracket:
//!
//! 1. **Candidate-pair vs all-pairs overlap parity**: the composed overlap
//!    enumeration (P1 at `2 r_max`, per-pair contact filter) matches an
//!    exhaustive `O(N^2)` overlap scan exactly on random 2D and 3D periodic
//!    polydisperse configurations (pair SET equality; the search is
//!    independent while the minimum-image metric formula is shared — the
//!    gate targets the composed search and filter, not floating-point
//!    noise).
//! 2. **Translation invariance**: a generic shift (then wrap) preserves the
//!    overlap pair set and the label-free component summary; an exact
//!    lattice-vector shift of a dyadic-quantized configuration (bit-stable
//!    under wrap) reproduces the identical pair sequence.
//! 3. **Particle-permutation invariance**: relabeling particles (positions
//!    and radii permuted together) relabels pair identities but leaves
//!    connectivity identical.
//! 4. **Periodic-representative invariance**: shifting individual particles
//!    by integer box vectors (choosing different periodic representatives)
//!    changes nothing.
//! 5. **Zero-radius limit**: as `r -> 0+` (radii are validated positive;
//!    the limit is taken at the smallest radius whose `(L/2r)^D` cell grid
//!    allocates) no pair survives except exactly coincident particles —
//!    exact, because fixture distances are `0` or `>= 1` against a contact
//!    cutoff of `1e-2`.
//! 6. **Fully-overlapping limit**: at the largest radius the minimum-image
//!    ceiling permits — contact cutoff exactly half the shortest side — a
//!    1D ring is the complete graph `K_N` (all `N(N-1)/2` pairs, one
//!    component), and a dense 2D grid collapses into a single component.
//!    (A truly complete 2D overlap graph is unreachable under minimum
//!    image: the box diagonal exceeds half the shortest side.)
//! 7. **Poisson particle-count statistics**: the documented chunked-Knuth
//!    count draw has mean and variance consistent with `lambda` through
//!    proper CI gates, both directly and through the adapter's
//!    `ActiveVertexCount` observable.
//! 8. **Monodisperse 2D disk threshold bracket**: the percolation
//!    probability (a component whose minimum-image-unwrapped x-extent
//!    exceeds 80% of the side) crosses 1/2 between `eta_c/1.15` and
//!    `1.15 eta_c` with `eta_c = 1.12808737(6)` for disks (Mertens &
//!    Moore, PRE 86, 061109 (2012), archived in `docs/ref/`) — a
//!    statistical bracket with real error bars, never an exact claim.
//! 9. **Equilibrium-measurement path**: a frozen caller-supplied
//!    configuration analyzed identically across repeated sweeps (no
//!    resampling) yields zero-variance summaries, through both the adapter
//!    and the zero-copy configuration entry point.
//! 10. **Determinism**: identical seeds reproduce bitwise-identical
//!     estimates.

use std::time::Instant;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze_with_labels, overlap_pairs, overlap_pairs_from_configuration, BooleanError,
    BooleanPercolationMC, ComponentSummary, ComponentWorkspace, ContinuumError, EdgeActivity,
    EdgeId, GraphView, ObservablePlan, OrthorhombicCell, ParticleConfiguration, PointProcess,
    Radii, StaticObservable, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Critical 2D disk filling factor of Mertens & Moore, PRE 86, 061109
/// (2012), Table I ("our work"): `eta_c = 1.12808737(6)`. Used only as a
/// statistical bracket — never labeled exact.
const DISK_ETA_C: f64 = 1.128_087_37;

fn flat<const D: usize>(positions: &[[f64; D]]) -> Vec<f64> {
    positions.iter().flat_map(|p| p.iter().copied()).collect()
}

fn decode<const D: usize>(flat: &[f64]) -> Vec<[f64; D]> {
    (0..flat.len() / D)
        .map(|particle| {
            let mut position = [0.0; D];
            position.copy_from_slice(&flat[particle * D..(particle + 1) * D]);
            position
        })
        .collect()
}

fn random_positions<const D: usize>(n: usize, lengths: [f64; D], seed: u64) -> Vec<[f64; D]> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut position = [0.0; D];
            for axis in 0..D {
                position[axis] = rng.random_range(0.0..lengths[axis]);
            }
            position
        })
        .collect()
}

/// Positions quantized to a dyadic grid (bit-stable under lattice shifts
/// followed by wrapping).
fn quantized_positions<const D: usize>(n: usize, lengths: [f64; D], seed: u64) -> Vec<[f64; D]> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut position = [0.0; D];
            for axis in 0..D {
                position[axis] = (rng.random_range(0u64..512) as f64) * lengths[axis] / 512.0;
            }
            position
        })
        .collect()
}

fn species_radii(n: usize, species: &[f64]) -> Radii {
    Radii::per_particle((0..n).map(|index| species[index % species.len()]))
        .expect("positive species radii")
}

fn sorted_pairs(pairs: &[[usize; 2]]) -> Vec<[usize; 2]> {
    let mut pairs = pairs.to_vec();
    pairs.sort_unstable();
    pairs
}

/// Exhaustive all-pairs overlap reference with an independently coded
/// minimum image (the metric formula is shared with the production cell;
/// the search and the per-pair filter composition are what this gate
/// checks).
fn bruteforce_overlap<const D: usize>(
    positions: &[[f64; D]],
    lengths: &[f64; D],
    radii: &[f64],
) -> Vec<[usize; 2]> {
    let mut pairs = Vec::new();
    for first in 0..positions.len() {
        for second in (first + 1)..positions.len() {
            let mut distance_squared = 0.0;
            for axis in 0..D {
                let raw = positions[second][axis] - positions[first][axis];
                let image = raw - lengths[axis] * (raw * lengths[axis].recip()).round();
                distance_squared += image * image;
            }
            let contact = radii[first] + radii[second];
            if distance_squared <= contact * contact {
                pairs.push([first, second]);
            }
        }
    }
    pairs
}

fn radii_snapshot(radii: &Radii, particles: usize) -> Vec<f64> {
    (0..particles).map(|index| radii.radius(index)).collect()
}

/// Guard: no fixture pair may sit on its contact distance (the filter's
/// boundary is inclusive; an invariance transform must not flip membership
/// through floating-point noise).
fn guard_off_contact<const D: usize>(positions: &[[f64; D]], lengths: &[f64; D], radii: &[f64]) {
    for first in 0..positions.len() {
        for second in (first + 1)..positions.len() {
            let mut distance_squared = 0.0;
            for axis in 0..D {
                let raw = positions[second][axis] - positions[first][axis];
                let image = raw - lengths[axis] * (raw * lengths[axis].recip()).round();
                distance_squared += image * image;
            }
            let contact = radii[first] + radii[second];
            assert!(
                (distance_squared.sqrt() - contact).abs() > 1.0e-9,
                "fixture pair ({first},{second}) sits on its contact distance"
            );
        }
    }
}

/// All-active activity masks mirroring the adapter's internal convention.
struct AllGrains(usize);

impl VertexActivity for AllGrains {
    fn vertex_count(&self) -> usize {
        self.0
    }

    fn vertex_active(&self, _vertex: VertexId) -> bool {
        true
    }
}

struct AllOverlaps(usize);

impl EdgeActivity for AllOverlaps {
    fn edge_count(&self) -> usize {
        self.0
    }

    fn edge_active(&self, _edge: EdgeId) -> bool {
        true
    }
}

/// The label-free summary fields (invariances compare these; only the
/// `largest_component` identity may legitimately change under a
/// relabeling).
fn invariant_fields(value: &ComponentSummary) -> (usize, usize, usize, usize, u128) {
    (
        value.active_vertex_count,
        value.active_edge_count,
        value.component_count,
        value.largest_component_size,
        value.raw_second_moment,
    )
}

fn overlap_summary<const D: usize>(
    positions: &[[f64; D]],
    lengths: [f64; D],
    radii: &Radii,
) -> ComponentSummary {
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let graph = overlap_pairs(&flat(positions), &cell, radii).unwrap();
    let mut workspace = ComponentWorkspace::new();
    let result = analyze_with_labels(
        &graph,
        &AllGrains(graph.vertex_count()),
        &AllOverlaps(graph.pair_count()),
        &[],
        &mut workspace,
    )
    .unwrap();
    result.summary()
}

#[test]
fn overlap_pairs_match_bruteforce_in_2d_and_3d_polydisperse() {
    for seed in [0x21u64, 0x43, 0x87] {
        // 2D: species radii 0.4/0.65/0.9 -> contact cutoff 1.8 <= 2.5.
        let lengths = [6.0, 5.0];
        let positions = random_positions(150, lengths, seed);
        let radii = species_radii(150, &[0.4, 0.65, 0.9]);
        let cell = OrthorhombicCell::new(lengths).unwrap();
        let graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();
        let expected = bruteforce_overlap(&positions, &lengths, &radii_snapshot(&radii, 150));
        assert!(
            !expected.is_empty(),
            "2D fixture should overlap (seed {seed:#x})"
        );
        assert_eq!(sorted_pairs(graph.pairs()), expected, "2D seed {seed:#x}");

        // 3D: species radii 0.3/0.5/0.55 -> contact cutoff 1.1 <= 2.0.
        let lengths = [5.0, 4.5, 4.0];
        let positions = random_positions(200, lengths, seed);
        let radii = species_radii(200, &[0.3, 0.5, 0.55]);
        let cell = OrthorhombicCell::new(lengths).unwrap();
        let graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();
        let expected = bruteforce_overlap(&positions, &lengths, &radii_snapshot(&radii, 200));
        assert!(
            !expected.is_empty(),
            "3D fixture should overlap (seed {seed:#x})"
        );
        assert_eq!(sorted_pairs(graph.pairs()), expected, "3D seed {seed:#x}");
    }
}

#[test]
fn translation_invariance_generic_and_lattice_shifts() {
    let lengths = [8.0, 8.0];

    // Generic (non-lattice) shift on a random fixture: relative geometry
    // unchanged, so the pair set and every label-free summary field are
    // invariant even though particles cross cell-list boundaries.
    let positions = random_positions(180, lengths, 0x5a);
    let radii = species_radii(180, &[0.35, 0.6]);
    guard_off_contact(&positions, &lengths, &radii_snapshot(&radii, 180));
    let reference = overlap_summary(&positions, lengths, &radii);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let reference_graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();

    let shift = [0.37, 2.71];
    let shifted: Vec<[f64; 2]> = positions
        .iter()
        .map(|p| {
            [
                (p[0] + shift[0]).rem_euclid(lengths[0]),
                (p[1] + shift[1]).rem_euclid(lengths[1]),
            ]
        })
        .collect();
    let shifted_graph = overlap_pairs(&flat(&shifted), &cell, &radii).unwrap();
    assert_eq!(
        sorted_pairs(shifted_graph.pairs()),
        sorted_pairs(reference_graph.pairs())
    );
    assert_eq!(
        invariant_fields(&overlap_summary(&shifted, lengths, &radii)),
        invariant_fields(&reference)
    );

    // Exact lattice shift on a dyadic-quantized fixture: wrapped
    // coordinates are bit-identical, so the canonical pair sequence itself
    // reproduces.
    let quantized = quantized_positions(200, lengths, 0x7f);
    let radii = species_radii(200, &[0.35, 0.6]);
    guard_off_contact(&quantized, &lengths, &radii_snapshot(&radii, 200));
    let reference_graph = overlap_pairs(&flat(&quantized), &cell, &radii).unwrap();
    let lattice_shifted: Vec<[f64; 2]> = quantized
        .iter()
        .map(|p| [(p[0] + 8.0).rem_euclid(8.0), (p[1] - 16.0).rem_euclid(8.0)])
        .collect();
    assert_eq!(
        lattice_shifted, quantized,
        "dyadic quantization wraps exactly"
    );
    let lattice_graph = overlap_pairs(&flat(&lattice_shifted), &cell, &radii).unwrap();
    assert_eq!(lattice_graph.pairs(), reference_graph.pairs());
}

#[test]
fn particle_permutation_invariance() {
    let lengths = [7.0, 6.0];
    let positions = random_positions(160, lengths, 0x6c);
    let radii = species_radii(160, &[0.4, 0.55, 0.8]);
    guard_off_contact(&positions, &lengths, &radii_snapshot(&radii, 160));

    let reference = overlap_summary(&positions, lengths, &radii);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let reference_graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();

    // Fisher-Yates permutation of the particle labels; positions and radii
    // travel together. `permutation[new] = old`, so the old->new map is its
    // inverse.
    let mut permutation: Vec<usize> = (0..160).collect();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x9e);
    for index in (1..permutation.len()).rev() {
        let swap_with = rng.random_range(0..=index);
        permutation.swap(index, swap_with);
    }
    let mut inverse = vec![0usize; permutation.len()];
    for (new, &old) in permutation.iter().enumerate() {
        inverse[old] = new;
    }
    let permuted_positions: Vec<[f64; 2]> =
        permutation.iter().map(|&index| positions[index]).collect();
    let permuted_radii = Radii::per_particle(permutation.iter().map(|&index| radii.radius(index)))
        .expect("permutation of valid radii");
    let permuted_graph = overlap_pairs(&flat(&permuted_positions), &cell, &permuted_radii).unwrap();

    // Pair identities relabel; the connectivity does not.
    let mapped: Vec<[usize; 2]> = reference_graph
        .pairs()
        .iter()
        .map(|&[left, right]| {
            let mut pair = [inverse[left], inverse[right]];
            pair.sort_unstable();
            pair
        })
        .collect();
    assert_eq!(sorted_pairs(permuted_graph.pairs()), sorted_pairs(&mapped));
    assert_eq!(
        invariant_fields(&overlap_summary(
            &permuted_positions,
            lengths,
            &permuted_radii
        )),
        invariant_fields(&reference)
    );
}

#[test]
fn periodic_representative_invariance() {
    let lengths = [6.0, 5.0];
    let positions = random_positions(140, lengths, 0x7b);
    let radii = species_radii(140, &[0.45, 0.7]);
    guard_off_contact(&positions, &lengths, &radii_snapshot(&radii, 140));

    let reference = overlap_summary(&positions, lengths, &radii);
    let cell = OrthorhombicCell::new(lengths).unwrap();
    let reference_graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();

    // Shift individual particles by (distinct) integer multiples of the
    // per-axis BOX LENGTHS: the same periodic configuration under different
    // representatives.
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x1d);
    let rewrapped: Vec<[f64; 2]> = positions
        .iter()
        .map(|p| {
            [
                p[0] + lengths[0] * (rng.random_range(-3i64..=3) as f64),
                p[1] + lengths[1] * (rng.random_range(-3i64..=3) as f64),
            ]
        })
        .collect();
    let rewrapped_graph = overlap_pairs(&flat(&rewrapped), &cell, &radii).unwrap();
    assert_eq!(
        sorted_pairs(rewrapped_graph.pairs()),
        sorted_pairs(reference_graph.pairs())
    );
    assert_eq!(
        invariant_fields(&overlap_summary(&rewrapped, lengths, &radii)),
        invariant_fields(&reference)
    );
}

#[test]
fn zero_radius_limit_keeps_only_coincident_pairs_exactly() {
    // Radii are validated strictly positive, so the zero-radius limit is
    // taken at the smallest radius whose (L/2r)^D cell grid allocates
    // comfortably (r = 5e-3 in a 6x6 box -> 600x600 cells). The limit
    // statement is exact for this fixture: distinct particles sit at
    // distances >= 1 against a contact cutoff of 1e-2, so the inclusive
    // comparison `d^2 <= (2r)^2` is decided exactly in floating point.
    let cell = OrthorhombicCell::new([6.0, 6.0]).unwrap();
    let mut positions = Vec::new();
    for row in 0..5 {
        for column in 0..5 {
            positions.push([0.5 + column as f64, 0.5 + row as f64]);
        }
    }
    let radii = Radii::monodisperse(5.0e-3).unwrap();
    let graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();
    assert_eq!(graph.pair_count(), 0, "distinct particles cannot overlap");

    // Two particles at the identical position: exactly one pair (d = 0).
    positions[7] = positions[6];
    let graph = overlap_pairs(&flat(&positions), &cell, &radii).unwrap();
    assert_eq!(graph.pairs(), &[[6usize, 7usize]]);

    // The same limit through per-particle radii.
    let tiny = Radii::per_particle([5.0e-3; 25]).unwrap();
    let graph = overlap_pairs(&flat(&positions), &cell, &tiny).unwrap();
    assert_eq!(graph.pairs(), &[[6usize, 7usize]]);
}

#[test]
fn fully_overlapping_limits() {
    // 1D ring at the largest permitted cutoff (exactly L/2): every pair's
    // minimum-image distance is <= L/2, so the graph is the complete K_6 —
    // truly fully overlapping, with the exact pair count 15.
    let ring = OrthorhombicCell::new([8.0]).unwrap();
    let positions: Vec<[f64; 1]> = [0.5, 1.5, 3.25, 4.0, 6.75, 7.5]
        .iter()
        .map(|&x| [x])
        .collect();
    let radii = Radii::monodisperse(2.0).unwrap();
    assert_eq!(radii.contact_cutoff(), 4.0);
    let graph = overlap_pairs(&flat(&positions), &ring, &radii).unwrap();
    assert_eq!(graph.pair_count(), 15, "complete graph K_6");
    assert_eq!(
        sorted_pairs(graph.pairs()),
        sorted_pairs(&bruteforce_overlap(
            &positions,
            &[8.0],
            &radii_snapshot(&radii, 6)
        ))
    );
    let summary = overlap_summary(&positions, [8.0], &radii);
    assert_eq!(summary.component_count, 1);
    assert_eq!(summary.largest_component_size, 6);
    assert_eq!(summary.raw_second_moment, 36);

    // 2D: the minimum-image ceiling caps the contact cutoff at L/2, below
    // the box diagonal, so a complete graph is unreachable — but
    // connectivity saturates: an 8x8 grid at r = L/4 collapses into one
    // component.
    let grid: Vec<[f64; 2]> = (0..8)
        .flat_map(|row| (0..8).map(move |column| [0.5 + column as f64, 0.5 + row as f64]))
        .collect();
    let radii = Radii::monodisperse(2.0).unwrap();
    let summary = overlap_summary(&grid, [8.0, 8.0], &radii);
    assert_eq!(summary.component_count, 1, "saturated connectivity");
    assert_eq!(summary.largest_component_size, 64);
    assert_eq!(summary.raw_second_moment, 64 * 64);
}

#[test]
fn poisson_count_draw_matches_mean_and_variance() {
    // Direct count-draw statistics through the documented chunked-Knuth
    // contract (the second case exercises the chunked path, mean > 512).
    // Gates: |mean - lambda| <= z sqrt(lambda/n) and |s^2 - lambda| <=
    // z sqrt((lambda + 2 lambda^2)/n) — the asymptotic normal CI of the
    // unbiased sample variance under Poisson (mu4 - sigma^4 = lambda +
    // 2 lambda^2).
    for (lambda, draws) in [(37.7_f64, 6000_usize), (1500.0, 800)] {
        let process = PointProcess::poisson(lambda).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x50_0f);
        let counts: Vec<f64> = (0..draws)
            .map(|_| process.draw_count(&mut rng) as f64)
            .collect();
        let n = draws as f64;
        let mean = counts.iter().sum::<f64>() / n;
        let variance = counts.iter().map(|c| (c - mean) * (c - mean)).sum::<f64>() / (n - 1.0);
        let z = 4.0;
        let mean_tolerance = z * (lambda / n).sqrt();
        let variance_tolerance = z * ((lambda + 2.0 * lambda * lambda) / n).sqrt();
        eprintln!(
            "BOOLEAN_POISSON lambda={lambda} draws={draws} mean={mean:.4} \
             (tol {mean_tolerance:.4}) variance={variance:.4} (tol {variance_tolerance:.4})"
        );
        assert!(
            (mean - lambda).abs() <= mean_tolerance,
            "lambda {lambda}: mean {mean} vs tolerance {mean_tolerance}"
        );
        assert!(
            (variance - lambda).abs() <= variance_tolerance,
            "lambda {lambda}: variance {variance} vs tolerance {variance_tolerance}"
        );
    }
}

#[test]
fn adapter_poisson_count_statistics_through_observable() {
    // The same statistics through the full adapter path: ActiveVertexCount
    // IS the per-realization particle count.
    let lambda = 40.0;
    let sweeps = 3000;
    let mut adapter = BooleanPercolationMC::<2>::sampling(
        OrthorhombicCell::new([16.0, 16.0]).unwrap(),
        Radii::monodisperse(0.25).unwrap(),
        PointProcess::poisson(lambda).unwrap(),
        ObservablePlan::new([StaticObservable::ActiveVertexCount], Vec::new()),
    )
    .unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0xC_0F), 0, 1);
    for _ in 0..sweeps {
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
    }
    let values = ctx.finalize_measurements();
    let estimate = &values["ActiveVertexCount"];
    let n = sweeps as f64;
    // Sample variance of the per-realization counts from the recorded
    // stderr (binsize 1): s^2 = stderr^2 * n.
    let mean = estimate.mean;
    let variance = estimate.stderr * estimate.stderr * n;
    let z = 4.0;
    eprintln!(
        "BOOLEAN_POISSON_ADAPTER lambda={lambda} sweeps={sweeps} mean={mean:.4} \
         variance={variance:.4} n_bins={}",
        estimate.n_bins
    );
    assert_eq!(estimate.n_bins, sweeps);
    assert!(
        (mean - lambda).abs() <= z * (lambda / n).sqrt(),
        "mean {mean}"
    );
    assert!(
        (variance - lambda).abs() <= z * ((lambda + 2.0 * lambda * lambda) / n).sqrt(),
        "variance {variance}"
    );
}

#[test]
fn frozen_configuration_is_measured_identically_every_sweep() {
    // The equilibrium-chain contract: an existing simulation's accepted
    // coordinates are a measurement target, never resampled. Repeated
    // sweeps must analyze the identical configuration and record identical
    // values (zero variance), through both the adapter and the zero-copy
    // configuration entry point.
    let lengths = [10.0, 10.0];
    let positions = random_positions(200, lengths, 0xF_F2);
    let radii = species_radii(200, &[0.3, 0.6]);
    let coordinates = flat(&positions);
    let cell = OrthorhombicCell::new(lengths).unwrap();

    let mut adapter = BooleanPercolationMC::<2>::measuring(
        cell,
        coordinates.clone(),
        radii.clone(),
        ObservablePlan::all(Vec::new()),
    )
    .unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(3), 0, 1);
    for _ in 0..12 {
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        // No resampling: the frozen buffer never changes.
        assert_eq!(adapter.coordinates(), coordinates.as_slice());
    }
    let values = ctx.finalize_measurements();
    for name in [
        "ActiveVertexCount",
        "ActiveEdgeCount",
        "ComponentCount",
        "LargestSize",
        "GiantFraction",
        "RawSecondMoment",
    ] {
        let estimate = &values[name];
        assert_eq!(estimate.n_bins, 12, "{name}");
        // Identical samples: the bin variance is exactly zero for integer
        // observables and at most one ulp of bin-mean rounding for the
        // rational GiantFraction (sum of 12 identical f64s / 12).
        assert!(
            estimate.stderr.abs() < 1.0e-12,
            "{name} must be identical every sweep (stderr {})",
            estimate.stderr
        );
    }
    // A nontrivial frozen configuration: overlapping, multi-component.
    assert!(values["ActiveEdgeCount"].mean > 0.0);
    assert!(values["ComponentCount"].mean >= 1.0);

    // The zero-copy configuration path enumerates the identical overlap set.
    let configuration = ParticleConfiguration::new(
        positions,
        vec![0; coordinates.len() / 2],
        OrthorhombicCell::new(lengths).unwrap(),
    )
    .unwrap();
    let via_configuration = overlap_pairs_from_configuration(&configuration, &radii).unwrap();
    let via_flat = overlap_pairs(&coordinates, &cell, &radii).unwrap();
    assert_eq!(via_configuration.pairs(), via_flat.pairs());
}

#[test]
fn sampling_is_bitwise_deterministic_per_seed() {
    let run = |seed: u64| {
        let mut adapter = BooleanPercolationMC::<3>::sampling(
            OrthorhombicCell::new([6.0, 5.5, 5.0]).unwrap(),
            species_radii(300, &[0.3, 0.5]),
            PointProcess::binomial(300).unwrap(),
            ObservablePlan::all(Vec::new()),
        )
        .unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
        for _ in 0..60 {
            adapter.sweep(&mut ctx);
            adapter.measure(&mut ctx);
        }
        ctx.finalize_measurements()
    };
    let first = run(0x0e0f);
    let second = run(0x0e0f);
    for name in [
        "ActiveVertexCount",
        "ActiveEdgeCount",
        "ComponentCount",
        "LargestSize",
        "GiantFraction",
        "RawSecondMoment",
        "FiniteClusterSusceptibility",
    ] {
        assert_eq!(first[name].mean, second[name].mean, "{name}");
        assert_eq!(first[name].n_bins, second[name].n_bins, "{name}");
    }
    let other = run(0x0e10);
    assert_ne!(first["ActiveEdgeCount"].mean, other["ActiveEdgeCount"].mean);
}

#[test]
fn typed_rejections_reach_integration_callers() {
    let cell = OrthorhombicCell::new([4.0, 4.0]).unwrap();
    let radii = Radii::monodisperse(0.5).unwrap();
    // Empty coordinates: the P1 rejection wraps unchanged.
    assert_eq!(
        overlap_pairs(&[], &cell, &radii),
        Err(BooleanError::Continuum(ContinuumError::EmptyCoordinates {
            dimension: 2
        }))
    );
    // Radius length mismatch against the coordinate count.
    let radii = Radii::per_particle([0.5; 3]).unwrap();
    assert_eq!(
        overlap_pairs(&[0.0; 4], &cell, &radii),
        Err(BooleanError::RadiiLengthMismatch {
            expected: 2,
            actual: 3
        })
    );
    // Contact cutoff beyond the minimum-image maximum.
    let radii = Radii::monodisperse(1.3).unwrap();
    assert!(matches!(
        overlap_pairs(&[0.5, 0.5, 1.0, 1.0], &cell, &radii),
        Err(BooleanError::Continuum(
            ContinuumError::CutoffTooLarge { .. }
        ))
    ));
}

/// Wilson score interval for a binomial proportion (the gate's proper CI).
fn wilson_bounds(successes: usize, trials: usize, z: f64) -> (f64, f64) {
    let n = trials as f64;
    let p = successes as f64 / n;
    let z_squared = z * z;
    let denominator = 1.0 + z_squared / n;
    let center = (p + z_squared / (2.0 * n)) / denominator;
    let half = z * (p * (1.0 - p) / n + z_squared / (4.0 * n * n)).sqrt() / denominator;
    (center - half, center + half)
}

/// Spanning criterion of the threshold gate: some component's UNWRAPPED
/// x-extent exceeds 80% of the side — a left-right crossing of the central
/// strip, measured correctly on the torus.
///
/// Raw x-extents are meaningless on a periodic axis (a two-particle cluster
/// straddling the seam has raw extent ~L), so each component is unfolded by
/// a BFS over its overlap edges, accumulating minimum-image x-displacements
/// into unwrapped coordinates `x_v + n_v L`. A component crosses iff its
/// unwrapped extent `max - min > 0.8 L`; any x-wrapping component
/// necessarily has unwrapped extent >= L and qualifies. Connectivity comes
/// from the production overlap enumeration; only the unfolding walk and the
/// geometric judgment are test-side.
fn spans_x_unwrapped(pairs: &[[usize; 2]], positions: &[[f64; 2]], side: f64) -> bool {
    let n = positions.len();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &[left, right] in pairs {
        adjacency[left].push(right);
        adjacency[right].push(left);
    }
    let mut unwrapped = vec![f64::NAN; n];
    for seed in 0..n {
        if !unwrapped[seed].is_nan() {
            continue;
        }
        unwrapped[seed] = positions[seed][0];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(seed);
        let mut low = positions[seed][0];
        let mut high = positions[seed][0];
        while let Some(vertex) = queue.pop_front() {
            for &neighbor in &adjacency[vertex] {
                if !unwrapped[neighbor].is_nan() {
                    continue;
                }
                let raw = positions[neighbor][0] - positions[vertex][0];
                let image = raw - side * (raw / side).round();
                let proposed = unwrapped[vertex] + image;
                // Snap to the representative nearest the proposal.
                let periods = ((proposed - positions[neighbor][0]) / side).round();
                unwrapped[neighbor] = positions[neighbor][0] + periods * side;
                low = low.min(unwrapped[neighbor]);
                high = high.max(unwrapped[neighbor]);
                queue.push_back(neighbor);
            }
        }
        if high - low > 0.8 * side {
            return true;
        }
    }
    false
}

#[test]
fn monodisperse_disk_threshold_statistical_bracket() {
    // Roadmap gate: the measured percolation probability crosses 1/2
    // between eta_c/1.15 and 1.15 eta_c on N = 5e4-1e5 particles. eta is
    // the 2D disk filling factor N pi r^2 / A; the anchor is Mertens &
    // Moore's Table I estimate 1.12808737(6) — a statistical bracket with
    // real error bars, never an exact threshold.
    const PARTICLES: usize = 50_000;
    const SIDE: f64 = 400.0;
    const REALIZATIONS: usize = 12;
    const Z: f64 = 2.0; // ~95% Wilson bounds

    let run_side = |factor: f64, seed: u64| {
        let eta = DISK_ETA_C * factor;
        let radius = (eta * SIDE * SIDE / (std::f64::consts::PI * PARTICLES as f64)).sqrt();
        let mut adapter = BooleanPercolationMC::<2>::sampling(
            OrthorhombicCell::new([SIDE, SIDE]).unwrap(),
            Radii::monodisperse(radius).unwrap(),
            PointProcess::binomial(PARTICLES).unwrap(),
            ObservablePlan::new([StaticObservable::GiantFraction], Vec::new()),
        )
        .unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
        let mut crossings = 0;
        let mut largest_fraction_sum = 0.0_f64;
        let start = Instant::now();
        for _ in 0..REALIZATIONS {
            // Draw the realization through the production sampling path;
            // the spanning judgment runs on the production overlap pairs.
            adapter.sweep(&mut ctx);
            let graph =
                overlap_pairs(adapter.coordinates(), adapter.cell(), adapter.radii()).unwrap();
            let positions = decode::<2>(adapter.coordinates());
            if spans_x_unwrapped(graph.pairs(), &positions, SIDE) {
                crossings += 1;
            }
            // Giant fraction from a plain component pass (record-side
            // observable check without a second labeled analysis).
            let summary = cmc_rs::analyze(
                &graph,
                &AllGrains(graph.vertex_count()),
                &AllOverlaps(graph.pair_count()),
                &[],
                &mut ComponentWorkspace::new(),
            )
            .unwrap()
            .summary();
            largest_fraction_sum +=
                summary.largest_component_size as f64 / graph.vertex_count() as f64;
        }
        eprintln!(
            "BOOLEAN_THRESHOLD eta={eta:.5} (factor {factor:.3} of eta_c={DISK_ETA_C}) \
             r={radius:.5} N={PARTICLES} realizations={REALIZATIONS} \
             crossings={crossings} P={:.4} giant_fraction_mean={:.4} \
             elapsed={:.1}s",
            crossings as f64 / REALIZATIONS as f64,
            largest_fraction_sum / REALIZATIONS as f64,
            start.elapsed().as_secs_f64(),
        );
        crossings
    };

    let low = run_side(1.0 / 1.15, 0xB_1D);
    let high = run_side(1.15, 0xB_1E);
    let (_, low_upper) = wilson_bounds(low, REALIZATIONS, Z);
    let (high_lower, _) = wilson_bounds(high, REALIZATIONS, Z);
    assert!(
        low_upper < 0.5,
        "P(eta_c/1.15) Wilson upper {low_upper:.4} must sit below 1/2"
    );
    assert!(
        high_lower > 0.5,
        "P(1.15 eta_c) Wilson lower {high_lower:.4} must sit above 1/2"
    );
    eprintln!(
        "BOOLEAN_THRESHOLD_BRACKET Wilson95: low in [0, {low_upper:.3}], \
         high in [{high_lower:.3}, 1] — the crossing of 1/2 lies inside \
         (eta_c/1.15, 1.15 eta_c); eta_c = {DISK_ETA_C} is a statistical \
         bracket anchor (Mertens & Moore 2012), never exact"
    );
}
