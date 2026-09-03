//! Exact-enumeration and reference-implementation validation for percolation.
//!
//! Three independent layers, each covering site, bond and mixed site-bond
//! modes:
//!
//! 1. **Hand-derived closed forms.** The 2x2 open square is small enough to
//!    enumerate exhaustively against moments derived by hand; the open chain
//!    has exact crossing forms `P = p^L` (site), `P = p^(L-1)` (bond) and
//!    `P = p_s^L * p_b^(L-1)` (mixed), and mixed 2x2 satisfies
//!    `P = 2 p_s^2 p_b - p_s^4 p_b^2` (only an active horizontal bond can
//!    cross, and the two rows coincide only when everything is open).
//! 2. **Independent algorithm cross-check.** `cluster_stats` (union find) is
//!    compared configuration-by-configuration against a flood-fill reference
//!    sharing no algorithmic path — exhaustively on small lattices, on seeded
//!    random configurations for larger ones, across every lattice family
//!    (chain, square, cubic, triangular, honeycomb, kagome, random graphs).
//! 3. **Structural identities.** Mixed percolation at `p_site = 1` must
//!    reproduce pure bond statistics exactly, and at `p_bond = 1` pure site
//!    statistics; crossing probability is monotone in each occupation
//!    probability separately; fixed seeds reproduce bitwise.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square, build_triangular,
    cluster_stats, Bond, BondType, ClusterStats, CsrLattice, OccupancyState, PercolationMC,
    PercolationMode,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const L: usize = 2;
const FROM: [usize; 2] = [0, 2];
const TO: [usize; 2] = [1, 3];
const ALL_MODES: [PercolationMode; 3] = [
    PercolationMode::Site,
    PercolationMode::Bond,
    PercolationMode::SiteBond,
];

/// Exhaustive enumeration is capped at 2^14 configurations; anything larger
/// is covered by seeded random sampling.
const MAX_EXHAUSTIVE_BITS: usize = 14;
/// Seeded random configurations tried when enumeration is infeasible.
const SAMPLED_CONFIGS: usize = 4_000;

/// Number of sampled elements (sites, bonds or both) for `mode`.
fn element_count(lattice: &CsrLattice, mode: PercolationMode) -> usize {
    match mode {
        PercolationMode::Site => lattice.n_sites,
        PercolationMode::Bond => lattice.n_edges(),
        PercolationMode::SiteBond => lattice.n_sites + lattice.n_edges(),
    }
}

/// Write mask bits into the occupancy arrays sampled by `mode`. Bit layout:
/// pure modes pack their sampled array into the low bits; the mixed mode
/// packs sites low and bonds above them.
fn fill_config(occupancy: &mut OccupancyState, mode: PercolationMode, mask: usize) {
    if mode.samples_sites() {
        for (bit, open) in occupancy.site_open.iter_mut().enumerate() {
            *open = (mask >> bit) & 1 == 1;
        }
    }
    if mode.samples_bonds() {
        let offset = match mode {
            PercolationMode::SiteBond => occupancy.site_open.len(),
            _ => 0,
        };
        for (bit, open) in occupancy.bond_open.iter_mut().enumerate() {
            *open = (mask >> (offset + bit)) & 1 == 1;
        }
    }
}

/// Occupied-element counts `(sites, bonds)` encoded in `mask`, mirroring the
/// bit layout of [`fill_config`].
fn mask_populations(lattice: &CsrLattice, mode: PercolationMode, mask: usize) -> (usize, usize) {
    match mode {
        PercolationMode::Site => (mask.count_ones() as usize, 0),
        PercolationMode::Bond => (0, mask.count_ones() as usize),
        PercolationMode::SiteBond => {
            let site_mask = (1usize << lattice.n_sites) - 1;
            (
                (mask & site_mask).count_ones() as usize,
                (mask >> lattice.n_sites).count_ones() as usize,
            )
        }
    }
}

/// Configuration weight at `(p_site, p_bond)`: the product over sampled
/// element kinds of `p^k (1-p)^(n-k)`.
fn config_weight(
    lattice: &CsrLattice,
    mode: PercolationMode,
    open_sites: usize,
    open_bonds: usize,
    p_site: f64,
    p_bond: f64,
) -> f64 {
    let (ns, ne) = (lattice.n_sites, lattice.n_edges());
    let site_factor =
        p_site.powi(open_sites as i32) * (1.0 - p_site).powi((ns - open_sites) as i32);
    let bond_factor =
        p_bond.powi(open_bonds as i32) * (1.0 - p_bond).powi((ne - open_bonds) as i32);
    match mode {
        PercolationMode::Site => site_factor,
        PercolationMode::Bond => bond_factor,
        PercolationMode::SiteBond => site_factor * bond_factor,
    }
}

/// Enumerate (or, for large lattices, densely sample) occupancy
/// configurations, returning probability-weighted `<max>`, `<sum(s^2)>`,
/// `<n_clusters>` and P(spanning) at `(p_site, p_bond)`. Pure modes use only
/// the probability of their sampled kind.
fn enumerate_moments(
    lattice: &CsrLattice,
    mode: PercolationMode,
    p_site: f64,
    p_bond: f64,
    from: &[usize],
    to: &[usize],
) -> (f64, f64, f64, f64) {
    let n_bits = element_count(lattice, mode);
    let mut occupancy = OccupancyState::new(lattice, mode);
    let (mut mean_max, mut mean_s2, mut mean_n, mut p_span) = (0.0, 0.0, 0.0, 0.0);
    for mask in 0..(1usize << n_bits) {
        fill_config(&mut occupancy, mode, mask);
        let (open_sites, open_bonds) = mask_populations(lattice, mode, mask);
        let weight = config_weight(lattice, mode, open_sites, open_bonds, p_site, p_bond);
        let stats = cluster_stats(lattice, &occupancy, from, to);
        mean_max += weight * stats.max_size as f64;
        mean_s2 += weight * stats.second_moment as f64;
        mean_n += weight * stats.n_clusters as f64;
        p_span += weight * f64::from(u8::from(stats.spanning));
    }
    (mean_max, mean_s2, mean_n, p_span)
}

/// Independent cluster statistics: flood fill with an explicit stack instead
/// of union find, so the reference shares no algorithmic path with
/// [`cmc_rs::cluster_stats`].
fn reference_stats(
    lattice: &CsrLattice,
    occupancy: &OccupancyState,
    from: &[usize],
    to: &[usize],
) -> ClusterStats {
    let edge_active = |edge_id: usize| match occupancy.mode {
        PercolationMode::Site => {
            let edge = &lattice.edges[edge_id];
            occupancy.site_open[edge.source] && occupancy.site_open[edge.target]
        }
        PercolationMode::Bond => occupancy.bond_open[edge_id],
        PercolationMode::SiteBond => {
            let edge = &lattice.edges[edge_id];
            occupancy.bond_open[edge_id]
                && occupancy.site_open[edge.source]
                && occupancy.site_open[edge.target]
        }
    };
    let tracked = |site: usize| match occupancy.mode {
        PercolationMode::Bond => true,
        _ => occupancy.site_open[site],
    };

    let mut seen = vec![false; lattice.n_sites];
    let mut stats = ClusterStats {
        max_size: 0,
        second_moment: 0,
        n_clusters: 0,
        spanning: false,
    };
    for seed in 0..lattice.n_sites {
        if seen[seed] || !tracked(seed) {
            continue;
        }
        seen[seed] = true;
        let mut stack = vec![seed];
        let mut size = 0usize;
        let mut touches_from = false;
        let mut touches_to = false;
        while let Some(site) = stack.pop() {
            size += 1;
            touches_from |= from.contains(&site);
            touches_to |= to.contains(&site);
            for (i, &edge_id) in lattice.edge_ids(site).iter().enumerate() {
                if !edge_active(edge_id) {
                    continue;
                }
                let neighbor = lattice.neighbors(site)[i];
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        stats.n_clusters += 1;
        stats.max_size = stats.max_size.max(size);
        stats.second_moment += (size as u64) * (size as u64);
        stats.spanning |= touches_from && touches_to;
    }
    stats
}

/// Compare `cluster_stats` against the flood-fill reference for every mode:
/// on every configuration (small lattices) or on seeded random configurations.
fn verify_against_reference(lattice: &CsrLattice, label: &str, from: &[usize], to: &[usize]) {
    for mode in ALL_MODES {
        let n_bits = element_count(lattice, mode);
        let exhaustive = n_bits <= MAX_EXHAUSTIVE_BITS;
        let trials = if exhaustive {
            1usize << n_bits
        } else {
            SAMPLED_CONFIGS
        };
        let mut occupancy = OccupancyState::new(lattice, mode);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);
        for trial in 0..trials {
            if exhaustive {
                fill_config(&mut occupancy, mode, trial);
            } else {
                occupancy.resample(0.5, 0.5, &mut rng);
            }
            let expected = reference_stats(lattice, &occupancy, from, to);
            let actual = cluster_stats(lattice, &occupancy, from, to);
            assert_eq!(
                actual, expected,
                "{label}/{mode:?} mismatch on configuration {trial}"
            );
        }
    }
}

/// Seeded Erdős–Rényi-style graph: the arbitrary-topology production path.
fn random_graph() -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
    let n = 10;
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            if rng.random::<f64>() < 0.35 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            }
        }
    }
    CsrLattice::try_from_edges(n, edges).expect("random graph is valid")
}

#[test]
fn exact_enumeration_matches_hand_derived_moments() {
    let lattice = build_square(L, L, false);
    // Site mode: (30/16, 76/16, 17/16, 7/16) at p = 1/2.
    let (max, s2, n, span) =
        enumerate_moments(&lattice, PercolationMode::Site, 0.5, 0.5, &FROM, &TO);
    assert_close(max, 1.875, 1e-12, "site <max>");
    assert_close(s2, 4.75, 1e-12, "site <s^2>");
    assert_close(n, 17.0 / 16.0, 1e-12, "site <n_clusters>");
    assert_close(span, 0.4375, 1e-12, "site P(span)");

    // Bond mode: (45/16, 164/16, 33/16, 12/16) at p = 1/2.
    let (max, s2, n, span) =
        enumerate_moments(&lattice, PercolationMode::Bond, 0.5, 0.5, &FROM, &TO);
    assert_close(max, 2.8125, 1e-12, "bond <max>");
    assert_close(s2, 10.25, 1e-12, "bond <s^2>");
    assert_close(n, 33.0 / 16.0, 1e-12, "bond <n_clusters>");
    assert_close(span, 0.75, 1e-12, "bond P(span)");
}

#[test]
fn exact_spanning_probability_follows_closed_form() {
    let lattice = build_square(L, L, false);
    // Site mode: spanning configs are the two rows (2 of size 2), all four
    // triples and the full config: P = 2p^2(1-p)^2 + 4p^3(1-p) + p^4.
    for p in [0.2, 0.44, 0.5927, 0.8] {
        let (_, _, _, span) = enumerate_moments(&lattice, PercolationMode::Site, p, p, &FROM, &TO);
        let closed_form =
            2.0 * p * p * (1.0 - p) * (1.0 - p) + 4.0 * p.powi(3) * (1.0 - p) + p.powi(4);
        assert_close(span, closed_form, 1e-12, "site spanning polynomial");
    }
}

/// Mixed 2x2 closed form: a crossing needs an active horizontal bond, so
/// P(span) = 2 p_s^2 p_b - p_s^4 p_b^2 — the two rows are the only routes
/// and coincide exactly when all four sites and both bonds are open.
#[test]
fn site_bond_spanning_matches_hand_derived_closed_form() {
    let lattice = build_square(L, L, false);
    for (p_s, p_b) in [(0.6, 0.7), (0.9, 0.4), (1.0, 0.5), (0.5, 1.0), (0.3, 0.3)] {
        let (_, _, _, span) =
            enumerate_moments(&lattice, PercolationMode::SiteBond, p_s, p_b, &FROM, &TO);
        let closed_form = 2.0 * p_s * p_s * p_b - p_s.powi(4) * p_b * p_b;
        assert_close(span, closed_form, 1e-12, "site-bond spanning polynomial");
    }
}

/// Degenerate mixed ensembles must reproduce the pure modes exactly:
/// `p_site = 1` reduces to pure bond, `p_bond = 1` to pure site.
#[test]
fn mixed_mode_reduces_to_pure_modes_at_degenerate_probabilities() {
    let lattice = build_square(L, L, false);
    for p in [0.2, 0.5, 0.8] {
        let pure_bond = enumerate_moments(&lattice, PercolationMode::Bond, p, p, &FROM, &TO);
        let mixed_as_bond =
            enumerate_moments(&lattice, PercolationMode::SiteBond, 1.0, p, &FROM, &TO);
        assert_close(mixed_as_bond.0, pure_bond.0, 1e-12, "p_site=1 <max>");
        assert_close(mixed_as_bond.1, pure_bond.1, 1e-12, "p_site=1 <s^2>");
        assert_close(mixed_as_bond.2, pure_bond.2, 1e-12, "p_site=1 <n>");
        assert_close(mixed_as_bond.3, pure_bond.3, 1e-12, "p_site=1 P(span)");

        let pure_site = enumerate_moments(&lattice, PercolationMode::Site, p, p, &FROM, &TO);
        let mixed_as_site =
            enumerate_moments(&lattice, PercolationMode::SiteBond, p, 1.0, &FROM, &TO);
        assert_close(mixed_as_site.0, pure_site.0, 1e-12, "p_bond=1 <max>");
        assert_close(mixed_as_site.1, pure_site.1, 1e-12, "p_bond=1 <s^2>");
        assert_close(mixed_as_site.2, pure_site.2, 1e-12, "p_bond=1 <n>");
        assert_close(mixed_as_site.3, pure_site.3, 1e-12, "p_bond=1 P(span)");
    }
}

/// 1D exact solution: on an open chain spanning requires every site (site
/// mode), every bond (bond mode) or both (mixed mode) to be occupied.
#[test]
fn chain_crossing_follows_exact_closed_forms() {
    let length = 8;
    let lattice = build_chain(length, false);
    let ends = (&[0][..], &[length - 1][..]);
    for (p_s, p_b) in [(0.3, 0.6), (0.6, 0.6), (0.9, 0.9)] {
        let (_, _, _, span) =
            enumerate_moments(&lattice, PercolationMode::Site, p_s, p_b, ends.0, ends.1);
        assert_close(
            span,
            p_s.powi(length as i32),
            1e-12,
            "chain site P(span) = p^L",
        );

        let (_, _, _, span) =
            enumerate_moments(&lattice, PercolationMode::Bond, p_s, p_b, ends.0, ends.1);
        assert_close(
            span,
            p_b.powi(length as i32 - 1),
            1e-12,
            "chain bond P(span) = p^(L-1)",
        );

        let (_, _, _, span) = enumerate_moments(
            &lattice,
            PercolationMode::SiteBond,
            p_s,
            p_b,
            ends.0,
            ends.1,
        );
        let closed_form = p_s.powi(length as i32) * p_b.powi(length as i32 - 1);
        assert_close(
            span,
            closed_form,
            1e-12,
            "chain mixed P(span) = p_s^L p_b^(L-1)",
        );
    }
}

/// Union-find vs flood-fill reference across every lattice family, an
/// arbitrary random graph and all three modes — exhaustive where feasible,
/// seeded random configurations otherwise.
#[test]
fn cluster_stats_matches_independent_bfs_reference() {
    let mut lattices: Vec<(&str, CsrLattice, Vec<usize>, Vec<usize>)> = Vec::new();
    lattices.push(("chain-8", build_chain(8, false), vec![0], vec![7]));
    lattices.push((
        "square-3x3",
        build_square(3, 3, false),
        vec![0, 3, 6],
        vec![2, 5, 8],
    ));
    // Row-major strides (1, 2, 2): x-planes are the even/odd site indices.
    lattices.push((
        "cubic-2x2x2",
        build_hypercubic(
            &[2, 2, 2],
            &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
            false,
        ),
        vec![0, 2, 4, 6],
        vec![1, 3, 5, 7],
    ));
    // Periodic families and the random graph: crossing between an arbitrary
    // site and the last one (topology coverage, not a geometric crossing).
    let triangular = build_triangular(2, 2);
    let triangular_end = triangular.n_sites - 1;
    lattices.push(("triangular-2x2", triangular, vec![0], vec![triangular_end]));
    let honeycomb = build_honeycomb(2, 2);
    let honeycomb_end = honeycomb.n_sites - 1;
    lattices.push(("honeycomb-2x2", honeycomb, vec![0], vec![honeycomb_end]));
    let kagome = build_kagome(2, 2);
    let kagome_end = kagome.n_sites - 1;
    lattices.push(("kagome-2x2", kagome, vec![0], vec![kagome_end]));
    let random = random_graph();
    let random_end = random.n_sites - 1;
    lattices.push(("random-graph", random, vec![0], vec![random_end]));

    for (label, lattice, from, to) in &lattices {
        assert!(
            to.iter().max().copied() < Some(lattice.n_sites),
            "{label}: spanning sets must be in range"
        );
        verify_against_reference(lattice, label, from, to);
    }
}

/// Thermodynamic sanity: the crossing probability never decreases with `p`.
#[test]
fn crossing_probability_is_monotone_in_p() {
    let side = 8;
    let lattice = build_square(side, side, false);
    let from: Vec<usize> = (0..side).map(|row| row * side).collect();
    let to: Vec<usize> = (0..side).map(|row| row * side + side - 1).collect();
    let trials = 20_000;
    for mode in [PercolationMode::Site, PercolationMode::Bond] {
        let mut occupancy = OccupancyState::new(&lattice, mode);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut previous = 0.0;
        for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] {
            let mut spans = 0usize;
            for _ in 0..trials {
                occupancy.resample(p, p, &mut rng);
                spans += usize::from(cluster_stats(&lattice, &occupancy, &from, &to).spanning);
            }
            let current = spans as f64 / trials as f64;
            assert!(
                current >= previous - 0.01,
                "{mode:?}: P(span) decreased from {previous:.4} to {current:.4} at p = {p}"
            );
            previous = current;
        }
        assert!(previous > 0.95, "{mode:?}: P(span) at p = 0.9 should be ~1");
    }
}

/// Mixed-mode monotonicity: the crossing probability never decreases in
/// `p_site` (fixed `p_bond`) nor in `p_bond` (fixed `p_site`).
#[test]
fn mixed_crossing_is_monotone_in_each_probability() {
    let side = 8;
    let lattice = build_square(side, side, false);
    let from: Vec<usize> = (0..side).map(|row| row * side).collect();
    let to: Vec<usize> = (0..side).map(|row| row * side + side - 1).collect();
    let trials = 10_000;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(321);
    let mut occupancy = OccupancyState::new(&lattice, PercolationMode::SiteBond);
    // p_fixed = 0.9 keeps the fixed coordinate deep supercritical, so the
    // varied-at-0.9 endpoint (0.9, 0.9) must saturate the crossing.
    for (p_fixed, varied) in [(0.9, "p_site"), (0.9, "p_bond")] {
        let mut previous = 0.0;
        for p in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let (p_s, p_b) = if varied == "p_site" {
                (p, p_fixed)
            } else {
                (p_fixed, p)
            };
            let mut spans = 0usize;
            for _ in 0..trials {
                occupancy.resample(p_s, p_b, &mut rng);
                spans += usize::from(cluster_stats(&lattice, &occupancy, &from, &to).spanning);
            }
            let current = spans as f64 / trials as f64;
            assert!(
                current >= previous - 0.015,
                "{varied}: P(span) decreased from {previous:.4} to {current:.4} at p = {p}"
            );
            previous = current;
        }
        assert!(previous > 0.95, "{varied}: P(span) at p = 0.9 should be ~1");
    }
}

/// Same base seed → bitwise identical results; different seed → different
/// stream. Guards accidental RNG reuse or thread-order dependence, for pure
/// and mixed modes alike.
#[test]
fn scheduler_runs_reproduce_for_fixed_seeds() {
    let params = |mode: &str| {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("mode", mode);
        if mode == "site-bond" {
            params.set("p_site", 0.6);
            params.set("p_bond", 0.5);
        } else {
            params.set("p", 0.5);
        }
        params
    };
    let config = |seed: u64| RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: 20_000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let mean = |params: &Params, seed: u64| {
        Scheduler::new(RayonBackend::new(1), config(seed))
            .run_one::<PercolationMC>(params)
            .get("Spanning")
            .expect("Spanning measured")
            .mean
    };
    for mode in ["site", "site-bond"] {
        let params = params(mode);
        let first = mean(&params, 424_242);
        assert_eq!(
            first,
            mean(&params, 424_242),
            "{mode}: same seed must reproduce"
        );
        assert_ne!(
            first,
            mean(&params, 424_243),
            "{mode}: different seed must change the stream"
        );
    }
}

/// End-to-end scheduler run: i.i.d. samples reproduce the exact 2x2 site
/// moments within 4 standard errors.
#[test]
fn scheduler_run_matches_exact_moments() {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", L);
    params.set("Ly", L);
    params.set("p", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: 200_000,
        binsize: 100,
        base_seed: 2026,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);

    let lattice = build_square(L, L, false);
    let (exact_max, exact_s2, exact_n, exact_span) =
        enumerate_moments(&lattice, PercolationMode::Site, 0.5, 0.5, &FROM, &TO);
    for (name, exact) in [
        ("MaxCluster", exact_max),
        ("SecondMoment", exact_s2),
        ("NClusters", exact_n),
        ("Spanning", exact_span),
    ] {
        let estimate = results
            .get(name)
            .unwrap_or_else(|| panic!("missing {name}"));
        let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
        assert!(
            z.abs() < 4.0,
            "{name}: mean {} vs exact {exact}, z = {z:.2}",
            estimate.mean
        );
    }
}

/// 1D closed form through the full scheduler stack (adapter, RNG context,
/// binning) for all three modes.
#[test]
fn scheduler_chain_crossing_matches_closed_form() {
    let run = |mode: &str, p_site: f64, p_bond: f64| {
        let mut params = Params::new();
        params.set("lattice_type", "chain");
        params.set("L", 6);
        params.set("mode", mode);
        if mode == "site-bond" {
            params.set("p_site", p_site);
            params.set("p_bond", p_bond);
        } else {
            params.set("p", p_site);
        }
        let config = RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: 100_000,
            binsize: 100,
            base_seed: 77,
            ..Default::default()
        };
        Scheduler::new(RayonBackend::new(1), config)
            .run_one::<PercolationMC>(&params)
            .get("Spanning")
            .expect("Spanning measured")
            .clone()
    };
    let check = |label: &str, estimate: &carlo_rs::Estimate, exact: f64| {
        let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
        assert!(
            z.abs() < 4.0,
            "{label}: {} vs {exact}, z = {z:.2}",
            estimate.mean
        );
    };

    let site = run("site", 0.8, 0.8);
    check("chain site", &site, 0.8_f64.powi(6));
    let bond = run("bond", 0.9, 0.9);
    check("chain bond", &bond, 0.9_f64.powi(5));
    let mixed = run("site-bond", 0.9, 0.8);
    let mixed_exact = 0.9_f64.powi(6) * 0.8_f64.powi(5);
    check("chain site-bond", &mixed, mixed_exact);
}

/// Long-run crossing check: at the square-lattice bond threshold p_c = 1/2
/// (exact by self-duality) the crossing probability of a large L tends to 1/2.
#[test]
#[ignore = "long: bond percolation threshold crossing check (~20 s)"]
fn bond_crossing_at_critical_probability_tends_to_half() {
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", 32);
    params.set("Ly", 32);
    params.set("mode", "bond");
    params.set("p", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 0,
        measurement_sweeps: 200_000,
        binsize: 100,
        base_seed: 7,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<PercolationMC>(&params);
    let spanning = results.get("Spanning").expect("Spanning measured");
    assert!(
        (spanning.mean - 0.5).abs() < 0.06,
        "P(cross) at p_c = {} deviates too far from 1/2",
        spanning.mean
    );
}

/// 3D production check: on a 16^3 cubic lattice well below (0.12) and well
/// above (0.40) the bond threshold p_c ≈ 0.2488 the front-to-back crossing
/// probability must be driven to 0 and 1 respectively.
#[test]
#[ignore = "long: 3D cubic crossing brackets the bond critical region (~10 s)"]
fn cubic_bond_crossing_brackets_the_critical_region() {
    let (lx, ly, lz) = (16, 16, 16);
    // Row-major strides: site = x + y*lx + z*lx*ly.
    let plane = |x: usize| -> String {
        let mut sites = Vec::new();
        for z in 0..lz {
            for y in 0..ly {
                sites.push((x + y * lx + z * lx * ly).to_string());
            }
        }
        sites.join(",")
    };
    let run = |p: f64| {
        let mut params = Params::new();
        params.set("lattice_type", "cubic");
        params.set("Lx", lx);
        params.set("Ly", ly);
        params.set("Lz", lz);
        params.set("pbc", false);
        params.set("mode", "bond");
        params.set("p", p);
        params.set("spanning_from", plane(0));
        params.set("spanning_to", plane(lx - 1));
        let config = RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: 15_000,
            binsize: 100,
            base_seed: 11,
            ..Default::default()
        };
        Scheduler::new(RayonBackend::new(1), config)
            .run_one::<PercolationMC>(&params)
            .get("Spanning")
            .expect("Spanning measured")
            .clone()
    };
    let below = run(0.12);
    assert!(
        below.mean < 0.05,
        "P(cross) at p = 0.12 < p_c should vanish, got {}",
        below.mean
    );
    let above = run(0.40);
    assert!(
        above.mean > 0.95,
        "P(cross) at p = 0.40 > p_c should saturate, got {}",
        above.mean
    );
}

fn assert_close(value: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (value - expected).abs() < tolerance,
        "{label}: {value} vs expected {expected}"
    );
}
