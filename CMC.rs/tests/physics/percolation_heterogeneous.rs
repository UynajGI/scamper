//! Exact and deterministic validation for heterogeneous (N1) percolation laws.
//!
//! The [`StaticLaw`] heterogeneous variants put one validated
//! [`cmc_rs::Probability`] on every vertex and/or every physical edge. This
//! file checks the roadmap §9 N1 deterministic gates through the owned
//! `StaticPercolationMC` adapter (scheduled with `Run::from_parts`, because
//! the `FromParams` schema stays uniform-only):
//!
//! 1. **Exact non-identically-distributed enumeration.** Small graphs are
//!    enumerated configuration-by-configuration with per-element weights
//!    `prod_i p_i^{x_i} (1 - p_i)^{1 - x_i}` (a Poisson-binomial-style
//!    product, not a binomial `p^k (1-p)^(n-k)` factor); scheduler statistics
//!    for distinct `p_i` match within |z| < 4. Hand-derived heterogeneous
//!    chain closed forms and endpoint-deterministic configurations are exact.
//! 2. **Uniform reduction.** A heterogeneous field of identical `p` values
//!    and the uniform law consume RNG in the same order, so identical seeds
//!    reproduce identical `Results` means bitwise, for site, bond, and mixed
//!    laws alike.
//! 3. **Mixed reduction identities.** A heterogeneous mixed law with every
//!    edge probability 1 reduces to the heterogeneous site law, and with
//!    every site probability 1 to the heterogeneous bond law — bitwise at
//!    identical seeds.
//! 4. **Endpoint exactness.** Fields of exact 0/1 probabilities produce
//!    deterministic configurations whose observables are hand-checkable.

use carlo_rs::{Context, Run, RunConfig, RunId, TaskId};
use cmc_rs::{
    analyze, build_chain, build_square, BoundaryQuery, ComponentWorkspace, CsrLattice, GraphView,
    ObservablePlan, Probability, StaticLaw, StaticPercolationMC, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const FROM2: [usize; 2] = [0, 2];
const TO2: [usize; 2] = [1, 3];

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("test probability")
}

fn probabilities(values: &[f64]) -> Vec<Probability> {
    values.iter().map(|&value| probability(value)).collect()
}

/// Adapter observation names compared by the bitwise gates.
const OBSERVED: [&str; 6] = [
    "ActiveVertexCount",
    "ActiveEdgeCount",
    "ComponentCount",
    "LargestSize",
    "RawSecondMoment",
    "BoundaryCrossing0",
];

fn square_plan(lattice: &CsrLattice) -> ObservablePlan {
    ObservablePlan::all(vec![BoundaryQuery::new(lattice, &FROM2, &TO2).unwrap()])
}

/// Run one adapter through `Run::from_parts` and return the `Results`.
///
/// Heterogeneous laws have no `FromParams` schema, so scheduling composes a
/// seeded `Context` directly; identical `seed` values reproduce identical
/// RNG streams across laws with the same draw order.
fn run_adapter(
    lattice: &CsrLattice,
    law: StaticLaw,
    plan: &ObservablePlan,
    seed: u64,
    sweeps: u64,
    binsize: usize,
) -> carlo_rs::Results {
    let adapter = StaticPercolationMC::new(lattice.clone(), law, plan.clone()).unwrap();
    let context = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, binsize);
    let mut run = Run::from_parts(
        context,
        adapter,
        TaskId::new(0),
        RunId::new(0),
        RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: sweeps,
            binsize,
            base_seed: seed,
            ..Default::default()
        },
    );
    run.run(sweeps);
    run.finalize(seed)
}

/// Weight of one configuration under per-element heterogeneous probabilities:
/// `prod p_i^{x_i} (1 - p_i)^{1 - x_i}` over the sampled elements.
fn heterogeneous_weight(mask: usize, values: &[f64], offset: usize) -> f64 {
    values
        .iter()
        .enumerate()
        .map(|(bit, &p)| {
            if (mask >> (offset + bit)) & 1 == 1 {
                p
            } else {
                1.0 - p
            }
        })
        .product()
}

/// Occupancy state for the enumeration reference. Like the uniform
/// enumeration tests in `percolation.rs`, the reference drives the
/// mode-free production `analyze` with mode-specific activity semantics
/// (site: edges join two open sites; bond: open bonds; mixed: both) and
/// differs from the adapter only in the per-element weights.
struct EnumState {
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl EnumState {
    fn site_active(&self, edge: usize, lattice: &CsrLattice) -> bool {
        let bond = &lattice.edges[edge];
        self.bonds[edge] && self.sites[bond.source] && self.sites[bond.target]
    }
}

/// Enumerate `(P(spanning), <LargestSize>, <sum s^2>)` of one heterogeneous
/// mode on `lattice`: site (`p_site` only, bonds implicitly all-active when
/// both endpoints are), bond (`p_bond` only), or mixed (both).
fn enumerate_heterogeneous(
    lattice: &CsrLattice,
    site_p: &[f64],
    bond_p: &[f64],
) -> (f64, f64, f64) {
    let samples_sites = !site_p.is_empty();
    let samples_bonds = !bond_p.is_empty();
    let n_bits = site_p.len() + bond_p.len();
    let queries = [BoundaryQuery::new(lattice, &FROM2, &TO2).unwrap()];
    let mut workspace = ComponentWorkspace::new();
    let mut state = EnumState {
        sites: vec![!samples_sites; lattice.vertex_count()],
        bonds: vec![!samples_bonds; lattice.edge_count()],
    };
    let (mut p_span, mut mean_max, mut mean_s2) = (0.0, 0.0, 0.0);
    for mask in 0..(1usize << n_bits) {
        if samples_sites {
            for (bit, open) in state.sites.iter_mut().enumerate() {
                *open = (mask >> bit) & 1 == 1;
            }
        }
        if samples_bonds {
            let offset = site_p.len();
            for (bit, open) in state.bonds.iter_mut().enumerate() {
                *open = (mask >> (offset + bit)) & 1 == 1;
            }
        }
        let weight = heterogeneous_weight(mask, site_p, 0)
            * heterogeneous_weight(mask, bond_p, site_p.len());
        // Mode-free production analyzer with the same activity semantics the
        // adapter applies: site mode activates every edge between open sites,
        // bond mode every open bond, mixed both conditions.
        struct Activity<'a> {
            state: &'a EnumState,
            lattice: &'a CsrLattice,
            samples_sites: bool,
            samples_bonds: bool,
        }
        impl<'a> cmc_rs::VertexActivity for Activity<'a> {
            fn vertex_count(&self) -> usize {
                self.lattice.vertex_count()
            }
            fn vertex_active(&self, vertex: cmc_rs::VertexId) -> bool {
                self.state.sites[vertex.index()] || !self.samples_sites
            }
        }
        impl<'a> cmc_rs::EdgeActivity for Activity<'a> {
            fn edge_count(&self) -> usize {
                self.lattice.edge_count()
            }
            fn edge_active(&self, edge: cmc_rs::EdgeId) -> bool {
                if !self.samples_bonds {
                    return true;
                }
                self.state.site_active(edge.index(), self.lattice)
            }
        }
        let activity = Activity {
            state: &state,
            lattice,
            samples_sites,
            samples_bonds,
        };
        let result = analyze(lattice, &activity, &activity, &queries, &mut workspace).unwrap();
        let summary = result.summary();
        p_span += weight * f64::from(u8::from(result.query_outcomes()[0]));
        mean_max += weight * summary.largest_component_size as f64;
        mean_s2 += weight * summary.raw_second_moment as f64;
    }
    (p_span, mean_max, mean_s2)
}

fn assert_close(value: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (value - expected).abs() < tolerance,
        "{label}: {value} vs expected {expected}"
    );
}

/// The enumerated weights must form a probability distribution.
#[test]
fn heterogeneous_enumeration_weights_sum_to_one() {
    let site_p = [0.1, 0.9, 0.4, 0.65];
    let bond_p = [0.2, 0.75, 0.5, 0.35];
    for (site, bond) in [
        (site_p.to_vec(), Vec::new()),
        (Vec::new(), bond_p.to_vec()),
        (site_p.to_vec(), bond_p.to_vec()),
    ] {
        let n_bits = site.len() + bond.len();
        let mut total = 0.0;
        for mask in 0..(1usize << n_bits) {
            total += heterogeneous_weight(mask, &site, 0)
                * heterogeneous_weight(mask, &bond, site.len());
        }
        assert_close(total, 1.0, 1e-12, "mixed weight normalization");
    }
}

/// Hand-derived heterogeneous chain closed forms through the scheduler: an
/// open chain spans only when every site (site mode), every bond (bond mode)
/// or both (mixed) is open, so `P = prod p_i` element-wise.
#[test]
fn heterogeneous_chain_crossing_follows_element_product() {
    let length = 6usize;
    let lattice = build_chain(length, false);
    let plan = ObservablePlan::all(vec![
        BoundaryQuery::new(&lattice, &[0], &[length - 1]).unwrap()
    ]);
    let site_p = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4];
    let bond_p = [0.75, 0.65, 0.55, 0.45, 0.35];
    let law = StaticLaw::site_heterogeneous(probabilities(&site_p));
    let results = run_adapter(&lattice, law, &plan, 31, 200_000, 100);
    let exact = site_p.iter().product::<f64>();
    let estimate = results.get("BoundaryCrossing0").unwrap();
    let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
    assert!(
        z.abs() < 4.0,
        "site: {} vs {exact}, z = {z:.2}",
        estimate.mean
    );

    let law = StaticLaw::bond_heterogeneous(probabilities(&bond_p));
    let results = run_adapter(&lattice, law, &plan, 32, 200_000, 100);
    let exact = bond_p.iter().product::<f64>();
    let estimate = results.get("BoundaryCrossing0").unwrap();
    let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
    assert!(
        z.abs() < 4.0,
        "bond: {} vs {exact}, z = {z:.2}",
        estimate.mean
    );

    let law = StaticLaw::mixed_heterogeneous(probabilities(&site_p), probabilities(&bond_p));
    let results = run_adapter(&lattice, law, &plan, 33, 200_000, 100);
    let exact = site_p.iter().product::<f64>() * bond_p.iter().product::<f64>();
    let estimate = results.get("BoundaryCrossing0").unwrap();
    let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
    assert!(
        z.abs() < 4.0,
        "mixed: {} vs {exact}, z = {z:.2}",
        estimate.mean
    );
}

/// Exact non-identically-distributed enumeration on the 2x2 open square with
/// four distinct probabilities per domain; scheduler statistics for spanning,
/// largest size, and raw second moment match within |z| < 4 for every mode.
#[test]
fn heterogeneous_scheduler_matches_exact_enumeration() {
    let lattice = build_square(2, 2, false);
    let plan = square_plan(&lattice);
    let site_p = [0.1, 0.9, 0.4, 0.65];
    let bond_p = [0.2, 0.75, 0.5, 0.35];

    let cases: [(&str, Vec<f64>, Vec<f64>); 3] = [
        ("site", site_p.to_vec(), Vec::new()),
        ("bond", Vec::new(), bond_p.to_vec()),
        ("mixed", site_p.to_vec(), bond_p.to_vec()),
    ];
    for (index, (mode, site, bond)) in cases.iter().enumerate() {
        let (exact_span, exact_max, exact_s2) = enumerate_heterogeneous(&lattice, site, bond);
        let law = match *mode {
            "site" => StaticLaw::site_heterogeneous(probabilities(site)),
            "bond" => StaticLaw::bond_heterogeneous(probabilities(bond)),
            _ => StaticLaw::mixed_heterogeneous(probabilities(site), probabilities(bond)),
        };
        for run_index in 0..4u64 {
            let results = run_adapter(
                &lattice,
                law.clone(),
                &plan,
                0x4845_5452_0000u64 + ((index * 16) as u64 + run_index) * 0x9E37,
                200_000,
                100,
            );
            for (name, exact) in [
                ("BoundaryCrossing0", exact_span),
                ("LargestSize", exact_max),
                ("RawSecondMoment", exact_s2),
            ] {
                let estimate = results.get(name).unwrap();
                let z = (estimate.mean - exact) / estimate.stderr.max(1e-12);
                assert!(
                    z.abs() < 4.0,
                    "{mode}/{name} run {run_index}: {} vs {exact}, z = {z:.2}",
                    estimate.mean
                );
            }
        }
    }
}

/// Uniform reduction, bitwise: a heterogeneous field of identical `p` values
/// consumes RNG exactly like the uniform law (same per-element draw order),
/// so identical seeds reproduce identical observable means bit-for-bit. This
/// holds for the exact endpoints as well (neither path consumes RNG).
#[test]
fn heterogeneous_uniform_reduction_is_bitwise_for_identical_seeds() {
    let lattice = build_square(3, 3, false);
    let from: Vec<usize> = (0..3).map(|row| row * 3).collect();
    let to: Vec<usize> = (0..3).map(|row| row * 3 + 2).collect();
    let plan = ObservablePlan::all(vec![BoundaryQuery::new(&lattice, &from, &to).unwrap()]);
    let vertices = lattice.vertex_count();
    let edges = lattice.edge_count();

    for p in [0.0, 0.3, 0.5, 1.0] {
        let uniform = StaticLaw::site(p).unwrap();
        let heterogeneous = StaticLaw::site_heterogeneous(vec![probability(p); vertices]);
        let left = run_adapter(&lattice, uniform, &plan, 0x554E_4946, 40_000, 100);
        let right = run_adapter(&lattice, heterogeneous, &plan, 0x554E_4946, 40_000, 100);
        for name in OBSERVED {
            assert_eq!(
                left.get(name).unwrap().mean.to_bits(),
                right.get(name).unwrap().mean.to_bits(),
                "site uniform reduction at p = {p} differs for {name}"
            );
        }

        let uniform = StaticLaw::bond(p).unwrap();
        let heterogeneous = StaticLaw::bond_heterogeneous(vec![probability(p); edges]);
        let left = run_adapter(&lattice, uniform, &plan, 0x554E_4946, 40_000, 100);
        let right = run_adapter(&lattice, heterogeneous, &plan, 0x554E_4946, 40_000, 100);
        for name in OBSERVED {
            assert_eq!(
                left.get(name).unwrap().mean.to_bits(),
                right.get(name).unwrap().mean.to_bits(),
                "bond uniform reduction at p = {p} differs for {name}"
            );
        }

        let uniform = StaticLaw::mixed(p, 1.0 - p).unwrap();
        let heterogeneous = StaticLaw::mixed_heterogeneous(
            vec![probability(p); vertices],
            vec![probability(1.0 - p); edges],
        );
        let left = run_adapter(&lattice, uniform, &plan, 0x554E_4946, 40_000, 100);
        let right = run_adapter(&lattice, heterogeneous, &plan, 0x554E_4946, 40_000, 100);
        for name in OBSERVED {
            assert_eq!(
                left.get(name).unwrap().mean.to_bits(),
                right.get(name).unwrap().mean.to_bits(),
                "mixed uniform reduction at p = {p} differs for {name}"
            );
        }
    }
}

/// Mixed reduction identities, bitwise: a heterogeneous mixed law whose edge
/// field is all-ones reduces to the heterogeneous site law, and whose site
/// field is all-ones to the heterogeneous bond law. Exact endpoint elements
/// skip RNG individually, so the draw order matches and identical seeds
/// reproduce identical results.
#[test]
fn heterogeneous_mixed_reduces_to_pure_laws_bitwise() {
    let lattice = build_square(3, 3, false);
    let from: Vec<usize> = (0..3).map(|row| row * 3).collect();
    let to: Vec<usize> = (0..3).map(|row| row * 3 + 2).collect();
    let plan = ObservablePlan::all(vec![BoundaryQuery::new(&lattice, &from, &to).unwrap()]);
    let site_p = [0.15, 0.85, 0.4, 0.6, 0.25, 0.75, 0.5, 0.35, 0.65];
    let bond_p = [
        0.2, 0.45, 0.7, 0.3, 0.55, 0.8, 0.4, 0.65, 0.15, 0.5, 0.35, 0.6,
    ];

    let site = StaticLaw::site_heterogeneous(probabilities(&site_p));
    let mixed_as_site = StaticLaw::mixed_heterogeneous(
        probabilities(&site_p),
        vec![probability(1.0); lattice.edge_count()],
    );
    let left = run_adapter(&lattice, site, &plan, 0x4D49_5845, 40_000, 100);
    let right = run_adapter(&lattice, mixed_as_site, &plan, 0x4D49_5845, 40_000, 100);
    for name in OBSERVED {
        assert_eq!(
            left.get(name).unwrap().mean.to_bits(),
            right.get(name).unwrap().mean.to_bits(),
            "p_bond = 1 identity differs for {name}"
        );
    }

    let bond = StaticLaw::bond_heterogeneous(probabilities(&bond_p));
    let mixed_as_bond = StaticLaw::mixed_heterogeneous(
        vec![probability(1.0); lattice.vertex_count()],
        probabilities(&bond_p),
    );
    let left = run_adapter(&lattice, bond, &plan, 0x4D49_5845, 40_000, 100);
    let right = run_adapter(&lattice, mixed_as_bond, &plan, 0x4D49_5845, 40_000, 100);
    for name in OBSERVED {
        assert_eq!(
            left.get(name).unwrap().mean.to_bits(),
            right.get(name).unwrap().mean.to_bits(),
            "p_site = 1 identity differs for {name}"
        );
    }
}

/// Exact endpoint fields are deterministic: every sweep redraws the same
/// configuration, so the ensemble means equal one hand-checked analysis.
#[test]
fn heterogeneous_endpoint_fields_are_deterministic() {
    let lattice = build_chain(6, false);
    let plan = ObservablePlan::all(vec![BoundaryQuery::new(&lattice, &[0], &[5]).unwrap()]);
    // Sites [1 0 1 0 1 0]: three isolated occupied sites, no span.
    let law = StaticLaw::site_heterogeneous(probabilities(&[1.0, 0.0, 1.0, 0.0, 1.0, 0.0]));
    let results = run_adapter(&lattice, law, &plan, 7, 100, 100);
    assert_eq!(results.get("ActiveVertexCount").unwrap().mean, 3.0);
    assert_eq!(results.get("ActiveEdgeCount").unwrap().mean, 0.0);
    assert_eq!(results.get("ComponentCount").unwrap().mean, 3.0);
    assert_eq!(results.get("LargestSize").unwrap().mean, 1.0);
    assert_eq!(results.get("RawSecondMoment").unwrap().mean, 3.0);
    assert_eq!(results.get("BoundaryCrossing0").unwrap().mean, 0.0);

    // Bonds [1 1 0 1 1]: components {0,1,2} and {3,4,5}, no span.
    let law = StaticLaw::bond_heterogeneous(probabilities(&[1.0, 1.0, 0.0, 1.0, 1.0]));
    let results = run_adapter(&lattice, law, &plan, 8, 100, 100);
    assert_eq!(results.get("ActiveVertexCount").unwrap().mean, 6.0);
    assert_eq!(results.get("ActiveEdgeCount").unwrap().mean, 4.0);
    assert_eq!(results.get("ComponentCount").unwrap().mean, 2.0);
    assert_eq!(results.get("LargestSize").unwrap().mean, 3.0);
    assert_eq!(results.get("RawSecondMoment").unwrap().mean, 18.0);
    assert_eq!(results.get("BoundaryCrossing0").unwrap().mean, 0.0);

    // Mixed sites [1 0 1 1 1 1] x bonds [1 1 0 1 1]: site 1 is closed and
    // bond (2,3) is closed, so the occupied sites form {0}, {2} and {3,4,5};
    // only bonds (3,4) and (4,5) join open sites.
    let law = StaticLaw::mixed_heterogeneous(
        probabilities(&[1.0, 0.0, 1.0, 1.0, 1.0, 1.0]),
        probabilities(&[1.0, 1.0, 0.0, 1.0, 1.0]),
    );
    let results = run_adapter(&lattice, law, &plan, 9, 100, 100);
    assert_eq!(results.get("ActiveVertexCount").unwrap().mean, 5.0);
    assert_eq!(results.get("ActiveEdgeCount").unwrap().mean, 2.0);
    assert_eq!(results.get("ComponentCount").unwrap().mean, 3.0);
    assert_eq!(results.get("LargestSize").unwrap().mean, 3.0);
    assert_eq!(results.get("RawSecondMoment").unwrap().mean, 11.0);
    assert_eq!(results.get("BoundaryCrossing0").unwrap().mean, 0.0);
}

/// Adapter-level typed rejection of heterogeneous domain mismatches and of
/// unvalidated probability values (the unit-level twin lives next to the
/// adapter; this pins the same behavior through the scheduling surface).
#[test]
fn heterogeneous_invalid_domains_are_rejected_at_construction() {
    let lattice = build_square(2, 2, false);
    let plan = square_plan(&lattice);
    // Wrong vertex-count field.
    let law = StaticLaw::site_heterogeneous(probabilities(&[0.5; 3]));
    assert!(
        StaticPercolationMC::new(lattice.clone(), law, plan.clone()).is_err(),
        "short site field must be rejected"
    );
    // Wrong edge-count field.
    let law = StaticLaw::bond_heterogeneous(probabilities(&[0.5; 5]));
    assert!(
        StaticPercolationMC::new(lattice.clone(), law, plan.clone()).is_err(),
        "long bond field must be rejected"
    );
    // Mixed fields validated one domain at a time.
    let law = StaticLaw::mixed_heterogeneous(probabilities(&[0.5; 4]), Vec::new());
    assert!(
        StaticPercolationMC::new(lattice.clone(), law, plan).is_err(),
        "empty mixed bond field must be rejected"
    );
    // Unvalidated values cannot be expressed: the constructors take
    // `Probability`, and every f64 entry point validates at `Probability::new`.
    for value in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        assert!(Probability::new(value).is_err(), "accepted {value}");
    }
}
