//! X2: quantitative two-layer cascade bracket and the cascade adapter.
//!
//! Anchor: Buldyrev et al., Nature 464, 1025 (2010) and the archived Saberi
//! review §V. For two fully interdependent Poisson random graphs with mean
//! degree `c`, the mutual giant component emerges at `c* ~ 2.4554`. This is
//! a finite-size bracket: the test asserts subcritical collapse and
//! supercritical survival, checks that the configuration-model layers
//! realize the requested mean degree, and brackets the 0.2-fraction crossing
//! around the literature value. The 0.2 level is a finite-size proxy: for a
//! first-order transition the crossing approaches `c*` from below as `N`
//! grows, so this is not a four-digit reproduction.

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    Bond, BondType, CsrLattice, MultiplexCascadeMC, MultiplexError, MultiplexNetwork,
    UndirectedGraphView,
};
use rand::{Rng, RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// One Poisson degree draw (Knuth) with mean `mean`, clamped to `n - 1`.
fn poisson_degree(rng: &mut impl Rng, mean: f64, n: usize) -> usize {
    let limit = (-mean).exp();
    let mut product = 1.0;
    let mut degree = 0usize;
    loop {
        product *= rng.random::<f64>();
        if product <= limit {
            return degree.min(n - 1);
        }
        degree += 1;
    }
}

/// Configuration-model Poisson random graph (self-loops and parallel edges
/// are allowed, as in the standard multigraph ensemble).
fn er_layer(n: usize, mean: f64, rng: &mut impl Rng) -> CsrLattice {
    let mut stubs = Vec::new();
    for vertex in 0..n {
        let degree = poisson_degree(rng, mean, n);
        stubs.extend(std::iter::repeat_n(vertex, degree));
    }
    if stubs.len() % 2 == 1 {
        stubs.pop();
    }
    for index in (1..stubs.len()).rev() {
        let swap = rng.random_range(0..=index);
        stubs.swap(index, swap);
    }
    let bonds = stubs
        .as_chunks::<2>()
        .0
        .iter()
        .map(|pair| Bond::new(pair[0], pair[1], BondType::Generic, 1.0))
        .collect::<Vec<_>>();
    CsrLattice::from_edges(n, bonds)
}

/// Mutual fraction and the realized mean degree averaged over both layers
/// (`(E_a + E_b) / N`; each layer alone has `2 E / N`).
fn cascade_fraction(vertices: usize, mean: f64, seed: u64) -> (f64, f64) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let layer_a = er_layer(vertices, mean, &mut rng);
    let layer_b = er_layer(vertices, mean, &mut rng);
    let realised = (layer_a.edge_count() + layer_b.edge_count()) as f64 / vertices as f64;
    let network = MultiplexNetwork::new(layer_a, layer_b).expect("shared vertex set");
    (
        network.mutual_component().expect("cascade").fraction(),
        realised,
    )
}

#[test]
fn two_layer_er_cascade_brackets_the_buldyrev_threshold() {
    let vertices = 1024;
    let seeds = 10u64;
    let mut rows = Vec::new();
    let mut mean_degree = 1.6;
    while mean_degree <= 3.4 + 1e-12 {
        let mut sum = 0.0;
        let mut realised_sum = 0.0;
        for seed in 0..seeds {
            let (fraction, realised) = cascade_fraction(
                vertices,
                mean_degree,
                0x4137_0001 ^ (mean_degree.to_bits()) ^ seed.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            );
            assert!(
                (realised - mean_degree).abs() < 0.25,
                "configuration model realized mean degree {realised:.3}, requested {:.3}",
                mean_degree
            );
            sum += fraction;
            realised_sum += realised;
        }
        let fraction = sum / seeds as f64;
        eprintln!(
            "X2 mean degree {mean_degree:.2} (realized {:.3}): mutual fraction {fraction:.4}",
            realised_sum / seeds as f64
        );
        rows.push((mean_degree, fraction));
        mean_degree += 0.1;
    }

    let (_, first) = rows[0];
    assert!(first < 0.01, "subcritical c=1.6 has a giant: {first:.4}");
    let (_, last) = *rows.last().expect("non-empty ladder");
    assert!(last > 0.2, "supercritical c=3.4 has no giant: {last:.4}");

    let crossing = rows
        .windows(2)
        .find_map(|pair| {
            let (c0, f0) = pair[0];
            let (c1, f1) = pair[1];
            (f0 < 0.2 && f1 >= 0.2).then(|| c0 + (0.2 - f0) / (f1 - f0) * (c1 - c0))
        })
        .expect("the mutual fraction crosses 0.2 inside the ladder");
    eprintln!("X2 0.2-fraction crossing: c={crossing:.3} (Buldyrev c*=2.4554)");
    assert!(
        (crossing - 2.4554).abs() < 0.5,
        "crossing c={crossing:.3} is outside the 0.5 bracket around 2.4554"
    );
    assert!(
        crossing > 2.0 && crossing < 3.2,
        "crossing c={crossing:.3} is outside the finite-size window"
    );
}

#[test]
fn cascade_adapter_records_survival_and_mutual_fraction() {
    // Layer A is a 3-path, layer B a 2-path: the mutual giant settles on
    // {0, 1}, i.e. 2/3 of the shared vertex set.
    let layer_a = CsrLattice::from_edges(
        3,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
        ],
    );
    let layer_b = CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
    let network = MultiplexNetwork::new(layer_a, layer_b).expect("shared vertex set");
    assert_eq!(network.vertices(), 3);

    let mut adapter = MultiplexCascadeMC::new(network, 1.0).expect("adapter");
    let mut context =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4137_0002), 0, 1);
    for _ in 0..4 {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let results = context.finalize_measurements();
    assert_eq!(results["MultiplexMutualFraction"].mean, 2.0 / 3.0);
    assert_eq!(results["MultiplexMutualVertices"].mean, 2.0);
    assert_eq!(results["MultiplexInitialVertices"].mean, 3.0);
    assert!(results["MultiplexCascadeIterations"].mean >= 2.0);

    let mut dead = MultiplexCascadeMC::new(
        MultiplexNetwork::new(
            CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]),
            CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]),
        )
        .expect("shared vertex set"),
        0.0,
    )
    .expect("adapter");
    let mut context =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4137_0003), 0, 1);
    dead.sweep(&mut context);
    dead.measure(&mut context);
    let results = context.finalize_measurements();
    assert_eq!(results["MultiplexMutualFraction"].mean, 0.0);
    assert_eq!(results["MultiplexInitialVertices"].mean, 0.0);

    let error = MultiplexCascadeMC::new(
        MultiplexNetwork::new(
            CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]),
            CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]),
        )
        .expect("shared vertex set"),
        1.5,
    )
    .err()
    .expect("survival outside [0, 1] is rejected");
    assert!(matches!(error, MultiplexError::Probability(source) if source.value() == 1.5));
}

#[test]
fn initial_activity_must_cover_the_shared_vertex_set() {
    let layer = || {
        CsrLattice::from_edges(
            3,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
            ],
        )
    };
    let network = MultiplexNetwork::new(layer(), layer()).expect("shared vertex set");
    assert!(matches!(
        network.mutual_component_from(&[true, true]),
        Err(MultiplexError::ActivityDomain {
            expected: 3,
            actual: 2
        })
    ));
    // A single surviving connected pair fixes the result immediately.
    let mutual = network
        .mutual_component_from(&[false, true, true])
        .expect("cascade");
    assert_eq!(mutual.active(), &[false, true, true]);
    assert_eq!(mutual.fraction(), 2.0 / 3.0);
}
