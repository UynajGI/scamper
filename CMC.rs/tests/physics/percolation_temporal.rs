//! X4: contact-process activity on a random regular graph.
//!
//! The contact process is the absorbing-state directed-percolation model
//! (Hinrichsen 2000, archived). On a locally tree-like `k`-regular graph the
//! tree/branching mean-field threshold is `lambda_c = 1/(k-1)`, i.e. 0.5 for
//! `k = 3`, but that is a reference line, not a finite-size certification.
//! On a finite graph the process always eventually absorbs; a finite-time
//! density is not a quasi-stationary order parameter. These tests therefore
//! check finite-time activity ordering and extinction against an independent
//! first-reaction simulator, and report the 0.05-density crossing only as a
//! finite-horizon observable. A longer run (T=100) extinguishes even at
//! lambda=0.55 for N<=2000, so the crossing is not a thermodynamic threshold.

use cmc_rs::{
    contact_ensemble, Bond, BondType, ContactProcessModel, CsrLattice, GraphView,
    UndirectedGraphView,
};
use rand::{Rng, RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Configuration-model random `degree`-regular multigraph.
fn regular_graph(vertices: usize, degree: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut stubs = Vec::with_capacity(vertices * degree);
    for vertex in 0..vertices {
        stubs.extend(std::iter::repeat_n(vertex, degree));
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
    CsrLattice::from_edges(vertices, bonds)
}

/// Independent first-reaction contact process: each possible event gets its
/// own exponential clock and the minimum fires. This is algorithmically
/// distinct from the direct inverse-CDF kernel in the library.
fn first_reaction(
    graph: &CsrLattice,
    lambda: f64,
    initial: &[bool],
    duration: f64,
    rng: &mut impl Rng,
) -> usize {
    let mut active = initial.to_vec();
    let mut time = 0.0;
    let mut events = Vec::new();
    loop {
        events.clear();
        for (vertex, &is_active) in active.iter().enumerate() {
            if is_active {
                events.push((false, vertex, 1.0_f64)); // recovery
            }
        }
        for edge in 0..graph.edge_count() {
            let id = graph.edge_id(edge).expect("dense edge");
            let [left, right] = graph.edge_endpoints(id);
            match (active[left.index()], active[right.index()]) {
                (true, false) => events.push((true, right.index(), lambda)),
                (false, true) => events.push((true, left.index(), lambda)),
                _ => {}
            }
        }
        let total: f64 = events.iter().map(|(_, _, rate)| *rate).sum();
        if total == 0.0 {
            break;
        }
        let mut min_delay = f64::INFINITY;
        let mut chosen = None;
        for &(infect, vertex, rate) in &events {
            let uniform = rng.random::<f64>().max(f64::MIN_POSITIVE);
            let delay = -uniform.ln() / rate;
            if delay < min_delay {
                min_delay = delay;
                chosen = Some((infect, vertex));
            }
        }
        if time + min_delay >= duration {
            break;
        }
        time += min_delay;
        let (infect, vertex) = chosen.expect("non-empty event set");
        active[vertex] = infect;
    }
    active.iter().filter(|&&value| value).count()
}

#[test]
fn independent_first_reaction_agrees_on_extinction_and_activity() {
    let vertices = 80;
    let graph = regular_graph(vertices, 3, 0x4138_0005);
    for (lambda, expect_active) in [(0.2_f64, false), (0.9_f64, true)] {
        let replicas = 6;
        let duration = 12.0;
        let mut direct_survivors = 0;
        let mut reference_survivors = 0;
        for replica in 0..replicas {
            let seed = 0x4138_0006 ^ lambda.to_bits() ^ (replica as u64) << 32;
            let model = ContactProcessModel::new(graph.clone(), lambda, 1.0).expect("model");
            let ensemble = contact_ensemble(&model, 0.05, 1, duration, seed).expect("ensemble");
            direct_survivors += ensemble.survived;

            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ 0xDEAD_BEEF);
            let mut initial = vec![false; vertices];
            for slot in &mut initial {
                *slot = rng.random::<f64>() < 0.05;
            }
            let active = first_reaction(&graph, lambda, &initial, duration, &mut rng);
            reference_survivors += usize::from(active > 0);
        }
        eprintln!(
            "X4 reference lambda={lambda:.2}: direct {direct_survivors}/{replicas}, first-reaction {reference_survivors}/{replicas}"
        );
        if expect_active {
            assert!(
                direct_survivors > 0 && reference_survivors > 0,
                "lambda={lambda} should show finite-time activity in both simulators"
            );
        } else {
            assert_eq!(direct_survivors, 0, "lambda={lambda} should be extinct");
            assert_eq!(reference_survivors, 0, "lambda={lambda} should be extinct");
        }
    }
}

#[test]
fn contact_process_activity_bracket_is_finite_time() {
    let vertices = 200;
    let degree = 3;
    let duration = 15.0;
    let initial_density = 0.05;
    let replicas = 4;
    let graphs = 2;
    let mut rows = Vec::new();
    let mut lambda: f64 = 0.2;
    while lambda <= 1.0 + 1e-12 {
        let mut density_sum = 0.0;
        let mut survived = 0usize;
        let mut total = 0usize;
        for graph_index in 0..graphs {
            let graph = regular_graph(
                vertices,
                degree,
                0x4138_0001 ^ lambda.to_bits() ^ ((graph_index as u64) << 32),
            );
            let model = ContactProcessModel::new(graph.clone(), lambda, 1.0).expect("model");
            assert_eq!(graph.vertex_count(), vertices);
            assert_eq!(graph.edge_count(), vertices * degree / 2);
            let ensemble = contact_ensemble(
                &model,
                initial_density,
                replicas,
                duration,
                0x4138_0002 ^ lambda.to_bits() ^ ((graph_index as u64) << 48),
            )
            .expect("ensemble");
            density_sum += ensemble.mean_density * replicas as f64;
            survived += ensemble.survived;
            total += replicas;
        }
        let density = density_sum / total as f64;
        let survival = survived as f64 / total as f64;
        eprintln!("X4 lambda={lambda:.2}: survival {survival:.3}, mean density {density:.4}");
        rows.push((lambda, density));
        lambda += 0.1;
    }

    let (_, first) = rows[0];
    assert!(
        first < 0.01,
        "subcritical lambda=0.2 has density {first:.4}"
    );
    let (_, last) = *rows.last().expect("non-empty ladder");
    assert!(last > 0.1, "supercritical lambda=1.0 has density {last:.4}");

    let crossing = rows
        .windows(2)
        .find_map(|pair| {
            let (l0, d0) = pair[0];
            let (l1, d1) = pair[1];
            (d0 < 0.05 && d1 >= 0.05).then(|| l0 + (0.05 - d0) / (d1 - d0) * (l1 - l0))
        })
        .expect("density crosses 0.05 inside the ladder");
    eprintln!(
        "X4 finite-time 0.05-density crossing: lambda={crossing:.3} (tree reference 0.5; NOT a threshold)"
    );
    // This is a finite-horizon observable, not a quasi-stationary threshold:
    // at T=100 even lambda=0.55 extinguishes on N<=2000. Keep only the
    // ordering/support checks that are meaningful at this horizon.
    assert!(
        crossing > 0.3 && crossing < 0.7,
        "finite-time crossing {crossing:.3} is outside the activity window"
    );
    let middle = rows[rows.len() / 2].1;
    assert!(
        last > middle,
        "lambda=1.0 density {last:.4} should exceed the mid-ladder density {middle:.4}"
    );
}

#[test]
fn absorbing_runs_report_zero_density() {
    let graph = regular_graph(60, 3, 0x4138_0003);
    let model = ContactProcessModel::new(graph, 0.1, 1.0).expect("model");
    let ensemble = contact_ensemble(&model, 0.01, 8, 20.0, 0x4138_0004).expect("ensemble");
    eprintln!(
        "X4 absorbing control: surviving {}/{}",
        ensemble.survived, ensemble.replicas
    );
    assert_eq!(ensemble.survived, 0);
    assert_eq!(ensemble.survival_fraction, 0.0);
    assert_eq!(ensemble.mean_density, 0.0);
}

#[test]
fn long_horizon_control_does_not_certify_the_crossing() {
    // Deterministic negative control: the lambda=0.55 activity at T=15 is
    // transient on this graph/seed set and has died by T=100. This is not a
    // general extinction theorem; it prevents reading the finite-time
    // crossing as a quasi-stationary threshold.
    let graph = regular_graph(200, 3, 0x4138_0103);
    let model = ContactProcessModel::new(graph, 0.55, 1.0).expect("model");
    let ensemble = contact_ensemble(&model, 0.05, 4, 100.0, 0x4138_0104).expect("ensemble");
    eprintln!(
        "X4 long-horizon control lambda=0.55 T=100: surviving {}/{}",
        ensemble.survived, ensemble.replicas
    );
    assert_eq!(
        ensemble.survived, 0,
        "the finite-time activity must not be presented as a threshold"
    );
}
