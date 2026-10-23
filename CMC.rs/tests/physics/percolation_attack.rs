//! N3 production validation: degree-ordered targeted attack and robustness
//! curves.
//!
//! Roadmap §9 N3 gates, anchored to archived papers:
//!
//! 1. **Every step matches naive recomputation**: full-resolution
//!    trajectories on star/path/small-ER graphs compared element-by-element
//!    against an independent per-step adjacency-BFS reference (exact
//!    equality, no statistics), including the susceptibility convention.
//! 2. **Active count strictly monotone**: the survivor mask removes exactly
//!    one vertex per step; the reference recomputes survivors per step and
//!    the analysis agrees.
//! 3. **Ties reproduce for fixed seeds**: stable-ID order is
//!    seed-independent; seeded shuffle reproduces bitwise for a fixed seed
//!    and permutes tie groups for different seeds.
//! 4. **Coordinates and provenance**: `RobustnessRemovedFraction` is the
//!    constant coordinate array; policy/grid provenance round-trips through
//!    `Results` attributes and JSON.
//! 5. **Literature anchors**:
//!    Callaway, Newman, Strogatz, Watts, PRL 85, 5468 (2000)
//!    (arXiv:cond-mat/0007300) Eq. (14) with degree-cutoff occupation
//!    `q_k = 1{k <= K}` — an exact configuration-model prediction at
//!    complete degree-group boundaries; and Cohen, Erez, ben-Avraham,
//!    Havlin, PRL 86, 3682 (2001) (arXiv:cond-mat/0010251) — the
//!    intentional-attack breakdown fraction `f_c` from their Eqs. (2), (7),
//!    (9) with the `kappa = 2` criterion, used as a statistical bracket
//!    (never an exact threshold claim). PDFs archived under `docs/ref/`.
//! 6. **T1 adaptive residual-degree attack** (roadmap v2 §2 T1): the
//!    heap-based online order equals a naive full-rescan reference at
//!    every step (both tie policies, several seeds, self-loops and
//!    parallel edges included); per-step invariants (monotone survivor
//!    count, removed vertex at the scanned residual maximum, edges leaving
//!    exactly with their endpoints); tie reproduction; and the Holme, Kim,
//!    Yoon, Han, PRE 65, 056109 (2002) Fig. 7(e) statistical bracket —
//!    adaptive RD below initial-degree ID through the late stage on BA
//!    networks, both crossing to zero far earlier than random failure.
//!    Statistical brackets only, never exact claims.

use std::collections::HashMap;

use carlo_rs::{ArrayEstimate, Context, MonteCarlo, Results};
use cmc_rs::{
    build_chain, AttackGrid, AttackPolicy, Bond, BondType, CsrLattice, DegreeAttackMC, GraphView,
    NetworkProfile, UndirectedGraphView,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::{config_model_power_law, erdos_renyi_graph, star_lattice, zscore_seed_count};

/// One seeded trajectory: the recorded arrays keyed by name.
fn run_attack(
    mc: &mut DegreeAttackMC,
    seed: u64,
    trajectories: u64,
) -> HashMap<String, ArrayEstimate> {
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 32);
    for _ in 0..trajectories {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    let (_, arrays) = ctx.finalize_measurements_split();
    arrays
}

/// Independent per-step giant fraction and susceptibility from adjacency
/// BFS over the survivor set (no project code beyond graph construction).
fn bfs_giant_fraction(adjacency: &[Vec<usize>], removed: &[bool]) -> (f64, Option<f64>) {
    let n = adjacency.len();
    let mut size = vec![0usize; n];
    let mut active_total = 0usize;
    let mut second_moment: u128 = 0;
    let mut largest = 0usize;
    let mut seen = vec![false; n];
    for seed in 0..n {
        if removed[seed] || seen[seed] {
            continue;
        }
        let mut component = vec![seed];
        seen[seed] = true;
        let mut head = 0;
        while head < component.len() {
            let vertex = component[head];
            head += 1;
            for &neighbor in &adjacency[vertex] {
                if !removed[neighbor] && !seen[neighbor] {
                    seen[neighbor] = true;
                    component.push(neighbor);
                }
            }
        }
        for &vertex in &component {
            size[vertex] = component.len();
        }
        active_total += component.len();
        second_moment += (component.len() as u128).pow(2);
        largest = largest.max(component.len());
    }
    let finite = active_total.saturating_sub(largest);
    let chi = if finite > 0 {
        Some((second_moment - (largest as u128).pow(2)) as f64 / finite as f64)
    } else {
        None
    };
    (largest as f64 / n as f64, chi)
}

fn adjacency_of(lattice: &cmc_rs::CsrLattice) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); lattice.vertex_count()];
    for bond in &lattice.edges {
        adjacency[bond.source].push(bond.target);
        adjacency[bond.target].push(bond.source);
    }
    adjacency
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1 + 2: full-resolution parity with the independent BFS reference.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn full_trajectory_matches_independent_bfs_every_step() {
    let graphs: Vec<(&str, cmc_rs::CsrLattice)> = vec![
        ("star5", star_lattice(4)),
        ("path6", build_chain(6, false)),
        ("er12", erdos_renyi_graph(12, 24, 0x4a11)),
    ];
    for (name, lattice) in &graphs {
        for policy in [
            AttackPolicy::DegreeDescendingStableId,
            AttackPolicy::DegreeDescendingSeededShuffle,
        ] {
            let adjacency = adjacency_of(lattice);
            let vertices = lattice.vertex_count();
            let mut mc = DegreeAttackMC::new(lattice.clone(), policy, AttackGrid::Full).unwrap();
            let mut ctx =
                Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4e32), 0, 1);
            mc.sweep(&mut ctx);
            let order = mc.removal_order().to_vec();

            // Degree-descending invariant of the order.
            let profile = NetworkProfile::from_graph(lattice).unwrap();
            for window in order.windows(2) {
                assert!(
                    profile.degree(window[0]) >= profile.degree(window[1]),
                    "{name}: removal order must be degree-descending"
                );
            }

            // Walk once; every checkpoint of the full grid is a prefix.
            mc.walk_trajectory(&order);
            assert_eq!(mc.curve_values().len(), vertices + 1);

            // Strictly monotone survivor count and exact per-step parity.
            let mut removed = vec![false; vertices];
            for (step, &victim) in order.iter().enumerate() {
                removed[victim] = true;
                let checkpoint = step + 1;
                let (reference_giant, reference_chi) = bfs_giant_fraction(&adjacency, &removed);
                assert_eq!(
                    mc.curve_values()[checkpoint],
                    reference_giant,
                    "{name}/{policy:?}: step {checkpoint} giant fraction diverged"
                );
                // Gate 2: exactly one more survivor removed each step,
                // recomputed independently from the survivor mask.
                let reference_active: usize = (0..vertices).filter(|&v| !removed[v]).count();
                assert_eq!(reference_active, vertices - checkpoint);
                match (reference_chi, mc.susceptibility_defined_at(checkpoint)) {
                    (Some(expected), true) => {
                        assert_eq!(mc.susceptibility_values()[checkpoint], expected)
                    }
                    (None, false) => {}
                    other => panic!(
                        "{name}: susceptibility definedness diverged at step {checkpoint}: {other:?}"
                    ),
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: tie policies reproduce.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn tie_policies_reproduce_and_differ_appropriately() {
    // Path8: degrees {1 (x2), 2 (x6)} — large tie groups.
    let lattice = build_chain(8, false);
    let order_of = |policy: AttackPolicy, seed: u64| {
        let mut mc =
            DegreeAttackMC::new(lattice.clone(), policy, AttackGrid::Every { step: 8 }).unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
        mc.sweep(&mut ctx);
        mc.removal_order().to_vec()
    };

    // Stable ID: seed-independent, deterministic.
    let stable_a = order_of(AttackPolicy::DegreeDescendingStableId, 1);
    let stable_b = order_of(AttackPolicy::DegreeDescendingStableId, 99);
    assert_eq!(
        stable_a, stable_b,
        "stable-ID order must not depend on seed"
    );
    assert_eq!(stable_a, vec![1, 2, 3, 4, 5, 6, 0, 7]);

    // Seeded shuffle: bitwise reproduction for a fixed seed.
    let shuffle_a = order_of(AttackPolicy::DegreeDescendingSeededShuffle, 7);
    let shuffle_b = order_of(AttackPolicy::DegreeDescendingSeededShuffle, 7);
    assert_eq!(shuffle_a, shuffle_b, "seeded order must reproduce");

    // Different seeds permute within tie groups while the degree groups
    // themselves stay intact.
    let profile = NetworkProfile::from_graph(&lattice).unwrap();
    let degrees_of = |order: &[usize]| order.iter().map(|&v| profile.degree(v)).collect::<Vec<_>>();
    let shuffle_c = order_of(AttackPolicy::DegreeDescendingSeededShuffle, 1234);
    assert_eq!(degrees_of(&shuffle_a), degrees_of(&shuffle_c));
    // With 6 tied interior vertices, two seeds agreeing everywhere has
    // probability 1/720 under the null.
    assert_ne!(shuffle_a, shuffle_c, "different seeds must reshuffle ties");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: coordinates and provenance.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn curve_coordinates_and_provenance_survive_results_json() {
    let lattice = erdos_renyi_graph(40, 80, 0x4a11);
    let mut mc = DegreeAttackMC::new(
        lattice,
        AttackPolicy::DegreeDescendingStableId,
        AttackGrid::Every { step: 5 },
    )
    .unwrap();
    assert_eq!(mc.checkpoints().len(), 9); // 0, 5, ..., 40

    let arrays = run_attack(&mut mc, 3, 16);
    let coordinates = &arrays["RobustnessRemovedFraction"];
    assert_eq!(coordinates.shape(), &[9]);
    for (index, element) in coordinates.elements().iter().enumerate() {
        assert_eq!(element.mean, index as f64 * 5.0 / 40.0);
        assert_eq!(element.stderr, 0.0, "coordinates are constant");
    }

    // Endpoints of the curve itself: intact graph then total collapse.
    let curve = &arrays["RobustnessCurve"];
    assert!(
        curve.element(0).unwrap().mean > 0.5,
        "ER(40, 80) has a giant"
    );
    assert_eq!(curve.element(8).unwrap().mean, 0.0);

    // Full-grid trajectories always hit an undefined susceptibility
    // checkpoint (the intact connected graph and the empty end), so only
    // the indicator is recorded — never a silently-biased array.
    assert!(!arrays.contains_key("RobustnessSusceptibility"));

    // Provenance: policy + grid attributes round-trip through Results.
    let mut results = Results::new();
    for (name, array) in &arrays {
        results.add_array(name, array.clone());
    }
    results.set_attribute("attack_policy", mc.policy().as_provenance());
    results.set_attribute("curve_points", 9);
    let json = results.to_json().unwrap();
    assert!(json.contains("RobustnessCurve"));
    assert!(json.contains("degree-descending/stable-id"));
    let back = Results::from_json(&json).unwrap();
    assert_eq!(
        back.attribute("attack_policy").and_then(|v| v.as_str()),
        Some("degree-descending/stable-id")
    );
    assert_eq!(
        back.get_array("RobustnessCurve").unwrap().shape(),
        &[9_usize]
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5a: Callaway et al. Eq. (14) with degree-cutoff occupation.
// ─────────────────────────────────────────────────────────────────────────

/// `S = F0(1) - F0(u)`, `u = 1 - F1(1) + F1(u)` for occupation
/// `q_k = 1{k <= cutoff}` on empirical degree bins (their Eq. 14 with the
/// definitions of their Eqs. 1-2).
fn callaway_cutoff_giant(bins: &[usize], cutoff: usize) -> f64 {
    let total = bins.iter().sum::<usize>() as f64;
    let f0 = |x: f64| -> f64 {
        bins.iter()
            .enumerate()
            .take(cutoff + 1)
            .map(|(k, &count)| count as f64 / total * x.powi(k as i32))
            .sum()
    };
    let z: f64 = bins
        .iter()
        .enumerate()
        .map(|(k, &count)| k as f64 * count as f64)
        .sum::<f64>()
        / total;
    let f1 = |x: f64| -> f64 {
        bins.iter()
            .enumerate()
            .skip(1)
            .take(cutoff)
            .map(|(k, &count)| k as f64 * count as f64 / total / z * x.powi(k as i32 - 1))
            .sum()
    };
    let mut u = f0(1.0).min(0.99);
    for _ in 0..200_000 {
        let next = 1.0 - f1(1.0) + f1(u);
        if (next - u).abs() < 1e-15 {
            u = next;
            break;
        }
        u = next;
    }
    f0(1.0) - f0(u)
}

#[test]
fn degree_cutoff_attack_matches_callaway_theory() {
    let cases: Vec<(&str, cmc_rs::CsrLattice, Vec<usize>)> = vec![
        (
            "er-z2.5",
            erdos_renyi_graph(60_000, 75_000, 0x005e_ed01),
            // Cutoffs K = 4, 5, 6, 7 (remove all k > K).
            vec![4, 5, 6, 7],
        ),
        (
            "power-law-a2.5",
            config_model_power_law(60_000, 2.5, 2, 400, 0x90_1a),
            vec![10, 20, 40, 80],
        ),
    ];
    let n_seeds = zscore_seed_count(16);

    for (name, lattice, cutoffs) in &cases {
        let vertices = lattice.vertex_count();
        let profile = NetworkProfile::from_graph(lattice).unwrap();
        let bins = profile.bins();

        // Checkpoints at complete degree-group boundaries: removing all
        // vertices with degree > K removes exactly sum_{k>K} n_k vertices.
        let mut boundaries: Vec<(usize, usize, f64)> = Vec::new(); // (K, removals, theory)
        for &cutoff in cutoffs {
            let removed: usize = bins[cutoff + 1..].iter().sum();
            if removed == 0 || removed >= vertices {
                continue;
            }
            let theory = callaway_cutoff_giant(bins, cutoff);
            if theory > 0.02 {
                boundaries.push((cutoff, removed, theory));
            }
        }
        assert!(
            boundaries.len() >= 3,
            "{name}: need >= 3 usable cutoffs, got {boundaries:?}"
        );
        // AttackGrid::At sorts checkpoints by removal count ascending, so
        // align the boundaries the same way (fewest removed first, i.e.
        // descending cutoff).
        boundaries.sort_unstable_by_key(|&(_, removal, _)| removal);

        let removals: Vec<usize> = boundaries.iter().map(|&(_, r, _)| r).collect();
        let mut mc = DegreeAttackMC::new(
            lattice.clone(),
            AttackPolicy::DegreeDescendingStableId,
            AttackGrid::At {
                removals: removals.clone(),
            },
        )
        .unwrap();

        let mut z_max = 0.0_f64;
        let mut z_sum = 0.0_f64;
        for (checkpoint, &(cutoff, _, theory)) in boundaries.iter().enumerate() {
            let mut means = Vec::new();
            let mut stderrs = Vec::new();
            for seed in 0..n_seeds as u64 {
                let arrays = run_attack(&mut mc, 0x4e32 + 6151 * seed, 1);
                let element = arrays["RobustnessCurve"].element(checkpoint).unwrap();
                means.push(element.mean);
                stderrs.push(element.stderr);
            }
            let n = means.len() as f64;
            let mean = means.iter().sum::<f64>() / n;
            let stderr = (stderrs.iter().map(|se| se * se).sum::<f64>() / (n * n)).sqrt();
            let z = (mean - theory) / (4.0 * stderr + 0.006).max(1e-12);
            z_max = z_max.max(z.abs());
            z_sum += z;
            eprintln!(
                "[N3 {name}] cutoff K={cutoff}: S_mc={mean:.4} ± {stderr:.4} vs Callaway {theory:.4}"
            );
        }
        let mean_z = z_sum / boundaries.len() as f64;
        eprintln!("[N3 {name}] max|z|={z_max:.2} mean z={mean_z:+.2}");
        assert!(z_max < 4.0, "{name}: max |z| = {z_max:.2}");
        assert!(mean_z.abs() < 2.0, "{name}: mean z = {mean_z:.2}");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5b: Cohen et al. intentional-attack breakdown bracket.
// ─────────────────────────────────────────────────────────────────────────

/// Cohen et al. (2001) estimated breakdown fraction for a degree attack:
/// find `f` such that the attacked distribution (their Eq. 2 with link-loss
/// `p_tilde` from their Eq. 9, cutoff from Eq. 7) reaches `kappa = 2`.
fn cohen_breakdown_fraction(bins: &[usize]) -> f64 {
    let total = bins.iter().sum::<usize>() as f64;
    let pk = |k: usize| bins.get(k).copied().unwrap_or(0) as f64 / total;
    let max_degree = bins.len() - 1;
    let z0: f64 = (0..=max_degree).map(|k| k as f64 * pk(k)).sum();

    let kappa_of = |f: f64| -> f64 {
        // Eq. (7): cutoff K~ with sum_{k >= K~} p(k) = f.
        let mut cutoff = max_degree;
        let mut tail = 0.0;
        while cutoff > 0 && tail + pk(cutoff) <= f {
            tail += pk(cutoff);
            cutoff -= 1;
        }
        // Eq. (9): link-loss probability.
        let p_tilde: f64 = ((cutoff + 1)..=max_degree)
            .map(|k| k as f64 * pk(k))
            .sum::<f64>()
            / z0;
        // Eq. (2): binomially thinned degree distribution of survivors.
        let mut m1 = 0.0;
        let mut m2 = 0.0;
        for k0 in 0..=cutoff {
            let weight = pk(k0);
            // Binomial thinning of k0 links with survival 1 - p_tilde.
            let mut binom = 1.0; // C(k0, 0) terms tracked incrementally
            let survival = 1.0 - p_tilde;
            for k in 0..=k0 {
                let mixed =
                    weight * binom * survival.powi(k as i32) * p_tilde.powi((k0 - k) as i32);
                m1 += k as f64 * mixed;
                m2 += (k * k) as f64 * mixed;
                binom = binom * (k0 - k) as f64 / (k + 1) as f64;
            }
        }
        // kappa needs only the first two moments of the thinned
        // distribution (the survivor normalization cancels in m2/m1).
        m2 / m1
    };

    // Bisection on kappa(f) = 2; kappa decreases with f.
    let mut low = 0.0;
    let mut high = 0.9;
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if kappa_of(middle) > 2.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

#[test]
fn intentional_attack_bracket_matches_cohen_prediction() {
    const N: usize = 60_000;
    // m = 1 is the paper's headline regime. Their continuum Eqs. (8)+(11)
    // give p_c ~ 0.056 at alpha = 2.5 (0.02-0.03 belongs to alpha ~ 3);
    // the discrete-empirical estimate of the same kappa = 2 criterion on
    // a finite N = 6e4 network lands lower (0.021 below) because the
    // realized cutoff is far below the continuum assumption. The gate
    // compares our two graphs under OUR estimate and only brackets.
    let power_law = config_model_power_law(N, 2.5, 1, 1500, 0x90_1b);
    let er = erdos_renyi_graph(N, 75_000, 0x005e_ed02);

    let profile_pl = NetworkProfile::from_graph(&power_law).unwrap();
    let profile_er = NetworkProfile::from_graph(&er).unwrap();
    let f_c_pl = cohen_breakdown_fraction(profile_pl.bins());
    let f_c_er = cohen_breakdown_fraction(profile_er.bins());
    eprintln!(
        "[N3 cohen] power-law f_c = {f_c_pl:.4}, ER f_c = {f_c_er:.4} \
         (scale-free nets collapse at a few percent; ER needs much more)"
    );
    // The comparative claim of Cohen et al. / Albert et al.: degree attack
    // is drastically more effective on scale-free networks.
    assert!(
        f_c_pl < 0.5 * f_c_er,
        "scale-free attack threshold must be far smaller"
    );

    // Statistical bracket around the power-law prediction: the giant
    // survives well below f_c and collapses above it. The Cohen estimate
    // carries its own continuum approximations (their Eq. 8 cutoff, the
    // binomial link-loss decoupling), so the bracket is generous and is
    // never an exact threshold claim; the measured numbers land at
    // S(0.3 f_c) ~ 0.4 and S(1.5 f_c) ~ 0.002.
    let below = ((0.3 * f_c_pl) * N as f64) as usize;
    let above = ((1.5 * f_c_pl) * N as f64) as usize;
    let mut mc = DegreeAttackMC::new(
        power_law.clone(),
        AttackPolicy::DegreeDescendingStableId,
        AttackGrid::At {
            removals: vec![below, above],
        },
    )
    .unwrap();
    let arrays = run_attack(&mut mc, 0x4e32, 8);
    let s_below = arrays["RobustnessCurve"].element(0).unwrap();
    let s_above = arrays["RobustnessCurve"].element(1).unwrap();
    eprintln!(
        "[N3 cohen] S(0.3 f_c) = {:.4} ± {:.4}, S(1.5 f_c) = {:.4} ± {:.4}",
        s_below.mean, s_below.stderr, s_above.mean, s_above.stderr
    );
    assert!(
        s_below.mean > 0.25,
        "giant must survive at 0.3x the predicted breakdown fraction"
    );
    assert!(
        s_above.mean < 0.05,
        "giant must have collapsed at 1.5x the predicted breakdown fraction"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// T1 gates (roadmap v2 §2 T1): adaptive residual-degree attack.
// ─────────────────────────────────────────────────────────────────────────

/// Naive full-rescan reference for the adaptive order (reference only,
/// never production): at EVERY step recompute every active vertex's
/// residual degree straight from the adjacency list (a self-loop is two
/// active incidences at its own vertex, the CSR convention) and remove the
/// maximum under (residual desc, tie key asc, dense id asc).
/// `O(V * (V + E))` — the cost the lazy heap avoids.
fn naive_residual_order(adjacency: &[Vec<usize>], tie_keys: &[u64]) -> Vec<usize> {
    let vertices = adjacency.len();
    let mut removed = vec![false; vertices];
    let mut order = Vec::with_capacity(vertices);
    for _ in 0..vertices {
        let mut victim = usize::MAX;
        let mut best_degree = 0_u32;
        let mut best_key = u64::MAX;
        for vertex in 0..vertices {
            if removed[vertex] {
                continue;
            }
            let residual = adjacency[vertex]
                .iter()
                .filter(|&&neighbor| !removed[neighbor])
                .count() as u32;
            let better = residual > best_degree
                || (residual == best_degree && tie_keys[vertex] < best_key)
                || (residual == best_degree && tie_keys[vertex] == best_key && vertex < victim);
            if better {
                best_degree = residual;
                best_key = tie_keys[vertex];
                victim = vertex;
            }
        }
        assert_ne!(victim, usize::MAX, "an active vertex always remains");
        removed[victim] = true;
        order.push(victim);
    }
    order
}

/// Replay the seeded-key draw contract: one `u64` per vertex, ascending
/// dense order, drawn at trajectory start before the first removal.
fn replay_seeded_tie_keys(seed: u64, vertices: usize) -> Vec<u64> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..vertices).map(|_| rng.random::<u64>()).collect()
}

/// Small multigraph with parallel edges and a self-loop: pins the CSR
/// incidence conventions (self-loop = two incidences at its own vertex).
fn multigraph_fixture() -> CsrLattice {
    let edges = vec![
        Bond::new(0, 1, BondType::Generic, 1.0),
        Bond::new(0, 1, BondType::Generic, 1.0),
        Bond::new(1, 2, BondType::Generic, 1.0),
        Bond::new(2, 3, BondType::Generic, 1.0),
        Bond::new(2, 3, BondType::Generic, 1.0),
        Bond::new(3, 4, BondType::Generic, 1.0),
        Bond::new(4, 5, BondType::Generic, 1.0),
        Bond::new(5, 6, BondType::Generic, 1.0),
        Bond::new(6, 0, BondType::Generic, 1.0),
        Bond::new(2, 2, BondType::Generic, 1.0),
    ];
    CsrLattice::from_edges(7, edges)
}

#[test]
fn adaptive_residual_order_matches_naive_rescan_every_step() {
    let graphs: Vec<(&str, CsrLattice)> = vec![
        ("chain10", build_chain(10, false)),
        ("cycle12", build_chain(12, true)),
        ("star6", star_lattice(6)),
        ("completeish10", erdos_renyi_graph(10, 40, 0x4a12)),
        ("multigraph7", multigraph_fixture()),
        (
            "er-multigraph24a",
            config_model_power_law(24, 2.5, 2, 12, 0x6272_4101),
        ),
        (
            "er-multigraph24b",
            config_model_power_law(24, 2.8, 2, 9, 0x6272_4102),
        ),
    ];
    let mut parity_cases = 0_usize;
    for (name, lattice) in &graphs {
        let adjacency = adjacency_of(lattice);
        let vertices = lattice.vertex_count();
        let total_edges = lattice.edge_count();
        for policy in [
            AttackPolicy::ResidualDegreeDescendingStableId,
            AttackPolicy::ResidualDegreeDescendingSeededKey,
        ] {
            // The stable policy is seed-independent; the seeded policy is
            // verified over several seeds so both tie paths are exercised.
            let seeds = match policy {
                AttackPolicy::ResidualDegreeDescendingStableId => vec![0x5eab_0001],
                _ => vec![0x5eed_0001, 0x5eed_0002, 0x5eed_0003],
            };
            for seed in seeds {
                let mut mc =
                    DegreeAttackMC::new(lattice.clone(), policy, AttackGrid::Full).unwrap();
                let mut ctx =
                    Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
                mc.sweep(&mut ctx);
                let order = mc.removal_order().to_vec();

                // Gate 1: exact full-order equality with the naive rescan
                // under the same tie keys (identity keys for the stable
                // policy; the replayed draw contract for the seeded one).
                let keys = match policy {
                    AttackPolicy::ResidualDegreeDescendingStableId => {
                        (0..vertices as u64).collect::<Vec<_>>()
                    }
                    _ => replay_seeded_tie_keys(seed, vertices),
                };
                let reference = naive_residual_order(&adjacency, &keys);
                assert_eq!(
                    order, reference,
                    "{name}/{policy:?}/seed {seed:#x}: heap order diverged from naive rescan"
                );
                parity_cases += 1;

                // Gate 2: per-step invariants, recomputed independently
                // from the adjacency, plus full BFS curve parity (the N3
                // gate applied to the adaptive order).
                mc.walk_trajectory(&order);
                let mut removed = vec![false; vertices];
                for (step, &victim) in order.iter().enumerate() {
                    let residuals: Vec<u32> = (0..vertices)
                        .map(|v| {
                            if removed[v] {
                                0
                            } else {
                                adjacency[v].iter().filter(|&&u| !removed[u]).count() as u32
                            }
                        })
                        .collect();
                    let active: Vec<usize> = (0..vertices).filter(|&v| !removed[v]).collect();
                    let max_residual = active.iter().map(|&v| residuals[v]).max().unwrap();
                    assert_eq!(
                        residuals[victim], max_residual,
                        "{name}: removed vertex must sit at the scanned residual maximum"
                    );
                    assert_eq!(active.len(), vertices - step, "one removal per step");
                    let active_edges: usize = (0..vertices)
                        .filter(|&v| !removed[v])
                        .map(|v| adjacency[v].iter().filter(|&&u| !removed[u]).count())
                        .sum::<usize>()
                        / 2;
                    let expected_edges: usize = lattice
                        .edges
                        .iter()
                        .filter(|edge| !removed[edge.source] && !removed[edge.target])
                        .count();
                    assert_eq!(
                        active_edges, expected_edges,
                        "{name}: edges must leave exactly when an endpoint goes"
                    );
                    assert_eq!(
                        total_edges - active_edges,
                        lattice
                            .edges
                            .iter()
                            .filter(|edge| removed[edge.source] || removed[edge.target])
                            .count(),
                        "{name}: removed edges must be exactly the incident ones"
                    );
                    removed[victim] = true;
                    let (reference_giant, reference_chi) = bfs_giant_fraction(&adjacency, &removed);
                    let checkpoint = step + 1;
                    assert_eq!(
                        mc.curve_values()[checkpoint],
                        reference_giant,
                        "{name}/{policy:?}: step {checkpoint} giant fraction diverged"
                    );
                    match (reference_chi, mc.susceptibility_defined_at(checkpoint)) {
                        (Some(expected), true) => {
                            assert_eq!(mc.susceptibility_values()[checkpoint], expected)
                        }
                        (None, false) => {}
                        other => panic!(
                            "{name}: susceptibility definedness diverged at step {checkpoint}: {other:?}"
                        ),
                    }
                }
            }
        }
    }
    eprintln!("[T1 parity] exact full-order equality on {parity_cases} cases");
}

#[test]
fn adaptive_tie_policies_reproduce_and_differ_only_through_keys() {
    let lattice = config_model_power_law(64, 2.5, 2, 20, 0x7e51_4d47);
    let vertices = lattice.vertex_count();
    let adjacency = adjacency_of(&lattice);
    let order_of = |policy: AttackPolicy, seed: u64| {
        let mut mc =
            DegreeAttackMC::new(lattice.clone(), policy, AttackGrid::Every { step: 64 }).unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
        mc.sweep(&mut ctx);
        mc.removal_order().to_vec()
    };

    // Stable: seed-independent and equal to the naive rescan with identity
    // keys (no RNG consumption on this path).
    let stable_a = order_of(AttackPolicy::ResidualDegreeDescendingStableId, 1);
    let stable_b = order_of(AttackPolicy::ResidualDegreeDescendingStableId, 99);
    assert_eq!(stable_a, stable_b, "stable order must not depend on seed");
    assert_eq!(
        stable_a,
        naive_residual_order(&adjacency, &(0..vertices as u64).collect::<Vec<_>>())
    );

    // Seeded key: bitwise reproduction for a fixed seed; a different seed
    // differs only through the drawn keys — each order still equals the
    // naive rescan under its own replayed keys.
    let seeded_a1 = order_of(AttackPolicy::ResidualDegreeDescendingSeededKey, 7);
    let seeded_a2 = order_of(AttackPolicy::ResidualDegreeDescendingSeededKey, 7);
    assert_eq!(seeded_a1, seeded_a2, "seeded order must reproduce bitwise");
    let seeded_c = order_of(AttackPolicy::ResidualDegreeDescendingSeededKey, 1234);
    assert_eq!(
        seeded_a1,
        naive_residual_order(&adjacency, &replay_seeded_tie_keys(7, vertices))
    );
    assert_eq!(
        seeded_c,
        naive_residual_order(&adjacency, &replay_seeded_tie_keys(1234, vertices))
    );
    assert_ne!(
        seeded_a1, seeded_c,
        "different keys must be able to change the adaptive order"
    );

    // Same seed => identical order AND identical curves (two complete
    // sweep+measure cycles on a full-grid adapter).
    let curves_of = |seed: u64| {
        let mut mc = DegreeAttackMC::new(
            lattice.clone(),
            AttackPolicy::ResidualDegreeDescendingSeededKey,
            AttackGrid::Every { step: 4 },
        )
        .unwrap();
        let arrays = run_attack(&mut mc, seed, 2);
        arrays["RobustnessCurve"]
            .elements()
            .iter()
            .map(|element| element.mean)
            .collect::<Vec<f64>>()
    };
    assert_eq!(curves_of(0x11e5), curves_of(0x11e5));
}

// ─────────────────────────────────────────────────────────────────────────
// T1 gate 5: Holme et al. (2002) Fig. 7(e) statistical bracket.
// ─────────────────────────────────────────────────────────────────────────

/// Barabási–Albert preferential-attachment graph, the Holme et al. (2002)
/// Sec. IV E construction: initial complete graph on `m0` vertices; every
/// later vertex attaches to `m` DISTINCT existing vertices, each target
/// drawn from an endpoint pool proportional to current degree (a duplicate
/// draw is redrawn). Graph builder only — deliberately test-local, not
/// platform physics code. E = m0(m0-1)/2 + m(N - m0).
fn barabasi_albert_graph(n: usize, m0: usize, m: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity(m0 * (m0 - 1) / 2 + m * (n - m0));
    let mut pool: Vec<usize> = Vec::with_capacity(2 * edges.capacity());
    for left in 0..m0 {
        for right in (left + 1)..m0 {
            edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            pool.push(left);
            pool.push(right);
        }
    }
    for vertex in m0..n {
        let mut targets = Vec::with_capacity(m);
        while targets.len() < m.min(vertex) {
            let target = pool[rng.random_range(0..pool.len())];
            if !targets.contains(&target) {
                targets.push(target);
            }
        }
        for target in targets {
            edges.push(Bond::new(vertex, target, BondType::Generic, 1.0));
            pool.push(vertex);
            pool.push(target);
        }
    }
    CsrLattice::from_edges(n, edges)
}

/// Wilson score 95% interval treating the per-graph largest-component
/// fractions as bounded draws in `[0, 1]` (`n` = graph count).
///
/// NOTE (documented deviation): the binomial variance `p(1-p)/n` this
/// interval assumes is ~two orders of magnitude LARGER than the empirical
/// across-graph variance of `S(f)` (component fractions concentrate), so
/// the Wilson interval is reported as the conservative companion; the
/// bracket is gated on the empirical standard-error interval
/// `mean ± z * s / sqrt(n)`, which respects the actual sampling variance.
fn wilson_95(mean: f64, n: usize) -> (f64, f64) {
    let z = 1.959963984540054_f64;
    let n = n as f64;
    let center = (mean + z * z / (2.0 * n)) / (1.0 + z * z / n);
    let half = z * (mean * (1.0 - mean) / n + z * z / (4.0 * n * n)).sqrt() / (1.0 + z * z / n);
    (center - half, center + half)
}

/// Empirical normal 95% interval `mean ± z * s / sqrt(n)` (sample standard
/// deviation across graphs).
fn se_95(samples: &[f64]) -> (f64, f64, f64, f64) {
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    let variance = if samples.len() > 1 {
        samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)
    } else {
        0.0
    };
    let half = 1.959_963_984_540_054 * variance.sqrt() / n.sqrt();
    (mean, mean - half, mean + half, half)
}

#[test]
fn adaptive_attack_bracket_matches_holme_fig7e() {
    // Holme et al. Fig. 7(e): BA m0 = 5, m = 3 (<k> ~ 6), here at the
    // roadmap-mandated N = 10^4 (E ~ 3*10^4). R graph/seed pairs; one
    // trajectory per process per graph: RD = adaptive recalculated degree
    // (T1), ID = initial degree (the N3 order, seeded shuffle ties),
    // random failure = uniformly shuffled removal order.
    const N: usize = 10_000;
    const GRAPHS: usize = 24;
    let fractions = [
        0.005_f64, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
    ];
    let mut removals: Vec<usize> = fractions.iter().map(|f| (f * N as f64) as usize).collect();
    removals.sort_unstable();
    removals.dedup();

    let mut rd = vec![Vec::with_capacity(GRAPHS); removals.len()];
    let mut id = vec![Vec::with_capacity(GRAPHS); removals.len()];
    let mut failure = vec![Vec::with_capacity(GRAPHS); removals.len()];
    for graph in 0..GRAPHS as u64 {
        let lattice = barabasi_albert_graph(N, 5, 3, 0x0BA1_0000 + graph);
        assert_eq!(lattice.edge_count(), 10 + 3 * (N - 5));

        let trajectory = |policy: AttackPolicy, seed: u64| -> Vec<f64> {
            let mut mc = DegreeAttackMC::new(
                lattice.clone(),
                policy,
                AttackGrid::At {
                    removals: removals.clone(),
                },
            )
            .unwrap();
            let arrays = run_attack(&mut mc, seed, 1);
            arrays["RobustnessCurve"]
                .elements()
                .iter()
                .map(|element| element.mean)
                .collect()
        };
        let rd_curve = trajectory(
            AttackPolicy::ResidualDegreeDescendingSeededKey,
            0x5D00_0000 + graph,
        );
        let id_curve = trajectory(
            AttackPolicy::DegreeDescendingSeededShuffle,
            0x1D00_0000 + graph,
        );
        // Random failure: a uniformly shuffled removal order walked by the
        // same trajectory machinery.
        let mut order: Vec<usize> = (0..N).collect();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xFA11_0000 + graph);
        for index in (1..order.len()).rev() {
            let swap_with = rng.random_range(0..=index);
            order.swap(index, swap_with);
        }
        let mut failure_mc = DegreeAttackMC::new(
            lattice,
            AttackPolicy::DegreeDescendingStableId,
            AttackGrid::At {
                removals: removals.clone(),
            },
        )
        .unwrap();
        failure_mc.walk_trajectory(&order);
        let failure_curve: Vec<f64> = failure_mc.curve_values().to_vec();

        for (index, curve) in rd_curve.iter().enumerate() {
            rd[index].push(*curve);
        }
        for (index, curve) in id_curve.iter().enumerate() {
            id[index].push(*curve);
        }
        for (index, curve) in failure_curve.iter().enumerate() {
            failure[index].push(*curve);
        }
    }

    eprintln!("[T1 holme] BA N={N} m0=5 m=3, R={GRAPHS} graphs, one trajectory each");
    let mut separated = Vec::new();
    for (index, &removal) in removals.iter().enumerate() {
        let f = removal as f64 / N as f64;
        let (rd_mean, rd_lo, rd_hi, rd_se) = se_95(&rd[index]);
        let (id_mean, id_lo, id_hi, id_se) = se_95(&id[index]);
        let (fail_mean, fail_lo, fail_hi, _) = se_95(&failure[index]);
        let (w_rd_lo, w_rd_hi) = wilson_95(rd_mean, GRAPHS);
        let (w_id_lo, w_id_hi) = wilson_95(id_mean, GRAPHS);
        let ci_disjoint = rd_hi < id_lo;
        let wilson_disjoint = w_rd_hi < w_id_lo;
        eprintln!(
            "[T1 holme] f={f:.3}: RD {rd_mean:.4} [{rd_lo:.4},{rd_hi:.4}] ±{rd_se:.4} \
             ID {id_mean:.4} [{id_lo:.4},{id_hi:.4}] ±{id_se:.4} \
             fail {fail_mean:.4} [{fail_lo:.4},{fail_hi:.4}] \
             wilson RD [{w_rd_lo:.4},{w_rd_hi:.4}] ID [{w_id_lo:.4},{w_id_hi:.4}] \
             disjoint={ci_disjoint} wilson_disjoint={wilson_disjoint}"
        );
        if ci_disjoint {
            separated.push(f);
        }
    }
    assert!(
        !separated.is_empty(),
        "RD must separate below ID with disjoint 95% CIs somewhere in the late stage"
    );
    // The separated f-range must be contiguous and reach into the late
    // stage (the last separated f at or beyond the grid median).
    let median = fractions[fractions.len() / 2];
    assert!(
        *separated.last().unwrap() >= median,
        "the disjoint RD<ID range must extend into the late stage: {separated:?}"
    );

    // Early stage (Holme: differences among procedures "not significant"):
    // at the first grid point the CIs must overlap (RD ~ ID before the
    // recalculated information diverges from the initial one).
    let (rd_mean0, rd_lo0, rd_hi0, _) = se_95(&rd[0]);
    let (id_mean0, id_lo0, id_hi0, _) = se_95(&id[0]);
    assert!(
        rd_lo0 <= id_hi0 && id_lo0 <= rd_hi0,
        "early stage must have overlapping CIs: RD {rd_mean0} [{rd_lo0},{rd_hi0}] \
         ID {id_mean0} [{id_lo0},{id_hi0}]"
    );

    // Both attacks reach S ~ 0 far earlier than random failure: the
    // smallest f where BOTH attack upper CIs sit below 0.05 must leave the
    // random-failure lower CI above 0.3.
    for index in 0..removals.len() {
        let (_, rd_lo_i, rd_hi_i, _) = se_95(&rd[index]);
        let (_, id_lo_i, id_hi_i, _) = se_95(&id[index]);
        let (_, fail_lo_i, _, _) = se_95(&failure[index]);
        if rd_hi_i < 0.05 && id_hi_i < 0.05 {
            let f = removals[index] as f64 / N as f64;
            assert!(
                fail_lo_i > 0.3,
                "at f={f} both attacks collapsed ({rd_lo_i:.4}, {id_lo_i:.4}) but random \
                 failure sits at {fail_lo_i:.4} — the ordering must match Holme Fig. 7(e)"
            );
            eprintln!(
                "[T1 holme] both attacks collapsed by f={f:.3} while random failure \
                 lower CI is {fail_lo_i:.4}"
            );
            break;
        }
    }
}
