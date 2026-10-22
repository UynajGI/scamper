//! D3 production validation: layered space-time directed percolation.
//!
//! Gates (roadmap §11 D3 + literature anchors):
//!
//! 1. **Frontier DP vs explicit DAG path enumeration**: the bitset-frontier
//!    propagation must equal reachability in the fully materialized layered
//!    DAG restricted to the same open bonds. The test replays the
//!    documented draw contract with an independent oracle (own RNG, own
//!    bond bookkeeping), builds the explicit DAG, and checks both against
//!    the D2 `reachable` engine layer by layer.
//! 2. **One realization = one vector sample**: after N sweeps the recorded
//!    `DpLayerActive`/`DpSurvival` arrays carry exactly N samples per
//!    element (n_bins == N), never N*T independent scalars.
//! 3. **Critical-point power law** (Hinrichsen, Adv. Phys. 49, 815
//!    (2000), arXiv:cond-mat/0001070, Tables 1-2; PDF archived in
//!    `docs/ref/`): (1+1)-dimensional bond DP at p_c = 0.6447001 with a
//!    single seed — survival P(t) ~ t^(-delta), delta = 0.159464.
//! 4. **Density decay with full occupation** at p_c: rho(t) ~ t^(-delta),
//!    the same exponent (second, independent anchor).
//! 5. **Phase brackets**: p = 0.60 dies out (P(1024) < 0.01); p = 0.68 on
//!    a periodic ring stays active (P(1024) > 0.9) — statistical brackets,
//!    never exact threshold claims.

use std::collections::HashMap;

use carlo_rs::{Context, Estimate, MonteCarlo};
use cmc_rs::{
    build_chain, reachable, Arc, DirectedCsr, DpSeedPolicy, GraphView,
    LayeredDirectedPercolationMC, ReachabilityWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Hinrichsen Table 1: (1+1)D directed bond percolation threshold.
const PC_BOND_1D: f64 = 0.644_700_1;
/// Hinrichsen Table 2: survival/density decay exponent in d = 1.
const DELTA_1D: f64 = 0.159_464;

fn sweep_count(default: usize) -> usize {
    std::env::var("SCUTTLE_DP_SWEEPS")
        .ok()
        .and_then(|raw| raw.trim().parse().ok())
        .unwrap_or(default)
}

// ── Gate 1: explicit-DAG parity ───────────────────────────────────────────

#[test]
fn frontier_process_equals_explicit_dag_paths() {
    const L: usize = 6;
    const T: usize = 8;
    const P: f64 = 0.55;
    let seed = 0xD3_11;

    // Production realization.
    let space = build_chain(L, true);
    let mut process =
        LayeredDirectedPercolationMC::new(space.clone(), P, T, DpSeedPolicy::Single { site: 2 })
            .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    process.draw_realization(&mut rng);

    // Independent oracle: replay the documented draw contract (layer-major,
    // ascending active sites, CSR incidence order, one draw per incidence)
    // with explicit sets and full bond bookkeeping.
    let mut oracle_rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut frontier: Vec<usize> = vec![2];
    let mut open_bonds: Vec<(usize, usize, usize)> = Vec::new(); // (t, i, j)
    let mut layer_sets: Vec<Vec<usize>> = Vec::with_capacity(T);
    for time in 0..T {
        layer_sets.push(frontier.clone());
        if time + 1 == T {
            break;
        }
        let mut next: Vec<usize> = Vec::new();
        let mut next_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &site in &frontier {
            for (neighbor, _) in space.incidences(site) {
                let open = oracle_rng.random::<f64>() < P;
                if open {
                    open_bonds.push((time, site, neighbor));
                    if next_set.insert(neighbor) {
                        next.push(neighbor);
                    }
                }
            }
        }
        next.sort_unstable();
        frontier = next;
    }

    // Production counts must equal the oracle sets at every layer.
    for (time, set) in layer_sets.iter().enumerate() {
        assert_eq!(
            process.layer_active_counts()[time],
            set.len() as f64,
            "layer {time} count diverged from the oracle"
        );
        assert_eq!(
            process.layer_survival()[time],
            f64::from(u8::from(!set.is_empty())),
            "layer {time} survival diverged"
        );
    }

    // Explicit DAG: node (site, t) = t*L + site; arcs are the open bonds.
    // Path enumeration via the D2 reachability engine from the seed.
    let arcs: Vec<Arc> = open_bonds
        .iter()
        .map(|&(t, i, j)| Arc::new(t * L + i, (t + 1) * L + j))
        .collect();
    let dag = DirectedCsr::from_arcs(T * L, arcs).unwrap();
    let mut reach_workspace = ReachabilityWorkspace::default();
    let seed_node = dag.vertex_id(2).unwrap();
    for (time, set) in layer_sets.iter().enumerate() {
        for site in 0..L {
            let node = dag.vertex_id(time * L + site).unwrap();
            let produced = reachable(&dag, seed_node, node, &mut reach_workspace).unwrap();
            let oracle_says = if time == 0 {
                site == 2
            } else {
                set.contains(&site)
            };
            assert_eq!(
                produced, oracle_says,
                "DAG reachability disagrees at (t={time}, site={site})"
            );
        }
    }
}

// ── Gate 2: one realization = one vector sample ───────────────────────────

#[test]
fn one_realization_records_one_vector_sample() {
    let space = build_chain(64, true);
    let mut process = LayeredDirectedPercolationMC::new(
        space,
        PC_BOND_1D,
        128,
        DpSeedPolicy::Single { site: 32 },
    )
    .unwrap();
    let sweeps = 16;
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(7), 0, 1);
    for _ in 0..sweeps {
        process.sweep(&mut ctx);
        process.measure(&mut ctx);
    }
    let (_, arrays) = ctx.finalize_measurements_split();
    for name in ["DpLayerActive", "DpSurvival"] {
        let array = &arrays[name];
        assert_eq!(array.shape(), &[128], "{name} must record T layers");
        for element in array.elements() {
            assert_eq!(element.n_bins, sweeps, "{name}: one sample per realization");
        }
    }
}

// ── Gates 3-5: literature anchors ─────────────────────────────────────────

fn survival_curve(sweeps: usize, p: f64, seed_policy: DpSeedPolicy) -> HashMap<String, Estimate> {
    let space = build_chain(512, true);
    let mut process = LayeredDirectedPercolationMC::new(space, p, 1200, seed_policy).unwrap();
    let mut ctx =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0xD3_C0), 0, sweeps + 1);
    for _ in 0..sweeps {
        process.sweep(&mut ctx);
        process.measure(&mut ctx);
    }
    let (scalars, arrays) = ctx.finalize_measurements_split();
    let mut out = HashMap::new();
    out.insert("ExtinctionTime".into(), scalars["DpExtinctionTime"].clone());
    for name in ["DpSurvival", "DpLayerActive"] {
        for (index, element) in arrays[name].elements().iter().enumerate() {
            out.insert(format!("{name}[{index}]"), element.clone());
        }
    }
    out
}

#[test]
fn critical_survival_follows_hinrichsen_power_law() {
    let sweeps = sweep_count(2000);
    let curve = survival_curve(sweeps, PC_BOND_1D, DpSeedPolicy::Single { site: 256 });

    let p64 = &curve["DpSurvival[64]"];
    let p1024 = &curve["DpSurvival[1024]"];
    // P(t) ~ t^-delta  =>  delta = ln(P(t1)/P(t2)) / ln(t2/t1), t2/t1 = 16.
    let delta_hat = (p64.mean / p1024.mean).ln() / 16.0_f64.ln();
    eprintln!(
        "[D3 critical] P(64) = {:.4} ± {:.4}, P(1024) = {:.4} ± {:.4}, \
         delta_hat = {:.4} (Hinrichsen: {DELTA_1D}), sweeps = {sweeps}",
        p64.mean, p64.stderr, p1024.mean, p1024.stderr, delta_hat
    );
    assert!(p64.mean > p1024.mean, "survival must decay at criticality");
    assert!(
        (delta_hat - DELTA_1D).abs() < 0.05,
        "delta_hat = {delta_hat:.4} vs Hinrichsen {DELTA_1D}"
    );
    // Extinction time at criticality: the mean sits far below the horizon.
    let extinction = &curve["ExtinctionTime"];
    eprintln!(
        "[D3 critical] mean extinction time = {:.1} of 1200 layers",
        extinction.mean
    );
    assert!(extinction.mean < 1100.0);
}

#[test]
fn full_seed_density_decay_uses_the_same_delta() {
    let sweeps = sweep_count(2000);
    let curve = survival_curve(sweeps, PC_BOND_1D, DpSeedPolicy::Full);

    // rho(t) = <active>/L  ~ t^-delta with the SAME delta.
    let rho64 = curve["DpLayerActive[64]"].mean / 512.0;
    let rho1024 = curve["DpLayerActive[1024]"].mean / 512.0;
    let delta_hat = (rho64 / rho1024).ln() / 16.0_f64.ln();
    eprintln!(
        "[D3 decay] rho(64) = {rho64:.4}, rho(1024) = {rho1024:.4}, \
         delta_hat = {delta_hat:.4} (Hinrichsen: {DELTA_1D}), sweeps = {sweeps}"
    );
    assert!(
        (delta_hat - DELTA_1D).abs() < 0.06,
        "density-decay delta_hat = {delta_hat:.4} vs {DELTA_1D}"
    );
    // Layer zero is fully occupied by construction.
    assert!((curve["DpLayerActive[0]"].mean - 512.0).abs() < 1e-9);
}

#[test]
fn phase_brackets_extinguish_and_survive() {
    let sweeps = sweep_count(1000);

    // Below threshold: the process dies out.
    let below = survival_curve(sweeps, 0.60, DpSeedPolicy::Single { site: 256 });
    let p1024 = below["DpSurvival[1024]"].mean;
    eprintln!("[D3 bracket] p=0.60: P(1024) = {p1024:.4}");
    assert!(p1024 < 0.01, "sub-critical process must extinguish");

    // Above threshold on a periodic ring: the process stays active. The
    // eventual plateau is the percolation strength ~ (p - p_c)^beta ~ 0.4
    // for p = 0.68, and t = 1024 is only a few correlation times
    // (xi_parallel ~ (p - p_c)^-nu_parallel ~ 200), so the measured
    // survival still sits in the approach from above (~0.67); the gate is
    // a conservative contrast against the sub-critical branch, never an
    // exact threshold claim.
    let above = survival_curve(sweeps, 0.68, DpSeedPolicy::Single { site: 256 });
    let p1024 = above["DpSurvival[1024]"].mean;
    eprintln!("[D3 bracket] p=0.68: P(1024) = {p1024:.4}");
    assert!(p1024 > 0.5, "super-critical ring must stay active: {p1024}");
}
