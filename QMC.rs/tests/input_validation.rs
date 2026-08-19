//! Input-validation audit (criterion G) for all four QMC.rs solvers:
//! invalid input must produce an error, never a panic on user data and
//! never silent acceptance of garbage.
//!
//! Covers every scheduler-ready `FromParams` surface (LatticeSpinQmc,
//! ImpurityQmc, OccupationWorldlineQmc, LongitudinalSpinBosonClusterQmc)
//! plus the direct public constructors they build on. Two source-side holes
//! were found by this audit and fixed in `src/lattice/mc.rs` (2026-08-19):
//! an unknown `model` name used to fall through to zero couplings (silently
//! simulating free spins), and an explicitly empty `edges` list used to
//! compile a coupling-free model.

use carlo_rs::{CarloError, FromParams, Params};
use qmc_rs::{
    CavityMode, CsrGraph, EdgeCoupling, EdgeSpec, ImpurityQmc, LatticeSpinQmc,
    LongitudinalSpinBosonClusterQmc, LongitudinalSpinBosonModel, LongitudinalWorldline,
    OccupationSpinBosonModel, OccupationWorldlineQmc, OccupationWorldlineSampler, PowerLawBath,
    RetardedKernel, SingleModeBath, SpinModelBuilder, SpinSpace, TabulatedBath,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = rand_xoshiro::Xoshiro256PlusPlus;

fn seeded() -> Rng {
    Xoshiro256PlusPlus::seed_from_u64(1)
}

/// Replace one field of a valid base configuration and require rejection on
/// the full construction path (the user-visible gate; `validate_params` is
/// only a cheap beta pre-check in QMC.rs).
fn with_param<T: FromParams<Rng = Rng>>(base: &[(&str, &str)], field: &str, value: &str) {
    let mut params = Params::new();
    for (key, raw) in base {
        params.set(key, raw);
    }
    params.set(field, value);
    assert_build_rejects::<T>(&params, field);
}

/// Require the full construction path (which needs an RNG) to reject.
fn assert_build_rejects<T: FromParams<Rng = Rng>>(params: &Params, label: &str) {
    let mut rng = seeded();
    let result = T::from_params(params, &mut rng);
    match result {
        Err(CarloError::InvalidConfig { .. }) => {}
        Err(other) => panic!("{label}: expected InvalidConfig, got {other:?}"),
        Ok(_) => panic!("{label}: invalid configuration was silently accepted"),
    }
}

// ── Lattice directed-loop solver ───────────────────────────────────────────

const LATTICE_BASE: [(&str, &str); 2] = [("beta", "1.0"), ("L", "4")];

#[test]
fn lattice_rejects_invalid_beta_and_model_names() {
    for value in ["-0.5", "0", "NaN", "not-a-number"] {
        with_param::<LatticeSpinQmc>(&LATTICE_BASE, "beta", value);
    }
    // A typo in the model name must be an error, not free spins.
    for value in ["heisenberk", "ising", "TFIM "] {
        with_param::<LatticeSpinQmc>(&LATTICE_BASE, "model", value);
    }
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "topology", "penrose-tiles");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "gauge", "marx");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "scattering", "boltzmann");
}

#[test]
fn lattice_rejects_invalid_geometry() {
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "L", "0");
    let mut params = Params::new();
    params.set("beta", 1.0);
    params.set("topology", "hypercubic");
    params.set("dims", "0,2");
    assert_build_rejects::<LatticeSpinQmc>(&params, "dims=0");

    for (label, edges) in [
        ("empty edge list", ""),
        ("out-of-range endpoint", "0:5"),
        ("self loop", "2:2"),
        ("non-finite weight", "0:1:0:NaN"),
    ] {
        let mut params = Params::new();
        params.set("beta", 1.0);
        params.set("topology", "edges");
        params.set("n_sites", 3);
        params.set("edges", edges);
        assert_build_rejects::<LatticeSpinQmc>(&params, label);
    }
    let mut params = Params::new();
    params.set("beta", 1.0);
    params.set("topology", "edges");
    params.set("n_sites", 3);
    params.set("edges", "0:1, 0:1");
    assert_build_rejects::<LatticeSpinQmc>(&params, "duplicate edge");

    let mut params = Params::new();
    params.set("beta", 1.0);
    params.set("topology", "adjacency");
    params.set("adjacency", "1;");
    assert_build_rejects::<LatticeSpinQmc>(&params, "asymmetric adjacency");
}

#[test]
fn lattice_rejects_invalid_spin_and_coupling_data() {
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "two_s", "0");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "spin", "0.6");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "spin", "-0.5");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "spin", "NaN");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "spins_by_site", "0.5, 0.5");
    with_param::<LatticeSpinQmc>(&LATTICE_BASE, "two_s_by_site", "1, 1");

    // Non-finite couplings are caught on the construction path.
    for (label, field, value) in [
        ("J = NaN", "J", "NaN"),
        ("h_x = inf", "h_x", "inf"),
        ("h_z = -inf", "h_z", "-inf"),
        ("D = NaN", "D", "NaN"),
        ("bond_shift below diagonal", "bond_shift", "-1e9"),
    ] {
        let mut params = Params::new();
        params.set("beta", 1.0);
        params.set("L", 4);
        params.set(field, value);
        assert_build_rejects::<LatticeSpinQmc>(&params, label);
    }
    // Frustrated (non-stoquastic after gauge) geometry is rejected loudly.
    let mut params = Params::new();
    params.set("beta", 1.0);
    params.set("model", "heisenberg");
    params.set("topology", "edges");
    params.set("n_sites", 3);
    params.set("edges", "0:1, 1:2, 2:0");
    params.set("J", "1.0");
    assert_build_rejects::<LatticeSpinQmc>(&params, "frustrated triangle");

    // Sanity: the valid base configuration still builds.
    let mut params = Params::new();
    params.set("beta", 1.0);
    params.set("L", 4);
    let mut rng = seeded();
    assert!(LatticeSpinQmc::from_params(&params, &mut rng).is_ok());
}

// ── Wormhole (spin-boson) solver ───────────────────────────────────────────

const WORMHOLE_BASE: [(&str, &str); 4] = [
    ("beta", "4.0"),
    ("model", "rabi"),
    ("bath", "single"),
    ("omega0", "1.0"),
];

#[test]
fn wormhole_rejects_invalid_beta_model_and_bath() {
    for value in ["-1.0", "NaN", "zero"] {
        with_param::<ImpurityQmc>(&WORMHOLE_BASE, "beta", value);
    }
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "model", "kondo");
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "bath", "lorrentz");
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "omega0", "0");
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "omega0", "-1");
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "omega0", "NaN");
    let rw_base: [(&str, &str); 5] = [
        ("beta", "4.0"),
        ("model", "rw_crw"),
        ("bath", "single"),
        ("omega0", "1.0"),
        ("vertex_scale", "0.3"),
    ];
    with_param::<ImpurityQmc>(&rw_base, "coupling_normalization", "weber");
    with_param::<ImpurityQmc>(&WORMHOLE_BASE, "loop_start", "middle");

    let mut base: Vec<(&str, &str)> = WORMHOLE_BASE.to_vec();
    base[1] = ("model", "xxz");
    base[2] = ("bath", "powerlaw");
    for (field, value) in [("s", "0"), ("omega_c", "-2"), ("s", "NaN")] {
        with_param::<ImpurityQmc>(&base, field, value);
    }
}

#[test]
fn wormhole_rejects_malformed_tabulated_baths_and_couplings() {
    let mut base: Vec<(&str, &str)> = WORMHOLE_BASE.to_vec();
    base[2] = ("bath", "tabulated");
    for (label, omegas, weights) in [
        ("count mismatch", "1.0, 2.0", "1.0"),
        ("negative weight", "1.0", "-1.0"),
        ("non-positive frequency", "0.0", "1.0"),
        ("unparsable", "one, 2.0", "1.0, 1.0"),
        ("zero total mass", "1.0", "0.0"),
    ] {
        let mut params = Params::new();
        for (key, raw) in &base {
            params.set(key, raw);
        }
        params.set("lambda", "0.2");
        params.set("bath_omegas", omegas);
        params.set("bath_weights", weights);
        assert_build_rejects::<ImpurityQmc>(&params, label);
    }

    // A tabulated normalized bath without an explicit coupling is rejected.
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("model", "rabi");
    params.set("bath", "tabulated");
    params.set("bath_omegas", "1.0, 2.0");
    params.set("bath_weights", "1.0, 1.0");
    assert_build_rejects::<ImpurityQmc>(&params, "tabulated bath without lambda");

    // Negative or non-finite effective couplings are rejected.
    for (label, field, value) in [
        ("negative lambda", "lambda", "-0.1"),
        ("NaN lambda", "lambda", "NaN"),
        ("negative lambda_z", "lambda_z", "-1.0"),
        ("NaN C", "C", "NaN"),
    ] {
        let mut params = Params::new();
        for (key, raw) in &WORMHOLE_BASE {
            params.set(key, raw);
        }
        params.set("model", if field == "lambda" { "rabi" } else { "xxz" });
        params.set(field, value);
        assert_build_rejects::<ImpurityQmc>(&params, label);
    }

    // crw_ratio is an rw_crw parameter and must be non-negative there.
    let rw_base: [(&str, &str); 5] = [
        ("beta", "4.0"),
        ("model", "rw_crw"),
        ("bath", "single"),
        ("omega0", "1.0"),
        ("vertex_scale", "0.3"),
    ];
    with_param::<ImpurityQmc>(&rw_base, "crw_ratio", "-0.5");

    // Specifying both constant conventions is ambiguous and rejected.
    let mut params = Params::new();
    for (key, raw) in &WORMHOLE_BASE {
        params.set(key, raw);
    }
    params.set("model", "rw_crw");
    params.set("C", "0.5");
    params.set("diagonal_shift", "0.2");
    assert_build_rejects::<ImpurityQmc>(&params, "C and diagonal_shift together");
}

// ── Occupation (cavity-QED) solver ─────────────────────────────────────────

const OCCUPATION_BASE: [(&str, &str); 4] = [
    ("beta", "2.0"),
    ("kind", "rabi"),
    ("omega0", "1.0"),
    ("cutoff", "4"),
];

#[test]
fn occupation_rejects_invalid_parameters() {
    for value in ["0", "-2", "NaN"] {
        with_param::<OccupationWorldlineQmc>(&OCCUPATION_BASE, "beta", value);
    }
    with_param::<OccupationWorldlineQmc>(&OCCUPATION_BASE, "kind", "dicke");

    for (label, field, value) in [
        ("non-positive omega", "omega0", "0"),
        ("NaN coupling", "g", "NaN"),
        ("zero cutoff", "cutoff", "0"),
        ("one slice", "slices", "1"),
        ("out-of-basis initial state", "initial_state", "9999"),
    ] {
        let mut params = Params::new();
        for (key, raw) in &OCCUPATION_BASE {
            params.set(key, raw);
        }
        params.set(field, value);
        assert_build_rejects::<OccupationWorldlineQmc>(&params, label);
    }

    // Indexed-mode convention with a malformed mode.
    let mut params = Params::new();
    params.set("beta", 2.0);
    params.set("omega_0", "1.0");
    params.set("g_0", "0.2");
    params.set("cutoff_0", "4");
    params.set("omega_1", "-0.5");
    params.set("g_1", "0.2");
    params.set("cutoff_1", "4");
    assert_build_rejects::<OccupationWorldlineQmc>(&params, "negative omega_1");
}

// ── Cluster (longitudinal SB) solver ───────────────────────────────────────

const CLUSTER_BASE: [(&str, &str); 5] = [
    ("beta", "4.0"),
    ("bath", "single"),
    ("omega0", "1.0"),
    ("g", "0.3"),
    ("tunnelling", "0.5"),
];

#[test]
fn cluster_rejects_invalid_parameters() {
    for value in ["-1.0", "NaN", "0"] {
        with_param::<LongitudinalSpinBosonClusterQmc>(&CLUSTER_BASE, "beta", value);
    }
    with_param::<LongitudinalSpinBosonClusterQmc>(&CLUSTER_BASE, "bath", "vibronic");
    with_param::<LongitudinalSpinBosonClusterQmc>(&CLUSTER_BASE, "omega0", "0");

    for (label, field, value) in [
        ("negative tunnelling", "tunnelling", "-0.5"),
        ("NaN bias", "epsilon", "NaN"),
        ("negative lambda", "lambda", "-0.2"),
        ("NaN lambda", "lambda", "NaN"),
    ] {
        let mut params = Params::new();
        for (key, raw) in &CLUSTER_BASE {
            params.set(key, raw);
        }
        params.set(field, value);
        assert_build_rejects::<LongitudinalSpinBosonClusterQmc>(&params, label);
    }

    // Tabulated multi-mode bath without an explicit lambda is rejected.
    let mut params = Params::new();
    params.set("beta", 4.0);
    params.set("bath", "tabulated");
    params.set("bath_omegas", "0.9, 1.7");
    params.set("bath_weights", "1.0, 0.6");
    params.set("tunnelling", "0.8");
    assert_build_rejects::<LongitudinalSpinBosonClusterQmc>(&params, "tabulated without lambda");
}

// ── Direct public constructors ─────────────────────────────────────────────

#[test]
fn direct_constructors_reject_invalid_input() {
    // Local spin space.
    assert!(SpinSpace::uniform(2, 0).is_err());
    assert!(SpinSpace::site_resolved(vec![]).is_err());
    assert!(SpinSpace::site_resolved(vec![1, 0]).is_err());

    // Graphs.
    assert!(CsrGraph::from_edges(0, [EdgeSpec::new(0, 1)]).is_err());
    assert!(CsrGraph::from_edges(2, [EdgeSpec::new(1, 1)]).is_err());
    assert!(CsrGraph::from_edges(2, [EdgeSpec::typed(0, 1, 0, f64::NAN)]).is_err());
    assert!(CsrGraph::from_adjacency(&[vec![1], Vec::new()]).is_err());
    assert!(CsrGraph::hypercubic(&[0, 2], true).is_err());

    // Baths.
    assert!(SingleModeBath::new(0.0).is_err());
    assert!(SingleModeBath::new(f64::NAN).is_err());
    assert!(PowerLawBath::new(0.0, 1.0).is_err());
    assert!(PowerLawBath::new(1.0, 0.0).is_err());
    assert!(TabulatedBath::new(vec![], vec![]).is_err());
    assert!(TabulatedBath::new(vec![1.0], vec![1.0, 2.0]).is_err());
    assert!(TabulatedBath::new(vec![1.0], vec![-1.0]).is_err());
    assert!(TabulatedBath::new(vec![0.0], vec![1.0]).is_err());
    assert!(TabulatedBath::new(vec![1.0], vec![0.0]).is_err());

    // Cavity modes and occupation sampler.
    assert!(CavityMode::new(0.0, 1.0, 3).is_err());
    assert!(CavityMode::new(1.0, f64::NAN, 3).is_err());
    assert!(CavityMode::new(1.0, 1.0, 0).is_err());
    let model = OccupationSpinBosonModel::rabi(1.0, vec![CavityMode::new(1.0, 0.2, 3).unwrap()])
        .expect("valid model");
    assert!(OccupationWorldlineSampler::new(model.clone(), 0.0, 4, 0).is_err());
    assert!(OccupationWorldlineSampler::new(model.clone(), 1.0, 1, 0).is_err());
    assert!(OccupationWorldlineSampler::new(model, 1.0, 4, 999).is_err());
    assert!(OccupationSpinBosonModel::rabi(1.0, vec![]).is_err());

    // Longitudinal cluster model and worldline.
    let bath = qmc_rs::Bath::SingleMode(SingleModeBath::new(1.0).unwrap());
    assert!(
        LongitudinalSpinBosonModel::with_default_quadrature(bath.clone(), -0.1, 0.5, 0.0).is_err()
    );
    assert!(
        LongitudinalSpinBosonModel::with_default_quadrature(bath.clone(), 0.1, -0.5, 0.0).is_err()
    );
    assert!(
        LongitudinalSpinBosonModel::with_default_quadrature(bath.clone(), 0.1, 0.5, f64::NAN)
            .is_err()
    );
    assert!(RetardedKernel::with_default_quadrature(&bath, -1.0).is_err());
    assert!(LongitudinalWorldline::from_kinks(4.0, 1, vec![1.0]).is_err());
    assert!(LongitudinalWorldline::new(0.0, 1).is_err());
    assert!(LongitudinalWorldline::new(4.0, 0).is_err());

    // Spin-lattice builder: non-positive shift margin is rejected.
    let graph = CsrGraph::chain(3, false).unwrap();
    let space = SpinSpace::uniform(3, 1).unwrap();
    assert!(SpinModelBuilder::new(graph.clone(), space.clone())
        .uniform_edge(EdgeCoupling::heisenberg(1.0))
        .shift_margin(0.0)
        .build()
        .is_err());
    assert!(SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(1.0))
        .build()
        .is_ok());
}
