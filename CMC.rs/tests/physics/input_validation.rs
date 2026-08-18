//! Input-validation audit (criterion G): invalid input → error, never
//! silent garbage and never an unintended panic on user data.
//!
//! Every scheduler-ready solver's `FromParams` path (and the direct
//! constructors of the non-scheduler kernels) is fed a representative set of
//! malformed configurations. Each must return `CarloError::InvalidConfig`
//! (or the kernel's own error type), not panic and not silently accept.
//!
//! This file closes the per-solver rejection-test gap identified in the
//! 2026-08-19 production audit: source-side validation was already extensive
//! but only sparsely covered by tests, and the Kawasaki/BKL adapters parsed
//! `beta`/`J` without checks — a `J = "NaN"` parameter reached an
//! assert-backed constructor as a panic (fixed in `dynamics/mc.rs`).

use carlo_rs::{CarloError, FromParams, Params};
use cmc_rs::{
    ClassicalMC, ContinuousHeatBathCore, HardSphereEventChainMC, HeatBathCore, HeisenbergModel,
    HybridCore, IsingGraphWormMC, IsingModel, IsingWangLandau, KawasakiIsingMC, KineticIsingBklMC,
    LennardJonesMuVt, LennardJonesNpt, LennardJonesNvt, MetropolisCore, MicrocanonicalCore,
    MolecularMetropolisCore, MoleculeTopology, MultiSpinIsing, PottsModel, SWCore,
    WangLandauConfig, WolffCore, XYModel,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type MetroIsing = ClassicalMC<IsingModel, MetropolisCore>;
type WolffIsing = ClassicalMC<IsingModel, WolffCore>;
type SwIsing = ClassicalMC<IsingModel, SWCore>;
type HeatBathIsing = ClassicalMC<IsingModel, HeatBathCore>;
type ContinuousHeatBathXy = ClassicalMC<XYModel, ContinuousHeatBathCore>;
type MicrocanonicalHeisenberg = ClassicalMC<HeisenbergModel, MicrocanonicalCore>;
type HybridIsing = ClassicalMC<IsingModel, HybridCore<MetropolisCore, WolffCore>>;
type MetropolisPotts = ClassicalMC<PottsModel, MetropolisCore>;

fn assert_invalid_config(result: Result<(), CarloError>, label: &str) {
    match result {
        Err(CarloError::InvalidConfig { .. }) => {}
        Err(other) => panic!("{label}: expected InvalidConfig, got {other:?}"),
        Ok(()) => panic!("{label}: invalid configuration was silently accepted"),
    }
}

/// A valid base configuration with one parameter replaced by an invalid one.
fn with_param<T: FromParams>(base: &[(&str, &str)], field: &str, value: &str) {
    let mut params = Params::new();
    for (key, raw) in base {
        params.set(key, raw);
    }
    params.set(field, value);
    assert_invalid_config(T::validate_params(&params), field);
}

const ISING_BASE: [(&str, &str); 2] = [("L", "4"), ("beta", "0.5")];
const POTTS_BASE: [(&str, &str); 3] = [("L", "4"), ("beta", "0.5"), ("q", "3")];
const XY_BASE: [(&str, &str); 2] = [("L", "4"), ("beta", "0.5")];

// ── Lattice kernels sharing the ClassicalMC parameter contract ─────────

#[test]
fn lattice_kernels_reject_invalid_temperature() {
    for value in ["-0.5", "NaN", "not-a-number"] {
        with_param::<MetroIsing>(&ISING_BASE, "beta", value);
        with_param::<WolffIsing>(&ISING_BASE, "beta", value);
        with_param::<SwIsing>(&ISING_BASE, "beta", value);
        with_param::<HeatBathIsing>(&ISING_BASE, "beta", value);
        with_param::<ContinuousHeatBathXy>(&XY_BASE, "beta", value);
        with_param::<MicrocanonicalHeisenberg>(&XY_BASE, "beta", value);
        with_param::<HybridIsing>(&ISING_BASE, "beta", value);
        with_param::<MetropolisPotts>(&POTTS_BASE, "beta", value);
    }
}

#[test]
fn lattice_kernels_reject_unknown_lattices_and_bad_dimensions() {
    with_param::<MetroIsing>(&ISING_BASE, "lattice_type", "penrose-tiles");
    with_param::<WolffIsing>(&ISING_BASE, "lattice_type", "penrose-tiles");
    with_param::<SwIsing>(&ISING_BASE, "lattice_type", "penrose-tiles");
    with_param::<HeatBathIsing>(&ISING_BASE, "lattice_type", "penrose-tiles");
    with_param::<HybridIsing>(&ISING_BASE, "lattice_type", "penrose-tiles");
    with_param::<MetropolisPotts>(&POTTS_BASE, "lattice_type", "penrose-tiles");

    with_param::<MetroIsing>(&ISING_BASE, "L", "0");
    with_param::<MetroIsing>(&ISING_BASE, "L", "-3");
    with_param::<MetroIsing>(&ISING_BASE, "L", "four");
    // Open-boundary triangles are unsupported geometry, not silently wrong.
    let mut params = Params::new();
    params.set("lattice_type", "triangular");
    params.set("Lx", 4usize);
    params.set("Ly", 4usize);
    params.set("pbc", "false");
    assert_invalid_config(MetroIsing::validate_params(&params), "triangular OBC");
}

#[test]
fn lattice_kernels_reject_invalid_models_and_initial_states() {
    // Potts q < 2 and malformed q.
    with_param::<MetropolisPotts>(&POTTS_BASE, "q", "1");
    with_param::<MetropolisPotts>(&POTTS_BASE, "q", "zero");
    // Non-finite couplings.
    with_param::<MetroIsing>(&ISING_BASE, "J", "NaN");
    with_param::<MetropolisPotts>(&POTTS_BASE, "J", "inf");
    // Unknown initial states (from_params path; needs a real RNG).
    let mut params = Params::new();
    params.set("L", 4usize);
    params.set("beta", 0.5);
    params.set("initial_state", "lukewarm");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
    let result = MetroIsing::from_params(&params, &mut rng);
    assert!(
        matches!(result, Err(CarloError::InvalidConfig { .. })),
        "unknown initial_state must be rejected"
    );
}

#[test]
fn multi_spin_ising_rejects_invalid_parameters() {
    let base = [("L", "4"), ("beta", "0.5")];
    with_param::<MultiSpinIsing>(&base, "beta", "-1.0");
    with_param::<MultiSpinIsing>(&base, "J", "NaN");
    with_param::<MultiSpinIsing>(&base, "L", "0");
}

#[test]
fn wang_landau_params_reject_invalid_configurations() {
    let base = [("L", "8"), ("beta", "0.5")];
    with_param::<IsingWangLandau>(&base, "wl_minimum_visited_fraction", "1.5");
    with_param::<IsingWangLandau>(&base, "wl_minimum_visited_fraction", "-0.1");
    with_param::<IsingWangLandau>(&base, "wl_flatness", "1.5");
    with_param::<IsingWangLandau>(&base, "wl_final_log_f", "2.0");
    with_param::<IsingWangLandau>(&base, "wl_initial_log_f", "0.0");
    with_param::<IsingWangLandau>(&base, "wl_flatness_check_interval", "0");
    with_param::<IsingWangLandau>(&base, "beta", "-1.0");
    // The exact-axis reference implementation is limited to 24 sites.
    with_param::<IsingWangLandau>(&base, "L", "25");

    // Direct config validation rejects the same classes.
    let config = WangLandauConfig {
        minimum_visited_fraction: 2.0,
        ..Default::default()
    };
    assert!(config.validate().is_err());
    let config = WangLandauConfig {
        flatness: -0.5,
        ..Default::default()
    };
    assert!(config.validate().is_err());
}

#[test]
fn worm_params_reject_invalid_configurations() {
    let base = [("lattice_type", "square"), ("Lx", "4"), ("Ly", "4")];
    with_param::<IsingGraphWormMC>(&base, "beta", "-0.5");
    with_param::<IsingGraphWormMC>(&base, "J", "-1.0");
    with_param::<IsingGraphWormMC>(&base, "J", "NaN");
    with_param::<IsingGraphWormMC>(&base, "worm_close_probability", "0.0");
    with_param::<IsingGraphWormMC>(&base, "worm_close_probability", "1.0");
    with_param::<IsingGraphWormMC>(&base, "worm_fugacity", "0.0");
    with_param::<IsingGraphWormMC>(&base, "worm_fugacity", "-1.0");
    with_param::<IsingGraphWormMC>(&base, "lattice_type", "hyperbolic");
    // Both fugacity parameterizations at once are ambiguous.
    let mut params = Params::new();
    params.set("lattice_type", "square");
    params.set("Lx", 4usize);
    params.set("Ly", 4usize);
    params.set("worm_fugacity", "0.25");
    params.set("log_worm_fugacity", "0.0");
    assert_invalid_config(
        IsingGraphWormMC::validate_params(&params),
        "double fugacity",
    );
}

#[test]
fn kawasaki_and_bkl_params_reject_invalid_values_without_panicking() {
    // Before the 2026-08-19 fix, a non-finite J or negative beta reached
    // assert-backed constructors as a panic on user data.
    let base = [("L", "4"), ("beta", "0.5")];
    for value in ["-0.5", "NaN"] {
        with_param::<KawasakiIsingMC>(&base, "beta", value);
        with_param::<KineticIsingBklMC>(&base, "beta", value);
    }
    for value in ["NaN", "inf"] {
        with_param::<KawasakiIsingMC>(&base, "J", value);
        with_param::<KineticIsingBklMC>(&base, "J", value);
    }
    with_param::<KawasakiIsingMC>(&base, "up_fraction", "1.5");
    with_param::<KineticIsingBklMC>(&base, "up_fraction", "-0.1");
    with_param::<KineticIsingBklMC>(&base, "kinetic_rate", "arrhenius");
    with_param::<KineticIsingBklMC>(&base, "attempt_frequency", "0.0");
    with_param::<KineticIsingBklMC>(&base, "event_time_per_sweep", "0.0");
    with_param::<KawasakiIsingMC>(&base, "kawasaki_attempts_per_sweep", "-1");

    // from_params itself must reject (not panic on) user data.
    let mut params = Params::new();
    params.set("L", 4usize);
    params.set("J", "NaN");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(2);
    assert!(
        matches!(
            KawasakiIsingMC::from_params(&params, &mut rng),
            Err(CarloError::InvalidConfig { .. })
        ),
        "Kawasaki from_params must reject J = NaN with an error"
    );
    let mut params = Params::new();
    params.set("L", 4usize);
    params.set("beta", "-1.0");
    assert!(
        matches!(
            KineticIsingBklMC::from_params(&params, &mut rng),
            Err(CarloError::InvalidConfig { .. })
        ),
        "BKL from_params must reject negative beta with an error"
    );
}

#[test]
fn event_chain_params_reject_invalid_values() {
    // The event-chain adapter validates through from_params (grid geometry
    // and kernel construction are validated downstream, never asserted).
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
    let mut rejects = |parameters: &[(&str, &str)]| {
        let mut params = Params::new();
        params.set("n_particles", 4usize);
        params.set("box_length", "6.0");
        params.set("diameter", "1.0");
        for (key, value) in parameters {
            params.set(key, value);
        }
        matches!(
            HardSphereEventChainMC::<2>::from_params(&params, &mut rng),
            Err(CarloError::InvalidConfig { .. })
        )
    };
    assert!(rejects(&[("box_length", "0.0")]));
    assert!(rejects(&[("box_length", "-6.0")]));
    assert!(rejects(&[("n_particles", "0")]));
    assert!(rejects(&[("diameter", "0.0")]));
    assert!(rejects(&[("diameter", "-1.0")]));
    assert!(rejects(&[("chains_per_sweep", "0")]));
    assert!(rejects(&[("chain_length", "0.0")]));
}

#[test]
fn particle_ensembles_reject_invalid_values() {
    let nvt_base = [
        ("n_particles", "32"),
        ("density", "0.7"),
        ("beta", "1.0"),
        ("cutoff", "2.5"),
    ];
    with_param::<LennardJonesNvt<2>>(&nvt_base, "beta", "-1.0");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "density", "0.0");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "n_particles", "0");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "cutoff", "0.0");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "max_displacement", "-0.1");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "sigma", "0.0");
    with_param::<LennardJonesNvt<2>>(&nvt_base, "box_length", "1.0");

    let npt_base = [
        ("n_particles", "32"),
        ("density", "0.7"),
        ("beta", "1.0"),
        ("pressure", "1.0"),
    ];
    with_param::<LennardJonesNpt<2>>(&npt_base, "pressure", "NaN");
    with_param::<LennardJonesNpt<2>>(&npt_base, "max_log_volume_change", "0.0");
    with_param::<LennardJonesNpt<2>>(&npt_base, "target_acceptance", "1.5");

    let muvt_base = [("n_particles", "0"), ("density", "0.7"), ("beta", "1.0")];
    with_param::<LennardJonesMuVt<2>>(&muvt_base, "chemical_potential", "NaN");
    with_param::<LennardJonesMuVt<2>>(&muvt_base, "log_activity", "inf");
    with_param::<LennardJonesMuVt<2>>(&muvt_base, "beta", "-1.0");
}

#[test]
fn molecular_kernel_rejects_invalid_scales_and_topologies() {
    let topology = MoleculeTopology::new(2, vec![vec![0, 1]]).unwrap();

    // Move scales must be finite and positive (errors, not panics).
    assert!(MolecularMetropolisCore::<2>::new(topology.clone(), -0.1, 0.3).is_err());
    assert!(MolecularMetropolisCore::<2>::new(topology.clone(), 0.3, 0.0).is_err());
    assert!(MolecularMetropolisCore::<2>::new(topology.clone(), f64::NAN, 0.3).is_err());
    // Plane rotations need at least two dimensions.
    let chain_topology = MoleculeTopology::new(2, vec![vec![0], vec![1]]).unwrap();
    assert!(MolecularMetropolisCore::<1>::new(chain_topology, 0.3, 0.3).is_err());

    // Topology corruption is rejected loudly.
    assert!(MoleculeTopology::new(3, vec![vec![0, 5]]).is_err());
    assert!(MoleculeTopology::new(2, vec![vec![0], vec![0]]).is_err());
    assert!(MoleculeTopology::new(2, vec![Vec::new()]).is_err());
    // A topology/state particle-count mismatch is caught at sweep time by
    // validation, so mismatched construction must at least be detectable.
    let wrong = MoleculeTopology::new(3, vec![vec![0, 1]]).unwrap();
    assert!(wrong.validate_particle_count(2).is_err());

    // Valid scales and topology remain accepted.
    assert!(MolecularMetropolisCore::<2>::new(topology, 0.3, 0.3).is_ok());
}
