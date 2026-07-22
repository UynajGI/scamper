//! Project usage tests — verify that exported public APIs compose correctly
//! with the Carlo.rs scheduler and produce physically sane results.
//!
//! These tests are NOT exact physics validation (that's in tests/physics/).
//! They verify that the pieces plug together and produce reasonable output.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::*;

// ── HybridCore: Metropolis + Wolff alternating ────────────────────────────
//
// HybridCore does not implement Default (it has two type parameters), so it
// cannot be used through Scheduler::run_one. We test it at the algorithm level.

#[test]
fn hybrid_metropolis_wolff_alternates_correctly() {
    use cmc_rs::{Initializable, IsingModel};

    let lattice = build_square(8, 8, true);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 0.0, 0.44);

    // Hot start: random spins
    let mut init_rng = rand::rng();
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut init_rng);
        system.spin_at_mut(site, 1).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    let mut hybrid = HybridCore::new(MetropolisCore::new(), WolffCore::new());
    let mut rng = rand::rng();

    // Run 10 hybrid sweeps — each does 1 Metropolis + 1 Wolff
    for _ in 0..10 {
        hybrid.sweep(&mut system, &model, &mut rng);
    }

    // Energy cache should be consistent
    let e = system.energy_error(&model);
    assert!(e.abs() < 1e-10, "energy cache should be consistent");
}

#[test]
fn hybrid_with_repetitions_runs() {
    use cmc_rs::{Initializable, IsingModel};

    let lattice = build_square(4, 4, true);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 0.0, 1.0);

    // Hot start
    let mut init_rng = rand::rng();
    for site in 0..system.n_sites() {
        let spin = model.random_spin(&mut init_rng);
        system.spin_at_mut(site, 1).copy_from_slice(&spin);
    }
    system.recompute_energy(&model);

    // 3 Metropolis sweeps + 1 Wolff sweep per hybrid step
    let mut hybrid = HybridCore::new(MetropolisCore::new(), WolffCore::new()).repetitions(3, 1);
    let mut rng = rand::rng();

    hybrid.sweep(&mut system, &model, &mut rng);

    let e = system.energy_error(&model);
    assert!(e.abs() < 1e-10, "energy cache should be consistent");
}

// ── ContinuousHeatBathCore for O(N) models ────────────────────────────────

#[test]
fn continuous_heat_bath_xy_model_runs() {
    let mut params = Params::new();
    params.set("Lx", 6);
    params.set("Ly", 6);
    params.set("beta", 1.0);
    params.set("J", 1.0);

    type XYSim = ClassicalMC<XYModel, ContinuousHeatBathCore>;

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 99,
            ..Default::default()
        },
    )
    .run_one::<XYSim>(&params);

    let m = results.get("Magnetization").expect("Magnetization");
    assert!(m.mean >= 0.0 && m.mean <= 1.0, "magnetization in [0,1]");
}

#[test]
fn continuous_heat_bath_heisenberg_model_runs() {
    let mut params = Params::new();
    params.set("Lx", 4);
    params.set("Ly", 4);
    params.set("beta", 2.0);
    params.set("J", 1.0);

    type HeisSim = ClassicalMC<HeisenbergModel, ContinuousHeatBathCore>;

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 11,
            ..Default::default()
        },
    )
    .run_one::<HeisSim>(&params);

    assert!(results.get("Energy").is_some());
    let e = results.get("Energy").unwrap();
    assert!(e.mean < 0.0, "ferromagnetic energy should be negative");
}

// ── 3D lattice via build_hypercubic ───────────────────────────────────────

#[test]
fn three_d_cubic_ising_metropolis_runs() {
    let mut params = Params::new();
    params.set("Lx", 4);
    params.set("Ly", 4);
    params.set("Lz", 4);
    params.set("beta", 0.5);
    params.set("J", 1.0);
    params.set("lattice_type", "cubic");

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy");
    assert!(
        e.mean < 0.0,
        "3D ferromagnetic Ising energy should be negative"
    );
}

#[test]
fn three_d_cubic_ising_wolff_runs() {
    let mut params = Params::new();
    params.set("Lx", 4);
    params.set("Ly", 4);
    params.set("Lz", 4);
    params.set("beta", 0.8);
    params.set("J", 1.0);
    params.set("lattice_type", "cubic");

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 50,
            measurement_sweeps: 100,
            binsize: 25,
            base_seed: 33,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<IsingModel, WolffCore>>(&params);

    let e = results.get("Energy").expect("Energy observable missing");
    assert!(
        e.mean < 0.0,
        "3D ferromagnetic Ising at β=0.8 should have negative energy, got {}",
        e.mean
    );
    assert!(
        results.get("Magnetization").is_some(),
        "Magnetization observable missing"
    );
}

// ── Honeycomb lattice ─────────────────────────────────────────────────────

#[test]
fn honeycomb_ising_metropolis_runs() {
    let mut params = Params::new();
    params.set("Lx", 8);
    params.set("Ly", 8);
    params.set("beta", 1.0);
    params.set("J", 1.0);
    params.set("lattice_type", "honeycomb");

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 55,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let e = results.get("Energy").expect("Energy");
    assert!(
        e.mean < 0.0,
        "honeycomb ferromagnetic energy should be negative"
    );
}

// ── Scheduler-ready particle adapters ─────────────────────────────────────

#[test]
fn lennard_jones_nvt_runs_through_scheduler() {
    let mut params = Params::new();
    params.set("n_particles", 32);
    params.set("density", 0.5);
    params.set("beta", 1.0);
    params.set("cutoff", 2.5);
    params.set("max_displacement", 0.15);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<LennardJonesNvt<3>>(&params);

    let e = results.get("Energy").expect("NVT should produce Energy");
    assert!(
        e.mean.is_finite(),
        "NVT energy should be finite, got {}",
        e.mean
    );
}

#[test]
fn lennard_jones_npt_runs_through_scheduler() {
    let mut params = Params::new();
    params.set("n_particles", 16);
    params.set("density", 0.5);
    params.set("beta", 1.0);
    params.set("pressure", 0.5);
    params.set("cutoff", 2.5);
    params.set("max_displacement", 0.1);
    params.set("max_volume_scale", 0.05);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 50,
            measurement_sweeps: 100,
            binsize: 25,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<LennardJonesNpt<3>>(&params);

    let e = results.get("Energy").expect("NPT should produce Energy");
    assert!(
        e.mean.is_finite(),
        "NPT energy should be finite, got {}",
        e.mean
    );
}

#[test]
fn lennard_jones_muvt_runs_through_scheduler() {
    let mut params = Params::new();
    params.set("n_particles", 8);
    params.set("density", 0.3);
    params.set("beta", 2.0);
    params.set("cutoff", 2.5);
    params.set("max_displacement", 0.1);
    params.set("chemical_potential", -1.0);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 50,
            measurement_sweeps: 100,
            binsize: 25,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<LennardJonesMuVt<3>>(&params);

    let e = results.get("Energy").expect("μVT should produce Energy");
    assert!(
        e.mean.is_finite(),
        "μVT energy should be finite, got {}",
        e.mean
    );
}

// ── Scheduler-ready dynamics adapters ─────────────────────────────────────

#[test]
fn kawasaki_ising_runs_through_scheduler() {
    let mut params = Params::new();
    params.set("Lx", 8);
    params.set("Ly", 8);
    params.set("beta", 1.0);
    params.set("J", 1.0);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 77,
            ..Default::default()
        },
    )
    .run_one::<KawasakiIsingMC>(&params);

    let e = results
        .get("Energy")
        .expect("Kawasaki should produce Energy");
    assert!(
        e.mean < 0.0,
        "ferromagnetic Ising at β=1.0 should have negative energy, got {}",
        e.mean
    );
}

// ── MultiSpinIsing: 64-replica bit-packed Ising ───────────────────────────

#[test]
fn multi_spin_ising_runs_and_produces_observables() {
    let mut params = Params::new();
    params.set("Lx", 8);
    params.set("Ly", 8);
    params.set("beta", 0.44);
    params.set("J", 1.0);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 300,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<MultiSpinIsing>(&params);

    assert!(
        results.get("Energy").is_some(),
        "MultiSpin should produce Energy"
    );
    let e = results.get("Energy").unwrap();
    assert!(e.mean < 0.0, "ferromagnetic energy should be negative");
}

// ── Macrostate types construction ─────────────────────────────────────────

#[test]
fn macrostate_types_construct() {
    let _energy = EnergyMacrostate;
    let _mag = MagnetizationMacrostate;
    let _particle = ParticleNumberMacrostate;
}

// ── statistical_efficiency function ───────────────────────────────────────

#[test]
fn statistical_efficiency_on_uncorrelated_data() {
    // Alternating data has short autocorrelation time
    let data: Vec<f64> = (0..1000)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    let ess = statistical_efficiency(&data, 1.0);
    assert!(ess.effective_samples > 0.0);
    assert!(ess.integrated_autocorrelation_time >= 0.0);
}

// ── CanonicalReweighting::mean_observable ─────────────────────────────────

#[test]
fn reweight_mean_observable_matches_direct_computation() {
    let lattice = build_chain(4, true);
    let model = IsingModel::new(1.0);
    let dos = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = dos.axis().unwrap();
    let log_density = dos.log_density().unwrap();
    let beta = 0.5;

    let rw = canonical_reweight(&axis, &log_density, beta).unwrap();

    // Mean energy from reweighting should match direct enumeration
    use cmc_rs::Hamiltonian;
    let mut z = 0.0;
    let mut e_sum = 0.0;
    for mask in 0..1u32 << lattice.n_sites {
        let spins: Vec<f64> = (0..lattice.n_sites)
            .map(|s| if mask & (1 << s) == 0 { -1.0 } else { 1.0 })
            .collect();
        let energy = model.compute_total_energy(&spins, &lattice, 0.0);
        let w = (-beta * energy).exp();
        z += w;
        e_sum += w * energy;
    }
    let exact_e = e_sum / z;

    assert!(
        (rw.mean_energy() - exact_e).abs() < 1e-10,
        "reweighted E={} should match exact E={}",
        rw.mean_energy(),
        exact_e
    );
}
