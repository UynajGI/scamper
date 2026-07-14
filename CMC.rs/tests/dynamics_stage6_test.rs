use carlo_rs::{Context, FromParams, MonteCarlo, Params};
use cmc_rs::{
    build_chain, Algorithm, BklIsingKernel, DynamicsError, GillespieKernel, HardSphereEventChain,
    HardSphereEventChainMC, IsingModel, KawasakiCore, KineticIsingBklMC, KineticIsingModel,
    KineticRateLaw, OrthorhombicCell, ParticleConfiguration, RejectionFreeModel, System,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Clone)]
struct TwoEventModel;

impl RejectionFreeModel for TwoEventModel {
    type State = [u64; 2];
    type Patch = ();

    fn event_count(&self, _state: &Self::State) -> usize {
        2
    }

    fn event_rate(&self, _state: &Self::State, event: usize) -> Result<f64, DynamicsError> {
        Ok([1.0, 3.0][event])
    }

    fn prepare_event(
        &self,
        _state: &Self::State,
        _event: usize,
        _patch: &mut Self::Patch,
    ) -> Result<(), DynamicsError> {
        Ok(())
    }

    fn commit_event(&self, state: &mut Self::State, event: usize, _patch: &Self::Patch) {
        state[event] += 1;
    }

    fn validate_state(&self, _state: &Self::State) -> Result<(), DynamicsError> {
        Ok(())
    }
}

#[test]
fn direct_gillespie_selects_rates_and_exponential_waits() {
    let mut kernel = GillespieKernel::new(TwoEventModel, [0, 0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6100);
    let mut waiting_sum = 0.0;
    for _ in 0..100_000 {
        waiting_sum += kernel.step(&mut rng).unwrap().unwrap().delta_time;
    }
    let counts = kernel.state();
    let fraction_second = counts[1] as f64 / (counts[0] + counts[1]) as f64;
    assert!((fraction_second - 0.75).abs() < 0.01);
    assert!((waiting_sum / 100_000.0 - 0.25).abs() < 0.005);
}

#[test]
fn fixed_event_time_windows_do_not_overshoot_measurement_clock() {
    let mut kernel = GillespieKernel::new(TwoEventModel, [0, 0]).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6101);
    for window in 1..=100 {
        kernel.advance_by(0.125, &mut rng).unwrap();
        assert!((kernel.event_time() - window as f64 * 0.125).abs() < 1e-13);
    }
}

#[test]
fn kawasaki_exchange_conserves_signed_magnetization_and_energy_cache() {
    let lattice = build_chain(8, true);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 1.0, 0.6);
    system.spins[..4].fill(-1.0);
    system.recompute_energy(&model);
    let initial_magnetization: f64 = system.spins.iter().sum();
    let mut kernel = KawasakiCore::new(16).with_cache_audit_interval(1);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6102);
    for _ in 0..5_000 {
        kernel.sweep(&mut system, &model, &mut rng);
    }
    assert_eq!(system.spins.iter().sum::<f64>(), initial_magnetization);
    assert!(system.energy_error(&model).abs() < 1e-10);
    assert!(kernel.accepts() > 0);
}

#[test]
fn bkl_rate_cache_and_snapshot_preserve_exact_future_trajectory() {
    let lattice = build_chain(8, true);
    let mut state = System::new(lattice, 1, 1.0, 0.45);
    for (site, spin) in state.spins.iter_mut().enumerate() {
        *spin = if site % 3 == 0 { -1.0 } else { 1.0 };
    }
    let model = KineticIsingModel::new(1.0, KineticRateLaw::glauber(1.0).unwrap()).unwrap();
    let mut original = BklIsingKernel::new(model.clone(), state.clone(), 7).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6103);
    for _ in 0..500 {
        original.step(&mut rng).unwrap();
    }
    original.validate().unwrap();
    let snapshot = original.save_snapshot();
    let mut restored = BklIsingKernel::new(model, state, 7).unwrap();
    restored.load_snapshot(&snapshot).unwrap();
    let mut restored_rng = rng.clone();
    for _ in 0..2_000 {
        original.step(&mut rng).unwrap();
        restored.step(&mut restored_rng).unwrap();
    }
    assert_eq!(original.save_snapshot(), restored.save_snapshot());
}

#[test]
fn bkl_fixed_time_sampling_matches_exact_small_ising_energy() {
    let beta = 0.4;
    let lattice = build_chain(4, true);
    let state = System::new(lattice, 1, 1.0, beta);
    let model = KineticIsingModel::new(1.0, KineticRateLaw::glauber(1.0).unwrap()).unwrap();
    let mut kernel = BklIsingKernel::new(model, state, 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x6104);
    for _ in 0..5_000 {
        kernel.advance_by(0.2, &mut rng).unwrap();
    }
    let mut energy_sum = 0.0;
    let samples = 100_000;
    for _ in 0..samples {
        kernel.advance_by(0.2, &mut rng).unwrap();
        energy_sum += kernel.state().energy;
    }
    let measured = energy_sum / samples as f64;
    let exact = exact_ring_energy(beta, 4);
    assert!((measured - exact).abs() < 0.06, "{measured} vs {exact}");
}

#[test]
fn hard_sphere_event_chain_transfers_lifting_at_exact_collision() {
    let cell = OrthorhombicCell::new([10.0, 5.0]).unwrap();
    let configuration =
        ParticleConfiguration::new(vec![[1.0, 1.0], [3.0, 1.0]], vec![0, 0], cell).unwrap();
    let mut kernel = HardSphereEventChain::new(configuration, 1.0, 3.0, 1).unwrap();
    let outcome = kernel.step_with_lifting(0, 0, 1).unwrap();
    assert_eq!(outcome.collisions, 1);
    assert_eq!(outcome.final_particle, 1);
    assert!((kernel.configuration().position(0)[0] - 2.0).abs() < 1e-12);
    assert!((kernel.configuration().position(1)[0] - 5.0).abs() < 1e-12);
    kernel.validate().unwrap();
}

#[test]
fn hard_sphere_event_chain_wraps_and_snapshot_restores() {
    let cell = OrthorhombicCell::new([8.0, 8.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[1.0, 1.0], [3.0, 3.0], [5.0, 5.0]],
        vec![0, 0, 0],
        cell,
    )
    .unwrap();
    let mut original = HardSphereEventChain::new(configuration.clone(), 1.0, 11.0, 3).unwrap();
    original.step_with_lifting(0, 0, 1).unwrap();
    let snapshot = original.save_snapshot();
    let mut restored = HardSphereEventChain::new(configuration, 1.0, 11.0, 3).unwrap();
    restored.load_snapshot(&snapshot).unwrap();
    for chain in 0..50 {
        let particle = chain % 3;
        let axis = chain % 2;
        let direction = if chain % 2 == 0 { 1 } else { -1 };
        original
            .step_with_lifting(particle, axis, direction)
            .unwrap();
        restored
            .step_with_lifting(particle, axis, direction)
            .unwrap();
    }
    assert_eq!(original.save_snapshot(), restored.save_snapshot());
}

#[test]
fn scheduler_adapters_advance_distinct_clocks() {
    let mut kinetic_params = Params::new();
    kinetic_params.set("lattice_type", "chain");
    kinetic_params.set("L", 8);
    kinetic_params.set("beta", 0.4);
    kinetic_params.set("event_time_per_sweep", 0.25);
    let rng = Xoshiro256PlusPlus::seed_from_u64(0x6105);
    let mut context = Context::new(rng, 0);
    let mut kinetic = KineticIsingBklMC::from_params(&kinetic_params, &mut context.rng).unwrap();
    kinetic.sweep(&mut context);
    assert!((context.event_time() - 0.25).abs() < 1e-13);
    assert_eq!(context.attempted_updates(), context.accepted_moves());

    let mut event_chain_params = Params::new();
    event_chain_params.set("n_particles", 4);
    event_chain_params.set("box_length", 6.0);
    event_chain_params.set("diameter", 1.0);
    event_chain_params.set("chain_length", 2.0);
    event_chain_params.set("chains_per_sweep", 3);
    let rng = Xoshiro256PlusPlus::seed_from_u64(0x6106);
    let mut context = Context::new(rng, 0);
    let mut event_chain =
        HardSphereEventChainMC::<2>::from_params(&event_chain_params, &mut context.rng).unwrap();
    event_chain.sweep(&mut context);
    assert_eq!(context.attempted_updates(), 3);
    assert_eq!(context.accepted_moves(), 3);
    assert_eq!(context.event_time(), 0.0);
}

fn exact_ring_energy(beta: f64, sites: usize) -> f64 {
    let mut partition = 0.0;
    let mut energy_sum = 0.0;
    for mask in 0..1usize << sites {
        let mut energy = 0.0;
        for site in 0..sites {
            let left = if mask & (1 << site) == 0 { -1.0 } else { 1.0 };
            let right_site = (site + 1) % sites;
            let right = if mask & (1 << right_site) == 0 {
                -1.0
            } else {
                1.0
            };
            energy -= left * right;
        }
        let weight = (-beta * energy).exp();
        partition += weight;
        energy_sum += weight * energy;
    }
    energy_sum / partition
}
