use super::common::assert_close;
use cmc_rs::{HardSphereEventChain, KineticRateLaw, OrthorhombicCell, ParticleConfiguration};

#[test]
fn kinetic_ising_rates_obey_local_detailed_balance_over_extreme_deltas() {
    for beta in [0.0, 0.2, 1.0, 8.0] {
        for delta in [-50.0, -4.0, -0.25, 0.0, 0.25, 4.0, 50.0] {
            for law in [
                KineticRateLaw::glauber(1.7).unwrap(),
                KineticRateLaw::metropolis(1.7).unwrap(),
            ] {
                let forward = law.rate(beta, delta).unwrap();
                let reverse = law.rate(beta, -delta).unwrap();
                assert!(forward.is_finite() && forward >= 0.0);
                assert!(reverse.is_finite() && reverse >= 0.0);
                if forward > 0.0 && reverse > 0.0 {
                    assert_close((forward / reverse).ln(), -beta * delta, 2e-13);
                }
            }
        }
    }
}

#[test]
fn hard_sphere_exact_contacts_lift_at_zero_distance_instead_of_overlapping() {
    let cell = OrthorhombicCell::new([10.0, 4.0]).unwrap();
    let configuration = ParticleConfiguration::new(
        vec![[1.0, 1.0], [2.0, 1.0], [3.0, 1.0]],
        vec![0, 0, 0],
        cell,
    )
    .unwrap();
    let mut kernel = HardSphereEventChain::new(configuration, 1.0, 0.5, 1).unwrap();
    let outcome = kernel.step_with_lifting(0, 0, 1).unwrap();
    assert_eq!(outcome.collisions, 2);
    assert_eq!(outcome.final_particle, 2);
    assert_eq!(kernel.configuration().position(0), &[1.0, 1.0]);
    assert_eq!(kernel.configuration().position(1), &[2.0, 1.0]);
    assert_eq!(kernel.configuration().position(2), &[3.5, 1.0]);
    kernel.validate().unwrap();
}

#[derive(Clone)]
struct AbsorbingModel;
impl cmc_rs::RejectionFreeModel for AbsorbingModel {
    type State = ();
    type Patch = ();
    fn event_count(&self, _: &Self::State) -> usize {
        3
    }
    fn event_rate(&self, _: &Self::State, _: usize) -> Result<f64, cmc_rs::DynamicsError> {
        Ok(0.0)
    }
    fn prepare_event(
        &self,
        _: &Self::State,
        _: usize,
        _: &mut Self::Patch,
    ) -> Result<(), cmc_rs::DynamicsError> {
        Ok(())
    }
    fn commit_event(&self, _: &mut Self::State, _: usize, _: &Self::Patch) {}
    fn validate_state(&self, _: &Self::State) -> Result<(), cmc_rs::DynamicsError> {
        Ok(())
    }
}

#[test]
fn gillespie_absorbing_state_advances_observation_clock_without_fake_events() {
    use cmc_rs::GillespieKernel;
    use rand::SeedableRng;
    let mut kernel = GillespieKernel::new(AbsorbingModel, ()).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xAB50B);
    assert!(kernel.step(&mut rng).unwrap().is_none());
    assert_eq!(kernel.events(), 0);
    assert_eq!(kernel.advance_by(2.75, &mut rng).unwrap(), 0);
    assert_close(kernel.event_time(), 2.75, 1e-15);
}

#[test]
fn kawasaki_dynamics_conserves_signed_magnetization_exactly() {
    use cmc_rs::{build_chain, Algorithm, IsingModel, KawasakiCore, System};
    use rand::SeedableRng;
    let model = IsingModel::new(1.0);
    let mut state = System::new(build_chain(10, true), 1, 1.0, 0.7);
    state.spins[..4].fill(-1.0);
    state.recompute_energy(&model);
    let initial: f64 = state.spins.iter().sum();
    let mut kernel = KawasakiCore::new(20);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xCAAA5A);
    for _ in 0..1_000 {
        kernel.sweep(&mut state, &model, &mut rng);
    }
    assert_eq!(state.spins.iter().sum::<f64>(), initial);
    assert_close(state.energy_error(&model), 0.0, 2e-10);
}
