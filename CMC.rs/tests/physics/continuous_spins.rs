use super::common::assert_close;
use cmc_rs::{
    build_chain, Algorithm, Hamiltonian, MicrocanonicalCore, ONModel, SimulationPhase, System,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn on_pair_energy_is_dot_product_on_unit_sphere() {
    let lattice = build_chain(2, false);
    let model = ONModel::<3>::new(1.7);
    let spins = [1.0, 0.0, 0.0, 0.0, 0.6, 0.8];
    assert_close(
        model.compute_total_energy(&spins, &lattice, 99.0),
        0.0,
        1e-15,
    );
    let spins = [0.0, 0.6, 0.8, 0.0, 0.6, 0.8];
    assert_close(
        model.compute_total_energy(&spins, &lattice, 99.0),
        -1.7,
        1e-15,
    );
}

#[test]
fn over_relaxation_preserves_energy_and_unit_norm_to_roundoff() {
    let lattice = build_chain(9, true);
    let model = ONModel::<3>::new(0.91);
    let mut system = System::new(lattice, 3, 0.0, 0.7);
    for (site, spin) in system.spins.chunks_exact_mut(3).enumerate() {
        let theta = 0.37 * site as f64;
        let z = 0.2 * ((site * 3 + 1) as f64).sin();
        let radial = (1.0 - z * z).sqrt();
        spin.copy_from_slice(&[radial * theta.cos(), radial * theta.sin(), z]);
    }
    system.recompute_energy(&model);
    let energy = system.energy;
    let mut kernel = MicrocanonicalCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xC01D);
    for _ in 0..200 {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }
    assert_close(system.energy, energy, 2e-10);
    assert_close(system.energy_error(&model), 0.0, 2e-10);
    for spin in system.spins.chunks_exact(3) {
        assert_close(spin.iter().map(|x| x * x).sum(), 1.0, 3e-12);
    }
}
