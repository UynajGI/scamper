use super::common::assert_close;
use cmc_rs::{
    build_chain, Algorithm, ContinuousHeatBathCore, Hamiltonian, MicrocanonicalCore, ONModel,
    SimulationPhase, System,
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
    for (site, spin) in system.spins.as_chunks_mut::<3>().0.iter_mut().enumerate() {
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
    for spin in system.spins.as_chunks::<3>().0 {
        assert_close(spin.iter().map(|x| x * x).sum(), 1.0, 3e-12);
    }
}

#[test]
fn continuous_heat_bath_infinite_t_is_uniform_on_sphere() {
    // At β→0 (infinite T), O(3) spins should be uniformly distributed on S².
    // For uniform-on-sphere: ⟨s_x⟩ = ⟨s_y⟩ = ⟨s_z⟩ = 0, ⟨|s|²⟩ = 1/3 per component.
    let lattice = build_chain(8, true);
    let model = ONModel::<3>::new(0.01); // very weak coupling → ~infinite T
    let mut system = System::new(lattice, 3, 0.0, 0.001); // β=0.001 ≈ infinite T

    let mut kernel = ContinuousHeatBathCore::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xBEEF);

    // Thermalize
    for _ in 0..500 {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }

    // Measure ⟨s_α⟩ for each component α and ⟨s_α²⟩
    let n_samples = 5000;
    let mut sx_sum = 0.0_f64;
    let mut sy_sum = 0.0_f64;
    let mut sz_sum = 0.0_f64;
    let mut sx2_sum = 0.0_f64;
    let mut sy2_sum = 0.0_f64;
    let mut sz2_sum = 0.0_f64;
    let mut count = 0_u64;

    for _ in 0..n_samples {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        for spin in system.spins.as_chunks::<3>().0 {
            sx_sum += spin[0];
            sy_sum += spin[1];
            sz_sum += spin[2];
            sx2_sum += spin[0] * spin[0];
            sy2_sum += spin[1] * spin[1];
            sz2_sum += spin[2] * spin[2];
            count += 1;
        }
    }

    let n = count as f64;
    let sx_mean = sx_sum / n;
    let sy_mean = sy_sum / n;
    let sz_mean = sz_sum / n;
    let sx2_mean = sx2_sum / n;
    let sy2_mean = sy2_sum / n;
    let sz2_mean = sz2_sum / n;

    // Uniform on S²: ⟨s_α⟩ ≈ 0 (within stochastic noise ~ 1/√n ≈ 0.005)
    assert!(sx_mean.abs() < 0.03, "⟨s_x⟩ = {sx_mean:.4}, expected ~0");
    assert!(sy_mean.abs() < 0.03, "⟨s_y⟩ = {sy_mean:.4}, expected ~0");
    assert!(sz_mean.abs() < 0.03, "⟨s_z⟩ = {sz_mean:.4}, expected ~0");

    // ⟨s_α²⟩ ≈ 1/3 for each component
    assert!(
        (sx2_mean - 1.0 / 3.0).abs() < 0.02,
        "⟨s_x²⟩ = {sx2_mean:.4}, expected ~1/3"
    );
    assert!(
        (sy2_mean - 1.0 / 3.0).abs() < 0.02,
        "⟨s_y²⟩ = {sy2_mean:.4}, expected ~1/3"
    );
    assert!(
        (sz2_mean - 1.0 / 3.0).abs() < 0.02,
        "⟨s_z²⟩ = {sz2_mean:.4}, expected ~1/3"
    );
}
