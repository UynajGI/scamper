//! Integration tests for algorithm cores.

#[cfg(test)]
mod tests {
    use crate::algorithms::common::Algorithm;
    use crate::algorithms::metropolis::MetropolisCore;
    use crate::algorithms::microcanonical::MicrocanonicalCore;
    use crate::algorithms::swendsen_wang::SWCore;
    use crate::algorithms::wolff::WolffCore;
    use crate::core::cache::BatchEnergyPatch;
    use crate::core::r#move::BatchSpinMove;
    use crate::core::trial::TrialEvaluator;
    use crate::lattice::graph::{build_chain, Bond, BondType, CsrLattice};
    use crate::lattice::models::{IsingModel, PottsModel, XYModel};
    use crate::lattice::state::System;
    use rand::SeedableRng;

    fn rng() -> rand_xoshiro::Xoshiro256PlusPlus {
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42)
    }

    #[test]
    fn metropolis_cache_matches_exact_energy() {
        let model = IsingModel::new(1.0);
        let mut system = System::new(build_chain(8, true), 1, 1.0, 0.8);
        system.recompute_energy(&model);
        let mut algorithm = MetropolisCore::new();
        let mut random = rng();
        for _ in 0..100 {
            algorithm.sweep(&mut system, &model, &mut random);
        }
        assert!(system.energy_error(&model).abs() < 1e-10);
    }

    #[test]
    fn wolff_batch_cache_matches_exact_energy() {
        let model = XYModel::new(1.0);
        let mut system = System::new(build_chain(8, true), 2, 0.0, 1.0);
        for spin in system.spins.chunks_exact_mut(2) {
            spin[0] = 1.0;
        }
        system.recompute_energy(&model);
        let mut algorithm = WolffCore::new();
        let mut random = rng();
        for _ in 0..20 {
            algorithm.sweep(&mut system, &model, &mut random);
            assert!(system.energy_error(&model).abs() < 1e-10);
        }
    }

    #[test]
    fn sw_potts_assigns_valid_independent_states() {
        let model = PottsModel::new(1.0, 5);
        let mut system = System::new(build_chain(16, false), 1, 0.0, 0.0);
        system.recompute_energy(&model);
        let mut algorithm = SWCore::new();
        let mut random = rng();
        algorithm.sweep(&mut system, &model, &mut random);
        assert!(system.spins.iter().all(|spin| (0.0..5.0).contains(spin)));
        // beta=0 forms singleton clusters; independent assignments should
        // almost surely produce more than one state with this fixed seed.
        assert!(system.spins.windows(2).any(|pair| pair[0] != pair[1]));
        assert!(system.energy_error(&model).abs() < 1e-10);
    }

    #[test]
    fn microcanonical_preserves_and_tracks_energy() {
        let model = XYModel::new(1.0);
        let mut system = System::new(build_chain(4, true), 2, 0.0, 1.0);
        system
            .spins
            .copy_from_slice(&[1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0]);
        system.recompute_energy(&model);
        let before = system.energy;
        let mut algorithm = MicrocanonicalCore::new();
        algorithm.sweep(&mut system, &model, &mut rng());
        assert!((system.energy - before).abs() < 1e-10);
        assert!(system.energy_error(&model).abs() < 1e-10);
    }

    #[test]
    fn batch_delta_handles_parallel_edges_and_self_loops() {
        let lattice = CsrLattice::from_edges(
            2,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(0, 1, BondType::Generic, 0.5),
                Bond::new(1, 1, BondType::Generic, 0.25),
            ],
        );
        let model = IsingModel::new(1.0);
        let mut system = System::new(lattice, 1, 1.0, 1.0);
        system.recompute_energy(&model);
        let mut movement = BatchSpinMove::new(1);
        movement.push(0, &[-1.0]);
        movement.push(1, &[-1.0]);
        let mut patch = BatchEnergyPatch::default();
        system.evaluate_trial(&model, &movement, &mut patch);
        <System as TrialEvaluator<IsingModel, BatchSpinMove>>::commit_trial(
            &mut system,
            &movement,
            &patch,
        );
        assert!(system.energy_error(&model).abs() < 1e-12);
    }
}
