use qmc_rs::impurity::spin_boson::occupation::{
    CavityMode, OccupationSpinBosonModel, OccupationWorldlineSampler, SpinState,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn mixed_radix_basis_round_trips() {
    let model = OccupationSpinBosonModel::rabi(
        1.0,
        vec![
            CavityMode::new(0.8, 0.1, 3).unwrap(),
            CavityMode::new(1.2, 0.2, 4).unwrap(),
        ],
    )
    .unwrap();
    let basis = model.basis();
    for spin in [SpinState::Down, SpinState::Up] {
        for n0 in 0..3 {
            for n1 in 0..4 {
                let state = basis.encode(spin, &[n0, n1]).unwrap();
                assert_eq!(basis.spin(state), spin);
                assert_eq!(basis.occupation(state, 0), n0);
                assert_eq!(basis.occupation(state, 1), n1);
            }
        }
    }
}

#[test]
fn rabi_hamiltonian_has_correct_sqrt_n_matrix_elements() {
    let coupling = 0.3;
    let model =
        OccupationSpinBosonModel::rabi(1.1, vec![CavityMode::new(0.9, coupling, 4).unwrap()])
            .unwrap();
    let basis = model.basis();
    let h = model.hamiltonian();
    let down0 = basis.encode(SpinState::Down, &[0]).unwrap();
    let up1 = basis.encode(SpinState::Up, &[1]).unwrap();
    let down2 = basis.encode(SpinState::Down, &[2]).unwrap();
    let up3 = basis.encode(SpinState::Up, &[3]).unwrap();
    assert!((h[down0][up1] + coupling).abs() < 1e-14);
    assert!((h[down2][up3] + coupling * 3.0_f64.sqrt()).abs() < 1e-14);
}

#[test]
fn jaynes_cummings_conserves_total_excitation() {
    let model = OccupationSpinBosonModel::jaynes_cummings(
        1.0,
        vec![CavityMode::new(1.0, 0.25, 5).unwrap()],
    )
    .unwrap();
    let basis = model.basis();
    let h = model.hamiltonian();
    for (left, row) in h.iter().enumerate() {
        for (right, &element) in row.iter().enumerate() {
            if left == right || element.abs() < 1e-14 {
                continue;
            }
            let excitation = |state: usize| {
                basis.occupation(state, 0) + usize::from(basis.spin(state) == SpinState::Up)
            };
            assert_eq!(excitation(left), excitation(right));
        }
    }
}

#[test]
fn uncoupled_mode_matches_truncated_bose_distribution() {
    let beta: f64 = 2.0;
    let omega: f64 = 0.8;
    let cutoff = 8usize;
    let model =
        OccupationSpinBosonModel::rabi(1.3, vec![CavityMode::new(omega, 0.0, cutoff).unwrap()])
            .unwrap();
    let initial = model.basis().encode(SpinState::Down, &[0]).unwrap();
    let mut sampler = OccupationWorldlineSampler::new(model, beta, 12, initial).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0CC0_0001);
    for _ in 0..2_000 {
        sampler.sweep(&mut rng).unwrap();
    }
    let samples = 30_000usize;
    let mean = (0..samples)
        .map(|_| {
            sampler.sweep(&mut rng).unwrap();
            sampler.measure().unwrap().mode_occupations[0]
        })
        .sum::<f64>()
        / samples as f64;
    let weights = (0..cutoff)
        .map(|n| (-beta * omega * n as f64).exp())
        .collect::<Vec<_>>();
    let exact = weights
        .iter()
        .enumerate()
        .map(|(n, w)| n as f64 * w)
        .sum::<f64>()
        / weights.iter().sum::<f64>();
    assert!(
        (mean - exact).abs() < 0.025,
        "sampled={mean}, exact={exact}"
    );
}

#[test]
fn transfer_worldline_is_independent_of_slice_count() {
    let beta = 1.7;
    let model =
        OccupationSpinBosonModel::rabi(1.0, vec![CavityMode::new(0.9, 0.22, 5).unwrap()]).unwrap();
    let initial = model.basis().encode(SpinState::Down, &[0]).unwrap();
    let mut results = Vec::new();
    for (slices, seed) in [(4, 11), (11, 12)] {
        let mut sampler =
            OccupationWorldlineSampler::new(model.clone(), beta, slices, initial).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        for _ in 0..2_000 {
            sampler.sweep(&mut rng).unwrap();
        }
        let samples = 25_000usize;
        let energy = (0..samples)
            .map(|_| {
                sampler.sweep(&mut rng).unwrap();
                sampler.measure().unwrap().energy
            })
            .sum::<f64>()
            / samples as f64;
        results.push(energy);
    }
    assert!(
        (results[0] - results[1]).abs() < 0.04,
        "energies={results:?}"
    );
}

#[test]
fn smoke_explicit_solver_reports_requested_cavity_observables() {
    let model = OccupationSpinBosonModel::rabi(
        0.7,
        vec![
            CavityMode::new(1.0, 0.15, 4).unwrap(),
            CavityMode::new(1.4, 0.08, 3).unwrap(),
        ],
    )
    .unwrap();
    let initial = model.basis().encode(SpinState::Down, &[0, 0]).unwrap();
    let sampler = OccupationWorldlineSampler::new(model, 1.2, 6, initial).unwrap();
    let obs = sampler.measure().unwrap();
    assert_eq!(obs.mode_occupations.len(), 2);
    assert_eq!(obs.mode_number_squared.len(), 2);
    assert_eq!(obs.mode_factorial_moments.len(), 2);
    assert_eq!(obs.mode_g2_zero.len(), 2);
    assert_eq!(obs.mode_cross_correlations.len(), 4);
    assert_eq!(obs.spin_boson_covariance_z_n.len(), 2);
    assert!((0.5..=1.0).contains(&obs.reduced_spin_purity));
}
