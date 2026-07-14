use super::common::{assert_close, direct_ising_energy, enumerate_ising};
use cmc_rs::{
    Bond, BondType, ClassicalMC, CsrLattice, Hamiltonian, HeatBathable, IsingModel,
    LocalFieldModel, MetropolisCore, ONModel, PottsModel, StandardStrategy, System,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn weighted_multigraph() -> CsrLattice {
    CsrLattice::from_edges(
        3,
        vec![
            Bond::new(0, 0, BondType::Generic, 0.37),
            Bond::new(0, 1, BondType::Generic, 1.25),
            Bond::new(0, 1, BondType::Generic, -0.4),
            Bond::new(1, 2, BondType::Generic, 0.8),
            Bond::new(2, 0, BondType::Generic, -1.1),
        ],
    )
}

#[test]
fn ising_energy_and_every_single_flip_delta_match_direct_edge_sum() {
    let lattice = weighted_multigraph();
    let model = IsingModel::new(0.73);
    for spins in enumerate_ising(lattice.n_sites) {
        let exact = direct_ising_energy(&spins, &lattice, model.j);
        assert_close(
            model.compute_total_energy(&spins, &lattice, 9.0),
            exact,
            1e-14,
        );
        for site in 0..lattice.n_sites {
            let mut proposed = spins.clone();
            proposed[site] = -proposed[site];
            let exact_delta = direct_ising_energy(&proposed, &lattice, model.j) - exact;
            assert_close(
                model.delta_energy(&spins, &lattice, site, &[proposed[site]]),
                exact_delta,
                2e-14,
            );
        }
    }
}

#[test]
fn potts_energy_and_every_replacement_delta_match_definition() {
    let lattice = weighted_multigraph();
    let model = PottsModel::new(0.61, 3);
    for code in 0..3usize.pow(lattice.n_sites as u32) {
        let mut rest = code;
        let mut spins = vec![0.0; lattice.n_sites];
        for spin in &mut spins {
            *spin = (rest % 3) as f64;
            rest /= 3;
        }
        let exact: f64 = lattice
            .edges
            .iter()
            .map(|edge| {
                if spins[edge.source] == spins[edge.target] {
                    -model.j * edge.weight
                } else {
                    0.0
                }
            })
            .sum();
        assert_close(
            model.compute_total_energy(&spins, &lattice, 4.0),
            exact,
            1e-14,
        );
        for site in 0..lattice.n_sites {
            for state in 0..3 {
                let mut proposed = spins.clone();
                proposed[site] = state as f64;
                let proposed_exact: f64 = lattice
                    .edges
                    .iter()
                    .map(|edge| {
                        if proposed[edge.source] == proposed[edge.target] {
                            -model.j * edge.weight
                        } else {
                            0.0
                        }
                    })
                    .sum();
                assert_close(
                    model.delta_energy(&spins, &lattice, site, &[state as f64]),
                    proposed_exact - exact,
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn self_loops_do_not_change_discrete_heat_bath_conditionals() {
    let base = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 0.7)]);
    let with_loop = CsrLattice::from_edges(
        2,
        vec![
            Bond::new(0, 1, BondType::Generic, 0.7),
            Bond::new(0, 0, BondType::Generic, 50.0),
        ],
    );

    let ising = IsingModel::new(1.2);
    let ising_spins = [-1.0, 1.0];
    let mut left = Xoshiro256PlusPlus::seed_from_u64(11);
    let mut right = left.clone();
    for _ in 0..1_000 {
        assert_eq!(
            ising.heat_bath_sample_site(&ising_spins, &base, 0, 0.83, &mut left),
            ising.heat_bath_sample_site(&ising_spins, &with_loop, 0, 0.83, &mut right)
        );
    }

    let potts = PottsModel::new(0.9, 4);
    let potts_spins = [3.0, 2.0];
    let mut left = Xoshiro256PlusPlus::seed_from_u64(12);
    let mut right = left.clone();
    for _ in 0..1_000 {
        assert_eq!(
            potts.heat_bath_sample_site(&potts_spins, &base, 0, 1.1, &mut left),
            potts.heat_bath_sample_site(&potts_spins, &with_loop, 0, 1.1, &mut right)
        );
    }
}

#[test]
fn self_loops_do_not_enter_on_local_field() {
    let base = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 0.7)]);
    let with_loop = CsrLattice::from_edges(
        2,
        vec![
            Bond::new(0, 1, BondType::Generic, 0.7),
            Bond::new(0, 0, BondType::Generic, 100.0),
        ],
    );
    let model = ONModel::<3>::new(1.3);
    let spins = [1.0, 0.0, 0.0, 0.0, 0.6, 0.8];
    let mut field_without = [0.0; 3];
    let mut field_with = [0.0; 3];
    model.local_field(&spins, &base, 0, &mut field_without);
    model.local_field(&spins, &with_loop, 0, &mut field_with);
    assert_eq!(field_without, field_with);
    assert_close(field_with[1], 1.3 * 0.7 * 0.6, 1e-15);
    assert_close(field_with[2], 1.3 * 0.7 * 0.8, 1e-15);
}

#[test]
fn built_in_models_reject_states_outside_their_physical_manifolds() {
    let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);

    let mut ising = System::new(lattice.clone(), 1, 1.0, 0.4);
    ising.spins[0] = 0.25;
    ising.energy = 0.0;
    assert!(ising.validate(&IsingModel::new(1.0)).is_err());

    let mut potts = System::new(lattice.clone(), 1, 0.0, 0.4);
    potts.spins[1] = 1.5;
    potts.energy = 0.0;
    assert!(potts.validate(&PottsModel::new(1.0, 3)).is_err());
    potts.spins[1] = 3.0;
    assert!(potts.validate(&PottsModel::new(1.0, 3)).is_err());

    let mut on = System::new(lattice, 3, 0.0, 0.4);
    on.spins.copy_from_slice(&[1.0, 0.0, 0.0, 0.0, 2.0, 0.0]);
    on.energy = 0.0;
    assert!(on.validate(&ONModel::<3>::new(1.0)).is_err());
}

#[test]
fn checkpoint_rejects_model_invalid_but_finite_spins() {
    let lattice = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice, 1, 1.0, 0.4);
    system.recompute_energy(&model);
    let mut mc: ClassicalMC<IsingModel, MetropolisCore<StandardStrategy>> =
        ClassicalMC::new(system, model, MetropolisCore::default());
    let mut snapshot = mc.save_snapshot();
    snapshot["spins"][0] = serde_json::json!(0.125);
    assert!(mc.load_snapshot(&snapshot).is_err());
}

#[derive(Debug, Clone, Copy)]
struct NonConstantSelfLoop;
impl cmc_rs::PairInteraction for NonConstantSelfLoop {
    fn spin_dim(&self) -> usize {
        1
    }
    fn coupling(&self) -> f64 {
        1.0
    }
    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64 {
        bond.weight * (left[0] + right[0])
    }
}

#[test]
fn self_loop_local_energy_is_counted_once_even_if_csr_incidences_are_reordered() {
    let mut lattice = CsrLattice::from_edges(
        2,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(0, 0, BondType::Generic, 2.0),
            Bond::new(0, 1, BondType::Generic, 3.0),
        ],
    );
    // Valid but noncanonical row ordering: edge 0 and edge 1 self-loop
    // incidences are interleaved. Physics must not depend on CSR row order.
    let row = lattice.offsets[0]..lattice.offsets[1];
    let pairs: Vec<_> = row
        .clone()
        .map(|i| (lattice.neighbors[i], lattice.edge_ids[i]))
        .collect();
    let desired = [pairs[0], pairs[2], pairs[1], pairs[3], pairs[4]];
    for (offset, &(neighbor, edge_id)) in desired.iter().enumerate() {
        let i = lattice.offsets[0] + offset;
        lattice.neighbors[i] = neighbor;
        lattice.edge_ids[i] = edge_id;
    }
    lattice.validate().unwrap();
    let model = NonConstantSelfLoop;
    let spins = [1.0, 4.0];
    let delta = model.delta_energy(&spins, &lattice, 0, &[2.0]);
    let mut proposed = spins;
    proposed[0] = 2.0;
    let exact = model.compute_total_energy(&proposed, &lattice, 1.0)
        - model.compute_total_energy(&spins, &lattice, 1.0);
    assert_close(delta, exact, 1e-14);
}
