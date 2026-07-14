use cmc_rs::{CsrLattice, Hamiltonian, IsingModel};

pub fn assert_close(left: f64, right: f64, tolerance: f64) {
    assert!(
        (left - right).abs() <= tolerance,
        "{left:.17e} != {right:.17e}; |Δ|={:.3e}, tolerance={tolerance:.3e}",
        (left - right).abs()
    );
}

pub fn enumerate_ising(n_sites: usize) -> Vec<Vec<f64>> {
    (0..1usize << n_sites)
        .map(|mask| {
            (0..n_sites)
                .map(|site| if mask & (1 << site) == 0 { -1.0 } else { 1.0 })
                .collect()
        })
        .collect()
}

pub fn direct_ising_energy(spins: &[f64], lattice: &CsrLattice, coupling: f64) -> f64 {
    lattice
        .edges
        .iter()
        .map(|edge| -coupling * edge.weight * spins[edge.source] * spins[edge.target])
        .sum()
}

pub fn exact_ising_moments(lattice: &CsrLattice, coupling: f64, beta: f64) -> (f64, f64, f64, f64) {
    let model = IsingModel::new(coupling);
    let mut z = 0.0;
    let mut e = 0.0;
    let mut e2 = 0.0;
    let mut m2 = 0.0;
    for spins in enumerate_ising(lattice.n_sites) {
        let energy = model.compute_total_energy(&spins, lattice, 1.0);
        let magnetization = spins.iter().sum::<f64>() / lattice.n_sites as f64;
        let weight = (-beta * energy).exp();
        z += weight;
        e += weight * energy;
        e2 += weight * energy * energy;
        m2 += weight * magnetization * magnetization;
    }
    (z, e / z, e2 / z, m2 / z)
}
