use cmc_rs::{CsrLattice, Hamiltonian, IsingModel};

/// Environment variable overriding the seed count of the multi-seed
/// z-score tests (nightly z-score monitoring).
pub const ZSCORE_SEEDS_ENV: &str = "SCUTTLE_ZSCORE_SEEDS";

/// Upper bound accepted for [`ZSCORE_SEEDS_ENV`] — guards against typos that
/// would silently turn a nightly run into a multi-day marathon.
const ZSCORE_SEEDS_MAX: usize = 4096;

/// Seed count for the multi-seed z-score tests.
///
/// Reads `SCUTTLE_ZSCORE_SEEDS` (see [`ZSCORE_SEEDS_ENV`]):
/// - unset or empty → `default` (the documented per-test seed count; this is
///   the regular-CI path and must stay byte-for-byte identical);
/// - an integer in `1..=4096` → that count (nightly high-power monitoring);
/// - anything else → panic with a clear message, because a silently degraded
///   monitoring run is worse than a loudly failed one.
pub fn zscore_seed_count(default: usize) -> usize {
    let raw = match std::env::var(ZSCORE_SEEDS_ENV) {
        Ok(raw) => raw,
        Err(std::env::VarError::NotPresent) => return default,
        Err(std::env::VarError::NotUnicode(raw)) => {
            panic!("{ZSCORE_SEEDS_ENV} is not valid Unicode: {raw:?}")
        }
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return default;
    }
    match trimmed.parse::<usize>() {
        Ok(n) if (1..=ZSCORE_SEEDS_MAX).contains(&n) => {
            eprintln!("[zscore] {ZSCORE_SEEDS_ENV}={n} (default {default})");
            n
        }
        _ => panic!(
            "{ZSCORE_SEEDS_ENV}={raw:?}: expected an integer in 1..={ZSCORE_SEEDS_MAX} \
             (unset the variable to use the default of {default} seeds)"
        ),
    }
}

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
