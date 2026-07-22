//! Finite-temperature Rabi model QMC: ⟨n⟩ FSS study.
//!
//! Uses the OccupationWorldlineSampler from QMC.rs.
//! For each (η, β, λ), run MC and measure ⟨n⟩.
//! Plot ñ = η·⟨n⟩ vs λ for different η — look for crossings.
//!
//! Output: results/qmc_fss/beta_{β}.csv

use qmc_rs::impurity::spin_boson::occupation::{
    model::{CavityMode, OccupationSpinBosonModel},
    OccupationWorldlineSampler,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fs;
use std::io::Write;
use std::time::Instant;

const DELTA: f64 = 1.0;
const OMEGAS: [f64; 3] = [50.0, 200.0, 1000.0];
const BETAS: [f64; 4] = [16.0, 64.0, 256.0, 1024.0];
const CUTOFF: usize = 200;
const SLICES: usize = 8;

/// Two-region scan: dense near λ_c=0.5 and dense near observed crossing λ≈5-9
const LAMBDA_VALUES: [f64; 25] = [
    0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0,
    4.0, 4.5, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5,
    7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0,
];

const THERMAL_SWEEPS: u64 = 30_000;
const MEASURE_SWEEPS: u64 = 120_000;

fn run_mc(omega: f64, g: f64, beta: f64, seed: u64) -> (f64, f64, f64, f64) {
    let model = OccupationSpinBosonModel::rabi(
        DELTA,
        vec![CavityMode::new(omega, g, CUTOFF).unwrap()],
    )
    .unwrap();

    let mut sampler = OccupationWorldlineSampler::new(model, beta, SLICES, 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut n_sum = 0.0_f64;
    let mut n2_sum = 0.0_f64;
    let mut sz_sum = 0.0_f64;
    let mut e_sum = 0.0_f64;
    let mut samples = 0_u64;

    for sweep in 0..(THERMAL_SWEEPS + MEASURE_SWEEPS) {
        sampler.sweep(&mut rng).unwrap();
        if sweep >= THERMAL_SWEEPS {
            let obs = sampler.measure().unwrap();
            n_sum += obs.total_boson_number;
            n2_sum += obs.total_boson_number * obs.total_boson_number;
            sz_sum += obs.sigma_z;
            e_sum += obs.energy;
            samples += 1;
        }
    }

    let n = samples as f64;
    (
        n_sum / n,
        (n2_sum / n - (n_sum / n).powi(2)).sqrt(),
        sz_sum / n,
        e_sum / n,
    )
}

fn main() {
    let out_dir = "results/qmc_fss";
    fs::create_dir_all(out_dir).unwrap();

    let start = Instant::now();
    eprintln!("=== Rabi QMC FSS: ñ = η·⟨n⟩ ===");
    eprintln!("Delta={DELTA}, Omegas={OMEGAS:?}, Betas={BETAS:?}");
    eprintln!("Cutoff={CUTOFF}, Slices={SLICES}");
    eprintln!("λ points: {LAMBDA_VALUES:?}");
    eprintln!("Thermal={THERMAL_SWEEPS}, Measure={MEASURE_SWEEPS}\n");

    let lambda_grid: Vec<f64> = LAMBDA_VALUES.to_vec();
    let n_lambda = lambda_grid.len();

    let mut log = fs::File::create(format!("{out_dir}/run.log")).unwrap();

    for &beta in &BETAS {
        eprintln!("--- beta = {beta} ---");
        let t0 = Instant::now();

        let mut csv = String::from("lambda");
        for &omega in &OMEGAS {
            csv.push_str(&format!(",n_omega{omega:.0},n_err_omega{omega:.0},ntilde_omega{omega:.0}"));
        }
        csv.push('\n');

        for (li, &lambda) in lambda_grid.iter().enumerate() {
            csv.push_str(&format!("{lambda:.6}"));

            for (oi, &omega) in OMEGAS.iter().enumerate() {
                let eta = omega / DELTA;
                let g = lambda * (omega * DELTA).sqrt();
                let seed = (li as u64) * 1000 + (oi as u64) * 100 + beta as u64;

                let (n_mean, n_err, _sz, _e) = run_mc(omega, g, beta, seed);
                let n_tilde = eta * n_mean;

                csv.push_str(&format!(",{n_mean:.10},{n_err:.10},{n_tilde:.10}"));
            }
            csv.push('\n');

            if li % 5 == 0 {
                eprintln!(
                    "  λ={lambda:.3} ({li}/{n_lambda}, {:.0}s)",
                    t0.elapsed().as_secs_f64()
                );
            }
        }

        let csv_path = format!("{out_dir}/beta_{beta:.0}.csv");
        fs::write(&csv_path, &csv).unwrap();
        writeln!(log, "beta={beta}: saved {csv_path} ({:.0}s)", t0.elapsed().as_secs_f64()).unwrap();

        eprintln!("  beta={beta} done ({:.0}s)\n", t0.elapsed().as_secs_f64());
    }

    eprintln!("Total: {:.0}s", start.elapsed().as_secs_f64());
}
