//! Finite-temperature Rabi model QMC via wormhole solver.
//!
//! The wormhole solver uses retarded interaction (integrates out bosons),
//! so there is NO boson cutoff limitation. This allows simulation at
//! arbitrary coupling λ and η.
//!
//! The sampled σz in the rotated basis corresponds to physical σx (the
//! Z₂-odd Rabi QPT order parameter). We measure:
//!   ⟨σx⟩, ⟨σx²⟩, ⟨σx⁴⟩ → U4 = 1 - ⟨σx⁴⟩/(3⟨σx²⟩²)
//!
//! Output: results/qmc_wormhole/beta_{β}.csv

use qmc_rs::algorithm::UpdateSchedule;
use qmc_rs::impurity::spin_boson::bath::{Bath, SingleModeBath};
use qmc_rs::impurity::spin_boson::model::ImpurityModel;
use qmc_rs::impurity::spin_boson::observables::measure_observables;
use qmc_rs::impurity::spin_boson::wormhole::{
    configuration::WormholeConfiguration,
    updates::WormholeEngine,
};
use qmc_rs::QmcKernel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fs;
use std::io::Write;
use std::time::Instant;

const DELTA: f64 = 1.0;
const OMEGAS: [f64; 3] = [50.0, 200.0, 1000.0];
const BETAS: [f64; 1] = [1024.0];

const THERMAL_SWEEPS: u64 = 15_000;
const MEASURE_SWEEPS: u64 = 40_000;

fn run_wormhole_mc(
    omega: f64,
    g_coupling: f64,
    beta: f64,
    seed: u64,
) -> (f64, f64, f64, f64, f64, usize) {
    // Build the Rabi model via rotated impurity.
    // effective coupling lambda_eff = g²/ω (this is the retarded interaction weight).
    let lambda_eff = g_coupling * g_coupling / omega;
    let bath = Bath::SingleMode(SingleModeBath::new(omega).unwrap());
    // tunnelling = Δ (the spin splitting in the rotated basis)
    let model = ImpurityModel::rotated_impurity(bath, lambda_eff, DELTA, None).unwrap();
    let model_clone = model.clone();

    let schedule = UpdateSchedule::new(1, 1, 16);
    let mut config = WormholeConfiguration::new(beta, 1).unwrap();
    let mut engine = WormholeEngine::new(model, schedule);
    engine.set_transverse_bins(64);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut sx_sum = 0.0_f64;
    let mut sx2_sum = 0.0_f64;
    let mut sx4_sum = 0.0_f64;
    let mut order_sum = 0.0_f64;
    let mut samples = 0_u64;

    for sweep in 0..(THERMAL_SWEEPS + MEASURE_SWEEPS) {
        engine.sweep(&mut config, &mut rng).unwrap();
        if sweep >= THERMAL_SWEEPS {
            let obs = measure_observables(&config, &model_clone, 16, &mut rng).unwrap();
            // magnetization_sigma_z is the sampled σz = physical σx (rotated basis)
            let sx = obs.magnetization_sigma_z;
            sx_sum += sx;
            sx2_sum += sx * sx;
            sx4_sum += sx.powi(4);
            order_sum += obs.expansion_order;
            samples += 1;
        }
    }

    let n = samples as f64;
    (
        sx_sum / n,
        sx2_sum / n,
        sx4_sum / n,
        order_sum / n,
        (sx2_sum / n - (sx_sum / n).powi(2)).sqrt(),
        samples as usize,
    )
}

fn main() {
    let out_dir = "results/qmc_wormhole";
    fs::create_dir_all(out_dir).unwrap();

    // Uniform grid: Δλ=0.05 in [0.05, 3.0] → 60 points
    let lam_min = 0.05; let lam_max = 3.0; let n_pts = 60;
    let lambdas: Vec<f64> = (0..n_pts).map(|i| lam_min + (lam_max - lam_min) * i as f64 / (n_pts - 1) as f64).collect();

    let start = Instant::now();
    eprintln!("=== Rabi QMC FSS via Wormhole: U4(σx) ===");
    eprintln!("Delta={DELTA}, Omegas={OMEGAS:?}, Betas={BETAS:?}");
    eprintln!("λ: {n_pts} pts, Δλ={:.3}", (lam_max - lam_min) / (n_pts - 1) as f64);
    eprintln!("Thermal={THERMAL_SWEEPS}, Measure={MEASURE_SWEEPS}\n");

    let mut log = fs::File::create(format!("{out_dir}/run.log")).unwrap();

    for &beta in &BETAS {
        eprintln!("--- beta = {beta} ---");
        let t0 = Instant::now();

        let mut csv = String::from("lambda");
        for &omega in &OMEGAS {
            csv.push_str(&format!(",sx_omega{omega:.0},sx2_omega{omega:.0},sx4_omega{omega:.0},u4_omega{omega:.0},order_omega{omega:.0}"));
        }
        csv.push('\n');

        for (li, &lambda) in lambdas.iter().enumerate() {
            csv.push_str(&format!("{lambda:.6}"));

            for (oi, &omega) in OMEGAS.iter().enumerate() {
                let g = lambda * (omega * DELTA).sqrt();
                let seed = (li as u64) * 1000 + (oi as u64) * 100 + beta as u64;

                let (sx, sx2, sx4, order, _sx_err, _n) =
                    run_wormhole_mc(omega, g, beta, seed);

                let u4 = if sx2.abs() > 1e-14 {
                    1.0 - sx4 / (3.0 * sx2 * sx2)
                } else {
                    0.0
                };

                csv.push_str(&format!(",{sx:.10},{sx2:.10},{sx4:.10},{u4:.10},{order:.6}"));
            }
            csv.push('\n');

            if li % 10 == 0 {
                eprintln!(
                    "  λ={lambda:.3} ({li}/{n_pts}, {:.0}s)",
                    t0.elapsed().as_secs_f64()
                );
            }
        }

        let csv_path = format!("{out_dir}/beta_{beta:.0}.csv");
        fs::write(&csv_path, &csv).unwrap();
        writeln!(
            log,
            "beta={beta}: saved {csv_path} ({:.0}s)",
            t0.elapsed().as_secs_f64()
        )
        .unwrap();

        eprintln!("  beta={beta} done ({:.0}s)\n", t0.elapsed().as_secs_f64());
    }

    eprintln!("Total: {:.0}s", start.elapsed().as_secs_f64());
}
