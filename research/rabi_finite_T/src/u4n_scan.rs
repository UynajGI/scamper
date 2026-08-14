//! U4(n) Binder cumulant of photon number.
//! U4(n) = 1 - κ₄/(3σ⁴) where κ₄ is the fourth central moment.
//! λ spacing = 0.005, T=0 (β=1024).

use qmc_rs::impurity::spin_boson::occupation::{
    model::{CavityMode, OccupationSpinBosonModel},
    OccupationWorldlineSampler,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fs;
use std::time::Instant;

const DELTA: f64 = 1.0;
const OMEGAS: [f64; 4] = [10.0, 50.0, 200.0, 500.0];
const BETA: f64 = 1024.0;
const CUTOFF: usize = 80;
const SLICES: usize = 8;

const THERMAL_SWEEPS: u64 = 50_000;
const MEASURE_SWEEPS: u64 = 100_000;

fn run(omega: f64, g: f64, seed: u64) -> (f64, f64, f64, f64, f64, f64) {
    let model = OccupationSpinBosonModel::rabi(
        DELTA,
        vec![CavityMode::new(omega, g, CUTOFF).unwrap()],
    ).unwrap();
    let mut sampler = OccupationWorldlineSampler::new(model, BETA, SLICES, 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let (mut s1, mut s2, mut s3, mut s4, mut n) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0_u64);
    for sweep in 0..(THERMAL_SWEEPS + MEASURE_SWEEPS) {
        sampler.sweep(&mut rng).unwrap();
        if sweep >= THERMAL_SWEEPS {
            let val = sampler.measure().unwrap().total_boson_number;
            s1 += val; s2 += val*val; s3 += val*val*val; s4 += val*val*val*val; n += 1;
        }
    }
    let f = n as f64;
    let m1 = s1/f; let m2 = s2/f; let m3 = s3/f; let m4 = s4/f;
    let var = m2 - m1*m1;
    let c4 = m4 - 4.0*m3*m1 + 6.0*m2*m1*m1 - 3.0*m1.powi(4);
    let u4 = if var > 1e-30 { 1.0 - c4/(3.0*var*var) } else { 0.0 };
    (m1, m2, m3, m4, var, u4)
}

fn ed_u4n(omega: f64, g: f64) -> f64 {
    // Same ED as before, compute moments of n from ground state
    let dim = 2 * CUTOFF;
    let mut h = vec![vec![0.0_f64; dim]; dim];
    for n0 in 0..CUTOFF {
        h[2*n0][2*n0] = omega*n0 as f64 - 0.5*DELTA;
        h[2*n0+1][2*n0+1] = omega*n0 as f64 + 0.5*DELTA;
        if n0+1 < CUTOFF {
            let a = g*((n0+1) as f64).sqrt();
            h[2*n0+1][2*(n0+1)]+=a; h[2*(n0+1)][2*n0+1]+=a;
            h[2*n0][2*(n0+1)+1]+=a; h[2*(n0+1)+1][2*n0]+=a;
        }
    }
    let shift = h.iter().flat_map(|r|r.iter()).map(|x|x.abs()).fold(0.0,f64::max)+1.0;
    let mut b = vec![vec![0.0;dim];dim];
    for i in 0..dim { for j in 0..dim { b[i][j]=-h[i][j]; } b[i][i]+=shift; }
    let mut v = vec![0.0_f64;dim]; v[0]=1.0;
    for _ in 0..1000 {
        let mut bv = vec![0.0;dim];
        for i in 0..dim { for j in 0..dim { bv[i]+=b[i][j]*v[j]; } }
        let nr = bv.iter().map(|x|x*x).sum::<f64>().sqrt();
        if nr<1e-30 {break}
        for i in 0..dim { v[i]=bv[i]/nr; }
    }
    let (mut m1, mut m2, mut m3, mut m4)=(0.0_f64,0.0_f64,0.0_f64,0.0_f64);
    for n0 in 0..CUTOFF {
        let p = v[2*n0]*v[2*n0] + v[2*n0+1]*v[2*n0+1];
        let nv = n0 as f64;
        m1+=nv*p; m2+=nv*nv*p; m3+=nv*nv*nv*p; m4+=nv*nv*nv*nv*p;
    }
    let var = m2 - m1*m1;
    let c4 = m4 - 4.0*m3*m1 + 6.0*m2*m1*m1 - 3.0*m1.powi(4);
    if var>1e-30 { 1.0 - c4/(3.0*var*var) } else { 0.0 }
}

fn main() {
    let out_dir = "results/qmc_verify";
    fs::create_dir_all(out_dir).unwrap();
    let lam_min = 0.200; let lam_max = 0.800; let n_pts = 121; // Δλ = 0.005
    let lam: Vec<f64> = (0..n_pts).map(|i| lam_min + (lam_max-lam_min)*i as f64/(n_pts-1) as f64).collect();
    eprintln!("=== U4(n) scan: λ∈[{lam_min},{lam_max}] Δλ=0.005 ===");
    let start=Instant::now();
    let mut csv=String::from("lambda");
    for &omega in &OMEGAS {
        let eta=omega/DELTA;
        csv.push_str(&format!(",u4n_qmc_eta{eta:.0},u4n_ed_eta{eta:.0},n_mean_eta{eta:.0}"));
    }
    csv.push('\n');
    for (li,&l) in lam.iter().enumerate() {
        csv.push_str(&format!("{l:.6}"));
        for (oi,&omega) in OMEGAS.iter().enumerate() {
            let g=l*(omega*DELTA).sqrt();
            let seed=(li as u64)*1000+(oi as u64)*100+BETA as u64;
            let (m1,_,_,_,_,u4)=run(omega,g,seed);
            let u4_ed=ed_u4n(omega,g);
            csv.push_str(&format!(",{u4:.10e},{u4_ed:.10e},{m1:.10e}"));
        }
        csv.push('\n');
        let dt=start.elapsed().as_secs_f64();
        if li%20==0 { eprintln!("  λ={l:.3} [{dt:.0}s]"); }
    }
    let path=format!("{out_dir}/u4n.csv");
    fs::write(&path,&csv).unwrap();
    eprintln!("\n{path} ({:.0}s)", start.elapsed().as_secs_f64());
}
