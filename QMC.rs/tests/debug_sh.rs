use qmc_rs::{MonteCarlo, Context, HeisenbergModel, SSECore};
use qmc_rs::lattice::builders::build_chain;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::SeedableRng;

#[test]
fn debug_max_length() {
    let n_sites = 8;
    let beta = 4.0;
    let lattice = build_chain(n_sites, true);
    let model = HeisenbergModel::new(lattice, beta, 1.0);
    let core = SSECore::new(model);

    println!("N = {}", n_sites);
    println!("beta = {}", beta);
    println!("max_length = {}", core.engine.op_seq.max_length);
    println!("Expected: N*beta*2 + 100 = {}", (n_sites as f64 * beta * 2.0) as usize + 100);

    // Count diagonal vs off-diagonal
    let mut n_diag = 0;
    let mut n_offdiag = 0;
    let mut n_id = 0;
    for v in &core.engine.op_seq.vertices {
        if v.vertex_idx == 0 { n_id += 1; }
        else if v.vertex_idx <= 4 { n_diag += 1; }
        else { n_offdiag += 1; }
    }
    println!("Initial: n_diag={}, n_offdiag={}, n_id={}", n_diag, n_offdiag, n_id);
}

#[test]
fn debug_operator_types() {
    let n_sites = 8;
    let beta = 4.0;
    let lattice = build_chain(n_sites, true);
    let model = HeisenbergModel::new(lattice, beta, 1.0);
    let mut core = SSECore::new(model);

    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    for i in 0..100 {
        core.sweep(&mut ctx);
        if (i+1) % 20 == 0 {
            let mut n_diag = 0;
            let mut n_offdiag = 0;
            let mut n_id = 0;
            for v in &core.engine.op_seq.vertices {
                if v.vertex_idx == 0 { n_id += 1; }
                else if v.vertex_idx <= 4 { n_diag += 1; }
                else { n_offdiag += 1; }
            }
            let n = core.engine.op_seq.n_operators;
            println!("Sweep {}: n={} (diag={}, offdiag={}, id={})",
                i+1, n, n_diag, n_offdiag, n_id);
        }
    }
}
