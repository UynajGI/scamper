mod percolation_support;

use std::mem::size_of;

use allocation_counter::measure;
use cmc_rs::{Bond, OccupancyState, PercolationMode};
use percolation_support::{cases, sample_and_analyze, Case, P_BOND, P_SITE};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const ALLOCATION_SAMPLES: usize = 128;

fn estimated_owned_storage_bytes(case: &Case, occupancy: &OccupancyState) -> usize {
    case.lattice.offsets.capacity() * size_of::<usize>()
        + case.lattice.neighbors.capacity() * size_of::<usize>()
        + case.lattice.edge_ids.capacity() * size_of::<usize>()
        + case.lattice.edges.capacity() * size_of::<Bond>()
        + occupancy.site_open.capacity()
        + occupancy.bond_open.capacity()
}

fn report_vec_bool_probe() {
    for logical_len in [1usize, 63, 64, 65, 4096, 8064, 11_520, 16_384] {
        let mut retained = None;
        let info = measure(|| {
            let bits = vec![false; logical_len];
            retained = Some((bits.len(), bits.capacity()));
            std::hint::black_box(&bits);
        });
        let (len, logical_capacity) = retained.expect("probe must construct Vec<bool>");
        eprintln!(
            "VEC_BOOL_PROBE len={} logical_capacity={} byte_backed_estimate={} \
             allocator_count={} allocator_bytes={} live_count={} live_bytes={}",
            len,
            logical_capacity,
            logical_capacity,
            info.count_total,
            info.bytes_total,
            info.count_current,
            info.bytes_current,
        );
    }
}

fn report_case(case: &Case, mode: PercolationMode) {
    let mut occupancy = OccupancyState::new(&case.lattice, mode);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0041_4c4c_4f43);
    for _ in 0..8 {
        std::hint::black_box(sample_and_analyze(case, &mut occupancy, &mut rng));
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            std::hint::black_box(sample_and_analyze(case, &mut occupancy, &mut rng));
        }
    });
    let storage_bytes = estimated_owned_storage_bytes(case, &occupancy);
    eprintln!(
        "PERCOLATION_ALLOC topology={} mode={} V={} E={} p_site={} p_bond={} \
         site_len={} site_logical_capacity={} bond_len={} bond_logical_capacity={} \
         estimated_owned_storage_bytes={} estimated_bytes_per_vertex={:.3} \
         estimated_bytes_per_edge={:.3} measured_samples={} allocations_per_sample={:.3} \
         allocated_bytes_per_sample={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.as_label(),
        case.lattice.n_sites,
        case.lattice.n_edges(),
        P_SITE,
        P_BOND,
        occupancy.site_open.len(),
        occupancy.site_open.capacity(),
        occupancy.bond_open.len(),
        occupancy.bond_open.capacity(),
        storage_bytes,
        storage_bytes as f64 / case.lattice.n_sites as f64,
        storage_bytes as f64 / case.lattice.n_edges() as f64,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

fn main() {
    report_vec_bool_probe();
    for case in cases() {
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            report_case(&case, mode);
        }
    }
}
