mod percolation_support;

use std::mem::size_of;

use allocation_counter::measure;
use cmc_rs::{
    Bond, BondBernoulli, BorrowedUndirectedCsr, GraphView, MixedBernoulli, OccupancyState,
    PercolationMode, Probability, ProbabilityField, SiteBernoulli, StaticConfiguration,
    UndirectedGraphView,
};
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

fn report_topology_view_probe() {
    let case = cases()
        .into_iter()
        .find(|case| case.name == "square")
        .expect("benchmark cases must include square");
    let endpoints: Vec<_> = case
        .lattice
        .edges
        .iter()
        .map(|edge| [edge.source, edge.target])
        .collect();
    let mut constructed = None;
    let construction_info = measure(|| {
        constructed = Some(BorrowedUndirectedCsr::new(
            &case.lattice.offsets,
            &case.lattice.neighbors,
            &case.lattice.edge_ids,
            &endpoints,
        ));
    });
    let view = constructed
        .expect("measurement must run constructor")
        .expect("owned lattice must form a valid borrowed view");
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            let edge_sum = (0..view.edge_count()).fold(0usize, |sum, edge| {
                let [left, right] = view.edge_endpoints(view.edge_id(edge).unwrap());
                sum.wrapping_add(left.index()).wrapping_add(right.index())
            });
            let incidence_sum = view.vertex_ids().fold(0usize, |sum, vertex| {
                view.incidences(vertex).fold(sum, |sum, incidence| {
                    sum.wrapping_add(incidence.neighbor.index())
                        .wrapping_add(incidence.edge.index())
                })
            });
            std::hint::black_box((edge_sum, incidence_sum));
        }
    });
    let input_bytes = case.lattice.offsets.len() * size_of::<usize>()
        + case.lattice.neighbors.len() * size_of::<usize>()
        + case.lattice.edge_ids.len() * size_of::<usize>()
        + endpoints.len() * size_of::<[usize; 2]>();
    eprintln!(
        "TOPOLOGY_VIEW_ALLOC topology={} V={} E={} I={} view_header_bytes={} \
         borrowed_payload_bytes=0 validation_scratch_bytes={} validation_allocations={} \
         input_bytes={} input_bytes_per_vertex={:.3} input_bytes_per_edge={:.3} \
         measured_scans={} allocations_per_scan={:.3} allocated_bytes_per_scan={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        view.vertex_count(),
        view.edge_count(),
        case.lattice.n_bonds,
        size_of::<BorrowedUndirectedCsr<'_>>(),
        construction_info.bytes_total,
        construction_info.count_total,
        input_bytes,
        input_bytes as f64 / view.vertex_count() as f64,
        input_bytes as f64 / view.edge_count() as f64,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("probe probability")
}

fn report_law_case(case: &Case) {
    let vertex_count = case.lattice.vertex_count();
    let edge_count = case.lattice.edge_count();
    let vertex_probabilities = (0..vertex_count)
        .map(|index| probability(if index % 2 == 0 { 0.25 } else { 0.75 }))
        .collect::<Vec<_>>();
    let edge_probabilities = (0..edge_count)
        .map(|index| probability(if index % 3 == 0 { 0.2 } else { 0.6 }))
        .collect::<Vec<_>>();

    let site = SiteBernoulli::new(probability(0.5).into());
    let bond = BondBernoulli::new(probability(0.5).into());
    let mixed = MixedBernoulli::new(probability(0.5).into(), probability(0.5).into());
    let heterogeneous = MixedBernoulli::new(
        ProbabilityField::Borrowed(&vertex_probabilities),
        ProbabilityField::Borrowed(&edge_probabilities),
    );
    let endpoint = MixedBernoulli::new(probability(0.0).into(), probability(1.0).into());
    let endpoint_reverse = MixedBernoulli::new(probability(1.0).into(), probability(0.0).into());
    let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0046_325f_414c_4c4f);

    mixed
        .sample(&case.lattice, &mut configuration, &mut rng)
        .expect("probe domains match");
    endpoint
        .sample(&case.lattice, &mut configuration, &mut rng)
        .expect("probe domains match");

    for (name, samples) in [
        ("site-uniform", ALLOCATION_SAMPLES),
        ("bond-uniform", ALLOCATION_SAMPLES),
        ("mixed-uniform", ALLOCATION_SAMPLES),
        ("mixed-heterogeneous", ALLOCATION_SAMPLES),
        ("mixed-endpoint-switch", ALLOCATION_SAMPLES),
    ] {
        let info = measure(|| {
            for _ in 0..samples {
                let result = match name {
                    "site-uniform" => site.sample(&case.lattice, &mut configuration, &mut rng),
                    "bond-uniform" => bond.sample(&case.lattice, &mut configuration, &mut rng),
                    "mixed-uniform" => mixed.sample(&case.lattice, &mut configuration, &mut rng),
                    "mixed-heterogeneous" => {
                        heterogeneous.sample(&case.lattice, &mut configuration, &mut rng)
                    }
                    "mixed-endpoint-switch" => endpoint
                        .sample(&case.lattice, &mut configuration, &mut rng)
                        .and_then(|()| {
                            endpoint_reverse.sample(&case.lattice, &mut configuration, &mut rng)
                        }),
                    _ => unreachable!("fixed probe case"),
                };
                result.expect("probe domains match");
            }
        });
        let owned_mask_bytes = vertex_count + edge_count;
        eprintln!(
            "PERCOLATION_LAW_ALLOC topology={} law={} V={} E={} measured_iterations={} \
             owned_mask_bytes={} bytes_per_vertex={:.3} bytes_per_edge={:.3} \
             allocations_per_iteration={:.3} allocated_bytes_per_iteration={:.3} \
             peak_live_allocations={} peak_live_bytes={}",
            case.name,
            name,
            vertex_count,
            edge_count,
            samples,
            owned_mask_bytes,
            owned_mask_bytes as f64 / vertex_count.max(1) as f64,
            owned_mask_bytes as f64 / edge_count.max(1) as f64,
            info.count_total as f64 / samples as f64,
            info.bytes_total as f64 / samples as f64,
            info.count_max,
            info.bytes_max,
        );
    }
}

fn main() {
    report_vec_bool_probe();
    report_topology_view_probe();
    let benchmark_cases = cases();
    for case in &benchmark_cases {
        report_law_case(case);
    }
    for case in benchmark_cases {
        for mode in [
            PercolationMode::Site,
            PercolationMode::Bond,
            PercolationMode::SiteBond,
        ] {
            report_case(&case, mode);
        }
    }
}
