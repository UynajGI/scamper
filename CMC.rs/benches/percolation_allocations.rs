mod percolation_support;

use std::mem::size_of;

use allocation_counter::measure;
use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze, analyze_winding, overlap_pairs, reweight, reweight_heterogeneous, scan_bond,
    scan_site, AttackGrid, AttackPolicy, Bond, BondBernoulli, BooleanPercolationMC,
    BorrowedUndirectedCsr, BoundaryQuery, ComponentWorkspace, DegreeAttackMC, GraphView,
    HetReweightWorkspace, MixedBernoulli, NewmanZiffWorkspace, ObservablePlan, OrthorhombicCell,
    PeriodicEmbedding, PointProcess, Probability, ProbabilityField, Radii, RandomFailureMC,
    ReweightWorkspace, ScanObservable, SiteBernoulli, StaticConfiguration, StaticLaw,
    StaticObservable, StaticPercolationMC, UndirectedGraphView, WindingWorkspace,
};
use percolation_support::{
    cases, sample_and_analyze, torus_cases, Case, ReferenceConfiguration, ReferenceMode, TorusCase,
    P_BOND, P_SITE,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

const ALLOCATION_SAMPLES: usize = 128;
const ADAPTER_WARMUP: usize = 8;

fn estimated_owned_storage_bytes(case: &Case, occupancy: &ReferenceConfiguration) -> usize {
    case.lattice.offsets.capacity() * size_of::<usize>()
        + case.lattice.neighbors.capacity() * size_of::<usize>()
        + case.lattice.edge_ids.capacity() * size_of::<usize>()
        + case.lattice.edges.capacity() * size_of::<Bond>()
        + occupancy.storage_bytes()
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

fn report_case(case: &Case, mode: ReferenceMode) {
    let mut occupancy = ReferenceConfiguration::new(case, mode);
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
         estimated_owned_storage_bytes={} estimated_bytes_per_vertex={:.3} \
         estimated_bytes_per_edge={:.3} measured_samples={} allocations_per_sample={:.3} \
         allocated_bytes_per_sample={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.label(),
        case.lattice.n_sites,
        case.lattice.n_edges(),
        P_SITE,
        P_BOND,
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

/// Adapter probe: `StaticPercolationMC::sweep` + `measure` must be
/// allocation-free in steady state. The first `measure` of each observable
/// name allocates inside `Measurements` (HashMap insert, first bin), and every
/// `bin_capacity`-th sample allocates one completed bin, so the context uses a
/// binsize strictly larger than warmup plus measured samples and both are
/// exhausted before the measured window.
fn report_adapter_case(case: &Case, mode: ReferenceMode) {
    let law = match mode {
        ReferenceMode::Site => StaticLaw::site(P_SITE),
        ReferenceMode::Bond => StaticLaw::bond(P_BOND),
        ReferenceMode::Mixed => StaticLaw::mixed(P_SITE, P_BOND),
    }
    .expect("probe probability");
    let plan = ObservablePlan::all(vec![BoundaryQuery::new(
        &case.lattice,
        &case.from,
        &case.to,
    )
    .expect("probe boundary query must be valid")]);
    let mut adapter = StaticPercolationMC::new(case.lattice.clone(), law, plan)
        .expect("probe adapter must construct");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x4633_4144_4150),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    let observable_count = 7 + 1 + 1; // plan scalars + defined indicator + one query
    eprintln!(
        "PERCOLATION_ADAPTER_ALLOC topology={} law={} V={} E={} observables={} \
         measured_samples={} allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.label(),
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

fn report_scan_case(case: &Case, mode: ReferenceMode) {
    let probability = Probability::new(P_BOND).expect("probe probability");
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut reweight_workspace = ReweightWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5343_414e);
    for _ in 0..ADAPTER_WARMUP {
        run_scan_once(
            case,
            mode,
            &mut scan_workspace,
            &mut reweight_workspace,
            &mut rng,
            probability,
        );
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            run_scan_once(
                case,
                mode,
                &mut scan_workspace,
                &mut reweight_workspace,
                &mut rng,
                probability,
            );
        }
    });
    let scan_capacity = scan_workspace.capacity_audit();
    let reweight_capacity = reweight_workspace.capacity_audit();
    let total_bytes = scan_capacity
        .total_bytes()
        .saturating_add(reweight_capacity.total_bytes());
    let elements = case
        .lattice
        .vertex_count()
        .saturating_add(case.lattice.edge_count());
    eprintln!(
        "PERCOLATION_SCAN_ALLOC topology={} mode={} V={} E={} p={} measured_samples={} \
         workspace_parent_bytes={} workspace_size_bytes={} workspace_vertex_active_bytes={} \
         workspace_edge_counted_bytes={} workspace_order_bytes={} workspace_points_bytes={} \
         workspace_scan_total_capacity_bytes={} workspace_reweight_weights_bytes={} \
         workspace_total_capacity_bytes={} workspace_bytes_per_vertex={:.3} \
         workspace_bytes_per_element={:.3} allocations_per_sample={:.3} \
         allocated_bytes_per_sample={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.label(),
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        P_BOND,
        ALLOCATION_SAMPLES,
        scan_capacity.parent,
        scan_capacity.size,
        scan_capacity.vertex_active,
        scan_capacity.edge_counted,
        scan_capacity.order,
        scan_capacity.points,
        scan_capacity.total_bytes(),
        reweight_capacity.weights,
        total_bytes,
        total_bytes as f64 / case.lattice.vertex_count().max(1) as f64,
        total_bytes as f64 / elements.max(1) as f64,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

fn run_scan_once(
    case: &Case,
    mode: ReferenceMode,
    scan_workspace: &mut NewmanZiffWorkspace,
    reweight_workspace: &mut ReweightWorkspace,
    rng: &mut Xoshiro256PlusPlus,
    probability: Probability,
) {
    let curve = match mode {
        ReferenceMode::Site => {
            scan_site(&case.lattice, scan_workspace, rng).expect("probe site scan")
        }
        ReferenceMode::Bond => {
            scan_bond(&case.lattice, scan_workspace, rng).expect("probe bond scan")
        }
        ReferenceMode::Mixed => unreachable!("scan probe covers pure modes"),
    };
    let value = reweight(
        &curve,
        probability,
        ScanObservable::LargestComponentSize,
        reweight_workspace,
    )
    .expect("probe reweight");
    std::hint::black_box(value);
}

/// T2 probe: one heterogeneous reweighting realization (site/bond scan plus
/// one Poisson-binomial reweight at a per-element probability field) must be
/// allocation-free in steady state after one prepare of each workspace
/// (same policy as the kcore/bootstrap/incremental probes).
fn report_heterogeneous_scan_case(case: &Case, mode: ReferenceMode) {
    let trials = match mode {
        ReferenceMode::Site => case.lattice.vertex_count(),
        ReferenceMode::Bond => case.lattice.edge_count(),
        ReferenceMode::Mixed => unreachable!("scan probe covers pure modes"),
    };
    let values: Vec<Probability> = (0..trials)
        .map(|index| probability(if index % 2 == 0 { 0.35 } else { 0.75 }))
        .collect();
    let probability_field = ProbabilityField::Borrowed(&values);
    let mut scan_workspace = NewmanZiffWorkspace::new();
    let mut reweight_workspace = HetReweightWorkspace::new();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4845_5A31);
    let run_once = |scan_workspace: &mut NewmanZiffWorkspace,
                    reweight_workspace: &mut HetReweightWorkspace,
                    rng: &mut Xoshiro256PlusPlus| {
        let curve = match mode {
            ReferenceMode::Site => {
                scan_site(&case.lattice, scan_workspace, rng).expect("probe site scan")
            }
            ReferenceMode::Bond => {
                scan_bond(&case.lattice, scan_workspace, rng).expect("probe bond scan")
            }
            ReferenceMode::Mixed => unreachable!("scan probe covers pure modes"),
        };
        let result = reweight_heterogeneous(
            &curve,
            &probability_field,
            ScanObservable::LargestComponentSize,
            reweight_workspace,
        )
        .expect("probe heterogeneous reweight");
        std::hint::black_box(result);
    };
    for _ in 0..ADAPTER_WARMUP {
        run_once(&mut scan_workspace, &mut reweight_workspace, &mut rng);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            run_once(&mut scan_workspace, &mut reweight_workspace, &mut rng);
        }
    });
    assert_eq!(
        info.count_total, 0,
        "heterogeneous reweight realizations allocated {} times in steady state",
        info.count_total
    );
    let capacity = reweight_workspace.capacity_audit();
    eprintln!(
        "HET_SCAN_ALLOC topology={} mode={} V={} E={} field_entries={} measured_samples={} \
         workspace_log_weights_bytes={} workspace_log_integers_bytes={} \
         workspace_uniform_weights_bytes={} workspace_het_total_capacity_bytes={} \
         allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.label(),
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        trials,
        ALLOCATION_SAMPLES,
        capacity.log_weights,
        capacity.log_integers,
        capacity.uniform,
        capacity.total_bytes(),
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

fn report_component_case(case: &Case, mode: ReferenceMode) {
    let vertex_count = case.lattice.vertex_count();
    let edge_count = case.lattice.edge_count();
    let query = BoundaryQuery::new(&case.lattice, &case.from, &case.to)
        .expect("probe boundary query must be valid");
    let queries = [query];
    let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut prepared_workspace = None;
    let workspace_prepare_info = measure(|| {
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(vertex_count, queries.len(), true)
            .expect("probe workspace reservation");
        prepared_workspace = Some(workspace);
    });
    let mut workspace = prepared_workspace.expect("probe must retain prepared workspace");
    let workspace_capacity = workspace.capacity_audit();
    let site = SiteBernoulli::new(probability(P_SITE).into());
    let bond = BondBernoulli::new(probability(P_BOND).into());
    let mixed = MixedBernoulli::new(probability(P_SITE).into(), probability(P_BOND).into());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_414c_4c4f);
    let sample = |configuration: &mut StaticConfiguration, rng: &mut Xoshiro256PlusPlus| {
        match mode {
            ReferenceMode::Site => site.sample(&case.lattice, configuration, rng),
            ReferenceMode::Bond => bond.sample(&case.lattice, configuration, rng),
            ReferenceMode::Mixed => mixed.sample(&case.lattice, configuration, rng),
        }
        .expect("probe dimensions match");
    };
    for _ in 0..8 {
        sample(&mut configuration, &mut rng);
        std::hint::black_box(
            analyze(
                &case.lattice,
                &configuration,
                &configuration,
                &queries,
                &mut workspace,
            )
            .expect("probe analysis dimensions match"),
        );
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            sample(&mut configuration, &mut rng);
            std::hint::black_box(
                analyze(
                    &case.lattice,
                    &configuration,
                    &configuration,
                    &queries,
                    &mut workspace,
                )
                .expect("probe analysis dimensions match"),
            );
        }
    });
    eprintln!(
        "PERCOLATION_COMPONENT_ALLOC topology={} law={} V={} E={} measured_samples={} \
         workspace_parent_bytes={} workspace_size_bytes={} workspace_canonical_bytes={} \
         workspace_stamp_bytes={} workspace_query_output_bytes={} workspace_label_bytes={} \
         workspace_total_capacity_bytes={} workspace_bytes_per_vertex={:.3} \
         workspace_bytes_per_query={:.3} workspace_prepare_allocations={} \
         workspace_prepare_allocator_bytes={} allocations_per_sample={:.3} \
         allocated_bytes_per_sample={:.3} peak_live_allocations={} peak_live_bytes={} \
         f0_reference_allocations_per_sample=7 f0_reference_bytes_per_sample=86016",
        case.name,
        mode.label(),
        vertex_count,
        edge_count,
        ALLOCATION_SAMPLES,
        workspace_capacity.parent,
        workspace_capacity.size,
        workspace_capacity.canonical,
        workspace_capacity.query_root_stamps,
        workspace_capacity.query_outcomes,
        workspace_capacity.labels,
        workspace_capacity.total_bytes(),
        workspace_capacity.total_bytes() as f64 / vertex_count.max(1) as f64,
        workspace_capacity.total_bytes() as f64 / queries.len().max(1) as f64,
        workspace_prepare_info.count_total,
        workspace_prepare_info.bytes_total,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// Heterogeneous (N1) adapter probe: `StaticPercolationMC` sweep+measure with
/// per-vertex and per-edge probability fields must be allocation-free in
/// steady state, exactly like the uniform adapter above.
fn report_heterogeneous_adapter_case(case: &Case, mode: ReferenceMode) {
    let vertex_count = case.lattice.vertex_count();
    let edge_count = case.lattice.edge_count();
    let vertex_probabilities = (0..vertex_count)
        .map(|index| probability(if index % 2 == 0 { 0.25 } else { 0.75 }))
        .collect::<Vec<_>>();
    let edge_probabilities = (0..edge_count)
        .map(|index| probability(if index % 3 == 0 { 0.2 } else { 0.6 }))
        .collect::<Vec<_>>();
    let law = match mode {
        ReferenceMode::Site => StaticLaw::site_heterogeneous(vertex_probabilities),
        ReferenceMode::Bond => StaticLaw::bond_heterogeneous(edge_probabilities),
        ReferenceMode::Mixed => {
            StaticLaw::mixed_heterogeneous(vertex_probabilities, edge_probabilities)
        }
    };
    let plan = ObservablePlan::all(vec![BoundaryQuery::new(
        &case.lattice,
        &case.from,
        &case.to,
    )
    .expect("probe boundary query must be valid")]);
    let mut adapter = StaticPercolationMC::new(case.lattice.clone(), law, plan)
        .expect("probe adapter must construct");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x4845_5445_4150),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    let observable_count = 7 + 1 + 1; // plan scalars + defined indicator + one query
    eprintln!(
        "PERCOLATION_ADAPTER_ALLOC topology={} law=heterogeneous-{} V={} E={} observables={} \
         measured_samples={} allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        mode.label(),
        vertex_count,
        edge_count,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// W2 winding-workspace probe: law sampling plus `analyze_winding` on a
/// builder torus must be allocation-free in steady state after one
/// `WindingWorkspace::prepare`.
fn report_winding_case(torus_case: &TorusCase) {
    let vertex_count = torus_case.lattice.vertex_count();
    let edge_count = torus_case.lattice.edge_count();
    let dimensions = torus_case.embedding.periodic_dimensions();
    let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut prepared = None;
    let prepare_info = measure(|| {
        let mut workspace = WindingWorkspace::new();
        workspace
            .prepare(vertex_count, dimensions)
            .expect("probe winding workspace reservation");
        prepared = Some(workspace);
    });
    let mut workspace = prepared.expect("probe must retain prepared workspace");
    let bond = BondBernoulli::new(probability(P_BOND).into());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_5749_4e44);
    for _ in 0..8 {
        bond.sample(&torus_case.lattice, &mut configuration, &mut rng)
            .expect("probe dimensions match");
        std::hint::black_box(
            analyze_winding(
                &torus_case.lattice,
                &torus_case.embedding,
                &configuration,
                &configuration,
                &mut workspace,
            )
            .expect("probe winding analysis"),
        );
    }
    // Audit after warmup so lazily grown storage (component table) is
    // included in the retained-capacity accounting.
    let capacity = workspace.capacity_audit();
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            bond.sample(&torus_case.lattice, &mut configuration, &mut rng)
                .expect("probe dimensions match");
            std::hint::black_box(
                analyze_winding(
                    &torus_case.lattice,
                    &torus_case.embedding,
                    &configuration,
                    &configuration,
                    &mut workspace,
                )
                .expect("probe winding analysis"),
            );
        }
    });
    eprintln!(
        "WINDING_WORKSPACE_ALLOC topology={} D={} V={} E={} measured_samples={} \
         workspace_parent_bytes={} workspace_size_bytes={} workspace_canonical_bytes={} \
         workspace_potential_bytes={} workspace_basis_bytes={} workspace_basis_count_bytes={} \
         workspace_root_slot_bytes={} workspace_scratch_bytes={} workspace_chain_bytes={} \
         workspace_components_bytes={} workspace_total_capacity_bytes={} \
         workspace_bytes_per_vertex={:.3} workspace_prepare_allocations={} \
         workspace_prepare_allocator_bytes={} allocations_per_sample={:.3} \
         allocated_bytes_per_sample={:.3} peak_live_allocations={} peak_live_bytes={}",
        torus_case.name,
        dimensions,
        vertex_count,
        edge_count,
        ALLOCATION_SAMPLES,
        capacity.parent,
        capacity.size,
        capacity.canonical,
        capacity.potential,
        capacity.basis,
        capacity.basis_counts,
        capacity.root_slot,
        capacity.scratch,
        capacity.chain,
        capacity.components,
        capacity.total_bytes(),
        capacity.total_bytes() as f64 / vertex_count.max(1) as f64,
        prepare_info.count_total,
        prepare_info.bytes_total,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// W2 adapter probe: `StaticPercolationMC` sweep + measure with the winding
/// observables requested must be allocation-free in steady state, exactly
/// like the plain adapter above.
fn report_winding_adapter_case(torus_case: &TorusCase) {
    let law = StaticLaw::bond(P_BOND).expect("probe probability");
    let plan = ObservablePlan::new(
        [
            StaticObservable::ActiveVertexCount,
            StaticObservable::ActiveEdgeCount,
            StaticObservable::ComponentCount,
            StaticObservable::LargestSize,
            StaticObservable::GiantFraction,
            StaticObservable::RawSecondMoment,
            StaticObservable::FiniteClusterSusceptibility,
            StaticObservable::Wrapping,
            StaticObservable::WindingRank,
        ],
        Vec::new(),
    );
    let mut adapter = StaticPercolationMC::with_embedding(
        torus_case.lattice.clone(),
        torus_case.embedding.clone(),
        law,
        plan,
    )
    .expect("probe winding adapter must construct");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x4633_5741_5053),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    eprintln!(
        "WINDING_ADAPTER_ALLOC topology={} law=bond-uniform V={} E={} observables=9+1 \
         measured_samples={} allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        torus_case.name,
        torus_case.lattice.vertex_count(),
        torus_case.lattice.edge_count(),
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// N2 facade probe: `RandomFailureMC::sweep` + `measure` (labels pass plus
/// six degree arrays) must be allocation-free in steady state under the same
/// binsize discipline as the adapter probe.
fn report_random_failure_case(case: &Case) {
    let plan = ObservablePlan::all(Vec::new());
    let mut facade =
        RandomFailureMC::new(case.lattice.clone(), P_SITE, plan).expect("probe facade");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x4633_4144_4150),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        facade.sweep(&mut context);
        facade.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            facade.sweep(&mut context);
            facade.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    // 7 plan scalars + susceptibility indicator + 6 degree arrays.
    let observable_count = 7 + 1 + 6;
    eprintln!(
        "PERCOLATION_RANDOM_FAILURE_ALLOC topology={} V={} E={} observables={} \
         measured_samples={} allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// N3 adapter probe: `DegreeAttackMC::sweep` + `measure` (one full attack
/// trajectory: order construction plus a component analysis per
/// checkpoint) must be allocation-free in steady state; buffers are
/// retained after the first trajectory and the grid keeps ~17 checkpoints
/// so the probe stays fast.
fn report_degree_attack_case(case: &Case) {
    let step = (case.lattice.vertex_count() / 16).max(1);
    let mut adapter = DegreeAttackMC::new(
        case.lattice.clone(),
        AttackPolicy::DegreeDescendingSeededShuffle,
        AttackGrid::Every { step },
    )
    .expect("probe adapter");
    let checkpoints = adapter.checkpoints().len();
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0xa11a),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    // RobustnessCurve + RobustnessRemovedFraction + susceptibility
    // indicator + susceptibility array (defined for these mid-range grids).
    let observable_count = 4;
    eprintln!(
        "PERCOLATION_DEGREE_ATTACK_ALLOC topology={} V={} E={} checkpoints={} \
         observables={} measured_trajectories={} allocations_per_trajectory={:.3} \
         allocated_bytes_per_trajectory={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        checkpoints,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// D2 probe: repeated weak/strong/reachability analyses on prepared
/// workspaces must be allocation-free in steady state.
fn report_directed_analysis_case(case: &Case) {
    use cmc_rs::{
        analyze_strong_components, analyze_weak_components, reachable, Arc, DirectedCsr,
        DirectedGraphView, ReachabilityWorkspace, StrongComponentWorkspace, WeakComponentWorkspace,
    };
    let vertices = case.lattice.vertex_count();
    let edge_count = case.lattice.edge_count();
    // Orient every undirected bond forward with occasional back arcs.
    let mut arcs = Vec::with_capacity(edge_count + edge_count / 8);
    for edge in &case.lattice.edges {
        let (source, target) = if edge.source < edge.target {
            (edge.source, edge.target)
        } else {
            (edge.target, edge.source)
        };
        arcs.push(Arc::new(source, target));
        if source % 8 == 0 && source > 0 {
            arcs.push(Arc::new(target, source));
        }
    }
    let graph = DirectedCsr::from_arcs(vertices, arcs).expect("valid arcs");
    let mut weak = WeakComponentWorkspace::default();
    let mut strong = StrongComponentWorkspace::default();
    let mut reach = ReachabilityWorkspace::default();
    let goal = graph.vertex_id(vertices - 1).unwrap();
    for _ in 0..ADAPTER_WARMUP {
        analyze_weak_components(&graph, &mut weak).unwrap();
        analyze_strong_components(&graph, &mut strong).unwrap();
        reachable(&graph, graph.vertex_id(0).unwrap(), goal, &mut reach).unwrap();
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            analyze_weak_components(&graph, &mut weak).unwrap();
            analyze_strong_components(&graph, &mut strong).unwrap();
            reachable(&graph, graph.vertex_id(0).unwrap(), goal, &mut reach).unwrap();
            std::hint::black_box((&mut weak, &mut strong, &mut reach));
        }
    });
    eprintln!(
        "DIRECTED_ANALYSIS_ALLOC topology={} V={} A={} measured_rounds={} \
         allocations_per_round={:.3} allocated_bytes_per_round={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        case.name,
        vertices,
        graph.arc_count(),
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// D3 probe: one complete layered directed-percolation realization
/// (frontier bitset propagation over all layers) must be allocation-free
/// in steady state.
fn report_directed_process_case(case: &Case) {
    use cmc_rs::{DpSeedPolicy, LayeredDirectedPercolationMC};
    let width = case.lattice.vertex_count();
    let mut process = LayeredDirectedPercolationMC::new(
        case.lattice.clone(),
        0.6447001,
        256,
        DpSeedPolicy::Single { site: width / 2 },
    )
    .expect("probe process");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0xd3a1),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        process.sweep(&mut context);
        process.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            process.sweep(&mut context);
            process.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    // DpLayerActive + DpSurvival arrays + extinction scalar.
    let observable_count = 3;
    eprintln!(
        "DIRECTED_PROCESS_ALLOC topology={} V={} layers=256 observables={} \
         measured_realizations={} allocations_per_realization={:.3} \
         allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        width,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// E3 probe: one complete k-core realization (site-keep activity draw,
/// residual-degree peeling to completion, summary recording) must be
/// allocation-free in steady state.
fn report_kcore_case(case: &Case) {
    use cmc_rs::{CoreThresholds, KCoreMC};
    let vertices = case.lattice.vertex_count();
    let mut process = KCoreMC::new(
        case.lattice.clone(),
        P_SITE,
        CoreThresholds::uniform(3).expect("probe threshold"),
        256,
    )
    .expect("probe process");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0xE3A1),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        process.sweep(&mut context);
        process.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            process.sweep(&mut context);
            process.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    // KCoreSurvivorCount array + three scalars.
    let observable_count = 4;
    // Review follow-up 2026-10-23: assert like the bootstrap probe — a
    // silent regression to nonzero allocations must fail the probe, not
    // merely print.
    assert_eq!(
        info.count_total, 0,
        "kcore realization allocated {} times in steady state",
        info.count_total
    );
    eprintln!(
        "KCORE_PROCESS_ALLOC topology={} V={} E={} k=3 keep={} recorded_rounds=256 \
         observables={} measured_realizations={} allocations_per_realization={:.3} \
         allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        vertices,
        case.lattice.edge_count(),
        P_SITE,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// Adapter probe: `BootstrapPercolationMC::sweep` + `measure` must be
/// allocation-free in steady state (E1 buffers prepared once). The context
/// binsize exceeds warmup plus measured samples so no accumulator bin
/// completes inside the measured window; the warmup covers the first
/// `measure` of each observable name.
fn report_bootstrap_case(case: &Case) {
    use cmc_rs::{BootstrapDynamics, BootstrapPercolationMC, BootstrapSeeds, BootstrapThreshold};
    let vertices = case.lattice.vertex_count();
    let mut process = BootstrapPercolationMC::new(
        case.lattice.clone(),
        BootstrapThreshold::Homogeneous { value: 2 },
        BootstrapSeeds::Density { p: P_SITE },
        BootstrapDynamics::SynchronousRounds,
        32,
    )
    .expect("probe process");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x00E2),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        process.sweep(&mut context);
        process.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            process.sweep(&mut context);
            process.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    assert_eq!(
        info.count_total, 0,
        "bootstrap synchronous realizations must allocate nothing in steady state"
    );

    // Review follow-up 2026-10-23: assert the QueueClosure clock too —
    // it shares the same E1 buffers; this pins the symmetry claim.
    let mut process = BootstrapPercolationMC::new(
        case.lattice.clone(),
        BootstrapThreshold::Homogeneous { value: 2 },
        BootstrapSeeds::Density { p: P_SITE },
        BootstrapDynamics::QueueClosure,
        32,
    )
    .expect("probe process (queue clock)");
    for _ in 0..ADAPTER_WARMUP {
        process.sweep(&mut context);
        process.measure(&mut context);
    }
    let queue_info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            process.sweep(&mut context);
            process.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    assert_eq!(
        queue_info.count_total, 0,
        "bootstrap queue-closure realizations must allocate nothing in steady state"
    );

    // BootstrapActiveCount array + two always-on scalars + the clock
    // scalar.
    let observable_count = 4;
    eprintln!(
        "BOOTSTRAP_PROCESS_ALLOC topology={} V={} E={} r=2 seed_density={} steps=32 \
         observables={} measured_realizations={} allocations_per_realization={:.3} \
         allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        vertices,
        case.lattice.edge_count(),
        P_SITE,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// G1 probe: one complete incremental realization (reset, every edge
/// activated, every vertex activated in a fresh random order — V + E
/// accepted operations, each appending one trajectory entry) must be
/// allocation-free in steady state on a reused workspace.
fn report_incremental_case(case: &Case) {
    use cmc_rs::{GraphView, IncrementalUnionWorkspace, UndirectedGraphView};
    use rand::RngExt;
    let vertices = case.lattice.vertex_count();
    let edges = case.lattice.edge_count();
    let mut workspace = IncrementalUnionWorkspace::new();
    workspace
        .prepare(vertices, edges)
        .expect("probe workspace reservation");
    let mut order: Vec<usize> = (0..vertices).collect();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00C1_6731);
    let run_once = |workspace: &mut IncrementalUnionWorkspace,
                    order: &mut Vec<usize>,
                    rng: &mut Xoshiro256PlusPlus| {
        workspace.reset();
        for index in 0..edges {
            let edge = case.lattice.edge_id(index).expect("dense edge index");
            workspace
                .activate_edge(&case.lattice, edge, index)
                .expect("probe edge activation");
        }
        for index in (1..order.len()).rev() {
            let swap_with = rng.random_range(0..=index);
            order.swap(index, swap_with);
        }
        for (step, &vertex) in order.iter().enumerate() {
            let vertex_id = case.lattice.vertex_id(vertex).expect("dense vertex index");
            workspace
                .activate_vertex(&case.lattice, vertex_id, step)
                .expect("probe vertex activation");
        }
        std::hint::black_box(workspace.trajectory().summaries().last().copied());
    };
    for _ in 0..8 {
        run_once(&mut workspace, &mut order, &mut rng);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            run_once(&mut workspace, &mut order, &mut rng);
        }
    });
    // The substrate's core promise: zero allocations per operation after
    // warmup — a silent regression must fail the probe, not merely print
    // (same policy as the kcore/bootstrap probes).
    assert_eq!(
        info.count_total, 0,
        "incremental realizations allocated {} times in steady state",
        info.count_total
    );
    let ops = vertices + edges;
    eprintln!(
        "INCREMENTAL_PROCESS_ALLOC topology={} V={} E={} ops_per_realization={} \
         measured_realizations={} allocations_per_realization={:.3} \
         allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        vertices,
        edges,
        ops,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// G2 probe: one complete Achlioptas realization (vertex re-initialization +
/// E candidate steps: 2 draws, 2 `component_size` scorings, 1 union each,
/// running to complete occupation) must be allocation-free in steady state —
/// every buffer (pool, candidates, added flags, accepted order, curve) is
/// reserved at construction. Asserts zero (same policy as the
/// kcore/bootstrap/incremental probes).
fn report_achlioptas_case(case: &Case) {
    use cmc_rs::{
        AchlioptasGrid, AchlioptasProcessMC, AchlioptasRule, CandidateSampling, EdgeUniverse,
        TiePolicy,
    };
    let edges = case.lattice.edge_count();
    let mut process = AchlioptasProcessMC::new(
        case.lattice.clone(),
        AchlioptasRule::Product,
        CandidateSampling::WithoutReplacement,
        2,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
        AchlioptasGrid::Every {
            step: (edges / 64).max(1),
        },
    )
    .expect("probe process");
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00C2_414C);
    for _ in 0..ADAPTER_WARMUP {
        process.run_realization(&mut rng);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            process.run_realization(&mut rng);
            std::hint::black_box(process.largest_curve());
        }
    });
    assert_eq!(
        info.count_total, 0,
        "achlioptas realizations allocated {} times in steady state",
        info.count_total
    );
    eprintln!(
        "ACHILLOPTAS_PROCESS_ALLOC topology={} V={} E={} m=2 rule=product \
         sampling=without-replacement steps_per_realization={} \
         measured_realizations={} allocations_per_realization={:.3} \
         allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        case.lattice.vertex_count(),
        edges,
        edges,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

/// G3 probe: one complete invasion realization (frozen threshold draw,
/// binary-heap minimum-frontier growth to exhaustion, G1 substrate drive,
/// C0 history compression, measurement) must be allocation-free in steady
/// state — the heap, invasion state, trapping scratch, and substrate are
/// reserved at construction. Site, bond, and the trapping variant are all
/// asserted (same policy as the kcore/bootstrap/incremental/achlioptas
/// probes).
fn report_invasion_case(case: &Case) {
    use cmc_rs::{InvasionMode, InvasionPercolationMC, InvasionStop, InvasionTrapping};
    let vertices = case.lattice.vertex_count();
    let variants: [(&str, InvasionMode, InvasionTrapping); 3] = [
        ("site", InvasionMode::Site, InvasionTrapping::NonTrapping),
        ("bond", InvasionMode::Bond, InvasionTrapping::NonTrapping),
        (
            "site-trapping",
            InvasionMode::Site,
            InvasionTrapping::Trapping {
                exit: vec![vertices - 1],
            },
        ),
    ];
    for (label, mode, trapping) in variants {
        let mut process = InvasionPercolationMC::new(
            case.lattice.clone(),
            mode,
            vec![0],
            InvasionStop::GrowAll,
            trapping,
            vertices.saturating_sub(1).max(1),
        )
        .expect("probe process");
        let mut context = Context::new_with_binsize(
            Xoshiro256PlusPlus::seed_from_u64(0x00C3_1A5A),
            0,
            ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
        );
        for _ in 0..ADAPTER_WARMUP {
            process.sweep(&mut context);
            process.measure(&mut context);
        }
        let info = measure(|| {
            for _ in 0..ALLOCATION_SAMPLES {
                process.sweep(&mut context);
                process.measure(&mut context);
                std::hint::black_box(&mut context);
            }
        });
        assert_eq!(
            info.count_total, 0,
            "invasion ({label}) realizations allocated {} times in steady state",
            info.count_total
        );
        // InvasionThresholdOrder array + three scalars.
        let observable_count = 4;
        eprintln!(
            "INVASION_PROCESS_ALLOC topology={} variant={} V={} E={} observables={} \
             measured_realizations={} allocations_per_realization={:.3} \
             allocated_bytes_per_realization={:.3} peak_live_allocations={} peak_live_bytes={}",
            case.name,
            label,
            vertices,
            case.lattice.edge_count(),
            observable_count,
            ALLOCATION_SAMPLES,
            info.count_total as f64 / ALLOCATION_SAMPLES as f64,
            info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
            info.count_max,
            info.bytes_max,
        );
    }
}

/// P2 probe: one complete Boolean realization (binomial point-process draw +
/// overlap construction via P1 + F3 analysis + plan record) must perform
/// exactly the allocations of the documented P1 transient construction —
/// zero beyond it. The transient (coordinate copy, cell-list buckets, pair
/// list) allocates EVERY realization by design under the no-persistent-CSR
/// rule, so warmup cannot exclude it; instead the probe replays the same
/// construction in isolation over identical coordinate snapshots and asserts
/// allocation-count equality: every steady-state allocation is attributable
/// to the documented transient and the adapter adds none. The frozen
/// measuring mode gets the same equality assertion (its sweep is a
/// documented no-op). Binomial counts keep the realized N fixed; a Poisson
/// adapter differs only by rare tail excursions growing the retained buffer
/// once (documented on the adapter).
fn report_boolean_case() {
    let count = 2048_usize;
    let side = (count as f64).sqrt();
    let cell = OrthorhombicCell::new([side, side]).expect("probe cell");
    let plan = || ObservablePlan::all(Vec::new());
    for (dispersion, radii) in [
        (
            "monodisperse",
            Radii::monodisperse(0.47).expect("probe radius"),
        ),
        (
            "polydisperse",
            Radii::per_particle((0..count).map(|index| [0.25, 0.5, 0.75][index % 3]))
                .expect("probe radii"),
        ),
    ] {
        let mut adapter = BooleanPercolationMC::<2>::sampling(
            cell,
            radii.clone(),
            PointProcess::binomial(count).expect("positive count"),
            plan(),
        )
        .expect("probe sampling adapter");
        let mut context = Context::new_with_binsize(
            Xoshiro256PlusPlus::seed_from_u64(0x0000_B001),
            0,
            ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
        );
        for _ in 0..ADAPTER_WARMUP {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
        }
        // Snapshot every realized coordinate set into a preallocated arena
        // (copying allocates nothing) so the replay window sees identical
        // geometry — pair-list growth patterns included.
        let coordinates = count * 2;
        let mut arena = vec![0.0_f64; ALLOCATION_SAMPLES * coordinates];
        let info_full = measure(|| {
            for sample in 0..ALLOCATION_SAMPLES {
                adapter.sweep(&mut context);
                arena[sample * coordinates..(sample + 1) * coordinates]
                    .copy_from_slice(adapter.coordinates());
                adapter.measure(&mut context);
                std::hint::black_box(&context);
            }
        });
        let info_transient = measure(|| {
            for sample in 0..ALLOCATION_SAMPLES {
                std::hint::black_box(
                    overlap_pairs(
                        &arena[sample * coordinates..(sample + 1) * coordinates],
                        &cell,
                        &radii,
                    )
                    .expect("probe geometry"),
                );
            }
        });
        assert_eq!(
            info_full.count_total, info_transient.count_total,
            "{dispersion} sampling realizations allocated {} times vs {} for the isolated \
             P1 transient; the adapter must add zero allocations",
            info_full.count_total, info_transient.count_total
        );
        eprintln!(
            "BOOLEAN_PERCOLATION_ALLOC mode=sampling radii={dispersion} N={} eta~{:.3} \
             measured_realizations={} transient_allocations_per_realization={:.3} \
             transient_bytes_per_realization={:.3} adapter_added_allocations={} \
             peak_live_allocations={} peak_live_bytes={}",
            count,
            std::f64::consts::PI
                * (0..count)
                    .map(|index| radii.radius(index).powi(2))
                    .sum::<f64>()
                / (side * side),
            ALLOCATION_SAMPLES,
            info_transient.count_total as f64 / ALLOCATION_SAMPLES as f64,
            info_transient.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
            info_full.count_total - info_transient.count_total,
            info_full.count_max,
            info_full.bytes_max,
        );
    }

    // Frozen measuring adapter: the equilibrium-measurement path.
    let mut frozen = Vec::with_capacity(count * 2);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0000_B002);
    for _ in 0..count {
        for axis in 0..2 {
            frozen.push(rng.random_range(0.0..cell.lengths()[axis]));
        }
    }
    let radii = Radii::monodisperse(0.47).expect("probe radius");
    let mut adapter =
        BooleanPercolationMC::<2>::measuring(cell, frozen.clone(), radii.clone(), plan())
            .expect("probe measuring adapter");
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x0000_B003),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info_full = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&context);
        }
    });
    let info_transient = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            std::hint::black_box(overlap_pairs(&frozen, &cell, &radii).expect("probe geometry"));
        }
    });
    assert_eq!(
        info_full.count_total, info_transient.count_total,
        "frozen measurements allocated {} times vs {} for the isolated P1 transient",
        info_full.count_total, info_transient.count_total
    );
    eprintln!(
        "BOOLEAN_PERCOLATION_ALLOC mode=measuring-frozen radii=monodisperse N={} \
         measured_measurements={} transient_allocations_per_measurement={:.3} \
         adapter_added_allocations={} peak_live_allocations={} peak_live_bytes={}",
        count,
        ALLOCATION_SAMPLES,
        info_transient.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info_full.count_total - info_transient.count_total,
        info_full.count_max,
        info_full.bytes_max,
    );
}

/// T4 probe: repeated `analyze_winding` on a FIXED continuum pair graph and
/// `ContinuumEmbedding` with re-drawn (toggled) edge activity must be
/// allocation-free in steady state after one `WindingWorkspace::prepare`
/// (same policy as the W2 torus probe). Construction — P1 pair enumeration
/// plus the embedding derivation — is a per-realization transient like the
/// `PairGraph` itself (the no-persistent-CSR rule), measured here for the
/// record and never asserted zero.
fn report_continuum_winding_case() {
    use cmc_rs::{CandidatePairGraph, ContinuumEmbedding};
    let count = 2048_usize;
    let side = 64.0;
    let cell = OrthorhombicCell::new([side, side]).expect("probe cell");
    let cutoff = 1.0;
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5434_434f);
    let mut coordinates = Vec::with_capacity(count * 2);
    for _ in 0..count {
        coordinates.push(rng.random_range(0.0..side));
        coordinates.push(rng.random_range(0.0..side));
    }
    let mut prepared = None;
    let construction_info = measure(|| {
        let graph = CandidatePairGraph::new(&coordinates, &cell, cutoff)
            .expect("probe coordinates")
            .pairs()
            .expect("probe pairs");
        let embedding = ContinuumEmbedding::from_pairs(graph.pairs(), &coordinates, &cell, cutoff)
            .expect("probe embedding");
        prepared = Some((graph, embedding));
    });
    let (graph, embedding) = prepared.expect("probe must retain the fixed pair graph");
    let vertex_count = graph.vertex_count();
    let edge_count = graph.edge_count();

    let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut workspace = WindingWorkspace::new();
    workspace
        .prepare(vertex_count, 2)
        .expect("probe winding workspace reservation");
    let bond = BondBernoulli::new(probability(P_BOND).into());
    for _ in 0..8 {
        bond.sample(&graph, &mut configuration, &mut rng)
            .expect("probe dimensions match");
        std::hint::black_box(
            analyze_winding(
                &graph,
                &embedding,
                &configuration,
                &configuration,
                &mut workspace,
            )
            .expect("probe winding analysis"),
        );
    }
    let capacity = workspace.capacity_audit();
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            bond.sample(&graph, &mut configuration, &mut rng)
                .expect("probe dimensions match");
            std::hint::black_box(
                analyze_winding(
                    &graph,
                    &embedding,
                    &configuration,
                    &configuration,
                    &mut workspace,
                )
                .expect("probe winding analysis"),
            );
        }
    });
    assert_eq!(
        info.count_total, 0,
        "steady-state continuum winding analysis allocated {} times on a fixed pair graph \
         and embedding",
        info.count_total
    );
    eprintln!(
        "CONTINUUM_WINDING_ALLOC topology=poisson-2d N={} E={} cutoff={} measured_samples={} \
         construction_allocations_per_realization={} \
         construction_allocator_bytes_per_realization={} \
         workspace_total_capacity_bytes={} workspace_bytes_per_vertex={:.3} \
         allocations_per_sample={:.3} allocated_bytes_per_sample={:.3} \
         peak_live_allocations={} peak_live_bytes={}",
        vertex_count,
        edge_count,
        cutoff,
        ALLOCATION_SAMPLES,
        construction_info.count_total,
        construction_info.bytes_total,
        capacity.total_bytes(),
        capacity.total_bytes() as f64 / vertex_count.max(1) as f64,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}

fn main() {
    report_vec_bool_probe();
    report_topology_view_probe();
    report_boolean_case();
    report_continuum_winding_case();
    let benchmark_cases = cases();
    for case in &benchmark_cases {
        report_law_case(case);
        report_random_failure_case(case);
        report_degree_attack_case(case);
        report_residual_attack_case(case);
        report_directed_analysis_case(case);
        report_directed_process_case(case);
        report_kcore_case(case);
        report_bootstrap_case(case);
        report_incremental_case(case);
        report_achlioptas_case(case);
        report_invasion_case(case);
        for mode in ReferenceMode::ALL {
            report_component_case(case, mode);
        }
        for mode in ReferenceMode::ALL {
            report_adapter_case(case, mode);
        }
        for mode in [ReferenceMode::Site, ReferenceMode::Bond] {
            report_scan_case(case, mode);
        }
        for mode in [ReferenceMode::Site, ReferenceMode::Bond] {
            report_heterogeneous_scan_case(case, mode);
        }
        for mode in ReferenceMode::ALL {
            report_heterogeneous_adapter_case(case, mode);
        }
    }
    for case in benchmark_cases {
        for mode in ReferenceMode::ALL {
            report_case(&case, mode);
        }
    }
    for torus_case in torus_cases() {
        report_winding_case(&torus_case);
        report_winding_adapter_case(&torus_case);
    }
}

/// T1 probe: `DegreeAttackMC::sweep` + `measure` with the adaptive
/// residual-degree policy (seeded-key draw + lazy-heap order generation +
/// checkpoint walk) must be allocation-free in steady state, exactly like
/// the N3 original-degree probe above (same binsize discipline; the heap,
/// residual, active, and tie-key buffers are reserved at construction).
fn report_residual_attack_case(case: &Case) {
    let step = (case.lattice.vertex_count() / 16).max(1);
    let mut adapter = DegreeAttackMC::new(
        case.lattice.clone(),
        AttackPolicy::ResidualDegreeDescendingSeededKey,
        AttackGrid::Every { step },
    )
    .expect("probe adapter");
    let checkpoints = adapter.checkpoints().len();
    let mut context = Context::new_with_binsize(
        Xoshiro256PlusPlus::seed_from_u64(0x7e51),
        0,
        ADAPTER_WARMUP + ALLOCATION_SAMPLES + 1,
    );
    for _ in 0..ADAPTER_WARMUP {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
    }
    let info = measure(|| {
        for _ in 0..ALLOCATION_SAMPLES {
            adapter.sweep(&mut context);
            adapter.measure(&mut context);
            std::hint::black_box(&mut context);
        }
    });
    // Same observable set as the N3 probe: curve + coordinates +
    // susceptibility indicator + susceptibility array.
    let observable_count = 4;
    // Same policy as the kcore/bootstrap/incremental probes: a silent
    // regression to nonzero allocations must fail the probe.
    assert_eq!(
        info.count_total, 0,
        "residual-degree attack trajectory allocated {} times in steady state",
        info.count_total
    );
    eprintln!(
        "RESIDUAL_ATTACK_ALLOC topology={} V={} E={} checkpoints={} observables={} \
         measured_trajectories={} allocations_per_trajectory={:.3} \
         allocated_bytes_per_trajectory={:.3} peak_live_allocations={} peak_live_bytes={}",
        case.name,
        case.lattice.vertex_count(),
        case.lattice.edge_count(),
        checkpoints,
        observable_count,
        ALLOCATION_SAMPLES,
        info.count_total as f64 / ALLOCATION_SAMPLES as f64,
        info.bytes_total as f64 / ALLOCATION_SAMPLES as f64,
        info.count_max,
        info.bytes_max,
    );
}
