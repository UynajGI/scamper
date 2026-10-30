//! A6: explicit specs, parameter round-trip, and provenance.
//!
//! The spec path must reproduce the parameter preset bit for bit (same seed,
//! same draws), reach the non-schema laws through `Run::from_parts`, and
//! stamp a versioned provenance object that survives the Results JSON round
//! trip.

use carlo_rs::{
    Context, FromParams, Params, RayonBackend, Results, RngPhase, RngStreamKey, Run, RunConfig,
    RunId, Scheduler, TaskId,
};
use cmc_rs::{
    BoundarySets, EstimatorSpec, GraphView, LatticeKind, PercolationSpec, Probability, ProcessSpec,
    StaticLaw, StaticObservable, StaticPercolationMC, TopologySpec, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use serde_json::{json, Value};

fn params(entries: &[(&str, &str)]) -> Params {
    let mut params = Params::new();
    for (key, value) in entries {
        params.set(key, *value);
    }
    params
}

fn replay_spec(spec: &PercolationSpec, sweeps: u64, seed: u64) -> Results {
    let adapter = spec.build().expect("valid spec");
    let mut run = Run::from_parts(
        Context::new_with_binsize(
            RngStreamKey::new(seed)
                .with_phase(RngPhase::Initialization)
                .seeded(),
            0,
            1,
        ),
        adapter,
        TaskId::new(0),
        RunId::new(0),
        RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: sweeps,
            binsize: 1,
            base_seed: seed,
            ..Default::default()
        },
    );
    run.run(sweeps);
    run.finalize(seed)
}

fn assert_same_results(expected: &Results, actual: &Results) {
    let names = |results: &Results| {
        results
            .estimates()
            .keys()
            .filter(|name| !name.starts_with("_ll_"))
            .cloned()
            .collect::<std::collections::BTreeSet<_>>()
    };
    assert_eq!(names(expected), names(actual), "estimate name sets differ");
    for (name, left) in expected.estimates() {
        let right = actual.get(name).unwrap_or_else(|| panic!("missing {name}"));
        assert_eq!(left.mean, right.mean, "mean of {name}");
        assert_eq!(left.stderr, right.stderr, "stderr of {name}");
        assert_eq!(left.n_bins, right.n_bins, "bins of {name}");
    }
    assert_eq!(
        expected.array_estimates().len(),
        actual.array_estimates().len()
    );
    for (name, left) in expected.array_estimates() {
        let right = actual
            .get_array(name)
            .unwrap_or_else(|| panic!("missing array {name}"));
        assert_eq!(left.shape, right.shape, "shape of {name}");
        for (index, (left, right)) in left.elements.iter().zip(&right.elements).enumerate() {
            assert_eq!(left.mean, right.mean, "{name}[{index}] mean");
            assert_eq!(left.stderr, right.stderr, "{name}[{index}] stderr");
        }
    }
}

#[test]
fn spec_replays_the_parameter_preset_bitwise() {
    let cases = [
        params(&[
            ("lattice_type", "square"),
            ("Lx", "4"),
            ("Ly", "4"),
            ("p", "0.5"),
        ]),
        params(&[
            ("lattice_type", "square"),
            ("Lx", "3"),
            ("Ly", "3"),
            ("pbc", "true"),
            ("p", "0.5"),
            ("observables", "wrapping,winding-rank"),
        ]),
    ];
    for case in cases {
        let seed = 0x4136_0001;
        let sweeps = 24;
        let expected = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 0,
                measurement_sweeps: sweeps,
                binsize: 1,
                base_seed: seed,
                ..Default::default()
            },
        )
        .run_one::<StaticPercolationMC>(&case);
        let spec = PercolationSpec::from_params(&case).expect("spec");
        let actual = replay_spec(&spec, sweeps, seed);
        assert_same_results(&expected, &actual);
    }
}

#[test]
fn stamped_provenance_survives_the_results_json_round_trip() {
    let case = params(&[
        ("lattice_type", "square"),
        ("Lx", "3"),
        ("Ly", "3"),
        ("pbc", "true"),
        ("p", "0.5"),
        ("observables", "wrapping,winding-rank"),
    ]);
    let spec = PercolationSpec::from_params(&case).expect("spec");
    let mut results = replay_spec(&spec, 8, 0x4136_0002);
    spec.stamp(&mut results);
    let provenance = spec.provenance().to_json();
    assert_eq!(
        results.metadata().attribute("percolation"),
        Some(&provenance)
    );

    let serialized: Value =
        serde_json::from_str(&results.to_json().expect("json")).expect("valid json");
    assert_eq!(
        serialized["metadata"]["attributes"]["percolation"],
        provenance
    );
    assert_eq!(
        serialized["metadata"]["attributes"]["percolation"],
        json!({
            "schema": "cmc.percolation.provenance",
            "version": 1,
            "topology": "square/3x3/pbc",
            "law": "site/p=0.5",
            "process": "independent-occupation/v1",
            "observables": ["Wrapping", "WindingRank"],
            "boundary_queries": 1,
            "estimator": "sample-moments/v1",
        })
    );
}

#[test]
fn non_schema_law_runs_through_from_parts_and_reports_its_spec() {
    let topology = TopologySpec::new(LatticeKind::Square { lx: 2, ly: 2 }, false);
    let law = StaticLaw::site_heterogeneous(
        [0.1, 0.9, 0.4, 0.65].map(|value| Probability::new(value).expect("p")),
    );
    let spec = PercolationSpec {
        topology,
        law,
        process: ProcessSpec::IndependentOccupation,
        observables: vec![StaticObservable::LargestSize],
        boundaries: Some(BoundarySets::new(vec![0], vec![1])),
        estimator: EstimatorSpec::SampleMoments,
    };
    assert_eq!(spec.provenance().law(), "site-heterogeneous/vertices=4");
    assert!(matches!(
        spec.to_params(),
        Err(cmc_rs::SpecError::LawNotParamsExpressible)
    ));

    let mut results = replay_spec(&spec, 8, 0x4136_0003);
    assert!(results.get("LargestSize").is_some());
    spec.stamp(&mut results);
    assert_eq!(
        results
            .metadata()
            .attribute("percolation")
            .expect("stamped")["law"],
        "site-heterogeneous/vertices=4"
    );
}

#[test]
fn spec_build_exposes_the_same_topology_as_the_preset() {
    let case = params(&[
        ("lattice_type", "square"),
        ("Lx", "5"),
        ("Ly", "3"),
        ("mode", "bond"),
        ("p", "0.4"),
    ]);
    let spec = PercolationSpec::from_params(&case).expect("spec");
    let built = spec.build().expect("adapter");
    let preset = StaticPercolationMC::from_params(&case, &mut Xoshiro256PlusPlus::seed_from_u64(1))
        .expect("preset");
    assert_eq!(built.lattice().vertex_count(), 15);
    assert_eq!(
        built.lattice().vertex_count(),
        preset.lattice().vertex_count()
    );
    assert_eq!(built.lattice().edge_count(), preset.lattice().edge_count());
    assert_eq!(built.law(), preset.law());
}
