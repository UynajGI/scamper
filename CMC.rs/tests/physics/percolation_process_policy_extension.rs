//! A4: external process-policy extension (scoring, thresholds).
//!
//! Each policy here is implemented outside the library crate, so the erased
//! enum variants are not the only way to drive the algorithms. The built-in
//! equivalent of every custom policy must reproduce the reference order
//! exactly.

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    build_chain, build_square, kcore_decompose, AchlioptasGrid, AchlioptasProcessMC,
    AchlioptasRule, AchlioptasScore, AchlioptasTie, BootstrapDynamics, BootstrapError,
    BootstrapPercolationMC, BootstrapSeedRule, BootstrapSeeds, BootstrapThreshold,
    BootstrapThresholdRule, CandidateSampling, CoreThresholdRule, CoreThresholds, CsrLattice,
    EdgeUniverse, GraphView, InvasionError, InvasionMode, InvasionPercolationMC, InvasionStop,
    InvasionStopReason, InvasionStopRule, InvasionTrapping, InvasionTrappingRule,
    InvasionWorkspace, KCoreError, KCoreMC, KCoreWorkspace, TiePolicy, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Product rule reimplemented downstream.
struct CustomProduct;

impl AchlioptasScore for CustomProduct {
    fn score(&self, size_a: usize, size_b: usize) -> u128 {
        (size_a as u128) * (size_b as u128)
    }

    fn name(&self) -> &str {
        "custom-product"
    }
}

/// Lowest-dense-index tie break reimplemented downstream.
struct CustomLowestIndex;

impl AchlioptasTie for CustomLowestIndex {
    fn beats(&self, candidate: usize, score: u128, best: usize, best_score: u128) -> bool {
        score < best_score || (score == best_score && candidate < best)
    }

    fn name(&self) -> &str {
        "custom-lowest-index"
    }
}

/// A tie policy that always keeps the latest draw (used to prove that a
/// custom tie actually changes the reference order).
struct LatestDraw;

impl AchlioptasTie for LatestDraw {
    fn beats(&self, _candidate: usize, score: u128, _best: usize, best_score: u128) -> bool {
        score <= best_score
    }

    fn name(&self) -> &str {
        "latest-draw"
    }
}

fn build_process<Rule: AchlioptasScore, Tie: AchlioptasTie>(
    lattice: &CsrLattice,
    rule: Rule,
    tie: Tie,
) -> AchlioptasProcessMC<Rule, Tie> {
    AchlioptasProcessMC::new(
        lattice.clone(),
        rule,
        CandidateSampling::WithoutReplacement,
        2,
        tie,
        EdgeUniverse::AllEdges,
        AchlioptasGrid::Every { step: 4 },
    )
    .expect("valid process")
}

#[test]
fn custom_achlioptas_policies_reproduce_the_preset_reference_order() {
    let lattice = build_square(4, 4, false);
    let mut preset = build_process(
        &lattice,
        AchlioptasRule::Product,
        TiePolicy::LowestEdgeIndex,
    );
    let mut custom = build_process(&lattice, CustomProduct, CustomLowestIndex);

    let seed = 0x4134_1001;
    preset.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    custom.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    assert_eq!(custom.rule().name(), "custom-product");
    assert_eq!(custom.tie().name(), "custom-lowest-index");
    assert_eq!(custom.accepted_edges(), preset.accepted_edges());
    assert_eq!(custom.largest_curve(), preset.largest_curve());

    let mut latest = build_process(&lattice, AchlioptasRule::Product, LatestDraw);
    latest.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    assert_ne!(
        latest.accepted_edges(),
        preset.accepted_edges(),
        "a different tie policy must change the reference order"
    );
}

/// Per-vertex bootstrap thresholds implemented downstream.
struct CustomThresholds(Vec<u32>);

impl BootstrapThresholdRule for CustomThresholds {
    fn validate(&self, vertices: usize) -> Result<(), BootstrapError> {
        if self.0.len() == vertices {
            Ok(())
        } else {
            Err(BootstrapError::ThresholdLengthMismatch {
                expected: vertices,
                actual: self.0.len(),
            })
        }
    }

    fn fill(&self, thresholds: &mut [u32]) {
        thresholds.copy_from_slice(&self.0);
    }

    fn name(&self) -> &str {
        "custom-per-vertex"
    }
}

#[test]
fn custom_bootstrap_thresholds_reproduce_the_preset_trajectory() {
    let graph = build_chain(6, false);
    let values = vec![2, 1, 1, 2, 3, 0];
    let seeds = BootstrapSeeds::Fixed {
        vertices: vec![0, 5],
    };
    let mut preset = BootstrapPercolationMC::new(
        graph.clone(),
        BootstrapThreshold::Heterogeneous {
            values: values.clone(),
        },
        seeds.clone(),
        BootstrapDynamics::QueueClosure,
        8,
    )
    .expect("preset adapter");
    let mut custom = BootstrapPercolationMC::new(
        graph,
        CustomThresholds(values.clone()),
        seeds,
        BootstrapDynamics::QueueClosure,
        8,
    )
    .expect("custom adapter");
    assert_eq!(custom.threshold().name(), "custom-per-vertex");

    preset.realize(&mut Xoshiro256PlusPlus::seed_from_u64(9));
    custom.realize(&mut Xoshiro256PlusPlus::seed_from_u64(9));
    assert_eq!(custom.last_outcome(), preset.last_outcome());
    assert_eq!(custom.last_history(), preset.last_history());

    // The adapter rejects a wrong-length custom field exactly like the enum.
    let error = BootstrapPercolationMC::new(
        build_chain(4, false),
        CustomThresholds(vec![1]),
        BootstrapSeeds::Fixed { vertices: vec![0] },
        BootstrapDynamics::QueueClosure,
        4,
    )
    .err()
    .expect("short custom field is rejected");
    assert!(matches!(
        error,
        BootstrapError::ThresholdLengthMismatch {
            expected: 4,
            actual: 1
        }
    ));
}

/// Per-vertex core thresholds implemented downstream.
struct CustomCore(Vec<u32>);

impl CoreThresholdRule for CustomCore {
    fn validate(&self, vertices: usize) -> Result<(), KCoreError> {
        if self.0.len() == vertices {
            Ok(())
        } else {
            Err(KCoreError::ThresholdLength {
                expected: vertices,
                actual: self.0.len(),
            })
        }
    }

    fn fill(&self, thresholds: &mut Vec<u32>, vertices: usize) {
        let _ = vertices;
        thresholds.extend_from_slice(&self.0);
    }

    fn name(&self) -> &str {
        "custom-per-vertex"
    }
}

#[test]
fn custom_core_thresholds_reproduce_the_preset_decomposition() {
    let graph = build_square(3, 3, false);
    let values = vec![2, 2, 1, 2, 1, 1, 1, 2, 2];
    let active = vec![true; graph.vertex_count()];
    let mut preset_workspace = KCoreWorkspace::default();
    let mut custom_workspace = KCoreWorkspace::default();
    let preset = kcore_decompose(
        &graph,
        &CoreThresholds::PerVertex(values.clone()),
        &active,
        &mut preset_workspace,
    )
    .expect("preset");
    let custom = kcore_decompose(
        &graph,
        &CustomCore(values.clone()),
        &active,
        &mut custom_workspace,
    )
    .expect("custom");
    assert_eq!(custom, preset);
    assert_eq!(
        custom_workspace.deletion_order(),
        preset_workspace.deletion_order()
    );
    assert_eq!(
        custom_workspace.round_survivors(),
        preset_workspace.round_survivors()
    );

    let mut preset_adapter = KCoreMC::new(
        graph.clone(),
        1.0,
        CoreThresholds::PerVertex(values.clone()),
        4,
    )
    .expect("preset adapter");
    let mut custom_adapter =
        KCoreMC::new(graph, 1.0, CustomCore(values), 4).expect("custom adapter");
    assert_eq!(custom_adapter.thresholds().name(), "custom-per-vertex");
    preset_adapter.draw_and_peel(&mut Xoshiro256PlusPlus::seed_from_u64(11));
    custom_adapter.draw_and_peel(&mut Xoshiro256PlusPlus::seed_from_u64(11));
    assert_eq!(
        custom_adapter.round_survivors(),
        preset_adapter.round_survivors()
    );
    assert_eq!(custom_adapter.summary(), preset_adapter.summary());
}

#[test]
fn custom_core_threshold_rejects_a_wrong_length_before_drawing() {
    let graph = build_chain(4, false);
    let error = KCoreMC::new(graph, 1.0, CustomCore(vec![1; 3]), 2)
        .err()
        .expect("short field is rejected");
    assert!(matches!(
        error,
        KCoreError::ThresholdLength {
            expected: 4,
            actual: 3
        }
    ));
}

#[test]
fn policy_traits_do_not_change_the_built_in_draw_contract() {
    // A sanity check that the generic adapter still exposes the same
    // context-driven sweep path (and the new generic parameters default to
    // the built-in enums at use sites).
    let graph = build_chain(5, false);
    let mut adapter: BootstrapPercolationMC<_> = BootstrapPercolationMC::new(
        graph,
        BootstrapThreshold::Homogeneous { value: 1 },
        BootstrapSeeds::Fixed { vertices: vec![0] },
        BootstrapDynamics::SynchronousRounds,
        6,
    )
    .expect("adapter");
    let mut context = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(3), 0, 1);
    adapter.sweep(&mut context);
    adapter.measure(&mut context);
    assert_eq!(adapter.last_history().len(), 6);
}

/// Step budget implemented downstream.
#[derive(Clone, Copy)]
struct CustomStepBudget(usize);

impl InvasionStopRule for CustomStepBudget {
    fn after_seeding(&self, _workspace: &InvasionWorkspace) -> Option<InvasionStopReason> {
        (self.0 == 0).then_some(InvasionStopReason::StepBudget)
    }

    fn after_step(
        &self,
        _workspace: &InvasionWorkspace,
        step: usize,
    ) -> Option<InvasionStopReason> {
        (step == self.0).then_some(InvasionStopReason::StepBudget)
    }

    fn validate(&self, _vertices: usize) -> Result<(), InvasionError> {
        Ok(())
    }

    fn name(&self) -> &str {
        "custom-step-budget"
    }
}

/// Non-trapping behaviour implemented downstream.
#[derive(Clone, Copy)]
struct CustomFree;

impl InvasionTrappingRule for CustomFree {
    fn validate(&self, _vertices: usize) -> Result<(), InvasionError> {
        Ok(())
    }

    fn name(&self) -> &str {
        "custom-free"
    }
}

/// Frontier filter that removes one selected vertex from the growth.
#[derive(Clone, Copy)]
struct BlockVertex(usize);

impl InvasionTrappingRule for BlockVertex {
    fn filters_frontier(&self) -> bool {
        true
    }

    fn invadable<G: UndirectedGraphView>(
        &self,
        _graph: &G,
        _workspace: &mut InvasionWorkspace,
        vertex: usize,
    ) -> bool {
        vertex != self.0
    }

    fn validate(&self, _vertices: usize) -> Result<(), InvasionError> {
        Ok(())
    }

    fn name(&self) -> &str {
        "block-vertex"
    }
}

#[test]
fn custom_invasion_policies_reproduce_the_preset_outcome() {
    let graph = build_chain(6, false);
    let seed = 0x4134_2001;
    let mut preset = InvasionPercolationMC::new(
        graph.clone(),
        InvasionMode::Site,
        vec![0],
        InvasionStop::Steps(3),
        InvasionTrapping::NonTrapping,
        6,
    )
    .expect("preset adapter");
    let mut custom = InvasionPercolationMC::new(
        graph.clone(),
        InvasionMode::Site,
        vec![0],
        CustomStepBudget(3),
        CustomFree,
        6,
    )
    .expect("custom adapter");
    assert_eq!(custom.stop().name(), "custom-step-budget");
    assert_eq!(custom.trapping().name(), "custom-free");

    preset.realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    custom.realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    assert_eq!(custom.last_outcome(), preset.last_outcome());
    assert_eq!(custom.last_order(), preset.last_order());
    assert_eq!(custom.last_history(), preset.last_history());
}

#[test]
fn custom_frontier_filter_changes_the_invasion_outcome() {
    let graph = build_chain(3, false);
    let mut filtered = InvasionPercolationMC::new(
        graph.clone(),
        InvasionMode::Site,
        vec![0],
        InvasionStop::GrowAll,
        BlockVertex(2),
        4,
    )
    .expect("filtered adapter");
    let mut free = InvasionPercolationMC::new(
        graph,
        InvasionMode::Site,
        vec![0],
        InvasionStop::GrowAll,
        InvasionTrapping::NonTrapping,
        4,
    )
    .expect("free adapter");
    let seed = 0x4134_2002;
    filtered.realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    free.realize(&mut Xoshiro256PlusPlus::seed_from_u64(seed));
    assert_eq!(filtered.last_outcome().trapped_sites, 1);
    assert_eq!(filtered.last_order(), &[0, 1]);
    assert_eq!(free.last_order(), &[0, 1, 2]);
}

/// Every-even-vertex seed law implemented downstream.
#[derive(Clone, Copy)]
struct EveryOther;

impl BootstrapSeedRule for EveryOther {
    fn validate(&self, _vertices: usize) -> Result<(), BootstrapError> {
        Ok(())
    }

    fn fill<R: rand::Rng + ?Sized>(&mut self, active: &mut [u8], _rng: &mut R) {
        for (index, slot) in active.iter_mut().enumerate() {
            *slot = u8::from(index % 2 == 0);
        }
    }

    fn name(&self) -> &str {
        "every-other"
    }
}

#[test]
fn custom_seed_rule_reproduces_the_preset_closure() {
    let graph = build_chain(6, false);
    let mut preset = BootstrapPercolationMC::new(
        graph.clone(),
        BootstrapThreshold::Homogeneous { value: 1 },
        BootstrapSeeds::Fixed {
            vertices: vec![0, 2, 4],
        },
        BootstrapDynamics::QueueClosure,
        7,
    )
    .expect("preset adapter");
    let mut custom = BootstrapPercolationMC::new(
        graph,
        BootstrapThreshold::Homogeneous { value: 1 },
        EveryOther,
        BootstrapDynamics::QueueClosure,
        7,
    )
    .expect("custom adapter");
    assert_eq!(custom.seeds().name(), "every-other");
    preset.realize(&mut Xoshiro256PlusPlus::seed_from_u64(13));
    custom.realize(&mut Xoshiro256PlusPlus::seed_from_u64(13));
    assert_eq!(custom.last_outcome(), preset.last_outcome());
    assert_eq!(custom.last_history(), preset.last_history());
}
