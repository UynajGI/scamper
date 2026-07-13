//! Frozen energy-bias Metropolis kernels for umbrella and multicanonical production.

use crate::algorithms::{Algorithm, SimulationPhase};
use crate::audit::{audit_lattice_cache, audit_macrostate_bin, should_audit_cache};
use crate::core::cache::EnergyPatch;
use crate::core::r#move::SiteSpinMove;
use crate::core::trial::TrialEvaluator;
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::generalized::{Histogram, LogBias, MacrostateAxis};
use crate::lattice::interaction::{Hamiltonian, Proposable};
use crate::lattice::proposal::{ProposalStrategy, StandardStrategy};
use crate::lattice::state::System;
use carlo_rs::accept_log_probability;
use rand::Rng;

/// Local Metropolis-Hastings kernel driven by a frozen total energy log weight.
#[derive(Debug, Clone)]
pub struct EnergyBiasCore<A, B, S = StandardStrategy> {
    axis: A,
    bias: B,
    strategy: S,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    histogram: Histogram,
    patch: EnergyPatch,
    energy_check_interval: u64,
    sweeps: u64,
    out_of_range_proposals: u64,
    last_visited_bin: Option<usize>,
}

impl<A, B> EnergyBiasCore<A, B, StandardStrategy>
where
    A: MacrostateAxis,
    B: LogBias,
{
    pub fn new(axis: A, bias: B) -> Self {
        Self::with_strategy(axis, bias, StandardStrategy::new())
    }
}

impl<A, B, S> EnergyBiasCore<A, B, S>
where
    A: MacrostateAxis,
    B: LogBias,
{
    pub fn with_strategy(axis: A, bias: B, strategy: S) -> Self {
        assert_eq!(
            axis.bins(),
            bias.bins(),
            "energy axis and bias must have the same bin count"
        );
        let histogram = Histogram::new(axis.bins()).expect("validated axis has at least one bin");
        Self {
            axis,
            bias,
            strategy,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            histogram,
            patch: EnergyPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
            out_of_range_proposals: 0,
            last_visited_bin: None,
        }
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }

    #[inline]
    pub const fn axis(&self) -> &A {
        &self.axis
    }

    #[inline]
    pub const fn bias(&self) -> &B {
        &self.bias
    }

    #[inline]
    pub const fn histogram(&self) -> &Histogram {
        &self.histogram
    }

    pub fn clear_histogram(&mut self) {
        self.histogram.clear();
    }

    /// Number of trial energies rejected because they were outside the fixed axis.
    #[inline]
    pub const fn out_of_range_proposals(&self) -> u64 {
        self.out_of_range_proposals
    }
}

impl<H, A, B, S> Algorithm<H> for EnergyBiasCore<A, B, S>
where
    H: Hamiltonian + Proposable,
    A: MacrostateAxis,
    B: LogBias,
    S: ProposalStrategy<H>,
{
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &H,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        let sites = self
            .order
            .prepare(system.n_sites(), self.visit_schedule, rng);
        for &site in sites {
            let old_bin = self
                .axis
                .bin(system.energy)
                .expect("accepted energy lies outside the generalized axis");
            let proposed_spin = self.strategy.propose(model, system, site, rng);
            let movement = SiteSpinMove::new(site, proposed_spin.spin);
            let delta = system.evaluate_trial(model, &movement, &mut self.patch);
            let new_bin = self.axis.bin(system.energy + delta.energy);
            if new_bin.is_none() {
                self.out_of_range_proposals = self.out_of_range_proposals.saturating_add(1);
            }
            let log_acceptance = new_bin.map_or(f64::NEG_INFINITY, |bin| {
                self.bias.log_weight_ratio(old_bin, bin) + proposed_spin.log_reverse_over_forward
            });
            let accepted = accept_log_probability(log_acceptance, rng);
            if accepted {
                <System as TrialEvaluator<H, SiteSpinMove>>::commit_trial(
                    system,
                    &movement,
                    &self.patch,
                );
            }
            self.strategy.record_result(accepted);
            let visited_bin = if accepted {
                new_bin.expect("accepted generalized trial has a bin")
            } else {
                old_bin
            };
            self.histogram.record(visited_bin);
            self.last_visited_bin = Some(visited_bin);
        }
        self.strategy.finish_sweep(phase.allows_adaptation());

        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_lattice_cache(system, model).expect("energy-bias cache audit failed");
            assert_eq!(
                self.histogram.bins(),
                self.axis.bins(),
                "energy-bias histogram/axis bin mismatch"
            );
            if let Some(bin) = self.last_visited_bin {
                audit_macrostate_bin(&self.axis, system.energy, bin)
                    .expect("energy-bias macrostate cache audit failed");
            }
        }
    }

    fn name(&self) -> &'static str {
        "Frozen generalized-energy Metropolis-Hastings"
    }
}
