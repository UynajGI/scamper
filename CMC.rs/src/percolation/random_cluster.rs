//! q = 1 random-cluster reduction (roadmap §15 R2).
//!
//! Exact statement validated by this stage, configuration by
//! configuration: **the q = 1 FK/random-cluster measure on a graph `G`
//! equals the product Bernoulli(p) measure on its bonds.**
//!
//! The random-cluster measure with edge parameter `p` and cluster weight
//! `q` on physical edges `E` is
//!
//! ```text
//! phi_{p,q}(A) = p^{|A|} (1-p)^{E-|A|} q^{k(A)} / Z(p, q),   A subset of E,
//! Z(p, q)      = sum_A p^{|A|} (1-p)^{E-|A|} q^{k(A)},
//! ```
//!
//! where `k(A)` is the number of connected components of `(V, A)` with
//! every vertex active — the bond-percolation convention of F3. At
//! `q = 1` the cluster term is exactly `1^{k(A)} = 1` for every
//! configuration, so it drops out of the weight entirely and
//!
//! ```text
//! Z(p, 1)   = sum_A p^{|A|} (1-p)^{E-|A|} = (p + (1-p))^E = 1,
//! phi_{p,1} = product over e of [ p   if e in A,
//!                                 (1-p) if e not in A ],
//! ```
//!
//! which is precisely the independent bond Bernoulli product measure
//! (with heterogeneous per-edge `p_e`, the product `prod_e p_e /
//! (1-p_e)` analogously). This is exact mathematics, not an
//! approximation: the collapse `q = 1 -> independent bonds` is the
//! founding observation of Fortuin and Kasteleyn (Physica 57, 536
//! (1972); pre-arXiv, so no PDF anchor is archived — the identity is
//! proved, not fitted). The random-cluster model itself originates
//! there.
//!
//! Because the `q = 1` bond draw *is* the Bernoulli product draw, this
//! adapter contains **no duplicated sampling law**: every realization is
//! drawn by delegating to the F2 [`BondBernoulli`] law (bitwise-identical
//! activity to `StaticLaw::bond` under the same RNG stream), and the
//! analysis runs through the R1 decoupled contract — the draw is copied
//! into a reused dense mask and analyzed via [`analyze_bonds`]
//! ([`ExternalBondConfiguration`]), a different activity-plumbing path
//! than `StaticPercolationMC`'s `StaticConfiguration`. Each realization is
//! exposed as R1's owned [`BondSample`] record through
//! [`RandomClusterMC::realization`]; that surface — Bernoulli delegate in,
//! `BondSample` + F3 `ObservablePlan` scalars out — is the contract R3
//! (q = 2 Edwards-Sokal/Swendsen-Wang) replaces the delegate behind while
//! keeping everything downstream unchanged.
//!
//! `q != 1` is a typed rejection ([`RandomClusterError::UnsupportedQ`]):
//! generic real-`q` samplers (Sweeny, Chayes-Machta) need a separate plan
//! per the roadmap, after R3.

use core::fmt;

use carlo_rs::{Context, MonteCarlo};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{
    analyze_bonds, BondBernoulli, BondSample, EdgeActivity, ObservablePlan, ObservablePlanError,
    Probability, ProbabilityError, ProbabilityField, SamplingError, StaticConfiguration,
};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Cluster-weight parameter `q` of the random-cluster measure.
///
/// R2 implements only `q = 1` (the exact Bernoulli collapse above), so the
/// value is pinned to `1.0` at construction: [`ClusterQ::q1`] builds the
/// one supported value directly and [`ClusterQ::new`] accepts exactly `1.0`,
/// returning the typed [`RandomClusterError::UnsupportedQ`] (which cites the
/// roadmap's separate-plan requirement) for anything else — including `NaN`.
/// The parameter is carried explicitly rather than hardcoded at each use so
/// R3 (q = 2) and any future real-`q` stage widen the domain in one place.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClusterQ {
    value: f64,
}

impl ClusterQ {
    /// The one supported cluster weight, `q = 1`.
    pub const fn q1() -> Self {
        Self { value: 1.0 }
    }

    /// Accept exactly `q = 1`; reject every other value (including `NaN`)
    /// with the roadmap-citing typed error.
    pub fn new(value: f64) -> Result<Self, RandomClusterError> {
        if value == 1.0 {
            Ok(Self { value })
        } else {
            Err(RandomClusterError::UnsupportedQ { value })
        }
    }

    /// The cluster weight (always `1.0` in R2).
    #[inline]
    pub const fn value(self) -> f64 {
        self.value
    }
}

impl fmt::Display for ClusterQ {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "q = {}", self.value)
    }
}

/// Owned scheduler adapter for the q = 1 random-cluster reduction.
///
/// One `sweep()` draws one complete bond realization by delegating to the
/// F2 [`BondBernoulli`] law (uniform `p`, or a validated per-edge
/// probability field through [`Self::new_heterogeneous`]) and mirrors it
/// into a reused dense mask; the following `measure()` analyzes that mask
/// through the R1 external-bond path and records the F3
/// [`ObservablePlan`] scalars. Samples are i.i.d. — set
/// `thermalization_sweeps = 0`. After `sweep()`, [`Self::mask`],
/// [`Self::configuration`], and [`Self::realization`] describe one and
/// the same realization.
pub struct RandomClusterMC {
    lattice: CsrLattice,
    cluster_q: ClusterQ,
    law: BondBernoulli<'static>,
    configuration: StaticConfiguration,
    mask: Vec<bool>,
    workspace: super::ComponentWorkspace,
    plan: ObservablePlan,
}

impl RandomClusterMC {
    /// Uniform-`p` q = 1 random-cluster adapter.
    ///
    /// Rejections mirror the static adapter's: invalid topology, winding
    /// observables (this adapter carries no embedding; winding an external
    /// bond configuration goes through the core `analyze_winding` with an
    /// `ExternalBondConfiguration` directly), moment plans that cannot
    /// convert exactly to `f64`, boundary-query domain mismatches, and an
    /// invalid `p`.
    pub fn new(
        lattice: CsrLattice,
        cluster_q: ClusterQ,
        bond_probability: f64,
        plan: ObservablePlan,
    ) -> Result<Self, RandomClusterError> {
        let probability = Probability::new(bond_probability)
            .map_err(|source| RandomClusterError::Probability("p", source))?;
        Self::from_field(
            lattice,
            cluster_q,
            ProbabilityField::Uniform(probability),
            plan,
        )
    }

    /// Heterogeneous q = 1 random-cluster adapter: one validated
    /// probability per physical edge (dense `EdgeId` order), exactly the
    /// F2 heterogeneous bond field. At `q = 1` the FK measure with
    /// per-edge `p_e` is the product `prod_e Bernoulli(p_e)`.
    pub fn new_heterogeneous(
        lattice: CsrLattice,
        cluster_q: ClusterQ,
        bond_probabilities: impl Into<Vec<Probability>>,
        plan: ObservablePlan,
    ) -> Result<Self, RandomClusterError> {
        Self::from_field(
            lattice,
            cluster_q,
            ProbabilityField::Owned(bond_probabilities.into()),
            plan,
        )
    }

    fn from_field(
        lattice: CsrLattice,
        cluster_q: ClusterQ,
        edges: ProbabilityField<'static>,
        plan: ObservablePlan,
    ) -> Result<Self, RandomClusterError> {
        lattice
            .validate()
            .map_err(|error| RandomClusterError::Topology(error.to_string()))?;
        if plan.requests_winding() {
            return Err(RandomClusterError::WindingUnsupported);
        }
        plan.validate(lattice.vertex_count())
            .map_err(RandomClusterError::Observable)?;
        for (index, query) in plan.boundary_queries().iter().enumerate() {
            if query.vertex_count() != lattice.vertex_count() {
                return Err(RandomClusterError::QueryDomain {
                    query: index,
                    expected: lattice.vertex_count(),
                    actual: query.vertex_count(),
                });
            }
        }
        // Uniform fields carry no length; heterogeneous fields must match
        // the topology's physical edge count exactly (an edgeless graph
        // with a uniform law is legal — `None != Some(0)` must not reject).
        if edges
            .len()
            .is_some_and(|actual| actual != lattice.edge_count())
        {
            return Err(RandomClusterError::LawDomain(
                SamplingError::ProbabilityFieldLength {
                    entity: "edge",
                    expected: lattice.edge_count(),
                    actual: edges.len().unwrap_or_default(),
                },
            ));
        }
        let configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
        let mask = vec![false; lattice.edge_count()];
        let mut workspace = super::ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| RandomClusterError::Capacity)?;
        Ok(Self {
            lattice,
            cluster_q,
            law: BondBernoulli::new(edges),
            configuration,
            mask,
            workspace,
            plan,
        })
    }

    /// The owned lattice under analysis.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The cluster-weight parameter (always `q = 1` in R2).
    pub const fn cluster_q(&self) -> ClusterQ {
        self.cluster_q
    }

    /// The F2 bond law every realization delegates to.
    pub const fn bond_law(&self) -> &BondBernoulli<'static> {
        &self.law
    }

    /// The observable plan recorded by `measure()`.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// The last drawn configuration (for reduction tests).
    pub const fn configuration(&self) -> &StaticConfiguration {
        &self.configuration
    }

    /// The dense mask view of the last drawn realization: one `bool` per
    /// physical edge in dense `EdgeId` order, bitwise-equal to the
    /// [`BondBernoulli`] draw and to the mask read back from
    /// [`Self::configuration`].
    #[inline]
    pub fn mask(&self) -> &[bool] {
        &self.mask
    }

    /// The last drawn realization as R1's owned open-edge record.
    ///
    /// This is the R3 contract surface: an Edwards-Sokal/Swendsen-Wang
    /// adapter exposes its actual FK bond draw as the same `BondSample`.
    /// Materialized on demand (`O(E)` plus one allocation) so the
    /// scheduler's sweep/measure path stays allocation-free.
    pub fn realization(&self) -> BondSample {
        BondSample::from_mask(&self.lattice, &self.mask)
            .expect("mask length pinned to the topology at construction")
    }

    /// Mirror the configuration's bond draw into the dense mask view.
    fn sync_mask(&mut self) {
        for (index, open) in self.mask.iter_mut().enumerate() {
            let edge = self.lattice.edge_id(index).expect("dense edge index");
            *open = self.configuration.edge_active(edge);
        }
    }
}

impl MonteCarlo for RandomClusterMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.law
            .sample(&self.lattice, &mut self.configuration, &mut ctx.rng)
            .expect("constructor validated the q = 1 law domain against the topology");
        self.sync_mask();
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let result = analyze_bonds(
            &self.lattice,
            &self.mask,
            self.plan.boundary_queries(),
            &mut self.workspace,
        )
        .expect("constructor prepared matching analysis domains");
        self.plan.record(
            self.lattice.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            None,
            ctx,
        );
    }

    fn name(&self) -> &'static str {
        "RandomClusterQ1"
    }
}

/// Validation failure while constructing a q = 1 random-cluster adapter.
#[derive(Debug)]
pub enum RandomClusterError {
    /// `q != 1`. At `q = 1` the FK measure equals independent bond
    /// Bernoulli (validated configuration by configuration in
    /// `tests/physics/random_cluster_q1.rs`); generic real-`q`
    /// random-cluster samplers (Sweeny, Chayes-Machta) need a separate
    /// plan after the R3 `q = 2` stage, per roadmap §15.
    UnsupportedQ {
        value: f64,
    },
    Probability(&'static str, ProbabilityError),
    Topology(String),
    LawDomain(SamplingError),
    Observable(ObservablePlanError),
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// The plan requests winding observables but this adapter carries no
    /// periodic embedding. Winding an external bond configuration is a
    /// core `analyze_winding` call, not an adapter feature.
    WindingUnsupported,
    Capacity,
}

impl fmt::Display for RandomClusterError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedQ { value } => write!(
                formatter,
                "random-cluster q must be exactly 1 in R2, got {value}: at q = 1 the FK \
                 measure is exactly independent bond Bernoulli, and generic real-q \
                 samplers (Sweeny, Chayes-Machta) need a separate plan (roadmap section \
                 15, after the R3 q = 2 stage)"
            ),
            Self::Probability(name, source) => write!(formatter, "invalid `{name}`: {source}"),
            Self::Topology(reason) => write!(formatter, "invalid topology: {reason}"),
            Self::LawDomain(source) => write!(formatter, "bond-law domain mismatch: {source}"),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} has V={actual}; adapter topology has V={expected}"
            ),
            Self::WindingUnsupported => write!(
                formatter,
                "winding observables are not part of the R2 adapter; call the core \
                 analyze_winding with an ExternalBondConfiguration instead"
            ),
            Self::Capacity => formatter.write_str("could not reserve component workspace"),
        }
    }
}

impl std::error::Error for RandomClusterError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, build_square, BoundaryQuery, StaticObservable, VertexActivity};
    use rand::SeedableRng;

    fn probability(value: f64) -> Probability {
        Probability::new(value).expect("test probability")
    }

    fn plan(lattice: &CsrLattice) -> ObservablePlan {
        let query = BoundaryQuery::new(lattice, &[0], &[lattice.vertex_count() - 1])
            .expect("test boundary query");
        ObservablePlan::all(vec![query])
    }

    fn context(seed: u64) -> Context<Xoshiro256PlusPlus> {
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1)
    }

    #[test]
    fn cluster_q_is_pinned_to_one_with_a_roadmap_cited_rejection() {
        assert_eq!(ClusterQ::q1().value(), 1.0);
        assert!(matches!(ClusterQ::new(1.0), Ok(q) if q == ClusterQ::q1()));
        for value in [0.0, 0.5, 2.0, -1.0, f64::NAN, f64::INFINITY] {
            let error = ClusterQ::new(value)
                .err()
                .unwrap_or_else(|| panic!("accepted q = {value}"));
            assert!(
                matches!(error, RandomClusterError::UnsupportedQ { value: v } if v.to_bits() == value.to_bits()),
                "wrong error for q = {value}: {error}"
            );
            let message = error.to_string();
            assert!(message.contains("exactly 1"), "{message}");
            assert!(message.contains("roadmap"), "{message}");
            assert!(message.contains("separate plan"), "{message}");
        }
    }

    #[test]
    fn constructor_rejections_mirror_the_static_adapter_vocabulary() {
        let lattice = build_square(3, 3, false);
        let q1 = ClusterQ::q1();
        let rejected = |result: Result<RandomClusterMC, RandomClusterError>| {
            result.err().expect("adapter must be rejected")
        };

        // Invalid p.
        assert!(matches!(
            rejected(RandomClusterMC::new(
                lattice.clone(),
                q1,
                1.5,
                plan(&lattice)
            )),
            RandomClusterError::Probability("p", _)
        ));
        // Winding observables need the core winding analyzer, not this adapter.
        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        assert!(matches!(
            rejected(RandomClusterMC::new(lattice.clone(), q1, 0.5, winding)),
            RandomClusterError::WindingUnsupported
        ));
        // Boundary query built for a different vertex count.
        let foreign = build_chain(5, false);
        let foreign_query =
            BoundaryQuery::new(&foreign, &[0], &[4]).expect("valid on the chain only");
        let foreign_plan =
            ObservablePlan::new([StaticObservable::LargestSize], vec![foreign_query]);
        assert!(matches!(
            rejected(RandomClusterMC::new(lattice.clone(), q1, 0.5, foreign_plan)),
            RandomClusterError::QueryDomain {
                query: 0,
                expected: 9,
                actual: 5
            }
        ));
        // Heterogeneous field length against the topology's edge count.
        let short = (0..5).map(|_| probability(0.5)).collect::<Vec<_>>();
        assert!(matches!(
            rejected(RandomClusterMC::new_heterogeneous(
                lattice.clone(),
                q1,
                short,
                plan(&lattice)
            )),
            RandomClusterError::LawDomain(SamplingError::ProbabilityFieldLength {
                entity: "edge",
                expected: 12,
                actual: 5
            })
        ));
        // A matching heterogeneous field constructs.
        let values = vec![probability(0.4); 12];
        let matching_plan = plan(&lattice);
        assert!(RandomClusterMC::new_heterogeneous(lattice, q1, values, matching_plan).is_ok());
    }

    #[test]
    fn sweep_leaves_mask_configuration_and_realization_bitwise_consistent() {
        let lattice = build_square(4, 5, false);
        let mut adapter = RandomClusterMC::new(
            lattice.clone(),
            ClusterQ::q1(),
            0.5,
            ObservablePlan::new([StaticObservable::GiantFraction], Vec::new()),
        )
        .expect("valid adapter");
        for seed in [1u64, 77, 2026] {
            adapter.sweep(&mut context(seed));
            let configuration_mask = (0..lattice.edge_count())
                .map(|index| {
                    adapter
                        .configuration()
                        .edge_active(lattice.edge_id(index).expect("dense edge index"))
                })
                .collect::<Vec<_>>();
            assert_eq!(adapter.mask(), configuration_mask.as_slice());
            assert_eq!(adapter.realization().to_mask(), configuration_mask);
            assert_eq!(
                adapter.realization().open_count(),
                configuration_mask.iter().filter(|&&open| open).count()
            );
            for vertex in lattice.vertex_ids() {
                assert!(adapter.configuration().vertex_active(vertex));
            }
        }
    }

    #[test]
    fn measure_records_endpoint_exact_scalars_through_the_external_path() {
        let lattice = build_chain(6, false);
        let mut adapter = RandomClusterMC::new(
            lattice.clone(),
            ClusterQ::q1(),
            1.0,
            ObservablePlan::new(
                [
                    StaticObservable::ComponentCount,
                    StaticObservable::GiantFraction,
                    StaticObservable::FiniteClusterSusceptibility,
                ],
                Vec::new(),
            ),
        )
        .expect("valid adapter");
        let mut ctx = context(9);
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        let values = ctx.finalize_measurements();
        assert_eq!(values["ComponentCount"].mean, 1.0);
        assert_eq!(values["GiantFraction"].mean, 1.0);
        assert_eq!(values["FiniteClusterSusceptibilityDefined"].mean, 0.0);

        let mut empty = RandomClusterMC::new(
            lattice,
            ClusterQ::q1(),
            0.0,
            ObservablePlan::new([StaticObservable::ComponentCount], Vec::new()),
        )
        .expect("valid adapter");
        let mut ctx = context(10);
        empty.sweep(&mut ctx);
        empty.measure(&mut ctx);
        assert_eq!(ctx.finalize_measurements()["ComponentCount"].mean, 6.0);
    }

    #[test]
    fn zero_edge_graph_samples_and_measures() {
        let lattice = crate::CsrLattice::from_edges(3, Vec::new());
        let mut adapter = RandomClusterMC::new(
            lattice,
            ClusterQ::q1(),
            0.5,
            ObservablePlan::new([StaticObservable::ComponentCount], Vec::new()),
        )
        .expect("edgeless graph is legal");
        let mut ctx = context(11);
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        assert!(adapter.mask().is_empty());
        assert_eq!(adapter.realization().open_count(), 0);
        assert_eq!(ctx.finalize_measurements()["ComponentCount"].mean, 3.0);
    }
}
