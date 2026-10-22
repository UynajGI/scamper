//! Owned facade adapter for random network failure (N2).
//!
//! Random failure removes each vertex independently with probability
//! `1 - keep`. The sampler is exactly the F2 uniform site Bernoulli law —
//! this adapter duplicates no sampling code — and every degree-conditioned
//! observable is conditioned on the **original** graph degree cached in the
//! [`NetworkProfile`], never on residual degree.
//!
//! Scalar component observables reuse the F4 [`ObservablePlan`] verbatim.
//! Degree-conditioned observables are recorded per sample as C0 element-
//! preserving arrays indexed by original degree `k = 0 ..= max_degree`:
//!
//! - `DegreeVertexCount[k]`: `n_k`, vertices with original degree `k`
//!   (constant across samples; makes empty bins machine-readable).
//! - `DegreeOccupied[k]`: 1 when `n_k > 0`, else 0 (constant).
//! - `ActiveFractionByDegree[k]`: active vertices of degree `k` over `n_k`.
//! - `GiantMembershipByDegree[k]`: vertices of degree `k` inside the
//!   canonical largest component over `n_k` (F4 tie policy: equal sizes
//!   resolve to the lowest component identity).
//! - `DegreeFirstMoment[k]`: `sum_{active v, deg k} s(v) / n_k`.
//! - `DegreeSecondMoment[k]`: `sum_{active v, deg k} s(v)^2 / n_k`.
//!
//! `s(v)` is the size of `v`'s component. The conditional mean component
//! size of a degree-`k` vertex, `S(k)` of Cohen et al., PRL 85, 4626
//! (2000), is estimated by the ratio of ensemble means
//! `DegreeFirstMoment / ActiveFractionByDegree` (not by a per-sample ratio,
//! which would be undefined when no degree-`k` vertex is active).
//!
//! Empty-bin semantics: degrees absent from the graph (`n_k = 0`) keep the
//! array shape stable and record exactly `0.0` in all four ratio arrays;
//! `DegreeOccupied[k] == 0` and `DegreeVertexCount[k] == 0` identify them.
//!
//! Physics anchors: giant fraction and threshold under uniform failure are
//! validated against Callaway, Newman, Strogatz, Watts, PRL 85, 5468
//! (2000) generating functions and the Cohen et al. `f_c = 1 - 1/(kappa-1)`
//! critical fraction; see `CMC.rs/PERCOLATION.md` (N2).

use carlo_rs::{Context, MonteCarlo};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::super::{
    analyze_with_labels, ComponentWorkspace, ObservablePlan, Probability, SiteBernoulli,
    StaticConfiguration,
};
use super::{NetworkPercolationError, NetworkProfile};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Facade over uniform site percolation adding degree-conditioned
/// observables for random-failure studies.
///
/// One `sweep()` draws one failure realization through the shared site law
/// (bitwise-identical activity to `StaticLaw::site(keep)` under the same
/// RNG stream); the following `measure()` analyzes it with component labels
/// and records the F4 scalar plan plus the degree arrays above. Configure
/// `thermalization_sweeps = 0`; samples are i.i.d.
pub struct RandomFailureMC {
    lattice: CsrLattice,
    profile: NetworkProfile,
    keep_probability: Probability,
    configuration: StaticConfiguration,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
    /// Reusable per-sample scratch; all lengths fixed at construction so
    /// the steady-state sample path allocates nothing.
    scratch: DegreeScratch,
}

/// Reusable per-sample degree counters and the shared value buffer.
///
/// One owned struct (rather than loose fields) lets [`RandomFailureMC::measure`]
/// borrow the scratch independently of the workspace that [`analyze_with_labels`]
/// borrows through its result. Public so callers holding their own component
/// labels (for example from a borrowed-graph analysis) can record the same
/// degree arrays through
/// [`RandomFailureMC::record_degree_observables`].
pub struct DegreeScratch {
    /// Component size per canonical vertex index (one pass per sample).
    root_size: Vec<u32>,
    active_count: Vec<u64>,
    giant_count: Vec<u64>,
    first_moment: Vec<u64>,
    second_moment: Vec<u128>,
    /// Refilled in place for each recorded array; capacity is retained.
    values: Vec<f64>,
}

impl DegreeScratch {
    fn new(vertices: usize, bins: usize) -> Self {
        Self {
            root_size: vec![0; vertices],
            active_count: vec![0; bins],
            giant_count: vec![0; bins],
            first_moment: vec![0; bins],
            second_moment: vec![0; bins],
            values: Vec::with_capacity(bins),
        }
    }
}

impl RandomFailureMC {
    /// Build the facade over an owned lattice with keep probability `keep`
    /// and an F4 scalar observable plan.
    ///
    /// Winding observables are rejected (no embedding exists for this
    /// facade), and `V^3` must fit the exact `f64` integer range because
    /// degree second moments sum `s(v)^2` over vertices
    /// (`sum_c s_c^3 <= V^3`).
    pub fn new(
        lattice: CsrLattice,
        keep: f64,
        plan: ObservablePlan,
    ) -> Result<Self, NetworkPercolationError> {
        let vertices = lattice.vertex_count();
        lattice
            .validate()
            .map_err(|error| NetworkPercolationError::Topology(error.to_string()))?;
        if plan.requests_winding() {
            return Err(NetworkPercolationError::WindingUnsupported);
        }
        plan.validate(vertices)
            .map_err(NetworkPercolationError::Observable)?;
        if (vertices as u128) * (vertices as u128) * (vertices as u128)
            > (1_u128 << f64::MANTISSA_DIGITS)
        {
            return Err(NetworkPercolationError::InexactDegreeMoment {
                total_vertices: vertices,
            });
        }
        for (index, query) in plan.boundary_queries().iter().enumerate() {
            if query.vertex_count() != vertices {
                return Err(NetworkPercolationError::QueryDomain {
                    query: index,
                    expected: vertices,
                    actual: query.vertex_count(),
                });
            }
        }
        let keep_probability = Probability::new(keep)
            .map_err(|source| NetworkPercolationError::Probability("keep", source))?;
        let profile =
            NetworkProfile::from_graph(&lattice).map_err(|_| NetworkPercolationError::Capacity)?;
        let configuration = StaticConfiguration::new(vertices, lattice.edge_count());
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(
                vertices,
                plan.boundary_queries().len(),
                /* with_labels = */ true,
            )
            .map_err(|_| NetworkPercolationError::Capacity)?;
        let scratch = DegreeScratch::new(vertices, profile.bins().len());
        Ok(Self {
            lattice,
            profile,
            keep_probability,
            configuration,
            workspace,
            plan,
            scratch,
        })
    }

    /// The owned lattice under analysis.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The cached original-degree profile.
    pub const fn profile(&self) -> &NetworkProfile {
        &self.profile
    }

    /// The site-keep probability of the failure process.
    pub const fn keep_probability(&self) -> f64 {
        self.keep_probability.get()
    }

    /// The scalar F4 observable plan recorded alongside the degree arrays.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// The current activity configuration (for reduction tests).
    pub const fn configuration(&self) -> &StaticConfiguration {
        &self.configuration
    }

    /// Record one sample's degree arrays into `ctx`.
    ///
    /// Exposed separately from `measure` so tests and external analyses can
    /// record degree observables for a configuration they built themselves.
    ///
    /// The value buffer is cleared and refilled in place: after
    /// construction this path performs no allocations.
    pub fn record_degree_observables(
        scratch: &mut DegreeScratch,
        profile: &NetworkProfile,
        labels: &[Option<crate::VertexId>],
        largest: Option<crate::VertexId>,
        ctx: &mut Context<Xoshiro256PlusPlus>,
    ) {
        let DegreeScratch {
            root_size,
            active_count,
            giant_count,
            first_moment,
            second_moment,
            values,
        } = scratch;
        root_size.fill(0);
        for label in labels.iter().flatten() {
            root_size[label.index()] += 1;
        }
        active_count.fill(0);
        giant_count.fill(0);
        first_moment.fill(0);
        second_moment.fill(0);
        for (index, label) in labels.iter().enumerate() {
            let Some(component) = label else {
                continue;
            };
            let bin = profile.degree(index) as usize;
            let size = u64::from(root_size[component.index()]);
            active_count[bin] += 1;
            first_moment[bin] += size;
            second_moment[bin] += u128::from(size) * u128::from(size);
            if largest == Some(*component) {
                giant_count[bin] += 1;
            }
        }

        let bins = profile.bins();
        let mut record = |name: &str, values: &[f64]| {
            ctx.measure_array(name, values)
                .expect("degree-bin shape is fixed by the profile");
        };

        values.clear();
        values.extend(
            bins.iter()
                .map(|&count| if count > 0 { count as f64 } else { 0.0 }),
        );
        record("DegreeVertexCount", values);
        values.clear();
        values.extend(bins.iter().map(|&count| f64::from(u8::from(count > 0))));
        record("DegreeOccupied", values);

        // Empty degree bins record exactly 0.0 — never NaN, never skipped —
        // so every sample keeps the full [0 ..= max_degree] shape.
        values.clear();
        values.extend(bins.iter().enumerate().map(|(bin, &count)| {
            if count > 0 {
                active_count[bin] as f64 / count as f64
            } else {
                0.0
            }
        }));
        record("ActiveFractionByDegree", values);
        values.clear();
        values.extend(bins.iter().enumerate().map(|(bin, &count)| {
            if count > 0 {
                giant_count[bin] as f64 / count as f64
            } else {
                0.0
            }
        }));
        record("GiantMembershipByDegree", values);
        values.clear();
        values.extend(bins.iter().enumerate().map(|(bin, &count)| {
            if count > 0 {
                first_moment[bin] as f64 / count as f64
            } else {
                0.0
            }
        }));
        record("DegreeFirstMoment", values);
        values.clear();
        values.extend(bins.iter().enumerate().map(|(bin, &count)| {
            if count > 0 {
                second_moment[bin] as f64 / count as f64
            } else {
                0.0
            }
        }));
        record("DegreeSecondMoment", values);
    }
}

impl MonteCarlo for RandomFailureMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Random failure is uniform site occupation: delegate to the F2 law
        // so the reduction to `StaticLaw::site(keep)` is by construction.
        SiteBernoulli::new(self.keep_probability.into())
            .sample(&self.lattice, &mut self.configuration, &mut ctx.rng)
            .expect("constructor validated domains");
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Disjoint field borrows: the analysis result borrows `workspace`
        // while the degree record path borrows `scratch` and `profile`.
        let Self {
            lattice,
            profile,
            configuration,
            workspace,
            plan,
            scratch,
            ..
        } = self;
        let result = analyze_with_labels(
            lattice,
            configuration,
            configuration,
            plan.boundary_queries(),
            workspace,
        )
        .expect("constructor prepared matching analysis domains");
        let vertices = lattice.vertex_count();
        let largest = result.summary().largest_component;
        let labels = result.labels().expect("labels were requested");
        plan.record(
            vertices,
            result.summary(),
            result.query_outcomes(),
            None,
            ctx,
        );
        Self::record_degree_observables(scratch, profile, labels, largest, ctx);
    }

    fn name(&self) -> &'static str {
        "RandomFailure"
    }
}
