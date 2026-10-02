//! Owned Carlo.rs adapter for independent static Bernoulli percolation.

use core::fmt;

use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{
    analyze, BondBernoulli, BoundaryQuery, ComponentWorkspace, MixedBernoulli, ObservablePlan,
    ObservablePlanError, Probability, SamplingError, SiteBernoulli, StaticConfiguration,
    StaticObservable,
};
use crate::classical_mc::{build_lattice_from_params, parse_bool, parse_param};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Uniform static law used only at the owned runtime-adapter boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StaticLaw {
    Site {
        probability: Probability,
    },
    Bond {
        probability: Probability,
    },
    Mixed {
        site_probability: Probability,
        bond_probability: Probability,
    },
}

impl StaticLaw {
    pub fn site(probability: f64) -> Result<Self, StaticPercolationError> {
        Ok(Self::Site {
            probability: Probability::new(probability)
                .map_err(|source| StaticPercolationError::Probability("p", source))?,
        })
    }

    pub fn bond(probability: f64) -> Result<Self, StaticPercolationError> {
        Ok(Self::Bond {
            probability: Probability::new(probability)
                .map_err(|source| StaticPercolationError::Probability("p", source))?,
        })
    }

    pub fn mixed(
        site_probability: f64,
        bond_probability: f64,
    ) -> Result<Self, StaticPercolationError> {
        Ok(Self::Mixed {
            site_probability: Probability::new(site_probability)
                .map_err(|source| StaticPercolationError::Probability("p_site", source))?,
            bond_probability: Probability::new(bond_probability)
                .map_err(|source| StaticPercolationError::Probability("p_bond", source))?,
        })
    }

    fn sample<R: Rng + ?Sized>(
        self,
        graph: &impl UndirectedGraphView,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        match self {
            Self::Site { probability } => {
                SiteBernoulli::new(probability.into()).sample(graph, configuration, rng)
            }
            Self::Bond { probability } => {
                BondBernoulli::new(probability.into()).sample(graph, configuration, rng)
            }
            Self::Mixed {
                site_probability,
                bond_probability,
            } => MixedBernoulli::new(site_probability.into(), bond_probability.into()).sample(
                graph,
                configuration,
                rng,
            ),
        }
    }
}

/// Owned convenience adapter for i.i.d. site, bond, or mixed percolation.
///
/// One `sweep()` draws one complete realization and the following `measure()`
/// analyzes it. Configure `thermalization_sweeps = 0`; the adapter cannot see
/// scheduler configuration and therefore cannot enforce this runtime setting.
/// Borrowed graphs and custom laws use the core APIs directly, optionally with
/// `Run::from_parts()`.
pub struct StaticPercolationMC {
    lattice: CsrLattice,
    law: StaticLaw,
    configuration: StaticConfiguration,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
}

impl StaticPercolationMC {
    pub fn new(
        lattice: CsrLattice,
        law: StaticLaw,
        plan: ObservablePlan,
    ) -> Result<Self, StaticPercolationError> {
        lattice
            .validate()
            .map_err(|error| StaticPercolationError::Topology(error.to_string()))?;
        plan.validate(lattice.vertex_count())
            .map_err(StaticPercolationError::Observable)?;
        for (index, query) in plan.boundary_queries().iter().enumerate() {
            if query.vertex_count() != lattice.vertex_count() {
                return Err(StaticPercolationError::QueryDomain {
                    query: index,
                    expected: lattice.vertex_count(),
                    actual: query.vertex_count(),
                });
            }
        }
        let configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| StaticPercolationError::Capacity)?;
        Ok(Self {
            lattice,
            law,
            configuration,
            workspace,
            plan,
        })
    }

    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    pub const fn law(&self) -> StaticLaw {
        self.law
    }

    pub const fn configuration(&self) -> &StaticConfiguration {
        &self.configuration
    }

    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }
}

impl MonteCarlo for StaticPercolationMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.law
            .sample(&self.lattice, &mut self.configuration, &mut ctx.rng)
            .expect("constructor validated static-law domains");
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let result = analyze(
            &self.lattice,
            &self.configuration,
            &self.configuration,
            self.plan.boundary_queries(),
            &mut self.workspace,
        )
        .expect("constructor prepared matching analysis domains");
        self.plan.record(
            self.lattice.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            ctx,
        );
    }

    fn name(&self) -> &'static str {
        "StaticPercolation"
    }
}

impl FromParams for StaticPercolationMC {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        parse_law(params)?;
        parse_observables(params)?;
        parse_index_list(params, "spanning_from")?;
        parse_index_list(params, "spanning_to")?;
        // Validate the full parameter-built topology and boundary sets here
        // as well, so `validate_params` accepts exactly what `from_params`
        // accepts (mirroring the `ClassicalMC` adapters).
        let pbc = parse_bool(params, "pbc", false)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        resolve_spanning_sets(params, &lattice)?;
        Ok(())
    }

    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let law = parse_law(params)?;
        let pbc = parse_bool(params, "pbc", false)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        let observables = parse_observables(params)?;
        let (from, to) = resolve_spanning_sets(params, &lattice)?;
        let query = BoundaryQuery::new(&lattice, &from, &to)
            .map_err(|error| invalid("spanning_from", error))?;
        let plan = ObservablePlan::new(observables, vec![query]);
        Self::new(lattice, law, plan).map_err(|error| invalid("percolation", error))
    }
}

fn parse_law(params: &Params) -> Result<StaticLaw, CarloError> {
    let mode = parse_param::<String>(params, "mode")?.unwrap_or_else(|| "site".to_string());
    let probability = |name: &str| -> Result<f64, CarloError> {
        parse_param::<f64>(params, name)?
            .ok_or_else(|| invalid(name, "missing required occupation probability"))
    };
    match mode.as_str() {
        "site" | "bond" => {
            if params.contains("p_site") || params.contains("p_bond") {
                return Err(invalid(
                    "p_site",
                    "pure modes accept only `p`, not `p_site` or `p_bond`",
                ));
            }
            let p = probability("p")?;
            if mode == "site" {
                StaticLaw::site(p)
            } else {
                StaticLaw::bond(p)
            }
            .map_err(|error| invalid("p", error))
        }
        "site-bond" => {
            if params.contains("p") {
                return Err(invalid(
                    "p",
                    "mixed mode accepts only `p_site` and `p_bond`",
                ));
            }
            StaticLaw::mixed(probability("p_site")?, probability("p_bond")?)
                .map_err(|error| invalid("p_site", error))
        }
        _ => Err(invalid("mode", "expected `site`, `bond`, or `site-bond`")),
    }
}

fn parse_observables(params: &Params) -> Result<Vec<StaticObservable>, CarloError> {
    let Some(raw) = parse_param::<String>(params, "observables")? else {
        return Ok(vec![
            StaticObservable::ActiveVertexCount,
            StaticObservable::ActiveEdgeCount,
            StaticObservable::ComponentCount,
            StaticObservable::LargestSize,
            StaticObservable::GiantFraction,
            StaticObservable::RawSecondMoment,
            StaticObservable::FiniteClusterSusceptibility,
        ]);
    };
    if raw.trim().is_empty() {
        return Err(invalid("observables", "selection must not be empty"));
    }
    raw.split(',')
        .map(|token| match token.trim() {
            "active-vertices" => Ok(StaticObservable::ActiveVertexCount),
            "active-edges" => Ok(StaticObservable::ActiveEdgeCount),
            "components" => Ok(StaticObservable::ComponentCount),
            "largest-size" => Ok(StaticObservable::LargestSize),
            "giant-fraction" => Ok(StaticObservable::GiantFraction),
            "raw-second-moment" => Ok(StaticObservable::RawSecondMoment),
            "finite-cluster-susceptibility" => Ok(StaticObservable::FiniteClusterSusceptibility),
            token => Err(invalid(
                "observables",
                format!("unknown observable `{token}`"),
            )),
        })
        .collect()
}

fn resolve_spanning_sets(
    params: &Params,
    lattice: &CsrLattice,
) -> Result<(Vec<usize>, Vec<usize>), CarloError> {
    match (
        parse_index_list(params, "spanning_from")?,
        parse_index_list(params, "spanning_to")?,
    ) {
        (Some(from), Some(to)) => return Ok((from, to)),
        (Some(_), None) | (None, Some(_)) => {
            return Err(invalid(
                "spanning_to",
                "`spanning_from` and `spanning_to` must be given together",
            ));
        }
        (None, None) => {}
    }
    let inferred = if params.contains("Lx") {
        "hypercubic"
    } else {
        "chain"
    };
    let lattice_type = params
        .get::<String>("lattice_type")
        .unwrap_or_else(|| inferred.to_string())
        .to_ascii_lowercase();
    match lattice_type.as_str() {
        "chain" if lattice.vertex_count() > 0 => Ok((vec![0], vec![lattice.vertex_count() - 1])),
        "square" => square_columns(params, lattice.vertex_count()),
        "hypercubic" if params.contains("Lx") && !params.contains("Lz") => {
            square_columns(params, lattice.vertex_count())
        }
        _ => Err(invalid(
            "spanning_from",
            format!("lattice type `{lattice_type}` requires explicit boundary sets"),
        )),
    }
}

fn square_columns(
    params: &Params,
    vertices: usize,
) -> Result<(Vec<usize>, Vec<usize>), CarloError> {
    let lx = parse_param::<usize>(params, "Lx")?.unwrap_or(4);
    let ly = parse_param::<usize>(params, "Ly")?.unwrap_or(lx);
    if lx == 0 || ly == 0 || lx.checked_mul(ly) != Some(vertices) {
        return Err(invalid("Lx", "square dimensions do not match topology"));
    }
    Ok((
        (0..ly).map(|row| row * lx).collect(),
        (0..ly).map(|row| row * lx + lx - 1).collect(),
    ))
}

fn parse_index_list(params: &Params, name: &str) -> Result<Option<Vec<usize>>, CarloError> {
    let Some(raw) = parse_param::<String>(params, name)? else {
        return Ok(None);
    };
    if raw.trim().is_empty() {
        return Err(invalid(name, "site list must not be empty"));
    }
    raw.split(',')
        .map(|token| {
            let token = token.trim();
            if token.is_empty() {
                return Err(invalid(name, "site list contains an empty entry"));
            }
            token
                .parse::<usize>()
                .map_err(|_| invalid(name, format!("cannot parse site index `{token}`")))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

fn invalid(field: impl Into<String>, reason: impl fmt::Display) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.to_string(),
    }
}

#[derive(Debug)]
pub enum StaticPercolationError {
    Probability(&'static str, super::ProbabilityError),
    Topology(String),
    Observable(ObservablePlanError),
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    Capacity,
}

impl fmt::Display for StaticPercolationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Probability(name, source) => write!(formatter, "invalid `{name}`: {source}"),
            Self::Topology(reason) => write!(formatter, "invalid topology: {reason}"),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} has V={actual}; adapter topology has V={expected}"
            ),
            Self::Capacity => formatter.write_str("could not reserve component workspace"),
        }
    }
}

impl std::error::Error for StaticPercolationError {}

#[cfg(test)]
mod tests {
    use super::*;
    use carlo_rs::{FromParams, RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

    fn rng() -> Xoshiro256PlusPlus {
        Xoshiro256PlusPlus::seed_from_u64(7)
    }

    fn square_params() -> Params {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 2);
        params.set("Ly", 2);
        params.set("p", 0.5);
        params
    }

    #[test]
    fn from_params_schema_and_structural_defaults() {
        let site = StaticPercolationMC::from_params(&square_params(), &mut rng()).unwrap();
        assert!(matches!(site.law(), StaticLaw::Site { .. }));
        assert_eq!(site.observable_plan().boundary_queries()[0].from().len(), 2);

        let mut bond = square_params();
        bond.set("mode", "bond");
        bond.set("spanning_from", "0,1");
        bond.set("spanning_to", "2,3");
        assert!(matches!(
            StaticPercolationMC::from_params(&bond, &mut rng())
                .unwrap()
                .law(),
            StaticLaw::Bond { .. }
        ));

        let mut mixed = square_params();
        mixed.set("mode", "site-bond");
        mixed.set("p_site", 0.7);
        mixed.set("p_bond", 0.4);
        assert!(mixed.get::<f64>("p").is_some());
        assert!(StaticPercolationMC::validate_params(&mixed).is_err());
        let mut mixed = Params::new();
        mixed.set("lattice_type", "square");
        mixed.set("Lx", 2);
        mixed.set("Ly", 2);
        mixed.set("mode", "site-bond");
        mixed.set("p_site", 0.7);
        mixed.set("p_bond", 0.4);
        assert!(matches!(
            StaticPercolationMC::from_params(&mixed, &mut rng())
                .unwrap()
                .law(),
            StaticLaw::Mixed { .. }
        ));
    }

    #[test]
    fn invalid_params_are_rejected_comprehensively() {
        let invalid = |mut params: Params, key: &str, value: &str| {
            params.set(key, value);
            assert!(StaticPercolationMC::from_params(&params, &mut rng()).is_err());
        };
        invalid(square_params(), "mode", "edge");
        invalid(square_params(), "p", "NaN");
        invalid(square_params(), "p", "1.1");
        invalid(square_params(), "p_site", "0.5");
        invalid(square_params(), "spanning_from", "0,");
        invalid(square_params(), "spanning_from", "0");
        let mut out_of_range = square_params();
        out_of_range.set("spanning_from", "0");
        out_of_range.set("spanning_to", "4");
        assert!(StaticPercolationMC::from_params(&out_of_range, &mut rng()).is_err());
        invalid(square_params(), "observables", "largest");

        let mut mixed_missing = Params::new();
        mixed_missing.set("mode", "site-bond");
        mixed_missing.set("p_site", 0.5);
        assert!(StaticPercolationMC::validate_params(&mixed_missing).is_err());

        let mut triangular = Params::new();
        triangular.set("lattice_type", "triangular");
        triangular.set("Lx", 2);
        triangular.set("Ly", 2);
        triangular.set("p", 0.5);
        assert!(StaticPercolationMC::from_params(&triangular, &mut rng()).is_err());
        // `validate_params` accepts exactly what `from_params` accepts, so the
        // topology and boundary defaults are checked there too (triangular has
        // no structural spanning default and no explicit sets were given).
        assert!(StaticPercolationMC::validate_params(&triangular).is_err());
        assert!(StaticPercolationMC::validate_params(&square_params()).is_ok());
    }

    #[test]
    fn scheduler_results_and_json_follow_final_schema() {
        let mut params = square_params();
        params.set(
            "observables",
            "active-vertices,largest-size,giant-fraction,finite-cluster-susceptibility",
        );
        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 0,
                measurement_sweeps: 64,
                binsize: 1,
                base_seed: 9,
                ..Default::default()
            },
        )
        .run_one::<StaticPercolationMC>(&params);
        for name in [
            "ActiveVertexCount",
            "LargestSize",
            "GiantFraction",
            "FiniteClusterSusceptibility",
            "FiniteClusterSusceptibilityDefined",
            "BoundaryCrossing0",
        ] {
            assert!(results.get(name).is_some(), "missing {name}");
        }
        assert!(results.get("RawSecondMoment").is_none());
        assert_eq!(results.metadata().thermalization_sweeps, 0);
        assert_eq!(results.metadata().measurement_sweeps, 64);
        let json: serde_json::Value = serde_json::from_str(&results.to_json().unwrap()).unwrap();
        assert!(json["observables"]["LargestSize"].is_object());
        assert_eq!(json["metadata"]["measurement_sweeps"], 64);
    }

    #[test]
    fn undefined_susceptibility_records_indicator_for_every_sample_only() {
        let mut params = Params::new();
        params.set("lattice_type", "chain");
        params.set("L", 1);
        params.set("p", 1.0);
        params.set("observables", "finite-cluster-susceptibility");
        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 0,
                measurement_sweeps: 7,
                binsize: 1,
                ..Default::default()
            },
        )
        .run_one::<StaticPercolationMC>(&params);
        assert!(results.get("FiniteClusterSusceptibility").is_none());
        let defined = results
            .get("FiniteClusterSusceptibilityDefined")
            .expect("defined indicator");
        assert_eq!(defined.mean, 0.0);
        assert_eq!(defined.n_bins, 7);
    }
}
