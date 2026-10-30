//! Explicit run specification, parameter round-trip, and provenance.
//!
//! The owned adapter has three ways in: a bespoke constructor, `FromParams`
//! (only the static preset), and [`Run::from_parts`](carlo_rs::Run::from_parts)
//! for laws the schema cannot express. [`PercolationSpec`] is the explicit
//! middle ground for the independent-occupation family: every field that
//! affects the run is present in the value (there are no hidden defaults once
//! a spec exists), the Params preset converts in both directions, and the
//! spec emits a versioned provenance object that can be stamped onto
//! [`Results`].
//!
//! Scope: site/bond/mixed [`StaticLaw`] draws composed with the
//! independent-occupation process, built-in [`StaticObservable`]s, and the
//! report-time sample-moment estimator. Custom laws, custom observables, and
//! the other process families stay on their dedicated constructors (or
//! `Run::from_parts`) until a later stage widens the schema.

use core::fmt;

use carlo_rs::{Params, Results};
use serde_json::{json, Value};

use super::carlo::{parse_law, parse_observables, resolve_spanning_sets};
use super::{
    BoundaryQuery, BoundaryQueryError, ObservablePlan, PercolationAdapter, StaticLaw,
    StaticObservable, StaticPercolationError,
};
use crate::classical_mc::{parse_bool, parse_param};
use crate::{
    build_chain, build_chain_with_embedding, build_honeycomb, build_hypercubic,
    build_hypercubic_with_embedding, build_kagome, build_square, build_square_with_embedding,
    build_triangular_with_embedding, BondType, CsrLattice, EmbeddingError, LatticeEmbedding,
};

/// Lattice shape with every dimension resolved (no parameter defaults).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LatticeKind {
    /// One-dimensional chain of `length` vertices.
    Chain { length: usize },
    /// Two-dimensional square lattice `lx * ly`.
    Square { lx: usize, ly: usize },
    /// Three-dimensional cubic lattice `lx * ly * lz`.
    Cubic { lx: usize, ly: usize, lz: usize },
    /// Two-dimensional triangular torus; periodic only.
    Triangular { lx: usize, ly: usize },
    /// Two-dimensional honeycomb torus; periodic only, even `lx`.
    Honeycomb { lx: usize, ly: usize },
    /// Two-dimensional kagome torus; periodic only.
    Kagome { lx: usize, ly: usize },
}

impl LatticeKind {
    /// Canonical `lattice_type` parameter name.
    pub const fn param_name(self) -> &'static str {
        match self {
            Self::Chain { .. } => "chain",
            Self::Square { .. } => "square",
            Self::Cubic { .. } => "cubic",
            Self::Triangular { .. } => "triangular",
            Self::Honeycomb { .. } => "honeycomb",
            Self::Kagome { .. } => "kagome",
        }
    }

    /// Whether this shape requires periodic boundaries.
    pub const fn periodic_only(self) -> bool {
        matches!(
            self,
            Self::Triangular { .. } | Self::Honeycomb { .. } | Self::Kagome { .. }
        )
    }
}

/// Resolved topology: lattice shape plus an explicit boundary choice.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TopologySpec {
    lattice: LatticeKind,
    periodic: bool,
}

impl TopologySpec {
    pub const fn new(lattice: LatticeKind, periodic: bool) -> Self {
        Self { lattice, periodic }
    }

    /// Resolve the standard lattice parameters into an explicit topology.
    ///
    /// Mirrors `build_lattice_from_params` (including its documented defaults,
    /// which become explicit here) without building the graph.
    pub fn from_params(params: &Params) -> Result<Self, SpecError> {
        let periodic = parse_bool(params, "pbc", false).map_err(SpecError::Params)?;
        let inferred = if params.contains("Lx") {
            "hypercubic"
        } else {
            "chain"
        };
        let lattice_type = params
            .get::<String>("lattice_type")
            .unwrap_or_else(|| inferred.to_string())
            .to_ascii_lowercase();
        let lattice = match lattice_type.as_str() {
            "chain" => {
                let length = if params.contains("L") {
                    dimension(params, "L", 10)?
                } else {
                    dimension(params, "Lx", 10)?
                };
                LatticeKind::Chain { length }
            }
            "square" => {
                let lx = dimension(params, "Lx", 4)?;
                let ly = dimension(params, "Ly", lx)?;
                LatticeKind::Square { lx, ly }
            }
            "cubic" => {
                let lx = dimension(params, "Lx", 4)?;
                let ly = dimension(params, "Ly", lx)?;
                let lz = dimension(params, "Lz", lx)?;
                LatticeKind::Cubic { lx, ly, lz }
            }
            "hypercubic" => {
                if params.contains("Lx") {
                    let lx = dimension(params, "Lx", 4)?;
                    let ly = dimension(params, "Ly", lx)?;
                    if params.contains("Lz") {
                        let lz = dimension(params, "Lz", lx)?;
                        LatticeKind::Cubic { lx, ly, lz }
                    } else {
                        LatticeKind::Square { lx, ly }
                    }
                } else {
                    LatticeKind::Chain {
                        length: dimension(params, "L", 10)?,
                    }
                }
            }
            "triangular" => {
                let lx = dimension(params, "Lx", 4)?;
                let ly = dimension(params, "Ly", lx)?;
                LatticeKind::Triangular { lx, ly }
            }
            "honeycomb" => {
                let lx = dimension(params, "Lx", 4)?;
                let ly = dimension(params, "Ly", lx)?;
                LatticeKind::Honeycomb { lx, ly }
            }
            "kagome" => {
                let lx = dimension(params, "Lx", 2)?;
                let ly = dimension(params, "Ly", lx)?;
                LatticeKind::Kagome { lx, ly }
            }
            _ => {
                return Err(SpecError::UnknownLattice {
                    lattice: lattice_type,
                })
            }
        };
        Ok(Self { lattice, periodic })
    }

    /// Explicit parameters that reproduce this topology.
    pub fn to_params(&self) -> Params {
        let mut params = Params::new();
        params.set("lattice_type", self.lattice.param_name());
        params.set("pbc", self.periodic);
        match self.lattice {
            LatticeKind::Chain { length } => params.set("L", length),
            LatticeKind::Square { lx, ly } => {
                params.set("Lx", lx);
                params.set("Ly", ly);
            }
            LatticeKind::Cubic { lx, ly, lz } => {
                params.set("Lx", lx);
                params.set("Ly", ly);
                params.set("Lz", lz);
            }
            LatticeKind::Triangular { lx, ly }
            | LatticeKind::Honeycomb { lx, ly }
            | LatticeKind::Kagome { lx, ly } => {
                params.set("Lx", lx);
                params.set("Ly", ly);
            }
        }
        params
    }

    pub const fn lattice(&self) -> LatticeKind {
        self.lattice
    }

    pub const fn periodic(&self) -> bool {
        self.periodic
    }

    /// Build the graph (and its periodic embedding where the family has one).
    ///
    /// Embedding availability exactly mirrors the parameter preset: chain,
    /// square, cubic, and triangular tori carry embeddings; honeycomb and
    /// kagome tori do not.
    pub fn build(&self) -> Result<(CsrLattice, Option<LatticeEmbedding>), SpecError> {
        if self.lattice.periodic_only() && !self.periodic {
            return Err(SpecError::PeriodicRequired {
                lattice: self.lattice.param_name(),
            });
        }
        match (self.lattice, self.periodic) {
            (LatticeKind::Chain { length }, true) => build_chain_with_embedding(length, true)
                .map(|(lattice, embedding)| (lattice, Some(embedding)))
                .map_err(SpecError::Embedding),
            (LatticeKind::Square { lx, ly }, true) => build_square_with_embedding(lx, ly, true)
                .map(|(lattice, embedding)| (lattice, Some(embedding)))
                .map_err(SpecError::Embedding),
            (LatticeKind::Cubic { lx, ly, lz }, true) => build_hypercubic_with_embedding(
                &[lx, ly, lz],
                &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                true,
            )
            .map(|(lattice, embedding)| (lattice, Some(embedding)))
            .map_err(SpecError::Embedding),
            (LatticeKind::Triangular { lx, ly }, true) => {
                if lx < 2 || ly < 2 {
                    return Err(SpecError::InvalidTopology {
                        reason: "triangular lattice requires both dimensions >= 2".to_string(),
                    });
                }
                build_triangular_with_embedding(lx, ly)
                    .map(|(lattice, embedding)| (lattice, Some(embedding)))
                    .map_err(SpecError::Embedding)
            }
            (LatticeKind::Honeycomb { lx, ly }, true) => {
                if lx < 2 || ly < 2 || !lx.is_multiple_of(2) {
                    return Err(SpecError::InvalidTopology {
                        reason: "honeycomb lattice requires even Lx >= 2 and Ly >= 2".to_string(),
                    });
                }
                Ok((build_honeycomb(lx, ly), None))
            }
            (LatticeKind::Kagome { lx, ly }, true) => {
                if lx < 2 || ly < 2 {
                    return Err(SpecError::InvalidTopology {
                        reason: "kagome lattice requires both dimensions >= 2".to_string(),
                    });
                }
                Ok((build_kagome(lx, ly), None))
            }
            (LatticeKind::Chain { length }, false) => Ok((build_chain(length, false), None)),
            (LatticeKind::Square { lx, ly }, false) => Ok((build_square(lx, ly, false), None)),
            (LatticeKind::Cubic { lx, ly, lz }, false) => Ok((
                build_hypercubic(
                    &[lx, ly, lz],
                    &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                    false,
                ),
                None,
            )),
            _ => Err(SpecError::PeriodicRequired {
                lattice: self.lattice.param_name(),
            }),
        }
    }

    /// Canonical provenance label, for example `square/16x16/pbc`.
    pub fn label(&self) -> String {
        let shape = match self.lattice {
            LatticeKind::Chain { length } => format!("{}", length),
            LatticeKind::Square { lx, ly }
            | LatticeKind::Triangular { lx, ly }
            | LatticeKind::Honeycomb { lx, ly }
            | LatticeKind::Kagome { lx, ly } => format!("{lx}x{ly}"),
            LatticeKind::Cubic { lx, ly, lz } => format!("{lx}x{ly}x{lz}"),
        };
        let boundary = if self.periodic { "pbc" } else { "open" };
        format!("{}/{shape}/{boundary}", self.lattice.param_name())
    }
}

fn dimension(params: &Params, name: &'static str, default: usize) -> Result<usize, SpecError> {
    let value = parse_param::<usize>(params, name)
        .map_err(SpecError::Params)?
        .unwrap_or(default);
    if value == 0 {
        return Err(SpecError::InvalidDimension { field: name });
    }
    Ok(value)
}

/// The sampling process the spec composes. One variant today: the preset
/// draws independent configurations, so thermalization is meaningless.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProcessSpec {
    /// One independent realization per sweep ([`StaticLaw`] draw).
    IndependentOccupation,
}

impl ProcessSpec {
    pub const fn label(self) -> &'static str {
        match self {
            Self::IndependentOccupation => "independent-occupation/v1",
        }
    }
}

/// The report-time estimator of recorded sample moments.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EstimatorSpec {
    /// Carlo.rs sample estimator: per-bin mean, then mean and standard
    /// error over bins (the v2 default; binsize stays in `RunConfig`).
    SampleMoments,
}

impl EstimatorSpec {
    pub const fn label(self) -> &'static str {
        match self {
            Self::SampleMoments => "sample-moments/v1",
        }
    }
}

/// An explicit boundary-crossing query between two dense vertex sets.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BoundarySets {
    from: Vec<usize>,
    to: Vec<usize>,
}

impl BoundarySets {
    pub fn new(from: Vec<usize>, to: Vec<usize>) -> Self {
        Self { from, to }
    }

    pub fn from(&self) -> &[usize] {
        &self.from
    }

    pub fn to(&self) -> &[usize] {
        &self.to
    }
}

/// A fully resolved independent-occupation percolation run.
///
/// Every field is explicit; construct with [`Self::from_params`] to resolve
/// the documented parameter defaults, or with a struct literal for a
/// programmatic preset. [`Self::build`] validates the combination and hands
/// back the scheduler-ready adapter.
///
/// ```
/// use carlo_rs::Params;
/// use cmc_rs::{GraphView, PercolationSpec};
///
/// let mut params = Params::new();
/// params.set("lattice_type", "square");
/// params.set("Lx", 4);
/// params.set("Ly", 4);
/// params.set("p", 0.5);
///
/// let spec = PercolationSpec::from_params(&params)?;
/// let adapter = spec.build()?;
/// assert_eq!(adapter.lattice().vertex_count(), 16);
/// assert_eq!(spec.provenance().topology(), "square/4x4/open");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct PercolationSpec {
    pub topology: TopologySpec,
    pub law: StaticLaw,
    pub process: ProcessSpec,
    pub observables: Vec<StaticObservable>,
    /// `None` records no boundary query; `Some` records one query.
    pub boundaries: Option<BoundarySets>,
    pub estimator: EstimatorSpec,
}

impl PercolationSpec {
    /// Resolve the standard static-preset parameters into an explicit spec.
    ///
    /// The two documented defaults (`observables` = the v2 built-in list,
    /// structural `spanning_from`/`spanning_to`) become explicit values, so
    /// the returned spec is self-contained.
    pub fn from_params(params: &Params) -> Result<Self, SpecError> {
        let topology = TopologySpec::from_params(params)?;
        let law = parse_law(params).map_err(SpecError::Params)?;
        let observables = parse_observables(params).map_err(SpecError::Params)?;
        let (lattice, _) = topology.build()?;
        let (from, to) = resolve_spanning_sets(params, &lattice).map_err(SpecError::Params)?;
        Ok(Self {
            topology,
            law,
            process: ProcessSpec::IndependentOccupation,
            observables,
            boundaries: Some(BoundarySets::new(from, to)),
            estimator: EstimatorSpec::SampleMoments,
        })
    }

    /// Explicit parameters for the schema-expressible subset of this spec.
    ///
    /// Heterogeneous laws, an empty observable list, and a missing or empty
    /// boundary query are programmatic-only shapes; they are typed rejections
    /// here rather than silently dropped.
    pub fn to_params(&self) -> Result<Params, SpecError> {
        let mut params = self.topology.to_params();
        match &self.law {
            StaticLaw::Site { probability } => {
                params.set("mode", "site");
                params.set("p", probability.get());
            }
            StaticLaw::Bond { probability } => {
                params.set("mode", "bond");
                params.set("p", probability.get());
            }
            StaticLaw::Mixed {
                site_probability,
                bond_probability,
            } => {
                params.set("mode", "site-bond");
                params.set("p_site", site_probability.get());
                params.set("p_bond", bond_probability.get());
            }
            StaticLaw::SiteHeterogeneous { .. }
            | StaticLaw::BondHeterogeneous { .. }
            | StaticLaw::MixedHeterogeneous { .. } => {
                return Err(SpecError::LawNotParamsExpressible);
            }
        }
        if self.observables.is_empty() {
            return Err(SpecError::EmptyObservables);
        }
        let observables = self
            .observables
            .iter()
            .map(|observable| observable_token(*observable))
            .collect::<Vec<_>>()
            .join(",");
        params.set("observables", observables);
        let Some(boundaries) = &self.boundaries else {
            return Err(SpecError::NoBoundaryQuery);
        };
        if boundaries.from.is_empty() || boundaries.to.is_empty() {
            return Err(SpecError::EmptyBoundarySets);
        }
        params.set(
            "spanning_from",
            boundaries
                .from
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(","),
        );
        params.set(
            "spanning_to",
            boundaries
                .to
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(","),
        );
        Ok(params)
    }

    /// Validate the spec and build its scheduler-ready adapter.
    ///
    /// Winding observables require a topology with a periodic embedding; the
    /// returned adapter carries that embedding only when the plan requests
    /// winding (exactly the parameter-preset behavior).
    pub fn build(&self) -> Result<PercolationAdapter<StaticLaw>, SpecError> {
        let (lattice, embedding) = self.topology.build()?;
        let queries = match &self.boundaries {
            Some(boundaries) => {
                vec![
                    BoundaryQuery::new(&lattice, &boundaries.from, &boundaries.to)
                        .map_err(SpecError::Query)?,
                ]
            }
            None => Vec::new(),
        };
        let plan = ObservablePlan::new(self.observables.iter().copied(), queries);
        let adapter = match embedding.filter(|_| plan.requests_winding()) {
            Some(embedding) => {
                PercolationAdapter::with_embedding(lattice, embedding, self.law.clone(), plan)
            }
            None => PercolationAdapter::new(lattice, self.law.clone(), plan),
        };
        adapter.map_err(SpecError::Adapter)
    }

    /// Versioned provenance for this spec.
    pub fn provenance(&self) -> PercolationProvenance {
        PercolationProvenance {
            version: 1,
            topology: self.topology.label(),
            law: law_label(&self.law),
            process: self.process.label().to_string(),
            observables: self
                .observables
                .iter()
                .map(|observable| observable.name().to_string())
                .collect(),
            boundary_queries: usize::from(self.boundaries.is_some()),
            estimator: self.estimator.label().to_string(),
        }
    }

    /// Stamp the provenance object onto finalized results under `percolation`.
    pub fn stamp(&self, results: &mut Results) {
        results.set_attribute("percolation", self.provenance().to_json());
    }
}

/// Versioned provenance attached to results.
#[derive(Clone, Debug, PartialEq)]
pub struct PercolationProvenance {
    version: u32,
    topology: String,
    law: String,
    process: String,
    observables: Vec<String>,
    boundary_queries: usize,
    estimator: String,
}

impl PercolationProvenance {
    pub const fn version(&self) -> u32 {
        self.version
    }

    pub fn topology(&self) -> &str {
        &self.topology
    }

    pub fn law(&self) -> &str {
        &self.law
    }

    pub fn process(&self) -> &str {
        &self.process
    }

    pub fn observables(&self) -> &[String] {
        &self.observables
    }

    pub const fn boundary_queries(&self) -> usize {
        self.boundary_queries
    }

    pub fn estimator(&self) -> &str {
        &self.estimator
    }

    /// JSON object recorded in `Results` metadata.
    pub fn to_json(&self) -> Value {
        json!({
            "schema": "cmc.percolation.provenance",
            "version": self.version,
            "topology": self.topology,
            "law": self.law,
            "process": self.process,
            "observables": self.observables,
            "boundary_queries": self.boundary_queries,
            "estimator": self.estimator,
        })
    }
}

fn observable_token(observable: StaticObservable) -> &'static str {
    match observable {
        StaticObservable::ActiveVertexCount => "active-vertices",
        StaticObservable::ActiveEdgeCount => "active-edges",
        StaticObservable::ComponentCount => "components",
        StaticObservable::LargestSize => "largest-size",
        StaticObservable::GiantFraction => "giant-fraction",
        StaticObservable::RawSecondMoment => "raw-second-moment",
        StaticObservable::FiniteClusterSusceptibility => "finite-cluster-susceptibility",
        StaticObservable::Wrapping => "wrapping",
        StaticObservable::WindingRank => "winding-rank",
    }
}

fn law_label(law: &StaticLaw) -> String {
    match law {
        StaticLaw::Site { probability } => format!("site/p={}", probability.get()),
        StaticLaw::Bond { probability } => format!("bond/p={}", probability.get()),
        StaticLaw::Mixed {
            site_probability,
            bond_probability,
        } => format!(
            "mixed/p_site={},p_bond={}",
            site_probability.get(),
            bond_probability.get()
        ),
        StaticLaw::SiteHeterogeneous { probabilities } => {
            format!("site-heterogeneous/vertices={}", probabilities.len())
        }
        StaticLaw::BondHeterogeneous { probabilities } => {
            format!("bond-heterogeneous/edges={}", probabilities.len())
        }
        StaticLaw::MixedHeterogeneous {
            site_probabilities,
            bond_probabilities,
        } => format!(
            "mixed-heterogeneous/vertices={},edges={}",
            site_probabilities.len(),
            bond_probabilities.len()
        ),
    }
}

/// Typed rejections for spec resolution, parameter round-trip, and building.
#[derive(Debug)]
pub enum SpecError {
    /// A wrapped parameter error from the preset parser.
    Params(carlo_rs::CarloError),
    /// Unknown `lattice_type` string.
    UnknownLattice { lattice: String },
    /// A zero dimension parameter.
    InvalidDimension { field: &'static str },
    /// A topology constraint that the shape itself cannot satisfy.
    InvalidTopology { reason: String },
    /// Honeycomb/kagome/triangular need periodic boundaries.
    PeriodicRequired { lattice: &'static str },
    /// The embedding builder rejected the topology.
    Embedding(EmbeddingError),
    /// A boundary vertex index is outside the built topology.
    Query(BoundaryQueryError),
    /// The adapter rejected the law/plan/embedding combination.
    Adapter(StaticPercolationError),
    /// Heterogeneous fields have no parameter representation.
    LawNotParamsExpressible,
    /// The Params preset requires at least one observable.
    EmptyObservables,
    /// A spec without a boundary query has no Params representation.
    NoBoundaryQuery,
    /// Empty boundary sets have no parameter representation.
    EmptyBoundarySets,
}

impl fmt::Display for SpecError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Params(source) => write!(formatter, "invalid parameters: {source}"),
            Self::UnknownLattice { lattice } => write!(
                formatter,
                "unknown lattice `{lattice}`; expected chain, square, cubic, hypercubic, \
                 triangular, honeycomb, or kagome"
            ),
            Self::InvalidDimension { field } => {
                write!(formatter, "`{field}` must be positive")
            }
            Self::InvalidTopology { reason } => {
                write!(formatter, "invalid topology: {reason}")
            }
            Self::PeriodicRequired { lattice } => write!(
                formatter,
                "lattice `{lattice}` requires periodic boundaries"
            ),
            Self::Embedding(source) => write!(formatter, "invalid periodic embedding: {source}"),
            Self::Query(source) => write!(formatter, "invalid boundary sets: {source}"),
            Self::Adapter(source) => write!(formatter, "invalid percolation run: {source}"),
            Self::LawNotParamsExpressible => formatter
                .write_str("heterogeneous law fields cannot be expressed as run parameters"),
            Self::EmptyObservables => {
                formatter.write_str("a parameter-expressible run requires at least one observable")
            }
            Self::NoBoundaryQuery => {
                formatter.write_str("a parameter-expressible run requires one boundary query")
            }
            Self::EmptyBoundarySets => {
                formatter.write_str("empty boundary sets cannot be expressed as run parameters")
            }
        }
    }
}

impl std::error::Error for SpecError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Params(source) => Some(source),
            Self::Embedding(source) => Some(source),
            Self::Query(source) => Some(source),
            Self::Adapter(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GraphView, Probability, StaticPercolationMC, UndirectedGraphView};
    use carlo_rs::FromParams;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    fn params(entries: &[(&str, &str)]) -> Params {
        let mut params = Params::new();
        for (key, value) in entries {
            params.set(key, *value);
        }
        params
    }

    fn preset_cases() -> Vec<Params> {
        vec![
            params(&[("p", "0.5")]),
            params(&[
                ("lattice_type", "chain"),
                ("L", "5"),
                ("pbc", "true"),
                ("p", "0.4"),
            ]),
            params(&[("lattice_type", "square"), ("p", "0.5")]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "3"),
                ("Ly", "5"),
                ("p", "0.5"),
            ]),
            params(&[
                ("lattice_type", "cubic"),
                ("Lx", "2"),
                ("Ly", "3"),
                ("Lz", "4"),
                ("p", "0.5"),
                ("spanning_from", "0,1"),
                ("spanning_to", "2,3"),
            ]),
            params(&[
                ("lattice_type", "hypercubic"),
                ("Lx", "3"),
                ("Ly", "4"),
                ("p", "0.5"),
            ]),
            params(&[
                ("lattice_type", "hypercubic"),
                ("L", "6"),
                ("p", "0.5"),
                ("spanning_from", "0"),
                ("spanning_to", "5"),
            ]),
            params(&[
                ("lattice_type", "triangular"),
                ("Lx", "3"),
                ("Ly", "3"),
                ("pbc", "true"),
                ("p", "0.5"),
                ("spanning_from", "0,1,2"),
                ("spanning_to", "6,7,8"),
            ]),
            params(&[
                ("lattice_type", "honeycomb"),
                ("Lx", "4"),
                ("Ly", "2"),
                ("pbc", "true"),
                ("p", "0.5"),
                ("spanning_from", "0"),
                ("spanning_to", "1"),
            ]),
            params(&[
                ("lattice_type", "kagome"),
                ("Lx", "3"),
                ("Ly", "3"),
                ("pbc", "true"),
                ("p", "0.5"),
                ("spanning_from", "0"),
                ("spanning_to", "1"),
            ]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "4"),
                ("Ly", "4"),
                ("mode", "bond"),
                ("p", "0.3"),
            ]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "4"),
                ("Ly", "4"),
                ("mode", "site-bond"),
                ("p_site", "0.7"),
                ("p_bond", "0.4"),
            ]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "4"),
                ("Ly", "4"),
                ("p", "0.5"),
                ("observables", "largest-size,giant-fraction"),
            ]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "4"),
                ("Ly", "4"),
                ("p", "0.5"),
                ("spanning_from", "0,1,2,3"),
                ("spanning_to", "12,13,14,15"),
            ]),
            params(&[
                ("lattice_type", "square"),
                ("Lx", "3"),
                ("Ly", "3"),
                ("p", "0.5"),
                ("pbc", "true"),
                ("observables", "wrapping,winding-rank"),
            ]),
        ]
    }

    #[test]
    fn parameter_round_trip_is_stable_for_every_preset_shape() {
        for cases in preset_cases() {
            let spec = PercolationSpec::from_params(&cases)
                .unwrap_or_else(|error| panic!("{cases:?}: {error}"));
            let explicit = spec.to_params().expect("schema-expressible spec");
            let again = PercolationSpec::from_params(&explicit)
                .unwrap_or_else(|error| panic!("{explicit:?}: {error}"));
            assert_eq!(spec, again, "round trip changed the spec: {cases:?}");
            // A second serialization is a fixed point.
            assert_eq!(explicit, again.to_params().expect("still expressible"));
        }
    }

    #[test]
    fn spec_build_matches_the_parameter_preset() {
        for cases in preset_cases() {
            let spec = PercolationSpec::from_params(&cases).expect("spec");
            let built = spec.build().expect("built adapter");
            let preset =
                StaticPercolationMC::from_params(&cases, &mut Xoshiro256PlusPlus::seed_from_u64(1))
                    .expect("preset adapter");
            assert_eq!(
                built.lattice().vertex_count(),
                preset.lattice().vertex_count(),
                "{cases:?}"
            );
            assert_eq!(
                built.lattice().edge_count(),
                preset.lattice().edge_count(),
                "{cases:?}"
            );
            assert_eq!(built.law(), preset.law(), "{cases:?}");
            assert_eq!(
                built.observable_plan().observables(),
                preset.observable_plan().observables(),
                "{cases:?}"
            );
            assert_eq!(
                built.observable_plan().boundary_queries(),
                preset.observable_plan().boundary_queries(),
                "{cases:?}"
            );
            assert_eq!(
                built.embedding().is_some(),
                preset.embedding().is_some(),
                "{cases:?}"
            );
        }
    }

    #[test]
    fn provenance_is_versioned_and_round_trips_through_results_json() {
        let cases = params(&[
            ("lattice_type", "square"),
            ("Lx", "3"),
            ("Ly", "3"),
            ("pbc", "true"),
            ("p", "0.5"),
            ("observables", "wrapping,winding-rank"),
        ]);
        let spec = PercolationSpec::from_params(&cases).expect("spec");
        let provenance = spec.provenance();
        let json = provenance.to_json();
        assert_eq!(json["schema"], "cmc.percolation.provenance");
        assert_eq!(json["version"], 1);
        assert_eq!(json["topology"], "square/3x3/pbc");
        assert_eq!(json["law"], "site/p=0.5");
        assert_eq!(json["process"], "independent-occupation/v1");
        assert_eq!(json["observables"], json!(["Wrapping", "WindingRank"]));
        assert_eq!(json["boundary_queries"], 1);
        assert_eq!(json["estimator"], "sample-moments/v1");

        let mut results = Results::new();
        spec.stamp(&mut results);
        assert_eq!(results.metadata().attribute("percolation"), Some(&json));
        let serialized: Value =
            serde_json::from_str(&results.to_json().expect("results json")).expect("json");
        assert_eq!(serialized["metadata"]["attributes"]["percolation"], json);
    }

    #[test]
    fn parameter_resolution_rejections_are_typed() {
        assert!(matches!(
            PercolationSpec::from_params(&params(&[("lattice_type", "hexagon"), ("p", "0.5")])),
            Err(SpecError::UnknownLattice { .. })
        ));
        assert!(matches!(
            PercolationSpec::from_params(&params(&[
                ("lattice_type", "chain"),
                ("L", "0"),
                ("p", "0.5")
            ])),
            Err(SpecError::InvalidDimension { field: "L" })
        ));
        assert!(matches!(
            PercolationSpec::from_params(&params(&[("p", "NaN")])),
            Err(SpecError::Params(_))
        ));
        assert!(matches!(
            PercolationSpec::from_params(&params(&[
                ("lattice_type", "triangular"),
                ("Lx", "3"),
                ("Ly", "3"),
                ("pbc", "true"),
                ("p", "0.5")
            ])),
            Err(SpecError::Params(_))
        ));
        assert!(matches!(
            PercolationSpec::from_params(&params(&[
                ("lattice_type", "honeycomb"),
                ("Lx", "3"),
                ("Ly", "2"),
                ("pbc", "true"),
                ("p", "0.5")
            ])),
            Err(SpecError::InvalidTopology { .. })
        ));
        assert!(matches!(
            PercolationSpec::from_params(&params(&[
                ("lattice_type", "kagome"),
                ("Lx", "2"),
                ("Ly", "2"),
                ("p", "0.5")
            ])),
            Err(SpecError::PeriodicRequired { lattice: "kagome" })
        ));
    }

    #[test]
    fn programmatic_shape_rejections_are_typed() {
        let topology = TopologySpec::new(LatticeKind::Square { lx: 2, ly: 2 }, false);
        let law = StaticLaw::site(0.5).expect("law");
        let spec = |observables: Vec<StaticObservable>, boundaries: Option<BoundarySets>| {
            PercolationSpec {
                topology,
                law: law.clone(),
                process: ProcessSpec::IndependentOccupation,
                observables,
                boundaries,
                estimator: EstimatorSpec::SampleMoments,
            }
        };

        assert!(matches!(
            spec(
                vec![StaticObservable::LargestSize],
                Some(BoundarySets::new(vec![9], vec![0]))
            )
            .build(),
            Err(SpecError::Query(_))
        ));
        assert!(matches!(
            spec(vec![StaticObservable::Wrapping], None).build(),
            Err(SpecError::Adapter(
                StaticPercolationError::WindingWithoutEmbedding
            ))
        ));
        assert!(matches!(
            spec(Vec::new(), Some(BoundarySets::new(vec![0], vec![1]))).to_params(),
            Err(SpecError::EmptyObservables)
        ));
        assert!(matches!(
            spec(vec![StaticObservable::LargestSize], None).to_params(),
            Err(SpecError::NoBoundaryQuery)
        ));
        assert!(matches!(
            spec(
                vec![StaticObservable::LargestSize],
                Some(BoundarySets::new(Vec::new(), vec![0]))
            )
            .to_params(),
            Err(SpecError::EmptyBoundarySets)
        ));

        let heterogeneous = PercolationSpec {
            topology,
            law: StaticLaw::site_heterogeneous(vec![Probability::new(0.5).expect("p"); 3]),
            process: ProcessSpec::IndependentOccupation,
            observables: vec![StaticObservable::LargestSize],
            boundaries: Some(BoundarySets::new(vec![0], vec![1])),
            estimator: EstimatorSpec::SampleMoments,
        };
        assert!(matches!(
            heterogeneous.to_params(),
            Err(SpecError::LawNotParamsExpressible)
        ));
        assert!(matches!(
            heterogeneous.build(),
            Err(SpecError::Adapter(StaticPercolationError::LawDomain(_)))
        ));
    }

    #[test]
    fn heterogeneous_specs_build_with_exact_field_lengths() {
        let topology = TopologySpec::new(LatticeKind::Square { lx: 2, ly: 2 }, false);
        let spec = PercolationSpec {
            topology,
            law: StaticLaw::site_heterogeneous(vec![Probability::new(0.5).expect("p"); 4]),
            process: ProcessSpec::IndependentOccupation,
            observables: vec![StaticObservable::LargestSize],
            boundaries: Some(BoundarySets::new(vec![0], vec![1])),
            estimator: EstimatorSpec::SampleMoments,
        };
        let adapter = spec.build().expect("matching heterogeneous field");
        assert_eq!(adapter.lattice().vertex_count(), 4);
        assert_eq!(spec.provenance().law(), "site-heterogeneous/vertices=4");
    }
}
