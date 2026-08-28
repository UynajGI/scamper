//! Sparse operator catalogs for sign-problem-free quantum-spin models.
//!
//! A model is compiled into positive one-site and two-site matrix elements of
//! `K = C - H`. The continuous-time engine only sees these operator terms.
//! Physical model helpers (Heisenberg, XXZ, XYZ, transverse-field Ising) are
//! therefore thin catalog builders rather than separate Monte Carlo codes.

use std::collections::{HashMap, VecDeque};

use rand::Rng;
use rand::RngExt;

use crate::graph::CsrGraph;
use crate::local_space::{BasisState, LocalHilbertSpace, SpinSpace};

use super::error::LatticeQmcError;
use super::scattering::{ScatteringPolicy, ScatteringTable};
use super::vertex::VertexKind;

const COUPLING_TOLERANCE: f64 = 1.0e-14;

/// Physical anisotropic exchange on one graph edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeCoupling {
    /// Coefficient of `Sx_i Sx_j` in the physical Hamiltonian.
    pub j_x: f64,
    /// Coefficient of `Sy_i Sy_j`.
    pub j_y: f64,
    /// Coefficient of `Sz_i Sz_j`.
    pub j_z: f64,
    /// Optional explicit positive shift `C_b`. Automatic shifts are safer.
    pub shift: Option<f64>,
}

impl EdgeCoupling {
    /// Isotropic Heisenberg exchange.
    pub const fn heisenberg(j: f64) -> Self {
        Self {
            j_x: j,
            j_y: j,
            j_z: j,
            shift: None,
        }
    }

    /// XXZ exchange.
    pub const fn xxz(j_xy: f64, j_z: f64) -> Self {
        Self {
            j_x: j_xy,
            j_y: j_xy,
            j_z,
            shift: None,
        }
    }

    /// Fully anisotropic XYZ exchange.
    pub const fn xyz(j_x: f64, j_y: f64, j_z: f64) -> Self {
        Self {
            j_x,
            j_y,
            j_z,
            shift: None,
        }
    }

    fn is_zero(self) -> bool {
        self.j_x.abs() <= COUPLING_TOLERANCE
            && self.j_y.abs() <= COUPLING_TOLERANCE
            && self.j_z.abs() <= COUPLING_TOLERANCE
    }
}

impl Default for EdgeCoupling {
    fn default() -> Self {
        Self::xyz(0.0, 0.0, 0.0)
    }
}

/// Physical one-site terms `D Sz^2 - h_z Sz - h_x Sx`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteCoupling {
    /// Transverse field in `-h_x Sx`.
    pub h_x: f64,
    /// Longitudinal field in `-h_z Sz`.
    pub h_z: f64,
    /// Single-ion anisotropy in `D Sz^2`.
    pub single_ion: f64,
    /// Optional explicit positive shift.
    pub shift: Option<f64>,
}

impl SiteCoupling {
    /// Construct site terms.
    pub const fn new(h_x: f64, h_z: f64, single_ion: f64) -> Self {
        Self {
            h_x,
            h_z,
            single_ion,
            shift: None,
        }
    }

    fn is_zero(self) -> bool {
        self.h_x.abs() <= COUPLING_TOLERANCE
            && self.h_z.abs() <= COUPLING_TOLERANCE
            && self.single_ion.abs() <= COUPLING_TOLERANCE
    }
}

impl Default for SiteCoupling {
    fn default() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

/// Handling of basis phases used to make off-diagonal matrix elements positive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaugePolicy {
    /// Solve the graph-wide `Z2` Marshall gauge when possible.
    Auto,
    /// Keep the input `Sz` product basis and reject negative matrix elements.
    Identity,
}

/// Location of one local operator term.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TermLocation {
    /// Unique graph edge index.
    Edge(usize),
    /// Site index.
    Site(usize),
}

/// One sparse positive local operator `K_t`.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatorTerm {
    label: String,
    location: TermLocation,
    sites: Box<[usize]>,
    kinds: Vec<VertexKind>,
    diagonal_lookup: HashMap<Vec<BasisState>, usize>,
    kind_lookup: HashMap<Vec<BasisState>, usize>,
    scattering: ScatteringTable,
    proposal_weight: f64,
    shift: f64,
}

impl OperatorTerm {
    fn new(
        label: impl Into<String>,
        location: TermLocation,
        sites: Vec<usize>,
        kinds: Vec<VertexKind>,
        space: &SpinSpace,
        scattering_policy: ScatteringPolicy,
        shift: f64,
    ) -> Result<Self, LatticeQmcError> {
        let dimensions = sites.iter().map(|&site| space.dimension(site)).collect();
        Self::from_sparse(
            label,
            location,
            sites,
            dimensions,
            kinds,
            scattering_policy,
            shift,
        )
    }

    /// Construct a custom positive one- or two-site sparse operator.
    ///
    /// This constructor is local-space agnostic and is the extension point for
    /// future bosonic and fermionic model catalogs. Fermionic callers must not
    /// use the positive engine until their global sign policy is implemented.
    #[allow(clippy::too_many_arguments)]
    pub fn from_sparse(
        label: impl Into<String>,
        location: TermLocation,
        sites: Vec<usize>,
        local_dimensions: Vec<usize>,
        kinds: Vec<VertexKind>,
        scattering_policy: ScatteringPolicy,
        shift: f64,
    ) -> Result<Self, LatticeQmcError> {
        if sites.is_empty() || sites.len() > 2 {
            return Err(LatticeQmcError::InvalidModel(
                "the current lattice backend supports one- and two-site terms".into(),
            ));
        }
        if local_dimensions.len() != sites.len()
            || local_dimensions.iter().any(|&dimension| dimension < 2)
        {
            return Err(LatticeQmcError::InvalidModel(
                "one local dimension is required for every term site".into(),
            ));
        }
        let leg_dimensions: Vec<usize> = local_dimensions
            .iter()
            .flat_map(|&dimension| [dimension, dimension])
            .collect();
        let mut diagonal_lookup = HashMap::new();
        let mut kind_lookup = HashMap::new();
        let mut proposal_weight = 0.0_f64;
        for (kind_id, kind) in kinds.iter().enumerate() {
            if kind.legs().len() != 2 * sites.len() {
                return Err(LatticeQmcError::InvalidModel(format!(
                    "term kind `{}` has the wrong leg count",
                    kind.name()
                )));
            }
            if kind_lookup.insert(kind.legs().to_vec(), kind_id).is_some() {
                return Err(LatticeQmcError::InvalidModel(format!(
                    "term `{}` has duplicate matrix elements",
                    kind.name()
                )));
            }
            if kind.is_diagonal() {
                let states: Vec<_> = kind
                    .legs()
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|pair| pair[0])
                    .collect();
                diagonal_lookup.insert(states, kind_id);
                proposal_weight = proposal_weight.max(kind.weight());
            }
        }
        if diagonal_lookup.is_empty() || proposal_weight <= 0.0 {
            return Err(LatticeQmcError::InvalidModel(
                "every continuous-time term needs a positive diagonal seed".into(),
            ));
        }
        let scattering = ScatteringTable::build(&kinds, &leg_dimensions, scattering_policy)?;
        Ok(Self {
            label: label.into(),
            location,
            sites: sites.into_boxed_slice(),
            kinds,
            diagonal_lookup,
            kind_lookup,
            scattering,
            proposal_weight,
            shift,
        })
    }

    /// Human-readable term label.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Graph location.
    pub fn location(&self) -> TermLocation {
        self.location
    }

    /// Sites acted on by the term.
    pub fn sites(&self) -> &[usize] {
        &self.sites
    }

    /// Local matrix-element catalog.
    pub fn kinds(&self) -> &[VertexKind] {
        &self.kinds
    }

    /// One matrix-element kind.
    pub fn kind(&self, kind: usize) -> &VertexKind {
        &self.kinds[kind]
    }

    /// Find a diagonal seed for current local states.
    pub fn diagonal_kind(&self, states: &[BasisState]) -> Option<usize> {
        self.diagonal_lookup.get(states).copied()
    }

    /// Find a kind by its full leg pattern.
    pub fn kind_for_legs(&self, legs: &[BasisState]) -> Option<usize> {
        self.kind_lookup.get(legs).copied()
    }

    /// Directed-loop scattering table.
    pub fn scattering(&self) -> &ScatteringTable {
        &self.scattering
    }

    /// Importance weight used to propose this term.
    pub fn proposal_weight(&self) -> f64 {
        self.proposal_weight
    }

    /// Constant shift contributed by this term.
    pub fn shift(&self) -> f64 {
        self.shift
    }
}

/// Positive sparse-operator model consumed by the continuous-time engine.
///
/// The trait deliberately contains no spin algebra. A future bosonic or
/// fermionic backend can provide another local space and operator catalog;
/// fermions additionally need a sign/determinant policy before they can use a
/// positive-weight engine.
pub trait PositiveOperatorModel: Send + Sync {
    /// Local Hilbert-space implementation.
    type Space: LocalHilbertSpace;

    /// Model label.
    fn name(&self) -> &str;
    /// Arbitrary graph topology.
    fn graph(&self) -> &CsrGraph;
    /// Local Hilbert space.
    fn space(&self) -> &Self::Space;
    /// Sparse local terms.
    fn terms(&self) -> &[OperatorTerm];
    /// One term.
    fn term(&self, term: usize) -> &OperatorTerm;
    /// Number of terms.
    fn term_count(&self) -> usize;
    /// Constant shift in `H = C - sum K_t`.
    fn constant_shift(&self) -> f64;
    /// Proposal probability for one term.
    fn term_probability(&self, term: usize) -> f64;
    /// Draw one term from the model proposal distribution.
    fn sample_term<R: Rng + ?Sized>(&self, rng: &mut R) -> Option<usize>;
}

/// Compiled sign-problem-free spin-lattice model.
#[derive(Debug, Clone, PartialEq)]
pub struct SpinLatticeModel {
    name: String,
    graph: CsrGraph,
    space: SpinSpace,
    terms: Vec<OperatorTerm>,
    proposal_cdf: Vec<f64>,
    proposal_probabilities: Vec<f64>,
    constant_shift: f64,
    gauge: Vec<i8>,
}

impl SpinLatticeModel {
    /// Isotropic Heisenberg helper.
    pub fn heisenberg(graph: CsrGraph, two_s: u16, j: f64) -> Result<Self, LatticeQmcError> {
        let space = SpinSpace::uniform(graph.site_count(), two_s)?;
        SpinModelBuilder::new(graph, space)
            .name("Heisenberg")
            .uniform_edge(EdgeCoupling::heisenberg(j))
            .build()
    }

    /// XXZ helper.
    pub fn xxz(graph: CsrGraph, two_s: u16, j_xy: f64, j_z: f64) -> Result<Self, LatticeQmcError> {
        let space = SpinSpace::uniform(graph.site_count(), two_s)?;
        SpinModelBuilder::new(graph, space)
            .name("XXZ")
            .uniform_edge(EdgeCoupling::xxz(j_xy, j_z))
            .build()
    }

    /// XYZ helper.
    pub fn xyz(
        graph: CsrGraph,
        two_s: u16,
        coupling: EdgeCoupling,
    ) -> Result<Self, LatticeQmcError> {
        let space = SpinSpace::uniform(graph.site_count(), two_s)?;
        SpinModelBuilder::new(graph, space)
            .name("XYZ")
            .uniform_edge(coupling)
            .build()
    }

    /// Transverse-field Ising helper.
    pub fn transverse_field_ising(
        graph: CsrGraph,
        two_s: u16,
        j_z: f64,
        h_x: f64,
    ) -> Result<Self, LatticeQmcError> {
        let space = SpinSpace::uniform(graph.site_count(), two_s)?;
        SpinModelBuilder::new(graph, space)
            .name("TransverseFieldIsing")
            .uniform_edge(EdgeCoupling::xxz(0.0, j_z))
            .uniform_site(SiteCoupling::new(h_x, 0.0, 0.0))
            .build()
    }

    /// Model name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Underlying arbitrary graph.
    pub fn graph(&self) -> &CsrGraph {
        &self.graph
    }

    /// Site-resolved spin space.
    pub fn space(&self) -> &SpinSpace {
        &self.space
    }

    /// Sparse local terms.
    pub fn terms(&self) -> &[OperatorTerm] {
        &self.terms
    }

    /// One term.
    pub fn term(&self, term: usize) -> &OperatorTerm {
        &self.terms[term]
    }

    /// Number of terms.
    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    /// Total `C` in `H = C - sum_t K_t`.
    pub fn constant_shift(&self) -> f64 {
        self.constant_shift
    }

    /// Marshall-gauge phase on each site.
    pub fn gauge(&self) -> &[i8] {
        &self.gauge
    }

    /// Proposal probability for a term.
    pub fn term_probability(&self, term: usize) -> f64 {
        self.proposal_probabilities[term]
    }

    /// Importance-sample one operator term.
    pub fn sample_term<R: Rng + ?Sized>(&self, rng: &mut R) -> Option<usize> {
        if self.terms.is_empty() {
            return None;
        }
        let u = rng.random::<f64>();
        let index = self.proposal_cdf.partition_point(|&value| value <= u);
        Some(index.min(self.terms.len() - 1))
    }
}

impl PositiveOperatorModel for SpinLatticeModel {
    type Space = SpinSpace;

    fn name(&self) -> &str {
        SpinLatticeModel::name(self)
    }

    fn graph(&self) -> &CsrGraph {
        SpinLatticeModel::graph(self)
    }

    fn space(&self) -> &Self::Space {
        SpinLatticeModel::space(self)
    }

    fn terms(&self) -> &[OperatorTerm] {
        SpinLatticeModel::terms(self)
    }

    fn term(&self, term: usize) -> &OperatorTerm {
        SpinLatticeModel::term(self, term)
    }

    fn term_count(&self) -> usize {
        SpinLatticeModel::term_count(self)
    }

    fn constant_shift(&self) -> f64 {
        SpinLatticeModel::constant_shift(self)
    }

    fn term_probability(&self, term: usize) -> f64 {
        SpinLatticeModel::term_probability(self, term)
    }

    fn sample_term<R: Rng + ?Sized>(&self, rng: &mut R) -> Option<usize> {
        SpinLatticeModel::sample_term(self, rng)
    }
}

/// Builder for site-dependent arbitrary-spin Hamiltonians.
#[derive(Debug, Clone)]
pub struct SpinModelBuilder {
    name: String,
    graph: CsrGraph,
    space: SpinSpace,
    edge_couplings: Vec<EdgeCoupling>,
    site_couplings: Vec<SiteCoupling>,
    gauge_policy: GaugePolicy,
    scattering_policy: ScatteringPolicy,
    shift_margin: f64,
}

impl SpinModelBuilder {
    /// Start from an arbitrary graph and site-resolved spin space.
    pub fn new(graph: CsrGraph, space: SpinSpace) -> Self {
        let edge_count = graph.edge_count();
        let site_count = graph.site_count();
        Self {
            name: "CustomSpinLattice".into(),
            graph,
            space,
            edge_couplings: vec![EdgeCoupling::default(); edge_count],
            site_couplings: vec![SiteCoupling::default(); site_count],
            gauge_policy: GaugePolicy::Auto,
            scattering_policy: ScatteringPolicy::LowBounce,
            shift_margin: 1.0e-10,
        }
    }

    /// Set a model label.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Apply one exchange coupling to every graph edge.
    pub fn uniform_edge(mut self, coupling: EdgeCoupling) -> Self {
        self.edge_couplings.fill(coupling);
        self
    }

    /// Set edge-resolved exchange couplings.
    pub fn edge_couplings(mut self, couplings: Vec<EdgeCoupling>) -> Result<Self, LatticeQmcError> {
        if couplings.len() != self.graph.edge_count() {
            return Err(LatticeQmcError::parameter(
                "edge_couplings",
                format!(
                    "expected {}, got {}",
                    self.graph.edge_count(),
                    couplings.len()
                ),
            ));
        }
        self.edge_couplings = couplings;
        Ok(self)
    }

    /// Apply one-site terms uniformly.
    pub fn uniform_site(mut self, coupling: SiteCoupling) -> Self {
        self.site_couplings.fill(coupling);
        self
    }

    /// Set site-resolved one-site terms.
    pub fn site_couplings(mut self, couplings: Vec<SiteCoupling>) -> Result<Self, LatticeQmcError> {
        if couplings.len() != self.graph.site_count() {
            return Err(LatticeQmcError::parameter(
                "site_couplings",
                format!(
                    "expected {}, got {}",
                    self.graph.site_count(),
                    couplings.len()
                ),
            ));
        }
        self.site_couplings = couplings;
        Ok(self)
    }

    /// Select the sign/gauge policy.
    pub fn gauge_policy(mut self, policy: GaugePolicy) -> Self {
        self.gauge_policy = policy;
        self
    }

    /// Select local scattering construction.
    pub fn scattering_policy(mut self, policy: ScatteringPolicy) -> Self {
        self.scattering_policy = policy;
        self
    }

    /// Set the strictly-positive automatic shift margin.
    pub fn shift_margin(mut self, margin: f64) -> Self {
        self.shift_margin = margin;
        self
    }

    /// Validate stoquasticity, solve the gauge, and compile sparse operators.
    pub fn build(self) -> Result<SpinLatticeModel, LatticeQmcError> {
        self.space.require_site_count(self.graph.site_count())?;
        if !self.shift_margin.is_finite() || self.shift_margin <= 0.0 {
            return Err(LatticeQmcError::parameter(
                "shift_margin",
                "must be finite and positive",
            ));
        }
        validate_finite_couplings(&self.edge_couplings, &self.site_couplings)?;
        let gauge = solve_gauge(
            &self.graph,
            &self.edge_couplings,
            &self.site_couplings,
            self.gauge_policy,
        )?;
        let mut terms = Vec::new();
        for (edge_id, coupling) in self.edge_couplings.iter().copied().enumerate() {
            if coupling.is_zero() {
                continue;
            }
            terms.push(build_edge_term(
                edge_id,
                coupling,
                &self.graph,
                &self.space,
                &gauge,
                self.scattering_policy,
                self.shift_margin,
            )?);
        }
        for (site, coupling) in self.site_couplings.iter().copied().enumerate() {
            if coupling.is_zero() {
                continue;
            }
            terms.push(build_site_term(
                site,
                coupling,
                &self.space,
                &gauge,
                self.scattering_policy,
                self.shift_margin,
            )?);
        }
        let constant_shift = terms.iter().map(OperatorTerm::shift).sum();
        let total_proposal_weight: f64 = terms.iter().map(OperatorTerm::proposal_weight).sum();
        let mut proposal_probabilities = Vec::with_capacity(terms.len());
        let mut proposal_cdf = Vec::with_capacity(terms.len());
        let mut cumulative = 0.0;
        for term in &terms {
            let probability = term.proposal_weight() / total_proposal_weight;
            proposal_probabilities.push(probability);
            cumulative += probability;
            proposal_cdf.push(cumulative);
        }
        if let Some(last) = proposal_cdf.last_mut() {
            *last = 1.0;
        }
        Ok(SpinLatticeModel {
            name: self.name,
            graph: self.graph,
            space: self.space,
            terms,
            proposal_cdf,
            proposal_probabilities,
            constant_shift,
            gauge,
        })
    }
}

fn validate_finite_couplings(
    edges: &[EdgeCoupling],
    sites: &[SiteCoupling],
) -> Result<(), LatticeQmcError> {
    for coupling in edges {
        for value in [coupling.j_x, coupling.j_y, coupling.j_z] {
            if !value.is_finite() {
                return Err(LatticeQmcError::InvalidModel(
                    "edge couplings must be finite".into(),
                ));
            }
        }
    }
    for coupling in sites {
        for value in [coupling.h_x, coupling.h_z, coupling.single_ion] {
            if !value.is_finite() {
                return Err(LatticeQmcError::InvalidModel(
                    "site couplings must be finite".into(),
                ));
            }
        }
    }
    Ok(())
}

fn sign(value: f64) -> i8 {
    if value >= 0.0 {
        1
    } else {
        -1
    }
}

fn solve_gauge(
    graph: &CsrGraph,
    edge_couplings: &[EdgeCoupling],
    site_couplings: &[SiteCoupling],
    policy: GaugePolicy,
) -> Result<Vec<i8>, LatticeQmcError> {
    let mut edge_constraints = vec![None; graph.edge_count()];
    for (edge_id, coupling) in edge_couplings.iter().copied().enumerate() {
        let scale = graph.edge(edge_id).weight;
        let exchange = 0.25 * scale * (coupling.j_x + coupling.j_y);
        let pair = 0.25 * scale * (coupling.j_x - coupling.j_y);
        let mut required = None;
        for coefficient in [exchange, pair] {
            if coefficient.abs() <= COUPLING_TOLERANCE {
                continue;
            }
            let constraint = -sign(coefficient);
            if required.is_some_and(|old| old != constraint) {
                return Err(LatticeQmcError::InvalidModel(format!(
                    "edge {edge_id} has exchange and pair-flip amplitudes with incompatible signs"
                )));
            }
            required = Some(constraint);
        }
        edge_constraints[edge_id] = required;
    }
    let site_constraints: Vec<Option<i8>> = site_couplings
        .iter()
        .map(|coupling| (coupling.h_x.abs() > COUPLING_TOLERANCE).then(|| sign(coupling.h_x)))
        .collect();

    if policy == GaugePolicy::Identity {
        let gauge = vec![1_i8; graph.site_count()];
        validate_gauge(graph, &gauge, &edge_constraints, &site_constraints)?;
        return Ok(gauge);
    }

    let mut relative = vec![0_i8; graph.site_count()];
    let mut gauge = vec![1_i8; graph.site_count()];
    for root in 0..graph.site_count() {
        if relative[root] != 0 {
            continue;
        }
        relative[root] = 1;
        let mut component = Vec::new();
        let mut queue = VecDeque::from([root]);
        while let Some(site) = queue.pop_front() {
            component.push(site);
            for neighbor in graph.neighbors(site) {
                let Some(constraint) = edge_constraints[neighbor.edge] else {
                    continue;
                };
                let expected = constraint * relative[site];
                if relative[neighbor.site] == 0 {
                    relative[neighbor.site] = expected;
                    queue.push_back(neighbor.site);
                } else if relative[neighbor.site] != expected {
                    return Err(LatticeQmcError::InvalidModel(
                        "off-diagonal signs are frustrated on the supplied graph".into(),
                    ));
                }
            }
        }
        let mut root_sign = None;
        for &site in &component {
            if let Some(required_site_phase) = site_constraints[site] {
                let required_root = required_site_phase * relative[site];
                if root_sign.is_some_and(|old| old != required_root) {
                    return Err(LatticeQmcError::InvalidModel(
                        "transverse fields conflict with the Marshall gauge".into(),
                    ));
                }
                root_sign = Some(required_root);
            }
        }
        let root_sign = root_sign.unwrap_or(1);
        for site in component {
            gauge[site] = root_sign * relative[site];
        }
    }
    validate_gauge(graph, &gauge, &edge_constraints, &site_constraints)?;
    Ok(gauge)
}

fn validate_gauge(
    graph: &CsrGraph,
    gauge: &[i8],
    edge_constraints: &[Option<i8>],
    site_constraints: &[Option<i8>],
) -> Result<(), LatticeQmcError> {
    for (site, required) in site_constraints.iter().enumerate() {
        if required.is_some_and(|phase| gauge[site] != phase) {
            return Err(LatticeQmcError::InvalidModel(format!(
                "site {site} requires gauge phase {required:?}"
            )));
        }
    }
    for (edge_id, required) in edge_constraints.iter().enumerate() {
        if let Some(required_product) = required {
            let edge = graph.edge(edge_id);
            if gauge[edge.source] * gauge[edge.target] != *required_product {
                return Err(LatticeQmcError::InvalidModel(format!(
                    "edge {edge_id} requires gauge product {required_product}"
                )));
            }
        }
    }
    Ok(())
}

fn positive_shift(
    maximum_diagonal: f64,
    requested: Option<f64>,
    margin: f64,
) -> Result<f64, LatticeQmcError> {
    let automatic = maximum_diagonal + margin * maximum_diagonal.abs().max(1.0);
    let shift = requested.unwrap_or(automatic);
    if !shift.is_finite() || shift <= maximum_diagonal {
        return Err(LatticeQmcError::InvalidModel(format!(
            "operator shift {shift} must exceed the maximum diagonal energy {maximum_diagonal}"
        )));
    }
    Ok(shift)
}

#[allow(clippy::too_many_arguments)]
fn build_edge_term(
    edge_id: usize,
    coupling: EdgeCoupling,
    graph: &CsrGraph,
    space: &SpinSpace,
    gauge: &[i8],
    scattering_policy: ScatteringPolicy,
    shift_margin: f64,
) -> Result<OperatorTerm, LatticeQmcError> {
    let edge = graph.edge(edge_id);
    let site_i = edge.source;
    let site_j = edge.target;
    let scale = edge.weight;
    let phase = f64::from(gauge[site_i] * gauge[site_j]);
    let exchange = -0.25 * scale * (coupling.j_x + coupling.j_y) * phase;
    let pair = -0.25 * scale * (coupling.j_x - coupling.j_y) * phase;
    if exchange < -COUPLING_TOLERANCE || pair < -COUPLING_TOLERANCE {
        return Err(LatticeQmcError::InvalidModel(format!(
            "edge {edge_id} remains non-stoquastic after gauge transformation"
        )));
    }
    let dimension_i = space.dimension(site_i);
    let dimension_j = space.dimension(site_j);
    let mut maximum_diagonal = f64::NEG_INFINITY;
    for raw_i in 0..dimension_i {
        let state_i = raw_i as BasisState;
        for raw_j in 0..dimension_j {
            let state_j = raw_j as BasisState;
            maximum_diagonal = maximum_diagonal
                .max(scale * coupling.j_z * space.m(site_i, state_i) * space.m(site_j, state_j));
        }
    }
    let shift = positive_shift(maximum_diagonal, coupling.shift, shift_margin)?;
    let mut kinds = Vec::new();
    for raw_i in 0..dimension_i {
        let state_i = raw_i as BasisState;
        for raw_j in 0..dimension_j {
            let state_j = raw_j as BasisState;
            let m_i = space.m(site_i, state_i);
            let m_j = space.m(site_j, state_j);
            let diagonal_weight = shift - scale * coupling.j_z * m_i * m_j;
            kinds.push(VertexKind::new(
                format!("edge-{edge_id}-diag-{state_i}-{state_j}"),
                vec![state_i, state_i, state_j, state_j],
                diagonal_weight,
            )?);

            if exchange > COUPLING_TOLERANCE {
                if let (Some(raise_i), Some(lower_j)) = (
                    space.raising_amplitude(site_i, state_i),
                    space.lowering_amplitude(site_j, state_j),
                ) {
                    kinds.push(VertexKind::new(
                        format!("edge-{edge_id}-exchange-plus-minus-{state_i}-{state_j}"),
                        vec![state_i, state_i + 1, state_j, state_j - 1],
                        exchange * raise_i * lower_j,
                    )?);
                }
                if let (Some(lower_i), Some(raise_j)) = (
                    space.lowering_amplitude(site_i, state_i),
                    space.raising_amplitude(site_j, state_j),
                ) {
                    kinds.push(VertexKind::new(
                        format!("edge-{edge_id}-exchange-minus-plus-{state_i}-{state_j}"),
                        vec![state_i, state_i - 1, state_j, state_j + 1],
                        exchange * lower_i * raise_j,
                    )?);
                }
            }
            if pair > COUPLING_TOLERANCE {
                if let (Some(raise_i), Some(raise_j)) = (
                    space.raising_amplitude(site_i, state_i),
                    space.raising_amplitude(site_j, state_j),
                ) {
                    kinds.push(VertexKind::new(
                        format!("edge-{edge_id}-pair-plus-plus-{state_i}-{state_j}"),
                        vec![state_i, state_i + 1, state_j, state_j + 1],
                        pair * raise_i * raise_j,
                    )?);
                }
                if let (Some(lower_i), Some(lower_j)) = (
                    space.lowering_amplitude(site_i, state_i),
                    space.lowering_amplitude(site_j, state_j),
                ) {
                    kinds.push(VertexKind::new(
                        format!("edge-{edge_id}-pair-minus-minus-{state_i}-{state_j}"),
                        vec![state_i, state_i - 1, state_j, state_j - 1],
                        pair * lower_i * lower_j,
                    )?);
                }
            }
        }
    }
    OperatorTerm::new(
        format!("edge-{edge_id}"),
        TermLocation::Edge(edge_id),
        vec![site_i, site_j],
        kinds,
        space,
        scattering_policy,
        shift,
    )
}

fn build_site_term(
    site: usize,
    coupling: SiteCoupling,
    space: &SpinSpace,
    gauge: &[i8],
    scattering_policy: ScatteringPolicy,
    shift_margin: f64,
) -> Result<OperatorTerm, LatticeQmcError> {
    let dimension = space.dimension(site);
    let mut maximum_diagonal = f64::NEG_INFINITY;
    for raw_state in 0..dimension {
        let state = raw_state as BasisState;
        let m = space.m(site, state);
        maximum_diagonal = maximum_diagonal.max(coupling.single_ion * m * m - coupling.h_z * m);
    }
    let shift = positive_shift(maximum_diagonal, coupling.shift, shift_margin)?;
    let transverse = 0.5 * coupling.h_x * f64::from(gauge[site]);
    if transverse < -COUPLING_TOLERANCE {
        return Err(LatticeQmcError::InvalidModel(format!(
            "site {site} transverse field remains non-stoquastic"
        )));
    }
    let mut kinds = Vec::new();
    for raw_state in 0..dimension {
        let state = raw_state as BasisState;
        let m = space.m(site, state);
        kinds.push(VertexKind::new(
            format!("site-{site}-diag-{state}"),
            vec![state, state],
            shift - coupling.single_ion * m * m + coupling.h_z * m,
        )?);
        if transverse > COUPLING_TOLERANCE {
            if let Some(amplitude) = space.raising_amplitude(site, state) {
                kinds.push(VertexKind::new(
                    format!("site-{site}-raise-{state}"),
                    vec![state, state + 1],
                    transverse * amplitude,
                )?);
            }
            if let Some(amplitude) = space.lowering_amplitude(site, state) {
                kinds.push(VertexKind::new(
                    format!("site-{site}-lower-{state}"),
                    vec![state, state - 1],
                    transverse * amplitude,
                )?);
            }
        }
    }
    OperatorTerm::new(
        format!("site-{site}"),
        TermLocation::Site(site),
        vec![site],
        kinds,
        space,
        scattering_policy,
        shift,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{CsrGraph, EdgeSpec};

    #[test]
    fn spin_two_xxz_catalog_is_positive() {
        let graph = CsrGraph::chain(4, true).expect("graph");
        let model = SpinLatticeModel::xxz(graph, 4, 1.0, 0.7).expect("model");
        assert!(!model.terms().is_empty());
        for term in model.terms() {
            for kind in term.kinds() {
                assert!(kind.weight() > 0.0);
            }
            assert!(term.scattering().diagnostics().max_detailed_balance_error < 1.0e-10);
        }
    }

    #[test]
    fn antiferromagnet_rejects_odd_cycle() {
        let graph = CsrGraph::from_edges(
            3,
            [
                EdgeSpec::new(0, 1),
                EdgeSpec::new(1, 2),
                EdgeSpec::new(2, 0),
            ],
        )
        .expect("triangle");
        let error = SpinLatticeModel::heisenberg(graph, 1, 1.0).expect_err("sign problem");
        assert!(error.to_string().contains("frustrated"));
    }

    #[test]
    fn transverse_field_ising_has_site_vertices() {
        let graph = CsrGraph::chain(3, false).expect("graph");
        let model = SpinLatticeModel::transverse_field_ising(graph, 1, -1.0, 0.5).expect("model");
        assert!(model
            .terms()
            .iter()
            .any(|term| matches!(term.location(), TermLocation::Site(_))));
    }
}
