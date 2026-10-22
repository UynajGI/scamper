//! D3: layered space-time directed percolation (absorbing-state process).
//!
//! One sweep generates one complete space-time realization over `T`
//! explicit causal layers and records it as ONE C0 vector sample — never
//! as T independent measurements. State lives in two frontier bitsets
//! (`u64` words), so memory scales with layer width, not space-time
//! volume; per-layer activity counts and survival indicators are the only
//! history retained.
//!
//! Physics: directed bond percolation. Layer `t+1` site `j` is active iff
//! some spatial neighbor `i` is active in layer `t` and the directed bond
//! `i -> j` is open at time `t`; bonds are independent Bernoulli(p) per
//! (time step, directed incidence). Anchored to Hinrichsen, Adv. Phys. 49,
//! 815 (2000) (arXiv:cond-mat/0001070, archived in `docs/ref/`):
//! (1+1)-dimensional bond DP threshold `p_c = 0.6447001(1)` (Table 1) and
//! survival/density decay exponent `delta = 0.159464(6)` (Table 2).
//!
//! RNG draw contract (bitwise reproducibility, and the contract the
//! explicit-DAG parity gate replays): per layer, iterate active sites in
//! ascending dense index order; for each, iterate its spatial incidences
//! in CSR order and draw one `f64` per incidence; the bond is open when
//! the draw is `< p`. Bonds out of inactive sites are never drawn (their
//! openness cannot influence the process law).

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::percolation::Probability;
use crate::{CsrLattice, GraphView, UndirectedGraphView, VertexId};

/// Initial condition of a space-time realization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DpSeedPolicy {
    /// A single active seed at the given dense site index (survival
    /// experiments: `P(t) ~ t^-delta` at criticality).
    Single { site: usize },
    /// Every site active in layer 0 (density-decay experiments:
    /// `rho(t) ~ t^-delta` at criticality).
    Full,
}

/// Owned adapter producing layered directed-percolation realizations.
///
/// `layers` counts the recorded layers including `t = 0`. Configure
/// `thermalization_sweeps = 0`; realizations are i.i.d.
pub struct LayeredDirectedPercolationMC {
    space: CsrLattice,
    probability: Probability,
    layers: usize,
    seed: DpSeedPolicy,
    /// Double-buffered frontier bitsets (ceil(L/64) words each).
    frontier: Vec<u64>,
    next: Vec<u64>,
    /// History of the last realization: active count per layer.
    layer_active: Vec<f64>,
    /// Survival indicator per layer of the last realization.
    layer_survival: Vec<f64>,
    /// First layer with zero activity (`layers` when the realization
    /// survives to the end).
    extinction_time: usize,
}

impl LayeredDirectedPercolationMC {
    /// Build the process over an owned spatial lattice.
    pub fn new(
        space: CsrLattice,
        keep: f64,
        layers: usize,
        seed: DpSeedPolicy,
    ) -> Result<Self, DpProcessError> {
        space
            .validate()
            .map_err(|error| DpProcessError::Topology(error.to_string()))?;
        if layers == 0 {
            return Err(DpProcessError::NoLayers);
        }
        if let DpSeedPolicy::Single { site } = seed {
            if site >= space.vertex_count() {
                return Err(DpProcessError::SeedDomain {
                    site,
                    vertices: space.vertex_count(),
                });
            }
        }
        let probability = Probability::new(keep).map_err(DpProcessError::Probability)?;
        let words = space.vertex_count().div_ceil(64);
        Ok(Self {
            space,
            probability,
            layers,
            seed,
            frontier: vec![0; words],
            next: vec![0; words],
            layer_active: vec![0.0; layers],
            layer_survival: vec![0.0; layers],
            extinction_time: layers,
        })
    }

    /// The spatial lattice the process runs on.
    pub const fn space(&self) -> &CsrLattice {
        &self.space
    }

    /// Open-bond probability.
    pub const fn probability(&self) -> f64 {
        self.probability.get()
    }

    /// Recorded layer count (including `t = 0`).
    pub const fn layer_count(&self) -> usize {
        self.layers
    }

    /// The initial condition.
    pub const fn seed_policy(&self) -> DpSeedPolicy {
        self.seed
    }

    /// Active-site count per layer of the last realization.
    pub fn layer_active_counts(&self) -> &[f64] {
        &self.layer_active
    }

    /// Survival indicator per layer of the last realization.
    pub fn layer_survival(&self) -> &[f64] {
        &self.layer_survival
    }

    /// First layer with zero activity (`layers` when it survives).
    pub const fn extinction_time(&self) -> usize {
        self.extinction_time
    }

    /// Draw one complete realization into the history buffers.
    ///
    /// Exposed for tests replaying the documented draw contract against an
    /// explicit layered DAG; `measure` records what it produced.
    pub fn draw_realization<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let probability = self.probability.get();
        self.frontier.fill(0);
        match self.seed {
            DpSeedPolicy::Single { site } => self.frontier[site / 64] |= 1 << (site % 64),
            DpSeedPolicy::Full => self.frontier.fill(!0u64),
        }
        if self.seed == DpSeedPolicy::Full {
            // Mask trailing bits beyond the vertex count in the last word.
            let vertices = self.space.vertex_count();
            let remainder = vertices % 64;
            if remainder != 0 {
                let last = self.frontier.len() - 1;
                self.frontier[last] &= (1u64 << remainder) - 1;
            }
        }
        self.extinction_time = self.layers;
        let mut extinct = false;
        for time in 0..self.layers {
            let mut active = 0usize;
            let mut any = false;
            for word in &self.frontier {
                active += word.count_ones() as usize;
                any |= *word != 0;
            }
            self.layer_active[time] = active as f64;
            self.layer_survival[time] = f64::from(u8::from(any));
            if !any && !extinct {
                extinct = true;
                self.extinction_time = time;
            }
            if time + 1 == self.layers {
                break;
            }
            // Propagate: ascending active sites, CSR incidence order, one
            // Bernoulli draw per incidence (the documented contract).
            self.next.fill(0);
            for (word_index, &word) in self.frontier.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let bit = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let site = word_index * 64 + bit;
                    let identity = VertexId::from_valid_index(site);
                    for incidence in UndirectedGraphView::incidences(&self.space, identity) {
                        if rng.random::<f64>() < probability {
                            let neighbor = incidence.neighbor.index();
                            self.next[neighbor / 64] |= 1 << (neighbor % 64);
                        }
                    }
                }
            }
            std::mem::swap(&mut self.frontier, &mut self.next);
        }
    }
}

impl MonteCarlo for LayeredDirectedPercolationMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.draw_realization(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure_array("DpLayerActive", &self.layer_active)
            .expect("layer count is fixed at construction");
        ctx.measure_array("DpSurvival", &self.layer_survival)
            .expect("layer count is fixed at construction");
        ctx.measure("DpExtinctionTime", self.extinction_time as f64);
    }

    fn name(&self) -> &'static str {
        "LayeredDirectedPercolation"
    }
}

/// Typed rejections for the layered process.
#[derive(Debug)]
pub enum DpProcessError {
    /// Zero requested layers.
    NoLayers,
    /// Keep probability outside `[0, 1]`.
    Probability(crate::percolation::ProbabilityError),
    /// Single-seed site outside the spatial lattice.
    SeedDomain { site: usize, vertices: usize },
    /// The spatial lattice failed validation.
    Topology(String),
}

impl core::fmt::Display for DpProcessError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoLayers => write!(formatter, "layered process needs at least one layer"),
            Self::Probability(source) => write!(formatter, "invalid bond probability: {source}"),
            Self::SeedDomain { site, vertices } => write!(
                formatter,
                "seed site {site} outside the spatial lattice (0..{vertices})"
            ),
            Self::Topology(reason) => write!(formatter, "invalid spatial lattice: {reason}"),
        }
    }
}

impl std::error::Error for DpProcessError {}
