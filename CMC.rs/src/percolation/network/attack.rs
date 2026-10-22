//! Owned adapter for degree-ordered targeted attack trajectories (N3).
//!
//! The first production attack removes vertices in descending **original**
//! degree (cached in the [`NetworkProfile`]); ties resolve by an explicit
//! policy — ascending dense vertex ID (deterministic, seed-independent) or a
//! seeded uniform shuffle of each tie group (one realization per seed).
//! Adaptive residual-degree attack is a different process and is out of
//! scope here by roadmap decision.
//!
//! One `sweep()` samples one complete removal trajectory (the tie shuffle
//! is the only randomness); the following `measure()` walks the trajectory,
//! analyzing the survivor configuration at every checkpoint, and records
//! the robustness curve as ONE C0 element-preserving array sample — one
//! trajectory is one realization, never a collection of independent scalar
//! samples:
//!
//! - `RobustnessCurve[i]`: `S_max / V` (fixed denominator `V`, the F4
//!   `GiantFraction` convention) at checkpoint `i`.
//! - `RobustnessSusceptibility[i]`: the F4 finite-cluster susceptibility
//!   (canonical-largest excluded, active-denominator) at checkpoint `i`.
//!   Recorded only when defined at **every** checkpoint of the trajectory;
//!   the scalar `RobustnessSusceptibilityDefined` indicator (0/1) is always
//!   recorded, mirroring the F4 defined-indicator convention.
//! - `RobustnessRemovedFraction[i]`: removed fraction `r_i / V` at
//!   checkpoint `i` — the curve coordinates, constant across samples, so a
//!   merged result row is self-describing.
//!
//! Provenance: callers record the attack policy and grid on the run's
//! `Results` metadata, e.g.
//! `results.set_attribute("attack_policy", mc.policy().as_provenance())`
//! (C0 generic attributes; Carlo.rs never interprets them).
//!
//! Physics anchors: checkpoint values at complete degree-group boundaries
//! are validated against Callaway et al., PRL 85, 5468 (2000) Eq. (14)
//! (`S = F0(1) - F0(u)`, `u = 1 - F1(1) + F1(u)` for degree-dependent
//! occupation `q_k`), and the attack threshold against Cohen et al., PRL
//! 86, 3682 (2001); both PDFs are archived in `docs/ref/`. Evidence in
//! `CMC.rs/PERCOLATION.md` (N3).

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::super::{analyze, ComponentWorkspace, StaticConfiguration};
use super::{NetworkPercolationError, NetworkProfile};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Removal-order policy for the degree attack.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AttackPolicy {
    /// Descending original degree; equal degrees resolve by ascending dense
    /// vertex ID. Deterministic: identical order for every seed.
    DegreeDescendingStableId,
    /// Descending original degree; each maximal run of equal-degree vertices
    /// is uniformly shuffled with the run RNG. One realization per seed;
    /// identical order for identical seeds.
    DegreeDescendingSeededShuffle,
}

impl AttackPolicy {
    /// Stable provenance string for `Results` metadata attributes.
    pub const fn as_provenance(self) -> &'static str {
        match self {
            Self::DegreeDescendingStableId => "degree-descending/stable-id",
            Self::DegreeDescendingSeededShuffle => "degree-descending/seeded-shuffle",
        }
    }
}

/// Where along the removal trajectory a checkpoint is recorded.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AttackGrid {
    /// A checkpoint after every removal (plus the intact graph): `V + 1`
    /// points. Exact but `O(V * (V + E))` per trajectory.
    Full,
    /// Checkpoints at `0, step, 2*step, ...` removals, always ending at `V`
    /// (the final interval may be shorter than `step`).
    Every { step: usize },
    /// Checkpoints at exactly these removal counts — no endpoints are
    /// forced, so a grid can stay inside the region where every checkpoint
    /// has a defined susceptibility. Useful for theory comparisons at
    /// degree-group boundaries. Values are deduplicated and sorted; every
    /// value must lie within `0..=V` and the grid must be non-empty.
    At { removals: Vec<usize> },
}

impl AttackGrid {
    /// Checkpoint removal counts for a graph with `vertices` vertices.
    fn checkpoints(&self, vertices: usize) -> Result<Vec<usize>, NetworkPercolationError> {
        match self {
            Self::Full => Ok((0..=vertices).collect()),
            Self::Every { step } => {
                if *step == 0 {
                    return Err(NetworkPercolationError::InvalidGrid("step must be >= 1"));
                }
                let mut points: Vec<usize> = (0..=vertices).step_by(*step).collect();
                if points.last() != Some(&vertices) {
                    points.push(vertices);
                }
                Ok(points)
            }
            Self::At { removals } => {
                if removals.is_empty() {
                    return Err(NetworkPercolationError::InvalidGrid(
                        "at least one checkpoint is required",
                    ));
                }
                let mut points = removals.clone();
                points.sort_unstable();
                points.dedup();
                if let Some(&removal) = points.last() {
                    if removal > vertices {
                        return Err(NetworkPercolationError::GridDomain { removal, vertices });
                    }
                }
                Ok(points)
            }
        }
    }
}

/// Owned adapter producing robustness curves for degree-ordered attacks.
///
/// Configure `thermalization_sweeps = 0`; each sweep is one complete
/// trajectory realization. There is no boundary query or scalar plan here:
/// the trajectory arrays above are the adapter's entire output.
pub struct DegreeAttackMC {
    lattice: CsrLattice,
    profile: NetworkProfile,
    policy: AttackPolicy,
    checkpoints: Vec<usize>,
    /// Dense vertex indices in removal order; rebuilt every sweep.
    order: Vec<usize>,
    configuration: StaticConfiguration,
    workspace: ComponentWorkspace,
    /// Curve values refilled every measure; coordinates are fixed.
    curve: Vec<f64>,
    susceptibility: Vec<f64>,
    coordinates: Vec<f64>,
    /// Whether the last trajectory had a defined susceptibility at every
    /// checkpoint (the whole-array defined indicator).
    all_susceptibility_defined: bool,
    /// Per-checkpoint susceptibility definedness of the last trajectory.
    susceptibility_defined: Vec<bool>,
}

impl DegreeAttackMC {
    /// Build the adapter over an owned lattice.
    ///
    /// `V^2` must fit the exact `f64` integer range because the
    /// susceptibility curve squares component sizes (the F4 moment gate).
    pub fn new(
        lattice: CsrLattice,
        policy: AttackPolicy,
        grid: AttackGrid,
    ) -> Result<Self, NetworkPercolationError> {
        let vertices = lattice.vertex_count();
        lattice
            .validate()
            .map_err(|error| NetworkPercolationError::Topology(error.to_string()))?;
        if vertices == 0 {
            return Err(NetworkPercolationError::EmptyGraph);
        }
        if (vertices as u128) * (vertices as u128) > (1_u128 << f64::MANTISSA_DIGITS) {
            return Err(NetworkPercolationError::Observable(
                super::super::ObservablePlanError::InexactMomentConversion {
                    total_vertices: vertices,
                },
            ));
        }
        let profile =
            NetworkProfile::from_graph(&lattice).map_err(|_| NetworkPercolationError::Capacity)?;
        let checkpoints = grid.checkpoints(vertices)?;
        let configuration = StaticConfiguration::new(vertices, lattice.edge_count());
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(vertices, 0, false)
            .map_err(|_| NetworkPercolationError::Capacity)?;
        let coordinates = checkpoints
            .iter()
            .map(|&removal| removal as f64 / vertices as f64)
            .collect();
        Ok(Self {
            lattice,
            profile,
            policy,
            checkpoints,
            order: Vec::new(),
            configuration,
            workspace,
            curve: Vec::new(),
            susceptibility: Vec::new(),
            coordinates,
            all_susceptibility_defined: true,
            susceptibility_defined: Vec::new(),
        })
    }

    /// The owned lattice under attack.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The original-degree profile the order is built from.
    pub const fn profile(&self) -> &NetworkProfile {
        &self.profile
    }

    /// The removal-order policy.
    pub const fn policy(&self) -> AttackPolicy {
        self.policy
    }

    /// Removal counts recorded per trajectory (the grid).
    pub fn checkpoints(&self) -> &[usize] {
        &self.checkpoints
    }

    /// Giant-fraction values of the last walked trajectory.
    pub fn curve_values(&self) -> &[f64] {
        &self.curve
    }

    /// Susceptibility values of the last walked trajectory (entries at
    /// checkpoints where the value is undefined read `0.0`; consult
    /// [`Self::susceptibility_defined_at`]).
    pub fn susceptibility_values(&self) -> &[f64] {
        &self.susceptibility
    }

    /// Whether the susceptibility was defined at checkpoint `index` of the
    /// last walked trajectory.
    pub fn susceptibility_defined_at(&self, index: usize) -> bool {
        self.susceptibility_defined
            .get(index)
            .copied()
            .unwrap_or(false)
    }

    /// The removal-order of the last sweep, as dense vertex indices.
    ///
    /// Deterministic ties (`StableId`) produce the same order for every
    /// seed; seeded shuffle reproduces for a fixed seed.
    pub fn removal_order(&self) -> &[usize] {
        &self.order
    }

    /// Build the removal order for one realization (tie randomness only).
    fn sample_order<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        self.order.clear();
        self.order.extend(0..self.lattice.vertex_count());
        self.order.sort_unstable_by(|&left, &right| {
            let left_degree = self.profile.degree(left);
            let right_degree = self.profile.degree(right);
            right_degree
                .cmp(&left_degree)
                .then_with(|| left.cmp(&right))
        });
        if matches!(self.policy, AttackPolicy::DegreeDescendingSeededShuffle) {
            // Uniformly shuffle every maximal equal-degree run in place.
            let mut start = 0;
            while start < self.order.len() {
                let degree = self.profile.degree(self.order[start]);
                let mut end = start + 1;
                while end < self.order.len() && self.profile.degree(self.order[end]) == degree {
                    end += 1;
                }
                // Fisher-Yates over order[start..end]; random_range already
                // returns an absolute index in [position, end).
                for position in start..end {
                    let pick = rng.random_range(position..end);
                    self.order.swap(position, pick);
                }
                start = end;
            }
        }
    }

    /// Walk the trajectory once, filling the curve buffers.
    ///
    /// Deterministic given `removal_order`; exposed for tests and for
    /// callers running trajectories without a Carlo.rs context.
    pub fn walk_trajectory(&mut self, removal_order: &[usize]) {
        let vertices = self.lattice.vertex_count();
        debug_assert_eq!(removal_order.len(), vertices);
        self.curve.clear();
        self.susceptibility.clear();
        self.susceptibility_defined.clear();
        let mut all_defined = true;
        let mut removed_so_far = 0usize;
        // Everyone starts active; edges are irrelevant to a vertex-removal
        // process and stay active so `analyze` uses every survivor edge.
        self.configuration.vertex_mask_mut().fill(true);
        self.configuration.edge_mask_mut().fill(true);
        for &checkpoint in &self.checkpoints {
            while removed_so_far < checkpoint {
                let victim = removal_order[removed_so_far];
                self.configuration.vertex_mask_mut().set(victim, false);
                removed_so_far += 1;
            }
            let summary = analyze(
                &self.lattice,
                &self.configuration,
                &self.configuration,
                &[],
                &mut self.workspace,
            )
            .expect("constructor prepared matching analysis domains")
            .summary();
            self.curve
                .push(summary.largest_component_size as f64 / vertices as f64);
            let finite = summary
                .active_vertex_count
                .saturating_sub(summary.largest_component_size);
            let defined = finite > 0;
            all_defined &= defined;
            if defined {
                let largest_squared = (summary.largest_component_size as u128).pow(2);
                self.susceptibility
                    .push((summary.raw_second_moment - largest_squared) as f64 / finite as f64);
            } else {
                self.susceptibility.push(0.0);
            }
            self.susceptibility_defined.push(defined);
        }
        self.all_susceptibility_defined = all_defined;
    }
}

impl MonteCarlo for DegreeAttackMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sample_order(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let order = std::mem::take(&mut self.order);
        self.walk_trajectory(&order);
        self.order = order;
        ctx.measure_array("RobustnessCurve", &self.curve)
            .expect("checkpoint count is fixed by the grid");
        ctx.measure_array("RobustnessRemovedFraction", &self.coordinates)
            .expect("checkpoint count is fixed by the grid");
        ctx.measure(
            "RobustnessSusceptibilityDefined",
            f64::from(u8::from(self.all_susceptibility_defined)),
        );
        if self.all_susceptibility_defined {
            ctx.measure_array("RobustnessSusceptibility", &self.susceptibility)
                .expect("checkpoint count is fixed by the grid");
        }
    }

    fn name(&self) -> &'static str {
        "DegreeAttack"
    }
}
