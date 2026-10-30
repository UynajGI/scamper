//! L2: chemical-distance observables on long-range-percolation boxes
//! (roadmap v2 §4, stage L2), over the L1 edge law
//! ([`LongRangeBox`], module `percolation::long_range`).
//!
//! # Observables
//!
//! The chemical distance `D(x, y)` is the unweighted graph distance on the
//! occupied edge set of one L1 realization (edge weights do not exist at
//! this layer — weighted shortest paths are the v2 F1/F2 track).
//!
//! - **Distance profile** ([`ChemicalDistanceWorkspace::profile`]): one
//!   `f64` per box vertex, entry `y` = `min_{x in sources} D(x, y)` from a
//!   caller-chosen non-empty source set (one vertex for point-to-point
//!   distances, several for a multi-source spread). One realization is ONE
//!   C0 vector sample of fixed shape `vertices` — never a scalar distance
//!   alone. Unreachable targets carry the F1 sentinel `f64::INFINITY`
//!   (zero-reset value of the distance array; reachability is the derived
//!   `is_finite` test).
//! - **Box diameter**: `max over all pairs (including the diagonal) of the
//!   finite chemical distance` — the maximum eccentricity over all
//!   vertices. Whenever the box graph is connected this is exactly the
//!   graph diameter; on a disconnected realization it is the largest
//!   within-component distance (isolates contribute 0), documented rather
//!   than rejected. Always defined, always finite.
//!
//! # Occupancy: the L1 draw contract replayed in-layer
//!
//! The occupied edge set is produced by replaying the L1 contract in this
//! layer (never by re-running the L1 substrate): one uniform per pair in
//! the documented ascending pair order — saturated and near-zero pairs
//! included, nothing skipped — with the pair at index k occupied iff
//! `uniforms[k] < LongRangeBox::pair_probability(r, beta, s)` (the pinned
//! bitwise expression, here evaluated once per ℓ1 distance `r` into a
//! table — the same deterministic function, the same bits). For identical
//! uniform fields the occupied set is bitwise the one `long_range_run`
//! exposes through `occupied_pairs()`; that parity is a validation gate.
//!
//! # Regime semantics (anchor, corrected)
//!
//! Anchor (archived
//! `docs/ref/md/Biskup_2003_long-range_chemical_distance.md`): Biskup,
//! *On the scaling of the chemical distance in long-range percolation*,
//! Ann. Probab. 32, 2938 (2004). The d = 1 regime boundaries fall at
//! s = 1 and s = 2, and the statements differ by observable — this layer
//! keeps them strictly apart:
//!
//! - **s in (1, 2)** (Biskup d < s < 2d): the sharp statement is the
//!   TYPICAL distance `D(x,y) = (log|x-y|)^(Delta + o(1))` for x, y in
//!   the infinite cluster, with `Delta^-1 = log2(2/s)` (abstract and
//!   §1.4, ~lines 24-28, 246-254). Bracketed through the profile.
//! - **s = 1** (Biskup s = d): Eq. (1.9) gives the DIAMETER
//!   `D_N = Θ(log N / log log N)`; the typical distance is "at most
//!   logarithmic" and the lower bound is explicitly missing (~lines
//!   235-245). Only the upper envelope is bracketed, never a lower claim.
//! - **s < 1** (Biskup s < d): Eq. (1.8) is the stochastic-dimension
//!   SUP statement `sup_{x,y in Z^d} D(x,y) = d/(d-s)` almost surely — a
//!   constant, not a power of N — with "a similar asymptotic statement"
//!   for the N → ∞ box diameter `D_N` (~lines 222-234, Theorem 4.1 of
//!   Benjamini-Berger). The DIAMETER is bracketed against that bounded
//!   constant scale; finite boxes may sit above the infinite-lattice sup
//!   because restricting to the box removes shortcut vertices outside it.
//!   The typical distance is smaller and unproved there, so no
//!   typical-distance gate runs at s <= 1.
//! - **s > 2** (Biskup s > 2d): Biskup expects the box DIAMETER `D_N` to
//!   scale linearly with N, proved in d = 1 (Benjamini-Berger; ~lines
//!   269-273). The gate here brackets the TYPICAL distance — the
//!   observable this layer profiles — as linear with a positive time
//!   constant, consistent with but not literally the quoted `D_N`
//!   statement.
//! - **s = 2** (Biskup s = 2d): "not very well understood" (~line 255) —
//!   deliberately NOT bracketed; no gate at s = 2.
//!
//! Biskup's theorems describe the infinite supercritical cluster on Z^d;
//! finite-box measurements can only be statistical brackets of these
//! asymptotics, and the validated brackets (with their finite-size
//! caveats) are recorded in PERCOLATION.md.
//!
//! # Connectivity parameterization
//!
//! At `beta >= 1` every r = 1 pair is surely occupied (L1 saturation), so
//! in d = 1 the nearest-neighbor chain spans the box: the box graph is
//! connected with probability 1 and the empirical disconnection rate is 0
//! by construction. All typical-distance gates run at `beta >= 1`
//! (disconnected pairs would pollute means; no conditioning needed). At
//! `beta < 1` disconnections occur and profile entries carry INFINITY —
//! the estimator below conditions per entry on reachability and reports
//! the reachable fraction, which is also the correct conditioning for
//! Biskup's `x, y in C_infinity` statement.
//!
//! # Cost model and the diameter cap
//!
//! The draw contract forces one full pair scan per realization:
//! `O(pairs)` occupancy decisions (`O(vertices^2)` — quadratic in the box
//! regardless of the realized edge count), then `O(vertices + m)` per
//! BFS. The exact diameter needs BFS from EVERY vertex, `O(vertices *
//! (vertices + m))` — `O(vertices^3)` in the complete-graph worst case, so
//! [`DIAMETER_VERTEX_CAP`] is a hard guard-rail rather than a feasibility
//! guarantee (the sparse regimes tested here are fast); the single-source
//! profile path has no cap beyond memory. The workspace reserves the
//! complete-graph worst case as CAPACITY only (`2 * pairs` adjacency
//! slots — untouched pages cost nothing) and leaves the neighbor buffer
//! empty, so after `prepare` a realization performs zero allocations and
//! the pre-run accessors report an empty adjacency.

use core::fmt;
use std::collections::TryReserveError;

use rand::{Rng, RngExt};

use super::estimator::Moments;
use super::long_range::LongRangeBox;
use super::process::ProcessError;

/// Largest box the all-vertex diameter pass accepts (`2^12` vertices).
///
/// A hard guard-rail: the pass is `O(vertices * (vertices + m))` and
/// `O(vertices^3)` in the complete-graph worst case, so this is not a
/// feasibility guarantee. The single-source profile path is not capped.
pub const DIAMETER_VERTEX_CAP: usize = 4096;

/// Whether a run also computes the box diameter (exact, BFS from every
/// vertex — see [`DIAMETER_VERTEX_CAP`]).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum DiameterTracking {
    /// Profile only: one multi-source BFS.
    #[default]
    Off,
    /// Profile plus the exact diameter: `vertices` BFS passes.
    AllVertexBfs,
}

/// Reusable chemical-distance workspace: the occupied-edge adjacency in
/// CSR form, the per-distance probability table, and the BFS buffers.
///
/// `prepare` reserves every buffer for the complete-graph worst case
/// (`2 * pairs` adjacency slots as capacity only — untouched pages cost
/// nothing) and is idempotent; a correctly sized workspace never
/// allocates in steady state. The neighbor buffer stays empty until a
/// run fills it, so the adjacency accessors report an empty graph before
/// the first run. The profile of the last run lives here (invalidated by
/// the next).
#[derive(Debug, Default)]
pub struct ChemicalDistanceWorkspace {
    /// CSR offsets, `vertices + 1` entries.
    offsets: Vec<usize>,
    /// CSR neighbors, `2 * pairs` slots at the worst case.
    neighbors: Vec<usize>,
    /// Fill cursors of the second occupancy pass.
    cursor: Vec<usize>,
    /// `pair_probability(r)` per ℓ1 distance (index 0 unused, set 0.0).
    probability: Vec<f64>,
    /// BFS queue.
    queue: Vec<usize>,
    /// Distance profile of the last run (one `f64` per vertex).
    profile: Vec<f64>,
    /// Integer distances of the all-vertex diameter passes.
    eccentric_scratch: Vec<usize>,
}

impl ChemicalDistanceWorkspace {
    /// An unprepared workspace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve every buffer for a box of `vertices` vertices, `pairs`
    /// unordered pairs, and ℓ1 distances `0..=distance_range`, resetting
    /// all run state. The neighbor buffer is reserved to `2 * pairs`
    /// slots but left EMPTY (a run resizes it to the realized edge
    /// count), so the pre-run adjacency is empty. Re-preparing with
    /// different sizes is supported and leaves no stale state.
    ///
    /// Validation of the derived lengths is done before any mutation;
    /// a capacity failure during reservation can leave earlier buffers
    /// re-prepared (no result is ever returned from a failed prepare).
    pub fn prepare(
        &mut self,
        vertices: usize,
        pairs: usize,
        distance_range: usize,
    ) -> Result<(), ChemicalDistanceError> {
        let offset_len = vertices
            .checked_add(1)
            .ok_or(ChemicalDistanceError::DimensionOverflow { field: "offsets" })?;
        let probability_len =
            distance_range
                .checked_add(1)
                .ok_or(ChemicalDistanceError::DimensionOverflow {
                    field: "probabilities",
                })?;
        let neighbor_slots = 2_usize
            .checked_mul(pairs)
            .ok_or(ChemicalDistanceError::DimensionOverflow { field: "neighbors" })?;
        reserve_exact(&mut self.offsets, offset_len, 0, "offsets")?;
        reserve_capacity(&mut self.neighbors, neighbor_slots, "neighbors")?;
        reserve_exact(&mut self.cursor, vertices, 0, "cursors")?;
        reserve_exact(&mut self.probability, probability_len, 0.0, "probabilities")?;
        reserve_exact(&mut self.queue, vertices, 0, "queue")?;
        reserve_exact(&mut self.profile, vertices, f64::INFINITY, "profile")?;
        reserve_exact(
            &mut self.eccentric_scratch,
            vertices,
            usize::MAX,
            "eccentricities",
        )?;
        Ok(())
    }

    /// Drop the last run's adjacency and profile, retaining all capacity.
    pub fn reset(&mut self) {
        self.offsets.clear();
        self.neighbors.clear();
        self.cursor.clear();
        self.probability.clear();
        self.queue.clear();
        self.profile.clear();
        self.eccentric_scratch.clear();
    }

    /// The distance profile of the last run: one `f64` per box vertex,
    /// `min_{x in sources} D(x, y)` at entry `y`, `f64::INFINITY` where
    /// unreachable. All INFINITY before the first run; empty after
    /// [`Self::reset`].
    #[inline]
    pub fn profile(&self) -> &[f64] {
        &self.profile
    }

    /// CSR offsets of the last run's occupied adjacency (`vertices + 1`
    /// entries; all zeros before the first run, empty after
    /// [`Self::reset`]).
    #[inline]
    pub fn adjacency_offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// CSR neighbor slots of the last run's occupied adjacency (degree
    /// blocks given by [`Self::adjacency_offsets`]; empty before the
    /// first run).
    #[inline]
    pub fn adjacency_neighbors(&self) -> &[usize] {
        &self.neighbors
    }

    /// Occupied edge count of the last run (zero before the first run).
    #[inline]
    pub fn occupied_edge_count(&self) -> usize {
        self.neighbors.len() / 2
    }
}

/// Reserve `buffer` to exactly hold `size` entries of `fill`, resetting
/// it. Capacity only ever grows: re-preparing a larger size reallocates
/// once, while a same-size re-prepare reuses the capacity (contents are
/// still reset).
fn reserve_exact<T>(
    buffer: &mut Vec<T>,
    size: usize,
    fill: T,
    name: &'static str,
) -> Result<(), ChemicalDistanceError>
where
    T: Clone,
{
    if buffer.capacity() < size {
        buffer.clear();
        buffer
            .try_reserve_exact(size)
            .map_err(|source| ChemicalDistanceError::capacity(name, source))?;
    }
    buffer.clear();
    buffer.resize(size, fill);
    Ok(())
}

/// Ensure `buffer` can hold `size` entries and leave it EMPTY (capacity
/// only; no reallocation when capacity already suffices).
fn reserve_capacity<T>(
    buffer: &mut Vec<T>,
    size: usize,
    name: &'static str,
) -> Result<(), ChemicalDistanceError> {
    if buffer.capacity() < size {
        buffer.clear();
        buffer
            .try_reserve_exact(size)
            .map_err(|source| ChemicalDistanceError::capacity(name, source))?;
    }
    buffer.clear();
    Ok(())
}

/// Scalar outcome of one chemical-distance realization.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ChemicalDistanceOutcome {
    /// Occupied pairs of the replayed L1 realization.
    pub occupied_edges: usize,
    /// Vertices reachable from the source set (finite profile entries).
    pub source_component_size: usize,
    /// Profile entries carrying the INFINITY sentinel.
    pub unreachable_count: usize,
    /// The exact box diameter (maximum finite distance over all pairs);
    /// `None` unless [`DiameterTracking::AllVertexBfs`] was requested.
    pub diameter: Option<usize>,
}

/// Run one chemical-distance realization on an L1 box (the reusable,
/// RNG-free core).
///
/// `uniforms` holds one `f64` per pair in the L1 consideration order
/// (ascending pair index, every entry finite and in `[0, 1]`); the
/// occupied set is the L1 draw-contract replay described in the module
/// documentation, materialized in `workspace` as a CSR adjacency. The
/// multi-source BFS over `sources` (non-empty, in-range; duplicates are
/// harmless — the seeded queue ignores a repeat) writes the profile;
/// `tracking` adds the exact all-vertex diameter. The workspace is
/// re-prepared idempotently; every validation rejection above leaves it
/// untouched, while a capacity failure during preparation can leave it
/// partially re-prepared (no result is returned in that case).
pub fn chemical_distance_run(
    box_: &LongRangeBox,
    beta: f64,
    s: f64,
    uniforms: &[f64],
    sources: &[usize],
    tracking: DiameterTracking,
    workspace: &mut ChemicalDistanceWorkspace,
) -> Result<ChemicalDistanceOutcome, ChemicalDistanceError> {
    validate_parameters(beta, s)?;
    let vertices = box_.vertices();
    if sources.is_empty() {
        return Err(ChemicalDistanceError::NoSources);
    }
    for (index, &source) in sources.iter().enumerate() {
        if source >= vertices {
            return Err(ChemicalDistanceError::SourceOutOfRange {
                index,
                source,
                vertices,
            });
        }
    }
    if tracking == DiameterTracking::AllVertexBfs && vertices > DIAMETER_VERTEX_CAP {
        return Err(ChemicalDistanceError::DiameterCapExceeded {
            vertices,
            cap: DIAMETER_VERTEX_CAP,
        });
    }
    if uniforms.len() != box_.pairs() {
        return Err(ChemicalDistanceError::UniformLengthMismatch {
            expected: box_.pairs(),
            actual: uniforms.len(),
        });
    }
    for (index, &value) in uniforms.iter().enumerate() {
        if !value.is_finite() || value < 0.0 || value > 1.0 {
            return Err(ChemicalDistanceError::UniformInvalid { index, value });
        }
    }
    let distance_range: usize = box_.sides().iter().map(|&side| side - 1).sum();
    workspace.prepare(vertices, box_.pairs(), distance_range)?;

    // Probability per ℓ1 distance — the pinned L1 expression evaluated
    // once per distance (bit-identical to the per-pair evaluation).
    for distance in 1..=distance_range {
        workspace.probability[distance] = LongRangeBox::pair_probability(distance, beta, s);
    }

    // Occupancy pass 1: degrees. The pair scan is the L1 consideration
    // order (lexicographic (x, y), x < y — ascending pair index).
    let mut pair = 0_usize;
    let mut occupied_edges = 0_usize;
    for first in 0..vertices {
        for second in (first + 1)..vertices {
            if uniforms[pair] < workspace.probability[box_.distance(first, second)] {
                workspace.offsets[first + 1] += 1;
                workspace.offsets[second + 1] += 1;
                occupied_edges += 1;
            }
            pair += 1;
        }
    }
    debug_assert_eq!(pair, box_.pairs());
    for index in 0..vertices {
        workspace.offsets[index + 1] += workspace.offsets[index];
    }
    workspace.neighbors.resize(2 * occupied_edges, usize::MAX);

    // Occupancy pass 2: fill the CSR blocks.
    workspace
        .cursor
        .copy_from_slice(&workspace.offsets[..vertices]);
    let mut pair = 0_usize;
    for first in 0..vertices {
        for second in (first + 1)..vertices {
            if uniforms[pair] < workspace.probability[box_.distance(first, second)] {
                workspace.neighbors[workspace.cursor[first]] = second;
                workspace.cursor[first] += 1;
                workspace.neighbors[workspace.cursor[second]] = first;
                workspace.cursor[second] += 1;
            }
            pair += 1;
        }
    }
    debug_assert_eq!(pair, box_.pairs());

    // Multi-source BFS: the profile.
    workspace.profile.fill(f64::INFINITY);
    workspace.queue.clear();
    for &source in sources {
        if workspace.profile[source] != 0.0 {
            workspace.profile[source] = 0.0;
            workspace.queue.push(source);
        }
    }
    let mut head = 0_usize;
    while head < workspace.queue.len() {
        let vertex = workspace.queue[head];
        head += 1;
        let distance = workspace.profile[vertex] + 1.0;
        for &neighbor in
            &workspace.neighbors[workspace.offsets[vertex]..workspace.offsets[vertex + 1]]
        {
            if workspace.profile[neighbor] == f64::INFINITY {
                workspace.profile[neighbor] = distance;
                workspace.queue.push(neighbor);
            }
        }
    }
    let unreachable_count = workspace
        .profile
        .iter()
        .filter(|entry| !entry.is_finite())
        .count();

    // The exact diameter: eccentricity of every vertex, maximum finite
    // distance over all pairs (the diagonal contributes 0, so the value
    // is defined even on an edgeless box).
    let diameter = match tracking {
        DiameterTracking::Off => None,
        DiameterTracking::AllVertexBfs => {
            let mut diameter = 0_usize;
            for source in 0..vertices {
                workspace.eccentric_scratch.fill(usize::MAX);
                workspace.eccentric_scratch[source] = 0;
                workspace.queue.clear();
                workspace.queue.push(source);
                let mut head = 0_usize;
                let mut eccentricity = 0_usize;
                while head < workspace.queue.len() {
                    let vertex = workspace.queue[head];
                    head += 1;
                    let distance = workspace.eccentric_scratch[vertex] + 1;
                    for &neighbor in &workspace.neighbors
                        [workspace.offsets[vertex]..workspace.offsets[vertex + 1]]
                    {
                        if workspace.eccentric_scratch[neighbor] == usize::MAX {
                            workspace.eccentric_scratch[neighbor] = distance;
                            eccentricity = eccentricity.max(distance);
                            workspace.queue.push(neighbor);
                        }
                    }
                }
                diameter = diameter.max(eccentricity);
            }
            Some(diameter)
        }
    };

    Ok(ChemicalDistanceOutcome {
        occupied_edges,
        source_component_size: vertices - unreachable_count,
        unreachable_count,
        diameter,
    })
}

/// Shared (beta, s) domain validation, mirroring the L1 core: both must
/// be finite and non-negative (zero is legal but degenerate there).
fn validate_parameters(beta: f64, s: f64) -> Result<(), ChemicalDistanceError> {
    if !beta.is_finite() || beta < 0.0 {
        return Err(ChemicalDistanceError::BetaInvalid { value: beta });
    }
    if !s.is_finite() || s < 0.0 {
        return Err(ChemicalDistanceError::ExponentInvalid { value: s });
    }
    Ok(())
}

/// Reusable realization sampler: the uniform field, the workspace, and
/// the profile of the last realization.
///
/// `new` reserves every buffer once (the complete-graph worst case, see
/// the module documentation); [`Self::realize`] then performs zero
/// steady-state allocations. The draw is one `f64` per pair in ascending
/// pair order — bit-identical to the stream of `LongRangePercolationMC`
/// at the same seed, so both layers realize the same occupied graph.
/// The returned profile of each call is invalidated by the next.
#[derive(Debug)]
pub struct ChemicalDistanceSampler {
    box_: LongRangeBox,
    beta: f64,
    s: f64,
    sources: Vec<usize>,
    uniforms: Vec<f64>,
    workspace: ChemicalDistanceWorkspace,
    outcome: ChemicalDistanceOutcome,
}

impl ChemicalDistanceSampler {
    /// Build the sampler. Every domain is validated here so `realize`
    /// cannot fail: `(beta, s)` as in L1, `sources` non-empty and
    /// in-range (duplicates allowed, harmless).
    pub fn new(
        box_: LongRangeBox,
        beta: f64,
        s: f64,
        sources: &[usize],
    ) -> Result<Self, ChemicalDistanceError> {
        validate_parameters(beta, s)?;
        if sources.is_empty() {
            return Err(ChemicalDistanceError::NoSources);
        }
        for (index, &source) in sources.iter().enumerate() {
            if source >= box_.vertices() {
                return Err(ChemicalDistanceError::SourceOutOfRange {
                    index,
                    source,
                    vertices: box_.vertices(),
                });
            }
        }
        let distance_range: usize = box_.sides().iter().map(|&side| side - 1).sum();
        let mut workspace = ChemicalDistanceWorkspace::new();
        workspace.prepare(box_.vertices(), box_.pairs(), distance_range)?;
        let mut uniforms = Vec::new();
        reserve_exact(&mut uniforms, box_.pairs(), 0.0, "uniforms")?;
        Ok(Self {
            box_,
            beta,
            s,
            sources: sources.to_vec(),
            uniforms,
            workspace,
            outcome: ChemicalDistanceOutcome::default(),
        })
    }

    /// Draw the uniform field (the L1 pair-order contract) and run one
    /// profile realization from the configured sources; the returned
    /// slice is the profile (one `f64` per vertex).
    pub fn realize<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
    ) -> Result<&[f64], ChemicalDistanceError> {
        self.outcome = self.run(rng, DiameterTracking::Off)?;
        Ok(self.workspace.profile())
    }

    /// [`Self::realize`] plus the exact all-vertex diameter (subject to
    /// [`DIAMETER_VERTEX_CAP`]); the profile stays available through
    /// [`Self::profile`].
    pub fn realize_with_diameter<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
    ) -> Result<ChemicalDistanceOutcome, ChemicalDistanceError> {
        self.outcome = self.run(rng, DiameterTracking::AllVertexBfs)?;
        Ok(self.outcome)
    }

    /// The shared draw + core call of both realize paths.
    fn run<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
        tracking: DiameterTracking,
    ) -> Result<ChemicalDistanceOutcome, ChemicalDistanceError> {
        for slot in self.uniforms.iter_mut() {
            *slot = rng.random::<f64>();
        }
        let Self {
            box_,
            beta,
            s,
            sources,
            uniforms,
            workspace,
            ..
        } = self;
        chemical_distance_run(box_, *beta, *s, uniforms, sources, tracking, workspace)
    }

    /// The box the sampler runs on.
    pub const fn geometry(&self) -> &LongRangeBox {
        &self.box_
    }

    /// The intensity `beta`.
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// The decay exponent `s`.
    pub const fn exponent(&self) -> f64 {
        self.s
    }

    /// The configured source set.
    #[inline]
    pub fn sources(&self) -> &[usize] {
        &self.sources
    }

    /// Outcome of the last realization.
    pub const fn last_outcome(&self) -> ChemicalDistanceOutcome {
        self.outcome
    }

    /// The distance profile of the last realization (all
    /// [`f64::INFINITY`] before the first).
    #[inline]
    pub fn profile(&self) -> &[f64] {
        self.workspace.profile()
    }

    /// The uniform field of the last realization (one `f64` per pair in
    /// consideration order) — the replay handle for the draw contract.
    #[inline]
    pub fn uniforms(&self) -> &[f64] {
        &self.uniforms
    }

    /// The workspace of the last realization (adjacency and profile).
    pub const fn workspace(&self) -> &ChemicalDistanceWorkspace {
        &self.workspace
    }
}

/// Accumulates independent chemical-distance realizations into sample
/// statistics.
///
/// Profile entries are conditioned per entry on reachability: INFINITY
/// entries are not accumulated, the per-entry reachable fraction is
/// reported (at `beta >= 1` in d = 1 it is identically 1), and the mean
/// of an entry is taken over the realizations where it is finite — the
/// conditioning of Biskup's `x, y in C_infinity` statement. Diameters
/// are always finite (their definition never produces INFINITY) and are
/// accumulated unconditionally. All buffers are reserved at `new`, so
/// neither `add` allocates.
#[derive(Clone, Debug)]
pub struct ChemicalDistanceEstimator {
    vertices: usize,
    profile_samples: usize,
    profile: Vec<Moments>,
    diameter: Moments,
}

impl ChemicalDistanceEstimator {
    /// Reserve the accumulators for profiles of `vertices` entries
    /// (`>= 1`).
    pub fn new(vertices: usize) -> Result<Self, ChemicalDistanceError> {
        if vertices == 0 {
            return Err(ChemicalDistanceError::NoVertices);
        }
        let mut profile = Vec::new();
        reserve_exact(&mut profile, vertices, Moments::new(), "profile moments")?;
        Ok(Self {
            vertices,
            profile_samples: 0,
            profile,
            diameter: Moments::new(),
        })
    }

    /// The reserved profile length.
    #[inline]
    pub fn vertices(&self) -> usize {
        self.vertices
    }

    /// Profiles accumulated so far.
    #[inline]
    pub const fn profile_samples(&self) -> usize {
        self.profile_samples
    }

    /// Diameters accumulated so far.
    #[inline]
    pub fn diameter_samples(&self) -> usize {
        self.diameter.samples()
    }

    /// Add one realization's profile. Length mismatches are typed
    /// rejections; nothing is added on rejection.
    pub fn add_profile(&mut self, profile: &[f64]) -> Result<(), ChemicalDistanceError> {
        if profile.len() != self.vertices {
            return Err(ChemicalDistanceError::ProfileLength {
                expected: self.vertices,
                actual: profile.len(),
            });
        }
        self.profile_samples += 1;
        for (moments, &value) in self.profile.iter_mut().zip(profile) {
            if value.is_finite() {
                moments.add(value);
            }
        }
        Ok(())
    }

    /// Add one realization's diameter (always finite by definition).
    /// The square is taken in `f64`, so no `usize` overflow is possible
    /// for any caller-supplied value.
    pub fn add_diameter(&mut self, diameter: usize) {
        self.diameter.add(diameter as f64);
    }

    /// Sample statistics, or `None` before any profile and any pair of
    /// diameters. Entries never reached in any realization have NaN mean
    /// and 0.0 error (the error is only defined from two reachable
    /// samples); entries reached once have zero error. Allocates the
    /// estimate's vectors (report-time, never per realization).
    pub fn estimate(&self) -> Option<ChemicalDistanceEstimate> {
        if self.profile_samples == 0 && self.diameter.samples() < 2 {
            return None;
        }
        let profile_mean = self
            .profile
            .iter()
            .map(|moments| moments.mean().unwrap_or(f64::NAN))
            .collect();
        let profile_standard_error = self
            .profile
            .iter()
            .map(|moments| {
                if moments.samples() < 2 {
                    0.0
                } else {
                    moments
                        .standard_error()
                        .expect("at least two reachable samples")
                }
            })
            .collect();
        let reachable_fraction = self
            .profile
            .iter()
            .map(|moments| {
                if self.profile_samples > 0 {
                    moments.samples() as f64 / self.profile_samples as f64
                } else {
                    f64::NAN
                }
            })
            .collect();
        let (diameter_mean, diameter_standard_error) = if self.diameter.samples() >= 2 {
            (
                self.diameter.mean().expect("at least two diameters"),
                self.diameter
                    .standard_error()
                    .expect("at least two diameters"),
            )
        } else {
            (f64::NAN, f64::NAN)
        };
        Some(ChemicalDistanceEstimate {
            profile_samples: self.profile_samples,
            profile_mean,
            profile_standard_error,
            reachable_fraction,
            diameter_samples: self.diameter.samples(),
            diameter_mean,
            diameter_standard_error,
        })
    }
}

/// Sample statistics over independent chemical-distance realizations.
#[derive(Clone, Debug, PartialEq)]
pub struct ChemicalDistanceEstimate {
    /// Profiles accumulated.
    pub profile_samples: usize,
    /// Per-vertex mean over the realizations where the entry was finite
    /// (NaN where never reached — see the estimator documentation).
    pub profile_mean: Vec<f64>,
    /// Standard error of [`Self::profile_mean`] (0.0 below two reachable
    /// samples).
    pub profile_standard_error: Vec<f64>,
    /// Fraction of realizations where each entry was finite (identically
    /// 1 at `beta >= 1` in d = 1).
    pub reachable_fraction: Vec<f64>,
    /// Diameters accumulated.
    pub diameter_samples: usize,
    /// Sample mean of the diameter (NaN with fewer than two samples).
    pub diameter_mean: f64,
    /// Standard error of [`Self::diameter_mean`] (NaN with fewer than
    /// two samples).
    pub diameter_standard_error: f64,
}

/// Typed rejections of the chemical-distance layer.
#[derive(Debug)]
pub enum ChemicalDistanceError {
    /// A negative or non-finite intensity.
    BetaInvalid { value: f64 },
    /// A negative or non-finite decay exponent.
    ExponentInvalid { value: f64 },
    /// An empty source set.
    NoSources,
    /// A zero-length profile (the estimator needs at least one entry).
    NoVertices,
    /// A source outside `0..vertices`.
    SourceOutOfRange {
        index: usize,
        source: usize,
        vertices: usize,
    },
    /// A uniform field with the wrong length.
    UniformLengthMismatch { expected: usize, actual: usize },
    /// A uniform outside `[0, 1]` or not finite.
    UniformInvalid { index: usize, value: f64 },
    /// An all-vertex diameter request on a box above the cost cap.
    DiameterCapExceeded { vertices: usize, cap: usize },
    /// A profile whose length differs from the estimator's vertex count.
    ProfileLength { expected: usize, actual: usize },
    /// A derived workspace length overflows `usize`.
    DimensionOverflow { field: &'static str },
    /// A buffer reservation failure.
    Capacity {
        buffer: &'static str,
        source: ProcessError,
    },
}

impl ChemicalDistanceError {
    /// Wrap a reservation failure with the buffer it happened on.
    fn capacity(buffer: &'static str, source: TryReserveError) -> Self {
        Self::Capacity {
            buffer,
            source: ProcessError::from(source),
        }
    }
}

impl fmt::Display for ChemicalDistanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BetaInvalid { value } => write!(
                formatter,
                "beta {value} is negative or not finite (zero is legal but degenerate)"
            ),
            Self::ExponentInvalid { value } => write!(
                formatter,
                "exponent s {value} is negative or not finite (zero is legal but degenerate)"
            ),
            Self::NoSources => write!(formatter, "at least one source vertex is required"),
            Self::NoVertices => write!(formatter, "the estimator needs at least one vertex"),
            Self::SourceOutOfRange {
                index,
                source,
                vertices,
            } => write!(
                formatter,
                "source {source} at position {index} is outside 0..{vertices}"
            ),
            Self::UniformLengthMismatch { expected, actual } => write!(
                formatter,
                "uniform field has {actual} entries, expected one per pair ({expected})"
            ),
            Self::UniformInvalid { index, value } => write!(
                formatter,
                "uniform {value} at pair {index} is outside [0, 1] or not finite"
            ),
            Self::DiameterCapExceeded { vertices, cap } => write!(
                formatter,
                "all-vertex diameter needs {vertices} BFS passes, above the cap {cap} \
                 (single-source profiles are not capped)"
            ),
            Self::ProfileLength { expected, actual } => write!(
                formatter,
                "profile has {actual} entries, expected one per vertex ({expected})"
            ),
            Self::DimensionOverflow { field } => write!(
                formatter,
                "the {field} buffer length overflows usize for the requested dimensions"
            ),
            Self::Capacity { buffer, source } => {
                write!(
                    formatter,
                    "could not prepare the {buffer} buffer of the chemical-distance workspace: {source}"
                )
            }
        }
    }
}

impl std::error::Error for ChemicalDistanceError {}
