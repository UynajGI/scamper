//! Generic arbitrary-spin directed-loop scattering.
//!
//! The table is built from a sparse positive local operator catalog. An
//! extended local state is `(vertex kind, entrance leg, delta)` where
//! `delta=±1` is the worm discontinuity in the local basis. Compatible states
//! are connected when changing the entrance and one exit leg maps one allowed
//! matrix element into another. Symmetric path weights enforce local detailed
//! balance exactly.

use std::collections::{HashMap, HashSet};

use rand::Rng;
use rand::RngExt;

use crate::local_space::BasisState;

use super::error::LatticeQmcError;
use super::vertex::VertexKind;

/// Scattering-table construction strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatteringPolicy {
    /// Greedy symmetric residual flow with low bounce.
    LowBounce,
    /// Symmetric fixed proposal followed by Metropolis acceptance.
    Metropolis,
}

/// One local directed-loop outcome.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatteringChoice {
    /// Matrix-element kind after the local update.
    pub new_kind: usize,
    /// Exit leg.
    pub exit_leg: usize,
    /// Discontinuity carried to the linked next leg.
    pub next_delta: i8,
    /// Row-normalized probability.
    pub probability: f64,
}

/// Construction diagnostics.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ScatteringDiagnostics {
    /// Maximum row-normalization error.
    pub max_row_error: f64,
    /// Maximum local detailed-balance residual.
    pub max_detailed_balance_error: f64,
    /// Mean bounce probability over valid extended states.
    pub mean_bounce_probability: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ExtendedState {
    kind: usize,
    entrance: usize,
    delta: i8,
}

/// Precomputed local scattering table for one operator term.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatteringTable {
    rows: Vec<Vec<ScatteringChoice>>,
    lookup: Vec<Option<usize>>,
    states: Vec<ExtendedState>,
    leg_count: usize,
    diagnostics: ScatteringDiagnostics,
    policy: ScatteringPolicy,
}

impl ScatteringTable {
    /// Build from positive local matrix elements and per-leg dimensions.
    pub fn build(
        kinds: &[VertexKind],
        leg_dimensions: &[usize],
        policy: ScatteringPolicy,
    ) -> Result<Self, LatticeQmcError> {
        if kinds.is_empty() {
            return Err(LatticeQmcError::InvalidModel(
                "an operator term needs at least one matrix element".into(),
            ));
        }
        let leg_count = kinds[0].legs().len();
        if leg_count != leg_dimensions.len() || leg_count == 0 {
            return Err(LatticeQmcError::InvalidModel(
                "leg dimensions do not match the operator arity".into(),
            ));
        }
        let mut kind_lookup = HashMap::<Vec<BasisState>, usize>::new();
        for (kind_id, kind) in kinds.iter().enumerate() {
            if kind.legs().len() != leg_count {
                return Err(LatticeQmcError::InvalidModel(
                    "all matrix elements of a term must have the same arity".into(),
                ));
            }
            for (leg, &state) in kind.legs().iter().enumerate() {
                if usize::from(state) >= leg_dimensions[leg] {
                    return Err(LatticeQmcError::InvalidModel(format!(
                        "state {state} exceeds dimension {} on leg {leg}",
                        leg_dimensions[leg]
                    )));
                }
            }
            if kind_lookup.insert(kind.legs().to_vec(), kind_id).is_some() {
                return Err(LatticeQmcError::InvalidModel(
                    "duplicate local leg pattern in one operator term".into(),
                ));
            }
        }

        let raw_count = kinds.len() * leg_count * 2;
        let mut lookup = vec![None; raw_count];
        let mut states = Vec::new();
        for (kind, vertex_kind) in kinds.iter().enumerate() {
            for entrance in 0..leg_count {
                for delta in [-1_i8, 1_i8] {
                    if shifted(vertex_kind.state(entrance), delta, leg_dimensions[entrance])
                        .is_some()
                    {
                        let compact = states.len();
                        lookup[raw_index(kind, entrance, delta, leg_count)] = Some(compact);
                        states.push(ExtendedState {
                            kind,
                            entrance,
                            delta,
                        });
                    }
                }
            }
        }

        let edges = compatible_edges(
            kinds,
            leg_dimensions,
            &kind_lookup,
            &lookup,
            &states,
            leg_count,
        );
        let rows = match policy {
            ScatteringPolicy::LowBounce => low_bounce_rows(kinds, &states, &edges),
            ScatteringPolicy::Metropolis => metropolis_rows(
                kinds,
                leg_dimensions,
                &kind_lookup,
                &lookup,
                &states,
                leg_count,
            ),
        };
        let diagnostics = diagnostics(kinds, &states, &rows);
        Ok(Self {
            rows,
            lookup,
            states,
            leg_count,
            diagnostics,
            policy,
        })
    }

    /// Whether an entrance discontinuity is represented by the table.
    pub fn has_row(&self, kind: usize, entrance: usize, delta: i8) -> bool {
        self.lookup
            .get(raw_index(kind, entrance, delta, self.leg_count))
            .and_then(|row| *row)
            .is_some()
    }

    /// Draw one local scattering outcome.
    pub fn sample<R: Rng + ?Sized>(
        &self,
        kind: usize,
        entrance: usize,
        delta: i8,
        rng: &mut R,
    ) -> Option<ScatteringChoice> {
        let row_id = self.lookup[raw_index(kind, entrance, delta, self.leg_count)]?;
        let row = &self.rows[row_id];
        let u = rng.random::<f64>();
        let mut cumulative = 0.0;
        for choice in row {
            cumulative += choice.probability;
            if u < cumulative {
                return Some(*choice);
            }
        }
        row.last().copied()
    }

    /// Number of legs in this term.
    pub fn leg_count(&self) -> usize {
        self.leg_count
    }

    /// Table diagnostics.
    pub fn diagnostics(&self) -> ScatteringDiagnostics {
        self.diagnostics
    }

    /// Construction policy.
    pub fn policy(&self) -> ScatteringPolicy {
        self.policy
    }
}

fn raw_index(kind: usize, entrance: usize, delta: i8, leg_count: usize) -> usize {
    let delta_slot = usize::from(delta > 0);
    2 * (kind * leg_count + entrance) + delta_slot
}

fn shifted(state: BasisState, delta: i8, dimension: usize) -> Option<BasisState> {
    match delta {
        1 if usize::from(state) + 1 < dimension => Some(state + 1),
        -1 if state > 0 => Some(state - 1),
        _ => None,
    }
}

fn transformed_kind(
    kinds: &[VertexKind],
    leg_dimensions: &[usize],
    kind_lookup: &HashMap<Vec<BasisState>, usize>,
    source: ExtendedState,
    exit: usize,
    exit_delta: i8,
) -> Option<usize> {
    let mut legs = kinds[source.kind].legs().to_vec();
    legs[source.entrance] = shifted(
        legs[source.entrance],
        source.delta,
        leg_dimensions[source.entrance],
    )?;
    legs[exit] = shifted(legs[exit], exit_delta, leg_dimensions[exit])?;
    kind_lookup.get(&legs).copied()
}

fn compatible_edges(
    kinds: &[VertexKind],
    leg_dimensions: &[usize],
    kind_lookup: &HashMap<Vec<BasisState>, usize>,
    lookup: &[Option<usize>],
    states: &[ExtendedState],
    leg_count: usize,
) -> Vec<(usize, usize)> {
    let mut edges = HashSet::new();
    for (source_id, &source) in states.iter().enumerate() {
        for exit in 0..leg_count {
            for exit_delta in [-1_i8, 1_i8] {
                let Some(new_kind) =
                    transformed_kind(kinds, leg_dimensions, kind_lookup, source, exit, exit_delta)
                else {
                    continue;
                };
                let reverse_delta = -exit_delta;
                let Some(target_id) = lookup[raw_index(new_kind, exit, reverse_delta, leg_count)]
                else {
                    continue;
                };
                if source_id != target_id {
                    edges.insert(if source_id < target_id {
                        (source_id, target_id)
                    } else {
                        (target_id, source_id)
                    });
                }
            }
        }
    }
    let mut edges: Vec<_> = edges.into_iter().collect();
    edges.sort_unstable();
    edges
}

fn low_bounce_rows(
    kinds: &[VertexKind],
    states: &[ExtendedState],
    edges: &[(usize, usize)],
) -> Vec<Vec<ScatteringChoice>> {
    let weights: Vec<f64> = states
        .iter()
        .map(|state| kinds[state.kind].weight())
        .collect();
    let mut residual = weights.clone();
    let mut flow = vec![HashMap::<usize, f64>::new(); states.len()];
    let scale = weights.iter().copied().fold(1.0_f64, f64::max);
    let tolerance = 64.0 * f64::EPSILON * scale;

    loop {
        let mut best = None;
        for &(left, right) in edges {
            let transferable = residual[left].min(residual[right]);
            if transferable <= tolerance {
                continue;
            }
            let score = (transferable, residual[left] + residual[right]);
            if best
                .as_ref()
                .is_none_or(|(_, _, best_score)| score > *best_score)
            {
                best = Some((left, right, score));
            }
        }
        let Some((left, right, (transferable, _))) = best else {
            break;
        };
        *flow[left].entry(right).or_insert(0.0) += transferable;
        *flow[right].entry(left).or_insert(0.0) += transferable;
        residual[left] = (residual[left] - transferable).max(0.0);
        residual[right] = (residual[right] - transferable).max(0.0);
    }

    states
        .iter()
        .enumerate()
        .map(|(source_id, source)| {
            let mut row = Vec::new();
            for (&target_id, &path_weight) in &flow[source_id] {
                let target = states[target_id];
                row.push(ScatteringChoice {
                    new_kind: target.kind,
                    exit_leg: target.entrance,
                    next_delta: -target.delta,
                    probability: path_weight / weights[source_id],
                });
            }
            row.sort_by_key(|choice| (choice.new_kind, choice.exit_leg, choice.next_delta));
            if residual[source_id] > tolerance || row.is_empty() {
                // A bounce leaves the local vertex unchanged and sends the
                // head back through the link it entered from. For a
                // multi-level local space the worm charge is therefore
                // preserved. Using `-source.delta` is only accidentally safe
                // for a binary spin-1/2 flip and can create an illegal
                // raising/lowering operation at the linked leg for Spin-S.
                row.push(ScatteringChoice {
                    new_kind: source.kind,
                    exit_leg: source.entrance,
                    next_delta: source.delta,
                    probability: residual[source_id].max(0.0) / weights[source_id],
                });
            }
            normalize_row(&mut row);
            row
        })
        .collect()
}

fn metropolis_rows(
    kinds: &[VertexKind],
    leg_dimensions: &[usize],
    kind_lookup: &HashMap<Vec<BasisState>, usize>,
    lookup: &[Option<usize>],
    states: &[ExtendedState],
    leg_count: usize,
) -> Vec<Vec<ScatteringChoice>> {
    let proposal_count = (2 * leg_count) as f64;
    states
        .iter()
        .map(|&source| {
            let source_weight = kinds[source.kind].weight();
            let mut probabilities = HashMap::<(usize, usize, i8), f64>::new();
            let mut bounce = 0.0;
            for exit in 0..leg_count {
                for exit_delta in [-1_i8, 1_i8] {
                    let Some(new_kind) = transformed_kind(
                        kinds,
                        leg_dimensions,
                        kind_lookup,
                        source,
                        exit,
                        exit_delta,
                    ) else {
                        bounce += 1.0 / proposal_count;
                        continue;
                    };
                    let target_raw = raw_index(new_kind, exit, -exit_delta, leg_count);
                    let Some(target_id) = lookup[target_raw] else {
                        bounce += 1.0 / proposal_count;
                        continue;
                    };
                    let target = states[target_id];
                    let accept = (kinds[target.kind].weight() / source_weight).min(1.0);
                    if target == source {
                        bounce += 1.0 / proposal_count;
                    } else {
                        *probabilities
                            .entry((target.kind, target.entrance, -target.delta))
                            .or_insert(0.0) += accept / proposal_count;
                        bounce += (1.0 - accept) / proposal_count;
                    }
                }
            }
            let mut row: Vec<_> = probabilities
                .into_iter()
                .map(
                    |((new_kind, exit_leg, next_delta), probability)| ScatteringChoice {
                        new_kind,
                        exit_leg,
                        next_delta,
                        probability,
                    },
                )
                .collect();
            row.sort_by_key(|choice| (choice.new_kind, choice.exit_leg, choice.next_delta));
            row.push(ScatteringChoice {
                new_kind: source.kind,
                exit_leg: source.entrance,
                next_delta: source.delta,
                probability: bounce,
            });
            normalize_row(&mut row);
            row
        })
        .collect()
}

fn normalize_row(row: &mut [ScatteringChoice]) {
    let sum: f64 = row.iter().map(|choice| choice.probability).sum();
    if sum > 0.0 {
        for choice in row {
            choice.probability /= sum;
        }
    }
}

fn diagnostics(
    kinds: &[VertexKind],
    states: &[ExtendedState],
    rows: &[Vec<ScatteringChoice>],
) -> ScatteringDiagnostics {
    let max_row_error = rows
        .iter()
        .map(|row| {
            let sum: f64 = row.iter().map(|choice| choice.probability).sum();
            (sum - 1.0).abs()
        })
        .fold(0.0_f64, f64::max);
    let mut max_detailed_balance_error = 0.0_f64;
    let mut bounce_sum = 0.0;
    for (source_id, source) in states.iter().enumerate() {
        let source_weight = kinds[source.kind].weight();
        for choice in &rows[source_id] {
            let is_bounce = choice.new_kind == source.kind
                && choice.exit_leg == source.entrance
                && choice.next_delta == source.delta;
            if is_bounce {
                // Identity paths obey detailed balance with themselves. Their
                // continuation charge follows the global bounce convention,
                // rather than the reverse-edge encoding used by non-bounces.
                bounce_sum += choice.probability;
                continue;
            }
            let target = ExtendedState {
                kind: choice.new_kind,
                entrance: choice.exit_leg,
                delta: -choice.next_delta,
            };
            let Some(target_id) = states.iter().position(|state| *state == target) else {
                max_detailed_balance_error = f64::INFINITY;
                continue;
            };
            let reverse_probability = rows[target_id]
                .iter()
                .filter(|reverse| {
                    reverse.new_kind == source.kind
                        && reverse.exit_leg == source.entrance
                        && reverse.next_delta == -source.delta
                })
                .map(|reverse| reverse.probability)
                .sum::<f64>();
            let target_weight = kinds[target.kind].weight();
            max_detailed_balance_error = max_detailed_balance_error.max(
                (source_weight * choice.probability - target_weight * reverse_probability).abs(),
            );
        }
    }
    ScatteringDiagnostics {
        max_row_error,
        max_detailed_balance_error,
        mean_bounce_probability: bounce_sum / states.len().max(1) as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generic_spin_one_table_obeys_balance() {
        let mut kinds = Vec::new();
        for state in 0..3_u16 {
            kinds.push(
                VertexKind::new(format!("diag-{state}"), vec![state, state], 1.0).expect("kind"),
            );
        }
        kinds.push(VertexKind::new("raise", vec![0, 1], 0.7).expect("kind"));
        kinds.push(VertexKind::new("lower", vec![1, 0], 0.7).expect("kind"));
        let table =
            ScatteringTable::build(&kinds, &[3, 3], ScatteringPolicy::LowBounce).expect("table");
        let diagnostics = table.diagnostics();
        assert!(diagnostics.max_row_error < 1.0e-12);
        assert!(diagnostics.max_detailed_balance_error < 1.0e-12);
    }
}
