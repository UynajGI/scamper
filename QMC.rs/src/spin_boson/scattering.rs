//! Local directed-loop scattering rules.
//!
//! A local extended state is `(vertex_kind, entrance_leg)` with statistical
//! weight equal to the corresponding vertex weight.  A valid non-bounce edge
//! joins two extended states when flipping the two associated legs maps one
//! vertex kind into the other.  Directed-loop detailed balance is satisfied by
//! assigning a symmetric non-negative path weight to every undirected edge.
//!
//! The default [`ScatteringPolicy::LowBounce`] solver greedily saturates the
//! largest compatible residual pair.  It preserves exact local detailed
//! balance and is substantially less bouncy than a generic Metropolis table.
//! It is a correctness-first graph solver, not a proof of the global LP
//! minimum for every exotic catalog.  [`ScatteringPolicy::Metropolis`] remains
//! available as a simple reference implementation.

use rand::Rng;
use rand::RngExt;

use super::vertex::{VertexKind, LEGS_PER_VERTEX};

/// Strategy used to construct local directed-loop rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatteringPolicy {
    /// Symmetric residual-flow allocation with low bounce.
    LowBounce,
    /// Symmetric four-exit proposal followed by Metropolis accept/reject.
    Metropolis,
}

/// One outcome in a precomputed local scattering row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatteringChoice {
    /// Vertex kind after flipping the entrance and accepted exit legs.
    pub new_kind: usize,
    /// Exit leg. Equal to the entrance leg for a bounce.
    pub exit_leg: usize,
    /// Row-normalized probability.
    pub probability: f64,
}

/// Diagnostics for a local scattering table.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ScatteringDiagnostics {
    /// Maximum probability-row normalization error.
    pub max_row_error: f64,
    /// Maximum local detailed-balance residual.
    pub max_detailed_balance_error: f64,
    /// Mean bounce probability over all extended local states.
    pub mean_bounce_probability: f64,
}

/// Exact local directed-loop scattering table.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatteringTable {
    rows: Vec<Vec<ScatteringChoice>>,
    diagnostics: ScatteringDiagnostics,
    policy: ScatteringPolicy,
}

impl ScatteringTable {
    /// Build a table using the selected policy.
    pub fn build(kinds: &[VertexKind], policy: ScatteringPolicy) -> Self {
        match policy {
            ScatteringPolicy::LowBounce => Self::low_bounce(kinds),
            ScatteringPolicy::Metropolis => Self::metropolis(kinds),
        }
    }

    /// Build an exact low-bounce table from symmetric local path weights.
    ///
    /// Each extended state starts with residual weight `W_v`. At every step,
    /// the compatible pair with the largest transferable residual is assigned
    /// a symmetric flow. Any final residual is a bounce. Symmetric flow gives
    /// `W_e P(e->e') = W_e' P(e'->e)` by construction.
    pub fn low_bounce(kinds: &[VertexKind]) -> Self {
        let state_count = kinds.len() * LEGS_PER_VERTEX;
        let weights: Vec<f64> = (0..state_count)
            .map(|state| kinds[state / LEGS_PER_VERTEX].weight())
            .collect();
        let mut residual = weights.clone();
        let mut flow = vec![vec![0.0_f64; state_count]; state_count];
        let edges = compatible_edges(kinds);
        let scale = weights.iter().copied().fold(1.0_f64, f64::max);
        let tolerance = 64.0 * f64::EPSILON * scale;

        loop {
            let mut best: Option<(usize, usize, f64, f64)> = None;
            for &(left, right) in &edges {
                let transferable = residual[left].min(residual[right]);
                if transferable <= tolerance {
                    continue;
                }
                let combined = residual[left] + residual[right];
                let replace = match best {
                    None => true,
                    Some((_, _, best_transferable, best_combined)) => {
                        match transferable.total_cmp(&best_transferable) {
                            std::cmp::Ordering::Greater => true,
                            std::cmp::Ordering::Equal => combined > best_combined,
                            std::cmp::Ordering::Less => false,
                        }
                    }
                };
                if replace {
                    best = Some((left, right, transferable, combined));
                }
            }
            let Some((left, right, transferable, _)) = best else {
                break;
            };
            flow[left][right] += transferable;
            flow[right][left] += transferable;
            residual[left] = (residual[left] - transferable).max(0.0);
            residual[right] = (residual[right] - transferable).max(0.0);
        }

        let mut rows = Vec::with_capacity(state_count);
        for (state, state_flow) in flow.iter().enumerate() {
            let mut row = Vec::new();
            for (target, path_weight) in state_flow.iter().copied().enumerate() {
                if path_weight > tolerance {
                    row.push(ScatteringChoice {
                        new_kind: target / LEGS_PER_VERTEX,
                        exit_leg: target % LEGS_PER_VERTEX,
                        probability: path_weight / weights[state],
                    });
                }
            }
            if residual[state] > tolerance || row.is_empty() {
                row.push(ScatteringChoice {
                    new_kind: state / LEGS_PER_VERTEX,
                    exit_leg: state % LEGS_PER_VERTEX,
                    probability: residual[state].max(0.0) / weights[state],
                });
            }
            normalize_row(&mut row);
            rows.push(row);
        }
        Self::from_rows(kinds, rows, ScatteringPolicy::LowBounce)
    }

    /// Build the generic symmetric-proposal Metropolis reference table.
    pub fn metropolis(kinds: &[VertexKind]) -> Self {
        let mut rows = Vec::with_capacity(kinds.len() * LEGS_PER_VERTEX);
        for (kind_id, kind) in kinds.iter().enumerate() {
            for entrance in 0..LEGS_PER_VERTEX {
                let mut accepted = Vec::new();
                let mut bounce_probability = 0.0;
                for exit in 0..LEGS_PER_VERTEX {
                    if let Some(new_kind) = kind_after_flips(kinds, kind_id, entrance, exit) {
                        let ratio = kinds[new_kind].weight() / kind.weight();
                        let accept = ratio.min(1.0);
                        let probability = accept / LEGS_PER_VERTEX as f64;
                        if exit == entrance {
                            bounce_probability += probability;
                        } else if probability > 0.0 {
                            accepted.push(ScatteringChoice {
                                new_kind,
                                exit_leg: exit,
                                probability,
                            });
                        }
                        bounce_probability += (1.0 - accept) / LEGS_PER_VERTEX as f64;
                    } else {
                        bounce_probability += 1.0 / LEGS_PER_VERTEX as f64;
                    }
                }
                accepted.push(ScatteringChoice {
                    new_kind: kind_id,
                    exit_leg: entrance,
                    probability: bounce_probability,
                });
                normalize_row(&mut accepted);
                rows.push(accepted);
            }
        }
        Self::from_rows(kinds, rows, ScatteringPolicy::Metropolis)
    }

    fn from_rows(
        kinds: &[VertexKind],
        rows: Vec<Vec<ScatteringChoice>>,
        policy: ScatteringPolicy,
    ) -> Self {
        let max_row_error = rows
            .iter()
            .map(|row| {
                let sum: f64 = row.iter().map(|choice| choice.probability).sum();
                (sum - 1.0).abs()
            })
            .fold(0.0_f64, f64::max);
        let max_detailed_balance_error = detailed_balance_error(kinds, &rows);
        let bounce_sum: f64 = rows
            .iter()
            .enumerate()
            .map(|(state, row)| {
                let kind = state / LEGS_PER_VERTEX;
                let entrance = state % LEGS_PER_VERTEX;
                transition_probability(row, kind, entrance)
            })
            .sum();
        let n_rows = rows.len().max(1) as f64;
        Self {
            rows,
            diagnostics: ScatteringDiagnostics {
                max_row_error,
                max_detailed_balance_error,
                mean_bounce_probability: bounce_sum / n_rows,
            },
            policy,
        }
    }

    /// Row for `(kind, entrance)`.
    pub fn row(&self, kind: usize, entrance: usize) -> &[ScatteringChoice] {
        &self.rows[kind * LEGS_PER_VERTEX + entrance]
    }

    /// Draw one local scattering outcome.
    pub fn sample<R: Rng + ?Sized>(
        &self,
        kind: usize,
        entrance: usize,
        rng: &mut R,
    ) -> ScatteringChoice {
        let row = self.row(kind, entrance);
        let u = rng.random::<f64>();
        let mut cumulative = 0.0;
        for choice in row {
            cumulative += choice.probability;
            if u < cumulative {
                return *choice;
            }
        }
        row[row.len() - 1]
    }

    /// Construction diagnostics.
    pub fn diagnostics(&self) -> ScatteringDiagnostics {
        self.diagnostics
    }

    /// Construction policy.
    pub fn policy(&self) -> ScatteringPolicy {
        self.policy
    }
}

/// Find the local kind produced by flipping two legs.
pub fn kind_after_flips(
    kinds: &[VertexKind],
    kind: usize,
    entrance: usize,
    exit: usize,
) -> Option<usize> {
    let mut legs = *kinds[kind].legs();
    legs[entrance] *= -1;
    legs[exit] *= -1;
    kinds.iter().position(|candidate| candidate.legs() == &legs)
}

fn compatible_edges(kinds: &[VertexKind]) -> Vec<(usize, usize)> {
    let state_count = kinds.len() * LEGS_PER_VERTEX;
    let mut edges = Vec::new();
    for state in 0..state_count {
        let kind = state / LEGS_PER_VERTEX;
        let entrance = state % LEGS_PER_VERTEX;
        for exit in 0..LEGS_PER_VERTEX {
            let Some(new_kind) = kind_after_flips(kinds, kind, entrance, exit) else {
                continue;
            };
            let target = new_kind * LEGS_PER_VERTEX + exit;
            if target <= state {
                continue;
            }
            if kind_after_flips(kinds, new_kind, exit, entrance) == Some(kind) {
                edges.push((state, target));
            }
        }
    }
    edges
}

fn normalize_row(row: &mut [ScatteringChoice]) {
    let sum: f64 = row.iter().map(|choice| choice.probability).sum();
    if sum > 0.0 && (sum - 1.0).abs() > f64::EPSILON {
        for choice in row {
            choice.probability /= sum;
        }
    }
}

fn transition_probability(row: &[ScatteringChoice], new_kind: usize, exit: usize) -> f64 {
    row.iter()
        .filter(|choice| choice.new_kind == new_kind && choice.exit_leg == exit)
        .map(|choice| choice.probability)
        .sum()
}

fn detailed_balance_error(kinds: &[VertexKind], rows: &[Vec<ScatteringChoice>]) -> f64 {
    let mut maximum = 0.0_f64;
    for (kind, vertex) in kinds.iter().enumerate() {
        for entrance in 0..LEGS_PER_VERTEX {
            let row = &rows[kind * LEGS_PER_VERTEX + entrance];
            for choice in row {
                if choice.exit_leg == entrance && choice.new_kind == kind {
                    continue;
                }
                let forward = vertex.weight() * choice.probability;
                let reverse_row = &rows[choice.new_kind * LEGS_PER_VERTEX + choice.exit_leg];
                let reverse_probability = transition_probability(reverse_row, kind, entrance);
                let reverse = kinds[choice.new_kind].weight() * reverse_probability;
                maximum = maximum.max((forward - reverse).abs());
            }
        }
    }
    maximum
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spin_boson::vertex::VertexKind;

    fn test_kinds() -> Vec<VertexKind> {
        vec![
            VertexKind::new("down", [-1, -1, -1, -1], 2.0, true).expect("vertex"),
            VertexKind::new("mixed", [-1, 1, 1, -1], 0.5, false).expect("vertex"),
            VertexKind::new("up", [1, 1, 1, 1], 1.0, true).expect("vertex"),
        ]
    }

    #[test]
    fn both_policies_obey_local_detailed_balance() {
        for policy in [ScatteringPolicy::LowBounce, ScatteringPolicy::Metropolis] {
            let table = ScatteringTable::build(&test_kinds(), policy);
            assert!(table.diagnostics().max_row_error < 1.0e-12);
            assert!(table.diagnostics().max_detailed_balance_error < 1.0e-12);
        }
    }

    #[test]
    fn low_bounce_is_not_worse_for_reference_catalog() {
        let low = ScatteringTable::low_bounce(&test_kinds());
        let metro = ScatteringTable::metropolis(&test_kinds());
        assert!(
            low.diagnostics().mean_bounce_probability
                <= metro.diagnostics().mean_bounce_probability
        );
    }
}
