//! Local directed-loop scattering rules.
//!
//! A local extended state is `(vertex_kind, entrance_leg)` with statistical
//! weight equal to the corresponding vertex weight. A valid non-bounce edge
//! joins two extended states when flipping the two associated legs maps one
//! vertex kind into the other. Directed-loop detailed balance is satisfied by
//! assigning a symmetric non-negative path weight to every undirected edge.
//!
//! [`ScatteringPolicy::LowBounce`] decomposes the extended-state graph into
//! connected components. Complete components contain at most four states for a
//! four-leg spin-1/2 vertex, so their directed-loop equations have an analytic
//! minimum-bounce solution. Unexpected non-complete components use a symmetric
//! Metropolis flow as a defensive fallback.

use std::collections::{HashMap, VecDeque};

use rand::Rng;
use rand::RngExt;

use crate::impurity::core::local_hilbert::Spin;
use crate::impurity::core::operators::{VertexKind, LEGS_PER_VERTEX};
use crate::impurity::ImpurityError;

/// Strategy used to construct local directed-loop rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatteringPolicy {
    /// Analytic minimum-bounce flow on complete connected components.
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
    pub fn build(kinds: &[VertexKind], policy: ScatteringPolicy) -> Result<Self, ImpurityError> {
        validate_catalog(kinds)?;
        match policy {
            ScatteringPolicy::LowBounce => Self::low_bounce(kinds),
            ScatteringPolicy::Metropolis => Self::metropolis(kinds),
        }
    }

    /// Build the analytic minimum-bounce table where the local graph permits it.
    pub fn low_bounce(kinds: &[VertexKind]) -> Result<Self, ImpurityError> {
        validate_catalog(kinds)?;
        let state_count = kinds.len() * LEGS_PER_VERTEX;
        let weights: Vec<f64> = (0..state_count)
            .map(|state| kinds[state / LEGS_PER_VERTEX].weight())
            .collect();
        let adjacency = compatibility_graph(kinds);
        let components = connected_components(&adjacency);
        let mut rows = vec![Vec::new(); state_count];

        for component in components {
            let complete = component.iter().enumerate().all(|(left_pos, &left)| {
                component
                    .iter()
                    .enumerate()
                    .all(|(right_pos, &right)| left_pos == right_pos || adjacency[left][right])
            });
            let component_weights: Vec<f64> =
                component.iter().map(|&state| weights[state]).collect();
            let flows = if complete && component.len() <= LEGS_PER_VERTEX {
                minimum_bounce_complete_graph(&component_weights)
            } else {
                symmetric_metropolis_flows(&component_weights, &component, &adjacency)
            };
            append_component_rows(&mut rows, &component, &weights, &adjacency, &flows)?;
        }

        Self::from_rows(kinds, rows, ScatteringPolicy::LowBounce)
    }

    /// Build the generic symmetric-proposal Metropolis reference table.
    pub fn metropolis(kinds: &[VertexKind]) -> Result<Self, ImpurityError> {
        validate_catalog(kinds)?;
        let mut rows = Vec::with_capacity(kinds.len() * LEGS_PER_VERTEX);
        for (kind_id, kind) in kinds.iter().enumerate() {
            for entrance in 0..LEGS_PER_VERTEX {
                let mut accepted = Vec::new();
                let mut bounce_probability = 0.0;
                for exit in 0..LEGS_PER_VERTEX {
                    if let Some(new_kind) = kind_after_flips(kinds, kind_id, entrance, exit) {
                        let log_acceptance = kinds[new_kind].weight().ln() - kind.weight().ln();
                        let accept = log_acceptance.min(0.0).exp();
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
                rows.push(accepted);
            }
        }
        Self::from_rows(kinds, rows, ScatteringPolicy::Metropolis)
    }

    fn from_rows(
        kinds: &[VertexKind],
        rows: Vec<Vec<ScatteringChoice>>,
        policy: ScatteringPolicy,
    ) -> Result<Self, ImpurityError> {
        let expected_rows = kinds.len() * LEGS_PER_VERTEX;
        if rows.len() != expected_rows {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "scattering table has {} rows, expected {expected_rows}",
                rows.len()
            )));
        }

        let scale = kinds.iter().map(VertexKind::weight).fold(1.0_f64, f64::max);
        let row_tolerance = 1.0e-10;
        let detailed_balance_tolerance = 1.0e-10 * scale.max(1.0);
        let mut max_row_error = 0.0_f64;

        for (state, row) in rows.iter().enumerate() {
            if row.is_empty() {
                return Err(ImpurityError::InvalidConfiguration(format!(
                    "empty scattering row for extended state {state}"
                )));
            }
            for choice in row {
                if choice.new_kind >= kinds.len() || choice.exit_leg >= LEGS_PER_VERTEX {
                    return Err(ImpurityError::InvalidConfiguration(format!(
                        "out-of-range scattering target in row {state}"
                    )));
                }
                if !choice.probability.is_finite() || choice.probability < 0.0 {
                    return Err(ImpurityError::InvalidConfiguration(format!(
                        "invalid scattering probability in row {state}"
                    )));
                }
            }
            let sum: f64 = row.iter().map(|choice| choice.probability).sum();
            let error = (sum - 1.0).abs();
            max_row_error = max_row_error.max(error);
            if error > row_tolerance {
                return Err(ImpurityError::InvalidConfiguration(format!(
                    "scattering row {state} sums to {sum}, not one"
                )));
            }
        }

        let max_detailed_balance_error = detailed_balance_error(kinds, &rows);
        if max_detailed_balance_error > detailed_balance_tolerance {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "local detailed-balance residual {max_detailed_balance_error} exceeds \
                 {detailed_balance_tolerance}"
            )));
        }

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
        Ok(Self {
            rows,
            diagnostics: ScatteringDiagnostics {
                max_row_error,
                max_detailed_balance_error,
                mean_bounce_probability: bounce_sum / n_rows,
            },
            policy,
        })
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

fn validate_catalog(kinds: &[VertexKind]) -> Result<(), ImpurityError> {
    if kinds.is_empty() {
        return Err(ImpurityError::parameter(
            "vertex catalog",
            "at least one positive vertex kind is required",
        ));
    }
    let mut patterns: HashMap<[Spin; LEGS_PER_VERTEX], usize> = HashMap::new();
    for (kind_id, kind) in kinds.iter().enumerate() {
        if let Some(previous) = patterns.insert(*kind.legs(), kind_id) {
            return Err(ImpurityError::parameter(
                "vertex catalog",
                format!(
                    "duplicate leg pattern for kinds {previous} and {kind_id}: {:?}",
                    kind.legs()
                ),
            ));
        }
    }
    Ok(())
}

fn compatibility_graph(kinds: &[VertexKind]) -> Vec<Vec<bool>> {
    let state_count = kinds.len() * LEGS_PER_VERTEX;
    let mut adjacency = vec![vec![false; state_count]; state_count];
    #[allow(clippy::needless_range_loop)]
    for state in 0..state_count {
        let kind = state / LEGS_PER_VERTEX;
        let entrance = state % LEGS_PER_VERTEX;
        for exit in 0..LEGS_PER_VERTEX {
            if exit == entrance {
                continue;
            }
            let Some(new_kind) = kind_after_flips(kinds, kind, entrance, exit) else {
                continue;
            };
            let target = new_kind * LEGS_PER_VERTEX + exit;
            if kind_after_flips(kinds, new_kind, exit, entrance) == Some(kind) {
                adjacency[state][target] = true;
                adjacency[target][state] = true;
            }
        }
    }
    adjacency
}

fn connected_components(adjacency: &[Vec<bool>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; adjacency.len()];
    let mut components = Vec::new();
    for root in 0..adjacency.len() {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        let mut queue = VecDeque::from([root]);
        let mut component = Vec::new();
        while let Some(state) = queue.pop_front() {
            component.push(state);
            for (target, connected) in adjacency[state].iter().copied().enumerate() {
                if connected && !seen[target] {
                    seen[target] = true;
                    queue.push_back(target);
                }
            }
        }
        components.push(component);
    }
    components
}

fn append_component_rows(
    rows: &mut [Vec<ScatteringChoice>],
    component: &[usize],
    weights: &[f64],
    adjacency: &[Vec<bool>],
    flows: &[Vec<f64>],
) -> Result<(), ImpurityError> {
    for (left_pos, &state) in component.iter().enumerate() {
        let weight = weights[state];
        let tolerance = 1.0e-10 * weight.max(f64::MIN_POSITIVE);
        let row_sum: f64 = flows[left_pos].iter().sum();
        if (row_sum - weight).abs() > tolerance {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "directed-loop path weights sum to {row_sum}, expected {weight}"
            )));
        }
        for (right_pos, &target) in component.iter().enumerate() {
            let path_weight = flows[left_pos][right_pos];
            if !path_weight.is_finite() || path_weight < -tolerance {
                return Err(ImpurityError::InvalidConfiguration(
                    "non-finite or negative directed-loop path weight".into(),
                ));
            }
            if path_weight <= 0.0 {
                continue;
            }
            if state != target && !adjacency[state][target] {
                return Err(ImpurityError::InvalidConfiguration(
                    "directed-loop flow assigned to a forbidden transition".into(),
                ));
            }
            rows[state].push(ScatteringChoice {
                new_kind: target / LEGS_PER_VERTEX,
                exit_leg: if state == target {
                    state % LEGS_PER_VERTEX
                } else {
                    target % LEGS_PER_VERTEX
                },
                probability: path_weight / weight,
            });
        }
    }
    Ok(())
}

fn minimum_bounce_complete_graph(weights: &[f64]) -> Vec<Vec<f64>> {
    let n = weights.len();
    let mut flows = vec![vec![0.0; n]; n];
    if n == 1 {
        flows[0][0] = weights[0];
        return flows;
    }

    let total: f64 = weights.iter().sum();
    let (maximum_index, &maximum_weight) = weights
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1))
        .expect("non-empty complete component");
    let other_weight = total - maximum_weight;
    let tolerance = 1.0e-14 * total.max(1.0);
    if maximum_weight > other_weight + tolerance {
        for (index, &weight) in weights.iter().enumerate() {
            if index != maximum_index {
                set_symmetric(&mut flows, maximum_index, index, weight);
            }
        }
        flows[maximum_index][maximum_index] = maximum_weight - other_weight;
        return flows;
    }

    match n {
        2 => {
            let transferred = weights[0].min(weights[1]);
            set_symmetric(&mut flows, 0, 1, transferred);
            flows[0][0] = weights[0] - transferred;
            flows[1][1] = weights[1] - transferred;
        }
        3 => {
            let [first, second, third] = [weights[0], weights[1], weights[2]];
            set_symmetric(&mut flows, 0, 1, 0.5 * (first + second - third));
            set_symmetric(&mut flows, 0, 2, 0.5 * (first + third - second));
            set_symmetric(&mut flows, 1, 2, 0.5 * (second + third - first));
        }
        4 => {
            let [first, second, third, fourth] = [weights[0], weights[1], weights[2], weights[3]];
            let x = (0.5 * (first + second - third - fourth)).max(0.0);
            let y = (0.5 * (first - second + third - fourth)).max(0.0);
            set_symmetric(&mut flows, 0, 1, x);
            set_symmetric(&mut flows, 0, 2, y);
            set_symmetric(&mut flows, 0, 3, first - x - y);
            set_symmetric(
                &mut flows,
                1,
                2,
                0.5 * (first + second + third - fourth) - x - y,
            );
            set_symmetric(
                &mut flows,
                1,
                3,
                0.5 * (-first + second - third + fourth) + y,
            );
            set_symmetric(
                &mut flows,
                2,
                3,
                0.5 * (-first - second + third + fourth) + x,
            );
            if !flow_matrix_is_valid(&flows, weights) {
                return complete_graph_metropolis(weights);
            }
        }
        _ => return complete_graph_metropolis(weights),
    }

    for row in &mut flows {
        for value in row {
            if *value < 0.0 && *value > -tolerance {
                *value = 0.0;
            }
        }
    }
    flows
}

fn complete_graph_metropolis(weights: &[f64]) -> Vec<Vec<f64>> {
    let n = weights.len();
    let mut flows = vec![vec![0.0; n]; n];
    let degree = (n - 1).max(1) as f64;
    for left in 0..n {
        for right in (left + 1)..n {
            set_symmetric(
                &mut flows,
                left,
                right,
                weights[left].min(weights[right]) / degree,
            );
        }
    }
    for (index, &weight) in weights.iter().enumerate() {
        flows[index][index] = weight - flows[index].iter().sum::<f64>();
    }
    flows
}

fn symmetric_metropolis_flows(
    weights: &[f64],
    component: &[usize],
    adjacency: &[Vec<bool>],
) -> Vec<Vec<f64>> {
    let n = weights.len();
    let max_degree = component
        .iter()
        .map(|&state| {
            component
                .iter()
                .filter(|&&target| adjacency[state][target])
                .count()
        })
        .max()
        .unwrap_or(1)
        .max(1) as f64;
    let mut flows = vec![vec![0.0; n]; n];
    for left in 0..n {
        for right in (left + 1)..n {
            if adjacency[component[left]][component[right]] {
                set_symmetric(
                    &mut flows,
                    left,
                    right,
                    weights[left].min(weights[right]) / max_degree,
                );
            }
        }
    }
    for (index, &weight) in weights.iter().enumerate() {
        flows[index][index] = weight - flows[index].iter().sum::<f64>();
    }
    flows
}

fn set_symmetric(flows: &mut [Vec<f64>], left: usize, right: usize, value: f64) {
    flows[left][right] = value;
    flows[right][left] = value;
}

fn flow_matrix_is_valid(flows: &[Vec<f64>], weights: &[f64]) -> bool {
    let scale = weights.iter().copied().fold(1.0_f64, f64::max);
    let tolerance = 1.0e-10 * scale;
    flows.iter().zip(weights).all(|(row, &weight)| {
        row.iter()
            .all(|value| value.is_finite() && *value >= -tolerance)
            && (row.iter().sum::<f64>() - weight).abs() <= tolerance
    })
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
            let table = ScatteringTable::build(&test_kinds(), policy).expect("table");
            assert!(table.diagnostics().max_row_error < 1.0e-12);
            assert!(table.diagnostics().max_detailed_balance_error < 1.0e-12);
        }
    }

    #[test]
    fn low_bounce_is_not_worse_for_reference_catalog() {
        let low = ScatteringTable::low_bounce(&test_kinds()).expect("low-bounce table");
        let metro = ScatteringTable::metropolis(&test_kinds()).expect("Metropolis table");
        assert!(
            low.diagnostics().mean_bounce_probability
                <= metro.diagnostics().mean_bounce_probability
        );
    }

    #[test]
    fn equal_three_state_component_has_zero_bounce() {
        let flows = minimum_bounce_complete_graph(&[1.0, 1.0, 1.0]);
        for (index, row) in flows.iter().enumerate() {
            assert!(row[index].abs() < 1.0e-14);
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1.0e-14);
        }
    }

    #[test]
    fn near_zero_weight_four_state_component_remains_bounce_free() {
        let weights = [16.0 * f64::EPSILON, 0.05, 0.06, 0.012];
        let flows = minimum_bounce_complete_graph(&weights);
        for (index, row) in flows.iter().enumerate() {
            assert!(row[index].abs() < 1.0e-14);
            assert!((row.iter().sum::<f64>() - weights[index]).abs() < 1.0e-13);
        }
    }

    #[test]
    fn dominant_weight_has_the_theoretical_minimum_bounce() {
        let weights = [5.0, 1.0, 1.5, 0.5];
        let flows = minimum_bounce_complete_graph(&weights);
        let total_bounce: f64 = flows.iter().enumerate().map(|(i, row)| row[i]).sum();
        let expected = 2.0 * 5.0 - weights.iter().sum::<f64>();
        assert!((total_bounce - expected).abs() < 1.0e-14);
    }

    #[test]
    fn duplicate_patterns_are_rejected() {
        let kinds = vec![
            VertexKind::new("first", [1, 1, 1, 1], 1.0, true).expect("vertex"),
            VertexKind::new("second", [1, 1, 1, 1], 2.0, true).expect("vertex"),
        ];
        assert!(ScatteringTable::build(&kinds, ScatteringPolicy::LowBounce).is_err());
    }
}
