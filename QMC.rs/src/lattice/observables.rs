//! Estimators for continuous-time spin-lattice configurations.

use super::configuration::{LatticeConfiguration, WorldlineIndex};
use super::error::LatticeQmcError;
use super::model::{SpinLatticeModel, TermLocation};

/// Common production observables for a spin-lattice run.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LatticeObservables {
    /// Uniform time-averaged magnetization per site.
    pub magnetization_z: f64,
    /// Absolute uniform magnetization.
    pub abs_magnetization_z: f64,
    /// Squared uniform magnetization.
    pub magnetization_z_squared: f64,
    /// Fourth moment.
    pub magnetization_z_fourth: f64,
    /// Gauge-staggered magnetization per site.
    pub staggered_magnetization_z: f64,
    /// Squared gauge-staggered magnetization.
    pub staggered_magnetization_z_squared: f64,
    /// Raw static uniform susceptibility `beta*N*<m_bar^2>`.
    pub susceptibility_z_raw: f64,
    /// Raw static gauge-staggered susceptibility.
    pub staggered_susceptibility_z_raw: f64,
    /// Exact expansion-order energy estimator.
    pub energy_total: f64,
    /// Energy per site.
    pub energy_per_site: f64,
    /// Total expansion order.
    pub expansion_order: f64,
    /// Diagonal expansion order.
    pub diagonal_order: f64,
    /// Off-diagonal expansion order.
    pub offdiagonal_order: f64,
    /// Vertices per site and unit imaginary time.
    pub vertex_density: f64,
    /// Exact imaginary-time average of `Sz_i Sz_j` over graph edges.
    pub nearest_neighbor_sz_correlation: f64,
    /// Fraction of sampled vertices on graph edges.
    pub edge_vertex_fraction: f64,
    /// Configuration sign. Positive sparse-operator models have sign one.
    pub average_sign: f64,
}

/// Measure common diagonal and expansion estimators.
pub fn measure_observables(
    configuration: &LatticeConfiguration,
    model: &SpinLatticeModel,
) -> Result<LatticeObservables, LatticeQmcError> {
    let index = WorldlineIndex::build(configuration, model)?;
    let site_magnetizations = site_time_averaged_magnetizations(configuration, model, &index);
    let n_sites = model.graph().site_count();
    let n_sites_f = n_sites as f64;
    let magnetization_z = site_magnetizations.iter().sum::<f64>() / n_sites_f;
    let staggered_magnetization_z = site_magnetizations
        .iter()
        .zip(model.gauge())
        .map(|(&magnetization, &phase)| magnetization * f64::from(phase))
        .sum::<f64>()
        / n_sites_f;
    let magnetization_z_squared = magnetization_z * magnetization_z;
    let staggered_magnetization_z_squared = staggered_magnetization_z * staggered_magnetization_z;
    let expansion_order = configuration.expansion_order() as f64;
    let diagonal_order = configuration.diagonal_order(model) as f64;
    let offdiagonal_order = expansion_order - diagonal_order;
    let energy_total = model.constant_shift() - expansion_order / configuration.beta();
    let edge_vertices = configuration
        .vertices()
        .iter()
        .filter(|vertex| matches!(model.term(vertex.term).location(), TermLocation::Edge(_)))
        .count() as f64;
    Ok(LatticeObservables {
        magnetization_z,
        abs_magnetization_z: magnetization_z.abs(),
        magnetization_z_squared,
        magnetization_z_fourth: magnetization_z_squared * magnetization_z_squared,
        staggered_magnetization_z,
        staggered_magnetization_z_squared,
        susceptibility_z_raw: configuration.beta() * n_sites_f * magnetization_z_squared,
        staggered_susceptibility_z_raw: configuration.beta()
            * n_sites_f
            * staggered_magnetization_z_squared,
        energy_total,
        energy_per_site: energy_total / n_sites_f,
        expansion_order,
        diagonal_order,
        offdiagonal_order,
        vertex_density: expansion_order / (configuration.beta() * n_sites_f),
        nearest_neighbor_sz_correlation: average_edge_correlation(configuration, model, &index),
        edge_vertex_fraction: if configuration.expansion_order() == 0 {
            0.0
        } else {
            edge_vertices / expansion_order
        },
        average_sign: 1.0,
    })
}

/// Exact imaginary-time average of `Sz` on every site.
pub fn site_time_averaged_magnetizations(
    configuration: &LatticeConfiguration,
    model: &SpinLatticeModel,
    index: &WorldlineIndex,
) -> Vec<f64> {
    (0..model.graph().site_count())
        .map(|site| {
            let events = index.events(site);
            let mut state = configuration.initial_states()[site];
            let mut previous_tau = 0.0;
            let mut integral = 0.0;
            for event in events {
                integral += (event.tau - previous_tau) * model.space().m(site, state);
                state = index.state_on_leg(configuration, model, event.outgoing_leg);
                previous_tau = event.tau;
            }
            integral += (configuration.beta() - previous_tau) * model.space().m(site, state);
            integral / configuration.beta()
        })
        .collect()
}

fn average_edge_correlation(
    configuration: &LatticeConfiguration,
    model: &SpinLatticeModel,
    index: &WorldlineIndex,
) -> f64 {
    if model.graph().edge_count() == 0 {
        return 0.0;
    }
    model
        .graph()
        .edges()
        .iter()
        .map(|edge| pair_time_average(configuration, model, index, edge.source, edge.target))
        .sum::<f64>()
        / model.graph().edge_count() as f64
}

fn pair_time_average(
    configuration: &LatticeConfiguration,
    model: &SpinLatticeModel,
    index: &WorldlineIndex,
    site_i: usize,
    site_j: usize,
) -> f64 {
    let events_i = index.events(site_i);
    let events_j = index.events(site_j);
    let mut position_i = 0;
    let mut position_j = 0;
    let mut state_i = configuration.initial_states()[site_i];
    let mut state_j = configuration.initial_states()[site_j];
    let mut previous_tau = 0.0;
    let mut integral = 0.0;

    while position_i < events_i.len() || position_j < events_j.len() {
        let next_i = events_i
            .get(position_i)
            .map_or(f64::INFINITY, |event| event.tau);
        let next_j = events_j
            .get(position_j)
            .map_or(f64::INFINITY, |event| event.tau);
        let next_tau = next_i.min(next_j);
        integral += (next_tau - previous_tau)
            * model.space().m(site_i, state_i)
            * model.space().m(site_j, state_j);

        while position_i < events_i.len() && events_i[position_i].tau.total_cmp(&next_tau).is_eq() {
            state_i = index.state_on_leg(configuration, model, events_i[position_i].outgoing_leg);
            position_i += 1;
        }
        while position_j < events_j.len() && events_j[position_j].tau.total_cmp(&next_tau).is_eq() {
            state_j = index.state_on_leg(configuration, model, events_j[position_j].outgoing_leg);
            position_j += 1;
        }
        previous_tau = next_tau;
    }
    integral += (configuration.beta() - previous_tau)
        * model.space().m(site_i, state_i)
        * model.space().m(site_j, state_j);
    integral / configuration.beta()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::CsrGraph;
    use crate::lattice::configuration::LatticeConfiguration;
    use crate::lattice::model::SpinLatticeModel;

    #[test]
    fn empty_spin_half_product_state_has_exact_magnetization() {
        let graph = CsrGraph::chain(4, true).expect("graph");
        let model = SpinLatticeModel::heisenberg(graph, 1, -1.0).expect("model");
        let configuration =
            LatticeConfiguration::new(3.0, vec![1, 1, 1, 1], &model).expect("configuration");
        let observables = measure_observables(&configuration, &model).expect("observables");
        assert!((observables.magnetization_z - 0.5).abs() < 1.0e-12);
        assert_eq!(observables.expansion_order, 0.0);
    }
}
