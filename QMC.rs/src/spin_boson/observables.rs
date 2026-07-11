//! Worldline estimators for the spin-boson impurity.

use rand::Rng;
use rand::RngExt;

use super::configuration::{event_spins, WorldlineIndex, WormholeConfiguration};
use super::error::SpinBosonError;
use super::model::SpinBosonModel;

/// One set of scalar measurements from a wormhole configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinBosonObservables {
    /// Imaginary-time averaged `sigma_z`.
    pub magnetization_sigma_z: f64,
    /// Imaginary-time averaged `S_z = sigma_z/2`.
    pub magnetization_s_z: f64,
    /// Squared averaged `sigma_z`.
    pub magnetization_sigma_z_squared: f64,
    /// Fourth power of averaged `sigma_z`.
    pub magnetization_sigma_z_fourth: f64,
    /// Static longitudinal susceptibility sample `beta * m_Sz^2`.
    pub susceptibility_z: f64,
    /// `sigma_z(beta/2) sigma_z(0)` averaged over random origins.
    pub correlation_sigma_z_half: f64,
    /// Corresponding `S_z` correlation.
    pub correlation_s_z_half: f64,
    /// Total expansion order.
    pub expansion_order: f64,
    /// Diagonal expansion order.
    pub diagonal_order: f64,
    /// Off-diagonal expansion order.
    pub offdiagonal_order: f64,
    /// Interaction-expansion energy estimator `-n/beta`.
    ///
    /// This includes arbitrary constant shifts introduced to keep vertex
    /// weights positive.  Model-specific physical energies should subtract
    /// those known shifts in a derived estimator.
    pub shifted_interaction_energy: f64,
}

/// Measure all built-in scalar observables.
pub fn measure_observables<R: Rng + ?Sized>(
    configuration: &WormholeConfiguration,
    model: &SpinBosonModel,
    correlation_samples: usize,
    rng: &mut R,
) -> Result<SpinBosonObservables, SpinBosonError> {
    let index = WorldlineIndex::build(configuration, model)?;
    let magnetization_sigma_z = integrated_sigma_z(configuration, model, &index);
    let magnetization_s_z = 0.5 * magnetization_sigma_z;
    let correlation_sigma_z_half = correlation_sigma_z(
        configuration,
        model,
        &index,
        0.5 * configuration.beta(),
        correlation_samples,
        rng,
    );
    let expansion_order = configuration.expansion_order() as f64;
    Ok(SpinBosonObservables {
        magnetization_sigma_z,
        magnetization_s_z,
        magnetization_sigma_z_squared: magnetization_sigma_z * magnetization_sigma_z,
        magnetization_sigma_z_fourth: magnetization_sigma_z.powi(4),
        susceptibility_z: configuration.beta() * magnetization_s_z * magnetization_s_z,
        correlation_sigma_z_half,
        correlation_s_z_half: 0.25 * correlation_sigma_z_half,
        expansion_order,
        diagonal_order: configuration.diagonal_order(model) as f64,
        offdiagonal_order: configuration.offdiagonal_order(model) as f64,
        shifted_interaction_energy: -expansion_order / configuration.beta(),
    })
}

/// Imaginary-time average of `sigma_z`.
pub fn integrated_sigma_z(
    configuration: &WormholeConfiguration,
    model: &SpinBosonModel,
    index: &WorldlineIndex,
) -> f64 {
    if index.events().is_empty() {
        return f64::from(configuration.empty_spin());
    }
    let first = index.events()[0];
    let mut spin = event_spins(configuration, model, first).0;
    let mut previous = 0.0;
    let mut total = 0.0;
    for event in index.events() {
        total += f64::from(spin) * (event.time - previous);
        spin = event_spins(configuration, model, *event).1;
        previous = event.time;
    }
    total += f64::from(spin) * (configuration.beta() - previous);
    total / configuration.beta()
}

/// Random-origin estimator of the longitudinal imaginary-time correlation.
pub fn correlation_sigma_z<R: Rng + ?Sized>(
    configuration: &WormholeConfiguration,
    model: &SpinBosonModel,
    index: &WorldlineIndex,
    delta_tau: f64,
    samples: usize,
    rng: &mut R,
) -> f64 {
    let sample_count = samples.max(1);
    let mut total = 0.0;
    for _ in 0..sample_count {
        let tau = rng.random::<f64>() * configuration.beta();
        let left = index.spin_at(configuration, model, tau);
        let right = index.spin_at(configuration, model, tau + delta_tau);
        total += f64::from(left * right);
    }
    total / sample_count as f64
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use crate::spin_boson::bath::{Bath, SingleModeBath};
    use crate::spin_boson::model::SpinBosonModel;

    use super::*;

    #[test]
    fn empty_worldline_estimators_are_exact() {
        let model = SpinBosonModel::jaynes_cummings(
            Bath::SingleMode(SingleModeBath::new(1.0).expect("mode")),
            0.2,
            0.0,
            None,
        )
        .expect("model");
        let configuration = WormholeConfiguration::new(4.0, -1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
        let observables =
            measure_observables(&configuration, &model, 8, &mut rng).expect("observables");
        assert!((observables.magnetization_sigma_z + 1.0).abs() < f64::EPSILON);
        assert!((observables.correlation_sigma_z_half - 1.0).abs() < f64::EPSILON);
        assert!(observables.expansion_order.abs() < f64::EPSILON);
    }
}
