//! Metropolis-Hastings algorithm — generic over model and proposal strategy.

use crate::algorithms::proposal_strategy::ProposalStrategy;
use crate::algorithms::standard_strategy::StandardStrategy;
use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Metropolis-Hastings algorithm with configurable proposal strategy.
pub struct MetropolisCore<MC: ModelMC, S: ProposalStrategy<MC> = StandardStrategy> {
    model: MC,
    strategy: S,
    snapshot_interval: Option<u64>,
}

impl<MC: ModelMC> MetropolisCore<MC, StandardStrategy> {
    /// Create with standard (model-native) proposal strategy.
    pub fn new(model: MC) -> Self {
        MetropolisCore {
            model,
            strategy: StandardStrategy::new(),
            snapshot_interval: None,
        }
    }
}

impl<MC: ModelMC, S: ProposalStrategy<MC>> MetropolisCore<MC, S> {
    /// Create with a custom proposal strategy.
    pub fn with_strategy(model: MC, strategy: S) -> Self {
        MetropolisCore {
            model,
            strategy,
            snapshot_interval: None,
        }
    }

    /// Set the snapshot recording interval (in sweeps).
    pub fn with_snapshot_interval(mut self, interval: u64) -> Self {
        self.snapshot_interval = Some(interval);
        self
    }

    pub fn model(&self) -> &MC {
        &self.model
    }

    pub fn model_mut(&mut self) -> &mut MC {
        &mut self.model
    }
}

impl<MC: ModelMC, S: ProposalStrategy<MC>> MonteCarlo for MetropolisCore<MC, S> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        let dim = self.model.spin_dim();
        for _ in 0..n {
            let site = ctx.rng.random_range(0..n);
            let (old_spin, new_spin) = self.strategy.propose_flip(&self.model, site, &mut ctx.rng);
            let de = self
                .strategy
                .compute_delta_e(&self.model, site, &old_spin, &new_spin);
            if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
                let spins = self.model.spins_mut();
                for d in 0..dim {
                    spins[site * dim + d] = new_spin[d];
                }
            }
        }
        self.strategy.adapt_after_sweep(&mut self.model);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        let magnetization = self.model.magnetization();
        ctx.measure("Energy", energy);
        ctx.measure("Energy_Squared", energy * energy);
        ctx.measure("Magnetization", magnetization);
        ctx.measure("Magnetization_Squared", magnetization * magnetization);

        if let Some(interval) = self.snapshot_interval {
            if ctx.sweep_count() % interval == 0 {
                ctx.measure_array("Snapshot", &self.model.snapshot());
            }
        }
    }
}

impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams
    for MetropolisCore<MC, StandardStrategy>
{
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        Ok(MetropolisCore::new(mc))
    }
}

/// OPSS variant: MetropolisCore with OPSSStrategy for any O(N) model.
impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams
    for MetropolisCore<MC, crate::algorithms::OPSSStrategy>
{
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        let sigma = params.get::<f64>("opss_sigma").unwrap_or(60.0);
        Ok(MetropolisCore::with_strategy(
            mc,
            crate::algorithms::OPSSStrategy::new(sigma),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use crate::models::IsingModel;
    use rand::SeedableRng;

    #[test]
    fn test_metropolis_sweep() {
        let lattice = build_chain(8, true);
        let mut model = IsingModel::new(lattice, 1.0, 1.0);
        // Start from a disordered state (alternating spins) -- higher energy
        for (i, s) in model.spins_mut().iter_mut().enumerate() {
            *s = if i % 2 == 0 { 1.0 } else { -1.0 };
        }
        let mut core = MetropolisCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        let initial_energy = core.model.total_energy();

        for _ in 0..100 {
            core.sweep(&mut ctx);
        }

        let final_energy = core.model.total_energy();
        // At low temperature, energy should decrease toward ground state
        assert!(final_energy < initial_energy);
    }

    #[test]
    fn test_metropolis_high_temp_disorder() {
        let lattice = build_chain(16, true);
        let model = IsingModel::new(lattice, 0.01, 1.0); // Very high temperature
        let mut core = MetropolisCore::new(model);
        let rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut ctx = Context::new(rng, 100);

        // At high temperature, most flips should be accepted
        let initial_spins = core.model.spins().to_vec();
        for _ in 0..1000 {
            core.sweep(&mut ctx);
        }
        let final_spins = core.model.spins().to_vec();

        // At high T, should have significant spin flips
        let n_flipped: usize = initial_spins
            .iter()
            .zip(final_spins.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert!(n_flipped > 0);
    }
}
