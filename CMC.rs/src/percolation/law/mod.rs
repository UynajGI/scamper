mod bernoulli;

pub use bernoulli::{BondBernoulli, MixedBernoulli, SamplingError, SiteBernoulli};

use rand::Rng;

use super::StaticConfiguration;
use crate::UndirectedGraphView;

/// Open extension point for a site/edge occupation law.
///
/// A law validates its own fields against the topology
/// ([`Self::validate_domain`]) and then draws one complete realization into a
/// reused [`StaticConfiguration`]. The trait is generic over the graph view so
/// the same law serves owned and borrowed topologies; it is statically
/// dispatched through the adapter's law parameter.
///
/// `&mut self` is deliberate: a correlated law may cache a reusable field
/// that is drawn per realization. Purely functional laws ignore the
/// mutability.
///
/// Continuum and long-range percolation are *processes*, not occupation
/// laws: their "edges" are discovered during the draw, so they implement the
/// process policy traits rather than this trait.
pub trait OccupationLaw<G: UndirectedGraphView> {
    /// Law-specific rejection type, surfaced by the adapter at construction.
    type Error;

    /// Validate every law-owned field length against the topology.
    ///
    /// Called once by the adapter constructor, before any draw. The default
    /// accepts arbitrary topology sizes for laws without sized fields.
    fn validate_domain(&self, _vertices: usize, _edges: usize) -> Result<(), Self::Error> {
        Ok(())
    }

    /// Draw one independent realization into `configuration`.
    ///
    /// Domain checks must run before the first RNG draw so a rejected call
    /// leaves both the configuration and the RNG stream untouched.
    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), Self::Error>;
}

impl<G: UndirectedGraphView> OccupationLaw<G> for SiteBernoulli<'_> {
    type Error = SamplingError;

    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        SiteBernoulli::sample(self, graph, configuration, rng)
    }
}

impl<G: UndirectedGraphView> OccupationLaw<G> for BondBernoulli<'_> {
    type Error = SamplingError;

    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        BondBernoulli::sample(self, graph, configuration, rng)
    }
}

impl<G: UndirectedGraphView> OccupationLaw<G> for MixedBernoulli<'_> {
    type Error = SamplingError;

    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        MixedBernoulli::sample(self, graph, configuration, rng)
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    use super::*;
    use crate::{
        build_chain, build_square, EdgeActivity, GraphView, Probability, ProbabilityField,
        StaticLaw, VertexActivity,
    };

    fn masks(
        graph: &impl UndirectedGraphView,
        configuration: &StaticConfiguration,
    ) -> (Vec<bool>, Vec<bool>) {
        let vertices = (0..graph.vertex_count())
            .map(|index| {
                configuration.vertex_active(graph.vertex_id(index).expect("dense vertex index"))
            })
            .collect();
        let edges = (0..graph.edge_count())
            .map(|index| configuration.edge_active(graph.edge_id(index).expect("dense edge index")))
            .collect();
        (vertices, edges)
    }

    fn draw_with_trait<L: OccupationLaw<crate::CsrLattice>>(
        law: &mut L,
        graph: &crate::CsrLattice,
        seed: u64,
    ) -> StaticConfiguration {
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        law.sample(graph, &mut configuration, &mut rng)
            .unwrap_or_else(|_| panic!("valid law must sample"));
        configuration
    }

    #[test]
    fn bernoulli_trait_path_is_bitwise_identical_to_the_inherent_draw() {
        let graph = build_chain(37, false);
        let probability = Probability::new(0.31).expect("probability");
        let mut rng_inherent = Xoshiro256PlusPlus::seed_from_u64(0x4132_0001);
        let mut inherent = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut law = SiteBernoulli::new(ProbabilityField::Uniform(probability));
        law.sample(&graph, &mut inherent, &mut rng_inherent)
            .expect("inherent draw");
        let via_trait = draw_with_trait(&mut law, &graph, 0x4132_0001);
        assert_eq!(masks(&graph, &inherent), masks(&graph, &via_trait));
    }

    #[test]
    fn static_law_trait_matches_the_fundamental_laws_bitwise() {
        let graph = build_chain(37, false);
        let probability = Probability::new(0.31).expect("probability");

        let mut static_law = StaticLaw::site(0.31).expect("law");
        let static_masks = masks(
            &graph,
            &draw_with_trait(&mut static_law, &graph, 0x4132_0002),
        );

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4132_0002);
        let mut inherent = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        SiteBernoulli::new(ProbabilityField::Uniform(probability))
            .sample(&graph, &mut inherent, &mut rng)
            .expect("inherent draw");
        assert_eq!(static_masks, masks(&graph, &inherent));
    }

    #[test]
    fn correlated_law_trait_matches_the_inherent_draw() {
        let graph = build_square(8, 8, true);
        let mut law = crate::CorrelatedSiteLaw::new(&[8, 8], 0.5, 0.4).expect("correlated law");

        let mut inherent = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4132_0003);
        law.sample(&graph, &mut inherent, &mut rng).expect("draw");

        let via_trait = draw_with_trait(&mut law, &graph, 0x4132_0003);
        assert_eq!(masks(&graph, &inherent), masks(&graph, &via_trait));
    }

    #[test]
    fn domain_rejections_precede_draws_and_name_the_field() {
        let graph = build_chain(12, false);
        let short = vec![Probability::new(0.5).expect("probability"); 11];
        let law = StaticLaw::site_heterogeneous(short);
        assert!(matches!(
            OccupationLaw::<crate::CsrLattice>::validate_domain(
                &law,
                graph.vertex_count(),
                graph.edge_count()
            ),
            Err(SamplingError::ProbabilityFieldLength {
                entity: "vertex",
                expected: 12,
                actual: 11,
            })
        ));

        let mut short_law = law;
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4132_0004);
        assert!(OccupationLaw::<crate::CsrLattice>::sample(
            &mut short_law,
            &graph,
            &mut configuration,
            &mut rng
        )
        .is_err());
        // A rejected draw must not advance the RNG: the stream still matches
        // a fresh generator at the same seed.
        let mut expected = Xoshiro256PlusPlus::seed_from_u64(0x4132_0004);
        assert_eq!(rng.next_u64(), expected.next_u64());
    }
}
