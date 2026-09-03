use core::fmt;

use rand::{Rng, RngExt};

use crate::UndirectedGraphView;

use super::super::activity::ActivityMask;
use super::super::{ProbabilityField, StaticConfiguration};

/// Domain-size mismatch detected before a Bernoulli sample mutates configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SamplingError {
    ConfigurationSize {
        expected_vertices: usize,
        actual_vertices: usize,
        expected_edges: usize,
        actual_edges: usize,
    },
    ProbabilityFieldLength {
        entity: &'static str,
        expected: usize,
        actual: usize,
    },
}

impl fmt::Display for SamplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ConfigurationSize {
                expected_vertices,
                actual_vertices,
                expected_edges,
                actual_edges,
            } => write!(
                formatter,
                "configuration has V={actual_vertices}, E={actual_edges}; graph requires V={expected_vertices}, E={expected_edges}"
            ),
            Self::ProbabilityFieldLength {
                entity,
                expected,
                actual,
            } => write!(
                formatter,
                "{entity} probability field has length {actual}; expected {expected}"
            ),
        }
    }
}

impl std::error::Error for SamplingError {}

/// Independent Bernoulli vertex activity. Every physical edge is active.
#[derive(Clone, Debug)]
pub struct SiteBernoulli<'a> {
    vertices: ProbabilityField<'a>,
}

impl<'a> SiteBernoulli<'a> {
    pub const fn new(vertices: ProbabilityField<'a>) -> Self {
        Self { vertices }
    }

    pub fn sample<G: UndirectedGraphView, R: Rng + ?Sized>(
        &self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        validate_configuration(graph, configuration)?;
        validate_field("vertex", &self.vertices, graph.vertex_count())?;
        sample_mask(&self.vertices, configuration.vertex_mask_mut(), rng);
        configuration.edge_mask_mut().fill(true);
        Ok(())
    }
}

/// Independent Bernoulli physical-edge activity. Every vertex is active.
#[derive(Clone, Debug)]
pub struct BondBernoulli<'a> {
    edges: ProbabilityField<'a>,
}

impl<'a> BondBernoulli<'a> {
    pub const fn new(edges: ProbabilityField<'a>) -> Self {
        Self { edges }
    }

    pub fn sample<G: UndirectedGraphView, R: Rng + ?Sized>(
        &self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        validate_configuration(graph, configuration)?;
        validate_field("edge", &self.edges, graph.edge_count())?;
        configuration.vertex_mask_mut().fill(true);
        sample_mask(&self.edges, configuration.edge_mask_mut(), rng);
        Ok(())
    }
}

/// Independent Bernoulli activity for both vertices and physical edges.
#[derive(Clone, Debug)]
pub struct MixedBernoulli<'a> {
    vertices: ProbabilityField<'a>,
    edges: ProbabilityField<'a>,
}

impl<'a> MixedBernoulli<'a> {
    pub const fn new(vertices: ProbabilityField<'a>, edges: ProbabilityField<'a>) -> Self {
        Self { vertices, edges }
    }

    pub fn sample<G: UndirectedGraphView, R: Rng + ?Sized>(
        &self,
        graph: &G,
        configuration: &mut StaticConfiguration,
        rng: &mut R,
    ) -> Result<(), SamplingError> {
        validate_configuration(graph, configuration)?;
        validate_field("vertex", &self.vertices, graph.vertex_count())?;
        validate_field("edge", &self.edges, graph.edge_count())?;
        sample_mask(&self.vertices, configuration.vertex_mask_mut(), rng);
        sample_mask(&self.edges, configuration.edge_mask_mut(), rng);
        Ok(())
    }
}

fn validate_configuration<G: UndirectedGraphView>(
    graph: &G,
    configuration: &StaticConfiguration,
) -> Result<(), SamplingError> {
    let expected_vertices = graph.vertex_count();
    let expected_edges = graph.edge_count();
    let actual_vertices = configuration.vertex_count();
    let actual_edges = configuration.edge_count();
    if (expected_vertices, expected_edges) == (actual_vertices, actual_edges) {
        Ok(())
    } else {
        Err(SamplingError::ConfigurationSize {
            expected_vertices,
            actual_vertices,
            expected_edges,
            actual_edges,
        })
    }
}

fn validate_field(
    entity: &'static str,
    field: &ProbabilityField<'_>,
    expected: usize,
) -> Result<(), SamplingError> {
    match field.len() {
        None => Ok(()),
        Some(actual) if actual == expected => Ok(()),
        Some(actual) => Err(SamplingError::ProbabilityFieldLength {
            entity,
            expected,
            actual,
        }),
    }
}

fn sample_mask<R: Rng + ?Sized>(
    field: &ProbabilityField<'_>,
    mask: &mut ActivityMask,
    rng: &mut R,
) {
    if let Some(probability) = field.uniform() {
        match probability.get() {
            0.0 => mask.fill(false),
            1.0 => mask.fill(true),
            probability => {
                for active in mask.dense_mut() {
                    *active = u8::from(rng.random::<f64>() < probability);
                }
            }
        }
        return;
    }

    for (index, active) in mask.dense_mut().iter_mut().enumerate() {
        let probability = field.at(index).get();
        *active = match probability {
            0.0 => 0,
            1.0 => 1,
            probability => u8::from(rng.random::<f64>() < probability),
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, BorrowedUndirectedCsr, EdgeActivity, Probability, VertexActivity};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    fn probability(value: f64) -> Probability {
        Probability::new(value).expect("test probability")
    }

    fn flags<G: UndirectedGraphView>(
        graph: &G,
        configuration: &StaticConfiguration,
    ) -> (Vec<bool>, Vec<bool>) {
        let vertices = graph
            .vertex_ids()
            .map(|vertex| configuration.vertex_active(vertex))
            .collect();
        let edges = (0..graph.edge_count())
            .map(|index| configuration.edge_active(graph.edge_id(index).unwrap()))
            .collect();
        (vertices, edges)
    }

    #[test]
    fn endpoints_use_complete_all_none_semantics() {
        let graph = build_chain(4, false);
        let mut configuration = StaticConfiguration::new(4, 3);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);

        SiteBernoulli::new(probability(0.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(
            flags(&graph, &configuration),
            (vec![false; 4], vec![true; 3])
        );

        BondBernoulli::new(probability(0.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(
            flags(&graph, &configuration),
            (vec![true; 4], vec![false; 3])
        );

        MixedBernoulli::new(probability(1.0).into(), probability(1.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(
            flags(&graph, &configuration),
            (vec![true; 4], vec![true; 3])
        );
    }

    #[test]
    fn fixed_seed_is_bitwise_reproducible() {
        let graph = build_chain(32, false);
        let law = MixedBernoulli::new(probability(0.37).into(), probability(0.61).into());
        let mut left = StaticConfiguration::new(32, 31);
        let mut right = StaticConfiguration::new(32, 31);
        let mut left_rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut right_rng = Xoshiro256PlusPlus::seed_from_u64(42);
        law.sample(&graph, &mut left, &mut left_rng).unwrap();
        law.sample(&graph, &mut right, &mut right_rng).unwrap();
        assert_eq!(flags(&graph, &left), flags(&graph, &right));
    }

    #[test]
    fn pure_laws_clear_stale_irrelevant_domains() {
        let graph = build_chain(8, false);
        let mut configuration = StaticConfiguration::new(8, 7);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
        MixedBernoulli::new(probability(0.0).into(), probability(0.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        SiteBernoulli::new(probability(1.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(
            flags(&graph, &configuration),
            (vec![true; 8], vec![true; 7])
        );

        MixedBernoulli::new(probability(0.0).into(), probability(0.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        BondBernoulli::new(probability(1.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(
            flags(&graph, &configuration),
            (vec![true; 8], vec![true; 7])
        );
    }

    #[test]
    fn mixed_reduces_to_pure_laws_with_matching_draw_order() {
        let graph = build_chain(24, false);
        let mut pure = StaticConfiguration::new(24, 23);
        let mut mixed = StaticConfiguration::new(24, 23);
        let mut pure_rng = Xoshiro256PlusPlus::seed_from_u64(8);
        let mut mixed_rng = Xoshiro256PlusPlus::seed_from_u64(8);
        SiteBernoulli::new(probability(0.43).into())
            .sample(&graph, &mut pure, &mut pure_rng)
            .unwrap();
        MixedBernoulli::new(probability(0.43).into(), probability(1.0).into())
            .sample(&graph, &mut mixed, &mut mixed_rng)
            .unwrap();
        assert_eq!(flags(&graph, &pure), flags(&graph, &mixed));

        let mut pure_rng = Xoshiro256PlusPlus::seed_from_u64(9);
        let mut mixed_rng = Xoshiro256PlusPlus::seed_from_u64(9);
        BondBernoulli::new(probability(0.57).into())
            .sample(&graph, &mut pure, &mut pure_rng)
            .unwrap();
        MixedBernoulli::new(probability(1.0).into(), probability(0.57).into())
            .sample(&graph, &mut mixed, &mut mixed_rng)
            .unwrap();
        assert_eq!(flags(&graph, &pure), flags(&graph, &mixed));
    }

    #[test]
    fn heterogeneous_fields_match_manual_reference() {
        let graph = build_chain(5, false);
        let vertex_values = [0.0, 0.2, 1.0, 0.7, 0.4].map(probability);
        let edge_values = [0.1, 1.0, 0.0, 0.8].map(probability);
        let law = MixedBernoulli::new(
            ProbabilityField::Borrowed(&vertex_values),
            ProbabilityField::Owned(edge_values.to_vec()),
        );
        let mut configuration = StaticConfiguration::new(5, 4);
        let mut actual_rng = Xoshiro256PlusPlus::seed_from_u64(17);
        law.sample(&graph, &mut configuration, &mut actual_rng)
            .unwrap();

        let mut reference_rng = Xoshiro256PlusPlus::seed_from_u64(17);
        let expected_vertices = vertex_values
            .iter()
            .map(|probability| match probability.get() {
                0.0 => false,
                1.0 => true,
                p => reference_rng.random::<f64>() < p,
            })
            .collect::<Vec<_>>();
        let expected_edges = edge_values
            .iter()
            .map(|probability| match probability.get() {
                0.0 => false,
                1.0 => true,
                p => reference_rng.random::<f64>() < p,
            })
            .collect::<Vec<_>>();
        assert_eq!(
            flags(&graph, &configuration),
            (expected_vertices, expected_edges)
        );
    }

    #[test]
    fn validation_failure_is_atomic_and_does_not_consume_rng() {
        let graph = build_chain(8, false);
        let mut configuration = StaticConfiguration::new(8, 7);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4154_4f4d_4943);
        MixedBernoulli::new(probability(0.4).into(), probability(0.6).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        let before = flags(&graph, &configuration);
        let vertices = vec![probability(0.5); 8];
        let edges = vec![probability(0.5); 6];
        let invalid = MixedBernoulli::new(
            ProbabilityField::Owned(vertices),
            ProbabilityField::Owned(edges),
        );
        let mut untouched_rng = rng.clone();
        assert_eq!(
            invalid.sample(&graph, &mut configuration, &mut rng),
            Err(SamplingError::ProbabilityFieldLength {
                entity: "edge",
                expected: 7,
                actual: 6,
            })
        );
        assert_eq!(flags(&graph, &configuration), before);
        for _ in 0..16 {
            assert_eq!(rng.random::<u64>(), untouched_rng.random::<u64>());
        }
    }
    #[test]
    fn uniform_endpoint_laws_do_not_consume_rng() {
        let graph = build_chain(8, false);
        for &p in &[0.0, 1.0] {
            let mut configuration = StaticConfiguration::new(8, 7);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5349_5445 + p as u64);
            let mut untouched_rng = rng.clone();
            SiteBernoulli::new(probability(p).into())
                .sample(&graph, &mut configuration, &mut rng)
                .unwrap();
            for _ in 0..16 {
                assert_eq!(rng.random::<u64>(), untouched_rng.random::<u64>());
            }
            let mut configuration = StaticConfiguration::new(8, 7);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x424f_4e44 + p as u64);
            let mut untouched_rng = rng.clone();
            BondBernoulli::new(probability(p).into())
                .sample(&graph, &mut configuration, &mut rng)
                .unwrap();
            for _ in 0..16 {
                assert_eq!(rng.random::<u64>(), untouched_rng.random::<u64>());
            }
        }
        for &(p_vertex, p_edge) in &[(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)] {
            let mut configuration = StaticConfiguration::new(8, 7);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                0x4d49_5845_4400 + (p_vertex as u64) * 2 + p_edge as u64,
            );
            let mut untouched_rng = rng.clone();
            MixedBernoulli::new(probability(p_vertex).into(), probability(p_edge).into())
                .sample(&graph, &mut configuration, &mut rng)
                .unwrap();
            for _ in 0..16 {
                assert_eq!(rng.random::<u64>(), untouched_rng.random::<u64>());
            }
        }
    }
    #[test]
    fn length_and_configuration_errors_precede_mutation() {
        let graph = build_chain(3, false);
        let probabilities = [probability(0.5)];
        let law = SiteBernoulli::new(ProbabilityField::Borrowed(&probabilities));
        let mut configuration = StaticConfiguration::new(3, 2);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        assert_eq!(
            law.sample(&graph, &mut configuration, &mut rng),
            Err(SamplingError::ProbabilityFieldLength {
                entity: "vertex",
                expected: 3,
                actual: 1,
            })
        );
        assert_eq!(
            flags(&graph, &configuration),
            (vec![false; 3], vec![false; 2])
        );

        configuration.resize(2, 2);
        let uniform = SiteBernoulli::new(probability(0.5).into());
        assert!(matches!(
            uniform.sample(&graph, &mut configuration, &mut rng),
            Err(SamplingError::ConfigurationSize { .. })
        ));
    }

    #[test]
    fn zero_size_graph_samples_without_rng_need() {
        let offsets = [0];
        let graph = BorrowedUndirectedCsr::new(&offsets, &[], &[], &[]).unwrap();
        let mut configuration = StaticConfiguration::new(0, 0);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        MixedBernoulli::new(probability(0.5).into(), probability(0.5).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(flags(&graph, &configuration), (Vec::new(), Vec::new()));
    }

    #[test]
    fn fresh_pure_laws_materialize_only_the_sampled_domain() {
        let graph = build_chain(16, false);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(6);

        let mut site_configuration = StaticConfiguration::new(16, 15);
        SiteBernoulli::new(probability(0.5).into())
            .sample(&graph, &mut site_configuration, &mut rng)
            .unwrap();
        assert_eq!(site_configuration.owned_mask_bytes(), 16);

        let mut bond_configuration = StaticConfiguration::new(16, 15);
        BondBernoulli::new(probability(0.5).into())
            .sample(&graph, &mut bond_configuration, &mut rng)
            .unwrap();
        assert_eq!(bond_configuration.owned_mask_bytes(), 15);
    }

    #[test]
    fn dense_capacity_survives_endpoint_samples() {
        let graph = build_chain(16, false);
        let mut configuration = StaticConfiguration::new(16, 15);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
        MixedBernoulli::new(probability(0.5).into(), probability(0.5).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        let bytes = configuration.owned_mask_bytes();
        MixedBernoulli::new(probability(0.0).into(), probability(1.0).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(configuration.owned_mask_bytes(), bytes);
        MixedBernoulli::new(probability(0.5).into(), probability(0.5).into())
            .sample(&graph, &mut configuration, &mut rng)
            .unwrap();
        assert_eq!(configuration.owned_mask_bytes(), bytes);
    }
}
