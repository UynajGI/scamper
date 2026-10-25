//! A2: external occupation-law extension.
//!
//! This test lives outside the library crate, so implementing
//! [`OccupationLaw`] here (using only public setters on
//! [`StaticConfiguration`]) proves a downstream process can supply its own
//! realization without editing CMC.rs. The generic adapter that drives such a
//! law is A3; here the law feeds the public analyzer directly.

use cmc_rs::{
    analyze, build_chain, ComponentWorkspace, CsrLattice, GraphView, OccupationLaw, SamplingError,
    StaticConfiguration, UndirectedGraphView,
};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Even-indexed sites occupied, every physical edge open.
#[derive(Debug, Default)]
struct AlternateSites;

impl OccupationLaw<CsrLattice> for AlternateSites {
    type Error = SamplingError;

    fn sample<R: Rng + ?Sized>(
        &mut self,
        graph: &CsrLattice,
        configuration: &mut StaticConfiguration,
        _rng: &mut R,
    ) -> Result<(), SamplingError> {
        for index in 0..graph.vertex_count() {
            let vertex = graph.vertex_id(index).expect("dense vertex index");
            configuration.set_vertex_active(vertex, index % 2 == 0);
        }
        for index in 0..graph.edge_count() {
            let edge = graph.edge_id(index).expect("dense edge index");
            configuration.set_edge_active(edge, true);
        }
        Ok(())
    }
}

#[test]
fn external_law_writes_a_configuration_the_public_analyzer_accepts() {
    let graph = build_chain(6, false);
    let mut law = AlternateSites;
    let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4133_0001);

    <AlternateSites as OccupationLaw<CsrLattice>>::validate_domain(
        &law,
        graph.vertex_count(),
        graph.edge_count(),
    )
    .expect("default domain check");
    law.sample(&graph, &mut configuration, &mut rng)
        .expect("infallible custom draw");

    let mut workspace = ComponentWorkspace::new();
    let result = analyze(&graph, &configuration, &configuration, &[], &mut workspace)
        .expect("analysis domains match the configuration");

    // Sites 0, 2, 4 are isolated: every chain edge touches an inactive site.
    assert_eq!(result.summary().active_vertex_count, 3);
    assert_eq!(result.summary().active_edge_count, 0);
    assert_eq!(result.summary().component_count, 3);
    assert_eq!(result.summary().largest_component_size, 1);
    assert_eq!(result.summary().raw_second_moment, 3);
}
