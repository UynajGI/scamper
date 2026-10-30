//! A4: external attack-order policy.
//!
//! This test lives outside the library crate, so implementing
//! [`AttackOrder`] here proves a downstream policy can drive the preset
//! adapter without editing CMC.rs. The built-in preset path is unchanged and
//! remains covered by the N3/T1 parity tests.

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{build_chain, AttackGrid, AttackOrder, AttackOrderContext, DegreeAttackMC, GraphView};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Descending dense ID: a maximally simple custom order with no RNG use.
#[derive(Debug, Default)]
struct ReverseId;

impl AttackOrder for ReverseId {
    fn name(&self) -> &str {
        "reverse-id"
    }

    fn order<R: rand::Rng + ?Sized>(&mut self, context: &mut AttackOrderContext<'_>, _rng: &mut R) {
        let vertices = context.graph().vertex_count();
        let order = context.order_mut();
        order.clear();
        order.extend((0..vertices).rev());
    }
}

#[test]
fn custom_policy_drives_the_preset_adapter() {
    let lattice = build_chain(6, false);
    let vertices = lattice.vertex_count();
    let mut adapter =
        DegreeAttackMC::new(lattice, ReverseId, AttackGrid::Every { step: 1 }).expect("adapter");
    assert_eq!(adapter.order_policy().name(), "reverse-id");

    let mut context =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x4134_0001), 0, 1);
    adapter.sweep(&mut context);
    assert_eq!(adapter.removal_order(), &[5, 4, 3, 2, 1, 0]);
    adapter.measure(&mut context);
    assert_eq!(adapter.curve_values().len(), vertices + 1);
    assert_eq!(adapter.curve_values()[0], 1.0);
    assert_eq!(adapter.curve_values()[vertices], 0.0);
}
