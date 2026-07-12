//! Magnetization observable.

use crate::lattice::interaction::{Hamiltonian, Measurable};
use crate::lattice::state::System;
use crate::observables::common::{MomentSpec, Observable};

pub struct Magnetization;

impl<H: Hamiltonian + Measurable> Observable<H> for Magnetization {
    fn name(&self) -> &str {
        "Magnetization"
    }

    fn measure(&self, system: &System, model: &H) -> f64 {
        model.magnetization(&system.spins)
    }

    fn moments(&self) -> &[MomentSpec] {
        static MOMENTS: [MomentSpec; 2] = [
            MomentSpec {
                name: "M2",
                order: 2,
            },
            MomentSpec {
                name: "M4",
                order: 4,
            },
        ];
        &MOMENTS
    }
}
