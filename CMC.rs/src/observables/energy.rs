//! Energy observables.

use crate::lattice::interaction::Hamiltonian;
use crate::lattice::state::System;
use crate::observables::common::{MomentSpec, Observable};

pub struct TotalEnergy;

impl<H: Hamiltonian> Observable<H> for TotalEnergy {
    fn name(&self) -> &str {
        "Energy"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy
    }

    fn moments(&self) -> &[MomentSpec] {
        static MOMENTS: [MomentSpec; 1] = [MomentSpec {
            name: "E2",
            order: 2,
        }];
        &MOMENTS
    }
}

pub struct EnergyPerSite;

impl<H: Hamiltonian> Observable<H> for EnergyPerSite {
    fn name(&self) -> &str {
        "EnergyPerSite"
    }

    fn measure(&self, system: &System, _model: &H) -> f64 {
        system.energy / system.n_sites() as f64
    }
}
