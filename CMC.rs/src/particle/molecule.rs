//! Rigid-molecule topology, whole-molecule translations and plane rotations.
//!
//! An optional one-body [`DipolarExternalField`] couples molecular dipoles to
//! a uniform external field: the term participates in every trial's
//! Metropolis-Hastings weight and is available for measurement through
//! [`MolecularMetropolisCore::external_field_energy`]. `ParticleSystem::energy`
//! remains the pair energy; the one-body term is tracked by the kernel.

use crate::algorithms::SimulationPhase;
use crate::audit::{audit_particle_cache, should_audit_cache};
use crate::core::acceptance::{AcceptanceRule, MetropolisHastingsAcceptance};
use crate::core::trial::{metropolis_hastings_step, ProposedMove, TrialEvaluator};
use crate::core::visit::{SiteOrder, VisitSchedule};
use crate::particle::{
    CanonicalParticleKernel, MoveMixture, PairPotential, ParticleAlgorithm, ParticleBatchMove,
    ParticleBatchPatch, ParticleConfiguration, ParticleError, ParticleSystem, SimulationCell,
};
use carlo_rs::accept_log_probability;
use rand::{Rng, RngExt};

/// Fixed grouping of atoms into rigid molecules.
///
/// The topology controls collective moves only. Pair-potential exclusions and
/// bonded terms remain the responsibility of a future force-field layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MoleculeTopology {
    molecules: Vec<Vec<usize>>,
    atom_to_molecule: Vec<Option<usize>>,
}

impl MoleculeTopology {
    pub fn new(particle_count: usize, molecules: Vec<Vec<usize>>) -> Result<Self, ParticleError> {
        let mut atom_to_molecule = vec![None; particle_count];
        for (molecule, atoms) in molecules.iter().enumerate() {
            if atoms.is_empty() {
                return Err(ParticleError::InvalidTopology(format!(
                    "molecule {molecule} contains no atoms"
                )));
            }
            for &atom in atoms {
                if atom >= particle_count {
                    return Err(ParticleError::InvalidTopology(format!(
                        "molecule {molecule} contains out-of-range atom {atom}"
                    )));
                }
                if atom_to_molecule[atom].replace(molecule).is_some() {
                    return Err(ParticleError::InvalidTopology(format!(
                        "atom {atom} belongs to multiple molecules"
                    )));
                }
            }
        }
        Ok(Self {
            molecules,
            atom_to_molecule,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.molecules.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.molecules.is_empty()
    }

    #[inline]
    pub fn atoms(&self, molecule: usize) -> &[usize] {
        &self.molecules[molecule]
    }

    #[inline]
    pub fn molecule_of(&self, atom: usize) -> Option<usize> {
        self.atom_to_molecule[atom]
    }

    pub fn validate_particle_count(&self, particle_count: usize) -> Result<(), ParticleError> {
        if self.atom_to_molecule.len() == particle_count {
            Ok(())
        } else {
            Err(ParticleError::InvalidTopology(format!(
                "topology was built for {} particles but state contains {particle_count}",
                self.atom_to_molecule.len()
            )))
        }
    }

    /// Number of atoms the topology was built for (molecule members and
    /// free atoms alike).
    #[inline]
    pub fn particle_count(&self) -> usize {
        self.atom_to_molecule.len()
    }
}

/// One-body coupling of rigid molecules to a uniform external field.
///
/// Every atom carries a point charge `q_a`; the one-body energy of one
/// molecule is `-E · mu` with the molecular dipole
/// `mu = Σ_a q_a (r_a - r_anchor)`, where `r_anchor` is the molecule's first
/// atom and displacements use the periodic minimum image. For a neutral
/// molecule (`Σ q_a = 0`) this is invariant under wrapping and under the
/// choice of anchor, so translation and rotation moves see a well-defined
/// energy difference. The kernel constructor therefore rejects non-neutral
/// molecules: a net charge would couple to the absolute (wrapped) position
/// and silently break both periodicity and detailed balance.
#[derive(Debug, Clone, PartialEq)]
pub struct DipolarExternalField<const D: usize> {
    field: [f64; D],
    charges: Vec<f64>,
}

impl<const D: usize> DipolarExternalField<D> {
    /// Define the field vector and the per-atom charge table.
    ///
    /// The table must contain one entry per particle of the target system;
    /// molecule neutrality is validated when the field is attached to a
    /// [`MolecularMetropolisCore`] (the topology is needed there).
    pub fn new(field: [f64; D], charges: Vec<f64>) -> Result<Self, ParticleError> {
        if D == 0 {
            return Err(ParticleError::ZeroDimension);
        }
        if charges.is_empty() {
            return Err(ParticleError::InvalidPotential(
                "external-field charge table must not be empty".to_string(),
            ));
        }
        for (axis, &component) in field.iter().enumerate() {
            if !component.is_finite() {
                return Err(ParticleError::InvalidPotential(format!(
                    "external-field component {axis} is non-finite"
                )));
            }
        }
        for (atom, &charge) in charges.iter().enumerate() {
            if !charge.is_finite() {
                return Err(ParticleError::InvalidPotential(format!(
                    "external-field charge of atom {atom} is non-finite"
                )));
            }
        }
        Ok(Self { field, charges })
    }

    /// Uniform field vector `E`.
    #[inline]
    pub const fn field(&self) -> &[f64; D] {
        &self.field
    }

    /// Per-atom point charges.
    #[inline]
    pub fn charges(&self) -> &[f64] {
        &self.charges
    }
}

/// Symmetric whole-molecule translation proposal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidMoleculeTranslation {
    max_displacement: f64,
}

impl RigidMoleculeTranslation {
    pub fn new(max_displacement: f64) -> Result<Self, ParticleError> {
        if !max_displacement.is_finite() || max_displacement <= 0.0 {
            return Err(ParticleError::InvalidMove(
                "molecular translation scale must be finite and positive".to_string(),
            ));
        }
        Ok(Self { max_displacement })
    }

    #[inline]
    pub const fn max_displacement(self) -> f64 {
        self.max_displacement
    }

    pub fn propose<const D: usize>(
        &self,
        configuration: &ParticleConfiguration<D>,
        topology: &MoleculeTopology,
        molecule: usize,
        rng: &mut impl Rng,
    ) -> ProposedMove<ParticleBatchMove<D>> {
        let atoms = topology.atoms(molecule);
        let mut displacement = [0.0; D];
        for component in &mut displacement {
            *component = rng.random_range(-self.max_displacement..self.max_displacement);
        }
        let mut positions = Vec::with_capacity(atoms.len());
        for &atom in atoms {
            let mut position = *configuration.position(atom);
            for axis in 0..D {
                position[axis] += displacement[axis];
            }
            configuration.cell().wrap(&mut position);
            positions.push(position);
        }
        ProposedMove::symmetric(
            ParticleBatchMove::new(atoms.to_vec(), positions)
                .expect("validated molecule topology must produce a valid batch move"),
        )
    }
}

/// Symmetric rigid rotation in a uniformly selected Cartesian coordinate plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidMoleculeRotation {
    max_angle: f64,
}

impl RigidMoleculeRotation {
    pub fn new(max_angle: f64) -> Result<Self, ParticleError> {
        if !max_angle.is_finite() || max_angle <= 0.0 {
            return Err(ParticleError::InvalidMove(
                "molecular rotation scale must be finite and positive".to_string(),
            ));
        }
        Ok(Self { max_angle })
    }

    #[inline]
    pub const fn max_angle(self) -> f64 {
        self.max_angle
    }

    pub fn propose<const D: usize>(
        &self,
        configuration: &ParticleConfiguration<D>,
        topology: &MoleculeTopology,
        molecule: usize,
        rng: &mut impl Rng,
    ) -> ProposedMove<ParticleBatchMove<D>> {
        assert!(D >= 2, "rigid rotation requires at least two dimensions");
        let atoms = topology.atoms(molecule);
        let anchor = *configuration.position(atoms[0]);
        let mut unwrapped = Vec::with_capacity(atoms.len());
        let mut center = [0.0; D];
        for &atom in atoms {
            let relative = configuration
                .cell()
                .displacement(&anchor, configuration.position(atom));
            let mut position = anchor;
            for axis in 0..D {
                position[axis] += relative[axis];
                center[axis] += position[axis];
            }
            unwrapped.push(position);
        }
        for coordinate in &mut center {
            *coordinate /= atoms.len() as f64;
        }

        let first_axis = rng.random_range(0..D);
        let mut second_axis = rng.random_range(0..D - 1);
        if second_axis >= first_axis {
            second_axis += 1;
        }
        let angle = rng.random_range(-self.max_angle..self.max_angle);
        let (sine, cosine) = angle.sin_cos();

        let mut positions = Vec::with_capacity(atoms.len());
        for mut position in unwrapped {
            let left = position[first_axis] - center[first_axis];
            let right = position[second_axis] - center[second_axis];
            position[first_axis] = center[first_axis] + cosine * left - sine * right;
            position[second_axis] = center[second_axis] + sine * left + cosine * right;
            configuration.cell().wrap(&mut position);
            positions.push(position);
        }
        ProposedMove::symmetric(
            ParticleBatchMove::new(atoms.to_vec(), positions)
                .expect("validated molecule topology must produce a valid batch move"),
        )
    }
}

/// Definition of a local torsion around a bonded axis in three dimensions.
///
/// The rotating atoms must exclude the two axis atoms. Bonded energies and
/// topology exclusions are intentionally left to a future force-field layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TorsionDefinition {
    axis_start: usize,
    axis_end: usize,
    rotating_atoms: Vec<usize>,
}

impl TorsionDefinition {
    pub fn new(
        particle_count: usize,
        axis_start: usize,
        axis_end: usize,
        rotating_atoms: Vec<usize>,
    ) -> Result<Self, ParticleError> {
        if axis_start >= particle_count || axis_end >= particle_count {
            return Err(ParticleError::InvalidTopology(
                "torsion axis contains an out-of-range atom".to_string(),
            ));
        }
        if axis_start == axis_end {
            return Err(ParticleError::InvalidTopology(
                "torsion axis atoms must be distinct".to_string(),
            ));
        }
        if rotating_atoms.is_empty() {
            return Err(ParticleError::InvalidTopology(
                "torsion must rotate at least one atom".to_string(),
            ));
        }
        let mut seen = std::collections::BTreeSet::new();
        for &atom in &rotating_atoms {
            if atom >= particle_count {
                return Err(ParticleError::InvalidTopology(format!(
                    "torsion contains out-of-range rotating atom {atom}"
                )));
            }
            if atom == axis_start || atom == axis_end {
                return Err(ParticleError::InvalidTopology(
                    "torsion rotating atoms must exclude the axis atoms".to_string(),
                ));
            }
            if !seen.insert(atom) {
                return Err(ParticleError::InvalidTopology(format!(
                    "torsion contains duplicate rotating atom {atom}"
                )));
            }
        }
        Ok(Self {
            axis_start,
            axis_end,
            rotating_atoms,
        })
    }

    #[inline]
    pub const fn axis(&self) -> (usize, usize) {
        (self.axis_start, self.axis_end)
    }

    #[inline]
    pub fn rotating_atoms(&self) -> &[usize] {
        &self.rotating_atoms
    }
}

/// Symmetric local torsion proposal using Rodrigues rotation in three dimensions.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionRotation {
    definition: TorsionDefinition,
    max_angle: f64,
}

impl TorsionRotation {
    pub fn new(definition: TorsionDefinition, max_angle: f64) -> Result<Self, ParticleError> {
        if !max_angle.is_finite() || max_angle <= 0.0 {
            return Err(ParticleError::InvalidMove(
                "torsion angle scale must be finite and positive".to_string(),
            ));
        }
        Ok(Self {
            definition,
            max_angle,
        })
    }

    #[inline]
    pub const fn definition(&self) -> &TorsionDefinition {
        &self.definition
    }

    #[inline]
    pub const fn max_angle(&self) -> f64 {
        self.max_angle
    }

    pub fn propose(
        &self,
        configuration: &ParticleConfiguration<3>,
        rng: &mut impl Rng,
    ) -> Result<ProposedMove<ParticleBatchMove<3>>, ParticleError> {
        let (axis_start, axis_end) = self.definition.axis();
        if axis_start >= configuration.len() || axis_end >= configuration.len() {
            return Err(ParticleError::InvalidTopology(
                "torsion definition does not match the particle state".to_string(),
            ));
        }
        let anchor = *configuration.position(axis_start);
        let axis = configuration
            .cell()
            .displacement(&anchor, configuration.position(axis_end));
        let norm_squared = axis.iter().map(|value| value * value).sum::<f64>();
        if !norm_squared.is_finite() || norm_squared <= f64::EPSILON {
            return Err(ParticleError::InvalidMove(
                "torsion axis has zero or non-finite length".to_string(),
            ));
        }
        let inverse_norm = norm_squared.sqrt().recip();
        let unit = [
            axis[0] * inverse_norm,
            axis[1] * inverse_norm,
            axis[2] * inverse_norm,
        ];
        let angle = rng.random_range(-self.max_angle..self.max_angle);
        let (sine, cosine) = angle.sin_cos();

        let mut positions = Vec::with_capacity(self.definition.rotating_atoms.len());
        for &atom in &self.definition.rotating_atoms {
            if atom >= configuration.len() {
                return Err(ParticleError::InvalidTopology(format!(
                    "torsion rotating atom {atom} is out of range for the state"
                )));
            }
            let relative = configuration
                .cell()
                .displacement(&anchor, configuration.position(atom));
            let dot = unit[0] * relative[0] + unit[1] * relative[1] + unit[2] * relative[2];
            let cross = [
                unit[1] * relative[2] - unit[2] * relative[1],
                unit[2] * relative[0] - unit[0] * relative[2],
                unit[0] * relative[1] - unit[1] * relative[0],
            ];
            let mut position = anchor;
            for axis in 0..3 {
                let rotated = relative[axis] * cosine
                    + cross[axis] * sine
                    + unit[axis] * dot * (1.0 - cosine);
                position[axis] += rotated;
            }
            configuration.cell().wrap(&mut position);
            positions.push(position);
        }
        Ok(ProposedMove::symmetric(ParticleBatchMove::new(
            self.definition.rotating_atoms.clone(),
            positions,
        )?))
    }
}

/// Collective move kinds used by [`MolecularMetropolisCore`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolecularMoveKind {
    Translation,
    Rotation,
}

/// One collective rigid move per molecule and sweep.
#[derive(Debug, Clone)]
pub struct MolecularMetropolisCore<const D: usize> {
    topology: MoleculeTopology,
    translation: RigidMoleculeTranslation,
    rotation: RigidMoleculeRotation,
    external_field: Option<DipolarExternalField<D>>,
    mixture: MoveMixture<MolecularMoveKind>,
    order: SiteOrder,
    visit_schedule: VisitSchedule,
    patch: ParticleBatchPatch,
    energy_check_interval: u64,
    sweeps: u64,
    last_phase: Option<SimulationPhase>,
}

impl<const D: usize> MolecularMetropolisCore<D> {
    pub fn new(
        topology: MoleculeTopology,
        max_displacement: f64,
        max_angle: f64,
    ) -> Result<Self, ParticleError> {
        if D < 2 {
            return Err(ParticleError::InvalidMove(
                "molecular rotation kernel requires D >= 2".to_string(),
            ));
        }
        let mut mixture = MoveMixture::new();
        mixture.add(MolecularMoveKind::Translation, 1.0)?;
        mixture.add(MolecularMoveKind::Rotation, 1.0)?;
        Ok(Self {
            topology,
            translation: RigidMoleculeTranslation::new(max_displacement)?,
            rotation: RigidMoleculeRotation::new(max_angle)?,
            external_field: None,
            mixture,
            order: SiteOrder::new(),
            visit_schedule: VisitSchedule::RandomPermutation,
            patch: ParticleBatchPatch::default(),
            energy_check_interval: 0,
            sweeps: 0,
            last_phase: None,
        })
    }

    #[inline]
    pub const fn topology(&self) -> &MoleculeTopology {
        &self.topology
    }

    #[inline]
    pub const fn move_mixture(&self) -> &MoveMixture<MolecularMoveKind> {
        &self.mixture
    }

    #[inline]
    pub fn move_mixture_mut(&mut self) -> &mut MoveMixture<MolecularMoveKind> {
        &mut self.mixture
    }

    pub fn with_visit_schedule(mut self, schedule: VisitSchedule) -> Self {
        self.visit_schedule = schedule;
        self
    }

    /// Periodically audit energy and packed cell membership; zero disables it.
    pub fn with_energy_check_interval(mut self, interval: u64) -> Self {
        self.energy_check_interval = interval;
        self
    }

    /// Attach a one-body dipolar external field to every trial move.
    ///
    /// The charge table must match the topology's particle count and every
    /// molecule must be neutral (|Σ q| within roundoff); both are rejected
    /// loudly here rather than producing wrap-dependent energies silently.
    pub fn with_external_field(
        mut self,
        field: DipolarExternalField<D>,
    ) -> Result<Self, ParticleError> {
        let particle_count = self.topology.particle_count();
        if field.charges().len() != particle_count {
            return Err(ParticleError::InvalidPotential(format!(
                "external-field charge table covers {} atoms but the topology has {particle_count}",
                field.charges().len()
            )));
        }
        let charge_scale = field
            .charges()
            .iter()
            .fold(0.0_f64, |scale, &charge| scale.max(charge.abs()));
        let neutrality_tolerance = 1e-9 * (1.0 + charge_scale);
        for molecule in 0..self.topology.len() {
            let total: f64 = self
                .topology
                .atoms(molecule)
                .iter()
                .map(|&atom| field.charges()[atom])
                .sum();
            if total.abs() > neutrality_tolerance {
                return Err(ParticleError::InvalidPotential(format!(
                    "molecule {molecule} carries net charge {total}: a dipolar external field \
                     requires neutral molecules under periodic boundaries"
                )));
            }
        }
        self.external_field = Some(field);
        Ok(self)
    }

    /// The attached one-body external field, if any.
    #[inline]
    pub const fn external_field(&self) -> Option<&DipolarExternalField<D>> {
        self.external_field.as_ref()
    }

    /// Total one-body field energy `-Σ_m E·mu_m` of a configuration.
    ///
    /// `None` when no field is attached. Each molecular dipole is measured
    /// through minimum-image displacements from the molecule's first atom,
    /// so the value is wrap-invariant for neutral molecules. Atoms that
    /// belong to no molecule are never moved by this kernel and are
    /// deliberately excluded: their charge would only add a constant.
    pub fn external_field_energy(&self, configuration: &ParticleConfiguration<D>) -> Option<f64> {
        let field = self.external_field.as_ref()?;
        let cell = configuration.cell();
        let mut energy = 0.0;
        for molecule in 0..self.topology.len() {
            let atoms = self.topology.atoms(molecule);
            let anchor = configuration.position(atoms[0]);
            let mut dipole = [0.0; D];
            for &atom in atoms {
                let shift = cell.displacement(anchor, configuration.position(atom));
                let charge = field.charges()[atom];
                for axis in 0..D {
                    dipole[axis] += charge * shift[axis];
                }
            }
            energy -= field
                .field()
                .iter()
                .zip(dipole)
                .map(|(e, mu)| e * mu)
                .sum::<f64>();
        }
        Some(energy)
    }
}

impl<const D: usize> CanonicalParticleKernel for MolecularMetropolisCore<D> {}

impl<const D: usize, P: PairPotential> ParticleAlgorithm<D, P> for MolecularMetropolisCore<D> {
    fn sweep_with_phase(
        &mut self,
        system: &mut ParticleSystem<D>,
        potential: &P,
        rng: &mut impl Rng,
        phase: SimulationPhase,
    ) {
        self.topology
            .validate_particle_count(system.len())
            .expect("molecular topology/state size mismatch");
        if self.last_phase != Some(phase) {
            if phase == SimulationPhase::Measurement {
                self.mixture.freeze();
            }
            self.last_phase = Some(phase);
        }

        let ensemble = system.canonical_ensemble();
        let acceptance = MetropolisHastingsAcceptance;
        let molecules = self
            .order
            .prepare(self.topology.len(), self.visit_schedule, rng);
        for &molecule in molecules {
            let kind = *self.mixture.select(rng);
            let proposal = match kind {
                MolecularMoveKind::Translation => {
                    self.translation
                        .propose(system.configuration(), &self.topology, molecule, rng)
                }
                MolecularMoveKind::Rotation => {
                    self.rotation
                        .propose(system.configuration(), &self.topology, molecule, rng)
                }
            };
            if let Some(field) = &self.external_field {
                // Combined pair + one-body transition weight. The pair part
                // comes from the transactional evaluator; the one-body part is
                // the field energy change of the moved atoms, computed with
                // minimum-image displacements so wrapped rotations are exact.
                let mut delta =
                    <ParticleSystem<D> as TrialEvaluator<P, ParticleBatchMove<D>>>::evaluate_trial(
                        system,
                        potential,
                        &proposal.movement,
                        &mut self.patch,
                    );
                let atoms = self.topology.atoms(molecule);
                let cell = system.configuration().cell();
                let mut dipole_change = [0.0; D];
                for (slot, &atom) in atoms.iter().enumerate() {
                    let charge = field.charges()[atom];
                    if charge == 0.0 {
                        continue;
                    }
                    let shift = cell.displacement(
                        system.configuration().position(atom),
                        &proposal.movement.positions()[slot],
                    );
                    for axis in 0..D {
                        dipole_change[axis] += charge * shift[axis];
                    }
                }
                let field_change: f64 = field
                    .field()
                    .iter()
                    .zip(dipole_change)
                    .map(|(e, mu)| e * mu)
                    .sum();
                delta.energy -= field_change;
                let log_acceptance =
                    acceptance.log_acceptance(&ensemble, &delta, proposal.log_reverse_over_forward);
                assert!(
                    !log_acceptance.is_nan(),
                    "molecular trial produced NaN log acceptance"
                );
                if accept_log_probability(log_acceptance, rng) {
                    <ParticleSystem<D> as TrialEvaluator<P, ParticleBatchMove<D>>>::commit_trial(
                        system,
                        &proposal.movement,
                        &self.patch,
                    );
                }
            } else {
                metropolis_hastings_step(
                    system,
                    potential,
                    &proposal,
                    &ensemble,
                    &acceptance,
                    &mut self.patch,
                    rng,
                );
            }
        }

        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.energy_check_interval) {
            audit_particle_cache(system, potential).expect("molecular particle cache audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Rigid-molecule Metropolis-Hastings"
    }
}
