//! Sign-free spin-boson impurity model catalogs for wormhole QMC.

use std::collections::HashMap;

use crate::impurity::core::kernel::{
    validate_sign_free_channels, KernelDirection, PairFlipGauge, SignFreeMetadata, SignFreeReport,
};
use crate::impurity::core::local_hilbert::Spin;
use crate::impurity::core::operators::{BasisTransform, VertexKind, LEGS_PER_VERTEX};
use crate::impurity::spin_boson::bath::Bath;
use crate::impurity::spin_boson::wormhole::scattering::{ScatteringPolicy, ScatteringTable};
use crate::impurity::ImpurityError;

/// Normalization convention for rotating- and counter-rotating amplitudes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CouplingNormalization {
    /// `r = 1`, `c = crw_ratio`.
    #[default]
    FixedRw,
    /// `r + c = 1`.
    FixedTotal,
    /// `r^2 + c^2 = 1`.
    FixedQuadratic,
}

impl CouplingNormalization {
    /// Dimensionless amplitudes in `rho = g (r sigma_- + c sigma_+)`.
    pub fn amplitudes(self, crw_ratio: f64) -> (f64, f64) {
        match self {
            Self::FixedRw => (1.0, crw_ratio),
            Self::FixedTotal => {
                let denominator = 1.0 + crw_ratio;
                (1.0 / denominator, crw_ratio / denominator)
            }
            Self::FixedQuadratic => {
                let denominator = (1.0 + crw_ratio * crw_ratio).sqrt();
                (1.0 / denominator, crw_ratio / denominator)
            }
        }
    }
}

/// Supported spin-boson Hamiltonian families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImpurityModelKind {
    JaynesCummings,
    RwCrw,
    Xxz,
    Xyz,
    RotatedImpurity,
    Custom,
}

/// One independently sampled retarded interaction type.
#[derive(Debug, Clone, PartialEq)]
pub struct InteractionChannel {
    name: String,
    bath: Bath,
    direction: KernelDirection,
    kinds: Vec<VertexKind>,
    diagonal_lookup: HashMap<(Spin, Spin), usize>,
    scattering: ScatteringTable,
    diagonal_shift: f64,
    sign_free: SignFreeMetadata,
}

impl InteractionChannel {
    /// Construct a backward-compatible custom channel. Pair-flip channels are
    /// marked as having an unspecified gauge and therefore cannot be composed
    /// with another interaction until metadata is supplied explicitly.
    pub fn new(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
    ) -> Result<Self, ImpurityError> {
        Self::with_scattering_policy(name, bath, direction, kinds, ScatteringPolicy::LowBounce)
    }

    pub fn with_scattering_policy(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
        scattering_policy: ScatteringPolicy,
    ) -> Result<Self, ImpurityError> {
        let pair_flip = catalog_has_pair_flips(&kinds);
        let metadata = SignFreeMetadata::new(
            BasisTransform::identity(),
            if pair_flip {
                PairFlipGauge::Unspecified
            } else {
                PairFlipGauge::NotPresent
            },
        );
        Self::with_metadata_and_policy(
            name,
            bath,
            direction,
            kinds,
            0.0,
            metadata,
            scattering_policy,
        )
    }

    /// Construct a channel with an explicit constant shift and sign-free basis.
    pub fn with_metadata(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
        diagonal_shift: f64,
        sign_free: SignFreeMetadata,
    ) -> Result<Self, ImpurityError> {
        Self::with_metadata_and_policy(
            name,
            bath,
            direction,
            kinds,
            diagonal_shift,
            sign_free,
            ScatteringPolicy::LowBounce,
        )
    }

    pub fn with_metadata_and_policy(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
        diagonal_shift: f64,
        sign_free: SignFreeMetadata,
        scattering_policy: ScatteringPolicy,
    ) -> Result<Self, ImpurityError> {
        if kinds.is_empty() {
            return Err(ImpurityError::parameter(
                "vertex catalog",
                "an interaction channel needs at least one vertex kind",
            ));
        }
        if !diagonal_shift.is_finite() {
            return Err(ImpurityError::parameter("C", "must be finite"));
        }

        let mut pattern_lookup: HashMap<[Spin; LEGS_PER_VERTEX], usize> = HashMap::new();
        let mut diagonal_lookup = HashMap::new();
        for (kind_id, kind) in kinds.iter().enumerate() {
            if let Some(previous) = pattern_lookup.insert(*kind.legs(), kind_id) {
                return Err(ImpurityError::parameter(
                    "vertex catalog",
                    format!(
                        "duplicate leg pattern for kinds {previous} and {kind_id}: {:?}",
                        kind.legs()
                    ),
                ));
            }
            if kind.is_diagonal() {
                let key = (kind.spin(0), kind.spin(2));
                if let Some(previous) = diagonal_lookup.insert(key, kind_id) {
                    return Err(ImpurityError::parameter(
                        "vertex catalog",
                        format!(
                            "duplicate diagonal seed for ({},{}): kinds {previous} and {kind_id}",
                            key.0, key.1
                        ),
                    ));
                }
            }
        }
        for spin_a in [-1, 1] {
            for spin_b in [-1, 1] {
                if !diagonal_lookup.contains_key(&(spin_a, spin_b)) {
                    return Err(ImpurityError::parameter(
                        "vertex catalog",
                        format!("missing diagonal seed for ({spin_a},{spin_b})"),
                    ));
                }
            }
        }
        let scattering = ScatteringTable::build(&kinds, scattering_policy)?;
        Ok(Self {
            name: name.into(),
            bath,
            direction,
            kinds,
            diagonal_lookup,
            scattering,
            diagonal_shift,
            sign_free,
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn bath(&self) -> &Bath {
        &self.bath
    }
    pub fn direction(&self) -> KernelDirection {
        self.direction
    }
    pub fn kinds(&self) -> &[VertexKind] {
        &self.kinds
    }
    pub fn kind(&self, kind: usize) -> &VertexKind {
        &self.kinds[kind]
    }
    pub fn diagonal_kind(&self, spin_a: Spin, spin_b: Spin) -> usize {
        self.diagonal_lookup[&(spin_a, spin_b)]
    }
    pub fn scattering(&self) -> &ScatteringTable {
        &self.scattering
    }
    pub fn diagonal_shift(&self) -> f64 {
        self.diagonal_shift
    }
    pub fn sign_free_metadata(&self) -> SignFreeMetadata {
        self.sign_free
    }
}

/// Complete single-spin impurity model sampled by the wormhole engine.
#[derive(Debug, Clone, PartialEq)]
pub struct ImpurityModel {
    kind: ImpurityModelKind,
    name: String,
    interactions: Vec<InteractionChannel>,
    sign_free: SignFreeReport,
}

impl ImpurityModel {
    pub fn from_interactions(
        name: impl Into<String>,
        interactions: Vec<InteractionChannel>,
    ) -> Result<Self, ImpurityError> {
        let sign_free = validate_model_channels(&interactions)?;
        Ok(Self {
            kind: ImpurityModelKind::Custom,
            name: name.into(),
            interactions,
            sign_free,
        })
    }

    pub fn jaynes_cummings(
        bath: Bath,
        lambda: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, ImpurityError> {
        validate_nonnegative("lambda", lambda)?;
        let offdiagonal = if lambda > 0.0 {
            vec![("Splus_A_Sminus_B", [-1, 1, 1, -1], lambda)]
        } else {
            Vec::new()
        };
        let (kinds, shift) = build_catalog(0.0, h_z, constant, &offdiagonal)?;
        let interaction = InteractionChannel::with_metadata(
            "jc",
            bath,
            KernelDirection::Directed,
            kinds,
            shift,
            SignFreeMetadata::default(),
        )?;
        Self::from_single(
            ImpurityModelKind::JaynesCummings,
            "JaynesCummings",
            interaction,
        )
    }

    /// Directed `rho^dagger(tau_a) rho(tau_b)` interaction with
    /// `rho = g (r S_- + c S_+)`.
    pub fn rw_crw(
        bath: Bath,
        vertex_scale: f64,
        crw_ratio: f64,
        tunnelling: f64,
        normalization: CouplingNormalization,
        constant: Option<f64>,
    ) -> Result<Self, ImpurityError> {
        validate_nonnegative("vertex_scale", vertex_scale)?;
        validate_nonnegative("crw_ratio", crw_ratio)?;
        if !tunnelling.is_finite() {
            return Err(ImpurityError::parameter("tunnelling", "must be finite"));
        }
        let (kinds, shift) =
            build_rw_crw_catalog(vertex_scale, crw_ratio, tunnelling, normalization, constant)?;
        let interaction = InteractionChannel::with_metadata(
            "rw_crw",
            bath,
            KernelDirection::Directed,
            kinds,
            shift,
            SignFreeMetadata::default(),
        )?;
        Self::from_single(ImpurityModelKind::RwCrw, "RwCrwImpurity", interaction)
    }

    pub fn xxz(
        bath: Bath,
        lambda_xy: f64,
        lambda_z: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, ImpurityError> {
        validate_nonnegative("lambda_xy", lambda_xy)?;
        validate_nonnegative("lambda_z", lambda_z)?;
        let exchange = 0.5 * lambda_xy;
        let mut offdiagonal = Vec::new();
        if exchange > 0.0 {
            offdiagonal.push(("Sminus_A_Splus_B", [1, -1, -1, 1], exchange));
            offdiagonal.push(("Splus_A_Sminus_B", [-1, 1, 1, -1], exchange));
        }
        let (kinds, shift) = build_catalog(lambda_z, h_z, constant, &offdiagonal)?;
        let interaction = InteractionChannel::with_metadata(
            "xxz",
            bath,
            KernelDirection::Symmetric,
            kinds,
            shift,
            SignFreeMetadata::default(),
        )?;
        Self::from_single(ImpurityModelKind::Xxz, "XxzImpurity", interaction)
    }

    pub fn xyz(
        bath: Bath,
        lambda_x: f64,
        lambda_y: f64,
        lambda_z: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, ImpurityError> {
        validate_nonnegative("lambda_x", lambda_x)?;
        validate_nonnegative("lambda_y", lambda_y)?;
        validate_nonnegative("lambda_z", lambda_z)?;
        let exchange = 0.25 * (lambda_x + lambda_y);
        let pair = 0.25 * (lambda_x - lambda_y).abs();
        let mut offdiagonal = Vec::new();
        if exchange > 0.0 {
            offdiagonal.push(("Sminus_A_Splus_B", [1, -1, -1, 1], exchange));
            offdiagonal.push(("Splus_A_Sminus_B", [-1, 1, 1, -1], exchange));
        }
        if pair > 0.0 {
            offdiagonal.push(("Splus_A_Splus_B", [-1, 1, -1, 1], pair));
            offdiagonal.push(("Sminus_A_Sminus_B", [1, -1, 1, -1], pair));
        }
        let (kinds, shift) = build_catalog(lambda_z, h_z, constant, &offdiagonal)?;
        let (basis, gauge) = if pair == 0.0 {
            (BasisTransform::identity(), PairFlipGauge::NotPresent)
        } else if lambda_x >= lambda_y {
            (BasisTransform::identity(), PairFlipGauge::Positive)
        } else {
            (BasisTransform::swap_xy_gauge(), PairFlipGauge::Negative)
        };
        let interaction = InteractionChannel::with_metadata(
            "xyz",
            bath,
            KernelDirection::Symmetric,
            kinds,
            shift,
            SignFreeMetadata::new(basis, gauge),
        )?;
        Self::from_single(ImpurityModelKind::Xyz, "XyzImpurity", interaction)
    }

    /// Original longitudinal Rabi/spin-boson model represented in a basis where
    /// sampled `S_z` is physical `S_x`.
    pub fn rotated_impurity(
        bath: Bath,
        lambda: f64,
        tunnelling: f64,
        constant: Option<f64>,
    ) -> Result<Self, ImpurityError> {
        validate_nonnegative("lambda", lambda)?;
        let exchange = 0.25 * lambda;
        let pair = 0.25 * lambda;
        let mut offdiagonal = Vec::new();
        if exchange > 0.0 {
            offdiagonal.push(("Sminus_A_Splus_B", [1, -1, -1, 1], exchange));
            offdiagonal.push(("Splus_A_Sminus_B", [-1, 1, 1, -1], exchange));
            offdiagonal.push(("Splus_A_Splus_B", [-1, 1, -1, 1], pair));
            offdiagonal.push(("Sminus_A_Sminus_B", [1, -1, 1, -1], pair));
        }
        let (kinds, shift) = build_catalog(0.0, tunnelling, constant, &offdiagonal)?;
        let metadata = SignFreeMetadata::new(
            BasisTransform::rotated_rabi(),
            if pair > 0.0 {
                PairFlipGauge::Positive
            } else {
                PairFlipGauge::NotPresent
            },
        );
        let interaction = InteractionChannel::with_metadata(
            "rotated_impurity",
            bath,
            KernelDirection::Symmetric,
            kinds,
            shift,
            metadata,
        )?;
        Self::from_single(
            ImpurityModelKind::RotatedImpurity,
            "RotatedImpurity",
            interaction,
        )
    }

    fn from_single(
        kind: ImpurityModelKind,
        name: impl Into<String>,
        interaction: InteractionChannel,
    ) -> Result<Self, ImpurityError> {
        let interactions = vec![interaction];
        let sign_free = validate_model_channels(&interactions)?;
        Ok(Self {
            kind,
            name: name.into(),
            interactions,
            sign_free,
        })
    }

    pub fn kind(&self) -> ImpurityModelKind {
        self.kind
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn interactions(&self) -> &[InteractionChannel] {
        &self.interactions
    }
    pub fn interaction(&self, interaction: usize) -> &InteractionChannel {
        &self.interactions[interaction]
    }
    pub fn interaction_count(&self) -> usize {
        self.interactions.len()
    }
    pub fn basis_transform(&self) -> BasisTransform {
        self.sign_free.basis
    }
    pub fn sign_free_report(&self) -> SignFreeReport {
        self.sign_free
    }
    pub fn total_diagonal_shift(&self) -> f64 {
        self.interactions
            .iter()
            .map(InteractionChannel::diagonal_shift)
            .sum()
    }
    /// Correct `-<n>/beta`, which is measured in the shifted expansion, back
    /// to the spin-plus-coupling energy by adding the model's constant shifts.
    pub fn corrected_spin_coupling_energy(&self, expansion_order: f64, beta: f64) -> f64 {
        -expansion_order / beta + self.total_diagonal_shift()
    }
}

fn validate_model_channels(
    interactions: &[InteractionChannel],
) -> Result<SignFreeReport, ImpurityError> {
    let metadata: Vec<_> = interactions
        .iter()
        .map(|channel| (channel.name(), channel.sign_free_metadata()))
        .collect();
    validate_sign_free_channels(&metadata)
}

fn catalog_has_pair_flips(kinds: &[VertexKind]) -> bool {
    kinds.iter().any(|kind| {
        !kind.is_diagonal()
            && kind.spin(0) != kind.spin(1)
            && kind.spin(2) != kind.spin(3)
            && kind.spin(0) == kind.spin(2)
            && kind.spin(1) == kind.spin(3)
    })
}

type OffDiagonalSpec<'a> = (&'a str, [Spin; LEGS_PER_VERTEX], f64);

fn build_catalog(
    lambda_z: f64,
    h_z: f64,
    constant: Option<f64>,
    offdiagonal: &[OffDiagonalSpec<'_>],
) -> Result<(Vec<VertexKind>, f64), ImpurityError> {
    let maximum_offdiagonal = offdiagonal
        .iter()
        .map(|(_, _, weight)| *weight)
        .fold(0.0_f64, f64::max);
    let minimum_base = [-1, 1]
        .into_iter()
        .flat_map(|spin_a| {
            [-1, 1].into_iter().map(move |spin_b| {
                0.25 * lambda_z * f64::from(spin_a * spin_b)
                    + 0.25 * h_z * f64::from(spin_a + spin_b)
            })
        })
        .fold(f64::INFINITY, f64::min);
    let automatic = -minimum_base
        + maximum_offdiagonal
            .max(0.25 * lambda_z)
            .max(0.5 * h_z.abs())
            .max(1.0e-8);
    let shift = constant.unwrap_or(automatic);
    if !shift.is_finite() {
        return Err(ImpurityError::parameter("C", "must be finite"));
    }

    let mut kinds = Vec::with_capacity(4 + offdiagonal.len());
    for spin_a in [-1, 1] {
        for spin_b in [-1, 1] {
            let weight = shift
                + 0.25 * lambda_z * f64::from(spin_a * spin_b)
                + 0.25 * h_z * f64::from(spin_a + spin_b);
            kinds.push(VertexKind::new(
                format!("diag_{spin_a:+}_{spin_b:+}"),
                [spin_a, spin_a, spin_b, spin_b],
                weight,
                true,
            )?);
        }
    }
    for (name, legs, weight) in offdiagonal {
        if *weight > 0.0 {
            kinds.push(VertexKind::new(*name, *legs, *weight, false)?);
        }
    }
    Ok((kinds, shift))
}

fn build_rw_crw_catalog(
    vertex_scale: f64,
    crw_ratio: f64,
    tunnelling: f64,
    normalization: CouplingNormalization,
    constant: Option<f64>,
) -> Result<(Vec<VertexKind>, f64), ImpurityError> {
    let diagonal_constant =
        constant.unwrap_or_else(|| 0.5 * tunnelling.abs() + 16.0 * f64::EPSILON);
    if !diagonal_constant.is_finite() {
        return Err(ImpurityError::parameter("C", "must be finite"));
    }

    let mut kinds = Vec::with_capacity(8);
    for spin_a in [-1, 1] {
        for spin_b in [-1, 1] {
            let weight = diagonal_constant + 0.25 * tunnelling * f64::from(spin_a + spin_b);
            kinds.push(VertexKind::new(
                format!("diag_{spin_a:+}_{spin_b:+}"),
                [spin_a, spin_a, spin_b, spin_b],
                weight,
                true,
            )?);
        }
    }

    let (rotating, counter_rotating) = normalization.amplitudes(crw_ratio);
    for spin_a in [-1, 1] {
        for spin_b in [-1, 1] {
            // rho^dagger at A = r S+ + c S-, rho at B = r S- + c S+.
            let amplitude_a = if spin_a == -1 {
                rotating
            } else {
                counter_rotating
            };
            let amplitude_b = if spin_b == 1 {
                rotating
            } else {
                counter_rotating
            };
            let weight = vertex_scale * amplitude_a * amplitude_b;
            if weight > 0.0 {
                kinds.push(VertexKind::new(
                    format!("rw_crw_{spin_a:+}_{spin_b:+}"),
                    [spin_a, -spin_a, spin_b, -spin_b],
                    weight,
                    false,
                )?);
            }
        }
    }
    Ok((kinds, diagonal_constant))
}

fn validate_nonnegative(field: &str, value: f64) -> Result<(), ImpurityError> {
    if !value.is_finite() || value < 0.0 {
        return Err(ImpurityError::parameter(
            field,
            format!("must be finite and non-negative, got {value}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impurity::spin_boson::bath::SingleModeBath;

    fn mode() -> Bath {
        Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"))
    }

    fn kind_weight(model: &ImpurityModel, legs: [Spin; LEGS_PER_VERTEX]) -> Option<f64> {
        model
            .interaction(0)
            .kinds()
            .iter()
            .find(|kind| kind.legs() == &legs)
            .map(VertexKind::weight)
    }

    #[test]
    fn pure_rotating_rw_crw_catalog_matches_jaynes_cummings() {
        let rw = ImpurityModel::rw_crw(
            mode(),
            0.4,
            0.0,
            0.1,
            CouplingNormalization::FixedRw,
            Some(0.6),
        )
        .expect("rw");
        let jc = ImpurityModel::jaynes_cummings(mode(), 0.4, 0.1, Some(0.6)).expect("jc");
        let rw_kinds = rw.interaction(0).kinds();
        let jc_kinds = jc.interaction(0).kinds();
        assert_eq!(rw_kinds.len(), jc_kinds.len());
        for (left, right) in rw_kinds.iter().zip(jc_kinds) {
            assert_eq!(left.legs(), right.legs());
            assert_eq!(left.is_diagonal(), right.is_diagonal());
            assert!((left.weight() - right.weight()).abs() < 1.0e-14);
        }
        assert_eq!(rw.interaction(0).direction(), KernelDirection::Directed);
    }

    #[test]
    fn rw_crw_weights_follow_rho_dagger_rho_orientation() {
        let scale = 0.7;
        let ratio = 0.2;
        let model = ImpurityModel::rw_crw(
            mode(),
            scale,
            ratio,
            0.0,
            CouplingNormalization::FixedRw,
            Some(0.5),
        )
        .expect("model");
        for spin_a in [-1, 1] {
            for spin_b in [-1, 1] {
                let first = if spin_a == -1 { 1.0 } else { ratio };
                let second = if spin_b == 1 { 1.0 } else { ratio };
                let expected = scale * first * second;
                let legs = [spin_a, -spin_a, spin_b, -spin_b];
                let actual = kind_weight(&model, legs).expect("off-diagonal kind");
                assert!((actual - expected).abs() < 1.0e-14);
            }
        }
    }

    #[test]
    fn fixed_total_equal_rw_crw_matches_rotated_rabi_weights() {
        let scale = 0.8;
        let tunnelling = 0.15;
        let constant = Some(0.6);
        let rw = ImpurityModel::rw_crw(
            mode(),
            scale,
            1.0,
            tunnelling,
            CouplingNormalization::FixedTotal,
            constant,
        )
        .expect("rw");
        let rabi =
            ImpurityModel::rotated_impurity(mode(), scale, tunnelling, constant).expect("rabi");
        for kind in rw.interaction(0).kinds() {
            let matching = kind_weight(&rabi, *kind.legs()).expect("matching kind");
            assert!((matching - kind.weight()).abs() < 1.0e-14);
        }
    }

    #[test]
    fn rotated_model_reports_physical_axis_map() {
        let model = ImpurityModel::rotated_impurity(mode(), 0.5, 0.2, Some(0.4)).expect("model");
        assert_eq!(model.basis_transform(), BasisTransform::rotated_rabi());
    }

    #[test]
    fn shifted_energy_is_corrected_by_catalog_constant() {
        let model = ImpurityModel::xxz(mode(), 0.2, 0.1, 0.0, Some(0.7)).expect("model");
        assert!((model.corrected_spin_coupling_energy(10.0, 5.0) + 1.3).abs() < 1e-14);
    }

    #[test]
    fn multi_channel_unspecified_pair_flip_gauge_is_rejected() {
        let (kinds, _) =
            build_catalog(0.0, 0.0, Some(1.0), &[("pair", [-1, 1, -1, 1], 0.2)]).expect("catalog");
        let pair = InteractionChannel::new("pair", mode(), KernelDirection::Symmetric, kinds)
            .expect("channel");
        let (diag_kinds, _) = build_catalog(0.0, 0.0, Some(1.0), &[]).expect("catalog");
        let diagonal =
            InteractionChannel::new("diagonal", mode(), KernelDirection::Symmetric, diag_kinds)
                .expect("channel");
        assert!(ImpurityModel::from_interactions("bad", vec![pair, diagonal]).is_err());
    }
}
