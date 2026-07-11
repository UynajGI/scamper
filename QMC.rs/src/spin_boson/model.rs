//! Sign-free impurity model catalogs for wormhole QMC.

use std::collections::HashMap;

use super::bath::{Bath, KernelDirection};
use super::error::SpinBosonError;
use super::scattering::{ScatteringPolicy, ScatteringTable};
use super::vertex::{Spin, VertexKind};

/// Supported impurity Hamiltonian families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinBosonModelKind {
    /// Rotating-wave coupling `a^dagger S_- + S_+ a`.
    JaynesCummings,
    /// U(1)-symmetric coordinate coupling with `lambda_x = lambda_y`.
    Xxz,
    /// Fully anisotropic coordinate coupling.
    Xyz,
    /// Original longitudinal spin-boson/Rabi model after a spin-axis rotation.
    RotatedSpinBoson,
    /// User-composed positive interaction channels.
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
}

impl InteractionChannel {
    /// Construct a channel and precompute its scattering table.
    pub fn new(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
    ) -> Result<Self, SpinBosonError> {
        Self::with_scattering_policy(name, bath, direction, kinds, ScatteringPolicy::LowBounce)
    }

    /// Construct a channel with an explicit local scattering policy.
    pub fn with_scattering_policy(
        name: impl Into<String>,
        bath: Bath,
        direction: KernelDirection,
        kinds: Vec<VertexKind>,
        scattering_policy: ScatteringPolicy,
    ) -> Result<Self, SpinBosonError> {
        if kinds.is_empty() {
            return Err(SpinBosonError::parameter(
                "vertex catalog",
                "an interaction channel needs at least one vertex kind",
            ));
        }
        let mut diagonal_lookup = HashMap::new();
        for (kind_id, kind) in kinds.iter().enumerate() {
            if kind.is_diagonal() {
                diagonal_lookup.insert((kind.spin(0), kind.spin(2)), kind_id);
            }
        }
        for spin_a in [-1, 1] {
            for spin_b in [-1, 1] {
                if !diagonal_lookup.contains_key(&(spin_a, spin_b)) {
                    return Err(SpinBosonError::parameter(
                        "vertex catalog",
                        format!("missing diagonal seed for ({spin_a},{spin_b})"),
                    ));
                }
            }
        }
        let scattering = ScatteringTable::build(&kinds, scattering_policy);
        Ok(Self {
            name: name.into(),
            bath,
            direction,
            kinds,
            diagonal_lookup,
            scattering,
        })
    }

    /// Channel name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Normalized bath shape.
    pub fn bath(&self) -> &Bath {
        &self.bath
    }

    /// Directed or symmetric kernel.
    pub fn direction(&self) -> KernelDirection {
        self.direction
    }

    /// Local vertex kinds.
    pub fn kinds(&self) -> &[VertexKind] {
        &self.kinds
    }

    /// One local vertex kind.
    pub fn kind(&self, kind: usize) -> &VertexKind {
        &self.kinds[kind]
    }

    /// Diagonal kind matching the two endpoint worldline spins.
    pub fn diagonal_kind(&self, spin_a: Spin, spin_b: Spin) -> usize {
        self.diagonal_lookup[&(spin_a, spin_b)]
    }

    /// Local scattering table.
    pub fn scattering(&self) -> &ScatteringTable {
        &self.scattering
    }
}

/// Complete single-impurity model sampled by the generic engine.
#[derive(Debug, Clone, PartialEq)]
pub struct SpinBosonModel {
    kind: SpinBosonModelKind,
    name: String,
    interactions: Vec<InteractionChannel>,
}

impl SpinBosonModel {
    /// Compose a custom sign-free impurity model from positive interaction channels.
    ///
    /// This is the extension point for additional retarded impurity models: the
    /// generic update engine only depends on the channel catalogs and does not
    /// need model-specific branches.
    pub fn from_interactions(
        name: impl Into<String>,
        interactions: Vec<InteractionChannel>,
    ) -> Result<Self, SpinBosonError> {
        if interactions.is_empty() {
            return Err(SpinBosonError::parameter(
                "interactions",
                "a model requires at least one interaction channel",
            ));
        }
        Ok(Self {
            kind: SpinBosonModelKind::Custom,
            name: name.into(),
            interactions,
        })
    }

    /// Jaynes-Cummings model with effective retarded weight
    /// `lambda = integral J(omega)/(pi omega) d omega`.
    pub fn jaynes_cummings(
        bath: Bath,
        lambda: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, SpinBosonError> {
        validate_nonnegative("lambda", lambda)?;
        let offdiagonal = if lambda > 0.0 {
            vec![("Splus_A_Sminus_B", [-1, 1, 1, -1], lambda)]
        } else {
            Vec::new()
        };
        let kinds = build_catalog(0.0, h_z, constant, &offdiagonal)?;
        let interaction = InteractionChannel::new("jc", bath, KernelDirection::Directed, kinds)?;
        Ok(Self {
            kind: SpinBosonModelKind::JaynesCummings,
            name: "JaynesCummings".into(),
            interactions: vec![interaction],
        })
    }

    /// U(1)-symmetric XXZ spin-boson model.
    ///
    /// `lambda_xy` and `lambda_z` are the normalized retarded couplings.  For
    /// Weber's power law they are `2 alpha_l omega_c / s`; for one coordinate
    /// mode they are `g_l^2 / omega_0`.
    pub fn xxz(
        bath: Bath,
        lambda_xy: f64,
        lambda_z: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, SpinBosonError> {
        validate_nonnegative("lambda_xy", lambda_xy)?;
        validate_nonnegative("lambda_z", lambda_z)?;
        let mut offdiagonal = Vec::new();
        let exchange = 0.5 * lambda_xy;
        if exchange > 0.0 {
            offdiagonal.push(("Sminus_A_Splus_B", [1, -1, -1, 1], exchange));
            offdiagonal.push(("Splus_A_Sminus_B", [-1, 1, 1, -1], exchange));
        }
        let kinds = build_catalog(lambda_z, h_z, constant, &offdiagonal)?;
        let interaction = InteractionChannel::new("xxz", bath, KernelDirection::Symmetric, kinds)?;
        Ok(Self {
            kind: SpinBosonModelKind::Xxz,
            name: "XxzSpinBoson".into(),
            interactions: vec![interaction],
        })
    }

    /// Fully anisotropic XYZ coordinate-coupled spin-boson model.
    ///
    /// The pair-flip coefficient is sampled with its absolute value.  If
    /// `lambda_x < lambda_y`, a global `z`-axis phase rotation exchanges the
    /// sign; closed partition-function configurations contain pair vertices in
    /// parity-compatible combinations, so the sign-free catalog is unchanged.
    pub fn xyz(
        bath: Bath,
        lambda_x: f64,
        lambda_y: f64,
        lambda_z: f64,
        h_z: f64,
        constant: Option<f64>,
    ) -> Result<Self, SpinBosonError> {
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
        let kinds = build_catalog(lambda_z, h_z, constant, &offdiagonal)?;
        let interaction = InteractionChannel::new("xyz", bath, KernelDirection::Symmetric, kinds)?;
        Ok(Self {
            kind: SpinBosonModelKind::Xyz,
            name: "XyzSpinBoson".into(),
            interactions: vec![interaction],
        })
    }

    /// Original spin-boson/Rabi model in the rotated basis where the bath
    /// couples to `S_x` and the tunnelling field becomes a diagonal `h_z`.
    pub fn rotated_spin_boson(
        bath: Bath,
        lambda: f64,
        tunnelling: f64,
        constant: Option<f64>,
    ) -> Result<Self, SpinBosonError> {
        let mut model = Self::xyz(bath, lambda, 0.0, 0.0, tunnelling, constant)?;
        model.kind = SpinBosonModelKind::RotatedSpinBoson;
        model.name = "RotatedSpinBoson".into();
        model.interactions[0].name = "rotated_spin_boson".into();
        Ok(model)
    }

    /// Model kind.
    pub fn kind(&self) -> SpinBosonModelKind {
        self.kind
    }

    /// Model name used by Carlo.rs metadata.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Independently sampled interaction channels.
    pub fn interactions(&self) -> &[InteractionChannel] {
        &self.interactions
    }

    /// One interaction channel.
    pub fn interaction(&self, interaction: usize) -> &InteractionChannel {
        &self.interactions[interaction]
    }

    /// Number of interaction channels eligible for diagonal insertion.
    pub fn interaction_count(&self) -> usize {
        self.interactions.len()
    }
}

type OffDiagonalSpec<'a> = (&'a str, [Spin; 4], f64);

fn build_catalog(
    lambda_z: f64,
    h_z: f64,
    constant: Option<f64>,
    offdiagonal: &[OffDiagonalSpec<'_>],
) -> Result<Vec<VertexKind>, SpinBosonError> {
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
        return Err(SpinBosonError::parameter("C", "must be finite"));
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
    Ok(kinds)
}

fn validate_nonnegative(field: &str, value: f64) -> Result<(), SpinBosonError> {
    if !value.is_finite() || value < 0.0 {
        return Err(SpinBosonError::parameter(
            field,
            format!("must be finite and non-negative, got {value}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spin_boson::bath::SingleModeBath;

    fn mode() -> Bath {
        Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"))
    }

    #[test]
    fn custom_model_accepts_positive_channels() {
        let kinds = build_catalog(0.0, 0.0, None, &[]).expect("catalog");
        let channel =
            InteractionChannel::new("identity", mode(), KernelDirection::Symmetric, kinds)
                .expect("channel");
        let model =
            SpinBosonModel::from_interactions("custom", vec![channel]).expect("custom model");
        assert_eq!(model.kind(), SpinBosonModelKind::Custom);
    }

    #[test]
    fn jc_has_one_flip_kind() {
        let model = SpinBosonModel::jaynes_cummings(mode(), 0.4, 0.2, None).expect("model");
        let offdiag = model
            .interaction(0)
            .kinds()
            .iter()
            .filter(|kind| !kind.is_diagonal());
        assert_eq!(offdiag.count(), 1);
    }

    #[test]
    fn xxz_has_two_exchange_kinds() {
        let model = SpinBosonModel::xxz(mode(), 0.4, 0.1, 0.0, None).expect("model");
        let offdiag = model
            .interaction(0)
            .kinds()
            .iter()
            .filter(|kind| !kind.is_diagonal());
        assert_eq!(offdiag.count(), 2);
    }

    #[test]
    fn xyz_has_pair_flips() {
        let model = SpinBosonModel::xyz(mode(), 0.5, 0.1, 0.2, 0.0, None).expect("model");
        let names: Vec<_> = model
            .interaction(0)
            .kinds()
            .iter()
            .map(VertexKind::name)
            .collect();
        assert!(names.contains(&"Splus_A_Splus_B"));
        assert!(names.contains(&"Sminus_A_Sminus_B"));
    }
}
