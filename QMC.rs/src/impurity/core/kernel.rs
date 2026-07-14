//! Retarded-kernel direction and sign-free composition metadata.

use crate::impurity::core::operators::BasisTransform;
use crate::impurity::ImpurityError;

/// Whether a retarded operator uses the directed propagator `D` or the
/// symmetrized propagator `D_+`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelDirection {
    /// Keep the orientation of imaginary-time propagation. Used by JC and
    /// rotating/counter-rotating channels.
    Directed,
    /// Sample the two orientations with equal probability. Used by Hermitian
    /// coordinate couplings such as XXZ and XYZ.
    Symmetric,
}

/// Gauge convention used to make an XYZ pair-flip matrix element non-negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairFlipGauge {
    /// The channel has no pair-flip vertices.
    NotPresent,
    /// The physical pair-flip coefficient is non-negative in the sampled basis.
    Positive,
    /// A global `pi/2` rotation around `Z` was used to absorb a negative sign.
    Negative,
    /// A custom channel contains pair flips but did not declare their gauge.
    Unspecified,
}

/// Sign-free assumptions attached to one interaction channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SignFreeMetadata {
    /// Basis in which all local matrix elements are sampled as non-negative.
    pub basis: BasisTransform,
    /// Pair-flip gauge, if pair flips are present.
    pub pair_flip_gauge: PairFlipGauge,
}

impl SignFreeMetadata {
    pub const fn new(basis: BasisTransform, pair_flip_gauge: PairFlipGauge) -> Self {
        Self {
            basis,
            pair_flip_gauge,
        }
    }

    pub const fn diagonal_or_exchange(basis: BasisTransform) -> Self {
        Self::new(basis, PairFlipGauge::NotPresent)
    }
}

impl Default for SignFreeMetadata {
    fn default() -> Self {
        Self::diagonal_or_exchange(BasisTransform::identity())
    }
}

/// Validated model-wide sign-free composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SignFreeReport {
    pub channel_count: usize,
    pub basis: BasisTransform,
    pub pair_flip_gauge: PairFlipGauge,
}

/// Validate that independently sampled channels use one common non-negative
/// basis and one compatible pair-flip gauge.
///
/// A single custom pair-flip channel may remain `Unspecified` for backward
/// compatibility.  Combining such a channel with another channel is rejected:
/// the relative sign between channels is then physical and cannot be inferred
/// from positive local weights alone.
pub fn validate_sign_free_channels(
    channels: &[(&str, SignFreeMetadata)],
) -> Result<SignFreeReport, ImpurityError> {
    let Some((_, first)) = channels.first().copied() else {
        return Err(ImpurityError::parameter(
            "interactions",
            "a model requires at least one interaction channel",
        ));
    };

    let mut gauge = PairFlipGauge::NotPresent;
    for (name, metadata) in channels {
        if metadata.basis != first.basis {
            return Err(ImpurityError::parameter(
                "interactions",
                format!(
                    "channel `{name}` uses a different sampled-to-physical basis; \
                     multiple channels must share one global basis rotation"
                ),
            ));
        }
        match metadata.pair_flip_gauge {
            PairFlipGauge::NotPresent => {}
            PairFlipGauge::Unspecified if channels.len() > 1 => {
                return Err(ImpurityError::parameter(
                    "interactions",
                    format!(
                        "channel `{name}` contains pair flips but does not declare their sign \
                         gauge; use InteractionChannel::with_metadata before composing channels"
                    ),
                ));
            }
            PairFlipGauge::Unspecified => gauge = PairFlipGauge::Unspecified,
            declared => match gauge {
                PairFlipGauge::NotPresent => gauge = declared,
                PairFlipGauge::Unspecified => gauge = declared,
                previous if previous == declared => {}
                previous => {
                    return Err(ImpurityError::parameter(
                        "interactions",
                        format!(
                            "incompatible pair-flip gauges across channels: {previous:?} and \
                             {declared:?}; no single global spin rotation makes all weights positive"
                        ),
                    ));
                }
            },
        }
    }

    Ok(SignFreeReport {
        channel_count: channels.len(),
        basis: first.basis,
        pair_flip_gauge: gauge,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incompatible_basis_is_rejected() {
        let result = validate_sign_free_channels(&[
            ("a", SignFreeMetadata::default()),
            (
                "b",
                SignFreeMetadata::diagonal_or_exchange(BasisTransform::rotated_rabi()),
            ),
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn conflicting_pair_flip_gauges_are_rejected() {
        let result = validate_sign_free_channels(&[
            (
                "positive",
                SignFreeMetadata::new(BasisTransform::identity(), PairFlipGauge::Positive),
            ),
            (
                "negative",
                SignFreeMetadata::new(BasisTransform::identity(), PairFlipGauge::Negative),
            ),
        ]);
        assert!(result.is_err());
    }
}
