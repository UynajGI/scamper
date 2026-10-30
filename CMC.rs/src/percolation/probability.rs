use core::fmt;

/// Error returned when constructing an occupation probability outside `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProbabilityError {
    value: f64,
}

impl ProbabilityError {
    #[inline]
    pub const fn value(self) -> f64 {
        self.value
    }
}

impl fmt::Display for ProbabilityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "probability must be finite within [0, 1], got {}",
            self.value
        )
    }
}

impl std::error::Error for ProbabilityError {}

/// A finite occupation probability in the closed interval `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Probability(f64);

impl Probability {
    pub fn new(value: f64) -> Result<Self, ProbabilityError> {
        if value.is_finite() && (0.0..=1.0).contains(&value) {
            Ok(Self(value))
        } else {
            Err(ProbabilityError { value })
        }
    }

    #[inline]
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// A validated uniform, borrowed heterogeneous, or owned heterogeneous field.
#[derive(Clone, Debug)]
pub enum ProbabilityField<'a> {
    Uniform(Probability),
    Borrowed(&'a [Probability]),
    Owned(Vec<Probability>),
}

impl ProbabilityField<'_> {
    /// Heterogeneous length, or `None` for a uniform field valid at any length.
    #[inline]
    pub fn len(&self) -> Option<usize> {
        match self {
            Self::Uniform(_) => None,
            Self::Borrowed(values) => Some(values.len()),
            Self::Owned(values) => Some(values.len()),
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == Some(0)
    }

    /// Probability at an index valid for the sampling domain.
    ///
    /// Heterogeneous fields panic for an invalid index; uniform fields are
    /// defined for every index and rely on the sampling domain for bounds.
    #[inline]
    pub fn at(&self, index: usize) -> Probability {
        match self {
            Self::Uniform(probability) => *probability,
            Self::Borrowed(values) => values[index],
            Self::Owned(values) => values[index],
        }
    }

    #[inline]
    pub const fn uniform(&self) -> Option<Probability> {
        match self {
            Self::Uniform(probability) => Some(*probability),
            Self::Borrowed(_) | Self::Owned(_) => None,
        }
    }
}

impl From<Probability> for ProbabilityField<'_> {
    fn from(probability: Probability) -> Self {
        Self::Uniform(probability)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probability(value: f64) -> Probability {
        Probability::new(value).expect("test probability")
    }

    #[test]
    fn rejects_invalid_values() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
            assert!(Probability::new(value).is_err(), "accepted {value}");
        }
        assert_eq!(Probability::new(0.0).unwrap().get(), 0.0);
        assert_eq!(Probability::new(1.0).unwrap().get(), 1.0);
    }

    #[test]
    fn fields_report_length_and_values() {
        let values = [probability(0.1), probability(0.9)];
        let uniform = ProbabilityField::Uniform(probability(0.5));
        let borrowed = ProbabilityField::Borrowed(&values);
        let owned = ProbabilityField::Owned(values.to_vec());

        assert_eq!(uniform.len(), None);
        assert_eq!(uniform.uniform(), Some(probability(0.5)));
        assert_eq!(borrowed.len(), Some(2));
        assert_eq!(owned.len(), Some(2));
        assert_eq!(borrowed.at(1), probability(0.9));
        assert_eq!(owned.at(0), probability(0.1));
        assert_eq!(borrowed.uniform(), None);
    }
}
