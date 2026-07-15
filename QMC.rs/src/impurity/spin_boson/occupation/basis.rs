//! Finite spin-boson occupation basis.

use crate::impurity::ImpurityError;

/// Spin label in the sampled sigma-z basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinState {
    Down,
    Up,
}

impl SpinState {
    #[inline]
    pub const fn sigma_z(self) -> f64 {
        match self {
            Self::Down => -1.0,
            Self::Up => 1.0,
        }
    }

    #[inline]
    pub const fn index(self) -> usize {
        match self {
            Self::Down => 0,
            Self::Up => 1,
        }
    }

    #[inline]
    pub const fn from_index(index: usize) -> Self {
        if index == 0 {
            Self::Down
        } else {
            Self::Up
        }
    }
}

/// Tensor-product basis of one spin and a finite number of truncated bosonic modes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OccupationBasis {
    cutoffs: Vec<usize>,
    strides: Vec<usize>,
    boson_dimension: usize,
    dimension: usize,
}

impl OccupationBasis {
    /// `cutoffs[m]` is the number of retained occupation states `0..cutoff`.
    pub fn new(cutoffs: Vec<usize>) -> Result<Self, ImpurityError> {
        if cutoffs.is_empty() {
            return Err(ImpurityError::parameter(
                "cutoffs",
                "at least one bosonic mode is required",
            ));
        }
        if cutoffs.contains(&0) {
            return Err(ImpurityError::parameter(
                "cutoffs",
                "every occupation cutoff must be positive",
            ));
        }
        let mut strides = Vec::with_capacity(cutoffs.len());
        let mut dimension = 1usize;
        for &cutoff in &cutoffs {
            strides.push(dimension);
            dimension = dimension.checked_mul(cutoff).ok_or_else(|| {
                ImpurityError::parameter("cutoffs", "occupation basis dimension overflows usize")
            })?;
        }
        let total = dimension.checked_mul(2).ok_or_else(|| {
            ImpurityError::parameter("cutoffs", "spin-boson basis dimension overflows usize")
        })?;
        Ok(Self {
            cutoffs,
            strides,
            boson_dimension: dimension,
            dimension: total,
        })
    }

    #[inline]
    pub const fn modes(&self) -> usize {
        self.cutoffs.len()
    }

    #[inline]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[inline]
    pub const fn boson_dimension(&self) -> usize {
        self.boson_dimension
    }

    #[inline]
    pub fn cutoff(&self, mode: usize) -> usize {
        self.cutoffs[mode]
    }

    #[inline]
    pub fn cutoffs(&self) -> &[usize] {
        &self.cutoffs
    }

    pub fn encode(&self, spin: SpinState, occupations: &[usize]) -> Result<usize, ImpurityError> {
        if occupations.len() != self.modes() {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "expected {} occupations, got {}",
                self.modes(),
                occupations.len()
            )));
        }
        let mut boson = 0usize;
        for (mode, &occupation) in occupations.iter().enumerate() {
            if occupation >= self.cutoffs[mode] {
                return Err(ImpurityError::InvalidConfiguration(format!(
                    "occupation {occupation} exceeds mode {mode} cutoff {}",
                    self.cutoffs[mode]
                )));
            }
            boson += occupation * self.strides[mode];
        }
        Ok(2 * boson + spin.index())
    }

    pub fn spin(&self, state: usize) -> SpinState {
        SpinState::from_index(state & 1)
    }

    pub fn occupation(&self, state: usize, mode: usize) -> usize {
        let boson = state / 2;
        (boson / self.strides[mode]) % self.cutoffs[mode]
    }

    pub fn decode_into(&self, state: usize, occupations: &mut Vec<usize>) -> SpinState {
        occupations.clear();
        occupations.reserve(self.modes());
        for mode in 0..self.modes() {
            occupations.push(self.occupation(state, mode));
        }
        self.spin(state)
    }

    pub(crate) fn shifted_state(&self, state: usize, mode: usize, delta: isize) -> Option<usize> {
        let occupation = self.occupation(state, mode) as isize;
        let shifted = occupation + delta;
        if shifted < 0 || shifted >= self.cutoffs[mode] as isize {
            return None;
        }
        let boson = state / 2;
        let next_boson = if delta > 0 {
            boson + self.strides[mode]
        } else {
            boson - self.strides[mode]
        };
        Some(2 * next_boson + (state & 1))
    }

    #[inline]
    pub(crate) const fn flipped_spin(&self, state: usize) -> usize {
        state ^ 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimension_is_two_times_product_of_cutoffs() {
        let basis = OccupationBasis::new(vec![3, 4]).expect("basis");
        assert_eq!(basis.modes(), 2);
        assert_eq!(basis.boson_dimension(), 12);
        assert_eq!(basis.dimension(), 24);
        assert_eq!(basis.cutoff(0), 3);
        assert_eq!(basis.cutoff(1), 4);
        assert_eq!(basis.cutoffs(), &[3, 4]);
    }

    #[test]
    fn single_mode_basis_has_two_spin_sectors() {
        let basis = OccupationBasis::new(vec![5]).expect("basis");
        assert_eq!(basis.boson_dimension(), 5);
        assert_eq!(basis.dimension(), 10);
    }

    #[test]
    fn empty_cutoffs_are_rejected() {
        assert!(OccupationBasis::new(vec![]).is_err());
    }

    #[test]
    fn zero_cutoff_is_rejected() {
        assert!(OccupationBasis::new(vec![3, 0]).is_err());
    }

    #[test]
    fn encode_decode_round_trips() {
        let basis = OccupationBasis::new(vec![3, 4]).expect("basis");
        for spin in [SpinState::Down, SpinState::Up] {
            for n0 in 0..basis.cutoff(0) {
                for n1 in 0..basis.cutoff(1) {
                    let state = basis.encode(spin, &[n0, n1]).expect("encode");
                    assert_eq!(basis.spin(state), spin);
                    assert_eq!(basis.occupation(state, 0), n0);
                    assert_eq!(basis.occupation(state, 1), n1);
                }
            }
        }
    }

    #[test]
    fn decode_into_returns_full_occupation_vector() {
        let basis = OccupationBasis::new(vec![2, 3, 4]).expect("basis");
        let state = basis.encode(SpinState::Up, &[1, 2, 3]).expect("encode");
        let mut occupations = Vec::new();
        let spin = basis.decode_into(state, &mut occupations);
        assert_eq!(spin, SpinState::Up);
        assert_eq!(occupations, &[1, 2, 3]);
    }

    #[test]
    fn encode_rejects_out_of_range_occupation() {
        let basis = OccupationBasis::new(vec![3]).expect("basis");
        assert!(basis.encode(SpinState::Up, &[3]).is_err());
    }

    #[test]
    fn encode_rejects_wrong_occupation_count() {
        let basis = OccupationBasis::new(vec![3, 4]).expect("basis");
        assert!(basis.encode(SpinState::Down, &[1]).is_err());
    }

    #[test]
    fn shifted_state_respects_cutoffs() {
        let basis = OccupationBasis::new(vec![3]).expect("basis");
        // |down, n=1> = state 2. Raising gives |down, n=2> = state 4.
        assert_eq!(basis.shifted_state(2, 0, 1), Some(4));
        // Lowering gives |down, n=0> = state 0.
        assert_eq!(basis.shifted_state(2, 0, -1), Some(0));
        // Raising past the cutoff returns None.
        let top = basis.encode(SpinState::Up, &[2]).expect("encode");
        assert_eq!(basis.shifted_state(top, 0, 1), None);
        // Lowering below zero returns None.
        assert_eq!(basis.shifted_state(0, 0, -1), None);
    }

    #[test]
    fn shifted_state_preserves_spin() {
        let basis = OccupationBasis::new(vec![4]).expect("basis");
        let up = basis.encode(SpinState::Up, &[1]).expect("encode");
        let raised = basis.shifted_state(up, 0, 1).expect("raise");
        assert_eq!(basis.spin(raised), SpinState::Up);
        assert_eq!(basis.occupation(raised, 0), 2);
    }

    #[test]
    fn flipped_spin_swaps_up_and_down() {
        let basis = OccupationBasis::new(vec![3]).expect("basis");
        let down0 = basis.encode(SpinState::Down, &[0]).expect("encode");
        let up0 = basis.encode(SpinState::Up, &[0]).expect("encode");
        assert_eq!(basis.flipped_spin(down0), up0);
        assert_eq!(basis.flipped_spin(up0), down0);
        assert_eq!(basis.flipped_spin(basis.flipped_spin(down0)), down0);
    }

    #[test]
    fn spin_state_arithmetic_is_consistent() {
        assert_eq!(SpinState::Down.sigma_z(), -1.0);
        assert_eq!(SpinState::Up.sigma_z(), 1.0);
        assert_eq!(SpinState::Down.index(), 0);
        assert_eq!(SpinState::Up.index(), 1);
        assert_eq!(SpinState::from_index(0), SpinState::Down);
        assert_eq!(SpinState::from_index(1), SpinState::Up);
    }
}
