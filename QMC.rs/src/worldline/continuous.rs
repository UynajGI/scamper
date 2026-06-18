use super::{State, Worldline};

/// Continuous-time worldline: kinks stored directly as `(tau, from, to)`.
///
/// All operations are O(log N) via binary search over kinks sorted by τ.
pub struct ContinuousWorldline {
    beta: f64,
    dim: u8,
    /// State at τ = 0 (same as state at τ → β⁻ due to periodicity).
    initial_state: State,
    /// Kinks sorted by τ ∈ (0, β). Each entry: at time τ, state flips from → to.
    kinks: Vec<(f64, State, State)>,
}

impl ContinuousWorldline {
    pub fn new(beta: f64, dim: u8, initial_state: State) -> Self {
        assert!(dim > 0);
        assert!(initial_state < dim);
        assert!(beta > 0.0);
        Self {
            beta,
            dim,
            initial_state,
            kinks: Vec::new(),
        }
    }

    /// Number of kinks strictly before `tau` (left-continuous convention:
    /// `state_at(τ)` returns the state just *before* the kink at τ).
    fn kink_before(&self, tau: f64) -> usize {
        self.kinks.partition_point(|k| k.0 < tau)
    }
}

impl Worldline for ContinuousWorldline {
    fn beta(&self) -> f64 {
        self.beta
    }

    fn dim(&self) -> u8 {
        self.dim
    }

    fn num_kinks(&self) -> usize {
        self.kinks.len()
    }

    fn state_at(&self, tau: f64) -> State {
        assert!(tau >= 0.0 && tau < self.beta);
        match self.kink_before(tau) {
            0 => self.initial_state,
            i => self.kinks[i - 1].2, // to_state of the previous kink
        }
    }

    fn for_each_kink(&self, mut f: impl FnMut(f64, State, State)) {
        for &(tau, from, to) in &self.kinks {
            f(tau, from, to);
        }
    }

    fn insert_kink(&mut self, tau: f64, to: State) {
        assert!(tau > 0.0 && tau < self.beta);
        assert!(to < self.dim);
        let from = self.state_at(tau);
        if from == to {
            return;
        }
        let idx = self.kink_before(tau);
        // Update the from_state of the next kink (if any) to reflect the flip
        if idx < self.kinks.len() {
            self.kinks[idx].1 = to;
        }
        self.kinks.insert(idx, (tau, from, to));
    }

    fn remove_kink(&mut self, idx: usize) {
        let (_tau, from, _to) = self.kinks.remove(idx);
        // The kink after the removed one now continues from `from`
        // (because the removed kink no longer flips state to `to`)
        if idx < self.kinks.len() {
            self.kinks[idx].1 = from;
        }
    }

    fn diagonal(&self) -> f64 {
        if self.kinks.is_empty() {
            return self.initial_state as f64;
        }
        let mut total = 0.0;
        let mut prev_tau = 0.0f64;
        let mut prev_state = self.initial_state as f64;
        for &(tau, _from, to) in &self.kinks {
            total += prev_state * (tau - prev_tau);
            prev_tau = tau;
            prev_state = to as f64;
        }
        // Final segment: from last kink to β
        total += prev_state * (self.beta - prev_tau);
        total / self.beta
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_worldline() {
        let wl = ContinuousWorldline::new(10.0, 3, 0);
        assert_eq!(wl.state_at(0.0), 0);
        assert_eq!(wl.state_at(5.0), 0);
        assert_eq!(wl.state_at(9.9), 0);
        assert_eq!(wl.num_kinks(), 0);
        assert_eq!(wl.diagonal(), 0.0);
    }

    #[test]
    fn test_insert_and_state_at() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 0);
        wl.insert_kink(3.0, 1);
        wl.insert_kink(7.0, 2);
        assert_eq!(wl.num_kinks(), 2);
        assert_eq!(wl.state_at(1.0), 0);
        assert_eq!(wl.state_at(3.0), 0); // left-continuous: before kink
        assert_eq!(wl.state_at(3.1), 1);
        assert_eq!(wl.state_at(5.0), 1);
        assert_eq!(wl.state_at(7.0), 1); // left-continuous: before second kink
        assert_eq!(wl.state_at(8.0), 2);
    }

    #[test]
    fn test_noop_same_state() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 0);
        wl.insert_kink(3.0, 1);
        wl.insert_kink(5.0, 1); // no-op: already state 1
        assert_eq!(wl.num_kinks(), 1);
    }

    #[test]
    fn test_remove_kink() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 0);
        wl.insert_kink(3.0, 1);
        wl.insert_kink(7.0, 2);
        wl.remove_kink(0); // remove the first kink
        assert_eq!(wl.num_kinks(), 1);
        assert_eq!(wl.state_at(1.0), 0);
        assert_eq!(wl.state_at(5.0), 0); // kink at 3 removed, back to 0
        assert_eq!(wl.state_at(8.0), 2);
    }

    #[test]
    fn test_for_each_kink() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 0);
        wl.insert_kink(3.0, 1);
        wl.insert_kink(7.0, 2);
        let mut kinks = Vec::new();
        wl.for_each_kink(|tau, from, to| kinks.push((tau, from, to)));
        assert_eq!(kinks, vec![(3.0, 0, 1), (7.0, 1, 2)]);
    }

    #[test]
    fn test_diagonal() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 1);
        // state 1 from 0 to 3, state 2 from 3 to 7, state 0 from 7 to 10
        wl.insert_kink(3.0, 2);
        wl.insert_kink(7.0, 0);
        let expected = (1.0 * 3.0 + 2.0 * 4.0 + 0.0 * 3.0) / 10.0; // 1.1
        assert!((wl.diagonal() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_insert_middle() {
        // Insert between existing kinks
        let mut wl = ContinuousWorldline::new(10.0, 4, 0);
        wl.insert_kink(2.0, 1);
        wl.insert_kink(8.0, 3);
        wl.insert_kink(5.0, 2);
        assert_eq!(wl.num_kinks(), 3);
        assert_eq!(wl.state_at(3.0), 1);
        assert_eq!(wl.state_at(6.0), 2);
        assert_eq!(wl.state_at(9.0), 3);
        assert_eq!(wl.state_at(0.5), 0);
    }
}
