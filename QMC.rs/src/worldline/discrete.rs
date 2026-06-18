use super::{State, Worldline};

/// Discrete-time worldline: fixed number M of Trotter slices of width Δτ = β/M.
///
/// `states[i]` is the state at τ = i·Δτ. Kinks are detected as adjacent slices
/// with differing states.
pub struct DiscreteWorldline {
    beta: f64,
    delta_tau: f64,
    dim: u8,
    /// states[i] = state on [i·Δτ, (i+1)·Δτ), i ∈ [0, M). Periodicity: state_M = state_0.
    states: Box<[State]>,
}

impl DiscreteWorldline {
    pub fn new(beta: f64, dim: u8, m: usize, initial_state: State) -> Self {
        assert!(dim > 0);
        assert!(initial_state < dim);
        assert!(beta > 0.0);
        assert!(m > 0);
        let delta_tau = beta / m as f64;
        Self {
            beta,
            delta_tau,
            dim,
            states: vec![initial_state; m].into_boxed_slice(),
        }
    }

    /// Number of time slices M.
    pub fn m(&self) -> usize {
        self.states.len()
    }

    pub fn delta_tau(&self) -> f64 {
        self.delta_tau
    }
}

impl Worldline for DiscreteWorldline {
    fn beta(&self) -> f64 {
        self.beta
    }

    fn dim(&self) -> u8 {
        self.dim
    }

    fn num_kinks(&self) -> usize {
        let m = self.states.len();
        let mut count = 0;
        for i in 0..m {
            if self.states[i] != self.states[(i + 1) % m] {
                count += 1;
            }
        }
        count
    }

    fn state_at(&self, tau: f64) -> State {
        assert!(tau >= 0.0 && tau < self.beta);
        let i = (tau / self.delta_tau) as usize;
        self.states[i]
    }

    fn for_each_kink(&self, mut f: impl FnMut(f64, State, State)) {
        let m = self.states.len();
        for i in 0..m {
            let cur = self.states[i];
            let next = self.states[(i + 1) % m];
            if cur != next {
                let tau = ((i + 1) as f64) * self.delta_tau;
                if tau < self.beta {
                    f(tau, cur, next);
                }
                // tau = β wraps to tau = 0, skip since it's redundant with the kink at β
            }
        }
    }

    fn insert_kink(&mut self, _tau: f64, _to: State) {
        unimplemented!("DiscreteWorldline insert_kink not yet implemented")
    }

    fn remove_kink(&mut self, _idx: usize) {
        unimplemented!("DiscreteWorldline remove_kink not yet implemented")
    }

    fn diagonal(&self) -> f64 {
        let sum: u64 = self.states.iter().map(|&s| s as u64).sum();
        sum as f64 / self.states.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_worldline() {
        let wl = DiscreteWorldline::new(10.0, 3, 10, 1);
        assert_eq!(wl.num_kinks(), 0);
        assert_eq!(wl.state_at(1.5), 1);
        assert_eq!(wl.diagonal(), 1.0);
    }

    #[test]
    fn test_manual_kinks() {
        // Manually construct a worldline with kinks via mutation would need insert_kink.
        // For now, verify the structure.
        let wl = DiscreteWorldline::new(10.0, 4, 100, 2);
        assert_eq!(wl.beta(), 10.0);
        assert_eq!(wl.dim(), 4);
        assert_eq!(wl.m(), 100);
        assert_eq!(wl.delta_tau(), 0.1);
    }

    #[test]
    fn test_state_at_mapping() {
        let wl = DiscreteWorldline::new(5.0, 3, 5, 0);
        // M=5, Δτ=1.0, τ=0,1,2,3,4 → slice 0,1,2,3,4
        assert_eq!(wl.state_at(0.5), 0);
        assert_eq!(wl.state_at(1.5), 0);
        assert_eq!(wl.state_at(4.5), 0);
    }
}
