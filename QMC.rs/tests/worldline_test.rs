use qmc_rs::worldline::{ContinuousWorldline, DiscreteWorldline, Worldline};

// ── helpers ──────────────────────────────────────────────────────────

/// Collect kinks via for_each_kink for assertion.
fn collect_kinks(wl: &impl Worldline) -> Vec<(f64, u8, u8)> {
    let mut v = Vec::with_capacity(wl.num_kinks());
    wl.for_each_kink(|tau, from, to| v.push((tau, from, to)));
    v
}

/// Verify internal consistency: for_each_kink count matches num_kinks, kinks
/// are sorted by τ, and state_at transitions match kink data.
fn assert_consistent<W: Worldline>(wl: &W) {
    let mut count = 0;
    wl.for_each_kink(|_, _, _| count += 1);
    assert_eq!(count, wl.num_kinks());

    let kinks = collect_kinks(wl);
    // Sorted by τ
    for w in kinks.windows(2) {
        assert!(w[0].0 < w[1].0, "kinks not sorted: {:?} >= {:?}", w[0], w[1]);
    }
    // Transitions match: from must equal state just before τ
    for &(tau, from, _to) in &kinks {
        if tau > 1e-12 {
            assert_eq!(wl.state_at(tau - 1e-12), from,
                "at τ={}, expected state before = {}", tau, from);
        }
    }
}

// ── ContinuousWorldline ──────────────────────────────────────────────

mod continuous {
    use super::*;

    #[test]
    fn spin_half_dim2() {
        let mut wl = ContinuousWorldline::new(10.0, 2, 0);
        assert_eq!(wl.dim(), 2);
        wl.insert_kink(2.0, 1);
        wl.insert_kink(5.0, 0);
        wl.insert_kink(8.0, 1);
        // 3 kinks → not periodic (odd kink count), but should be consistent
        assert_eq!(wl.num_kinks(), 3);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(0.0), 0);
        assert_eq!(wl.state_at(3.0), 1);
        assert_eq!(wl.state_at(6.0), 0);
        assert_eq!(wl.state_at(9.0), 1);
        let kinks = collect_kinks(&wl);
        assert_eq!(kinks, vec![(2.0, 0, 1), (5.0, 1, 0), (8.0, 0, 1)]);
    }

    #[test]
    fn spin_one_dim3() {
        let mut wl = ContinuousWorldline::new(12.0, 3, 1); // m=0
        assert_eq!(wl.dim(), 3);
        wl.insert_kink(3.0, 2);  // m=0 → m=+1
        wl.insert_kink(6.0, 0);  // m=+1 → m=-1
        wl.insert_kink(9.0, 1);  // m=-1 → m=0
        assert_eq!(wl.num_kinks(), 3);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(1.0), 1);
        assert_eq!(wl.state_at(4.0), 2);
        assert_eq!(wl.state_at(7.0), 0);
        assert_eq!(wl.state_at(10.0), 1);
    }

    #[test]
    fn spin_three_half_dim4() {
        let mut wl = ContinuousWorldline::new(8.0, 4, 0); // m=-3/2
        for m in 1u8..4 {
            wl.insert_kink(m as f64 * 2.0, m);
        }
        assert_eq!(wl.num_kinks(), 3);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(1.0), 0);
        assert_eq!(wl.state_at(3.0), 1);
        assert_eq!(wl.state_at(5.0), 2);
        assert_eq!(wl.state_at(7.0), 3);
    }

    #[test]
    fn spin_five_half_dim6() {
        let mut wl = ContinuousWorldline::new(10.0, 6, 3); // m=0
        wl.insert_kink(2.0, 5);  // → m=+5/2
        wl.insert_kink(4.0, 0);  // → m=-5/2
        wl.insert_kink(7.0, 3);  // → m=0
        assert_eq!(wl.num_kinks(), 3);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(1.0), 3);
        assert_eq!(wl.state_at(3.0), 5);
        assert_eq!(wl.state_at(5.0), 0);
        assert_eq!(wl.state_at(8.0), 3);
    }

    #[test]
    fn bosonic_dim10() {
        // n_max = 9 occupation states
        let mut wl = ContinuousWorldline::new(20.0, 10, 0);
        wl.insert_kink(2.0, 3);
        wl.insert_kink(5.0, 7);
        wl.insert_kink(10.0, 1);
        wl.insert_kink(15.0, 9);
        wl.insert_kink(18.0, 4);
        assert_eq!(wl.num_kinks(), 5);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(0.5), 0);
        assert_eq!(wl.state_at(3.0), 3);
        assert_eq!(wl.state_at(7.0), 7);
        assert_eq!(wl.state_at(12.0), 1);
        assert_eq!(wl.state_at(16.0), 9);
        assert_eq!(wl.state_at(19.0), 4);
    }

    #[test]
    fn insert_many_kinks_stress() {
        let beta = 101.0;
        let mut wl = ContinuousWorldline::new(beta, 5, 0);
        for i in 1..=100u8 {
            wl.insert_kink(i as f64, i % 5);
        }
        assert_eq!(wl.num_kinks(), 100);
        assert_consistent(&wl);
        // state after 50th kink: 50%5=0
        assert_eq!(wl.state_at(50.5), 0);
    }

    #[test]
    fn remove_all_kinks() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 1);
        wl.insert_kink(2.0, 0);
        wl.insert_kink(5.0, 2);
        wl.insert_kink(8.0, 1);
        for i in (0..3).rev() {
            wl.remove_kink(i);
        }
        assert_eq!(wl.num_kinks(), 0);
        assert_eq!(wl.state_at(3.0), 1);
        assert_eq!(wl.state_at(9.0), 1);
        assert_consistent(&wl);
    }

    #[test]
    fn remove_from_middle() {
        let mut wl = ContinuousWorldline::new(10.0, 4, 0);
        wl.insert_kink(2.0, 1); // kink 0
        wl.insert_kink(5.0, 2); // kink 1  ← remove this
        wl.insert_kink(8.0, 3); // kink 2
        wl.remove_kink(1);
        assert_eq!(wl.num_kinks(), 2);
        // After removal: state goes 0→1 at τ=2, then 1→3 at τ=8 (from was updated)
        assert_eq!(wl.state_at(3.0), 1);
        assert_eq!(wl.state_at(6.0), 1); // was 2, now 1
        assert_eq!(wl.state_at(9.0), 3);
        assert_consistent(&wl);
    }

    #[test]
    fn diagonal_spin_half() {
        let mut wl = ContinuousWorldline::new(10.0, 2, 0);
        // 0→1 at τ=4 → state 0 for [0,4), state 1 for [4,10)
        wl.insert_kink(4.0, 1);
        let expected = (0.0 * 4.0 + 1.0 * 6.0) / 10.0; // 0.6
        assert!((wl.diagonal() - expected).abs() < 1e-10);
    }

    #[test]
    fn noop_insert_same_state() {
        let mut wl = ContinuousWorldline::new(10.0, 3, 1);
        wl.insert_kink(3.0, 2);
        wl.insert_kink(5.0, 2); // no-op
        wl.insert_kink(5.5, 2); // no-op
        assert_eq!(wl.num_kinks(), 1);
        assert_consistent(&wl);
    }
}

// ── DiscreteWorldline ────────────────────────────────────────────────

mod discrete {
    use super::*;

    #[test]
    fn spin_half_dim2() {
        let wl = DiscreteWorldline::new(10.0, 2, 10, 0);
        assert_eq!(wl.beta(), 10.0);
        assert_eq!(wl.dim(), 2);
        assert_eq!(wl.m(), 10);
        assert_eq!(wl.delta_tau(), 1.0);
        assert_eq!(wl.num_kinks(), 0);
        assert_consistent(&wl);
        assert_eq!(wl.state_at(1.5), 0);
        assert_eq!(wl.state_at(5.5), 0);
    }

    #[test]
    fn spin_one_dim3() {
        let wl = DiscreteWorldline::new(5.0, 3, 50, 1);
        assert_eq!(wl.dim(), 3);
        assert_eq!(wl.state_at(0.3), 1);
        assert_eq!(wl.state_at(2.5), 1);
        assert_eq!(wl.state_at(4.9), 1);
        assert_consistent(&wl);
    }

    #[test]
    fn bosonic_dim10() {
        let wl = DiscreteWorldline::new(8.0, 10, 80, 5);
        assert_eq!(wl.dim(), 10);
        assert_eq!(wl.state_at(3.3), 5);
        assert_eq!(wl.diagonal(), 5.0);
        assert_consistent(&wl);
    }

    #[test]
    fn spin_five_half_dim6() {
        let wl = DiscreteWorldline::new(6.0, 6, 60, 3);
        assert_eq!(wl.dim(), 6);
        assert_consistent(&wl);
    }

    #[test]
    fn many_slices() {
        let wl = DiscreteWorldline::new(1.0, 5, 1000, 2);
        assert_eq!(wl.state_at(0.0005), 2);
        assert_eq!(wl.state_at(0.5), 2);
        assert_eq!(wl.state_at(0.9995), 2);
    }
}
