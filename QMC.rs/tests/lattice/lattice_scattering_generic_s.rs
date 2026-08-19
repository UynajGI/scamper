//! Generic-spin-S scattering-table identities (criterion A).
//!
//! The deterministic directed-loop balance validation previously covered only
//! S=1 (custom catalog) and S=2 at 1e-10 (XXZ). Production use spans arbitrary
//! integer and half-integer spins on every supported model family, so the
//! machine-precision row-normalization and local detailed-balance identities
//! are checked here for `S ∈ {1/2, 1, 3/2, 2, 5/2}` across the model catalog,
//! for **both** user-selectable scattering policies.
//!
//! The identities are purely algebraic (exact integer `2S` ladder algebra),
//! so any violation is a bug, not noise: tolerance 1e-12, no statistics.

use qmc_rs::{
    CsrGraph, EdgeCoupling, LocalHilbertSpace, ScatteringPolicy, SiteCoupling, SpinLatticeModel,
    SpinModelBuilder, SpinSpace,
};

const POLICIES: [ScatteringPolicy; 2] = [ScatteringPolicy::LowBounce, ScatteringPolicy::Metropolis];

/// Generic spins covered: 1/2, 1, 3/2, 2, 5/2 (integer `2S`).
const TWO_S: [u16; 5] = [1, 2, 3, 4, 5];

fn assert_table_balanced(label: &str, model: &SpinLatticeModel) {
    assert!(
        !model.terms().is_empty(),
        "{label}: model compiled to zero operator terms"
    );
    for term in model.terms() {
        let diagnostics = term.scattering().diagnostics();
        assert!(
            diagnostics.max_row_error < 1.0e-12,
            "{label} / term `{}`: row normalization error {} exceeds 1e-12",
            term.label(),
            diagnostics.max_row_error
        );
        assert!(
            diagnostics.max_detailed_balance_error < 1.0e-12,
            "{label} / term `{}`: local detailed-balance residual {} exceeds 1e-12",
            term.label(),
            diagnostics.max_detailed_balance_error
        );
        assert!(
            (0.0..=1.0).contains(&diagnostics.mean_bounce_probability),
            "{label} / term `{}`: mean bounce probability {} outside [0, 1]",
            term.label(),
            diagnostics.mean_bounce_probability
        );
    }
}

fn check_model(label: &str, build: impl Fn(ScatteringPolicy) -> SpinLatticeModel) {
    for policy in POLICIES {
        let model = build(policy);
        assert_table_balanced(&format!("{label} ({policy:?})"), &model);
    }
}

fn open_chain(sites: usize) -> CsrGraph {
    // Open chain: bipartite, so antiferromagnetic off-diagonal terms admit
    // the Marshall gauge at every spin magnitude.
    CsrGraph::chain(sites, false).expect("chain")
}

#[test]
fn heisenberg_scattering_balance_holds_for_generic_spin() {
    for &two_s in &TWO_S {
        for &j in &[1.0, -1.0] {
            check_model(
                &format!("Heisenberg S={}/2 J={j}", f64::from(two_s)),
                |policy| {
                    SpinModelBuilder::new(
                        open_chain(3),
                        SpinSpace::uniform(3, two_s).expect("space"),
                    )
                    .name("Heisenberg")
                    .uniform_edge(EdgeCoupling::heisenberg(j))
                    .scattering_policy(policy)
                    .build()
                    .expect("model")
                },
            );
        }
    }
}

#[test]
fn xxz_and_xy_scattering_balance_holds_for_generic_spin() {
    for &two_s in &TWO_S {
        check_model(&format!("XXZ S={}/2", f64::from(two_s)), |policy| {
            SpinModelBuilder::new(open_chain(3), SpinSpace::uniform(3, two_s).expect("space"))
                .name("XXZ")
                .uniform_edge(EdgeCoupling::xxz(0.8, 1.3))
                .scattering_policy(policy)
                .build()
                .expect("model")
        });
        check_model(&format!("XY S={}/2", f64::from(two_s)), |policy| {
            SpinModelBuilder::new(open_chain(4), SpinSpace::uniform(4, two_s).expect("space"))
                .name("XY")
                .uniform_edge(EdgeCoupling::xxz(1.1, 0.0))
                .scattering_policy(policy)
                .build()
                .expect("model")
        });
    }
}

#[test]
fn xyz_pair_flip_scattering_balance_holds_for_generic_spin() {
    // |Jx| > |Jy| keeps the exchange and pair-flip gauge constraints
    // compatible on a bipartite chain.
    for &two_s in &TWO_S {
        check_model(&format!("XYZ S={}/2", f64::from(two_s)), |policy| {
            SpinModelBuilder::new(open_chain(3), SpinSpace::uniform(3, two_s).expect("space"))
                .name("XYZ")
                .uniform_edge(EdgeCoupling::xyz(1.2, 0.5, 0.9))
                .scattering_policy(policy)
                .build()
                .expect("model")
        });
    }
}

#[test]
fn transverse_field_ising_scattering_balance_holds_for_generic_spin() {
    for &two_s in &TWO_S {
        check_model(&format!("TFIM S={}/2", f64::from(two_s)), |policy| {
            SpinModelBuilder::new(open_chain(3), SpinSpace::uniform(3, two_s).expect("space"))
                .name("TFIM")
                .uniform_edge(EdgeCoupling::xxz(0.0, -1.0))
                .uniform_site(SiteCoupling::new(0.6, 0.0, 0.0))
                .scattering_policy(policy)
                .build()
                .expect("model")
        });
    }
}

#[test]
fn single_ion_and_longitudinal_field_terms_balance_for_generic_spin() {
    // Site terms with single-ion anisotropy D and longitudinal field h_z do
    // not constrain the Marshall gauge, so they compose with exchange edges.
    for &two_s in &TWO_S {
        check_model(
            &format!("Heisenberg+D+h_z S={}/2", f64::from(two_s)),
            |policy| {
                SpinModelBuilder::new(open_chain(3), SpinSpace::uniform(3, two_s).expect("space"))
                    .name("AnisotropicSite")
                    .uniform_edge(EdgeCoupling::heisenberg(1.0))
                    .uniform_site(SiteCoupling::new(0.0, 0.25, -0.4))
                    .scattering_policy(policy)
                    .build()
                    .expect("model")
            },
        );
    }
}

#[test]
fn site_resolved_mixed_spin_model_balance_holds() {
    // Mixed magnitudes exercise different leg dimensions inside one term.
    let two_s_by_site = [1_u16, 4, 3, 5];
    let space = SpinSpace::site_resolved(two_s_by_site.to_vec()).expect("mixed spin space");
    check_model("mixed-spin Heisenberg", |policy| {
        SpinModelBuilder::new(open_chain(4), space.clone())
            .name("MixedSpin")
            .uniform_edge(EdgeCoupling::heisenberg(0.9))
            .scattering_policy(policy)
            .build()
            .expect("model")
    });
}

#[test]
fn generic_spin_ladder_algebra_satisfies_exact_sum_rules() {
    // Deterministic identity layer beneath the tables: the integer-2S ladder
    // amplitudes must reproduce the closed-form angular-momentum sum rules
    //   sum_m m = 0
    //   sum_m m^2 = S(S+1)(2S+1)/3
    //   sum_m [S(S+1) - m(m+1)] = (2/3) S(S+1)(2S+1)
    // at machine precision for every supported generic spin.
    for &two_s in &TWO_S {
        let space = SpinSpace::uniform(1, two_s).expect("space");
        let dimension = space.dimension(0);
        let s = space.spin(0);
        let mut sum_m = 0.0;
        let mut sum_m2 = 0.0;
        let mut sum_raising = 0.0;
        for state in 0..dimension {
            let m = space.m(0, state as u16);
            sum_m += m;
            sum_m2 += m * m;
            if let Some(amplitude) = space.raising_amplitude(0, state as u16) {
                sum_raising += amplitude * amplitude;
            }
        }
        assert!(
            sum_m.abs() < 1.0e-12,
            "S={s}: sum of m values is {sum_m}, expected 0"
        );
        let expected_m2 = s * (s + 1.0) * (2.0 * s + 1.0) / 3.0;
        assert!(
            (sum_m2 - expected_m2).abs() < 1.0e-12,
            "S={s}: sum m^2 = {sum_m2}, expected {expected_m2}"
        );
        let expected_raising = (2.0 / 3.0) * s * (s + 1.0) * (2.0 * s + 1.0);
        assert!(
            (sum_raising - expected_raising).abs() < 1.0e-12,
            "S={s}: ladder sum rule = {sum_raising}, expected {expected_raising}"
        );
        // Half-integer spins must carry no floating-point rounding convention.
        assert_eq!(space.m_twice(0, 0), -i32::from(two_s));
        assert_eq!(space.m_twice(0, two_s), i32::from(two_s));
    }
}
