//! Acceptance tests for T-02.1 — the shallow-water state type and the fixed
//! physical parameters of the 1.5-layer reduced-gravity model.
//!
//! Every expected value here comes from an independent source: the C-grid
//! staggering fixed in ADR-0003, the definitions in `CONTEXT.md`, or published
//! equatorial-Pacific values. None of them was produced by running this code.

use engine::{
    Field2D, Grid, OceanState, PhysicalParams, PhysicalParamsError, StateVector, H_STAGGERING,
    U_STAGGERING, V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981): a density contrast of about
/// 5 kg/m³ across the thermocline over a reference density near 1025 kg/m³.
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: a damping timescale of about two years, the
/// order used for the linear equatorial-wave problem.
const PACIFIC_DAMPING_PER_S: f64 = 1.0 / (2.0 * 365.0 * 86_400.0);

/// The parameter set the derived-scale tests are anchored on.
fn pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        PACIFIC_DAMPING_PER_S,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

#[test]
fn at_rest_allocates_each_variable_at_its_c_grid_staggering() {
    // ADR-0003 puts `h` at cell centers, `u` at east/west faces and `v` at
    // north/south faces, so the face fields carry one extra line of points
    // along the axis they are staggered on. The shapes come from the grid
    // rather than from index arithmetic here, per CODING_STANDARDS.md
    // § Scope guards.
    let grid = Grid::new(4, 3).expect("a 4x3 basin is non-degenerate");
    let state = OceanState::at_rest(grid);

    assert_eq!(
        (state.h().nx(), state.h().ny()),
        grid.field_shape(H_STAGGERING)
    );
    assert_eq!(
        (state.u().nx(), state.u().ny()),
        grid.field_shape(U_STAGGERING)
    );
    assert_eq!(
        (state.v().nx(), state.v().ny()),
        grid.field_shape(V_STAGGERING)
    );
    assert_eq!(state.grid(), grid);
}

#[test]
fn the_rest_state_is_an_undisturbed_thermocline_at_its_mean_depth() {
    // `h` is an anomaly, never a total depth (CONTEXT.md): at rest every
    // anomaly is exactly zero, and the total thermocline depth is `H + 0 = H`.
    // Exactly zero, not "small" — this is a definition, not a computation, so
    // no tolerance is warranted.
    let grid = Grid::new(3, 2).expect("a 3x2 basin is non-degenerate");
    let state = OceanState::at_rest(grid);
    let params = pacific_params();

    assert!(state.h().as_slice().iter().all(|&h_m| h_m == 0.0));
    assert!(state.u().as_slice().iter().all(|&u| u == 0.0));
    assert!(state.v().as_slice().iter().all(|&v| v == 0.0));

    assert_eq!(
        state.total_thermocline_depth_m(&params, 2, 1),
        Some(PACIFIC_MEAN_DEPTH_M)
    );
    assert_eq!(state.total_thermocline_depth_m(&params, 3, 1), None);
}

#[test]
fn total_thermocline_depth_adds_the_anomaly_to_the_mean_depth() {
    // A 20 m deeper-than-average thermocline over a 150 m mean layer is a
    // total depth of 170 m — the sum is exact in binary floating point, so it
    // is asserted exactly.
    let grid = Grid::new(2, 2).expect("a 2x2 basin is non-degenerate");
    let mut state = OceanState::at_rest(grid);
    *state
        .h_mut()
        .get_mut(1, 0)
        .expect("(1, 0) is a cell center") = 20.0;

    assert_eq!(
        state.total_thermocline_depth_m(&pacific_params(), 1, 0),
        Some(PACIFIC_MEAN_DEPTH_M + 20.0)
    );
}

#[test]
fn physical_parameters_are_stored_in_si_units_unchanged() {
    // The acceptance criterion is SI throughout, so a parameter must come back
    // out in the unit its name states — no hidden rescaling on the way in.
    let params = pacific_params();
    assert_eq!(
        params.reduced_gravity_m_per_s2(),
        PACIFIC_REDUCED_GRAVITY_M_PER_S2
    );
    assert_eq!(params.mean_depth_m(), PACIFIC_MEAN_DEPTH_M);
    assert_eq!(params.rayleigh_damping_per_s(), PACIFIC_DAMPING_PER_S);
    assert_eq!(
        params.beta_per_m_per_s(),
        engine::EQUATORIAL_BETA_PER_M_PER_S
    );
    assert_eq!(
        params.reference_density_kg_per_m3(),
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3
    );
}

#[test]
fn the_kelvin_wave_speed_matches_the_observed_equatorial_pacific_value() {
    // `c = √(g'·H)` (CONTEXT.md). With g' = 0.05 m/s² and H = 150 m that is
    // √7.5 = 2.7386 m/s, and the observed phase speed of the first baroclinic
    // equatorial Kelvin wave in the Pacific is about 2.7 m/s (Gill,
    // *Atmosphere–Ocean Dynamics*; Cane & Sarachik 1981).
    let speed_m_per_s = pacific_params().kelvin_wave_speed_m_per_s();

    // 5%: the published figure is quoted to two significant figures.
    assert!(
        (speed_m_per_s - 2.7).abs() < 0.05 * 2.7,
        "c = {speed_m_per_s} m/s is not the observed ~2.7 m/s"
    );
    // Against the definition itself the only error is one sqrt rounding, a few
    // ulps — 1e-12 relative is orders of magnitude above that and still far
    // tighter than any physical tolerance.
    let analytic = 7.5_f64.sqrt();
    assert!((speed_m_per_s - analytic).abs() < 1e-12 * analytic);
}

#[test]
fn the_equatorial_deformation_radius_matches_the_published_scale() {
    // `Le = √(c/β)` (CONTEXT.md). With c = 2.7386 m/s and β = 2.3×10⁻¹¹
    // m⁻¹s⁻¹ that is √(1.1907×10¹¹) = 3.45×10⁵ m, and the equatorial
    // deformation radius of the Pacific is quoted as about 350 km.
    let radius_m = pacific_params().equatorial_deformation_radius_m();

    // 10%: the published 350 km is a one-significant-figure scale estimate,
    // and it depends on which baroclinic mode is quoted.
    assert!(
        (radius_m - 350.0e3).abs() < 0.10 * 350.0e3,
        "Le = {radius_m} m is not the observed ~350 km"
    );
    let analytic = (7.5_f64.sqrt() / engine::EQUATORIAL_BETA_PER_M_PER_S).sqrt();
    assert!((radius_m - analytic).abs() < 1e-12 * analytic);
}

#[test]
fn unphysical_parameters_are_rejected_naming_the_value_and_the_bound() {
    // A scenario carrying a negative mean depth is invalid *input*, so it is a
    // `Result` with an actionable message, per CODING_STANDARDS.md
    // § Correctness and failure.
    let err = PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        -150.0,
        PACIFIC_DAMPING_PER_S,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect_err("a thermocline cannot sit at a negative mean depth");
    assert_eq!(
        err,
        PhysicalParamsError::NotPositive {
            parameter: "mean_depth_m",
            value: -150.0,
        }
    );
    let message = err.to_string();
    assert!(message.contains("mean_depth_m"), "{message}");
    assert!(message.contains("-150"), "{message}");

    for (reduced_gravity, mean_depth, damping, beta, density) in [
        (
            0.0,
            PACIFIC_MEAN_DEPTH_M,
            PACIFIC_DAMPING_PER_S,
            engine::EQUATORIAL_BETA_PER_M_PER_S,
            engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
        ),
        (
            f64::NAN,
            PACIFIC_MEAN_DEPTH_M,
            PACIFIC_DAMPING_PER_S,
            engine::EQUATORIAL_BETA_PER_M_PER_S,
            engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
        ),
        (
            PACIFIC_REDUCED_GRAVITY_M_PER_S2,
            PACIFIC_MEAN_DEPTH_M,
            PACIFIC_DAMPING_PER_S,
            0.0,
            engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
        ),
        (
            PACIFIC_REDUCED_GRAVITY_M_PER_S2,
            PACIFIC_MEAN_DEPTH_M,
            PACIFIC_DAMPING_PER_S,
            engine::EQUATORIAL_BETA_PER_M_PER_S,
            0.0,
        ),
    ] {
        assert!(
            PhysicalParams::new(reduced_gravity, mean_depth, damping, beta, density).is_err(),
            "accepted an unphysical parameter set"
        );
    }

    let err = PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        -1.0e-8,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect_err("negative damping would amplify rather than damp");
    assert_eq!(
        err,
        PhysicalParamsError::Negative {
            parameter: "rayleigh_damping_per_s",
            value: -1.0e-8,
        }
    );
}

#[test]
fn zero_damping_is_admissible_because_the_undamped_limit_is_a_validation_target() {
    // `01-scientific-model.md` validates energy conservation in the `r = 0`,
    // `τ = 0` limit, so `r = 0` must be a constructible parameter set.
    let params = PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        0.0,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the undamped limit is a physical configuration");
    assert_eq!(params.rayleigh_damping_per_s(), 0.0);
}

#[test]
fn state_vector_arithmetic_is_componentwise_over_all_three_variables() {
    // Hand-computed: starting from 1.0 everywhere and adding 2.0 × 3.0 gives
    // 7.0 in every point of every variable. Exact in binary floating point.
    let grid = Grid::new(2, 2).expect("a 2x2 basin is non-degenerate");
    let mut state = uniform_state(grid, 1.0);
    let other = uniform_state(grid, 3.0);

    state.add_scaled(2.0, &other);
    for values in [
        state.h().as_slice(),
        state.u().as_slice(),
        state.v().as_slice(),
    ] {
        assert!(values.iter().all(|&value| value == 7.0), "{values:?}");
    }

    state.assign(&other);
    for values in [
        state.h().as_slice(),
        state.u().as_slice(),
        state.v().as_slice(),
    ] {
        assert!(values.iter().all(|&value| value == 3.0), "{values:?}");
    }
}

#[test]
#[should_panic(expected = "grid")]
fn combining_states_over_different_basins_panics_rather_than_truncating() {
    // The shape of an `OceanState` is not in its type, so `StateVector`'s
    // contract requires a panic rather than a silent truncation: a mismatch
    // means the calling code is wrong.
    let mut state = OceanState::at_rest(Grid::new(2, 2).expect("non-degenerate"));
    let other = OceanState::at_rest(Grid::new(3, 2).expect("non-degenerate"));
    state.add_scaled(1.0, &other);
}

/// A state with every point of every variable set to `value` — a shape no
/// physical state has, used only to check the vector-space arithmetic.
fn uniform_state(grid: Grid, value: f64) -> OceanState {
    let mut state = OceanState::at_rest(grid);
    fill(state.h_mut(), value);
    fill(state.u_mut(), value);
    fill(state.v_mut(), value);
    state
}

fn fill(field: &mut Field2D<f64>, value: f64) {
    for slot in field.as_mut_slice() {
        *slot = value;
    }
}
