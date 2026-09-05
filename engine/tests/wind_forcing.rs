//! Acceptance tests for T-03.1 — the `WindStress` trait, the steady
//! trade-wind scenario, and the thermocline tilt it drives.
//!
//! Two things are under test. The first is the forcing itself: a
//! `SteadyTradeWinds` is an easterly stress that does not change with `x` or
//! `t` and decays away from the equator as a Gaussian, and sampling it onto
//! the Arakawa C-grid puts `τx` on the east/west faces and `τy` on the
//! north/south ones. Every expected value there is the closed form
//! `τ₀·exp(−(y/Ly)²)` written out from the module's own definition, evaluated
//! independently in the test.
//!
//! The second is the ticket's acceptance criterion: *steady forcing run to
//! equilibrium produces a thermocline that is deeper in the west than the
//! east, sanity-checked against the analytic tilt formula.*
//!
//! # The analytic tilt
//!
//! Take the equations of `docs/planning/01-scientific-model.md` in the limit
//! the uniform trade winds admit: a stress independent of `y` is balanced by a
//! thermocline independent of `y`, so `∂h/∂y = 0`, `τy = 0`, and the whole
//! problem collapses onto the zonal pair
//!
//! ```text
//! 0 = −g'·∂h/∂x + τx/(ρ₀·H) − r·u
//! 0 = −H·∂u/∂x                − r·h
//! ```
//!
//! in a closed basin `0 ≤ x ≤ L`, where no-normal-flow reads `u(0) = u(L) = 0`.
//! Eliminating `u` between the two gives `∂²h/∂x² = k²·h` with `k = r/c` and
//! `c = √(g'·H)` the Kelvin wave speed, and the two wall conditions pick out
//! the odd solution about mid-basin:
//!
//! ```text
//!            τx                     
//! h(x) = ───────────────── · sinh(k·(x − L/2))
//!        ρ₀·H·g'·k·cosh(kL/2)
//! ```
//!
//! This is the damped Stommel-type balance of Cane & Sarachik (1981) and of
//! Gill, *Atmosphere–Ocean Dynamics* § 11.7, and it reduces to the familiar
//! undamped tilt `h(x) = τx·(x − L/2)/(ρ₀·H·g')` as `k → 0`. With `τx < 0` it
//! is positive in the west and negative in the east: the deep-west, shallow-
//! east mean state of `CONTEXT.md`, *Thermocline tilt*.
//!
//! Two things the two-dimensional model does that this one-dimensional balance
//! does not, and how each is accounted for:
//!
//! - **Rotation.** The `f·v` coupling has no counterpart in the balance above,
//!   and it does not vanish on a beta-plane: `|f| = β·|y|` is zero only on the
//!   equator itself. Its influence therefore shrinks with the basin's
//!   meridional extent, and
//!   [`the_tilt_approaches_the_analytic_balance_as_the_basin_narrows`]
//!   is the convergence test that says so — the honest form of this check per
//!   CODING_STANDARDS.md § Tests, rather than one threshold on one basin.
//! - **Discretisation.** The centred second difference solves
//!   `2·(cosh(k̃·Δx) − 1) = (k·Δx)²` rather than `k̃ = k`, so the discrete decay
//!   rate is off by `k̃/k − 1 ≈ −(k·Δx)²/24`. That is what sets the fixed
//!   tolerance of the point check below.

use engine::{
    max_stable_dt, Basin, BasinError, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing,
    SteadyTradeWinds, WaveSpeed, WindStress, WindStressError, WindStressField, H_STAGGERING,
    U_STAGGERING, V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;
/// Reference seawater density `ρ₀`, in kg/m³ — `CONTEXT.md` and Gill, appendix 3.
const REFERENCE_DENSITY_KG_PER_M3: f64 = engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3;

/// Rayleigh damping `r` the tilt runs use, in s⁻¹: an `e`-folding time of about
/// 11.6 days. Far stronger than the equatorial Pacific's own damping, for the
/// reason `rayleigh_damping.rs` spells out — the basin has to reach its
/// equilibrium inside a run of CFL-admissible steps.
const STRONG_DAMPING_PER_S: f64 = 1.0e-6;

/// Zonal wind stress `τ₀` of the trade-wind scenarios, in Pa. Easterly
/// trade-wind stress is `τx < 0` (`CONTEXT.md`), and 0.05 Pa is the observed
/// scale of the equatorial Pacific's mean zonal stress.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*).
const BASIN_LX_M: f64 = 1.0e7;
/// Cell height of the test basins, in metres. Held fixed while the meridional
/// cell *count* varies, so that the convergence test changes the basin's
/// reach away from the equator — and therefore `|f| = β·|y|` — and nothing
/// else about the discretisation.
const BASIN_DY_M: f64 = 1.0e5;
/// Zonal cell count of the test basins.
const BASIN_NX: usize = 40;

/// One solar day, in seconds — the shortest interval the steadiness check
/// samples over.
const DAY_S: f64 = 86_400.0;
/// One tropical year, in seconds (365.24 days). The period T-03.2's seasonal
/// modulation will carry, and therefore the interval a *steady* scenario must
/// be unmoved by.
const YEAR_S: f64 = 365.24 * DAY_S;
/// Ten years, in seconds — several ENSO cycles, and long enough that a slow
/// drift in a supposedly steady field would show.
const DECADE_S: f64 = 10.0 * YEAR_S;

/// Meridional positions the trade-wind profile is probed at, in metres north
/// of the equator: on the equator, one test cell off it, a quarter of the way
/// to the pole, and far outside any basin this model uses.
const PROBE_LATITUDES_M: [f64; 4] = [0.0, BASIN_DY_M, -2.5e6, 1.0e7];

/// Meridional cell counts the tilt is measured on, in decreasing order. The
/// basin reaches `±ny·Δy/2`, so each entry halves the largest `|f| = β·|y|`
/// the run sees, and with it the rotation coupling the one-dimensional
/// analytic balance omits.
const NARROWING_BASIN_CELLS: [usize; 3] = [8, 4, 2];

/// Length of a run to equilibrium, in seconds — about 350 days.
///
/// Two timescales have to have expired: the basin's adjustment time
/// `L/c = 3.7×10⁶ s`, and the damping time `1/r = 10⁶ s`. This is eight of the
/// first and thirty of the second, so the slowest transient is down by
/// `exp(−30) ≈ 10⁻¹³` of its initial size — far below every tolerance below.
const RUN_TO_EQUILIBRIUM_S: f64 = 3.0e7;

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// tens of ulps of `f64` (ε ≈ 2.2×10⁻¹⁶) for the handful of operations per
/// point the expression costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

/// The equatorial deformation radius `Le = √(c/β)`, in metres — the
/// meridional scale over which equatorial waves decay away from the equator
/// (`CONTEXT.md`, *Equatorial deformation radius*), and the natural width for
/// a wind field meant to drive that waveguide. About 3.45×10⁵ m for the
/// parameters below.
fn equatorial_deformation_radius_m(params: PhysicalParams) -> f64 {
    (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt()
}

/// The equatorial-Pacific parameter set at a given Rayleigh damping.
fn pacific_params(rayleigh_damping_per_s: f64) -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        rayleigh_damping_per_s,
        BETA_PER_M_PER_S,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// A basin [`BASIN_LX_M`] wide and `ny` cells tall, centred on the equator.
fn equatorial_basin(ny: usize) -> Basin {
    let grid = Grid::new(BASIN_NX, ny).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / BASIN_NX as f64, BASIN_DY_M)
        .expect("a basin spanned by whole cells has positive spacing");
    Basin::centered_on_equator(grid, spacing)
}

// --- The trade-wind field itself. ---

#[test]
fn the_trade_winds_blow_easterly_on_the_equator() {
    // The sign convention of `CONTEXT.md`: the alizés blow from the east, so
    // `τx < 0`. On the equator the Gaussian is exactly one, so the stress is
    // the configured amplitude itself.
    let winds = SteadyTradeWinds::with_meridional_decay(
        TRADE_WIND_STRESS_PA,
        equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S)),
    )
    .expect("an easterly stress with a positive decay scale");
    let (tau_x_pa, _) = winds.stress(0.0, 0.0, 0.0);

    assert!(
        tau_x_pa < 0.0,
        "trade winds must be easterly, got {tau_x_pa}"
    );
    assert_eq!(tau_x_pa, TRADE_WIND_STRESS_PA);
}

#[test]
fn the_trade_winds_decay_away_from_the_equator_as_a_gaussian() {
    // Checked against `τ₀·exp(−(y/Ly)²)` evaluated here, not against the
    // module's own arithmetic: at one decay scale the stress is `τ₀/e`, at two
    // it is `τ₀/e⁴`, and it is symmetric about the equator because `y` enters
    // squared.
    let decay_scale_m = equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S));
    let winds = SteadyTradeWinds::with_meridional_decay(TRADE_WIND_STRESS_PA, decay_scale_m)
        .expect("an easterly stress with a positive decay scale");

    for (scales_off_equator, expected_fraction) in [
        (0.0, 1.0),
        (1.0, 1.0 / std::f64::consts::E),
        (2.0, (-4.0_f64).exp()),
        (3.0, (-9.0_f64).exp()),
    ] {
        let expected_pa = TRADE_WIND_STRESS_PA * expected_fraction;
        for sign in [-1.0, 1.0] {
            let y_m = sign * scales_off_equator * decay_scale_m;
            let (tau_x_pa, _) = winds.stress(0.0, y_m, 0.0);
            assert!(
                (tau_x_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
                "at y = {y_m} m the stress is {tau_x_pa} Pa, expected {expected_pa} Pa"
            );
        }
    }
}

#[test]
fn uniform_trade_winds_have_no_meridional_structure() {
    // The `Ly → ∞` limit: the same stress at the equator and a quarter of the
    // globe away from it.
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");

    for y_m in PROBE_LATITUDES_M {
        assert_eq!(winds.stress(0.0, y_m, 0.0).0, TRADE_WIND_STRESS_PA);
    }
}

#[test]
fn the_trade_winds_are_zonal_and_steady() {
    // "Steady" is the whole name of the scenario: the same stress at every
    // `x`, and at every `t` from the first step to ten years in. The seasonal
    // modulation is T-03.2's.
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let winds = SteadyTradeWinds::with_meridional_decay(
        TRADE_WIND_STRESS_PA,
        equatorial_deformation_radius_m(params),
    )
    .expect("an easterly stress with a positive decay scale");
    // Probed off the equator, where the profile actually varies, so a stress
    // that drifted with `x` or `t` could not hide behind the Gaussian's peak.
    let probe_y_m = BASIN_DY_M;
    let reference = winds.stress(0.0, probe_y_m, 0.0);

    for x_m in [0.0, BASIN_LX_M / 2.0, BASIN_LX_M] {
        for t_s in [0.0, DAY_S, YEAR_S, DECADE_S] {
            let stress = winds.stress(x_m, probe_y_m, t_s);
            assert_eq!(stress, reference, "at x = {x_m} m, t = {t_s} s");
            assert_eq!(stress.1, 0.0, "the alizés carry no meridional stress");
        }
    }
}

#[test]
fn a_wind_that_is_not_easterly_is_refused_by_name() {
    // Invalid scenario input is a `Result` naming the offending value, not a
    // panic and not a silently flipped sign (CODING_STANDARDS.md
    // § Correctness and failure).
    for value_pa in [0.0, 0.05, f64::NAN, f64::NEG_INFINITY] {
        let error = SteadyTradeWinds::uniform(value_pa)
            .expect_err("only a strictly negative stress is a trade wind");
        let WindStressError::NotEasterly {
            value_pa: rejected_pa,
        } = error
        else {
            panic!("expected the stress itself to be rejected, got {error}");
        };
        // Compared bitwise rather than with `==`, so that the NaN case checks
        // the value was carried through rather than trivially passing.
        assert_eq!(rejected_pa.to_bits(), value_pa.to_bits());
        let message = error.to_string();
        assert!(message.contains("strictly negative"), "{message}");
    }
}

#[test]
fn a_decay_scale_that_is_not_a_distance_is_refused_by_name() {
    let negative_scale_m = -equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S));
    for value_m in [0.0, negative_scale_m, f64::NAN, f64::INFINITY] {
        let error = SteadyTradeWinds::with_meridional_decay(TRADE_WIND_STRESS_PA, value_m)
            .expect_err("a decay scale must be a finite, positive distance");
        let WindStressError::ScaleNotPositive {
            parameter,
            value_m: rejected_m,
        } = error
        else {
            panic!("expected the decay scale to be rejected, got {error}");
        };
        assert_eq!(parameter, "meridional_decay_scale_m");
        assert_eq!(rejected_m.to_bits(), value_m.to_bits());
        let message = error.to_string();
        assert!(message.contains("meridional_decay_scale_m"), "{message}");
    }
}

#[test]
fn a_basin_at_a_position_that_is_not_a_position_is_refused_by_name() {
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[1]);
    let error = Basin::new(basin.grid(), basin.spacing(), f64::NAN, 0.0)
        .expect_err("a basin edge must be a finite position");
    let BasinError::NotFinite { parameter, value_m } = error;
    assert_eq!(parameter, "western_edge_x_m");
    assert!(
        value_m.is_nan(),
        "the rejected value must be carried through"
    );

    let error = Basin::new(basin.grid(), basin.spacing(), 0.0, f64::INFINITY)
        .expect_err("a basin edge must be a finite position");
    assert_eq!(
        error,
        BasinError::NotFinite {
            parameter: "southern_edge_y_m",
            value_m: f64::INFINITY,
        }
    );
}

// --- Sampling the trait onto the C-grid. ---

#[test]
fn the_basin_and_the_beta_plane_agree_on_where_each_row_sits() {
    // The forcing and the rotation must not disagree about which row is the
    // equator, or the wind would be centred on one latitude and the waveguide
    // on another.
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[0]);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());

    assert_eq!(basin.southern_edge_y_m(), plane.southern_edge_y_m());
    for staggering in [H_STAGGERING, U_STAGGERING, V_STAGGERING] {
        for j in 0..=basin.grid().ny() {
            assert_eq!(
                basin.y_of_row_m(staggering, j),
                plane.y_of_row_m(staggering, j),
                "row {j} at {staggering:?}"
            );
        }
    }
}

#[test]
fn sampling_puts_each_component_on_the_faces_its_equation_lives_on() {
    // `τx` accelerates `u`, which lives on the east/west faces; `τy`
    // accelerates `v`, on the north/south ones. Each interior face must carry
    // the trait's value at that face's own position, evaluated here from the
    // Gaussian rather than read back from the field.
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let decay_scale_m = equatorial_deformation_radius_m(params);
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[0]);
    let winds = SteadyTradeWinds::with_meridional_decay(TRADE_WIND_STRESS_PA, decay_scale_m)
        .expect("an easterly stress with a positive decay scale");

    let field = WindStressField::sampled(basin, &winds, 0.0);

    let (nx, ny) = (basin.grid().nx(), basin.grid().ny());
    assert_eq!(field.grid(), basin.grid());
    for j in 0..field.tau_x_pa().ny() {
        let y_m = basin.y_of_row_m(U_STAGGERING, j);
        let scaled = y_m / decay_scale_m;
        let expected_pa = TRADE_WIND_STRESS_PA * (-scaled * scaled).exp();
        for i in 1..nx {
            let sampled_pa = *field.tau_x_pa().get(i, j).expect("an interior face");
            assert!(
                (sampled_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
                "τx at face ({i}, {j}) is {sampled_pa} Pa, expected {expected_pa} Pa"
            );
        }
    }
    for j in 0..=ny {
        for i in 0..nx {
            assert_eq!(
                *field.tau_y_pa().get(i, j).expect("an in-bounds face"),
                0.0,
                "the alizés carry no meridional stress"
            );
        }
    }
}

#[test]
fn the_closed_basins_walls_carry_no_sampled_stress() {
    // The rule the `forcing` module header derives at length: a wall face has
    // water on one side only, so a sampled field leaves it at exactly zero. It
    // is what keeps a wind-driven closed basin closed until T-04.2 gives the
    // boundary a condition of its own — without it the coasts pass water and
    // the basin never tilts at all.
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[0]);
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");

    let field = WindStressField::sampled(basin, &winds, 0.0);

    let (nx, ny) = (basin.grid().nx(), basin.grid().ny());
    for j in 0..field.tau_x_pa().ny() {
        assert_eq!(*field.tau_x_pa().get(0, j).expect("the western wall"), 0.0);
        assert_eq!(*field.tau_x_pa().get(nx, j).expect("the eastern wall"), 0.0);
    }
    for i in 0..field.tau_y_pa().nx() {
        assert_eq!(*field.tau_y_pa().get(i, 0).expect("the southern wall"), 0.0);
        assert_eq!(
            *field.tau_y_pa().get(i, ny).expect("the northern wall"),
            0.0
        );
    }
}

#[test]
fn re_sampling_in_place_writes_every_point() {
    // A time-varying scenario re-samples one buffer per RK4 stage rather than
    // allocating a field per stage (CODING_STANDARDS.md § Performance), so no
    // point may survive from the previous contents.
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[0]);
    let winds = SteadyTradeWinds::with_meridional_decay(
        TRADE_WIND_STRESS_PA,
        equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S)),
    )
    .expect("an easterly stress with a positive decay scale");
    let expected = WindStressField::sampled(basin, &winds, 0.0);

    let mut field = WindStressField::uniform_including_walls(basin.grid(), 1.0, -1.0);
    field.sample(basin, &winds, 0.0);

    assert_eq!(field, expected);
}

#[test]
#[should_panic(expected = "wind stress field was built for")]
fn sampling_over_the_wrong_basin_panics() {
    // A shape mismatch means the calling code is wrong, which is what panics
    // are for (CODING_STANDARDS.md § Correctness and failure).
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");
    let mut field = WindStressField::calm(equatorial_basin(NARROWING_BASIN_CELLS[0]).grid());
    field.sample(equatorial_basin(NARROWING_BASIN_CELLS[1]), &winds, 0.0);
}

// --- The acceptance criterion: the thermocline tilts. ---

/// The steady thermocline anomaly of the one-dimensional damped balance, in
/// metres, at `x_m` metres east of the western boundary.
///
/// The closed form derived in this file's header;
/// `docs/planning/01-scientific-model.md` names the same balance as the
/// Sverdrup/Stommel-type validation target of Epic 07.
fn analytic_tilt_m(x_m: f64, basin: Basin, params: PhysicalParams, tau_x_pa: f64) -> f64 {
    let decay_per_m = params.rayleigh_damping_per_s() / params.kelvin_wave_speed_m_per_s();
    let basin_width_m = basin.zonal_extent_m();
    let amplitude_m = tau_x_pa
        / (params.reference_density_kg_per_m3()
            * params.mean_thermocline_depth_m()
            * params.reduced_gravity_m_per_s2()
            * decay_per_m
            * (decay_per_m * basin_width_m / 2.0).cosh());
    amplitude_m * (decay_per_m * (x_m - basin_width_m / 2.0)).sinh()
}

/// Run `basin` from rest under `winds` for `run_s` seconds, and return the
/// final state.
fn run_to_equilibrium(
    basin: Basin,
    params: PhysicalParams,
    winds: &impl WindStress,
    run_s: f64,
) -> OceanState {
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    let dt_s = max_stable_dt(basin.spacing(), wave_speed);
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());
    let mut solver = Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
        .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"));

    // Driven by the trait itself, re-sampled at every RK4 stage: this is the
    // path a scenario takes, and the one that makes the forcing plumbing of
    // this ticket load-bearing rather than decorative.
    let mut state = OceanState::at_rest(basin.grid());
    let steps = (run_s / dt_s).ceil() as usize;
    for step in 0..steps {
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, winds);
    }
    state
}

/// The meridional mean of `h` in the column `i`, in metres — the quantity the
/// `y`-independent analytic balance is a statement about.
fn column_mean_h_m(state: &OceanState, i: usize) -> f64 {
    let ny = state.grid().ny();
    (0..ny)
        .map(|j| *state.h().get(i, j).expect("an in-bounds cell"))
        .sum::<f64>()
        / ny as f64
}

/// The east–west thermocline tilt of `state`, in metres: how much deeper the
/// westernmost column is than the easternmost one. Positive is the observed
/// mean state (`CONTEXT.md`, *Thermocline tilt*).
fn tilt_m(state: &OceanState) -> f64 {
    column_mean_h_m(state, 0) - column_mean_h_m(state, state.grid().nx() - 1)
}

#[test]
fn steady_trade_winds_leave_the_thermocline_deeper_in_the_west() {
    // The ticket's acceptance criterion, on the narrowest of the test basins —
    // the one where the rotation coupling the analytic balance omits is
    // smallest. Both halves are checked: the sign, which is the criterion
    // itself, and the magnitude against the closed form.
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[NARROWING_BASIN_CELLS.len() - 1]);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");

    let state = run_to_equilibrium(basin, params, &winds, RUN_TO_EQUILIBRIUM_S);

    let west_m = column_mean_h_m(&state, 0);
    let east_m = column_mean_h_m(&state, basin.grid().nx() - 1);
    assert!(
        west_m > east_m,
        "easterly trade winds must deepen the west and shoal the east, got west {west_m} m \
         and east {east_m} m"
    );

    // The centred second difference of the discrete balance decays at `k̃`
    // rather than `k`, with `k̃/k − 1 ≈ −(k·Δx)²/24 = −3.5×10⁻⁴` for this
    // basin; 10⁻³ is that truncation rounded up to the next order of
    // magnitude. The rotation coupling is the other error source, and the
    // convergence test below is what pins it.
    const DISCRETE_DECAY_TOLERANCE: f64 = 1.0e-3;
    let x_of = |i: usize| basin.x_of_column_m(H_STAGGERING, i);
    let analytic_tilt_m = analytic_tilt_m(x_of(0), basin, params, TRADE_WIND_STRESS_PA)
        - analytic_tilt_m(
            x_of(basin.grid().nx() - 1),
            basin,
            params,
            TRADE_WIND_STRESS_PA,
        );
    let measured_tilt_m = tilt_m(&state);
    let relative_error = (measured_tilt_m - analytic_tilt_m).abs() / analytic_tilt_m.abs();
    assert!(
        relative_error <= DISCRETE_DECAY_TOLERANCE,
        "the tilt is {measured_tilt_m} m against an analytic {analytic_tilt_m} m, a relative \
         error of {relative_error}"
    );
}

#[test]
fn the_equilibrium_is_an_equilibrium() {
    // "Run to equilibrium" has to mean something: after the run, further steps
    // must not move the thermocline. The slowest transient decays like
    // `exp(−r·t)`, which after `RUN_TO_EQUILIBRIUM_S` is `exp(−30) ≈ 10⁻¹³` of
    // its initial size, so 10⁻¹⁰ of the tilt is a bound the state must already
    // be well inside.
    const RESIDUAL_DRIFT_TOLERANCE: f64 = 1.0e-10;
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[NARROWING_BASIN_CELLS.len() - 1]);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");

    let settled = run_to_equilibrium(basin, params, &winds, RUN_TO_EQUILIBRIUM_S);
    let settled_further = run_to_equilibrium(basin, params, &winds, 2.0 * RUN_TO_EQUILIBRIUM_S);

    let drift_m = (tilt_m(&settled_further) - tilt_m(&settled)).abs();
    assert!(
        drift_m <= RESIDUAL_DRIFT_TOLERANCE * tilt_m(&settled).abs(),
        "the tilt moved {drift_m} m over a second run of the same length; it has not settled"
    );
}

#[test]
fn the_tilt_approaches_the_analytic_balance_as_the_basin_narrows() {
    // Convergence rather than a point check (CODING_STANDARDS.md § Tests). The
    // one term the analytic balance drops is the rotation coupling `f·v`, and
    // `|f| = β·|y|` is bounded by `β·ny·Δy/2`: halving the basin's meridional
    // cell count halves the largest rotation rate the run sees, so the
    // discrepancy must shrink with it. Requiring each error to be at most half
    // the previous is the weakest statement of that — the coupling enters the
    // balance at least linearly in `|f|`.
    const REQUIRED_ERROR_REDUCTION: f64 = 0.5;
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let winds = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
        .expect("an easterly stress is a trade wind");

    let errors: Vec<f64> = NARROWING_BASIN_CELLS
        .iter()
        .map(|&ny| {
            let basin = equatorial_basin(ny);
            let state = run_to_equilibrium(basin, params, &winds, RUN_TO_EQUILIBRIUM_S);
            let x_of = |i: usize| basin.x_of_column_m(H_STAGGERING, i);
            let analytic_m = analytic_tilt_m(x_of(0), basin, params, TRADE_WIND_STRESS_PA)
                - analytic_tilt_m(
                    x_of(basin.grid().nx() - 1),
                    basin,
                    params,
                    TRADE_WIND_STRESS_PA,
                );
            (tilt_m(&state) - analytic_m).abs() / analytic_m.abs()
        })
        .collect();

    for (pair, ny) in errors.windows(2).zip(NARROWING_BASIN_CELLS.windows(2)) {
        assert!(
            pair[1] <= REQUIRED_ERROR_REDUCTION * pair[0],
            "halving the basin from {} to {} cells took the error from {} to {}, which is not \
             the reduction a rotation-driven discrepancy must show",
            ny[0],
            ny[1],
            pair[0],
            pair[1]
        );
    }
}

#[test]
fn trade_winds_that_decay_off_the_equator_still_tilt_the_thermocline() {
    // The other profile the ticket names. A stress that falls away from the
    // equator has no closed-form steady state — it drives a meridional
    // circulation the one-dimensional balance says nothing about — so this
    // checks the sign the criterion is about, on the widest of the test basins,
    // which reaches ±4×10⁵ m — a little over one deformation radius, so the
    // decay is visible across it.
    let basin = equatorial_basin(NARROWING_BASIN_CELLS[0]);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let winds = SteadyTradeWinds::with_meridional_decay(
        TRADE_WIND_STRESS_PA,
        equatorial_deformation_radius_m(params),
    )
    .expect("an easterly stress with a positive decay scale");

    let state = run_to_equilibrium(basin, params, &winds, RUN_TO_EQUILIBRIUM_S);

    let measured_tilt_m = tilt_m(&state);
    assert!(
        measured_tilt_m > 0.0,
        "easterly trade winds must deepen the west, got a tilt of {measured_tilt_m} m"
    );
}
