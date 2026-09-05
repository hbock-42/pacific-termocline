//! T-12.2 — the atmosphere's answer to the SST anomaly, and the proof that a
//! run with the feedback switched off is the run this project already
//! validated.
//!
//! The response under test is the statistical, Gill-type one of
//! `docs/planning/01-scientific-model.md` § *Phase 2*: the tropical atmosphere
//! adjusts to an equatorial heating anomaly far faster than the ocean adjusts
//! to the wind, so its wind anomaly is a *diagnostic* function of the SST
//! anomaly of the instant, and the equatorially trapped zonal wind of Gill's
//! (*Q. J. R. Meteorol. Soc.* 106, 1980) solution is a Gaussian about the
//! equator:
//!
//! ```text
//! τx'(x, y, t) = μ · ⟨T'⟩(t) · exp(−(y/L_a)²)      τy' = 0
//! ⟨T'⟩ = Σ T'ᵢⱼ·exp(−(yⱼ/L_a)²) / Σ exp(−(yⱼ/L_a)²)
//! ```
//!
//! Four strands are checked here, and each has an independent source:
//!
//! - **The pattern.** The formula above is written out from theory in this
//!   file and never asked of the engine.
//! - **The index.** Its weights are the same Gaussian, so an anomaly moved off
//!   the equator must weigh exactly `exp((y₀² − y₁²)/L_a²)` less — a closed
//!   form, not a measured ratio.
//! - **The sign.** `CONTEXT.md`, *Bjerknes feedback*: warmer eastern SST
//!   weakens the trade winds. A warm anomaly must make `τx` less easterly, and
//!   the thermocline flatter, than the same run without the feedback.
//! - **The regression.** The ticket's acceptance criterion: at zero feedback
//!   strength the run is the uncoupled one, bit for bit.

use engine::forcing::StageForcing;
use engine::sst::{SstParams, DEFAULT_SURFACE_DRAG_PER_S};
use engine::wind_response::{
    CoupledWind, SstWindResponse, WindResponseError, WindResponseParams,
    ATMOSPHERIC_GRAVITY_WAVE_SPEED_M_PER_S, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M,
};
use engine::{
    Basin, BetaPlane, CompositeWind, Grid, OceanState, PhysicalParams, Scenario, ScenarioConfig,
    ScenarioError, Solver, Spacing, SteadyTradeWinds, WindForcing, WindStress, WindStressField,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s² (Gill, *Atmosphere–Ocean Dynamics*, ch. 11).
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H`, in metres.
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;
/// Seconds in a day.
const SECONDS_PER_DAY: f64 = 86_400.0;
/// Rayleigh damping `r`, in s⁻¹: a 100-day decay.
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0 / (100.0 * SECONDS_PER_DAY);
/// Meridional gradient of the Coriolis parameter, in m⁻¹s⁻¹ (`CONTEXT.md`).
const BETA_PER_M_PER_S: f64 = 2.3e-11;

/// Mixed-layer depth `H_m`, in metres (Zebiak & Cane, *Mon. Wea. Rev.* 115,
/// 1987, § 2b).
const MIXED_LAYER_DEPTH_M: f64 = 50.0;
/// Zonal gradient of the mean SST, in K/m: 6 K over the basin's 15 000 km,
/// negative because the ocean cools eastward.
const MEAN_ZONAL_SST_GRADIENT_K_PER_M: f64 = -4.0e-7;
/// Sensitivity `γ = ∂T_sub/∂h` of the entrained water, in K/m (Zebiak & Cane
/// 1987, § 2c).
const SUBSURFACE_SENSITIVITY_K_PER_M: f64 = 0.1;
/// Thermal damping `ε_T`, in s⁻¹: a 125-day relaxation (Zebiak & Cane 1987,
/// § 2b).
const THERMAL_DAMPING_PER_S: f64 = 1.0 / (125.0 * SECONDS_PER_DAY);

/// Equatorial trade-wind stress, in Pa. Easterly, so negative
/// (`CONTEXT.md`, *Wind stress*).
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// Feedback strength `μ`, in Pa/K, for the tests that want the loop closed.
///
/// A 1 K warming relaxing the trades by 0.02 Pa — 40% of
/// [`TRADE_WIND_STRESS_PA`] — which is the order the observed regression of
/// equatorial zonal stress on the Niño-3 index carries.
const FEEDBACK_STRENGTH_PA_PER_K: f64 = 0.02;

fn physical_params() -> PhysicalParams {
    PhysicalParams::new(
        REDUCED_GRAVITY_M_PER_S2,
        MEAN_THERMOCLINE_DEPTH_M,
        RAYLEIGH_DAMPING_PER_S,
        BETA_PER_M_PER_S,
        SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("these are the standard equatorial-Pacific parameters")
}

fn sst_params() -> SstParams {
    SstParams::new(
        MIXED_LAYER_DEPTH_M,
        DEFAULT_SURFACE_DRAG_PER_S,
        MEAN_ZONAL_SST_GRADIENT_K_PER_M,
        SUBSURFACE_SENSITIVITY_K_PER_M,
        THERMAL_DAMPING_PER_S,
    )
    .expect("these are the standard Zebiak-Cane mixed-layer parameters")
}

fn response_params(feedback_strength_pa_per_k: f64) -> WindResponseParams {
    WindResponseParams::new(
        feedback_strength_pa_per_k,
        DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M,
    )
    .expect("a non-negative strength and the default scale are a valid response")
}

/// A basin `ny` rows tall spanning `meridional_extent_m`, centred on the
/// equator, with `ny` odd so that one row of cell centers lies exactly on it.
fn equatorial_basin(nx: usize, ny: usize, dx_m: f64, meridional_extent_m: f64) -> Basin {
    assert!(
        ny % 2 == 1,
        "an odd row count puts a center row on the equator"
    );
    let grid = Grid::new(nx, ny).expect("a basin has cells");
    let spacing = Spacing::new(dx_m, meridional_extent_m / ny as f64).expect("cells have width");
    Basin::centered_on_equator(grid, spacing)
}

/// The index of the cell-center row that lies on the equator.
const fn equator_row(ny: usize) -> usize {
    (ny - 1) / 2
}

/// The test basin: 24 columns by 21 rows over 4000 km of latitude, the shape
/// the T-12.1 suite used for the same kind of question.
fn test_basin() -> Basin {
    equatorial_basin(24, 21, 2.0e5, 4.0e6)
}

/// The Gill meridional structure of the atmospheric response, written out from
/// theory: `exp(−(y/L_a)²)`, the equatorially trapped zonal wind of Gill
/// (1980).
fn analytic_meridional_structure(y_m: f64, meridional_scale_m: f64) -> f64 {
    let scaled = y_m / meridional_scale_m;
    (-scaled * scaled).exp()
}

// ---------------------------------------------------------------------------
// The pattern and the index
// ---------------------------------------------------------------------------

#[test]
fn the_default_meridional_scale_is_the_atmospheric_equatorial_rossby_radius() {
    // Gill (1980) § 2: the equatorially trapped zonal wind of the Kelvin part
    // of the response falls off as `exp(−βy²/(2c_a))`, which is
    // `exp(−(y/L_a)²)` with `L_a = √(2·c_a/β)`. The constant the engine
    // carries is that number, quoted to two significant figures — which is the
    // whole tolerance here.
    let derived_m = (2.0 * ATMOSPHERIC_GRAVITY_WAVE_SPEED_M_PER_S / BETA_PER_M_PER_S).sqrt();
    let relative_error = (DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M - derived_m).abs() / derived_m;
    assert!(
        relative_error < 2.0e-2,
        "the default meridional scale is {DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M} m, but \
         √(2·c_a/β) is {derived_m} m"
    );
}

#[test]
fn a_uniform_anomaly_gives_the_gill_gaussian_patch() {
    // With `T'` the same everywhere, the weighted index is that value exactly
    // — the weights cancel — so the response is the closed form above with
    // `⟨T'⟩ = T₀`, and nothing about the grid enters.
    let basin = test_basin();
    let anomaly_k = 1.5;
    let params = response_params(FEEDBACK_STRENGTH_PA_PER_K);
    let mut response = SstWindResponse::new(basin, params);

    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(anomaly_k);
    response.observe(state.sst_anomaly_k().expect("a coupled state carries `T'`"));

    // A weighted mean of one repeated value is that value, to the rounding of
    // summing `nx · ny` terms: 504 additions at a relative error of at most
    // one machine epsilon each is under `1e-13`.
    assert!(
        (response.index_k() - anomaly_k).abs() < 1e-13 * anomaly_k,
        "the index of a uniform {anomaly_k} K anomaly is {}",
        response.index_k()
    );

    for row in 0..basin.grid().ny() {
        let y_m = basin.y_of_row_m(engine::H_STAGGERING, row);
        let expected_pa = FEEDBACK_STRENGTH_PA_PER_K
            * anomaly_k
            * analytic_meridional_structure(y_m, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M);
        let (tau_x_pa, tau_y_pa) = response.stress(0.0, y_m, 0.0);
        // The engine evaluates the same closed form, so the two agree to a
        // handful of ulps; `1e-15` relative is that, with room for the
        // exponential's own rounding.
        assert!(
            (tau_x_pa - expected_pa).abs() <= 1e-15 * expected_pa.abs().max(1e-3),
            "row {row} at y = {y_m} m answers {tau_x_pa} Pa, not {expected_pa} Pa"
        );
        assert_eq!(
            tau_y_pa, 0.0,
            "the statistical response is zonal only, but row {row} carries a meridional stress"
        );
    }
}

#[test]
fn the_index_weights_the_equator_the_way_the_pattern_does() {
    // The same anomaly, in one row and then in another, must weigh in the
    // ratio of the two Gaussian weights — the closed form below, which is what
    // makes the index the projection of `T'` onto the atmosphere's equatorial
    // mode rather than a plain basin average.
    let basin = test_basin();
    let grid = basin.grid();
    let equator = equator_row(grid.ny());
    let offset_row = equator + 4;
    let params = response_params(FEEDBACK_STRENGTH_PA_PER_K);

    let index_of = |row: usize| {
        let mut state = OceanState::at_rest_with_sst_anomaly(grid);
        let anomaly = state
            .sst_anomaly_k_mut()
            .expect("a coupled state carries `T'`");
        for column in 0..grid.nx() {
            *anomaly
                .get_mut(column, row)
                .expect("the row is inside the basin") = 1.0;
        }
        let mut response = SstWindResponse::new(basin, params);
        response.observe(anomaly);
        response.index_k()
    };

    let y_equator_m = basin.y_of_row_m(engine::H_STAGGERING, equator);
    let y_offset_m = basin.y_of_row_m(engine::H_STAGGERING, offset_row);
    let expected_ratio =
        analytic_meridional_structure(y_offset_m, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M)
            / analytic_meridional_structure(y_equator_m, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M);

    let ratio = index_of(offset_row) / index_of(equator);
    // Both indices are the same sum with one weight changed, so the ratio is
    // exact bar rounding; `1e-12` relative is the accumulated error of two
    // 504-term sums.
    assert!(
        (ratio - expected_ratio).abs() < 1e-12 * expected_ratio,
        "an anomaly {y_offset_m} m off the equator weighs {ratio} of one on it, not \
         {expected_ratio}"
    );
}

#[test]
fn the_response_is_calm_while_the_ocean_is_at_its_climatology() {
    // `T'` is an anomaly, so a mixed layer sitting at its climatological
    // temperature must leave the prescribed trades exactly as written.
    let basin = test_basin();
    let mut response = SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K));
    let state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    response.observe(state.sst_anomaly_k().expect("a coupled state carries `T'`"));

    assert_eq!(response.index_k(), 0.0);
    assert_eq!(response.stress(0.0, 0.0, 0.0), (0.0, 0.0));
}

#[test]
fn the_response_composes_with_the_prescribed_trades_by_addition() {
    // The T-03.3 pattern: forcings are superimposed, so the stress the ocean
    // feels is the pointwise sum of the trades and the atmosphere's answer.
    let basin = test_basin();
    let anomaly_k = 2.0;
    let params = response_params(FEEDBACK_STRENGTH_PA_PER_K);
    let mut response = SstWindResponse::new(basin, params);
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(anomaly_k);
    response.observe(state.sst_anomaly_k().expect("a coupled state carries `T'`"));

    let trades = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress");
    let mut coupled = CoupledWind::new(
        basin,
        CompositeWind::new().with(trades),
        SstWindResponse::new(basin, params),
    );
    let field = coupled.at(0.0, &state);

    for row in 0..field.tau_x_pa().ny() {
        let y_m = basin.y_of_row_m(engine::U_STAGGERING, row);
        let expected_pa = TRADE_WIND_STRESS_PA + response.stress(0.0, y_m, 0.0).0;
        for column in 0..field.tau_x_pa().nx() {
            let actual_pa = *field
                .tau_x_pa()
                .get(column, row)
                .expect("the loop bounds are the field's own");
            // The two are the same two numbers added in the same order, so
            // this is exact.
            assert_eq!(
                actual_pa, expected_pa,
                "the coupled stress at column {column}, row {row} is not the sum"
            );
        }
    }
}

#[test]
fn a_warm_anomaly_weakens_the_alizes() {
    // `CONTEXT.md`, *Bjerknes feedback*: warmer eastern SST → weaker trade
    // winds. Easterly is `τx < 0`, so "weaker" is "less negative".
    let basin = test_basin();
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(1.0);

    let mut coupled = CoupledWind::new(
        basin,
        CompositeWind::new()
            .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress")),
        SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K)),
    );
    let equator = equator_row(basin.grid().ny());
    let warm_pa = *coupled
        .at(0.0, &state)
        .tau_x_pa()
        .get(0, equator)
        .expect("the equator row is inside the basin");
    assert!(
        warm_pa > TRADE_WIND_STRESS_PA && warm_pa < 0.0,
        "a 1 K warm anomaly left the equatorial stress at {warm_pa} Pa, not between \
         {TRADE_WIND_STRESS_PA} Pa and calm"
    );

    // And the other way round: a cold anomaly strengthens them.
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(-1.0);
    let cold_pa = *coupled
        .at(0.0, &state)
        .tau_x_pa()
        .get(0, equator)
        .expect("the equator row is inside the basin");
    assert!(
        cold_pa < TRADE_WIND_STRESS_PA,
        "a 1 K cold anomaly left the equatorial stress at {cold_pa} Pa, which is no stronger \
         than the prescribed {TRADE_WIND_STRESS_PA} Pa"
    );
}

// ---------------------------------------------------------------------------
// The acceptance criterion: zero feedback is the model that was validated
// ---------------------------------------------------------------------------

/// Step a coupled basin `steps` times under the alizés and return its four
/// fields. `feedback_strength_pa_per_k` of `None` builds the T-12.1 run — the
/// prescribed forcing, with no response object anywhere in the loop.
fn stepped_fields(
    feedback_strength_pa_per_k: Option<f64>,
    steps: u64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let params = physical_params();
    let basin = test_basin();
    let plane = BetaPlane::of_basin(params, basin);
    let dt_s = 1800.0;
    let mut solver = Solver::coupled_to_sst(
        basin.grid(),
        basin.spacing(),
        params,
        plane,
        dt_s,
        sst_params(),
    )
    .expect("half an hour is well inside both timestep bounds");
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    // A thermocline anomaly to give the coupling something to feed on, exactly
    // as T-12.1's regression does.
    state.h_mut().as_mut_slice().fill(5.0);
    let wind = || {
        CompositeWind::new()
            .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress"))
    };

    match feedback_strength_pa_per_k {
        None => {
            let mut forcing = WindForcing::new(basin, wind());
            for step in 0..steps {
                solver.step_with_forcing(&mut state, step as f64 * dt_s, &mut forcing);
            }
        }
        Some(strength) => {
            let mut forcing = CoupledWind::new(
                basin,
                wind(),
                SstWindResponse::new(basin, response_params(strength)),
            );
            for step in 0..steps {
                solver.step_with_forcing(&mut state, step as f64 * dt_s, &mut forcing);
            }
        }
    }

    (
        state.h().as_slice().to_vec(),
        state.u().as_slice().to_vec(),
        state.v().as_slice().to_vec(),
        state
            .sst_anomaly_k()
            .expect("a coupled state carries `T'`")
            .as_slice()
            .to_vec(),
    )
}

#[test]
fn zero_feedback_strength_is_the_prescribed_wind_run_bit_for_bit() {
    // The ticket's acceptance criterion, in the same form T-12.1 asserted its
    // own: not "close", but identical. Anything else would mean the response
    // had re-associated the arithmetic of a model that is already validated.
    let prescribed = stepped_fields(None, 200);
    let unfed_back = stepped_fields(Some(0.0), 200);
    assert_eq!(prescribed.0, unfed_back.0, "the thermocline `h` differs");
    assert_eq!(prescribed.1, unfed_back.1, "the zonal current `u` differs");
    assert_eq!(
        prescribed.2, unfed_back.2,
        "the meridional current `v` differs"
    );
    assert_eq!(prescribed.3, unfed_back.3, "the SST anomaly `T'` differs");
}

#[test]
fn a_positive_feedback_strength_actually_changes_the_run() {
    // The companion the test above needs: "nothing changed" is only worth
    // asserting if the machinery was capable of changing something.
    let unfed_back = stepped_fields(Some(0.0), 200);
    let fed_back = stepped_fields(Some(FEEDBACK_STRENGTH_PA_PER_K), 200);
    assert_ne!(
        unfed_back.0, fed_back.0,
        "closing the loop left the thermocline untouched"
    );
}

#[test]
fn the_closed_loop_flattens_the_thermocline_a_warm_anomaly_made() {
    // The Bjerknes chain end to end, in the direction `CONTEXT.md` states it:
    // a warm mixed layer weakens the alizés, weaker alizés pile less water in
    // the west, and the thermocline tilt relaxes. The zonal tilt is what
    // `steady_wind_tilt.rs` validated for the prescribed case, so the check is
    // that closing the loop reduces it.
    let params = physical_params();
    let basin = test_basin();
    let grid = basin.grid();
    let plane = BetaPlane::of_basin(params, basin);
    let dt_s = 1800.0;
    let steps = 400;
    let equator = equator_row(grid.ny());

    let tilt_m = |feedback_strength_pa_per_k: f64| {
        let mut solver =
            Solver::coupled_to_sst(grid, basin.spacing(), params, plane, dt_s, sst_params())
                .expect("half an hour is well inside both timestep bounds");
        let mut state = OceanState::at_rest_with_sst_anomaly(grid);
        // A uniformly warm mixed layer, so the atmosphere has something to
        // answer from the first step.
        state
            .sst_anomaly_k_mut()
            .expect("a coupled state carries `T'`")
            .as_mut_slice()
            .fill(1.0);
        let mut forcing = CoupledWind::new(
            basin,
            CompositeWind::new()
                .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress")),
            SstWindResponse::new(basin, response_params(feedback_strength_pa_per_k)),
        );
        for step in 0..steps {
            solver.step_with_forcing(&mut state, step as f64 * dt_s, &mut forcing);
        }
        let west_m = *state
            .h()
            .get(0, equator)
            .expect("the western column is inside the basin");
        let east_m = *state
            .h()
            .get(grid.nx() - 1, equator)
            .expect("the eastern column is inside the basin");
        west_m - east_m
    };

    let prescribed_m = tilt_m(0.0);
    let coupled_m = tilt_m(FEEDBACK_STRENGTH_PA_PER_K);
    assert!(
        prescribed_m > 0.0,
        "the alizés should tilt the thermocline down to the west, but the tilt is \
         {prescribed_m} m"
    );
    assert!(
        coupled_m < prescribed_m,
        "closing the loop on a warm anomaly left the tilt at {coupled_m} m, no flatter than the \
         prescribed {prescribed_m} m"
    );
}

// ---------------------------------------------------------------------------
// The config parameter
// ---------------------------------------------------------------------------

#[test]
fn an_unphysical_response_is_refused_by_name() {
    assert_eq!(
        WindResponseParams::new(-1.0, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M),
        Err(WindResponseError::Negative {
            parameter: "wind_feedback_strength_pa_per_k",
            value: -1.0,
        })
    );
    assert_eq!(
        WindResponseParams::new(0.02, 0.0),
        Err(WindResponseError::NotPositive {
            parameter: "wind_response_meridional_scale_m",
            value: 0.0,
        })
    );
}

/// A scenario file with an optional `[sst]` section appended.
fn scenario_toml(sst_section: &str) -> String {
    format!(
        "[basin]\n\
         resolution_deg = 2.0\n\
         \n\
         [physics]\n\
         reduced_gravity_m_per_s2 = 0.05\n\
         mean_thermocline_depth_m = 150.0\n\
         rayleigh_damping_per_s = 1e-7\n\
         \n\
         [run]\n\
         dt_s = 1800.0\n\
         total_steps = 10\n\
         output_every_n_steps = 5\n\
         {sst_section}"
    )
}

/// The `[sst]` section of a coupled scenario, with `extra` lines appended.
fn sst_section(extra: &str) -> String {
    format!(
        "\n[sst]\n\
         mixed_layer_depth_m = 50.0\n\
         mean_zonal_sst_gradient_k_per_m = -4e-7\n\
         subsurface_temperature_sensitivity_k_per_m = 0.1\n\
         thermal_damping_per_s = 9.26e-8\n\
         {extra}"
    )
}

#[test]
fn a_coupled_scenario_that_says_nothing_about_the_wind_has_no_feedback() {
    // The T-12.1 scenario, unchanged: the section that switches the SST
    // equation on does not by itself close the loop, so every file written
    // before this ticket still describes the run it described.
    let scenario = Scenario::from_toml(&scenario_toml(&sst_section(""))).expect("a valid scenario");
    let response = scenario
        .wind_response_params()
        .expect("a coupled scenario carries a response");
    assert_eq!(response.feedback_strength_pa_per_k(), 0.0);
    assert_eq!(
        response.meridional_scale_m(),
        DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M
    );
}

#[test]
fn an_uncoupled_scenario_has_no_wind_response_at_all() {
    let scenario = Scenario::from_toml(&scenario_toml("")).expect("a valid scenario");
    assert!(scenario.wind_response_params().is_none());
}

#[test]
fn the_feedback_strength_is_a_config_parameter_and_round_trips_through_toml() {
    let source = scenario_toml(&sst_section(
        "wind_feedback_strength_pa_per_k = 0.02\n\
         wind_response_meridional_scale_m = 1.5e6\n",
    ));
    let config = ScenarioConfig::from_toml(&source).expect("a valid scenario");
    let scenario = config.build().expect("a runnable scenario");
    let response = scenario
        .wind_response_params()
        .expect("a coupled scenario carries a response");
    assert_eq!(response.feedback_strength_pa_per_k(), 0.02);
    assert_eq!(response.meridional_scale_m(), 1.5e6);

    let reparsed =
        ScenarioConfig::from_toml(&config.to_toml().expect("TOML can hold these numbers"))
            .expect("what the engine writes, the engine reads");
    assert_eq!(reparsed, config);
}

#[test]
fn a_negative_feedback_strength_is_refused_by_the_scenario_loader() {
    let source = scenario_toml(&sst_section("wind_feedback_strength_pa_per_k = -0.02\n"));
    let error = ScenarioConfig::from_toml(&source)
        .expect("the file parses")
        .build()
        .expect_err("an anti-Bjerknes atmosphere is not a scenario");
    assert!(
        matches!(error, ScenarioError::WindResponse(_)),
        "the loader reported {error:?} rather than a wind-response error"
    );
    assert!(
        error
            .to_string()
            .contains("wind_feedback_strength_pa_per_k"),
        "the message does not name the offending parameter: {error}"
    );
}

#[test]
fn a_feedback_asked_for_without_the_sst_coupling_is_a_field_the_format_rejects() {
    // The response reads `T'`, which only an `[sst]` scenario integrates, so
    // the parameters live in that section and nowhere else. A file that tries
    // to put them elsewhere is a file with a key the format does not define.
    let source = scenario_toml("\nwind_feedback_strength_pa_per_k = 0.02\n");
    assert!(ScenarioConfig::from_toml(&source).is_err());
}

// ---------------------------------------------------------------------------
// The whole run
// ---------------------------------------------------------------------------

#[test]
fn a_run_records_the_stress_the_ocean_actually_felt() {
    // The frames of a coupled run carry the total stress — the trades plus the
    // atmosphere's answer — because that is the forcing the step read. A run
    // whose frames showed the prescribed trades alone would be describing a
    // scenario that was not integrated.
    let basin = test_basin();
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(1.0);
    let mut forcing = CoupledWind::new(
        basin,
        CompositeWind::new()
            .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress")),
        SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K)),
    );
    let equator = equator_row(basin.grid().ny());
    let felt_pa = *StageForcing::at(&mut forcing, 0.0, &state)
        .tau_x_pa()
        .get(0, equator)
        .expect("the equator row is inside the basin");
    assert!(felt_pa > TRADE_WIND_STRESS_PA);

    // And a calm-anomaly stage sees the prescribed trades exactly, so the
    // total field is not carrying a stale answer from the warm one.
    let at_rest = OceanState::at_rest_with_sst_anomaly(basin.grid());
    let calm_pa = *StageForcing::at(&mut forcing, 0.0, &at_rest)
        .tau_x_pa()
        .get(0, equator)
        .expect("the equator row is inside the basin");
    assert_eq!(calm_pa, TRADE_WIND_STRESS_PA);
}

#[test]
fn the_stress_field_a_coupled_forcing_writes_covers_the_basin_it_was_built_for() {
    let basin = test_basin();
    let forcing = CoupledWind::new(
        basin,
        CompositeWind::new()
            .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress")),
        SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K)),
    );
    assert_eq!(forcing.basin().grid(), basin.grid());
}

#[test]
fn a_response_built_for_another_basin_is_a_bug_and_panics() {
    let basin = test_basin();
    let mut response = SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K));
    let elsewhere =
        OceanState::at_rest_with_sst_anomaly(Grid::new(8, 5).expect("a basin has cells"));
    let anomaly = elsewhere
        .sst_anomaly_k()
        .expect("a coupled state carries `T'`");
    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        response.observe(anomaly);
    }));
    assert!(
        panicked.is_err(),
        "observing an anomaly over another basin was accepted"
    );
}

#[test]
fn a_field_of_the_wrong_shape_never_reaches_the_index() {
    // The companion of the panic above: the field a coupled state carries is
    // the one the response was built for, so the ordinary path does not panic.
    let basin = test_basin();
    let mut response = SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K));
    let state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    response.observe(state.sst_anomaly_k().expect("a coupled state carries `T'`"));
    assert_eq!(response.index_k(), 0.0);
}

#[test]
fn the_response_is_a_wind_stress_like_any_other() {
    // The deliverable: it implements `WindStress`, so it stacks in a
    // `CompositeWind` beside the trades exactly as the burst of T-03.3 does.
    let basin = test_basin();
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries `T'`")
        .as_mut_slice()
        .fill(1.0);
    let mut response = SstWindResponse::new(basin, response_params(FEEDBACK_STRENGTH_PA_PER_K));
    response.observe(state.sst_anomaly_k().expect("a coupled state carries `T'`"));
    let alone = response.stress(0.0, 0.0, 0.0).0;

    let stacked = CompositeWind::new()
        .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress"))
        .with(response);
    let field = WindStressField::sampled(basin, &stacked, 0.0);
    let equator = equator_row(basin.grid().ny());
    assert_eq!(
        *field
            .tau_x_pa()
            .get(0, equator)
            .expect("the equator row is inside the basin"),
        TRADE_WIND_STRESS_PA + alone
    );
}
