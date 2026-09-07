//! Acceptance tests for T-08.5 — the engine's simulation API without a
//! filesystem.
//!
//! ADR-0012 puts the engine in the browser: the visualizer links it and
//! computes runs itself, so everything the browser needs — building a
//! `Scenario`, building a `Solver`, taking a step — has to be reachable with
//! the `fs` feature off. This file walks exactly that path and names nothing
//! the feature gates: no `std::fs`, no `std::path`, no `Scenario::load`, no
//! `RunWriter::create`, no `run_scenario`.
//!
//! It is a compile-time claim as much as a runtime one. What proves the gating
//! itself is `cargo build -p engine --no-default-features --target
//! wasm32-unknown-unknown`, which CI runs; what this file adds is that the
//! portable half still computes what it computed before it was a separately
//! compilable half — T-08.5 changes what is compiled, not what is computed.
//!
//! Nothing below is measured out of a run. The scenario's values are read off
//! the TOML text by eye, and the state after one step is the closed-form first
//! step of the momentum equation of ADR-0003.

use engine::{
    BetaPlane, OceanState, Scenario, Solver, WindForcing, WindStress,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3, U_STAGGERING,
};

/// A scenario small enough to step in a test, held as text rather than in a
/// file — the shape the browser has (ADR-0012).
///
/// A 40° × 10° box at 1°, which is 40 × 10 cells, under a steady easterly.
const SCENARIO_TOML: &str = r#"
[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -5.0
northern_latitude_deg = 5.0
resolution_deg = 1.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 4
output_every_n_steps = 2

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0
"#;

// The values of `SCENARIO_TOML`, read off it by eye.
const RESOLUTION_DEG: f64 = 1.0;
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.06;
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0e-7;
const DT_S: f64 = 3600.0;
const TOTAL_STEPS: u64 = 4;
const EQUATORIAL_ZONAL_STRESS_PA: f64 = -0.05;
const MERIDIONAL_DECAY_SCALE_M: f64 = 361_000.0;
/// 40° of longitude and 10° of latitude at 1° a cell.
const EXPECTED_NX: usize = 40;
const EXPECTED_NY: usize = 10;
/// A frame every two steps over four steps is the initial state and two more.
const EXPECTED_FRAME_COUNT: u64 = 3;
/// Two steps of an hour.
const EXPECTED_INTERVAL_S: f64 = 7200.0;

/// Relative tolerance on the current after one step, assembled from the terms
/// the closed form below drops.
///
/// From rest under a stress that does not vary zonally, the correction to
/// `u = X·dt` over the first step is the Rayleigh damping the wind works
/// against: integrating `du/dt = X − r·u` from zero gives
/// `u = (X/r)·(1 − e^{−r·dt}) = X·dt·(1 − r·dt/2 + …)`, whose leading term is
/// `r·dt/2 = 1.8e-4`. The Coriolis term enters at `(f·dt)²/2`, which for
/// `f = β·y` half a degree off the equator is `1e-5`. Twice the larger of the
/// two leaves room for the RK4 truncation of both without admitting a third
/// term.
const FIRST_STEP_RELATIVE_TOLERANCE: f64 = 1.0e-3;

/// A `Scenario` is constructible from TOML **text**, with no path anywhere.
#[test]
fn a_scenario_is_built_from_toml_text() {
    let scenario = Scenario::from_toml(SCENARIO_TOML).expect("the scenario text is valid");

    let bounds = scenario.bounds();
    assert_eq!(bounds.resolution_deg(), RESOLUTION_DEG);
    assert_eq!(scenario.basin().grid().nx(), EXPECTED_NX);
    assert_eq!(scenario.basin().grid().ny(), EXPECTED_NY);

    let params = scenario.physical_params();
    assert_eq!(params.reduced_gravity_m_per_s2(), REDUCED_GRAVITY_M_PER_S2);
    assert_eq!(params.mean_thermocline_depth_m(), MEAN_THERMOCLINE_DEPTH_M);
    assert_eq!(params.rayleigh_damping_per_s(), RAYLEIGH_DAMPING_PER_S);

    let schedule = scenario.output_schedule();
    assert_eq!(schedule.dt_s(), DT_S);
    assert_eq!(schedule.total_steps(), TOTAL_STEPS);
    assert_eq!(schedule.frame_count(), EXPECTED_FRAME_COUNT);
    assert_eq!(schedule.interval_s(), EXPECTED_INTERVAL_S);
}

/// The wind a scenario carries is sampled without a filesystem too — the
/// stress a browser frame records is the stress its step read (ADR-0012).
#[test]
fn the_scenario_wind_is_sampled_from_text_alone() {
    let scenario = Scenario::from_toml(SCENARIO_TOML).expect("the scenario text is valid");

    // On the equator the profile's Gaussian is one, so `τx = τ₀` exactly and
    // `τy` is zero: `τx(x, y, t) = τ₀·exp(−(y/Ly)²)`, `τy = 0`.
    let (zonal_pa, meridional_pa) = scenario.wind().stress(0.0, 0.0, 0.0);
    assert_eq!(zonal_pa, EQUATORIAL_ZONAL_STRESS_PA);
    assert_eq!(meridional_pa, 0.0);
}

/// The solver that scenario implies takes a step, and the step is the one the
/// equations say it is.
///
/// The browser's loop is this loop (ADR-0012, *Decision*): hold a `Scenario`,
/// build a `Solver`, step it. What is asserted is the closed-form first step
/// of `∂u/∂t = −g'·∂h/∂x + f·v − r·u + τx/(ρ₀·H)` from rest. The stress of
/// `SCENARIO_TOML` does not vary with longitude, so away from the walls the
/// pressure gradient is zero at the start of the step and stays zero through
/// it, and `v` starts at zero: the first step is the wind's impulse,
/// `u = X·dt` with `X = τx/(ρ₀·H)`, less the damping named in
/// `FIRST_STEP_RELATIVE_TOLERANCE`.
#[test]
fn a_solver_steps_a_scenario_that_never_saw_a_file() {
    let scenario = Scenario::from_toml(SCENARIO_TOML).expect("the scenario text is valid");
    let basin = scenario.basin();
    let params = scenario.physical_params();
    let schedule = scenario.output_schedule();

    let mut solver = Solver::new(
        basin.grid(),
        basin.spacing(),
        params,
        BetaPlane::of_basin(params, basin),
        schedule.dt_s(),
    )
    .expect("the scenario's timestep is inside the CFL bound");
    let mut state = OceanState::at_rest(basin.grid());
    let mut forcing = WindForcing::new(basin, scenario.wind());

    solver.step_with_forcing(&mut state, 0.0, &mut forcing);

    // Mid-basin, so the walls' no-normal-flow condition is many cells away,
    // and one row off the equator, which is the closest a 1° grid's cell
    // centers come to it.
    let column = EXPECTED_NX / 2;
    let row = EXPECTED_NY / 2;
    let scaled = basin.y_of_row_m(U_STAGGERING, row) / MERIDIONAL_DECAY_SCALE_M;
    let stress_pa = EQUATORIAL_ZONAL_STRESS_PA * (-scaled * scaled).exp();
    let expected_u_m_per_s =
        stress_pa * DT_S / (SEAWATER_REFERENCE_DENSITY_KG_PER_M3 * MEAN_THERMOCLINE_DEPTH_M);

    let u_m_per_s = *state
        .u()
        .get(column, row)
        .expect("the sampled cell is inside the basin");
    let relative_error = (u_m_per_s - expected_u_m_per_s).abs() / expected_u_m_per_s.abs();
    assert!(
        relative_error < FIRST_STEP_RELATIVE_TOLERANCE,
        "u after one step is {u_m_per_s} m/s, and the wind's impulse is \
         {expected_u_m_per_s} m/s: a relative gap of {relative_error}"
    );
}
