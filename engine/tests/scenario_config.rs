//! T-03.4 — a scenario is a TOML file, and this is what it must mean.
//!
//! Two acceptance criteria, and the tests below are grouped by them:
//!
//! 1. each of the three example configs loads and produces *the corresponding*
//!    [`WindStress`] implementation with the right parameters, and
//! 2. an invalid config — a bad grid size, an unknown forcing type — fails
//!    with a clear error rather than a panic.
//!
//! Expected values come from two independent sources, never from running the
//! loader: the parameters are read off the example files by eye and written
//! here as literals, and the stresses those parameters imply are computed from
//! the closed-form expressions in `CONTEXT.md` and
//! `docs/planning/01-scientific-model.md`.

use std::path::{Path, PathBuf};

use engine::scenario::{Scenario, ScenarioError, ScenarioWind};
use engine::{ScenarioConfig, WindStress, TROPICAL_YEAR_S};

/// Relative tolerance for a stress compared against a closed-form expression
/// evaluated here. Both sides are the same handful of `f64` multiplications
/// and one `exp`, differing only in the order they are written, so the gap is
/// a few units in the last place; 1e-12 is roughly 4500 ulp at these
/// magnitudes and still four orders of magnitude tighter than any physically
/// meaningful difference.
const STRESS_RELATIVE_TOLERANCE: f64 = 1e-12;

/// Parameters shared by all three example files, read off them by eye.
const NX: usize = 200;
const NY: usize = 60;
const DX_M: f64 = 50_000.0;
const DY_M: f64 = 50_000.0;
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.06;
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0e-7;
const DT_S: f64 = 3600.0;
const TOTAL_STEPS: u64 = 17_520;
const OUTPUT_EVERY_N_STEPS: u64 = 24;
const TRADE_STRESS_PA: f64 = -0.05;
const TRADE_DECAY_SCALE_M: f64 = 361_000.0;

/// The example scenarios the ticket asks for, one per scenario of
/// `docs/planning/01-scientific-model.md`.
const EXAMPLE_FILE_NAMES: [&str; 3] = [
    "steady-trades.toml",
    "seasonal-cycle.toml",
    "wind-burst.toml",
];

fn scenarios_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("scenarios")
}

fn load(file_name: &str) -> Scenario {
    Scenario::load(&scenarios_dir().join(file_name))
        .unwrap_or_else(|error| panic!("{file_name} should load: {error}"))
}

fn assert_close(actual: f64, expected: f64, what: &str) {
    let tolerance = expected.abs() * STRESS_RELATIVE_TOLERANCE;
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: expected {expected}, got {actual} (tolerance {tolerance})"
    );
}

/// The Gaussian of `forcing`: `exp(−(offset/scale)²)`, which is `1/e` at
/// `offset = scale`. Written out here so the tests do not borrow the
/// implementation's helper.
fn gaussian(offset: f64, scale: f64) -> f64 {
    (-(offset / scale) * (offset / scale)).exp()
}

fn assert_shared_basin_physics_and_run(scenario: &Scenario) {
    let basin = scenario.basin();
    assert_eq!(basin.grid().nx(), NX);
    assert_eq!(basin.grid().ny(), NY);
    assert_eq!(basin.spacing().dx_m(), DX_M);
    assert_eq!(basin.spacing().dy_m(), DY_M);
    // No `western_edge_x_m` or `southern_edge_y_m` in any example file, so the
    // basin is the one centred on the equator with its west wall at x = 0.
    assert_eq!(basin.western_edge_x_m(), 0.0);
    assert_eq!(basin.southern_edge_y_m(), -(NY as f64 * DY_M) / 2.0);

    let physics = scenario.physical_params();
    assert_eq!(physics.reduced_gravity_m_per_s2(), REDUCED_GRAVITY_M_PER_S2);
    assert_eq!(physics.mean_thermocline_depth_m(), MEAN_THERMOCLINE_DEPTH_M);
    assert_eq!(physics.rayleigh_damping_per_s(), RAYLEIGH_DAMPING_PER_S);
    // β and ρ₀ are absent from the example files, so they must come from the
    // named constants of `params`.
    assert_eq!(
        physics.beta_per_m_per_s(),
        engine::EQUATORIAL_BETA_PER_M_PER_S
    );
    assert_eq!(
        physics.reference_density_kg_per_m3(),
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3
    );

    let schedule = scenario.output_schedule();
    assert_eq!(schedule.dt_s(), DT_S);
    assert_eq!(schedule.total_steps(), TOTAL_STEPS);
    assert_eq!(
        schedule.interval_s(),
        OUTPUT_EVERY_N_STEPS as f64 * DT_S,
        "a frame every {OUTPUT_EVERY_N_STEPS} steps"
    );
}

// ---------------------------------------------------------------------------
// Criterion 1: each example config loads and produces the corresponding
// `WindStress` implementation with the right parameters.
// ---------------------------------------------------------------------------

#[test]
fn every_example_scenario_loads() {
    for file_name in EXAMPLE_FILE_NAMES {
        let path = scenarios_dir().join(file_name);
        if let Err(error) = Scenario::load(&path) {
            panic!("{file_name} should load, got: {error}");
        }
    }
}

#[test]
fn steady_trades_example_builds_steady_trade_winds() {
    let scenario = load("steady-trades.toml");
    assert_shared_basin_physics_and_run(&scenario);

    let [ScenarioWind::Steady(trades)] = scenario.winds() else {
        panic!(
            "the steady-trades scenario should carry exactly one SteadyTradeWinds, got {:?}",
            scenario.winds()
        );
    };
    assert_eq!(trades.equatorial_zonal_stress_pa(), TRADE_STRESS_PA);
    assert_eq!(
        trades.meridional_decay_scale_m(),
        Some(TRADE_DECAY_SCALE_M),
        "the file states a decay scale, so it must not become the structureless field"
    );

    // τx(x, y, t) = τ₀·exp(−(y/Ly)²), independent of x and t
    // (docs/planning/01-scientific-model.md, *Control scenario*). One decay
    // scale off the equator the stress is τ₀/e.
    let wind = scenario.wind();
    let (tau_x_pa, tau_y_pa) = wind.stress(0.0, TRADE_DECAY_SCALE_M, 0.0);
    assert_close(
        tau_x_pa,
        TRADE_STRESS_PA * std::f64::consts::E.recip(),
        "τx one decay scale north of the equator",
    );
    assert_eq!(tau_y_pa, 0.0, "the alizés are zonal");
}

#[test]
fn seasonal_cycle_example_builds_seasonal_trade_winds() {
    let scenario = load("seasonal-cycle.toml");
    assert_shared_basin_physics_and_run(&scenario);

    // Read off `scenarios/seasonal-cycle.toml`.
    const RELATIVE_AMPLITUDE: f64 = 0.2;
    const PEAK_TIME_S: f64 = 18_144_000.0;

    let [ScenarioWind::Seasonal(seasonal)] = scenario.winds() else {
        panic!(
            "the seasonal scenario should carry exactly one SeasonalTradeWinds, got {:?}",
            scenario.winds()
        );
    };
    assert_eq!(seasonal.relative_amplitude(), RELATIVE_AMPLITUDE);
    assert_eq!(seasonal.peak_time_s(), PEAK_TIME_S);
    assert_eq!(
        seasonal.steady().equatorial_zonal_stress_pa(),
        TRADE_STRESS_PA,
        "the season modulates the steady field the same file describes"
    );
    assert_eq!(
        seasonal.steady().meridional_decay_scale_m(),
        Some(TRADE_DECAY_SCALE_M)
    );

    // τx = τ₀·exp(−(y/Ly)²)·(1 + a·cos(2π(t − t_peak)/T_year)) (CONTEXT.md,
    // *Seasonal cycle*): strongest by a factor (1 + a) at the peak, weakest by
    // (1 − a) half a year later.
    let wind = scenario.wind();
    let (strongest_pa, _) = wind.stress(0.0, 0.0, PEAK_TIME_S);
    assert_close(
        strongest_pa,
        TRADE_STRESS_PA * (1.0 + RELATIVE_AMPLITUDE),
        "τx on the equator at the seasonal peak",
    );
    let (weakest_pa, _) = wind.stress(0.0, 0.0, PEAK_TIME_S + TROPICAL_YEAR_S / 2.0);
    assert_close(
        weakest_pa,
        TRADE_STRESS_PA * (1.0 - RELATIVE_AMPLITUDE),
        "τx on the equator half a year after the seasonal peak",
    );
}

#[test]
fn wind_burst_example_stacks_a_burst_on_the_trades() {
    let scenario = load("wind-burst.toml");
    assert_shared_basin_physics_and_run(&scenario);

    // Read off `scenarios/wind-burst.toml`.
    const PEAK_STRESS_PA: f64 = 0.04;
    const CENTER_X_M: f64 = 2_000_000.0;
    const ZONAL_SCALE_M: f64 = 1_000_000.0;
    const MERIDIONAL_SCALE_M: f64 = 361_000.0;
    const PEAK_TIME_S: f64 = 31_556_926.08;
    const DURATION_S: f64 = 864_000.0;

    let [ScenarioWind::Steady(trades), ScenarioWind::Burst(burst)] = scenario.winds() else {
        panic!(
            "the burst scenario should carry trades then a burst, in that order, got {:?}",
            scenario.winds()
        );
    };
    assert_eq!(trades.equatorial_zonal_stress_pa(), TRADE_STRESS_PA);
    assert_eq!(burst.peak_zonal_stress_pa(), PEAK_STRESS_PA);
    assert_eq!(burst.center_x_m(), CENTER_X_M);
    assert_eq!(burst.zonal_scale_m(), ZONAL_SCALE_M);
    assert_eq!(burst.meridional_scale_m(), MERIDIONAL_SCALE_M);
    assert_eq!(burst.peak_time_s(), PEAK_TIME_S);
    assert_eq!(burst.duration_s(), DURATION_S);

    // The equations are linear in the stress, so the two components add
    // (CONTEXT.md, *Westerly wind burst*): at the burst's centre in space and
    // time the total is τ₀ + τ_burst, still easterly for these numbers.
    let wind = scenario.wind();
    let (peak_total_pa, _) = wind.stress(CENTER_X_M, 0.0, PEAK_TIME_S);
    assert_close(
        peak_total_pa,
        TRADE_STRESS_PA + PEAK_STRESS_PA,
        "τx at the centre of the burst",
    );
    assert!(
        peak_total_pa < 0.0,
        "a 0.04 Pa burst does not reverse 0.05 Pa of trades"
    );

    // One e-folding away on each of the burst's three axes the anomaly is
    // τ_burst/e³, while the trades have only their own meridional decay.
    let x_m = CENTER_X_M + ZONAL_SCALE_M;
    let y_m = MERIDIONAL_SCALE_M;
    let t_s = PEAK_TIME_S + DURATION_S;
    let (off_peak_pa, _) = wind.stress(x_m, y_m, t_s);
    let expected_pa = TRADE_STRESS_PA * gaussian(y_m, TRADE_DECAY_SCALE_M)
        + PEAK_STRESS_PA * std::f64::consts::E.powi(-3);
    assert_close(off_peak_pa, expected_pa, "τx one e-folding off the burst");
}

#[test]
fn every_example_round_trips_through_toml() {
    // serde-based (de)serialization, not just deserialization: what the loader
    // read has to write back out as the same scenario, so a run can record the
    // scenario that produced it.
    for file_name in EXAMPLE_FILE_NAMES {
        let source = std::fs::read_to_string(scenarios_dir().join(file_name))
            .expect("the example file is in the repository");
        let config = ScenarioConfig::from_toml(&source)
            .unwrap_or_else(|error| panic!("{file_name}: {error}"));
        let written = config.to_toml().expect("a valid scenario serializes");
        let reparsed = ScenarioConfig::from_toml(&written)
            .unwrap_or_else(|error| panic!("{file_name}, written back: {error}"));
        assert_eq!(reparsed, config, "{file_name} should survive a round trip");
    }
}

// ---------------------------------------------------------------------------
// Criterion 2: invalid configs fail with a clear error, not a panic.
// ---------------------------------------------------------------------------

/// A valid config, as a template the failure cases mutate one line of.
const VALID_TOML: &str = r#"
[basin]
nx = 200
ny = 60
dx_m = 50000.0
dy_m = 50000.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 17520
output_every_n_steps = 24

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
"#;

fn error_from(toml: &str) -> ScenarioError {
    Scenario::from_toml(toml).expect_err("this scenario should be rejected")
}

#[test]
fn the_template_the_failure_cases_mutate_is_itself_valid() {
    // Otherwise a case below could pass for the wrong reason.
    assert!(Scenario::from_toml(VALID_TOML).is_ok());
}

#[test]
fn a_zero_grid_extent_is_rejected_by_name() {
    let error = error_from(&VALID_TOML.replace("ny = 60", "ny = 0"));
    let message = error.to_string();
    assert!(
        message.contains("ny") && message.contains('0'),
        "the message should name the offending extent, got: {message}"
    );
}

#[test]
fn a_negative_grid_extent_is_rejected_at_the_line_that_carries_it() {
    // A count of cells has no negative values to reject downstream, so this
    // one is refused while the file is still being parsed. The message still
    // has to point at the offending line rather than at the file.
    let error = error_from(&VALID_TOML.replace("nx = 200", "nx = -200"));
    let message = error.to_string();
    assert!(
        message.contains("nx = -200"),
        "the message should quote the line it rejected, got: {message}"
    );
}

#[test]
fn a_negative_cell_spacing_is_rejected_by_name() {
    let error = error_from(&VALID_TOML.replace("dx_m = 50000.0", "dx_m = -50000.0"));
    let message = error.to_string();
    assert!(
        message.contains("dx"),
        "the message should name the offending spacing, got: {message}"
    );
}

#[test]
fn an_unknown_forcing_type_is_rejected_by_name() {
    let error = error_from(&VALID_TOML.replace("steady_trade_winds", "hurricane"));
    let message = error.to_string();
    assert!(
        message.contains("hurricane"),
        "the message should name the forcing that was asked for, got: {message}"
    );
    assert!(
        message.contains("steady_trade_winds"),
        "the message should list the forcings that do exist, got: {message}"
    );
}

#[test]
fn a_misspelled_forcing_parameter_is_rejected_rather_than_ignored() {
    // Silently dropping `equatorial_zonal_stres_pa` would run a scenario
    // nobody asked for.
    let error =
        error_from(&VALID_TOML.replace("equatorial_zonal_stress_pa", "equatorial_zonal_stres_pa"));
    let message = error.to_string();
    assert!(
        message.contains("equatorial_zonal_stres_pa"),
        "the message should name the key it did not recognise, got: {message}"
    );
}

#[test]
fn a_westerly_trade_wind_is_rejected() {
    let error = error_from(&VALID_TOML.replace(
        "equatorial_zonal_stress_pa = -0.05",
        "equatorial_zonal_stress_pa = 0.05",
    ));
    let message = error.to_string();
    assert!(
        message.contains("0.05"),
        "the message should name the stress it rejected, got: {message}"
    );
}

#[test]
fn a_non_positive_reduced_gravity_is_rejected_by_name() {
    let error = error_from(&VALID_TOML.replace(
        "reduced_gravity_m_per_s2 = 0.06",
        "reduced_gravity_m_per_s2 = 0.0",
    ));
    let message = error.to_string();
    assert!(
        message.contains("reduced_gravity_m_per_s2"),
        "the message should name the offending parameter, got: {message}"
    );
}

#[test]
fn a_timestep_past_the_cfl_bound_is_rejected_before_the_run_starts() {
    // c = √(0.06 · 150) = 3.0 m/s and dx = dy = 50 km, so the bound of
    // `termocline-numerics` is 0.8 · dx / c ≈ 13 333 s (see its module
    // comment). A day-long step is an order of magnitude past it.
    let error = error_from(&VALID_TOML.replace("dt_s = 3600.0", "dt_s = 86400.0"));
    let message = error.to_string();
    assert!(
        message.contains("86400"),
        "the message should name the timestep asked for, got: {message}"
    );
    assert!(
        message.contains("CFL"),
        "the message should name the bound it violated, got: {message}"
    );
}

#[test]
fn a_malformed_file_is_reported_rather_than_panicked() {
    let error = error_from("[basin\nnx = 200");
    assert!(
        !error.to_string().is_empty(),
        "a parse failure still has to say something"
    );
}

#[test]
fn a_missing_section_is_reported_by_name() {
    let (before_run, _) = VALID_TOML
        .split_once("[run]")
        .expect("the template has a [run] section to drop");
    let error = error_from(before_run);
    let message = error.to_string();
    assert!(
        message.contains("missing field") && message.contains("run"),
        "the message should say which section is missing, got: {message}"
    );
}

#[test]
fn a_scenario_that_is_not_a_file_is_reported_with_its_path() {
    let missing = scenarios_dir().join("no-such-scenario.toml");
    let error = Scenario::load(&missing).expect_err("there is no such file");
    assert!(
        error.to_string().contains("no-such-scenario.toml"),
        "the message should name the file it could not read, got: {error}"
    );
}
