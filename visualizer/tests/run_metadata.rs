//! T-08.1 acceptance criterion: opening one of the example runs shows correct
//! header metadata — grid size, scenario name, frame count — in the UI.
//!
//! The rows the UI draws are built here, away from a window and a GPU, so the
//! claim "the shell shows the right numbers" is an assertion rather than a
//! screenshot. Drawing them is a one-liner over [`LoadedRun::metadata`]; what
//! could be wrong is the values, and they are checked below.
//!
//! The expected values come from `engine/scenarios/steady-trades.toml` and
//! `CONTEXT.md`, *Basin*: 120°E–80°W by 25°S–25°N at 0.5°, which is 320 × 100
//! cells. Nothing here was read back off a run this code produced.

use termocline_format::{
    frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader,
    FORMAT_VERSION,
};
use visualizer::{LoadedRun, RunBytes};

/// The basin of `CONTEXT.md`, *Basin*.
const PACIFIC: BasinExtent = BasinExtent::new(120.0, -80.0, -25.0, 25.0);
/// 50° of latitude and 160° of longitude at the scenario's 0.5° resolution.
const NX: usize = 320;
const NY: usize = 100;
/// `steady-trades.toml` writes a frame every 24 steps of an hour.
const FRAME_INTERVAL_S: f64 = 86_400.0;

/// The header `steady-trades.toml` produces, built from the scenario file
/// rather than from a run this code read back.
fn steady_trades_header(frame_count: u64) -> RunHeader {
    let grid = GridSpec::new(NX, NY, PACIFIC).expect("320 x 100 is a valid basin");
    RunHeader::new(
        grid,
        PhysicalParams {
            mean_depth_m: 150.0,
            reduced_gravity_m_per_s2: 0.06,
            beta_per_m_per_s: 2.3e-11,
            rayleigh_damping_per_s: 1.0e-7,
            reference_density_kg_per_m3: 1025.0,
        },
        "steady-trades",
        OutputTiming {
            frame_count,
            interval_s: FRAME_INTERVAL_S,
        },
    )
}

/// A run's two byte sources, as they would arrive from a file drop or a fetch.
fn run_bytes(header: &RunHeader) -> RunBytes {
    let grid = header.grid;
    let field = |variable| vec![0.0; grid.field_len(variable)];
    let mut frames = Vec::new();
    for index in 0..header.output.frame_count {
        #[allow(clippy::cast_precision_loss)]
        let t_s = index as f64 * header.output.interval_s;
        let frame = Frame::new(
            t_s,
            &grid,
            field(termocline_format::Variable::ThermoclineDepthAnomaly),
            field(termocline_format::Variable::ZonalCurrentAnomaly),
            field(termocline_format::Variable::MeridionalCurrentAnomaly),
            field(termocline_format::Variable::ZonalWindStress),
            field(termocline_format::Variable::MeridionalWindStress),
        )
        .expect("fields sized from the grid fit it");
        frames.extend(
            bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
        );
    }
    RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames,
    }
}

/// The value of the row labelled `label`, or a panic naming what was there.
fn row(run: &LoadedRun, label: &str) -> String {
    let rows = run.metadata();
    rows.iter()
        .find(|row| row.label == label)
        .unwrap_or_else(|| panic!("no row {label:?} among {:?}", rows))
        .value
        .clone()
}

#[test]
fn shows_the_scenario_name_grid_size_and_frame_count_of_a_loaded_run() {
    let header = steady_trades_header(11);
    let run = LoadedRun::from_bytes("run-demo", &run_bytes(&header)).expect("the run loads");

    assert_eq!(row(&run, "Scenario"), "steady-trades");
    assert_eq!(row(&run, "Grid"), "320 × 100 cells");
    assert_eq!(row(&run, "Frames"), "11");
}

#[test]
fn names_where_the_run_came_from() {
    // The shell loads from a directory, a pair of dropped files or a URL, and
    // says which — a run that loaded from the wrong place looks identical
    // otherwise.
    let run = LoadedRun::from_bytes("/tmp/run-demo", &run_bytes(&steady_trades_header(1)))
        .expect("loads");
    assert_eq!(run.source(), "/tmp/run-demo");
}

#[test]
fn shows_the_basin_the_run_covers() {
    // The extent is a west-to-east span across the antimeridian, not a min and
    // a max, so 120.0 to -80.0 reads as 120°E to 80°W and never as a
    // 200°-wide basin running the other way.
    let run = LoadedRun::from_bytes("run", &run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(row(&run, "Basin"), "120.0°E – 80.0°W, 25.0°S – 25.0°N");
}

#[test]
fn shows_the_output_cadence_and_the_model_time_it_spans() {
    // 11 frames a day apart span 10 intervals — ten days of model time, not
    // eleven; an off-by-one here would misdate every frame Epic 09 draws.
    let run = LoadedRun::from_bytes("run", &run_bytes(&steady_trades_header(11))).expect("loads");
    assert_eq!(row(&run, "Frame interval"), "86400 s (1.00 days)");
    assert_eq!(row(&run, "Model time"), "10.00 days");
}

#[test]
fn a_run_with_no_frames_spans_no_model_time() {
    let run = LoadedRun::from_bytes("run", &run_bytes(&steady_trades_header(0))).expect("loads");
    assert_eq!(row(&run, "Frames"), "0");
    assert_eq!(row(&run, "Model time"), "0.00 days");
}

#[test]
fn reports_the_variables_each_frame_carries() {
    let run = LoadedRun::from_bytes("run", &run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(row(&run, "Variables"), "h, u, v, tau_x, tau_y");
    assert_eq!(row(&run, "Format version"), FORMAT_VERSION.to_string());
}

#[test]
fn a_run_from_an_unknown_format_version_is_refused_rather_than_shown() {
    // The frames of a future version may have any layout, so a header this
    // build does not understand is an error the shell reports, not metadata it
    // shows with the wrong labels.
    let mut header = steady_trades_header(1);
    let bytes = run_bytes(&header);
    header.format_version = FORMAT_VERSION + 1;
    let bumped = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: bytes.frames,
    };
    let error = LoadedRun::from_bytes("run", &bumped).expect_err("the version is not readable");
    let message = error.to_string();
    assert!(
        message.contains(&(FORMAT_VERSION + 1).to_string()),
        "the message should name the version it found: {message}"
    );
}

#[test]
fn a_header_that_is_not_a_header_is_an_error_not_a_panic() {
    let bytes = RunBytes {
        header: b"{ not json".to_vec(),
        frames: Vec::new(),
    };
    LoadedRun::from_bytes("run", &bytes).expect_err("malformed input is refused");
}

#[test]
fn writes_small_physical_parameters_in_scientific_notation() {
    // β = 2.3e-11 m^-1 s^-1 written out in full is eleven leading zeros; a
    // reader checking the run was integrated with the parameters they meant
    // should not have to count them.
    let run = LoadedRun::from_bytes("run", &run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(row(&run, "Coriolis gradient β"), "2.3e-11 m^-1 s^-1");
    assert_eq!(row(&run, "Rayleigh damping r"), "1e-7 s^-1");
    // Values of everyday magnitude stay as they are written in the scenario.
    assert_eq!(row(&run, "Mean depth H"), "150 m");
    assert_eq!(row(&run, "Reduced gravity g'"), "0.06 m s^-2");
}
