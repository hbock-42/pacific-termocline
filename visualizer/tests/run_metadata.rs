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

mod common;

use common::{encoded_frames, run_bytes, steady_trades_header};
use termocline_format::FORMAT_VERSION;
use visualizer::{LoadedRun, RunBytes};

/// The value of the row labelled `label`, or a panic naming what was there.
fn row(run: &LoadedRun, label: &str) -> String {
    let rows = run.metadata();
    rows.iter()
        .find(|row| row.label == label)
        .unwrap_or_else(|| panic!("no row {label:?} among {rows:?}"))
        .value
        .clone()
}

#[test]
fn shows_the_scenario_name_grid_size_and_frame_count_of_a_loaded_run() {
    let header = steady_trades_header(11);
    let run = LoadedRun::from_bytes("run-demo", run_bytes(&header)).expect("the run loads");

    assert_eq!(row(&run, "Scenario"), "steady-trades");
    assert_eq!(row(&run, "Grid"), "320 × 100 cells");
    assert_eq!(row(&run, "Frames"), "11");
}

#[test]
fn names_where_the_run_came_from() {
    // The shell loads from a directory, a pair of dropped files or a URL, and
    // says which — a run that loaded from the wrong place looks identical
    // otherwise.
    let run =
        LoadedRun::from_bytes("/tmp/run-demo", run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(run.source(), "/tmp/run-demo");
}

#[test]
fn shows_the_basin_the_run_covers() {
    // The extent is a west-to-east span across the antimeridian, not a min and
    // a max, so 120.0 to -80.0 reads as 120°E to 80°W and never as a
    // 200°-wide basin running the other way.
    let run = LoadedRun::from_bytes("run", run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(row(&run, "Basin"), "120.0°E – 80.0°W, 25.0°S – 25.0°N");
}

#[test]
fn shows_the_output_cadence_and_the_model_time_it_spans() {
    // 11 frames a day apart span 10 intervals — ten days of model time, not
    // eleven; an off-by-one here would misdate every frame Epic 09 draws.
    let run = LoadedRun::from_bytes("run", run_bytes(&steady_trades_header(11))).expect("loads");
    assert_eq!(row(&run, "Frame interval"), "86400 s (1.00 days)");
    assert_eq!(row(&run, "Model time"), "10.00 days");
}

#[test]
fn a_run_with_no_frames_spans_no_model_time() {
    let run = LoadedRun::from_bytes("run", run_bytes(&steady_trades_header(0))).expect("loads");
    assert_eq!(row(&run, "Frames"), "0");
    assert_eq!(row(&run, "Model time"), "0.00 days");
}

#[test]
fn reports_the_variables_each_frame_carries() {
    let run = LoadedRun::from_bytes("run", run_bytes(&steady_trades_header(1))).expect("loads");
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
    let error = LoadedRun::from_bytes("run", bumped).expect_err("the version is not readable");
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
    LoadedRun::from_bytes("run", bytes).expect_err("malformed input is refused");
}

#[test]
fn writes_small_physical_parameters_in_scientific_notation() {
    // β = 2.3e-11 m^-1 s^-1 written out in full is eleven leading zeros; a
    // reader checking the run was integrated with the parameters they meant
    // should not have to count them.
    let run = LoadedRun::from_bytes("run", run_bytes(&steady_trades_header(1))).expect("loads");
    assert_eq!(row(&run, "Coriolis gradient β"), "2.3e-11 m^-1 s^-1");
    assert_eq!(row(&run, "Rayleigh damping r"), "1e-7 s^-1");
    // Values of everyday magnitude stay as they are written in the scenario.
    assert_eq!(row(&run, "Mean depth H"), "150 m");
    assert_eq!(row(&run, "Reduced gravity g'"), "0.06 m s^-2");
}

#[test]
fn a_frame_count_the_file_does_not_keep_is_refused_rather_than_shown() {
    // The header is a claim about the file beside it, and the panel reports it
    // as fact. A header promising three frames next to two — one run's header
    // beside another run's frames — must not read as a shorter run that
    // happens to be fine.
    let header = steady_trades_header(3);
    let short = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames(&header, 2),
    };
    let error = LoadedRun::from_bytes("run", short).expect_err("a short run is refused");
    assert!(
        error.to_string().contains('3'),
        "the message should name the count that was promised: {error}"
    );

    // And the same in the other direction: frames past the promised count are
    // data nothing would ever draw.
    let long = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames(&header, 4),
    };
    LoadedRun::from_bytes("run", long).expect_err("an overlong run is refused");
}
