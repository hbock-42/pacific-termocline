//! Loading a run from a directory: the native affordance ADR-0006 keeps, and
//! the one the acceptance criterion's example run arrives by on a desktop.
//!
//! The web path reaches [`LoadedRun`] with the same [`RunBytes`], so what is
//! left to check here is the part only a filesystem has: that the two files
//! are found by the names the format gives them, and that a directory missing
//! one of them says which.

#![cfg(not(target_arch = "wasm32"))]

use std::path::PathBuf;

use termocline_format::{
    frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader,
    Variable, FRAME_FILE_NAME, HEADER_FILE_NAME,
};
use visualizer::{read_run_directory, LoadedRun};

/// A one-frame run on a two-by-two basin, written to a directory of its own.
///
/// Small on purpose: what is under test is the directory, not the grid.
fn write_run(directory: &PathBuf) -> RunHeader {
    let extent = BasinExtent::new(120.0, -80.0, -25.0, 25.0);
    let grid = GridSpec::new(2, 2, extent).expect("2 x 2 is a valid basin");
    let header = RunHeader::new(
        grid,
        PhysicalParams {
            mean_depth_m: 150.0,
            reduced_gravity_m_per_s2: 0.06,
            beta_per_m_per_s: 2.3e-11,
            rayleigh_damping_per_s: 1.0e-7,
            reference_density_kg_per_m3: 1025.0,
        },
        "directory-run",
        OutputTiming {
            frame_count: 1,
            interval_s: 86_400.0,
        },
    );
    let field = |variable| vec![0.0; grid.field_len(variable)];
    let frame = Frame::new(
        0.0,
        &grid,
        field(Variable::ThermoclineDepthAnomaly),
        field(Variable::ZonalCurrentAnomaly),
        field(Variable::MeridionalCurrentAnomaly),
        field(Variable::ZonalWindStress),
        field(Variable::MeridionalWindStress),
    )
    .expect("fields sized from the grid fit it");

    std::fs::create_dir_all(directory).expect("the run directory is created");
    std::fs::write(
        directory.join(HEADER_FILE_NAME),
        serde_json::to_vec(&header).expect("a header serializes"),
    )
    .expect("the header is written");
    std::fs::write(
        directory.join(FRAME_FILE_NAME),
        bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
    )
    .expect("the frames are written");
    header
}

/// A directory of this test's own, named so two tests never share one.
fn scratch(name: &str) -> PathBuf {
    let directory = std::env::temp_dir().join(format!("termocline-viz-{name}"));
    let _ = std::fs::remove_dir_all(&directory);
    directory
}

#[test]
fn a_run_directory_loads_into_the_same_metadata_the_web_path_produces() {
    let directory = scratch("directory-loads");
    let header = write_run(&directory);

    let bytes = read_run_directory(&directory).expect("the run directory is read");
    let run =
        LoadedRun::from_bytes(directory.display().to_string(), &bytes).expect("the run loads");

    assert_eq!(run.header(), &header);
    let grid = run
        .metadata()
        .into_iter()
        .find(|row| row.label == "Grid")
        .expect("a grid row");
    assert_eq!(grid.value, "2 × 2 cells");
}

#[test]
fn a_directory_missing_one_of_the_two_files_names_the_one_it_is_missing() {
    let directory = scratch("directory-half-copied");
    write_run(&directory);
    std::fs::remove_file(directory.join(FRAME_FILE_NAME)).expect("the frames are removed");

    let error = read_run_directory(&directory).expect_err("a half-copied run is refused");
    assert!(
        error.contains(FRAME_FILE_NAME),
        "the message should name the missing file: {error}"
    );
}
