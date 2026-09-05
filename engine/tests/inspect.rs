//! Acceptance tests for T-06.4 — the `inspect` command.
//!
//! The criterion is that *the output matches the header written by T-05.2 for
//! a known test run*. So the run here is a real one — the T-02.5 solver
//! stepped forward and sampled by the T-05.2 writer — and the whole rendering
//! is compared against a summary written out by hand from the scenario
//! constants below. Nothing in the expected text was taken from running the
//! command: a field the renderer forgets, mislabels or prints in the wrong
//! unit fails the comparison.
//!
//! There is no tolerance anywhere in this file. Printing a header is a
//! transcription, not an approximation — a parameter that reached the terminal
//! merely "close" would misreport the run someone is sanity-checking — so
//! every number is compared as the exact text of its shortest round-tripping
//! form.
//!
//! Both sides of the command are covered: the library rendering, and the
//! binary a user actually types, which is run as a subprocess so the exit
//! status and the two streams are the real ones.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{self, Command};
use std::sync::atomic::{AtomicU32, Ordering};

use engine::{
    inspect_run, render_header, BetaPlane, Grid, OceanState, OutputSchedule, PhysicalParams,
    RunWriter, Solver, Spacing, WindStressField, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
use termocline_format::{BasinExtent, GridSpec, RunHeader, FRAME_FILE_NAME};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: an `e`-folding time of about 11.6 days, the
/// value the Epic 02 tests damp at.
const DAMPING_PER_S: f64 = 1.0e-6;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*).
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: an equatorial channel
/// reaching ±500 km. Different from [`BASIN_LX_M`] so an x/y swap cannot pass.
const BASIN_LY_M: f64 = 1.0e6;

/// Cells along x and y of the test basin. Different from one another so a
/// transposed grid line cannot pass.
const NX: usize = 6;
const NY: usize = 4;

/// Where the test basin sits on the globe: the equatorial Pacific of
/// `CONTEXT.md`, 120°E to 80°W, and a ±5° equatorial band. All four differ so
/// a boundary printed in the wrong slot cannot pass.
const WEST_DEG_EAST: f64 = 120.0;
const EAST_DEG_EAST: f64 = -80.0;
const SOUTH_DEG_NORTH: f64 = -5.0;
const NORTH_DEG_NORTH: f64 = 5.0;

/// Free text naming the scenario, as a header carries it.
const SCENARIO: &str = "steady trade winds over a resting basin";

/// The timestep the test run steps at, in seconds: 15 minutes, which
/// `a_run_is_written_inside_the_cfl_bound` checks is inside the gravity-wave
/// CFL maximum for this basin. A round number so the output cadence below is
/// exact in seconds and can be written out by hand.
const DT_S: f64 = 900.0;
/// Length of the test run, in solver steps.
const TOTAL_STEPS: u64 = 12;
/// Output cadence: one frame every fourth step, so the run is a decimated
/// series rather than every step.
const EVERY_N_STEPS: u64 = 4;
/// Frames the cadence above asks for: one at step 0 and one every fourth step
/// through step 12, i.e. steps 0, 4, 8 and 12.
const EXPECTED_FRAMES: u64 = 4;
/// Model time between those frames, in seconds: four steps of 15 minutes.
const EXPECTED_INTERVAL_S: f64 = 3600.0;

/// Zonal wind stress of the test run, in Pa. Easterly trade-wind stress is
/// `τx < 0` (`CONTEXT.md`).
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the test run, in Pa. Different in magnitude from
/// [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

/// The summary the header above must produce, written out by hand from the
/// scenario constants rather than from the command's output. `1` is the
/// current [`termocline_format::FORMAT_VERSION`], stated as a literal for the
/// same reason.
fn expected_summary() -> String {
    "\
format version: 1
scenario: steady trade winds over a resting basin
grid: 6 x 4 cells
basin extent: 120.0 to -80.0 degrees east, -5.0 to 5.0 degrees north
mean thermocline depth H = 150.0 m
reduced gravity g' = 0.05 m s^-2
beta = 2.3e-11 m^-1 s^-1
Rayleigh damping r = 1e-6 s^-1
reference density rho_0 = 1025.0 kg m^-3
frames: 4, one every 3600.0 s
variables: h [m], u [m s^-1], v [m s^-1], tau_x [N m^-2], tau_y [N m^-2]
"
    .to_owned()
}

/// The whole output of inspecting the run in `directory`: the directory it
/// was read from, then the summary of its header.
fn expected_output(directory: &Path) -> String {
    format!("run: {}\n{}", directory.display(), expected_summary())
}

fn params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        DAMPING_PER_S,
        EQUATORIAL_BETA_PER_M_PER_S,
        SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the Pacific parameter set is physical")
}

fn extent() -> BasinExtent {
    BasinExtent::new(
        WEST_DEG_EAST,
        EAST_DEG_EAST,
        SOUTH_DEG_NORTH,
        NORTH_DEG_NORTH,
    )
}

fn basin() -> (Grid, Spacing) {
    let grid = Grid::new(NX, NY).expect("the test basin has cells on both axes");
    let spacing = Spacing::new(BASIN_LX_M / NX as f64, BASIN_LY_M / NY as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

fn schedule() -> OutputSchedule {
    OutputSchedule::new(DT_S, TOTAL_STEPS, EVERY_N_STEPS)
        .expect("a positive timestep and a non-zero cadence are a valid schedule")
}

fn header() -> RunHeader {
    RunHeader::new(
        GridSpec::new(NX, NY, extent()).expect("the test basin has cells on both axes"),
        params().into(),
        SCENARIO,
        schedule().timing(),
    )
}

/// The known test run, written into `dir` by the T-05.2 writer: the T-02.5
/// solver stepped from a single deep spot in the west, sampled at the cadence
/// above.
fn write_run_to_files(dir: &Path) {
    let (grid, spacing) = basin();
    let schedule = schedule();
    let mut writer = RunWriter::create(dir, &header()).expect("the scratch directory is writable");

    let plane = BetaPlane::centered_on_equator(params(), spacing, grid);
    let mut solver = Solver::new(grid, spacing, params(), plane, DT_S)
        .expect("a 15-minute step is inside this basin's CFL bound");
    let wind = WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let mut state = OceanState::at_rest(grid);
    *state
        .h_mut()
        .get_mut(1, 1)
        .expect("(1, 1) is inside a 6 by 4 basin") = 20.0;

    for step in 0..=schedule.total_steps() {
        let t_s = step as f64 * DT_S;
        if schedule.writes_at_step(step) {
            writer
                .append(t_s, &state, &wind)
                .expect("the state covers the basin the header describes");
        }
        if step < schedule.total_steps() {
            solver.step(&mut state, t_s, |_| &wind);
        }
    }
    writer
        .finish()
        .expect("the run wrote every frame it promised");
}

/// A directory under the system temp directory, removed when the test ends.
struct ScratchDir {
    path: PathBuf,
}

impl ScratchDir {
    fn new(name: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("termocline-t064-{name}-{}-{unique}", process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("the system temp directory is writable");
        Self { path }
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for ScratchDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// The `inspect` subcommand of the engine binary, run on `directory`.
fn run_cli(directory: &Path) -> process::Output {
    Command::new(env!("CARGO_BIN_EXE_termocline"))
        .arg("inspect")
        .arg("--run")
        .arg(directory)
        .output()
        .expect("the engine binary is built before its integration tests run")
}

#[test]
fn the_run_is_written_inside_the_cfl_bound() {
    // The scenario's timestep is a hand-picked round number rather than a
    // fraction of the bound, so the test run is only a real run if that number
    // is stable. This is the guard on that constant, not a physics check.
    let (_, spacing) = basin();
    let wave_speed = engine::WaveSpeed::new(params().kelvin_wave_speed_m_per_s())
        .expect("`√(g'·H)` of a physical parameter set is positive");
    assert!(
        DT_S < engine::max_stable_dt(spacing, wave_speed),
        "the test run's timestep is outside the gravity-wave CFL bound"
    );
}

#[test]
fn the_rendered_header_matches_the_run_that_was_written() {
    assert_eq!(render_header(&header()), expected_summary());
}

#[test]
fn inspecting_a_run_directory_reports_its_header() {
    let scratch = ScratchDir::new("directory");
    write_run_to_files(scratch.path());

    let summary = inspect_run(scratch.path()).expect("the run directory holds a readable run");

    assert_eq!(summary, expected_output(scratch.path()));
}

#[test]
fn the_command_prints_the_summary_and_succeeds() {
    let scratch = ScratchDir::new("cli");
    write_run_to_files(scratch.path());

    let output = run_cli(scratch.path());

    assert!(
        output.status.success(),
        "inspecting a readable run failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).expect("the summary is UTF-8"),
        expected_output(scratch.path())
    );
}

#[test]
fn a_run_whose_frames_are_unreadable_still_reports_its_header() {
    // The point of the command is sanity-checking a run without the
    // visualizer, including a run a crash cut short. The header is on disk
    // from the run's first moment, so inspecting it must not depend on the
    // frames beside it decoding.
    let scratch = ScratchDir::new("truncated");
    write_run_to_files(scratch.path());
    fs::write(scratch.path().join(FRAME_FILE_NAME), b"not a frame")
        .expect("the scratch directory is writable");

    let summary = inspect_run(scratch.path()).expect("the header is readable on its own");

    assert_eq!(summary, expected_output(scratch.path()));
}

#[test]
fn a_run_missing_its_frame_file_is_an_actionable_error() {
    // The command reads runs through `RunReader`, which opens both of a run's
    // files (ADR-0004). A directory holding only a header is half a run, so it
    // is refused by name rather than reported as if it were whole.
    let scratch = ScratchDir::new("no-frames");
    write_run_to_files(scratch.path());
    let frames = scratch.path().join(FRAME_FILE_NAME);
    fs::remove_file(&frames).expect("the run wrote a frame file");

    let output = run_cli(scratch.path());
    let stderr = String::from_utf8(output.stderr).expect("the error message is UTF-8");

    assert!(
        !output.status.success(),
        "inspecting a run with no frame file reported success"
    );
    assert!(
        stderr.contains(&frames.display().to_string()),
        "the error does not name the file it could not open: {stderr}"
    );
}

#[test]
fn a_missing_run_directory_is_an_actionable_error() {
    let scratch = ScratchDir::new("missing");
    let absent = scratch.path().join("no-such-run");

    let output = run_cli(&absent);
    let stderr = String::from_utf8(output.stderr).expect("the error message is UTF-8");

    assert!(
        !output.status.success(),
        "inspecting a directory that does not exist reported success"
    );
    assert!(
        stderr.contains(&absent.display().to_string()),
        "the error does not name the run it could not read: {stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "an unreadable run panicked rather than returning an error: {stderr}"
    );
}

#[test]
fn the_test_run_writes_the_frames_the_summary_reports() {
    // The frame count and cadence in `expected_summary` are written out by
    // hand; this is the guard that the scenario above actually asks for them,
    // so a summary that merely echoed a different schedule could not pass.
    let timing = schedule().timing();
    assert_eq!(timing.frame_count, EXPECTED_FRAMES);
    assert_eq!(timing.interval_s.to_bits(), EXPECTED_INTERVAL_S.to_bits());
}
