//! The acceptance criterion of T-05.3 that spans both crates: *reading a run
//! produced by T-05.2 yields frames identical to what was written, in order*.
//!
//! So the run here is a real one — the T-02.5 solver stepped forward, sampled
//! by the T-05.2 writer at a cadence coarser than its timestep — and it is
//! read back through `RunReader` and compared against the states that went in.
//! The reader's own behaviour on hand-built runs, and its memory bound, are
//! tested inside `termocline-format`; what this file adds is that the writer
//! and the reader agree on the bytes between them.
//!
//! Nothing here has a tolerance. A round trip through a file is an identity,
//! not an approximation — a value that came back merely "close" would silently
//! corrupt an archived run — so every field is compared as an IEEE-754 bit
//! pattern.
//!
//! Per [ADR-0006] the same run is read twice: once from the two files on disk
//! through the `fs` convenience, and once from the two byte buffers a browser
//! would hold. The frames must be the same either way, because the source of
//! the bytes is not part of the format.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU32, Ordering};

use engine::{
    BetaPlane, Grid, OceanState, OutputSchedule, PhysicalParams, RunWriter, Solver, Spacing,
    WindStress,
};
use termocline_format::{
    BasinExtent, Frame, GridSpec, RunHeader, RunReadError, RunReader, Variable,
};

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
/// transposed field cannot pass.
const NX: usize = 6;
const NY: usize = 4;

/// Amplitude of the initial thermocline depth anomaly, in metres. A 20 m
/// departure is the scale of an observed equatorial Pacific anomaly.
const H_AMPLITUDE_M: f64 = 20.0;

/// Zonal wind stress of the test run, in Pa. Easterly trade-wind stress is
/// `τx < 0` (`CONTEXT.md`).
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the test run, in Pa. Different in magnitude from
/// [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

/// Length of the test run, in solver steps.
const TOTAL_STEPS: u64 = 12;
/// Output cadence: one frame every fourth step, so the run read back is a
/// decimated series rather than every step.
const EVERY_N_STEPS: u64 = 4;
/// Frames the cadence above asks for: one at step 0 and one every fourth step
/// through step 12, i.e. steps 0, 4, 8 and 12. Written out by hand rather than
/// taken from the code under test.
const EXPECTED_FRAMES: u64 = 4;

/// The equatorial Pacific basin of `CONTEXT.md`: 120°E–80°W, 25°S–25°N.
fn extent() -> BasinExtent {
    BasinExtent::new(120.0, -80.0, -25.0, 25.0)
}

fn params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        DAMPING_PER_S,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

fn basin() -> (Grid, Spacing) {
    let grid = Grid::new(NX, NY).expect("cell counts are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / NX as f64, BASIN_LY_M / NY as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

/// The timestep the run steps at, in seconds: half the gravity-wave CFL
/// maximum for this basin, so the run is comfortably inside both stability
/// bounds.
fn timestep_s(spacing: Spacing) -> f64 {
    let wave_speed = engine::WaveSpeed::new(params().kelvin_wave_speed_m_per_s())
        .expect("`√(g'·H)` of a physical parameter set is positive");
    engine::max_stable_dt(spacing, wave_speed) / 2.0
}

fn schedule(spacing: Spacing) -> OutputSchedule {
    OutputSchedule::new(timestep_s(spacing), TOTAL_STEPS, EVERY_N_STEPS)
        .expect("a positive timestep and a non-zero cadence are a valid schedule")
}

fn header(spacing: Spacing) -> RunHeader {
    RunHeader::new(
        GridSpec::new(NX, NY, extent()).expect("the test basin has cells on both axes"),
        params().into(),
        "steady alizes over a resting basin",
        schedule(spacing).timing(),
    )
}

/// The initial state: a single deep spot in the west, so the fields are
/// neither uniform nor symmetric and a mislabelled or transposed field cannot
/// survive the comparison.
fn initial_state(grid: Grid) -> OceanState {
    let mut state = OceanState::at_rest(grid);
    *state
        .h_mut()
        .get_mut(1, 1)
        .expect("(1, 1) is inside a 6 by 4 basin") = H_AMPLITUDE_M;
    state
}

/// The test run, written through `writer`. Returns the model times and states
/// that were written, in order, so the reader's output can be checked against
/// them.
fn write_run<W: std::io::Write>(
    writer: &mut RunWriter<W>,
    grid: Grid,
    spacing: Spacing,
) -> Vec<(f64, OceanState)> {
    let schedule = schedule(spacing);
    let dt_s = schedule.dt_s();
    let plane = BetaPlane::centered_on_equator(params(), spacing, grid);
    let mut solver =
        Solver::new(grid, spacing, params(), plane, dt_s).expect("half the CFL maximum is stable");
    let wind = WindStress::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let mut state = initial_state(grid);
    let mut written = Vec::new();
    for step in 0..=schedule.total_steps() {
        let t_s = step as f64 * dt_s;
        if schedule.writes_at_step(step) {
            writer
                .append(t_s, &state, &wind)
                .expect("the state covers the basin the header describes");
            written.push((t_s, state.clone()));
        }
        if step < schedule.total_steps() {
            solver.step(&mut state, t_s, |_| &wind);
        }
    }
    written
}

/// The run written into a scratch directory, with the states that went in.
fn write_run_to_files(dir: &Path) -> Vec<(f64, OceanState)> {
    let (grid, spacing) = basin();
    let mut writer =
        RunWriter::create(dir, &header(spacing)).expect("the scratch directory is writable");
    let written = write_run(&mut writer, grid, spacing);
    writer
        .finish()
        .expect("the run wrote every frame it promised");
    written
}

/// The same run written into memory: the two files as byte buffers, which is
/// how a run reaches a browser (ADR-0006).
fn write_run_to_memory() -> (Vec<u8>, Vec<u8>, Vec<(f64, OceanState)>) {
    let (grid, spacing) = basin();
    let mut header_bytes = Vec::new();
    let mut writer = RunWriter::new(&mut header_bytes, Vec::new(), &header(spacing))
        .expect("a vector never fails a write");
    let written = write_run(&mut writer, grid, spacing);
    let frame_bytes = writer
        .finish()
        .expect("the run wrote every frame it promised");
    (header_bytes, frame_bytes, written)
}

fn assert_bit_identical(context: &str, written: &[f64], read: &[f64]) {
    assert_eq!(written.len(), read.len(), "{context}: length changed");
    for (i, (w, r)) in written.iter().zip(read).enumerate() {
        assert_eq!(
            w.to_bits(),
            r.to_bits(),
            "{context}[{i}]: {w} came back as {r}"
        );
    }
}

/// Every frame of `read` against the states that were written, in order.
fn assert_run_round_trips(written: &[(f64, OceanState)], read: &[Frame]) {
    assert_eq!(
        read.len(),
        written.len(),
        "the reader returned a different number of frames than the run wrote"
    );
    for (index, (frame, (t_s, state))) in read.iter().zip(written).enumerate() {
        assert_eq!(
            frame.t_s().to_bits(),
            t_s.to_bits(),
            "frame {index}: model time changed"
        );
        assert_bit_identical(&format!("frame {index} h"), state.h().as_slice(), frame.h());
        assert_bit_identical(&format!("frame {index} u"), state.u().as_slice(), frame.u());
        assert_bit_identical(&format!("frame {index} v"), state.v().as_slice(), frame.v());
        // The forcing is uniform and constant over this run, so it is checked
        // against the stress the run was driven with rather than a buffer.
        for value in frame.tau_x() {
            assert_eq!(value.to_bits(), TRADE_WIND_STRESS_X_PA.to_bits());
        }
        for value in frame.tau_y() {
            assert_eq!(value.to_bits(), TRADE_WIND_STRESS_Y_PA.to_bits());
        }
    }
}

fn collect(reader: RunReader<impl std::io::Read>) -> Vec<Frame> {
    reader
        .collect::<Result<Vec<Frame>, RunReadError>>()
        .expect("a run this build wrote is a run this build reads")
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
            std::env::temp_dir().join(format!("termocline-t053-{name}-{}-{unique}", process::id()));
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

#[test]
fn a_run_written_to_files_reads_back_frame_for_frame() {
    let scratch = ScratchDir::new("files");
    let written = write_run_to_files(scratch.path());

    let reader = RunReader::open(scratch.path()).expect("the run directory holds a readable run");
    assert_eq!(reader.header().output.frame_count, EXPECTED_FRAMES);
    assert_run_round_trips(&written, &collect(reader));
}

#[test]
fn a_run_held_in_memory_reads_back_frame_for_frame() {
    // The browser's path (ADR-0006): two buffers, no filesystem anywhere.
    let (header_bytes, frame_bytes, written) = write_run_to_memory();

    let reader = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the buffers hold a readable run");
    assert_run_round_trips(&written, &collect(reader));
}

#[test]
fn a_run_reads_back_the_same_whether_it_came_from_a_file_or_a_buffer() {
    // The source of the bytes is not part of the format, so the two paths are
    // required to agree — not merely to work.
    let scratch = ScratchDir::new("agree");
    write_run_to_files(scratch.path());
    let (header_bytes, frame_bytes, _) = write_run_to_memory();

    let from_files = RunReader::open(scratch.path()).expect("the run directory is readable");
    let from_memory = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the buffers hold a readable run");
    assert_eq!(from_files.header(), from_memory.header());
    assert_eq!(collect(from_files), collect(from_memory));
}

#[test]
fn the_header_describes_the_run_before_any_frame_is_read() {
    // The deliverable's header accessor: a caller sizes its buffers and labels
    // its axes from the header, without paying for a frame.
    let (_, spacing) = basin();
    let (header_bytes, frame_bytes, _) = write_run_to_memory();

    let reader = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the buffers hold a readable run");
    let read = reader.header();

    assert_eq!(read, &header(spacing));
    assert_eq!(read.grid.nx(), NX);
    assert_eq!(read.grid.ny(), NY);
    assert_eq!(read.output.frame_count, EXPECTED_FRAMES);
    assert_eq!(
        read.output.interval_s,
        EVERY_N_STEPS as f64 * timestep_s(spacing)
    );
    assert_eq!(read.physical_params.mean_depth_m, PACIFIC_MEAN_DEPTH_M);
    let symbols: Vec<&str> = read.variables.iter().map(|v| v.symbol.as_str()).collect();
    assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y"]);
    assert_eq!(
        Variable::ALL.len(),
        symbols.len(),
        "the header lists every variable a frame carries"
    );
}

#[test]
fn a_run_whose_frame_file_was_lost_is_refused_by_name() {
    // A half-copied run directory is invalid input, not a broken invariant:
    // the reader says which file it could not open rather than panicking.
    let scratch = ScratchDir::new("lost-frames");
    write_run_to_files(scratch.path());
    fs::remove_file(scratch.path().join(engine::FRAME_FILE_NAME)).expect("the frame file exists");

    let error = RunReader::open(scratch.path()).expect_err("a run without frames is not a run");
    let message = error.to_string();
    assert!(matches!(error, RunReadError::Open { .. }), "{error:?}");
    assert!(
        message.contains(engine::FRAME_FILE_NAME),
        "the error names the file it could not open: {message}"
    );
}

#[test]
fn a_run_whose_header_was_lost_is_refused_by_name() {
    // The other half of the same mistake, and it must not be reported as a
    // frame failure: the two files are opened by name, so the error says which
    // one is missing.
    let scratch = ScratchDir::new("lost-header");
    write_run_to_files(scratch.path());
    fs::remove_file(scratch.path().join(engine::HEADER_FILE_NAME)).expect("the header exists");

    let error = RunReader::open(scratch.path()).expect_err("a run without a header is not a run");
    let message = error.to_string();
    assert!(matches!(error, RunReadError::Open { .. }), "{error:?}");
    assert!(
        message.contains(engine::HEADER_FILE_NAME),
        "the error names the file it could not open: {message}"
    );
}
