//! Acceptance tests for T-05.2 — the engine-side run writer.
//!
//! The criterion is that *a short test run produces a header file and a frame
//! file readable back through `termocline-format`, with the expected number of
//! frames for the configured output interval*. So the run here is a real one:
//! the T-02.5 solver stepping the full right-hand side, with the writer
//! sampling it at a cadence coarser than the timestep.
//!
//! Nothing here has a physical tolerance, because nothing here is physics.
//! Writing a run is an identity, not an approximation — a value that came back
//! merely "close" would silently corrupt an archived run — so the fields are
//! compared as IEEE-754 bit patterns and the frame count is compared against
//! the cadence arithmetic written out by hand.
//!
//! # Why the reader here is `serde_json` plus `bincode` rather than a `RunReader`
//!
//! T-05.3 is the reader ticket, and it is blocked by this one. What
//! `termocline-format` offers today is the `serde` types themselves, so
//! reading a run back means decoding those types out of the two files — which
//! is exactly what T-05.3 will wrap. Per [ADR-0006] that reader is defined
//! over a byte source rather than a path, so these tests read every run out of
//! a plain byte stream, front to back, with no seeking and no trailer: if the
//! format written here needed either, the web reader could not be built on it.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU32, Ordering};

use engine::{
    BetaPlane, Grid, OceanState, OutputSchedule, OutputScheduleError, PhysicalParams,
    RunWriteError, RunWriter, Solver, Spacing, WindStressField,
};
use termocline_format::{
    frame_encoding, BasinExtent, FormatError, Frame, GridSpec, RunHeader, Variable, FORMAT_VERSION,
    FRAME_FILE_NAME, HEADER_FILE_NAME,
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
/// reaching ±500 km, which leaves the gravity-wave CFL bound the binding one
/// (see `time_stepping.rs`). Different from [`BASIN_LX_M`] so an x/y swap
/// cannot pass.
const BASIN_LY_M: f64 = 1.0e6;

/// Cells along x and y of the test basin. Different from one another so a
/// transposed field cannot pass, and small enough to write out by hand.
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

/// Length of the short test run, in solver steps.
const TOTAL_STEPS: u64 = 12;
/// Output cadence: one frame every fourth step, so the run is decimated rather
/// than saved whole — which is the point of the ticket.
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

/// A basin of `nx` by `ny` cells spanning [`BASIN_LX_M`] by [`BASIN_LY_M`].
fn basin(nx: usize, ny: usize) -> (Grid, Spacing) {
    let grid = Grid::new(nx, ny).expect("cell counts are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

/// The timestep the test run steps at, in seconds: half the gravity-wave CFL
/// maximum for this basin, so the run is comfortably inside both stability
/// bounds and `Solver::new` has no reason to refuse it.
fn timestep_s(spacing: Spacing) -> f64 {
    let wave_speed = engine::WaveSpeed::new(params().kelvin_wave_speed_m_per_s())
        .expect("`√(g'·H)` of a physical parameter set is positive");
    engine::max_stable_dt(spacing, wave_speed) / 2.0
}

fn schedule(spacing: Spacing) -> OutputSchedule {
    OutputSchedule::new(timestep_s(spacing), TOTAL_STEPS, EVERY_N_STEPS)
        .expect("a positive timestep and a non-zero cadence are a valid schedule")
}

fn grid_spec(nx: usize, ny: usize) -> GridSpec {
    GridSpec::new(nx, ny, extent()).expect("the test basin has cells on both axes")
}

fn header(spacing: Spacing, nx: usize, ny: usize) -> RunHeader {
    RunHeader::new(
        grid_spec(nx, ny),
        params().into(),
        "steady alizes over a resting basin",
        schedule(spacing).timing(),
    )
}

/// The initial state of the test run: a single deep spot in the west, so the
/// fields are neither uniform nor symmetric and a mislabelled or transposed
/// field cannot survive the comparison.
fn initial_state(grid: Grid) -> OceanState {
    let mut state = OceanState::at_rest(grid);
    *state
        .h_mut()
        .get_mut(1, 1)
        .expect("(1, 1) is inside a 6 by 4 basin") = H_AMPLITUDE_M;
    state
}

/// The run of the acceptance criteria, written through `writer`: `TOTAL_STEPS`
/// steps of the T-02.5 solver, with a frame appended whenever the schedule
/// says so. Returns the states that were written, in order, so a reader's
/// output can be checked against them.
fn write_short_run<W: std::io::Write>(
    writer: &mut RunWriter<W>,
    grid: Grid,
    spacing: Spacing,
) -> Vec<(f64, OceanState)> {
    let schedule = schedule(spacing);
    let dt_s = schedule.dt_s();
    let plane = BetaPlane::centered_on_equator(params(), spacing, grid);
    let mut solver =
        Solver::new(grid, spacing, params(), plane, dt_s).expect("half the CFL maximum is stable");
    let wind = WindStressField::uniform_including_walls(
        grid,
        TRADE_WIND_STRESS_X_PA,
        TRADE_WIND_STRESS_Y_PA,
    );

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

/// A run written into a scratch directory, with the states that went in.
fn write_short_run_to_files(dir: &Path) -> Vec<(f64, OceanState)> {
    let (grid, spacing) = basin(NX, NY);
    let mut writer = RunWriter::create(dir, &header(spacing, NX, NY))
        .expect("the scratch directory is writable");
    let written = write_short_run(&mut writer, grid, spacing);
    writer
        .finish()
        .expect("the run wrote every frame it promised");
    written
}

/// The same run written into memory: the two files as byte buffers.
fn write_short_run_to_memory() -> (Vec<u8>, Vec<u8>) {
    let (grid, spacing) = basin(NX, NY);
    let mut header_bytes = Vec::new();
    let mut writer = RunWriter::new(&mut header_bytes, Vec::new(), &header(spacing, NX, NY))
        .expect("a vector never fails a write");
    write_short_run(&mut writer, grid, spacing);
    let frame_bytes = writer
        .finish()
        .expect("the run wrote every frame it promised");
    (header_bytes, frame_bytes)
}

/// Decode a header out of a byte source, the way a reader would.
fn read_header(bytes: &[u8]) -> RunHeader {
    serde_json::from_slice(bytes).expect("the writer wrote a JSON header")
}

/// Decode every frame of a run out of a byte stream, front to back, taking the
/// count from the header exactly as [ADR-0006]'s source-agnostic reader must:
/// no seek, no index, no trailer.
///
/// [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md
fn read_frames(header: &RunHeader, bytes: &[u8]) -> Vec<Frame> {
    let config = frame_encoding();
    let mut cursor = std::io::Cursor::new(bytes);
    let mut frames = Vec::new();
    for index in 0..header.output.frame_count {
        let frame: Frame = bincode::serde::decode_from_std_read(&mut cursor, config)
            .unwrap_or_else(|error| panic!("frame {index} decodes: {error}"));
        frame
            .validate(&header.grid)
            .expect("every frame fits the grid the header describes");
        frames.push(frame);
    }
    assert_eq!(
        cursor.position() as usize,
        bytes.len(),
        "the frame file holds exactly the frames the header promised, and nothing after them"
    );
    frames
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

/// A directory under the system temp directory, removed when the test ends.
struct ScratchDir {
    path: PathBuf,
}

impl ScratchDir {
    fn new(name: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("termocline-t052-{name}-{}-{unique}", process::id()));
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
fn a_short_run_produces_a_header_file_and_a_frame_file() {
    // The acceptance criterion, end to end: run the solver, write it out, and
    // read both files back through `termocline-format`'s types.
    let scratch = ScratchDir::new("short-run");
    let written = write_short_run_to_files(scratch.path());

    let header_bytes =
        fs::read(scratch.path().join(HEADER_FILE_NAME)).expect("the header file was written");
    let frame_bytes =
        fs::read(scratch.path().join(FRAME_FILE_NAME)).expect("the frame file was written");

    let header = read_header(&header_bytes);
    assert_eq!(header.format_version, FORMAT_VERSION);
    assert_eq!(header.grid.nx(), NX);
    assert_eq!(header.grid.ny(), NY);

    let frames = read_frames(&header, &frame_bytes);
    // Steps 0, 4, 8 and 12 of a twelve-step run at a cadence of four.
    assert_eq!(frames.len() as u64, EXPECTED_FRAMES);
    assert_eq!(written.len() as u64, EXPECTED_FRAMES);
}

#[test]
fn the_header_records_the_frame_count_the_output_interval_asks_for() {
    // A reader over a byte source has only the header to tell it how many
    // frames to expect, so the count it carries is load-bearing, not a note.
    let (_, spacing) = basin(NX, NY);
    let schedule = schedule(spacing);
    let timing = schedule.timing();

    assert_eq!(timing.frame_count, EXPECTED_FRAMES);
    // The output interval is the cadence in steps times the timestep, not the
    // timestep: this is the decimation the ticket exists for.
    assert_eq!(
        timing.interval_s.to_bits(),
        (EVERY_N_STEPS as f64 * schedule.dt_s()).to_bits()
    );
    assert!(
        timing.interval_s > schedule.dt_s(),
        "a decimated run writes less often than it steps"
    );
}

#[test]
fn the_output_interval_decides_which_steps_are_written() {
    // The cadence arithmetic, written out by hand: a run of `total` steps at a
    // cadence of `every` writes at steps 0, every, 2·every, ... up to and
    // including `total`, which is `total / every + 1` frames (integer
    // division, so a run whose length is not a whole number of intervals
    // simply stops at the last one that fits).
    for (total, every, expected_frames) in [
        (12_u64, 4_u64, 4_u64),
        (12, 1, 13),
        (12, 12, 2),
        (12, 5, 3), // steps 0, 5, 10 — step 12 is not a multiple of 5
        (0, 4, 1),  // a run of no steps still writes its initial state
    ] {
        let schedule = OutputSchedule::new(1.0, total, every)
            .expect("a positive timestep and a non-zero cadence are a valid schedule");
        assert_eq!(
            schedule.frame_count(),
            expected_frames,
            "{total} steps at a cadence of {every}"
        );

        let written: Vec<u64> = (0..=total)
            .filter(|s| schedule.writes_at_step(*s))
            .collect();
        let by_hand: Vec<u64> = (0..expected_frames).map(|k| k * every).collect();
        assert_eq!(written, by_hand, "{total} steps at a cadence of {every}");
    }
}

#[test]
fn every_frame_comes_back_exactly_as_it_was_written() {
    // Serialization is an identity, not an approximation, so the fields are
    // compared as bit patterns: a run that came back "close" would be a
    // corrupted archive.
    let scratch = ScratchDir::new("readback");
    let written = write_short_run_to_files(scratch.path());

    let header = read_header(
        &fs::read(scratch.path().join(HEADER_FILE_NAME)).expect("the header file was written"),
    );
    let frames = read_frames(
        &header,
        &fs::read(scratch.path().join(FRAME_FILE_NAME)).expect("the frame file was written"),
    );

    assert_eq!(frames.len(), written.len());
    for (index, (frame, (t_s, state))) in frames.iter().zip(&written).enumerate() {
        assert_eq!(
            frame.t_s().to_bits(),
            t_s.to_bits(),
            "frame {index} carries the model time it was written at"
        );
        assert_bit_identical(&format!("frame {index} h"), state.h().as_slice(), frame.h());
        assert_bit_identical(&format!("frame {index} u"), state.u().as_slice(), frame.u());
        assert_bit_identical(&format!("frame {index} v"), state.v().as_slice(), frame.v());
        // The forcing is constant in this run, so every frame carries the same
        // stress the solver was driven with — in the format's `N m^-2`, which
        // is the pascal the engine states it in.
        for tau in frame.tau_x() {
            assert_eq!(tau.to_bits(), TRADE_WIND_STRESS_X_PA.to_bits());
        }
        for tau in frame.tau_y() {
            assert_eq!(tau.to_bits(), TRADE_WIND_STRESS_Y_PA.to_bits());
        }
    }

    // Frames are in the order they were written, and model time advances by
    // exactly one output interval between them.
    let interval_s = header.output.interval_s;
    for (earlier, later) in frames.iter().zip(frames.iter().skip(1)) {
        let gap = later.t_s() - earlier.t_s();
        // Exact to a few ulps: both times are `step · dt` with `step` a whole
        // multiple of the cadence, so the difference is exact but for the
        // rounding of one subtraction (ε ≈ 2.2e-16).
        assert!(
            (gap - interval_s).abs() <= 8.0 * f64::EPSILON * interval_s,
            "frames are one output interval apart: {gap} vs {interval_s}"
        );
    }
}

#[test]
fn the_header_is_complete_before_the_first_frame_is_appended() {
    // The writer writes the header once, when the run opens (T-05.2), so a
    // reader watching a run in progress — or picking over one that crashed —
    // finds a whole header rather than a placeholder to be patched at the end.
    let scratch = ScratchDir::new("header-first");
    let (grid, spacing) = basin(NX, NY);
    let mut writer = RunWriter::create(scratch.path(), &header(spacing, NX, NY))
        .expect("the scratch directory is writable");

    let at_open =
        fs::read(scratch.path().join(HEADER_FILE_NAME)).expect("the header file was written");
    let parsed = read_header(&at_open);
    assert_eq!(parsed.output.frame_count, EXPECTED_FRAMES);
    assert_eq!(parsed, header(spacing, NX, NY));

    write_short_run(&mut writer, grid, spacing);
    writer
        .finish()
        .expect("the run wrote every frame it promised");

    let at_close =
        fs::read(scratch.path().join(HEADER_FILE_NAME)).expect("the header file is still there");
    assert_eq!(
        at_open, at_close,
        "the header is written once, not rewritten"
    );
}

#[test]
fn the_header_records_the_parameters_the_run_was_integrated_with() {
    // A header that disagreed with the run would mislabel every archived
    // frame; the units carried across are the SI ones both crates state.
    let (_, spacing) = basin(NX, NY);
    let recorded = header(spacing, NX, NY).physical_params;
    let engine_params = params();

    assert_eq!(
        recorded.mean_depth_m.to_bits(),
        engine_params.mean_thermocline_depth_m().to_bits()
    );
    assert_eq!(
        recorded.reduced_gravity_m_per_s2.to_bits(),
        engine_params.reduced_gravity_m_per_s2().to_bits()
    );
    assert_eq!(
        recorded.beta_per_m_per_s.to_bits(),
        engine_params.beta_per_m_per_s().to_bits()
    );
    assert_eq!(
        recorded.rayleigh_damping_per_s.to_bits(),
        engine_params.rayleigh_damping_per_s().to_bits()
    );
    assert_eq!(
        recorded.reference_density_kg_per_m3.to_bits(),
        engine_params.reference_density_kg_per_m3().to_bits()
    );
}

#[test]
fn a_run_written_to_memory_is_the_same_run_written_to_files() {
    // ADR-0006: the reader is defined over a byte source, so the writer must
    // be defined over a byte sink — a browser has no filesystem, and a run
    // held in memory has to be the same bytes as one on disk.
    let scratch = ScratchDir::new("byte-sink");
    write_short_run_to_files(scratch.path());
    let (header_bytes, frame_bytes) = write_short_run_to_memory();

    assert_eq!(
        fs::read(scratch.path().join(HEADER_FILE_NAME)).expect("the header file was written"),
        header_bytes
    );
    assert_eq!(
        fs::read(scratch.path().join(FRAME_FILE_NAME)).expect("the frame file was written"),
        frame_bytes
    );
}

#[test]
fn the_frame_file_is_a_plain_concatenation_of_encoded_frames() {
    // What makes the format readable from a byte source front to back: every
    // frame is the same size on this grid, they are written one after another,
    // and there is no index or trailer to seek to.
    //
    // The size is bincode's, derived from its specification rather than
    // measured: with the standard configuration a struct is its fields in
    // order, an `f64` is 8 fixed bytes, and a sequence is a variable-length
    // length prefix — one byte for a length below 251 — followed by its
    // elements. So a frame is `t` (8) plus, for each of the five fields, one
    // length byte plus 8 bytes per point: h over 6×4 cells, u and tau_x over
    // 7×4 faces, v and tau_y over 6×5 faces.
    let grid = grid_spec(NX, NY);
    let points: usize = Variable::ALL.iter().map(|v| grid.field_len(*v)).sum();
    // 8 bytes of `t`, one length byte per field, and 8 bytes per point.
    const T_BYTES: usize = 8;
    const LENGTH_PREFIX_BYTES: usize = 1;
    const F64_BYTES: usize = 8;
    let expected_frame_bytes =
        T_BYTES + Variable::ALL.len() * LENGTH_PREFIX_BYTES + points * F64_BYTES;
    assert_eq!(expected_frame_bytes, 8 + 5 + (24 + 28 + 30 + 28 + 30) * 8);

    let (header_bytes, frame_bytes) = write_short_run_to_memory();
    let header = read_header(&header_bytes);
    assert_eq!(
        frame_bytes.len(),
        expected_frame_bytes * header.output.frame_count as usize
    );
}

#[test]
fn the_same_scenario_written_twice_is_byte_identical() {
    // CODING_STANDARDS.md § Correctness and failure: identical scenario in,
    // byte-identical output.
    let first = write_short_run_to_memory();
    let second = write_short_run_to_memory();
    assert_eq!(first.0, second.0, "the headers differ");
    assert_eq!(first.1, second.1, "the frames differ");
}

#[test]
fn a_state_that_does_not_fit_the_header_grid_is_rejected_by_name() {
    // Invalid input returns a `Result` naming the offending value and the
    // bound it violated (CODING_STANDARDS.md), rather than writing a run whose
    // frames disagree with its own header.
    let (grid, spacing) = basin(NX, NY);
    // A header for a wider basin than the state the caller then hands over.
    let mut writer = RunWriter::new(Vec::new(), Vec::new(), &header(spacing, NX + 1, NY))
        .expect("a vector never fails a write");

    let wind = WindStressField::uniform_including_walls(
        grid,
        TRADE_WIND_STRESS_X_PA,
        TRADE_WIND_STRESS_Y_PA,
    );
    let error = writer
        .append(0.0, &OceanState::at_rest(grid), &wind)
        .expect_err("a 6-column state does not cover a 7-column basin");

    assert!(
        matches!(
            &error,
            RunWriteError::Frame(FormatError::FieldShape {
                variable: Variable::ThermoclineDepthAnomaly,
                expected,
                actual,
            }) if *expected == (NX + 1) * NY && *actual == NX * NY
        ),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.contains('h'), "{message}");
    assert!(message.contains("expected 28"), "{message}");
}

#[test]
fn appending_more_frames_than_the_header_promises_is_rejected() {
    // The header's frame count is what a byte-source reader trusts; a frame
    // past it would be a frame no reader ever sees, so the writer refuses it
    // rather than writing an unreadable tail.
    let (grid, spacing) = basin(NX, NY);
    let mut writer = RunWriter::new(Vec::new(), Vec::new(), &header(spacing, NX, NY))
        .expect("a vector never fails a write");
    let wind = WindStressField::uniform_including_walls(
        grid,
        TRADE_WIND_STRESS_X_PA,
        TRADE_WIND_STRESS_Y_PA,
    );
    let state = OceanState::at_rest(grid);

    for _ in 0..EXPECTED_FRAMES {
        writer
            .append(0.0, &state, &wind)
            .expect("the header promised this many frames");
    }
    let error = writer
        .append(0.0, &state, &wind)
        .expect_err("one frame past the count the header promised");

    assert!(
        matches!(error, RunWriteError::TooManyFrames { promised } if promised == EXPECTED_FRAMES),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.contains(&EXPECTED_FRAMES.to_string()), "{message}");
}

#[test]
fn finishing_a_run_short_of_its_frame_count_is_rejected() {
    // The mirror image: a run that stopped early leaves a header promising
    // frames the file does not hold, which a reader would read off the end of.
    let (grid, spacing) = basin(NX, NY);
    let mut writer = RunWriter::new(Vec::new(), Vec::new(), &header(spacing, NX, NY))
        .expect("a vector never fails a write");
    let wind = WindStressField::uniform_including_walls(
        grid,
        TRADE_WIND_STRESS_X_PA,
        TRADE_WIND_STRESS_Y_PA,
    );
    writer
        .append(0.0, &OceanState::at_rest(grid), &wind)
        .expect("the header promised four frames");

    let error = writer
        .finish()
        .expect_err("one frame written of the four the header promised");
    assert!(
        matches!(
            error,
            RunWriteError::MissingFrames { promised, written }
                if promised == EXPECTED_FRAMES && written == 1
        ),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.contains('4'), "{message}");
    assert!(message.contains('1'), "{message}");
}

#[test]
fn a_schedule_with_no_output_cadence_is_rejected() {
    // A cadence of zero steps has no meaning and would divide by zero on the
    // way to a frame count; a non-positive timestep is not a run.
    let error = OutputSchedule::new(1.0, TOTAL_STEPS, 0)
        .expect_err("a run cannot write a frame every zero steps");
    assert_eq!(error, OutputScheduleError::CadenceIsZero);
    assert!(error.to_string().contains("every_n_steps"), "{error}");

    for bad_dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = OutputSchedule::new(bad_dt, TOTAL_STEPS, EVERY_N_STEPS)
            .expect_err("a timestep must be a finite, positive duration");
        // Compared as bit patterns: NaN is one of the rejected values, and it
        // is not equal to itself.
        assert!(
            matches!(
                error,
                OutputScheduleError::TimestepNotPositive { dt_s }
                    if dt_s.to_bits() == bad_dt.to_bits()
            ),
            "{bad_dt} was rejected as something else"
        );
    }
}
