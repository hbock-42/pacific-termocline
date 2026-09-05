//! Acceptance tests for T-05.3 — the reader API.
//!
//! The reader is exercised here against runs assembled *by hand* out of
//! `serde_json` and `bincode`, so nothing in this file learns what the format
//! is by asking the reader. The companion test in `engine/tests/run_reader.rs`
//! closes the loop the acceptance criteria actually name — a run produced by
//! the T-05.2 writer, read back frame for frame — and the memory bound lives
//! in `reader_memory.rs`.
//!
//! Nothing here has a tolerance. Reading a run back is an identity, not an
//! approximation, so every field is compared as an IEEE-754 bit pattern.
//!
//! Every byte source below is an in-memory `Cursor`: per [ADR-0006] the reader
//! must work where there is no filesystem, and a test suite that only ever
//! hands it files would not notice the day it stopped doing so.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::io::Cursor;

use termocline_format::{
    frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader,
    RunReadError, RunReader, Variable, FORMAT_VERSION,
};

/// A basin small enough to write out by hand.
const NX: usize = 3;
const NY: usize = 2;

/// Frames in the fixture run. Small, and not a power of two, so an off-by-one
/// in the count is visible.
const FRAME_COUNT: u64 = 5;

/// Model time between frames, in seconds: one day of output cadence.
const INTERVAL_S: f64 = 86_400.0;

fn grid() -> GridSpec {
    // The equatorial Pacific basin of CONTEXT.md: 120°E-80°W, 25°S-25°N.
    GridSpec::new(NX, NY, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
        .expect("a 3x2 basin is a valid grid")
}

fn header(frame_count: u64) -> RunHeader {
    // Scenario values, not physical constants: they exist to be read back.
    let params = PhysicalParams {
        mean_depth_m: 150.0,
        reduced_gravity_m_per_s2: 0.05,
        beta_per_m_per_s: 2.28e-11,
        rayleigh_damping_per_s: 1.0e-6,
        reference_density_kg_per_m3: 1025.0,
    };
    RunHeader::new(
        grid(),
        params,
        "a run assembled by hand",
        OutputTiming {
            frame_count,
            interval_s: INTERVAL_S,
        },
    )
}

/// Values chosen to stress the decoder rather than to be physical: a negative
/// zero, a subnormal and the extremes of f64 all have bit patterns a lossy
/// decoder would not return unchanged. `index` shifts the pattern so two
/// frames, or two fields of one frame, are never the same buffer.
fn awkward_values(len: usize, index: usize) -> Vec<f64> {
    let stressors = [
        -0.0_f64,
        f64::MIN_POSITIVE / 2.0,
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        -1.234_567_890_123_456_7e-9,
    ];
    (0..len)
        .map(|i| stressors[(i + index) % stressors.len()])
        .collect()
}

/// Frame `index` of the fixture run, at model time `index * INTERVAL_S`.
fn frame(index: u64) -> Frame {
    let g = grid();
    let offset = index as usize;
    Frame::new(
        index as f64 * INTERVAL_S,
        &g,
        awkward_values(g.field_len(Variable::ThermoclineDepthAnomaly), offset),
        awkward_values(g.field_len(Variable::ZonalCurrentAnomaly), offset + 1),
        awkward_values(g.field_len(Variable::MeridionalCurrentAnomaly), offset + 2),
        awkward_values(g.field_len(Variable::ZonalWindStress), offset + 3),
        awkward_values(g.field_len(Variable::MeridionalWindStress), offset + 4),
    )
    .expect("every field is built at the length the grid asks for")
}

/// The two byte sources of a run: the JSON header, and `frame_count` frames
/// concatenated with nothing between them.
fn run_bytes(frame_count: u64) -> (Vec<u8>, Vec<u8>) {
    let header_bytes = serde_json::to_vec(&header(frame_count)).expect("the header serializes");
    let mut frame_bytes = Vec::new();
    for index in 0..frame_count {
        bincode::serde::encode_into_std_write(frame(index), &mut frame_bytes, frame_encoding())
            .expect("a vector never fails a write");
    }
    (header_bytes, frame_bytes)
}

fn reader(header_bytes: Vec<u8>, frame_bytes: Vec<u8>) -> RunReader<Cursor<Vec<u8>>> {
    RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the fixture run has a readable header")
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

fn assert_frames_identical(index: u64, written: &Frame, read: &Frame) {
    assert_eq!(
        written.t_s().to_bits(),
        read.t_s().to_bits(),
        "frame {index}: model time changed"
    );
    for variable in Variable::ALL {
        assert_bit_identical(
            &format!("frame {index} field {}", variable.symbol()),
            written.field(variable),
            read.field(variable),
        );
    }
}

#[test]
fn the_header_is_available_before_a_single_frame_is_read() {
    // The deliverable is an iterator *plus* a header accessor: a caller sizes
    // its buffers and labels its axes from the header, which means reading it
    // must not cost a frame.
    let (header_bytes, frame_bytes) = run_bytes(FRAME_COUNT);
    let reader = reader(header_bytes, frame_bytes);

    assert_eq!(reader.header(), &header(FRAME_COUNT));
    assert_eq!(reader.header().output.frame_count, FRAME_COUNT);
    assert_eq!(reader.remaining_frames(), FRAME_COUNT);
}

#[test]
fn every_frame_comes_back_exactly_as_it_was_written_in_order() {
    // The acceptance criterion: frames identical to what was written, in
    // order. The expected values are the fixture frames, built independently
    // of the reader.
    let (header_bytes, frame_bytes) = run_bytes(FRAME_COUNT);
    let read: Vec<Frame> = reader(header_bytes, frame_bytes)
        .collect::<Result<_, _>>()
        .expect("the fixture run reads back");

    assert_eq!(read.len(), FRAME_COUNT as usize);
    for (index, read) in read.iter().enumerate() {
        assert_frames_identical(index as u64, &frame(index as u64), read);
    }
}

#[test]
fn the_frame_count_shrinks_as_the_run_is_consumed() {
    // A caller that streams frames wants to know how many are left without
    // counting them itself; the count comes from the header, so it is known
    // before the frames are read.
    let (header_bytes, frame_bytes) = run_bytes(FRAME_COUNT);
    let mut reader = reader(header_bytes, frame_bytes);

    for expected_remaining in (0..FRAME_COUNT).rev() {
        let frame = reader
            .next()
            .expect("a promised frame is there")
            .expect("it decodes");
        assert_eq!(
            frame.t_s(),
            (FRAME_COUNT - 1 - expected_remaining) as f64 * INTERVAL_S
        );
        assert_eq!(reader.remaining_frames(), expected_remaining);
    }
    assert!(reader.next().is_none(), "the run ends after its last frame");
}

#[test]
fn a_run_of_one_frame_reads_back() {
    // The initial state alone is a legal run: the writer's schedule always
    // saves step 0, so `frame_count` is never zero and one is its floor.
    let (header_bytes, frame_bytes) = run_bytes(1);
    let read: Vec<Frame> = reader(header_bytes, frame_bytes)
        .collect::<Result<_, _>>()
        .expect("a one-frame run reads back");

    assert_eq!(read.len(), 1);
    assert_frames_identical(0, &frame(0), &read[0]);
}

#[test]
fn a_run_written_in_a_future_format_version_is_refused_by_name() {
    // The header carries a version so a reader can tell whether it understands
    // a file rather than guessing; decoding frames of an unknown layout would
    // produce plausible garbage.
    let (header_bytes, frame_bytes) = run_bytes(FRAME_COUNT);
    let json = String::from_utf8(header_bytes).expect("the header is UTF-8 JSON");
    let future = json.replace(
        &format!("\"format_version\":{FORMAT_VERSION}"),
        &format!("\"format_version\":{}", FORMAT_VERSION + 1),
    );
    assert_ne!(json, future, "the fixture must actually carry a version");

    let err = RunReader::new(Cursor::new(future.into_bytes()), Cursor::new(frame_bytes))
        .expect_err("a run from the future is not readable");
    assert!(
        matches!(
            err,
            RunReadError::UnsupportedVersion { found, supported }
                if found == FORMAT_VERSION + 1 && supported == FORMAT_VERSION
        ),
        "{err:?}"
    );
    let message = err.to_string();
    assert!(
        message.contains(&format!("{}", FORMAT_VERSION + 1)),
        "{message}"
    );
}

#[test]
fn a_run_whose_frames_stop_early_is_refused_rather_than_ending_quietly() {
    // A run cut short by a crash has a header promising more frames than the
    // file holds. Ending the iteration quietly would present a truncated run
    // as a complete one.
    let (header_bytes, frame_bytes) = run_bytes(FRAME_COUNT);
    let truncated = frame_bytes[..frame_bytes.len() / 2].to_vec();

    let mut reader = reader(header_bytes, truncated);
    let error = loop {
        match reader.next() {
            Some(Ok(_)) => {}
            Some(Err(error)) => break error,
            None => panic!("a truncated run must not end as if it were complete"),
        }
    };
    assert!(
        matches!(
            error,
            RunReadError::Truncated {
                promised: FRAME_COUNT,
                ..
            }
        ),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.contains("5"), "{message}");
}

#[test]
fn a_frame_file_longer_than_its_header_promises_is_refused() {
    // The mirror of the writer's `TooManyFrames`: bytes past the promised
    // count mean the header and the frames disagree, and the reader says so
    // rather than dropping data on the floor.
    let (header_bytes, mut frame_bytes) = run_bytes(FRAME_COUNT);
    bincode::serde::encode_into_std_write(frame(FRAME_COUNT), &mut frame_bytes, frame_encoding())
        .expect("a vector never fails a write");

    let mut reader = reader(header_bytes, frame_bytes);
    for _ in 0..FRAME_COUNT {
        reader
            .next()
            .expect("a promised frame is there")
            .expect("it decodes");
    }
    let error = reader
        .next()
        .expect("the extra frame is noticed")
        .expect_err("it is refused");
    assert!(
        matches!(error, RunReadError::TrailingBytes { promised } if promised == FRAME_COUNT),
        "{error:?}"
    );
}

#[test]
fn a_frame_that_does_not_fit_the_headers_grid_is_refused_by_name() {
    // The frame bytes do not record the grid, so a frame from another basin
    // decodes happily; only the header says how long each field should be.
    let (header_bytes, _) = run_bytes(1);
    let other = GridSpec::new(NX + 1, NY, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
        .expect("a 4x2 basin is a valid grid");
    let foreign = Frame::new(
        0.0,
        &other,
        vec![0.0; other.field_len(Variable::ThermoclineDepthAnomaly)],
        vec![0.0; other.field_len(Variable::ZonalCurrentAnomaly)],
        vec![0.0; other.field_len(Variable::MeridionalCurrentAnomaly)],
        vec![0.0; other.field_len(Variable::ZonalWindStress)],
        vec![0.0; other.field_len(Variable::MeridionalWindStress)],
    )
    .expect("the frame fits the grid it was built on");
    let mut frame_bytes = Vec::new();
    bincode::serde::encode_into_std_write(foreign, &mut frame_bytes, frame_encoding())
        .expect("a vector never fails a write");

    let error = reader(header_bytes, frame_bytes)
        .next()
        .expect("the frame is there")
        .expect_err("it does not fit the header's basin");
    let message = error.to_string();
    assert!(matches!(error, RunReadError::Frame(_)), "{error:?}");
    assert!(message.contains('h'), "{message}");
}

#[test]
fn a_header_that_is_not_json_is_refused() {
    // Invalid input returns a Result with an actionable message, rather than
    // panicking (CODING_STANDARDS.md).
    let error = RunReader::new(
        Cursor::new(b"not a header".to_vec()),
        Cursor::new(Vec::new()),
    )
    .expect_err("a run needs a header");
    assert!(matches!(error, RunReadError::Header(_)), "{error:?}");
}
