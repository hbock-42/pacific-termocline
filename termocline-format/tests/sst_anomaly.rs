//! The acceptance criteria of T-05.4, at the level of the format itself: the
//! mixed-layer SST anomaly `T'` travels in a frame, an absent one is absent
//! rather than zero, and a run written by format version 1 still reads.
//!
//! Nothing here has a tolerance. A frame is a transcription, not an
//! approximation — a format that returned "close enough" values would corrupt
//! every archived run — so the assertions compare exact IEEE-754 bit patterns.
//!
//! The version 1 fixture below is *not* produced by this build. Version 1's
//! frame layout is written out here as its own struct, from the specification
//! the old format had — model time and the five fields of the linear core, and
//! nothing after them — so the test reads bytes this build could not have
//! written, which is the only thing that proves an archive still opens.

use std::io::Cursor;

use serde::Serialize;
use termocline_format::{
    frame_encoding, BasinExtent, FormatError, Frame, GridSpec, OutputTiming, PhysicalParams,
    RunHeader, RunReadError, RunReader, Variable, FORMAT_VERSION,
};

/// A basin small enough to write out by hand, so every expected value below
/// comes from the fixture rather than from running the code.
const NX: usize = 3;
const NY: usize = 2;

/// Model time between the fixture's frames, in seconds: a day.
const INTERVAL_S: f64 = 86_400.0;

fn grid() -> GridSpec {
    // The equatorial Pacific basin of CONTEXT.md: 120°E-80°W, 25°S-25°N.
    GridSpec::new(NX, NY, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
        .expect("a 3x2 basin is a valid grid")
}

fn params() -> PhysicalParams {
    PhysicalParams {
        // The 1.5-layer equatorial Pacific of CONTEXT.md: c = sqrt(g'H) = 3 m/s.
        mean_depth_m: 150.0,
        reduced_gravity_m_per_s2: 0.06,
        beta_per_m_per_s: 2.3e-11,
        rayleigh_damping_per_s: 1.0e-7,
        reference_density_kg_per_m3: 1025.0,
    }
}

fn header(frame_count: u64) -> RunHeader {
    RunHeader::new(
        grid(),
        params(),
        "a run assembled by hand",
        OutputTiming {
            frame_count,
            interval_s: INTERVAL_S,
        },
    )
}

/// Values chosen to stress the encoder rather than to be physical: a negative
/// zero, a subnormal and the extremes of f64 all have bit patterns a lossy
/// encoder would not return unchanged. `offset` shifts the pattern so two
/// fields of one frame are never the same buffer.
fn awkward_values(len: usize, offset: usize) -> Vec<f64> {
    let stressors = [
        -0.0_f64,
        f64::MIN_POSITIVE / 2.0,
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        -1.234_567_890_123_456_7e-9,
    ];
    (0..len)
        .map(|i| stressors[(i + offset) % stressors.len()])
        .collect()
}

/// A frame of the linear core at model time `t_s`.
fn core_frame(t_s: f64) -> Frame {
    let g = grid();
    Frame::new(
        t_s,
        &g,
        awkward_values(g.field_len(Variable::ThermoclineDepthAnomaly), 0),
        awkward_values(g.field_len(Variable::ZonalCurrentAnomaly), 1),
        awkward_values(g.field_len(Variable::MeridionalCurrentAnomaly), 2),
        awkward_values(g.field_len(Variable::ZonalWindStress), 3),
        awkward_values(g.field_len(Variable::MeridionalWindStress), 4),
    )
    .expect("every field is built at the length the grid asks for")
}

/// The same frame with `T'` on it, as a coupled run writes one.
fn coupled_frame(t_s: f64) -> Frame {
    let g = grid();
    core_frame(t_s)
        .with_sst_anomaly(&g, awkward_values(g.field_len(Variable::SstAnomaly), 5))
        .expect("the anomaly is built at the length the grid asks for")
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

/// A run's two byte sources: the JSON header, and `frames` encoded one after
/// another with nothing between them.
fn run_bytes(header: &RunHeader, frames: &[Frame]) -> (Vec<u8>, Vec<u8>) {
    let header_bytes = serde_json::to_vec(header).expect("the header is serializable");
    let mut frame_bytes = Vec::new();
    for frame in frames {
        bincode::serde::encode_into_std_write(frame, &mut frame_bytes, frame_encoding())
            .expect("the frame encodes");
    }
    (header_bytes, frame_bytes)
}

fn read_back(header_bytes: Vec<u8>, frame_bytes: Vec<u8>) -> (RunHeader, Vec<Frame>) {
    let mut reader = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the buffers hold a readable run");
    let header = reader.header().clone();
    let frames = reader
        .by_ref()
        .collect::<Result<Vec<_>, _>>()
        .expect("every frame of the fixture run decodes");
    (header, frames)
}

#[test]
fn a_coupled_run_round_trips_the_sst_anomaly_bit_for_bit() {
    // The ticket's first criterion. `T'` is written, read back, and is the
    // same 24 bit patterns it went in as.
    let written: Vec<Frame> = (0..3)
        .map(|k| coupled_frame(k as f64 * INTERVAL_S))
        .collect();
    let header = header(written.len() as u64).with_sst_anomaly();
    let (header_bytes, frame_bytes) = run_bytes(&header, &written);

    let (read_header, read) = read_back(header_bytes, frame_bytes);
    assert_eq!(read_header, header);
    assert_eq!(read.len(), written.len());
    for (index, (written, read)) in written.iter().zip(&read).enumerate() {
        assert_eq!(written.t_s().to_bits(), read.t_s().to_bits());
        let written_sst = written
            .sst_anomaly_k()
            .expect("the fixture frame carries an anomaly");
        let read_sst = read
            .sst_anomaly_k()
            .expect("a coupled run's frame carries its anomaly");
        assert_bit_identical(&format!("frame {index} sst"), written_sst, read_sst);
        // And the linear core is untouched by the extra field beside it.
        assert_bit_identical(&format!("frame {index} h"), written.h(), read.h());
        assert_bit_identical(&format!("frame {index} u"), written.u(), read.u());
        assert_bit_identical(&format!("frame {index} v"), written.v(), read.v());
    }
}

#[test]
fn a_coupled_run_declares_the_anomaly_in_its_header_with_its_unit() {
    // A reader never guesses at the meaning of a field (ADR-0004), so the
    // variable list gains an entry rather than the frame gaining a field
    // nothing announced. The unit is kelvin: `T'` is a temperature difference.
    let coupled = header(1).with_sst_anomaly();
    let symbols: Vec<&str> = coupled
        .variables
        .iter()
        .map(|v| v.symbol.as_str())
        .collect();
    let units: Vec<&str> = coupled.variables.iter().map(|v| v.unit.as_str()).collect();

    assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y", "sst"]);
    assert_eq!(units, ["m", "m s^-1", "m s^-1", "N m^-2", "N m^-2", "K"]);
    assert!(coupled.carries(Variable::SstAnomaly));
    // The core's five entries are exactly an uncoupled run's, in order.
    assert_eq!(coupled.variables[..5], header(1).variables[..]);
}

#[test]
fn an_uncoupled_run_reports_its_missing_anomaly_as_absent_rather_than_zero() {
    // The ticket's trap. A basin of zeros would round-trip perfectly and would
    // claim the ocean sat at exactly its climatological temperature; absence
    // says the run has no `T'`, which is the true statement about it.
    let written: Vec<Frame> = (0..3).map(|k| core_frame(k as f64 * INTERVAL_S)).collect();
    let header = header(written.len() as u64);
    let (header_bytes, frame_bytes) = run_bytes(&header, &written);

    let (read_header, read) = read_back(header_bytes, frame_bytes);
    assert!(!read_header.carries(Variable::SstAnomaly));
    for (index, frame) in read.iter().enumerate() {
        assert_eq!(frame.sst_anomaly_k(), None, "frame {index}");
        assert_eq!(frame.field(Variable::SstAnomaly), None, "frame {index}");
        assert!(!frame.carries_sst_anomaly(), "frame {index}");
    }
}

#[test]
fn an_uncoupled_frame_pays_one_byte_for_the_field_it_does_not_have() {
    // "An uncoupled run does not pay for a field it does not have." The price
    // is bincode's one-byte `Option` tag, derived from its specification: a
    // `None` is the tag and nothing else. Measured against the version 1
    // frame, which is this frame with no place for `T'` at all, so the
    // difference is exactly what version 2 costs a run that has none.
    const OPTION_TAG_BYTES: usize = 1;
    const LENGTH_PREFIX_BYTES: usize = 1;
    const F64_BYTES: usize = 8;

    let (_, uncoupled) = run_bytes(&header(1), &[core_frame(0.0)]);
    let (_, version_1) = version_1_run(1);
    assert_eq!(uncoupled.len(), version_1.len() + OPTION_TAG_BYTES);

    // A field of zeros would have cost the tag, a length byte and one f64 a
    // cell — 26 bytes on this 3x2 basin against 1, and it would have been a
    // temperature nobody computed.
    let fabricated = core_frame(0.0)
        .with_sst_anomaly(&grid(), vec![0.0; NX * NY])
        .expect("a zero anomaly is the right length");
    let (_, fabricated) = run_bytes(&header(1).with_sst_anomaly(), &[fabricated]);
    assert_eq!(
        fabricated.len(),
        version_1.len() + OPTION_TAG_BYTES + LENGTH_PREFIX_BYTES + NX * NY * F64_BYTES
    );

    // And the anomaly a coupled run really has costs the same as the
    // fabricated one: absence is the only thing being saved.
    let (_, coupled) = run_bytes(&header(1).with_sst_anomaly(), &[coupled_frame(0.0)]);
    assert_eq!(coupled.len(), fabricated.len());
}

#[test]
fn a_frame_carrying_a_variable_its_header_never_declared_is_refused() {
    // The header's variable list is what a reader indexes a run by, so the
    // reader holds the frames to it as well as to the grid. A run offering a
    // `T'` its header never announced is a run whose frames and whose labels
    // disagree, and reading it would mean reporting a field under a heading
    // that does not exist.
    let (header_bytes, frame_bytes) = run_bytes(&header(1), &[coupled_frame(0.0)]);

    let mut reader = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the header itself is valid");
    let error = reader
        .next()
        .expect("the header promises a frame")
        .expect_err("an undeclared anomaly is not read as data");

    assert!(
        matches!(
            error,
            RunReadError::Frame(FormatError::UndeclaredVariable {
                variable: Variable::SstAnomaly,
                declared: false,
            })
        ),
        "{error:?}"
    );
    assert!(error.to_string().contains("sst"), "{error}");
}

#[test]
fn a_frame_missing_a_variable_its_header_declared_is_refused() {
    // The other direction, and the one the ticket's trap is about: a header
    // promising `T'` beside frames that have none would leave a reader with a
    // declared variable and nothing behind it — which is exactly where a
    // buffer of zeros would get invented to fill the hole.
    let (header_bytes, frame_bytes) = run_bytes(&header(1).with_sst_anomaly(), &[core_frame(0.0)]);

    let mut reader = RunReader::new(Cursor::new(header_bytes), Cursor::new(frame_bytes))
        .expect("the header itself is valid");
    let error = reader
        .next()
        .expect("the header promises a frame")
        .expect_err("a promised anomaly that is not there is not read as absent");

    assert!(
        matches!(
            error,
            RunReadError::Frame(FormatError::UndeclaredVariable {
                variable: Variable::SstAnomaly,
                declared: true,
            })
        ),
        "{error:?}"
    );
}

#[test]
fn an_anomaly_that_does_not_cover_the_basin_is_refused_by_name() {
    // Invalid input returns a Result naming the offending value and the bound
    // it violated (CODING_STANDARDS.md), rather than writing a short field.
    let error = core_frame(0.0)
        .with_sst_anomaly(&grid(), vec![0.0; NX * NY - 1])
        .expect_err("an anomaly one cell short does not cover the basin");

    assert_eq!(
        error,
        FormatError::FieldShape {
            variable: Variable::SstAnomaly,
            expected: NX * NY,
            actual: NX * NY - 1,
        }
    );
    assert!(error.to_string().contains("sst"), "{error}");
}

// --- Runs written before this change -----------------------------------

/// Format version 1's frame, as that version specified it: model time and the
/// five fields of the linear core, with nothing after them. Serialize only —
/// this is a fixture writer, and nothing in the build writes version 1.
#[derive(Serialize)]
struct FrameV1 {
    #[serde(rename = "t")]
    t_s: f64,
    h: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    tau_x: Vec<f64>,
    tau_y: Vec<f64>,
}

impl FrameV1 {
    /// The version 1 frame carrying the same values as [`core_frame`], so the
    /// two layouts can be compared value for value.
    fn of(t_s: f64) -> Self {
        let g = grid();
        Self {
            t_s,
            h: awkward_values(g.field_len(Variable::ThermoclineDepthAnomaly), 0),
            u: awkward_values(g.field_len(Variable::ZonalCurrentAnomaly), 1),
            v: awkward_values(g.field_len(Variable::MeridionalCurrentAnomaly), 2),
            tau_x: awkward_values(g.field_len(Variable::ZonalWindStress), 3),
            tau_y: awkward_values(g.field_len(Variable::MeridionalWindStress), 4),
        }
    }
}

/// A whole run as version 1 wrote it: a header stamped `1`, listing the five
/// variables of the linear core, beside frames in the version 1 layout.
fn version_1_run(frame_count: u64) -> (Vec<u8>, Vec<u8>) {
    let mut header = header(frame_count);
    header.format_version = 1;
    let header_bytes = serde_json::to_vec(&header).expect("the header is serializable");

    let mut frame_bytes = Vec::new();
    for index in 0..frame_count {
        bincode::serde::encode_into_std_write(
            FrameV1::of(index as f64 * INTERVAL_S),
            &mut frame_bytes,
            frame_encoding(),
        )
        .expect("the version 1 frame encodes");
    }
    (header_bytes, frame_bytes)
}

#[test]
fn a_run_written_before_this_change_still_reads() {
    // The ticket's second criterion, and ADR-0011's decision: a version 1 run
    // is a complete run of the linear core, so it opens rather than being
    // refused for a field it never had.
    const FRAME_COUNT: u64 = 3;
    let (header_bytes, frame_bytes) = version_1_run(FRAME_COUNT);

    let (read_header, frames) = read_back(header_bytes, frame_bytes);

    assert_eq!(read_header.format_version, 1);
    assert_ne!(read_header.format_version, FORMAT_VERSION);
    assert_eq!(frames.len() as u64, FRAME_COUNT);
    for (index, frame) in frames.iter().enumerate() {
        // The frames really were decoded one after another, not resynchronized
        // by luck: each carries its own model time and its own values.
        let expected = FrameV1::of(index as f64 * INTERVAL_S);
        assert_eq!(
            frame.t_s().to_bits(),
            expected.t_s.to_bits(),
            "frame {index}"
        );
        assert_bit_identical(&format!("frame {index} h"), &expected.h, frame.h());
        assert_bit_identical(&format!("frame {index} u"), &expected.u, frame.u());
        assert_bit_identical(&format!("frame {index} v"), &expected.v, frame.v());
        assert_bit_identical(
            &format!("frame {index} tau_x"),
            &expected.tau_x,
            frame.tau_x(),
        );
        assert_bit_identical(
            &format!("frame {index} tau_y"),
            &expected.tau_y,
            frame.tau_y(),
        );
    }
}

#[test]
fn a_run_written_before_this_change_has_no_sst_anomaly_rather_than_a_zero_one() {
    // The honest half of reading an old archive: version 1 had no `T'` at all,
    // so a reader is told there is none. Zeros would be a temperature this
    // build invented for a run that never computed one.
    let (header_bytes, frame_bytes) = version_1_run(2);
    let (read_header, frames) = read_back(header_bytes, frame_bytes);

    assert!(!read_header.carries(Variable::SstAnomaly));
    let symbols: Vec<&str> = read_header
        .variables
        .iter()
        .map(|v| v.symbol.as_str())
        .collect();
    assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y"]);
    for (index, frame) in frames.iter().enumerate() {
        assert_eq!(frame.sst_anomaly_k(), None, "frame {index}");
    }
}

#[test]
fn a_version_1_frame_can_be_fetched_at_random_out_of_a_run_held_whole() {
    // The scrubber's path into an old archive (T-08.3, ADR-0006): the reader's
    // random-access half decodes by the header's version too, so the version
    // that reaches a browser tab is the version it decodes with.
    const FRAME_COUNT: u64 = 4;
    let (header_bytes, frame_bytes) = version_1_run(FRAME_COUNT);
    let header: RunHeader =
        serde_json::from_slice(&header_bytes).expect("the fixture header is JSON");

    let mut offset = 0;
    for index in 0..FRAME_COUNT {
        let (frame, used) = termocline_format::decode_frame(&frame_bytes[offset..], &header)
            .expect("every version 1 frame decodes");
        assert_eq!(
            frame.t_s().to_bits(),
            (index as f64 * INTERVAL_S).to_bits(),
            "frame {index}"
        );
        assert_eq!(frame.sst_anomaly_k(), None, "frame {index}");
        offset += used;
    }
    assert_eq!(
        offset,
        frame_bytes.len(),
        "the reported lengths should account for every byte of the run"
    );
}

#[test]
fn a_version_0_run_is_refused_by_name() {
    // The range is bounded at both ends. Version 0 was never written, so there
    // is no layout for it and the reader says so rather than guessing.
    let mut header = header(1);
    header.format_version = 0;
    let header_bytes = serde_json::to_vec(&header).expect("the header is serializable");

    let error = RunReader::new(Cursor::new(header_bytes), Cursor::new(Vec::new()))
        .expect_err("version 0 is not a version this build reads");
    assert!(
        matches!(error, RunReadError::UnsupportedVersion { found: 0, .. }),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.contains('0'), "{message}");
}
