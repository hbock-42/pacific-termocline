//! The acceptance criteria of T-05.1: a header and a frame each survive a
//! serialize -> deserialize round trip losslessly.
//!
//! Nothing here has a tolerance. Serialization is an identity, not an
//! approximation: a format that returns "close enough" floats would silently
//! corrupt every archived run, so the assertions compare exact IEEE-754 bit
//! patterns rather than values within a bound.

use termocline_format::{
    BasinExtent, FormatError, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader, Variable,
    FORMAT_VERSION,
};

/// A basin small enough to write out by hand, so every expected value below
/// comes from the fixture rather than from running the code.
const NX: usize = 3;
const NY: usize = 2;

fn extent() -> BasinExtent {
    // The equatorial Pacific basin of CONTEXT.md: 120°E-80°W, 25°S-25°N.
    BasinExtent::new(120.0, -80.0, -25.0, 25.0)
}

fn grid() -> GridSpec {
    GridSpec::new(NX, NY, extent()).expect("a 3x2 basin is a valid grid")
}

fn header() -> RunHeader {
    // Scenario values, not physical constants: they exist to be round-tripped.
    // JSON carries floats as decimal text, so the fixture leans on the cases
    // where that is least obviously lossless — a negative zero, a subnormal,
    // and the extremes of f64 — rather than on tidy decimals.
    let params = PhysicalParams {
        mean_depth_m: -0.0,
        reduced_gravity_m_per_s2: f64::MIN_POSITIVE / 2.0,
        beta_per_m_per_s: f64::MAX,
        rayleigh_damping_per_s: f64::MIN,
        reference_density_kg_per_m3: 1.234_567_890_123_456_7e-9,
    };
    RunHeader::new(
        grid(),
        params,
        "steady alizes, one westerly wind burst",
        OutputTiming {
            frame_count: 7,
            interval_s: 86_400.0,
        },
    )
}

/// Values chosen to stress the encoder rather than to be physical: a negative
/// zero, a subnormal, and the extremes of f64 all have bit patterns that a
/// lossy or normalizing encoder would not return unchanged.
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

fn frame() -> Frame {
    let g = grid();
    Frame::new(
        1_234.5,
        &g,
        awkward_values(g.field_len(Variable::ThermoclineDepthAnomaly), 0),
        awkward_values(g.field_len(Variable::ZonalCurrentAnomaly), 1),
        awkward_values(g.field_len(Variable::MeridionalCurrentAnomaly), 2),
        awkward_values(g.field_len(Variable::ZonalWindStress), 3),
        awkward_values(g.field_len(Variable::MeridionalWindStress), 4),
    )
    .expect("every field is built at the length the grid asks for")
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

#[test]
fn header_survives_a_json_round_trip_bit_for_bit() {
    let original = header();
    let json = serde_json::to_string(&original).expect("the header is serializable");
    let decoded: RunHeader = serde_json::from_str(&json).expect("the header is deserializable");

    // Bit-for-bit, per the acceptance criteria: re-encoding the decoded header
    // must reproduce the original bytes exactly.
    let reencoded = serde_json::to_string(&decoded).expect("the header is serializable");
    assert_eq!(json, reencoded);
    assert_eq!(original, decoded);

    // `PartialEq` on f64 treats 0.0 and -0.0 as equal, so the params are also
    // compared as bit patterns: JSON's decimal text is where a sign or a low
    // bit would quietly go missing.
    let written = original.physical_params;
    let read = decoded.physical_params;
    assert_bit_identical(
        "physical_params",
        &[
            written.mean_depth_m,
            written.reduced_gravity_m_per_s2,
            written.beta_per_m_per_s,
            written.rayleigh_damping_per_s,
            written.reference_density_kg_per_m3,
        ],
        &[
            read.mean_depth_m,
            read.reduced_gravity_m_per_s2,
            read.beta_per_m_per_s,
            read.rayleigh_damping_per_s,
            read.reference_density_kg_per_m3,
        ],
    );
}

#[test]
fn a_written_header_is_self_describing() {
    // The point of a JSON header (ADR-0004) is that a reader never has to
    // guess: the version, the grid, the units and the output cadence are all
    // on the page.
    let json = serde_json::to_value(header()).expect("the header is serializable");
    assert_eq!(json["format_version"], FORMAT_VERSION);
    assert_eq!(json["grid"]["nx"], NX);
    assert_eq!(json["grid"]["ny"], NY);
    assert_eq!(json["output"]["frame_count"], 7);

    let units: Vec<&str> = json["variables"]
        .as_array()
        .expect("variables is a list")
        .iter()
        .map(|v| v["unit"].as_str().expect("each variable states its unit"))
        .collect();
    // Units of h, u, v, tau_x, tau_y from 01-scientific-model.md.
    assert_eq!(units, ["m", "m s^-1", "m s^-1", "N m^-2", "N m^-2"]);

    let symbols: Vec<&str> = json["variables"]
        .as_array()
        .expect("variables is a list")
        .iter()
        .map(|v| {
            v["symbol"]
                .as_str()
                .expect("each variable states its symbol")
        })
        .collect();
    assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y"]);
}

#[test]
fn a_frame_is_written_under_the_field_names_the_format_specifies() {
    // The frame's fields are named for the model's symbols, and its time is
    // `t` — the unit lives in the Rust name and in the header, not on the
    // wire, so a reader in another language finds what it was promised.
    let json = serde_json::to_value(frame()).expect("the frame is serializable");
    let mut keys: Vec<&str> = json
        .as_object()
        .expect("a frame is a record")
        .keys()
        .map(String::as_str)
        .collect();
    keys.sort_unstable();
    assert_eq!(keys, ["h", "sst", "t", "tau_x", "tau_y", "u", "v"]);
    // `sst` is on the page of every frame, and it is null on a frame of a run
    // that did not couple SST: the field says the run has no `T'`, rather than
    // being quietly dropped so that a reader cannot tell an uncoupled run from
    // a frame that lost a field (T-05.4).
    assert!(json["sst"].is_null());

    let coupled = frame()
        .with_sst_anomaly(&grid(), awkward_values(NX * NY, 5))
        .expect("the anomaly is built at the length the grid asks for");
    let json = serde_json::to_value(coupled).expect("the frame is serializable");
    assert_eq!(
        json["sst"]
            .as_array()
            .expect("a coupled frame's anomaly is a list")
            .len(),
        NX * NY
    );
}

#[test]
fn frame_survives_a_bincode_round_trip_value_for_value() {
    let original = frame();
    let config = bincode::config::standard();
    let bytes = bincode::serde::encode_to_vec(&original, config).expect("the frame encodes");
    let (decoded, consumed): (Frame, usize) =
        bincode::serde::decode_from_slice(&bytes, config).expect("the frame decodes");

    assert_eq!(
        consumed,
        bytes.len(),
        "the frame decoder read the whole blob"
    );
    assert_eq!(original.t_s().to_bits(), decoded.t_s().to_bits());
    assert_bit_identical("h", original.h(), decoded.h());
    assert_bit_identical("u", original.u(), decoded.u());
    assert_bit_identical("v", original.v(), decoded.v());
    assert_bit_identical("tau_x", original.tau_x(), decoded.tau_x());
    assert_bit_identical("tau_y", original.tau_y(), decoded.tau_y());
}

#[test]
fn a_frame_carries_one_value_per_point_of_each_staggered_field() {
    // The C-grid staggering of ADR-0003 gives the face fields one extra line
    // of points: u spans nx+1 by ny, v spans nx by ny+1, h nx by ny.
    let g = grid();
    assert_eq!(g.field_len(Variable::ThermoclineDepthAnomaly), NX * NY);
    assert_eq!(g.field_len(Variable::ZonalCurrentAnomaly), (NX + 1) * NY);
    assert_eq!(
        g.field_len(Variable::MeridionalCurrentAnomaly),
        NX * (NY + 1)
    );
    // Wind stress forces the currents, so each component sits where the
    // current it forces sits.
    assert_eq!(g.field_len(Variable::ZonalWindStress), (NX + 1) * NY);
    assert_eq!(g.field_len(Variable::MeridionalWindStress), NX * (NY + 1));
}

#[test]
fn a_frame_that_does_not_fit_the_grid_is_rejected_by_name() {
    // Invalid input returns a Result naming the offending value and the bound
    // it violated (CODING_STANDARDS.md), rather than writing a truncated run.
    let g = grid();
    let err = Frame::new(
        0.0,
        &g,
        vec![0.0; NX * NY],
        vec![0.0; NX * NY], // u needs (nx + 1) * ny
        vec![0.0; NX * (NY + 1)],
        vec![0.0; (NX + 1) * NY],
        vec![0.0; NX * (NY + 1)],
    )
    .expect_err("a cell-centred u field does not cover the eastern boundary");

    assert_eq!(
        err,
        FormatError::FieldShape {
            variable: Variable::ZonalCurrentAnomaly,
            expected: (NX + 1) * NY,
            actual: NX * NY,
        }
    );
    let message = err.to_string();
    assert!(message.contains('u'), "{message}");
    assert!(message.contains("expected 8"), "{message}");
}

#[test]
fn a_basin_with_no_cells_is_rejected() {
    let err = GridSpec::new(0, NY, extent()).expect_err("a basin needs at least one column");
    let message = err.to_string();
    assert!(message.contains("nx is 0"), "{message}");
}

#[test]
fn a_header_claiming_an_empty_basin_is_rejected_on_the_way_in() {
    // A truncated or hand-edited header is invalid input, not a broken
    // invariant, so reading one back returns an error rather than panicking
    // the first time something asks the grid for a field length.
    let json = serde_json::to_string(&header()).expect("the header is serializable");
    let corrupted = json.replace("\"nx\":3", "\"nx\":0");
    assert_ne!(
        json, corrupted,
        "the fixture must actually name a cell count"
    );

    let err = serde_json::from_str::<RunHeader>(&corrupted)
        .expect_err("a basin with no columns is not a grid");
    let message = err.to_string();
    assert!(message.contains("nx is 0"), "{message}");
}
