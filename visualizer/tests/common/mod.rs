//! The run fixture both integration tests build on.
//!
//! The visualizer must not link the engine (ADR-0001), so a run under test is
//! written here from `termocline_format` alone. The values come from
//! `engine/scenarios/steady-trades.toml` and `CONTEXT.md`, *Basin* — never
//! from reading back a run this code produced.

// Each integration test compiles this module separately and uses part of it,
// so what one test does not touch is dead code in that test's build.
#![allow(dead_code)]

use termocline_format::{
    frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader, Variable,
};
use visualizer::RunBytes;

/// The basin of `CONTEXT.md`, *Basin*: 120°E–80°W by 25°S–25°N.
pub const PACIFIC: BasinExtent = BasinExtent::new(120.0, -80.0, -25.0, 25.0);

/// 160° of longitude and 50° of latitude at `steady-trades.toml`'s 0.5°
/// resolution.
pub const NX: usize = 320;
pub const NY: usize = 100;

/// `steady-trades.toml` writes a frame every 24 steps of an hour.
pub const FRAME_INTERVAL_S: f64 = 86_400.0;

/// Frames `steady-trades.toml` writes: 17 520 steps of an hour with a frame
/// every 24 makes 730 daily frames, and the frame at t = 0 makes 731.
pub const STEADY_TRADES_FRAMES: u64 = 731;

/// The physical parameters of `steady-trades.toml`: g' = 0.06 m s^-2 and
/// H = 150 m give c = √(g'H) = 3.0 m s^-1, the observed first-baroclinic
/// Kelvin speed of the equatorial Pacific (CONTEXT.md).
pub const STEADY_TRADES_PARAMS: PhysicalParams = PhysicalParams {
    mean_depth_m: 150.0,
    reduced_gravity_m_per_s2: 0.06,
    beta_per_m_per_s: 2.3e-11,
    rayleigh_damping_per_s: 1.0e-7,
    reference_density_kg_per_m3: 1025.0,
};

/// The header `steady-trades.toml` produces for a run of `frame_count` frames.
pub fn steady_trades_header(frame_count: u64) -> RunHeader {
    header_on(
        GridSpec::new(NX, NY, PACIFIC).expect("320 x 100 is a valid basin"),
        "steady-trades",
        frame_count,
    )
}

/// The header of a run on `grid`, for a test that wants a smaller basin than
/// the scenario's.
pub fn header_on(grid: GridSpec, scenario: &str, frame_count: u64) -> RunHeader {
    RunHeader::new(
        grid,
        STEADY_TRADES_PARAMS,
        scenario,
        OutputTiming {
            frame_count,
            interval_s: FRAME_INTERVAL_S,
        },
    )
}

/// The two byte sources of the run `header` describes, as they would arrive
/// from a file drop or an HTTP fetch.
///
/// The fields are zero-filled: what is under test is the header and the shape
/// of the frames, never their values.
pub fn run_bytes(header: &RunHeader) -> RunBytes {
    RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames(header, header.output.frame_count),
    }
}

/// `count` encoded frames on `header`'s grid, spaced by its output interval.
///
/// Every field is zero: what these serve is a test of the header and of the
/// shape of the frames, never of their values.
pub fn encoded_frames(header: &RunHeader, count: u64) -> Vec<u8> {
    let zero_h = vec![0.0; header.grid.field_len(Variable::ThermoclineDepthAnomaly)];
    encoded_frames_with_h(header, count, |_| zero_h.clone())
}

/// `count` encoded frames on `header`'s grid whose thermocline depth anomaly
/// `h` is `h_m(index)`, in metres, and whose every other field is zero.
///
/// `h` alone because it is the field the basin map draws; the currents and the
/// wind stress only have to be the right shape for the frame to encode.
pub fn encoded_frames_with_h(
    header: &RunHeader,
    count: u64,
    h_m: impl Fn(u64) -> Vec<f64>,
) -> Vec<u8> {
    encoded_frames_with_fields(header, count, |index| FrameFields {
        h_m: h_m(index),
        ..FrameFields::calm(header)
    })
}

/// The fields of one frame that a test cares about: the anomaly the map draws
/// and the stress the overlay draws. The currents are not among them — nothing
/// in the visualizer reads `u` or `v` yet — so they are zero-filled below.
pub struct FrameFields {
    /// Thermocline depth anomaly `h`, in metres, at cell centres.
    pub h_m: Vec<f64>,
    /// Zonal wind stress `τx`, in pascals, at east/west faces.
    pub tau_x_pa: Vec<f64>,
    /// Meridional wind stress `τy`, in pascals, at north/south faces.
    pub tau_y_pa: Vec<f64>,
}

impl FrameFields {
    /// An ocean at rest under no wind, on `header`'s grid.
    pub fn calm(header: &RunHeader) -> Self {
        let field = |variable| vec![0.0; header.grid.field_len(variable)];
        Self {
            h_m: field(Variable::ThermoclineDepthAnomaly),
            tau_x_pa: field(Variable::ZonalWindStress),
            tau_y_pa: field(Variable::MeridionalWindStress),
        }
    }
}

/// `count` encoded frames on `header`'s grid, spaced by its output interval,
/// carrying the fields `fields(index)` gives for each.
pub fn encoded_frames_with_fields(
    header: &RunHeader,
    count: u64,
    fields: impl Fn(u64) -> FrameFields,
) -> Vec<u8> {
    let grid = header.grid;
    let zero = |variable| vec![0.0; grid.field_len(variable)];
    let mut frames = Vec::new();
    for index in 0..count {
        #[allow(clippy::cast_precision_loss)]
        let t_s = index as f64 * header.output.interval_s;
        let FrameFields {
            h_m,
            tau_x_pa,
            tau_y_pa,
        } = fields(index);
        let frame = Frame::new(
            t_s,
            &grid,
            h_m,
            zero(Variable::ZonalCurrentAnomaly),
            zero(Variable::MeridionalCurrentAnomaly),
            tau_x_pa,
            tau_y_pa,
        )
        .expect("fields sized from the grid fit it");
        frames.extend(
            bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
        );
    }
    frames
}
