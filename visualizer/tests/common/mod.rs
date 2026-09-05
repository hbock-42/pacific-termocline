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
pub fn encoded_frames(header: &RunHeader, count: u64) -> Vec<u8> {
    let grid = header.grid;
    let field = |variable| vec![0.0; grid.field_len(variable)];
    let mut frames = Vec::new();
    for index in 0..count {
        #[allow(clippy::cast_precision_loss)]
        let t_s = index as f64 * header.output.interval_s;
        let frame = Frame::new(
            t_s,
            &grid,
            field(Variable::ThermoclineDepthAnomaly),
            field(Variable::ZonalCurrentAnomaly),
            field(Variable::MeridionalCurrentAnomaly),
            field(Variable::ZonalWindStress),
            field(Variable::MeridionalWindStress),
        )
        .expect("fields sized from the grid fit it");
        frames.extend(
            bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
        );
    }
    frames
}
