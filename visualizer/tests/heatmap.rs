//! T-08.2 acceptance criterion: the rendered heatmap of a known run's
//! equilibrium frame shows a deeper thermocline in the west than in the east.
//!
//! The criterion is written as a visual smoke test, and a smoke test nobody
//! runs is not a criterion. So the pixels are asserted instead of looked at:
//! [`Heatmap`] is the same colour-mapped image the shell uploads to the GPU,
//! built without a window or a device, and what a reader would check by eye —
//! warm in the west, cool in the east, neutral where the anomaly changes sign
//! — is checked here by index.
//!
//! # Where the expected values come from
//!
//! The tilt is T-07.4's measured equilibrium of `engine/scenarios/steady-trades.toml`:
//! `h` = +38.2 m at the western wall and -28.2 m at the eastern, a 60.1 m
//! west-to-east drop. It is the *input* these tests draw, named by the
//! acceptance criterion itself ("matching T-07.4's known result"), not an
//! expected output read back off this code; what is asserted about it — the
//! colours — comes from ColorBrewer. `CONTEXT.md`, *Thermocline tilt*, is what
//! makes the sign right: sustained easterly stress piles the warm layer up in
//! the west, and a positive `h` is a deeper-than-average thermocline.
//!
//! The colours are ColorBrewer's 11-class `RdBu` diverging scheme
//! (Brewer, Harrower & Pennsylvania State University, <https://colorbrewer2.org>),
//! transcribed below from the published table and never read back off this
//! code.

mod common;

use common::{encoded_frames_with_h, steady_trades_header, NX, NY};
use termocline_format::RunHeader;
use visualizer::{DivergingScale, Heatmap, LoadedRun, RunBytes};

/// T-07.4's equilibrium `h` at the western wall of the steady-trades basin, in
/// metres. Positive is deeper than the mean depth `H` (`CONTEXT.md`).
const WESTERN_WALL_H_M: f64 = 38.2;
/// T-07.4's equilibrium `h` at the eastern wall, in metres.
const EASTERN_WALL_H_M: f64 = -28.2;

/// ColorBrewer 11-class `RdBu`, reversed so it runs from the most negative
/// anomaly to the most positive: blue (shallow) through near-white (no
/// anomaly) to red (deep).
const RD_BU_11: [[u8; 3]; 11] = [
    [5, 48, 97],
    [33, 102, 172],
    [67, 147, 195],
    [146, 197, 222],
    [209, 229, 240],
    [247, 247, 247],
    [253, 219, 199],
    [244, 165, 130],
    [214, 96, 77],
    [178, 24, 43],
    [103, 0, 31],
];

/// The neutral colour a zero anomaly is drawn in: the middle class of `RdBu`.
const NEUTRAL: [u8; 3] = RD_BU_11[5];
/// The colour of the deepest anomaly on the scale.
const DEEPEST: [u8; 3] = RD_BU_11[10];
/// The colour of the shallowest anomaly on the scale.
const SHALLOWEST: [u8; 3] = RD_BU_11[0];

/// A rounding tolerance of one 8-bit level per channel.
///
/// The ramp is interpolated in `f64` and quantized to `u8` once, so a colour
/// asked for at a class boundary can land a level either side of the published
/// value. Nothing about the test's meaning survives a wider bound: two
/// adjacent `RdBu` classes differ by tens of levels.
const CHANNEL_TOLERANCE: i32 = 1;

/// The fraction of the scale's half-range below which an anomaly is not
/// resolvable in eight bits.
///
/// The middle class of `RdBu` moves 48 levels of blue over a fifth of the
/// half-range, so one level is worth about `half_range / 240` metres and
/// anything smaller quantizes to the neutral colour by construction. One part
/// in 128 is the next round number clear of that. It is a bound on what the
/// display can show, not a tolerance on the physics.
const RESOLVABLE_FRACTION: f64 = 1.0 / 128.0;

/// Assert `actual` is `expected` to within [`CHANNEL_TOLERANCE`].
#[track_caller]
fn assert_close(actual: [u8; 3], expected: [u8; 3], what: &str) {
    for channel in 0..3 {
        let difference = i32::from(actual[channel]) - i32::from(expected[channel]);
        assert!(
            difference.abs() <= CHANNEL_TOLERANCE,
            "{what}: expected {expected:?}, got {actual:?}"
        );
    }
}

/// Whether a colour is on the warm (deep) half of the scale.
fn is_warm(rgb: [u8; 3]) -> bool {
    rgb[0] > rgb[2]
}

/// Whether a colour is on the cool (shallow) half of the scale.
fn is_cool(rgb: [u8; 3]) -> bool {
    rgb[2] > rgb[0]
}

/// The equilibrium tilt of T-07.4 as a field on a basin `nx` cells wide:
/// linear in `x` between the two wall values, uniform in `y`.
///
/// A linear profile is not what the solver produces in detail — the tilt has
/// structure near the equator — but the criterion is about its sign and its
/// direction, and those the straight line carries exactly.
fn tilt_field(nx: usize, ny: usize, west_m: f64, east_m: f64) -> Vec<f64> {
    let mut field = Vec::with_capacity(nx * ny);
    for _ in 0..ny {
        for i in 0..nx {
            #[allow(clippy::cast_precision_loss)]
            let fraction = i as f64 / (nx - 1) as f64;
            field.push(west_m + (east_m - west_m) * fraction);
        }
    }
    field
}

/// A one-frame run of `header`'s shape whose only frame carries `h_m`.
fn run_of_one_frame(header: &RunHeader, h_m: Vec<f64>) -> LoadedRun {
    let bytes = RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames_with_h(header, header.output.frame_count, |_| h_m.clone()),
    };
    LoadedRun::from_bytes("run-steady-trades", bytes).expect("the run loads")
}

/// The heatmap of the only frame of a one-frame run carrying `h_m`.
fn heatmap_of(header: &RunHeader, h_m: Vec<f64>) -> Heatmap {
    let run = run_of_one_frame(header, h_m);
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid")
}

#[test]
fn the_equilibrium_tilt_is_drawn_deep_in_the_west_and_shallow_in_the_east() {
    // The acceptance criterion itself. The run is steady-trades at its
    // measured equilibrium; the map must show the two ends of the basin in
    // opposite halves of the scale, with the deep end in the west.
    let header = steady_trades_header(1);
    let heatmap = heatmap_of(
        &header,
        tilt_field(NX, NY, WESTERN_WALL_H_M, EASTERN_WALL_H_M),
    );

    let row = heatmap.height() / 2;
    let west = heatmap.pixel(0, row).expect("the western column is drawn");
    let east = heatmap
        .pixel(heatmap.width() - 1, row)
        .expect("the eastern column is drawn");

    assert!(
        is_warm(west),
        "the western wall is the deepest anomaly in the basin and must be drawn warm: {west:?}"
    );
    assert!(
        is_cool(east),
        "the eastern wall is a shallow anomaly and must be drawn cool: {east:?}"
    );
    // +38.2 m is the largest magnitude in the field, so it saturates the
    // scale: the western wall is the warm end of `RdBu`, not merely warm.
    assert_close(west, DEEPEST, "the western wall");
}

#[test]
fn every_column_is_drawn_on_the_side_of_the_scale_its_anomaly_is_on() {
    // The claim a reader makes by eye is not just "the ends differ" but "the
    // colour changes sign where the anomaly does". Checking every column also
    // rules out a map that happens to be right at the two walls and scrambled
    // between them.
    let header = steady_trades_header(1);
    let field = tilt_field(NX, NY, WESTERN_WALL_H_M, EASTERN_WALL_H_M);
    let heatmap = heatmap_of(&header, field.clone());

    let row = heatmap.height() / 2;
    // Row `row` of the image is row `height - 1 - row` of the field: the map
    // is drawn north up and the field is stored south first.
    let field_row = &field[(heatmap.height() - 1 - row) * heatmap.width()..];
    let resolvable_m = heatmap.scale().half_range_m() * RESOLVABLE_FRACTION;
    for (x, &anomaly_m) in field_row.iter().take(heatmap.width()).enumerate() {
        // The handful of columns either side of the zero contour are neutral
        // because they are neutral, not because the map is wrong.
        if anomaly_m.abs() < resolvable_m {
            continue;
        }
        let rgb = heatmap.pixel(x, row).expect("every column is drawn");
        if anomaly_m > 0.0 {
            assert!(
                is_warm(rgb),
                "h = {anomaly_m} m at x = {x} is drawn {rgb:?}"
            );
        } else {
            assert!(
                is_cool(rgb),
                "h = {anomaly_m} m at x = {x} is drawn {rgb:?}"
            );
        }
    }
}

#[test]
fn a_tilt_the_other_way_round_is_drawn_the_other_way_round() {
    // A map that ignored x, or flipped it, would pass the test above on a
    // symmetric field. El Niño collapses the tilt and can reverse it, so the
    // reversed case is a state the visualizer will really be asked to draw.
    let header = steady_trades_header(1);
    let heatmap = heatmap_of(
        &header,
        tilt_field(NX, NY, EASTERN_WALL_H_M, WESTERN_WALL_H_M),
    );

    let row = heatmap.height() / 2;
    let west = heatmap.pixel(0, row).expect("the western column is drawn");
    let east = heatmap
        .pixel(heatmap.width() - 1, row)
        .expect("the eastern column is drawn");
    assert!(
        is_cool(west),
        "a reversed tilt is shallow in the west: {west:?}"
    );
    assert!(
        is_warm(east),
        "a reversed tilt is deep in the east: {east:?}"
    );
}

#[test]
fn the_map_is_drawn_north_up() {
    // The field is row-major with `j` increasing northward, and an image's
    // first row is its top one. A map drawn straight from the buffer would put
    // the southern hemisphere on top, which reads as a plausible basin.
    let header = steady_trades_header(1);
    let mut field = vec![0.0; NX * NY];
    // Northernmost row deep, southernmost shallow.
    for i in 0..NX {
        field[(NY - 1) * NX + i] = WESTERN_WALL_H_M;
        field[i] = EASTERN_WALL_H_M;
    }
    let heatmap = heatmap_of(&header, field);

    let top = heatmap.pixel(NX / 2, 0).expect("the top row is drawn");
    let bottom = heatmap
        .pixel(NX / 2, heatmap.height() - 1)
        .expect("the bottom row is drawn");
    assert!(
        is_warm(top),
        "the northern edge belongs at the top: {top:?}"
    );
    assert!(
        is_cool(bottom),
        "the southern edge belongs at the bottom: {bottom:?}"
    );
}

#[test]
fn the_map_carries_one_pixel_per_cell_of_the_basin() {
    // `h` sits at cell centres (ADR-0003), so its field is exactly nx by ny —
    // no boundary row or column, unlike `u` and `v`.
    let heatmap = heatmap_of(&steady_trades_header(1), vec![0.0; NX * NY]);
    assert_eq!(heatmap.width(), NX);
    assert_eq!(heatmap.height(), NY);
}

#[test]
fn the_scale_runs_from_the_shallowest_class_of_rd_bu_to_the_deepest() {
    // The eleven published classes, in order, at the eleven values that land
    // on them. `h` is signed, so the scale is diverging and symmetric about
    // zero: the same magnitude of anomaly is the same distance from neutral
    // whichever way it goes.
    let scale = DivergingScale::symmetric_over(&[-WESTERN_WALL_H_M, WESTERN_WALL_H_M]);
    assert_eq!(scale.half_range_m(), WESTERN_WALL_H_M);
    for (class, expected) in RD_BU_11.iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let fraction = class as f64 / (RD_BU_11.len() - 1) as f64;
        let value_m = WESTERN_WALL_H_M * (2.0 * fraction - 1.0);
        assert_close(scale.color(value_m), *expected, &format!("h = {value_m} m"));
    }
}

#[test]
fn no_anomaly_is_drawn_neutral() {
    // Zero must land exactly on the middle class, not near it: the neutral
    // colour is what makes the sign of the anomaly readable at a glance.
    let scale = DivergingScale::symmetric_over(&[-WESTERN_WALL_H_M, WESTERN_WALL_H_M]);
    assert_eq!(scale.color(0.0), NEUTRAL);
}

#[test]
fn a_run_at_rest_is_drawn_entirely_neutral() {
    // A field of zeros has no range to normalize by. Dividing by it would tint
    // the whole basin whichever way the arithmetic happened to fall.
    let heatmap = heatmap_of(&steady_trades_header(1), vec![0.0; NX * NY]);
    assert_eq!(heatmap.scale().half_range_m(), 0.0);
    for x in 0..heatmap.width() {
        assert_eq!(heatmap.pixel(x, 0), Some(NEUTRAL));
    }
}

#[test]
fn anomalies_past_the_ends_of_the_scale_are_drawn_at_the_ends() {
    let scale = DivergingScale::symmetric_over(&[-WESTERN_WALL_H_M, WESTERN_WALL_H_M]);
    assert_eq!(scale.color(WESTERN_WALL_H_M * 10.0), DEEPEST);
    assert_eq!(scale.color(-WESTERN_WALL_H_M * 10.0), SHALLOWEST);
}

#[test]
fn a_field_that_has_blown_up_is_drawn_as_missing_rather_than_as_an_anomaly() {
    // A NaN in `h` means the integration diverged. Drawn neutral it would read
    // as an undisturbed patch of ocean, which is the one reading that must not
    // be available; black is off the scale entirely.
    let mut field = vec![WESTERN_WALL_H_M; NX * NY];
    field[0] = f64::NAN;
    let heatmap = heatmap_of(&steady_trades_header(1), field);

    // The scale ignores the non-finite value rather than being destroyed by it.
    assert_eq!(heatmap.scale().half_range_m(), WESTERN_WALL_H_M);
    // Row 0 of the field is the southernmost, so it is the map's last row.
    let blown_up = heatmap
        .pixel(0, heatmap.height() - 1)
        .expect("the corner is drawn");
    assert_eq!(blown_up, [0, 0, 0]);
}

#[test]
fn the_frame_drawn_is_the_frame_asked_for() {
    // The heatmap is of one chosen frame index, so picking the wrong frame is
    // the failure that looks most like success.
    let header = steady_trades_header(3);
    let bytes = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames_with_h(&header, 3, |index| {
            #[allow(clippy::cast_precision_loss)]
            let level = index as f64;
            vec![level; NX * NY]
        }),
    };
    let run = LoadedRun::from_bytes("run", bytes).expect("the run loads");

    for index in 0..3 {
        let frame = run.frame(index).expect("the run has three frames");
        #[allow(clippy::cast_precision_loss)]
        let expected = index as f64;
        assert_eq!(frame.h()[0], expected);
        assert_eq!(frame.t_s(), expected * common::FRAME_INTERVAL_S);
    }
    assert!(run.frame(3).is_none(), "a fourth frame was never written");
}

#[test]
fn one_scale_covers_the_whole_run_so_a_decaying_anomaly_is_seen_to_decay() {
    // Rescaling per frame would renormalize every frame back to full
    // saturation, and a tilt collapsing into an El Niño — the thing this
    // visualizer exists to show — would look like a tilt that never moved.
    let header = steady_trades_header(2);
    let strong = vec![WESTERN_WALL_H_M; NX * NY];
    let weak = vec![WESTERN_WALL_H_M / 4.0; NX * NY];
    let bytes = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames_with_h(&header, 2, |index| {
            if index == 0 {
                strong.clone()
            } else {
                weak.clone()
            }
        }),
    };
    let run = LoadedRun::from_bytes("run", bytes).expect("the run loads");
    assert_eq!(run.anomaly_scale().half_range_m(), WESTERN_WALL_H_M);

    let map_of = |index| {
        let frame = run.frame(index).expect("the run has two frames");
        Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
            .expect("the frame fits its own grid")
    };
    assert_close(
        map_of(0).pixel(0, 0).expect("the corner is drawn"),
        DEEPEST,
        "the strong frame saturates the scale",
    );
    // A quarter of the way up the warm half of an eleven-class ramp is class
    // 6.25 — well short of the warm end, which is the whole point.
    let quarter = map_of(1).pixel(0, 0).expect("the corner is drawn");
    assert!(
        is_warm(quarter) && quarter != DEEPEST,
        "a quarter-strength anomaly must read as weaker, not as saturated: {quarter:?}"
    );
}

#[test]
fn the_colour_bar_is_the_scale_it_labels() {
    // The bar is what tells a reader which colour means what. Drawn from
    // anything but the same scale it would be a legend for another map.
    let scale = DivergingScale::symmetric_over(&[-WESTERN_WALL_H_M, WESTERN_WALL_H_M]);
    let samples = RD_BU_11.len();
    let bar = scale.bar_rgb(samples);
    assert_eq!(bar.len(), samples * 3);
    for (class, expected) in RD_BU_11.iter().enumerate() {
        let sample = [bar[class * 3], bar[class * 3 + 1], bar[class * 3 + 2]];
        assert_close(sample, *expected, &format!("bar sample {class}"));
    }
}
