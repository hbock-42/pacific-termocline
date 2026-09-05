//! T-09.3 acceptance criterion: the cross-section of a known equilibrium frame
//! matches the analytic tilt of T-07.4.
//!
//! The criterion is written as a visual match, and a visual match nobody runs
//! is not a criterion. So the line is asserted instead of looked at:
//! [`CrossSection`] is the same list of points the shell draws the polyline
//! from, built without a window or a device, and what a reader would check by
//! eye — high in the west, falling monotonically, crossing zero once, low in
//! the east — is checked here by index ([ADR-0006]).
//!
//! # Where the expected values come from
//!
//! The profile is T-07.4's measured equilibrium of
//! `engine/scenarios/steady-trades.toml`: `h` = +38.2 m at the western wall and
//! -28.2 m at the eastern, a 66.4 m west-to-east drop, linear between. It is
//! the *input* these tests draw, named by the acceptance criterion itself
//! ("matches the analytic tilt from T-07.4"), and the expected value at each
//! longitude is that straight line evaluated there — never a number read back
//! off this code. `CONTEXT.md`, *Thermocline tilt*, is what makes the sign
//! right: sustained easterly stress piles the warm layer up in the west, and a
//! positive `h` is a deeper-than-average thermocline.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod common;

use common::{encoded_frames_with_h, header_on, steady_trades_header, NX, NY, PACIFIC};
use termocline_format::{BasinExtent, GridSpec, RunHeader, Variable};
use visualizer::{CrossSection, LoadedRun, RunBytes};

/// T-07.4's equilibrium `h` at the western wall of the steady-trades basin, in
/// metres. Positive is deeper than the mean depth `H` (`CONTEXT.md`).
const WESTERN_WALL_H_M: f64 = 38.2;
/// T-07.4's equilibrium `h` at the eastern wall, in metres.
const EASTERN_WALL_H_M: f64 = -28.2;

/// How far a sampled anomaly may sit from the analytic line, in metres.
///
/// The section reads a value out of the frame and, where the equator falls
/// between two rows, averages two of them. That is a handful of flops on
/// quantities of order 40 m, so the round-off is of order `40 · 2.2e-16` ≈
/// 1e-14 m. A nanometre is six orders of magnitude clear of that and still
/// twelve orders below the metre-scale structure the criterion is about.
const SAMPLE_TOLERANCE_M: f64 = 1e-9;

/// How far a computed plot position may sit from the fraction the geometry
/// gives. Dimensionless, and round-off on quantities of order one.
const POSITION_TOLERANCE: f64 = 1e-12;

/// How far a computed longitude may sit from the one the basin's geometry
/// gives, in degrees. Round-off on degree-scale quantities, as above.
const LONGITUDE_TOLERANCE_DEG: f64 = 1e-9;

/// The half-degree resolution of `steady-trades.toml` (`engine/src/basin.rs`).
const RESOLUTION_DEG: f64 = 0.5;

/// T-07.4's equilibrium `h` at `fraction` of the way from the western wall to
/// the eastern one: the analytic tilt, evaluated independently of this crate.
fn analytic_tilt_m(fraction: f64) -> f64 {
    WESTERN_WALL_H_M + (EASTERN_WALL_H_M - WESTERN_WALL_H_M) * fraction
}

/// The equilibrium tilt of T-07.4 as a field on a basin `nx` by `ny` cells:
/// the analytic line in `x`, uniform in `y`.
///
/// The wall values sit at the walls, so the cell centres carry the line
/// sampled at `(i + 0.5)/nx` of the way across.
fn tilt_field(nx: usize, ny: usize) -> Vec<f64> {
    let mut field = Vec::with_capacity(nx * ny);
    for _ in 0..ny {
        for i in 0..nx {
            #[allow(clippy::cast_precision_loss)]
            let fraction = (i as f64 + 0.5) / nx as f64;
            field.push(analytic_tilt_m(fraction));
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

/// The cross-section of the only frame of a one-frame run carrying `h_m`.
fn section_of(header: &RunHeader, h_m: Vec<f64>) -> CrossSection {
    let run = run_of_one_frame(header, h_m);
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid")
}

#[test]
fn the_equilibrium_cross_section_matches_the_analytic_tilt() {
    // The acceptance criterion itself. The run is steady-trades at its
    // equilibrium, and the line drawn along the equator is compared point by
    // point with the analytic tilt of T-07.4.
    let header = steady_trades_header(1);
    let section = section_of(&header, tilt_field(NX, NY));
    let points = section.points();
    assert_eq!(points.len(), NX, "one point per cell of the basin");

    for (i, point) in points.iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let expected_m = analytic_tilt_m((i as f64 + 0.5) / NX as f64);
        assert!(
            (point.h_m() - expected_m).abs() < SAMPLE_TOLERANCE_M,
            "point {i}: expected {expected_m} m, got {} m",
            point.h_m()
        );
    }

    // What a reader checks by eye, in the order they check it.
    assert!(
        points[0].h_m() > 0.0,
        "the thermocline is deeper than average at the western wall"
    );
    assert!(
        points[NX - 1].h_m() < 0.0,
        "and shallower than average at the eastern one"
    );
    assert!(
        points.windows(2).all(|pair| pair[1].h_m() < pair[0].h_m()),
        "the tilt falls monotonically from west to east"
    );
    let crossings = points
        .windows(2)
        .filter(|pair| pair[0].h_m().signum() != pair[1].h_m().signum())
        .count();
    assert_eq!(crossings, 1, "and changes sign exactly once");
}

#[test]
fn the_section_runs_west_to_east_across_the_antimeridian() {
    let header = steady_trades_header(1);
    let section = section_of(&header, tilt_field(NX, NY));
    let points = section.points();

    // The basin spans 120°E to 80°W eastward across the antimeridian
    // (`CONTEXT.md`, *Basin*), so the first cell centre is half a cell east of
    // 120°E and the last is half a cell west of 80°W. Longitudes are written
    // in the same degrees-east convention `BasinExtent` uses, so the eastern
    // half of the basin is negative.
    let first_deg_east = PACIFIC.west_deg_east + RESOLUTION_DEG / 2.0;
    let last_deg_east = PACIFIC.east_deg_east - RESOLUTION_DEG / 2.0;
    assert!(
        (points[0].longitude_deg_east() - first_deg_east).abs() < LONGITUDE_TOLERANCE_DEG,
        "the first point sits at {first_deg_east}°E, not {}°E",
        points[0].longitude_deg_east()
    );
    assert!(
        (points[NX - 1].longitude_deg_east() - last_deg_east).abs() < LONGITUDE_TOLERANCE_DEG,
        "the last point sits at {last_deg_east}°E, not {}°E",
        points[NX - 1].longitude_deg_east()
    );

    // The axis position is monotonic even though the longitude label wraps
    // through the antimeridian: a plot whose x fell back by 360° in the middle
    // would draw the basin folded over itself.
    assert!(
        points
            .windows(2)
            .all(|pair| pair[1].x_fraction() > pair[0].x_fraction()),
        "the axis runs west to east without folding at the antimeridian"
    );
    assert!(points[0].x_fraction() > 0.0 && points[NX - 1].x_fraction() < 1.0);
}

/// A field whose `h` at every cell is the latitude of that cell's centre, in
/// metres per degree north: nothing varies with `x`, so what the section reads
/// says only which row (or rows) it read.
fn latitude_marked_field(grid: GridSpec) -> Vec<f64> {
    let extent = grid.extent();
    #[allow(clippy::cast_precision_loss)]
    let dy_deg = (extent.north_deg_north - extent.south_deg_north) / grid.ny() as f64;
    let mut field = Vec::with_capacity(grid.field_len(Variable::ThermoclineDepthAnomaly));
    for j in 0..grid.ny() {
        #[allow(clippy::cast_precision_loss)]
        let latitude_deg = (j as f64 + 0.5).mul_add(dy_deg, extent.south_deg_north);
        for _ in 0..grid.nx() {
            field.push(latitude_deg);
        }
    }
    field
}

#[test]
fn the_section_is_taken_at_the_equator_when_it_falls_between_two_rows() {
    // The scenario basin is 25°S–25°N in 100 half-degree cells, so no cell
    // centre is on the equator: the two nearest sit at 0.25°S and 0.25°N. A
    // section that took either one alone would be off the equator, where the
    // waveguide the whole model is about is not centred.
    let header = steady_trades_header(1);
    let grid = header.grid;
    let section = section_of(&header, latitude_marked_field(grid));
    assert_eq!(section.latitude_deg_north(), 0.0);
    assert_eq!(section.rows_averaged(), 2);
    for point in section.points() {
        assert!(
            point.h_m().abs() < SAMPLE_TOLERANCE_M,
            "a section on the equator reads 0, not {}",
            point.h_m()
        );
    }
}

#[test]
fn the_section_is_taken_at_the_equator_when_a_row_sits_on_it() {
    // 101 half-degree cells from 25.25°S to 25.25°N put a cell centre exactly
    // on the equator, and then there is nothing to average.
    let extent = BasinExtent::new(120.0, -80.0, -25.25, 25.25);
    let grid = GridSpec::new(NX, 101, extent).expect("320 x 101 is a valid basin");
    let header = header_on(grid, "equatorial-row", 1);
    let section = section_of(&header, latitude_marked_field(grid));
    assert_eq!(section.latitude_deg_north(), 0.0);
    assert_eq!(section.rows_averaged(), 1);
    for point in section.points() {
        assert!(
            point.h_m().abs() < SAMPLE_TOLERANCE_M,
            "a section on the equator reads 0, not {}",
            point.h_m()
        );
    }
}

#[test]
fn the_vertical_axis_is_the_runs_scale_so_a_collapsing_tilt_is_seen_to_collapse() {
    // Two frames: the equilibrium tilt, then a tenth of it. The second frame
    // must be drawn a tenth as tall, which is only true if the axis is the
    // run's and not each frame's — per-frame autoscaling would renormalize the
    // collapse away, which is the one thing this view exists to show.
    let header = steady_trades_header(2);
    let full = tilt_field(NX, NY);
    let collapsed: Vec<f64> = full.iter().map(|h_m| h_m / 10.0).collect();
    let bytes = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames_with_h(&header, 2, |index| {
            if index == 0 {
                full.clone()
            } else {
                collapsed.clone()
            }
        }),
    };
    let run = LoadedRun::from_bytes("collapsing", bytes).expect("the run loads");
    let scale = run.anomaly_scale();
    let section_of_frame = |index: u64| {
        let frame = run.frame(index).expect("the run holds this frame");
        CrossSection::of_frame(run.header().grid, &frame, scale)
            .expect("the frame fits its own grid")
    };

    // The run's scale reaches as far as the deepest anomaly anywhere in it,
    // which is the first cell centre of the uncollapsed frame.
    #[allow(clippy::cast_precision_loss)]
    let deepest_m = analytic_tilt_m(0.5 / NX as f64);
    assert!((scale.half_range_m() - deepest_m).abs() < SAMPLE_TOLERANCE_M);

    // Plot coordinates run from the top of the panel down, so a deeper-than-
    // average anomaly is above the zero line and a shallower one below it.
    let first = section_of_frame(0);
    let second = section_of_frame(1);
    let (_, y_full) = first
        .plot_position(&first.points()[0])
        .expect("a finite anomaly has a position");
    let (_, y_collapsed) = second
        .plot_position(&second.points()[0])
        .expect("a finite anomaly has a position");
    // The scale is symmetric about zero, so half-height is the zero line: a
    // value's height above it is proportional to the value.
    let above_zero = |y: f64| 0.5 - y;
    assert!(
        (above_zero(y_collapsed) * 10.0 - above_zero(y_full)).abs() < POSITION_TOLERANCE,
        "a tenth of the anomaly is drawn a tenth as far from the zero line"
    );
}

#[test]
fn an_anomaly_that_is_not_a_number_has_no_place_on_the_line() {
    // A `NaN` in `h` means the integration diverged. The line must break
    // there rather than be drawn through a value the run never produced.
    let header = steady_trades_header(1);
    let mut field = tilt_field(NX, NY);
    for row in 0..NY {
        field[row * NX + 7] = f64::NAN;
    }
    let section = section_of(&header, field);
    let points = section.points();
    assert!(section.plot_position(&points[7]).is_none());
    assert!(section.plot_position(&points[6]).is_some());
    assert!(section.plot_position(&points[8]).is_some());
}
