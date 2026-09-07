//! T-13.1 acceptance criteria: the equatorial side view draws the ocean, not a
//! line plot, and draws it at the thermocline's *total* depth.
//!
//! The view is the teaching answer to the same question the scientific views
//! answer: the basin map plots `h`, a depth anomaly, on a top-down map, so
//! reading it needs a reader who already knows the thermocline is an interface
//! and that red means deeper. This one draws the interface where it actually
//! is — `H + h` metres below the sea surface (`CONTEXT.md`, *Thermocline depth
//! anomaly*) — with warm water above it and cold water below.
//!
//! Like [`visualizer::CrossSection`], [`visualizer::Heatmap`] and
//! [`visualizer::WindOverlay`], none of it knows what a GPU is: a
//! [`SideView`] is a list of columns and the two coloured bands each one is
//! drawn as, on a unit rectangle. So the criteria — the depth of the
//! interface, its tilt, which layer is on top, whether the two layers can be
//! told apart — are asserted here rather than looked at ([ADR-0006]).
//!
//! # Where the expected values come from
//!
//! The equilibrium profile of `engine/scenarios/steady-trades.toml`, as
//! T-07.4 measured it and as T-13.1 restates it: `h` = +38.2 m at the western
//! wall, +8.4 m mid-basin and −28.2 m at the eastern wall (`tests/common`).
//! With that scenario's `H` = 150 m the interface therefore sits 188.2 m down
//! in the west, 158.4 m down mid-basin and 121.8 m down in the east — the
//! numbers this file checks, and none of them read back out of this crate.
//!
//! The profile *between* those three stations is a piecewise-linear stand-in,
//! and is not claimed to be T-07.4's: the measured profile is smooth and
//! curved. It carries exactly the properties these tests are about — the three
//! measured stations, a monotonic fall from west to east, one change of sign —
//! and `tests/cross_section.rs` holds the same data's curvature to the closed
//! form already.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod common;

use common::{
    encoded_frames_with_h, steady_trades_header, EASTERN_WALL_H_M, MID_BASIN_H_M, NX, NY,
    STEADY_TRADES_PARAMS, WESTERN_WALL_H_M,
};
use termocline_format::RunHeader;
use visualizer::{
    CrossSection, LoadedRun, RunBytes, SideView, SideViewColumn, COLD_LAYER_RGB, WARM_LAYER_RGB,
};

/// How far a computed position on the panel may sit from the fraction the
/// geometry gives. Dimensionless, and round-off on quantities of order one.
const POSITION_TOLERANCE: f64 = 1e-12;

/// How far the interface at the outermost cell centre may sit from the value
/// measured *at the wall*, in metres.
///
/// The measured stations are the walls and the middle of the basin; a cell
/// centre is half a cell in from a wall, which is `0.5/320` of the basin's
/// width. The steepest stretch of the measured profile is the eastern half, at
/// `(−28.2 − 8.4)/0.5` = 73.2 m per unit fraction, so half a cell carries the
/// interface `73.2 · 0.5/320` = 0.11 m away from the wall's own value. This
/// bound is that offset rounded up, and it is two orders of magnitude below
/// the 66 m of tilt the criterion is about.
const HALF_CELL_TOLERANCE_M: f64 = 0.15;

/// The smallest contrast ratio, by WCAG 2.1's definition, that counts as two
/// layers a reader can tell apart without a colour bar.
///
/// 3:1 is WCAG 2.1's own threshold for non-text contrast (success criterion
/// 1.4.11, <https://www.w3.org/TR/WCAG21/#non-text-contrast>): the bound
/// published for exactly this question — whether two adjacent blocks of colour
/// read as two things — rather than one chosen here.
const MIN_CONTRAST_RATIO: f64 = 3.0;

/// The mean thermocline depth `H` of `steady-trades.toml`, in metres.
const MEAN_DEPTH_M: f64 = STEADY_TRADES_PARAMS.mean_depth_m;

/// T-07.4's measured equilibrium `h` at `fraction` of the way east across the
/// basin, in metres: linear between the three measured stations.
fn measured_tilt_m(fraction: f64) -> f64 {
    if fraction <= 0.5 {
        WESTERN_WALL_H_M + (MID_BASIN_H_M - WESTERN_WALL_H_M) * (fraction / 0.5)
    } else {
        MID_BASIN_H_M + (EASTERN_WALL_H_M - MID_BASIN_H_M) * ((fraction - 0.5) / 0.5)
    }
}

/// The equilibrium tilt as a field on a basin `nx` by `ny` cells: the measured
/// profile in `x`, uniform in `y`.
fn tilt_field(nx: usize, ny: usize) -> Vec<f64> {
    let mut field = Vec::with_capacity(nx * ny);
    for _ in 0..ny {
        for i in 0..nx {
            #[allow(clippy::cast_precision_loss)]
            let fraction = (i as f64 + 0.5) / nx as f64;
            field.push(measured_tilt_m(fraction));
        }
    }
    field
}

/// A run of `header`'s shape whose frames carry the fields `h_m(index)` gives.
fn run_of(header: &RunHeader, h_m: impl Fn(u64) -> Vec<f64>) -> LoadedRun {
    let bytes = RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames_with_h(header, header.output.frame_count, h_m),
    };
    LoadedRun::from_bytes("run-steady-trades", bytes).expect("the run loads")
}

/// The side view of frame `index` of `run`.
fn side_view_of(run: &LoadedRun, index: u64) -> SideView {
    let frame = run.frame(index).expect("the run holds this frame");
    let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    SideView::of_section(&section, run.header().physical_params.mean_depth_m)
}

/// The side view of the only frame of a one-frame run carrying `h_m`.
fn side_view_of_field(h_m: Vec<f64>) -> SideView {
    let header = steady_trades_header(1);
    let run = run_of(&header, |_| h_m.clone());
    side_view_of(&run, 0)
}

/// The column nearest `fraction` of the way east across the basin.
///
/// The columns tile the basin west to east, one per cell, so the column at a
/// fraction is that fraction of the way along them — the same arithmetic the
/// shell uses to place one.
fn column_at(view: &SideView, fraction: f64) -> &SideViewColumn {
    let columns = view.columns();
    #[allow(clippy::cast_precision_loss, clippy::cast_sign_loss)]
    #[allow(clippy::cast_possible_truncation)]
    let index = ((fraction * columns.len() as f64) as usize).min(columns.len() - 1);
    &columns[index]
}

#[test]
fn the_interface_is_drawn_at_total_depth_not_at_an_anomaly() {
    // The first acceptance criterion, and the reason the ticket exists: `h` is
    // an anomaly about zero (`CONTEXT.md`), and an ocean drawn from it alone
    // would put the thermocline at the surface half the time. What is drawn is
    // `H + h`, so the interface sits where the model says it is.
    let view = side_view_of_field(tilt_field(NX, NY));

    // The three measured stations, checked against T-07.4's own numbers rather
    // than against the profile this test generated from them: `H` plus the
    // measured anomaly, at the cell nearest each station. The cell centre is
    // half a cell in from a wall, which is what [`HALF_CELL_TOLERANCE_M`]
    // allows for.
    for (fraction, measured_h_m) in [
        (0.0, WESTERN_WALL_H_M),
        (0.5, MID_BASIN_H_M),
        (1.0, EASTERN_WALL_H_M),
    ] {
        let column = column_at(&view, fraction);
        let expected_m = MEAN_DEPTH_M + measured_h_m;
        assert!(
            (column.interface_depth_m() - expected_m).abs() < HALF_CELL_TOLERANCE_M,
            "at {fraction} of the way east the interface is {} m down, not {expected_m} m",
            column.interface_depth_m()
        );
        // And it is emphatically not the anomaly: the whole point of the view.
        // The two are `H` apart, and `H` is 150 m against anomalies of tens.
        assert!(
            (column.interface_depth_m() - measured_h_m).abs() > MEAN_DEPTH_M / 2.0,
            "the interface is drawn at a total depth, not at the anomaly {measured_h_m} m"
        );
    }

    // The ticket's own arithmetic, stated as it states it: about 188 m in the
    // west and about 122 m in the east. Those are `H + h` rounded to the metre
    // from 188.2 and 121.8, so the bound is that rounding — half a metre — plus
    // the half-cell offset above.
    let rounding_m = 0.5 + HALF_CELL_TOLERANCE_M;
    let west_m = column_at(&view, 0.0).interface_depth_m();
    let east_m = column_at(&view, 1.0).interface_depth_m();
    assert!(
        (west_m - 188.0).abs() < rounding_m,
        "the western interface is about 188 m down, not {west_m} m"
    );
    assert!(
        (east_m - 122.0).abs() < rounding_m,
        "the eastern interface is about 122 m down, not {east_m} m"
    );
}

#[test]
fn the_interface_slopes_down_towards_the_west() {
    // The third acceptance criterion: run against the control scenario at
    // equilibrium the view shows a clearly deeper interface in the west and a
    // shallower one in the east — the ~60 m tilt of `CONTEXT.md`,
    // *Thermocline tilt*.
    let view = side_view_of_field(tilt_field(NX, NY));
    let columns = view.columns();
    assert_eq!(columns.len(), NX, "one column per cell of the basin");

    assert!(
        columns
            .windows(2)
            .all(|pair| pair[1].interface_depth_m() < pair[0].interface_depth_m()),
        "the interface shoals monotonically from west to east"
    );
    // The measured tilt is 66.4 m between the walls; the outermost cell
    // centres are half a cell in from them, so the drawn drop is a little
    // less. Anything above 60 m is the tilt the ticket names and nothing else
    // in the run is that large.
    let drop_m = columns[0].interface_depth_m() - columns[NX - 1].interface_depth_m();
    assert!(
        drop_m > 60.0,
        "the west-to-east drop is the ~60 m tilt, not {drop_m} m"
    );

    // Every interface is below the sea surface and above the bottom of the
    // panel: an ocean, not a chart that ran off its axis.
    assert!(columns.iter().all(|column| column.interface_depth_m() > 0.0
        && column.interface_depth_m() < view.deepest_drawn_depth_m()));
}

#[test]
fn depth_increases_downward_from_the_sea_surface() {
    // The second half of the second criterion. The panel's coordinates are the
    // convention every other view here uses — `y` down from the north-west
    // corner of a unit rectangle — so the sea surface is at 0 and a deeper
    // interface is drawn *lower* on the panel, which is the one thing a
    // side-on picture of the ocean has to get right.
    let view = side_view_of_field(tilt_field(NX, NY));
    let west = column_at(&view, 0.0);
    let east = column_at(&view, 1.0);
    let fraction_of =
        |column: &SideViewColumn| column.interface_fraction().expect("a finite interface");

    assert!(
        fraction_of(west) > fraction_of(east),
        "the deeper western interface is drawn further down the panel"
    );
    // And the fraction is the depth, measured from the surface: the panel's
    // top edge is 0 m and its bottom edge is the deepest depth drawn.
    for column in view.columns() {
        let expected = column.interface_depth_m() / view.deepest_drawn_depth_m();
        assert!(
            (fraction_of(column) - expected).abs() < POSITION_TOLERANCE,
            "the interface is drawn at its own share of the panel's depth"
        );
    }
}

#[test]
fn warm_water_is_drawn_above_the_interface_and_cold_water_below_it() {
    // The rest of the second criterion: the picture is two layers, and which
    // is which is not a matter of interpretation. The warm layer runs from the
    // sea surface down to the interface and the cold water fills everything
    // beneath it, at every longitude — the upper layer of `CONTEXT.md` sitting
    // on the abyss.
    let view = side_view_of_field(tilt_field(NX, NY));
    for column in view.columns() {
        let [warm, cold] = column.bands().expect("a finite interface is drawn");
        let interface = column
            .interface_fraction()
            .expect("a finite interface has a place");

        assert_eq!(warm.rgb, WARM_LAYER_RGB);
        assert_eq!(cold.rgb, COLD_LAYER_RGB);
        assert!(
            (warm.top_fraction - 0.0).abs() < POSITION_TOLERANCE,
            "the warm layer starts at the sea surface"
        );
        assert!((warm.bottom_fraction - interface).abs() < POSITION_TOLERANCE);
        assert!((cold.top_fraction - interface).abs() < POSITION_TOLERANCE);
        assert!(
            (cold.bottom_fraction - 1.0).abs() < POSITION_TOLERANCE,
            "and the cold water fills the panel beneath it"
        );
        assert!(
            warm.bottom_fraction > warm.top_fraction && cold.bottom_fraction > cold.top_fraction,
            "both layers have depth, so neither is a line"
        );
    }
}

/// The relative luminance of an sRGB colour, by the definition WCAG 2.1 gives
/// (<https://www.w3.org/TR/WCAG21/#dfn-relative-luminance>). Transcribed from
/// the specification rather than taken from this crate.
fn relative_luminance(rgb: [u8; 3]) -> f64 {
    let channel = |value: u8| {
        let srgb = f64::from(value) / 255.0;
        if srgb <= 0.040_45 {
            srgb / 12.92
        } else {
            ((srgb + 0.055) / 1.055).powf(2.4)
        }
    };
    0.0722_f64.mul_add(
        channel(rgb[2]),
        0.2126_f64.mul_add(channel(rgb[0]), 0.7152 * channel(rgb[1])),
    )
}

#[test]
fn the_two_layers_can_be_told_apart_without_a_colour_bar() {
    // The second criterion's own words. Two blocks of colour are told apart by
    // their contrast, and contrast that survives every dichromacy is contrast
    // in *luminance* rather than in hue — so the bound is WCAG's contrast
    // ratio, which is a function of luminance alone and therefore says nothing
    // about whether the reader can see red.
    let (warm, cold) = (
        relative_luminance(WARM_LAYER_RGB),
        relative_luminance(COLD_LAYER_RGB),
    );
    let ratio = (warm.max(cold) + 0.05) / (warm.min(cold) + 0.05);
    assert!(
        ratio >= MIN_CONTRAST_RATIO,
        "the warm and cold layers contrast at {ratio}:1, below the {MIN_CONTRAST_RATIO}:1 that \
         makes two blocks of colour read as two things"
    );
    // The warm layer is the lighter of the two: sunlit water above, dark
    // abyss beneath, which is the way round a reader already expects.
    assert!(
        warm > cold,
        "the warm upper layer is the lighter of the two"
    );
}

#[test]
fn the_depth_axis_is_the_runs_so_the_ocean_does_not_breathe_as_the_run_is_scrubbed() {
    // The axis is the run's, for the reason `crate::heatmap` gives for the
    // colour scale being the run's: a panel rescaled frame by frame would draw
    // every frame's interface at the same place, and a tilt that collapses over
    // a run — which is El Niño (`CONTEXT.md`, *ENSO*) — would look like a tilt
    // that never moved.
    let header = steady_trades_header(2);
    let full = tilt_field(NX, NY);
    let collapsed: Vec<f64> = full.iter().map(|h_m| h_m / 10.0).collect();
    let run = run_of(&header, |index| {
        if index == 0 {
            full.clone()
        } else {
            collapsed.clone()
        }
    });

    let first = side_view_of(&run, 0);
    let second = side_view_of(&run, 1);
    assert_eq!(
        first.deepest_drawn_depth_m(),
        second.deepest_drawn_depth_m(),
        "both frames are drawn on the same depth axis"
    );

    // And on that shared axis the collapsed frame is visibly flatter: the
    // west-to-east drop across the panel falls by the same factor of ten the
    // field did.
    let drop_of = |view: &SideView| {
        let columns = view.columns();
        columns[0].interface_fraction().expect("finite")
            - columns[NX - 1].interface_fraction().expect("finite")
    };
    assert!(
        (drop_of(&first) - drop_of(&second) * 10.0).abs() < POSITION_TOLERANCE,
        "a tenth of the tilt is drawn a tenth as steep"
    );

    // Every interface of every frame fits inside the panel, which is what the
    // run-wide axis is for.
    for view in [&first, &second] {
        assert!(view
            .columns()
            .iter()
            .all(|column| column.interface_fraction().is_some_and(|f| f < 1.0)));
    }
}

#[test]
fn an_interface_that_is_not_a_number_is_not_drawn() {
    // A `NaN` in `h` means the integration diverged there. The view draws no
    // water at all in that column rather than an ocean at some depth the run
    // never produced — the same refusal `CrossSection` makes by breaking its
    // line.
    let mut field = tilt_field(NX, NY);
    for row in 0..NY {
        field[row * NX + 7] = f64::NAN;
    }
    let view = side_view_of_field(field);
    let columns = view.columns();
    assert!(columns[7].bands().is_none());
    assert!(columns[7].interface_fraction().is_none());
    assert!(columns[6].bands().is_some());
    assert!(columns[8].bands().is_some());
}

#[test]
fn the_side_view_reads_the_same_ocean_the_cross_section_does() {
    // The view is not a second extraction of the equator: it is the section's
    // own data, drawn as an ocean. So the two agree column for column, and the
    // latitude the section was read at is the latitude the side view reports.
    let header = steady_trades_header(1);
    let run = run_of(&header, |_| tilt_field(NX, NY));
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    let view = SideView::of_section(&section, MEAN_DEPTH_M);

    assert_eq!(view.latitude_deg_north(), section.latitude_deg_north());
    assert_eq!(view.mean_depth_m(), MEAN_DEPTH_M);
    assert_eq!(view.columns().len(), section.points().len());
    // Column for column and in the section's own order, so a place on one view
    // is the same place on the other, and the depth is that point's `h` under
    // `H`.
    for (column, point) in view.columns().iter().zip(section.points()) {
        assert!((column.interface_depth_m() - (MEAN_DEPTH_M + point.h_m())).abs() < f64::EPSILON);
    }
}

#[test]
fn a_column_where_the_upper_layer_has_vanished_is_not_drawn() {
    // `h ≤ −H` means the upper layer has no thickness left. The model is
    // linear and assumes anomalies small against `H` (`docs/planning/
    // 01-scientific-model.md`), so this is the model broken rather than an
    // ocean it computed — and a column drawn full of cold water would be a
    // picture of an outcropped thermocline, which is a claim the run does not
    // support. It is drawn as no reading, exactly as a `NaN` is.
    let mut field = tilt_field(NX, NY);
    for row in 0..NY {
        // A metre past the vanishing point, so the test is about the state and
        // not about the boundary.
        field[row * NX + 11] = -MEAN_DEPTH_M - 1.0;
        // And exactly at it: an interface at the sea surface is the same state.
        field[row * NX + 12] = -MEAN_DEPTH_M;
    }
    let view = side_view_of_field(field);
    let columns = view.columns();
    for index in [11, 12] {
        assert!(
            columns[index].bands().is_none(),
            "an upper layer of no thickness is drawn as no reading, not as cold water"
        );
        assert!(columns[index].interface_fraction().is_none());
    }
    assert!(columns[10].bands().is_some());
    assert!(columns[13].bands().is_some());
}
