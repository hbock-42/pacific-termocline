//! T-09.4 acceptance criterion: a point near the eastern boundary of a
//! wind-burst run shows the delayed thermocline signal, arriving after the
//! western perturbation.
//!
//! The criterion is written as something a reader sees, and something nobody
//! runs is not a criterion. So the delay is measured instead of looked at:
//! [`PointSeries`] is the same list of samples the shell draws the chart from,
//! built without a window or a device, so "the signal arrives later in the
//! east" becomes an assertion about *when* the series peaks ([ADR-0006]).
//!
//! # What these tests do and do not pin
//!
//! They are tests of the **view**, not of the physics. That a Kelvin wave
//! crosses the basin at `c` is the engine's claim, and it is validated against
//! the engine in `engine/tests/kelvin_wave.rs`; the visualizer must not link
//! the engine at all (ADR-0001), so no run here comes out of the solver.
//!
//! What is pinned here is that the view *transports* the timing that is in the
//! run it was given: the frames carry a pulse whose crest reaches each cell at
//! a stated instant, and the series must peak at that instant, at that cell,
//! and not at another — which is what a view that transposed its axes, read
//! the wrong column, or dropped frames would fail. The fixture and the
//! expected delay share `c`, so the assertion is not evidence about `c`; it is
//! evidence that the arrival in the file survives the trip to the chart.
//!
//! # Where the expected values come from
//!
//! The fixture is the analytic signal, and the assertions are about the view.
//!
//! A westerly wind burst launches an equatorially trapped **Kelvin wave**,
//! which travels **eastward only** and is **non-dispersive** at
//! `c = √(g'·H)` (`CONTEXT.md`, *Kelvin wave*, *Kelvin wave speed*). With
//! `engine/scenarios/wind-burst.toml`'s `g' = 0.06 m s^-2` and `H = 150 m`
//! that is `c = √9 = 3.0 m s^-1` exactly — the observed first-baroclinic speed
//! of the equatorial Pacific. So a downwelling pulse released at the burst's
//! centre `x_b` at its peak time `t_peak` has its crest at
//!
//! ```text
//! x_crest(t) = x_b + c·(t − t_peak),
//! ```
//!
//! and a cell at `x` sees the crest at `t_peak + (x − x_b)/c`. The fixture
//! writes exactly that pulse — a Gaussian in `x − x_crest(t)` on the steady
//! trade-wind tilt of T-07.4 — and the tests assert that the series the view
//! builds carries the arrival times that closed form predicts. The travel time
//! is never read back out of the code that produced the frames: it is
//! `(x_east − x_west)/c`, arithmetic on the scenario file's own numbers.
//!
//! Every scenario constant used here — the burst's centre, width and peak
//! time, and the basin's bounds and resolution — is transcribed from
//! `engine/scenarios/wind-burst.toml`; the degrees-to-metres projection is
//! `engine/src/basin.rs`'s [`METRES_PER_DEGREE_OF_ARC`], transcribed for the
//! same reason (the visualizer must not link the engine, ADR-0001).
//!
//! # The SST anomaly, present and absent
//!
//! Frames may carry the Epic 12 mixed-layer SST anomaly `T'`, and an
//! uncoupled run's frames carry no such field at all (T-05.4,
//! `termocline_format::Frame`). Both cases are asserted, because the one way
//! this view could lie is by reporting an absent `T'` as 0 K — a claim that
//! the mixed layer sat at its climatological temperature, which an uncoupled
//! run never made.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod common;

use std::f64::consts::PI;

use common::{
    encoded_frames_with_fields, header_on, FrameFields, EASTERN_WALL_H_M, FRAME_INTERVAL_S, NX,
    STEADY_TRADES_FRAMES, WESTERN_WALL_H_M,
};
use termocline_format::{BasinExtent, GridSpec, RunHeader, Variable};
use visualizer::{BasinPoint, LoadedRun, PointSeries, RunBytes, SeriesSample};

/// Earth's mean radius, in metres, as `engine/src/basin.rs` states it.
const EARTH_MEAN_RADIUS_M: f64 = 6_371_008.8;

/// Metres per degree of arc at Earth's mean radius: `R·π/180`, the projection
/// `engine/src/basin.rs` uses to turn the basin's degrees into metres.
const METRES_PER_DEGREE_OF_ARC: f64 = EARTH_MEAN_RADIUS_M * PI / 180.0;

/// The zonal span of the basin, in degrees: 120°E to 80°W (`CONTEXT.md`,
/// *Basin*).
const BASIN_WIDTH_DEG: f64 = 160.0;

/// The basin's width, in metres. About 17 790 km, which is the "basin 17 800 km
/// wide" `wind-burst.toml` describes.
const BASIN_WIDTH_M: f64 = BASIN_WIDTH_DEG * METRES_PER_DEGREE_OF_ARC;

/// The Kelvin wave speed of `wind-burst.toml`, in metres per second:
/// `c = √(g'·H) = √(0.06 · 150) = 3.0` (`CONTEXT.md`, *Kelvin wave speed*).
const KELVIN_SPEED_M_PER_S: f64 = 3.0;

/// `wind-burst.toml`'s `center_x_m`: the burst sits 2 000 km east of the
/// western boundary, in the warm pool.
const BURST_CENTER_X_M: f64 = 2.0e6;

/// `wind-burst.toml`'s `zonal_scale_m`: the burst is 1 000 km wide, and so is
/// the pulse it launches.
const BURST_ZONAL_SCALE_M: f64 = 1.0e6;

/// `wind-burst.toml`'s `peak_time_s`: the burst fires one tropical year into
/// the run, once the trade-driven tilt has spun up.
const BURST_PEAK_TIME_S: f64 = 31_556_926.08;

/// The crest height of the downwelling pulse the fixture writes, in metres.
///
/// The fixture's own, not a measured value, and nothing below depends on it:
/// every assertion here is about *when* the pulse arrives, and the ones about
/// how far it moves the thermocline are stated as fractions of this.
const PULSE_AMPLITUDE_M: f64 = 12.0;

/// The trough depth of the upwelling pulse, in metres.
///
/// The mirror of [`PULSE_AMPLITUDE_M`]. The criterion asks for the
/// "shoaling/deepening" signal, and both halves happen: a westerly burst
/// launches a downwelling Kelvin wave, an easterly anomaly against the trades
/// an upwelling one, and a view that only got the deepening right would be
/// half a view.
const UPWELLING_AMPLITUDE_M: f64 = -12.0;

/// The equatorial band the fixture is written on: the scenario's zonal grid,
/// over the one degree of latitude straddling the equator.
///
/// The scenario basin is 320 × 100 cells, and 731 frames of it is a run far
/// larger than a test should build. Two rows at ±0.25° are where the
/// scenario's own equatorial rows sit — the axis of the waveguide a Kelvin
/// wave travels along — and a series is read at one cell, so the rows this
/// leaves out are not what is under test.
const WAVEGUIDE: BasinExtent = BasinExtent::new(120.0, -80.0, -0.5, 0.5);

/// Rows of the waveguide fixture.
const WAVEGUIDE_NY: usize = 2;

/// How far a sampled anomaly may sit from the analytic value it was written
/// from, in metres.
///
/// The series reads a value out of a decoded frame and does no arithmetic on
/// it, so the only error is the `bincode` round-trip, which is exact for
/// `f64`. A picometre is far below that and still fifteen orders below the
/// metre-scale structure the criterion is about.
const SAMPLE_TOLERANCE_M: f64 = 1e-12;

/// How far a measured arrival time may sit from the closed form's, in seconds.
///
/// One output interval. The series is sampled once a day, so the crest of a
/// pulse can only be located to the frame it is nearest — half an interval
/// either side — and the delay is a difference of two such locations, so its
/// error is bounded by a whole one. The travel time being measured is 61 days,
/// sixty times this bound.
const ARRIVAL_TOLERANCE_S: f64 = FRAME_INTERVAL_S;

/// How far a plot position may sit from the fraction the geometry gives.
/// Dimensionless, and round-off on quantities of order one.
const POSITION_TOLERANCE: f64 = 1e-12;

/// The steady thermocline tilt `fraction` of the way east across the basin, in
/// metres: T-07.4's measured equilibrium at the two walls with a straight line
/// between them.
///
/// The straight line is a stand-in, as `tests/cross_section.rs` says of the
/// same one — the measured profile is curved. It carries the only property
/// these tests need of it: it does not change with time, so a peak in time at
/// a fixed cell is the pulse and nothing else.
fn steady_tilt_m(fraction: f64) -> f64 {
    WESTERN_WALL_H_M + (EASTERN_WALL_H_M - WESTERN_WALL_H_M) * fraction
}

/// The analytic thermocline depth anomaly at `x_m` east of the western wall at
/// model time `t_s`, in metres: the steady tilt plus the eastward-travelling
/// Kelvin pulse of crest height `amplitude_m` the burst released.
fn analytic_h_m(x_m: f64, t_s: f64, amplitude_m: f64) -> f64 {
    let crest_x_m = KELVIN_SPEED_M_PER_S.mul_add(t_s - BURST_PEAK_TIME_S, BURST_CENTER_X_M);
    let offset = (x_m - crest_x_m) / BURST_ZONAL_SCALE_M;
    amplitude_m.mul_add((-offset * offset).exp(), steady_tilt_m(x_m / BASIN_WIDTH_M))
}

/// How far east of the western wall the centre of column `column` sits, in
/// metres.
fn cell_center_x_m(column: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let fraction = (column as f64 + 0.5) / NX as f64;
    fraction * BASIN_WIDTH_M
}

/// The waveguide grid the fixture is written on.
fn waveguide_grid() -> GridSpec {
    GridSpec::new(NX, WAVEGUIDE_NY, WAVEGUIDE).expect("320 x 2 is a valid basin")
}

/// The header of the wind-burst fixture run.
fn burst_header() -> RunHeader {
    header_on(waveguide_grid(), "wind-burst", STEADY_TRADES_FRAMES)
}

/// The wind-burst run itself: 731 daily frames of the analytic pulse on the
/// steady tilt, uniform across the two rows of the waveguide.
///
/// The wind stress is left calm. What the burst does to the ocean is written
/// into `h` directly, from the closed form, rather than integrated — this is a
/// test of a view, and a view that only agreed with the run when the run came
/// from the solver would be a test of the solver.
fn burst_run() -> LoadedRun {
    pulse_run(PULSE_AMPLITUDE_M)
}

/// The same run carrying a pulse of crest height `amplitude_m`: positive for
/// the downwelling wave a westerly burst launches, negative for the upwelling
/// one that shoals the thermocline as it passes.
fn pulse_run(amplitude_m: f64) -> LoadedRun {
    let header = burst_header();
    let grid = header.grid;
    let (width, _height) = grid
        .grid()
        .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
    let cells = grid.field_len(Variable::ThermoclineDepthAnomaly);
    let frames = encoded_frames_with_fields(&header, header.output.frame_count, |index| {
        #[allow(clippy::cast_precision_loss)]
        let t_s = index as f64 * FRAME_INTERVAL_S;
        FrameFields {
            h_m: (0..cells)
                .map(|offset| analytic_h_m(cell_center_x_m(offset % width), t_s, amplitude_m))
                .collect(),
            ..FrameFields::calm(&header)
        }
    });
    LoadedRun::from_bytes(
        "wind-burst",
        RunBytes {
            header: serde_json::to_vec(&header).expect("a header serializes"),
            frames,
        },
    )
    .expect("a run written from its own header loads")
}

/// The cell at `column` and `row`, picked the way a reader picks one: off the
/// map, by the fraction of it the click landed at.
///
/// Row 0 of the map is the northernmost (`visualizer::Heatmap`), so the row
/// asked for here — counted north from the southern wall, as the field counts
/// it — is `rows - 1 - row` rows down the map.
fn point_at(grid: GridSpec, column: usize, row: usize) -> BasinPoint {
    let (width, height) = grid
        .grid()
        .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
    #[allow(clippy::cast_precision_loss)]
    let east = (column as f64 + 0.5) / width as f64;
    #[allow(clippy::cast_precision_loss)]
    let down = ((height - 1 - row) as f64 + 0.5) / height as f64;
    BasinPoint::at_map_fraction(grid, east, down).expect("a fraction inside the map names a cell")
}

/// The model time at which `series` is at its deepest, in seconds.
///
/// The pulse is the only thing in the fixture that changes with time at a
/// fixed cell, so its crest is where the series peaks.
fn deepest_t_s(series: &PointSeries) -> f64 {
    series
        .samples()
        .iter()
        .max_by(|a, b| a.h_m().total_cmp(&b.h_m()))
        .expect("the series has samples")
        .t_s()
}

/// The sample of `series` nearest model time `t_s`.
fn sample_nearest(series: &PointSeries, t_s: f64) -> SeriesSample {
    *series
        .samples()
        .iter()
        .min_by(|a, b| (a.t_s() - t_s).abs().total_cmp(&(b.t_s() - t_s).abs()))
        .expect("the series has samples")
}

#[test]
fn a_click_on_the_map_names_the_cell_under_it() {
    let grid = waveguide_grid();
    // The north-west corner of the map is the first column of the northernmost
    // row, and the northernmost row is the *last* row of the field.
    let north_west = BasinPoint::at_map_fraction(grid, 0.0, 0.0).expect("the corner is on the map");
    assert_eq!(north_west.column(), 0);
    assert_eq!(north_west.row(), WAVEGUIDE_NY - 1);
    assert!(north_west.latitude_deg_north() > 0.0, "the north half");
    assert!(
        (north_west.longitude_deg_east() - (120.0 + BASIN_WIDTH_DEG / NX as f64 / 2.0)).abs()
            < 1e-9,
        "the centre of the westernmost column"
    );

    // And the south-east corner is the last column of the field's first row.
    let south_east =
        BasinPoint::at_map_fraction(grid, 0.999, 0.999).expect("the corner is on the map");
    assert_eq!(south_east.column(), NX - 1);
    assert_eq!(south_east.row(), 0);
    assert!(south_east.latitude_deg_north() < 0.0, "the south half");

    // A click off the map names nothing: the reader who missed it did not mean
    // the coast.
    assert_eq!(BasinPoint::at_map_fraction(grid, 1.0, 0.5), None);
    assert_eq!(BasinPoint::at_map_fraction(grid, -0.01, 0.5), None);
    assert_eq!(BasinPoint::at_map_fraction(grid, 0.5, f64::NAN), None);
}

#[test]
fn the_series_holds_one_sample_per_frame_of_the_run() {
    let run = burst_run();
    let series = PointSeries::at_point(&run, point_at(run.header().grid, NX - 1, 0))
        .expect("a cell of this run's own basin");
    assert_eq!(series.samples().len() as u64, STEADY_TRADES_FRAMES);
    for (index, sample) in series.samples().iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let expected_s = index as f64 * FRAME_INTERVAL_S;
        assert!((sample.t_s() - expected_s).abs() < f64::EPSILON * expected_s.max(1.0));
    }
}

#[test]
fn the_series_reads_the_cell_the_click_named() {
    let run = burst_run();
    let column = NX - 1;
    let series = PointSeries::at_point(&run, point_at(run.header().grid, column, 0))
        .expect("a cell of this run's own basin");
    let x_m = cell_center_x_m(column);
    for sample in series.samples() {
        assert!(
            (sample.h_m() - analytic_h_m(x_m, sample.t_s(), PULSE_AMPLITUDE_M)).abs()
                < SAMPLE_TOLERANCE_M,
            "the series carries the run's own values at that cell"
        );
    }
}

#[test]
fn an_eastern_point_sees_the_burst_a_basin_crossing_after_the_west() {
    // The acceptance criterion. A point near the eastern boundary deepens, and
    // it deepens later than the west by exactly the time a Kelvin wave takes to
    // cross the water between them.
    let run = burst_run();
    let grid = run.header().grid;
    // The burst sits 2 000 km east of the western wall; this is the column
    // whose centre is nearest it.
    let west_column = (0..NX)
        .min_by(|&a, &b| {
            (cell_center_x_m(a) - BURST_CENTER_X_M)
                .abs()
                .total_cmp(&(cell_center_x_m(b) - BURST_CENTER_X_M).abs())
        })
        .expect("the basin has columns");
    let east_column = NX - 1;

    let west = PointSeries::at_point(&run, point_at(grid, west_column, 0))
        .expect("a cell of this run's own basin");
    let east = PointSeries::at_point(&run, point_at(grid, east_column, 0))
        .expect("a cell of this run's own basin");

    let measured_delay_s = deepest_t_s(&east) - deepest_t_s(&west);
    assert!(
        measured_delay_s > 0.0,
        "the eastern signal arrives after the western perturbation, not with it"
    );
    // The closed form: a non-dispersive Kelvin wave covers the water between
    // the two cells at c, and nothing else in the fixture moves.
    let expected_delay_s =
        (cell_center_x_m(east_column) - cell_center_x_m(west_column)) / KELVIN_SPEED_M_PER_S;
    assert!(
        (measured_delay_s - expected_delay_s).abs() < ARRIVAL_TOLERANCE_S,
        "the eastern point peaks {measured_delay_s} s after the west; \
         a Kelvin wave crosses that water in {expected_delay_s} s"
    );
}

#[test]
fn the_eastern_point_is_undisturbed_until_the_wave_reaches_it() {
    // The other half of the criterion: the signal is *delayed*. At the instant
    // the burst fires in the west, the eastern point is still sitting on the
    // steady tilt; a basin crossing later it has deepened by the pulse.
    let run = burst_run();
    let east_column = NX - 1;
    let series = PointSeries::at_point(&run, point_at(run.header().grid, east_column, 0))
        .expect("a cell of this run's own basin");
    let x_m = cell_center_x_m(east_column);
    let tilt_m = steady_tilt_m(x_m / BASIN_WIDTH_M);

    let when_it_fires = sample_nearest(&series, BURST_PEAK_TIME_S);
    // The pulse is 15.8 e-folding widths away at that instant, so the Gaussian
    // there is exp(-249) — zero to any precision an f64 keeps. A millimetre is
    // an enormously looser bound than that and still far below the 12 m the
    // arrival is worth.
    assert!(
        (when_it_fires.h_m() - tilt_m).abs() < 1e-3,
        "the east has not moved when the burst fires in the west"
    );

    let crossing_s = (x_m - BURST_CENTER_X_M) / KELVIN_SPEED_M_PER_S;
    let on_arrival = sample_nearest(&series, BURST_PEAK_TIME_S + crossing_s);
    // Sampled daily against a pulse 3.9 days wide (its 1 000 km width at
    // 3 m s^-1), the nearest frame to the crest is within half a day of it, so
    // the crest is caught at exp(-(0.5/3.9)^2) = 98% of full height at worst.
    assert!(
        on_arrival.h_m() - tilt_m > 0.9 * PULSE_AMPLITUDE_M,
        "a basin crossing later the thermocline there has deepened"
    );
}

#[test]
fn an_upwelling_pulse_arrives_in_the_east_as_a_delayed_shoaling() {
    // The other half of "shoaling/deepening". The same geometry with the
    // pulse's sign reversed: the eastern point rises *shallower* than the
    // steady tilt, and it does so a basin crossing after the perturbation, not
    // with it.
    let run = pulse_run(UPWELLING_AMPLITUDE_M);
    let east_column = NX - 1;
    let series = PointSeries::at_point(&run, point_at(run.header().grid, east_column, 0))
        .expect("a cell of this run's own basin");
    let x_m = cell_center_x_m(east_column);
    let tilt_m = steady_tilt_m(x_m / BASIN_WIDTH_M);

    let shallowest = series
        .samples()
        .iter()
        .min_by(|a, b| a.h_m().total_cmp(&b.h_m()))
        .expect("the series has samples");
    // Caught within half a day of the crest of a pulse 3.9 days wide, as in
    // the deepening test above, so at worst 98% of full depth.
    assert!(
        tilt_m - shallowest.h_m() > 0.9 * UPWELLING_AMPLITUDE_M.abs(),
        "the thermocline there shoals as the wave passes"
    );
    let expected_s = BURST_PEAK_TIME_S + (x_m - BURST_CENTER_X_M) / KELVIN_SPEED_M_PER_S;
    assert!(
        shallowest.t_s() > BURST_PEAK_TIME_S,
        "the shoaling is delayed, not simultaneous with the burst"
    );
    assert!(
        (shallowest.t_s() - expected_s).abs() < ARRIVAL_TOLERANCE_S,
        "and it arrives a basin crossing later"
    );
}

#[test]
fn the_axis_reaches_the_series_own_extreme_and_puts_the_run_at_both_ends() {
    let run = burst_run();
    let series = PointSeries::at_point(&run, point_at(run.header().grid, NX - 1, 0))
        .expect("a cell of this run's own basin");
    let deepest = series
        .samples()
        .iter()
        .max_by(|a, b| a.h_m().total_cmp(&b.h_m()))
        .expect("the series has samples");
    let shallowest = series
        .samples()
        .iter()
        .min_by(|a, b| a.h_m().total_cmp(&b.h_m()))
        .expect("the series has samples");
    let (_east, deepest_down) = series
        .plot_position(deepest)
        .expect("a finite anomaly has a position");
    let (_east, shallowest_down) = series
        .plot_position(shallowest)
        .expect("a finite anomaly has a position");
    // `y` is down, so a deeper-than-average anomaly is drawn *above* a
    // shallower one — the same way round as the map's warm-for-deep colours.
    assert!(deepest_down < shallowest_down);
    // The axis is symmetric about zero and reaches this series' own largest
    // magnitude. At a cell near the eastern boundary that is the shallow
    // steady tilt, which therefore sits exactly on the bottom of the chart.
    assert!(shallowest.h_m() < 0.0, "the eastern thermocline is shallow");
    assert!((shallowest_down - 1.0).abs() < POSITION_TOLERANCE);

    let first = series.samples().first().expect("the series has samples");
    let last = series.samples().last().expect("the series has samples");
    assert!(series.time_fraction(first.t_s()).abs() < POSITION_TOLERANCE);
    assert!((series.time_fraction(last.t_s()) - 1.0).abs() < POSITION_TOLERANCE);
}

/// A small basin for the SST tests: nothing about `T'` needs 731 frames of the
/// scenario grid.
fn small_grid() -> GridSpec {
    GridSpec::new(4, 3, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
        .expect("a 4x3 basin is a valid grid")
}

/// A three-frame run on [`small_grid`], carrying `T'` if `sst_anomaly_k` gives
/// a field for a frame.
fn small_run(sst_anomaly_k: impl Fn(u64) -> Option<Vec<f64>>, couples_sst: bool) -> LoadedRun {
    let header = header_on(small_grid(), "sst", 3);
    let header = if couples_sst {
        header.with_sst_anomaly()
    } else {
        header
    };
    let frames =
        encoded_frames_with_fields(&header, header.output.frame_count, |index| FrameFields {
            sst_anomaly_k: sst_anomaly_k(index),
            ..FrameFields::calm(&header)
        });
    LoadedRun::from_bytes(
        "sst",
        RunBytes {
            header: serde_json::to_vec(&header).expect("a header serializes"),
            frames,
        },
    )
    .expect("a run written from its own header loads")
}

#[test]
fn a_point_picked_off_another_basin_is_refused_rather_than_read() {
    // A cell is an index pair *into a particular basin*. Asked of a run whose
    // basin is a different shape, the same pair names a different place
    // entirely, so it is refused rather than answered with a plausible-looking
    // series of somewhere else.
    let elsewhere = BasinPoint::at_map_fraction(small_grid(), 0.5, 0.5)
        .expect("the middle of a map is on the map");
    let run = burst_run();
    assert!(PointSeries::at_point(&run, elsewhere).is_err());
}

#[test]
fn an_uncoupled_run_reports_no_sst_anomaly_rather_than_zero() {
    // The honesty test. An uncoupled run's frames carry no `T'` at all, and
    // the series must say so — not report 0 K, which would claim the mixed
    // layer sat at its climatological temperature for the whole run.
    let run = small_run(|_| None, false);
    let series = PointSeries::at_point(&run, point_at(run.header().grid, 1, 1))
        .expect("a cell of this run's own basin");
    assert!(!series.carries_sst_anomaly());
    assert_eq!(series.sst_scale(), None);
    for sample in series.samples() {
        assert_eq!(sample.sst_anomaly_k(), None);
        assert_eq!(series.sst_plot_position(sample), None);
    }
}

#[test]
fn a_coupled_run_carries_its_sst_anomaly_at_the_chosen_cell() {
    // `T'` is cell-centred like `h`, so the cell the reader picked is the cell
    // the second line is read at. The field is the cell's own offset in kelvin
    // plus the frame's index, so a series that read the wrong cell or the wrong
    // frame would carry the wrong numbers rather than merely the wrong shape.
    let grid = small_grid();
    let cells = grid.field_len(Variable::SstAnomaly);
    let run = small_run(
        |index| {
            #[allow(clippy::cast_precision_loss)]
            Some(
                (0..cells)
                    .map(|offset| offset as f64 + index as f64 / 10.0)
                    .collect(),
            )
        },
        true,
    );
    let (column, row) = (1, 1);
    let series = PointSeries::at_point(&run, point_at(grid, column, row))
        .expect("a cell of this run's own basin");
    assert!(series.carries_sst_anomaly());
    let (width, _height) = grid.grid().field_shape(Variable::SstAnomaly.staggering());
    #[allow(clippy::cast_precision_loss)]
    let offset_k = (row * width + column) as f64;
    for (index, sample) in series.samples().iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let expected_k = offset_k + index as f64 / 10.0;
        assert_eq!(sample.sst_anomaly_k(), Some(expected_k));
        assert!(series.sst_plot_position(sample).is_some());
    }
    // The axis reaches the largest magnitude in the series, which is the last
    // frame's value at that cell.
    #[allow(clippy::cast_precision_loss)]
    let largest_k = offset_k + 2.0 / 10.0;
    assert_eq!(
        series
            .sst_scale()
            .expect("a coupled run has one")
            .half_range_k(),
        largest_k
    );
}
