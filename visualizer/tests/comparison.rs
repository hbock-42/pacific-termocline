//! T-09.5: what it means for two runs to be comparable, and the one colour
//! scale both of them are drawn on.
//!
//! The acceptance criterion of the ticket — the two panels stay frame-synced —
//! is asserted on the panel itself, in `src/app.rs`, because that is where the
//! frame chooser is. What is asserted here is the part that decides whether
//! the side-by-side says anything true: that neither run is renormalized onto
//! a scale of its own, and that two runs which are not comparable are refused
//! rather than drawn.
//!
//! # Where the expected values come from
//!
//! The tilt is T-07.4's measured equilibrium of `engine/scenarios/steady-trades.toml`:
//! `h` = +38.2 m at the western wall and -28.2 m at the eastern. It is an
//! *input* here, as it is in `tests/heatmap.rs`, never an expected output read
//! back off this code. The second run is that same tilt at half the amplitude,
//! which is what a comparison exists to show.
//!
//! The colours are ColorBrewer's 11-class `RdBu` diverging scheme (Brewer,
//! Harrower & Pennsylvania State University, <https://colorbrewer2.org>),
//! transcribed from the published table, and the value asserted at the half
//! amplitude is interpolated between two of its classes by hand — never read
//! back off the ramp this code builds.

mod common;

use common::{
    encoded_frames_with_h, header_on, EASTERN_WALL_H_M, FRAME_INTERVAL_S, NX, NY, PACIFIC,
    STEADY_TRADES_PARAMS, WESTERN_WALL_H_M,
};
use termocline_format::{BasinExtent, GridSpec, OutputTiming, RunHeader, Variable};
use visualizer::{Comparison, Difference, Heatmap, LoadedRun, Mismatch, RunBytes};

/// ColorBrewer 11-class `RdBu`, reversed so it runs from the most negative
/// anomaly to the most positive, as `tests/heatmap.rs` transcribes it.
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

/// The colour of the deepest anomaly a scale reaches.
const DEEPEST: [u8; 3] = RD_BU_11[10];

/// A rounding tolerance of one 8-bit level per channel.
///
/// The ramp interpolates in `f64` and quantizes to eight bits, so a colour
/// asserted from the published table may land one level either side of the
/// value computed by hand. Nothing wider: two levels would let a genuinely
/// different class pass.
const CHANNEL_TOLERANCE: u8 = 1;

/// The grid of the steady-trades basin, which both runs of a comparison share
/// unless a test is about them not sharing it.
fn steady_trades_grid() -> GridSpec {
    GridSpec::new(NX, NY, PACIFIC).expect("320 x 100 over the Pacific is a valid basin")
}

/// A run of `frame_count` frames whose every frame carries the steady-trades
/// tilt scaled by `amplitude`: `h` is `+38.2 · amplitude` m in the westernmost
/// column, `-28.2 · amplitude` m in the easternmost, and zero between.
///
/// `amplitude` is what a comparison is for: the same scenario forced harder or
/// damped further produces the same shape at a different size, and the point
/// of one scale across both panels is that the difference in size survives to
/// the screen.
fn tilted_run(header: &RunHeader, amplitude: f64) -> LoadedRun {
    let grid = header.grid;
    let (nx, ny) = (grid.nx(), grid.ny());
    let mut field = vec![0.0; grid.field_len(Variable::ThermoclineDepthAnomaly)];
    for j in 0..ny {
        field[j * nx] = WESTERN_WALL_H_M * amplitude;
        field[j * nx + nx - 1] = EASTERN_WALL_H_M * amplitude;
    }
    LoadedRun::from_bytes(
        header.scenario_description.clone(),
        RunBytes {
            header: serde_json::to_vec(header).expect("a header serializes"),
            frames: encoded_frames_with_h(header, header.output.frame_count, |_| field.clone()),
        },
    )
    .expect("a run written from its own header loads")
}

/// The header of a steady-trades run of `frame_count` frames, named `scenario`.
fn header(scenario: &str, frame_count: u64) -> RunHeader {
    header_on(steady_trades_grid(), scenario, frame_count)
}

/// The colour the westernmost cell of `run`'s first frame is drawn in, on
/// `scale`.
fn western_wall_color(run: &LoadedRun, scale: visualizer::DivergingScale) -> [u8; 3] {
    let frame = run.frame(0).expect("the run has a first frame");
    let heatmap = Heatmap::of_frame(run.header().grid, &frame, scale)
        .expect("a frame of a run fits that run's grid");
    heatmap
        .pixel(0, 0)
        .expect("the northwest corner is inside the basin")
}

/// Assert two colours agree channel for channel, within the quantization
/// tolerance.
fn assert_color_near(drawn: [u8; 3], expected: [u8; 3], what: &str) {
    for channel in 0..3 {
        assert!(
            drawn[channel].abs_diff(expected[channel]) <= CHANNEL_TOLERANCE,
            "{what}: drawn {drawn:?} against the expected {expected:?}"
        );
    }
}

#[test]
fn both_runs_are_drawn_on_one_scale_reaching_as_far_as_the_louder() {
    let full = tilted_run(&header("steady-trades", 3), 1.0);
    let half = tilted_run(&header("steady-trades-damped", 3), 0.5);
    let comparison = Comparison::of(&full, &half).expect("two runs on one grid are comparable");

    // The louder run's own scale, because it is the one that has to hold
    // everything both runs contain: +38.2 m is the largest magnitude in either.
    assert_eq!(comparison.scale().half_range_m(), WESTERN_WALL_H_M);
    // And it is the same scale on both sides — not each run's own, which is
    // the whole of what makes the two panels readable against each other.
    assert_eq!(comparison.scale(), full.anomaly_scale());
    assert_ne!(comparison.scale(), half.anomaly_scale());
}

#[test]
fn the_quieter_run_is_not_renormalized_up_to_the_louder_ones_colours() {
    let full = tilted_run(&header("steady-trades", 3), 1.0);
    let half = tilted_run(&header("steady-trades-damped", 3), 0.5);
    let comparison = Comparison::of(&full, &half).expect("two runs on one grid are comparable");

    // On its own scale the half-amplitude run's western wall is the deepest
    // colour the ramp has — indistinguishable from the full run's western
    // wall, which is the misreading this view exists to prevent.
    assert_color_near(
        western_wall_color(&half, half.anomaly_scale()),
        DEEPEST,
        "the quieter run on its own scale saturates",
    );

    // On the shared scale it is half way out: h/half_range = 0.5 puts it at
    // (0.5 + 1)/2 = 0.75 of the ramp, which is anchor 0.75 × 10 = 7.5 — half
    // way between RdBu classes 7 and 8. Interpolated by hand from the
    // published table: (244+214)/2 = 229, (165+96)/2 = 130.5, (130+77)/2 =
    // 103.5.
    assert_color_near(
        western_wall_color(&half, comparison.scale()),
        [229, 131, 104],
        "the quieter run on the shared scale reads as quieter",
    );
    // The louder run still reaches the end of the ramp, so the shared scale
    // costs it nothing.
    assert_color_near(
        western_wall_color(&full, comparison.scale()),
        DEEPEST,
        "the louder run on the shared scale still saturates",
    );
}

#[test]
fn each_runs_own_range_is_reported_beside_the_shared_one() {
    let full = tilted_run(&header("steady-trades", 3), 1.0);
    let half = tilted_run(&header("steady-trades-damped", 3), 0.5);
    let comparison = Comparison::of(&full, &half).expect("two runs on one grid are comparable");

    let ranges = comparison
        .differences()
        .into_iter()
        .find_map(|difference| match difference {
            Difference::AnomalyRange {
                left_half_range_m,
                right_half_range_m,
            } => Some((left_half_range_m, right_half_range_m)),
            _ => None,
        })
        .expect("two runs of different amplitudes differ in range");
    assert_eq!(ranges, (WESTERN_WALL_H_M, WESTERN_WALL_H_M * 0.5));
}

#[test]
fn runs_on_differently_shaped_grids_are_refused() {
    let left = header("steady-trades", 3);
    let mut right = header("steady-trades-coarse", 3);
    right.grid = GridSpec::new(NX / 2, NY / 2, PACIFIC).expect("160 x 50 is a valid basin");

    let refusal = Comparison::of(&tilted_run(&left, 1.0), &tilted_run(&right, 1.0))
        .expect_err("two resolutions of one basin are not the same cells");
    assert!(matches!(refusal, Mismatch::Grid { .. }));
    // The message names both sides: a refusal a reader cannot act on is worse
    // than the misleading picture it replaced.
    let message = refusal.to_string();
    assert!(message.contains("320"), "{message}");
    assert!(message.contains("160"), "{message}");
}

#[test]
fn runs_over_different_basins_are_refused() {
    let left = header("steady-trades", 3);
    let mut right = header("atlantic", 3);
    right.grid = GridSpec::new(NX, NY, BasinExtent::new(-60.0, 20.0, -25.0, 25.0))
        .expect("an Atlantic basin is a valid grid");

    let refusal = Comparison::of(&tilted_run(&left, 1.0), &tilted_run(&right, 1.0))
        .expect_err("the same cell count over another ocean is another place");
    assert!(matches!(refusal, Mismatch::Grid { .. }));
}

#[test]
fn runs_written_at_different_cadences_are_refused() {
    let left = header("steady-trades", 3);
    let mut right = header("steady-trades-hourly", 3);
    right.output = OutputTiming {
        frame_count: 3,
        // A frame an hour against a frame a day: frame 2 of one run and frame
        // 2 of the other are 23 hours of model time apart, so a synced index
        // is not a synced time.
        interval_s: 3_600.0,
    };

    let refusal = Comparison::of(&tilted_run(&left, 1.0), &tilted_run(&right, 1.0))
        .expect_err("a synced index only syncs time when the cadence matches");
    assert!(matches!(refusal, Mismatch::FrameInterval { .. }));
    let message = refusal.to_string();
    assert!(message.contains("86400"), "{message}");
    assert!(message.contains("3600"), "{message}");
}

#[test]
fn runs_of_different_lengths_are_compared_over_the_frames_they_share() {
    let long = tilted_run(&header("steady-trades", 731), 1.0);
    let short = tilted_run(&header("steady-trades-year", 366), 1.0);
    let comparison = Comparison::of(&long, &short).expect("a shorter run is still comparable");

    // 366 frames, not 731: past the end of the shorter run there is nothing to
    // put in the second panel, and holding its last frame while the other ran
    // on would draw a steady ocean the run never produced.
    assert_eq!(comparison.frame_count(), 366);
    let stated = comparison.differences().into_iter().any(|difference| {
        matches!(
            difference,
            Difference::Length {
                left_frames: 731,
                right_frames: 366,
                shared_frames: 366,
            }
        )
    });
    assert!(stated, "the shortened comparison is stated, not silent");
}

#[test]
fn a_coupled_run_and_an_uncoupled_one_are_compared_on_h_and_told_apart() {
    // The side-by-side draws the thermocline depth anomaly, and every run
    // carries it (`Variable::LINEAR_CORE`) — so a run that couples SST and one
    // that does not are comparable in exactly the field on screen. That the
    // one carries `T'` and the other does not is a difference between them,
    // not a reason to refuse.
    let uncoupled_header = header("steady-trades", 3);
    let coupled_header = header("steady-trades-coupled", 3).with_sst_anomaly();
    let uncoupled = tilted_run(&uncoupled_header, 1.0);
    let coupled = coupled_run(&coupled_header, 1.0);

    let comparison =
        Comparison::of(&uncoupled, &coupled).expect("both runs carry the field on screen");
    assert_eq!(comparison.frame_count(), 3);
    let stated = comparison.differences().into_iter().any(|difference| {
        matches!(
            difference,
            Difference::SstAnomaly {
                left: false,
                right: true
            }
        )
    });
    assert!(stated, "a coupled run beside an uncoupled one says so");
}

#[test]
fn differing_physical_parameters_are_named() {
    let left = header("steady-trades", 3);
    let mut right = header("steady-trades-stratified", 3);
    // A stronger stratification: c = √(g'H) rises with g', so the two runs are
    // of the same basin with a different wave speed — a comparison worth
    // drawing, and one whose cause a reader should be told.
    right.physical_params.reduced_gravity_m_per_s2 = 0.09;

    let (left, right) = (tilted_run(&left, 1.0), tilted_run(&right, 1.0));
    let comparison = Comparison::of(&left, &right)
        .expect("a different g' does not stop two runs being comparable");
    let named = comparison
        .differences()
        .into_iter()
        .find_map(|difference| match difference {
            Difference::PhysicalParam {
                name, left, right, ..
            } => Some((name, left, right)),
            _ => None,
        })
        .expect("the parameter that differs is named");
    assert_eq!(
        named,
        (
            "Reduced gravity g'",
            STEADY_TRADES_PARAMS.reduced_gravity_m_per_s2,
            0.09
        )
    );
}

#[test]
fn two_runs_of_one_scenario_differ_in_nothing_worth_stating() {
    let header = header("steady-trades", 3);
    let left = tilted_run(&header, 1.0);
    let right = tilted_run(&header, 1.0);
    let comparison = Comparison::of(&left, &right).expect("a run is comparable with itself");
    assert!(
        comparison.differences().is_empty(),
        "{:?}",
        comparison.differences()
    );
    assert_eq!(comparison.frame_count(), 3);
}

/// A coupled run: the tilt of [`tilted_run`], with a mixed-layer SST anomaly
/// alongside it so the frames match a `with_sst_anomaly` header.
fn coupled_run(header: &RunHeader, amplitude: f64) -> LoadedRun {
    use termocline_format::{frame_encoding, Frame};

    let grid = header.grid;
    let (nx, ny) = (grid.nx(), grid.ny());
    let mut h_m = vec![0.0; grid.field_len(Variable::ThermoclineDepthAnomaly)];
    for j in 0..ny {
        h_m[j * nx] = WESTERN_WALL_H_M * amplitude;
        h_m[j * nx + nx - 1] = EASTERN_WALL_H_M * amplitude;
    }
    let zero = |variable| vec![0.0; grid.field_len(variable)];
    let mut frames = Vec::new();
    for index in 0..header.output.frame_count {
        #[allow(clippy::cast_precision_loss)]
        let t_s = index as f64 * FRAME_INTERVAL_S;
        let frame = Frame::new(
            t_s,
            &grid,
            h_m.clone(),
            zero(Variable::ZonalCurrentAnomaly),
            zero(Variable::MeridionalCurrentAnomaly),
            zero(Variable::ZonalWindStress),
            zero(Variable::MeridionalWindStress),
        )
        .expect("fields sized from the grid fit it")
        .with_sst_anomaly(&grid, zero(Variable::SstAnomaly))
        .expect("an SST field sized from the grid fits it");
        frames.extend(
            bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
        );
    }
    LoadedRun::from_bytes(
        header.scenario_description.clone(),
        RunBytes {
            header: serde_json::to_vec(header).expect("a header serializes"),
            frames,
        },
    )
    .expect("a coupled run written from its own header loads")
}
