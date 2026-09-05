//! T-08.3 acceptance criterion: dragging the scrubber updates the displayed
//! frame with no perceptible lag on a moderately-sized test run.
//!
//! Two things have to hold for that, and both are asserted here without a
//! window, because what a reader would check by dragging is a property of
//! [`LoadedRun::frame`]:
//!
//! - **The frame shown is the frame chosen**, whichever order the frames are
//!   asked for in. A drag arrives as a jump, not a walk: the reader lets go of
//!   the scrubber at frame 600 having passed through none of the frames in
//!   between, and the panel draws whichever frame it lands on.
//! - **Reaching a frame costs the same wherever it is in the run.** Forward-
//!   only reading (`termocline_format::reader`) decodes every frame before the
//!   one asked for, so a drag across a 731-frame run gets steadily slower the
//!   further right it goes — which is exactly the lag the criterion rules out.
//!
//! # Where the expected values come from
//!
//! The run under test is written here from `termocline_format` alone, with a
//! thermocline depth anomaly `h` that is the frame's own index in metres. The
//! frame `n` should hold is then known before the code is asked for it, rather
//! than read back off it. Its length is `engine/scenarios/steady-trades.toml`'s
//! — 17 520 steps of an hour with a frame every 24 makes 730 daily frames, and
//! the frame at t = 0 makes 731.

mod common;

use std::hint::black_box;
use std::time::{Duration, Instant};

use common::{
    encoded_frames_with_h, header_on, steady_trades_header, FRAME_INTERVAL_S, PACIFIC,
    STEADY_TRADES_FRAMES as FRAMES,
};
use termocline_format::{GridSpec, RunHeader, Variable};
use visualizer::{Heatmap, LoadedRun, RunBytes};

/// The run the timing is taken on is long, so its basin is small: what makes a
/// drag slow is how many frames stand before the one asked for, and a frame of
/// the scenario's 320 × 100 basin is 1.3 MB — 731 of them is a gigabyte no
/// test should hold.
const TIMING_GRID: (usize, usize) = (16, 8);

/// A run of `count` frames on `header`'s grid whose `h` is everywhere the
/// frame's own index, in metres.
fn run_of(header: &RunHeader, count: u64) -> LoadedRun {
    let cells = header.grid.field_len(Variable::ThermoclineDepthAnomaly);
    let frames = encoded_frames_with_h(header, count, |index| {
        #[allow(clippy::cast_precision_loss)]
        let h_m = index as f64;
        vec![h_m; cells]
    });
    LoadedRun::from_bytes(
        "frame-access",
        RunBytes {
            header: serde_json::to_vec(header).expect("a header serializes"),
            frames,
        },
    )
    .expect("a run written from its own header loads")
}

/// The run the timing is taken on: 731 frames of a small basin.
fn long_run() -> LoadedRun {
    let (nx, ny) = TIMING_GRID;
    let grid = GridSpec::new(nx, ny, PACIFIC).expect("16 x 8 is a valid basin");
    run_of(&header_on(grid, "frame-access", FRAMES), FRAMES)
}

/// The shortest of `samples` attempts at fetching frame `index`.
///
/// The shortest rather than the mean: every sample does the same work, so what
/// separates them is scheduling noise the machine added, and the smallest
/// sample is the one with least of it.
fn shortest_fetch(run: &LoadedRun, index: u64, samples: u32) -> Duration {
    (0..samples)
        .map(|_| {
            let start = Instant::now();
            black_box(run.frame(black_box(index)));
            start.elapsed()
        })
        .min()
        .expect("at least one sample")
}

#[test]
fn a_frame_asked_for_out_of_order_is_the_frame_that_comes_back() {
    let run = long_run();
    // A drag's worth of jumps: to the end, back to the start, and around the
    // middle in neither direction.
    for index in [730, 0, 365, 729, 1, 366, 730, 0] {
        let frame = run
            .frame(index)
            .expect("the run holds every frame it counts");
        #[allow(clippy::cast_precision_loss)]
        let expected_h_m = index as f64;
        assert!(
            frame.h().iter().all(|h_m| *h_m == expected_h_m),
            "frame {index} should hold h = {expected_h_m} m everywhere",
        );
        #[allow(clippy::cast_precision_loss)]
        let expected_t_s = index as f64 * FRAME_INTERVAL_S;
        assert!(
            (frame.t_s() - expected_t_s).abs() < f64::EPSILON,
            "frame {index} should be at t = {expected_t_s} s, not {}",
            frame.t_s(),
        );
    }
}

#[test]
fn the_frames_of_the_scenario_grid_are_reachable_in_any_order() {
    let header = steady_trades_header(4);
    let run = run_of(&header, 4);
    for index in [3, 0, 2, 1] {
        #[allow(clippy::cast_precision_loss)]
        let expected_h_m = index as f64;
        let frame = run
            .frame(index)
            .expect("the run holds every frame it counts");
        assert_eq!(frame.h().first().copied(), Some(expected_h_m));
    }
}

#[test]
fn past_the_end_of_the_run_there_is_no_frame() {
    let run = long_run();
    assert!(run.frame(FRAMES).is_none());
    assert!(run.frame(u64::MAX).is_none());
}

#[test]
fn reaching_the_last_frame_costs_what_reaching_the_first_costs() {
    let run = long_run();
    // Enough samples that the shortest of them is dominated by the work and
    // not by whatever else the machine was doing.
    const SAMPLES: u32 = 64;
    let first = shortest_fetch(&run, 0, SAMPLES);
    let last = shortest_fetch(&run, FRAMES - 1, SAMPLES);

    // Reading forward from the start of the run, the last frame costs 731
    // decodes to the first frame's one. Fetching one frame wherever it sits
    // costs one decode either way, so the honest ratio is 1. The bound sits an
    // order of magnitude above that — timing noise, and the frame bytes one
    // end of the run may have left warmer in cache than the other — and an
    // order of magnitude below the 731 a forward-only walk would cost.
    const RATIO_BOUND: u32 = 20;
    assert!(
        last <= first * RATIO_BOUND,
        "the last of {FRAMES} frames took {last:?} against the first frame's {first:?}: \
         reaching a frame still costs more the further into the run it is",
    );
}

/// Below this, a response to a direct manipulation feels instantaneous rather
/// than laggy: 0.1 s, the first of Nielsen's three response-time limits
/// (Nielsen, *Usability Engineering*, 1993, ch. 5, after Miller 1968). The
/// acceptance criterion's "no perceptible lag" is that limit.
const PERCEPTIBLE_LAG: Duration = Duration::from_millis(100);

/// Frames of the scenario's own 320 × 100 basin a test can afford: each is
/// about 1.3 MB decoded.
const SCENARIO_GRID_FRAMES: u64 = 8;

#[test]
fn a_frame_of_the_scenario_basin_is_fetched_and_coloured_inside_that_limit() {
    let header = steady_trades_header(SCENARIO_GRID_FRAMES);
    let run = run_of(&header, SCENARIO_GRID_FRAMES);
    let scale = run.anomaly_scale();
    // What a drag does for each frame it lands on, bar the texture upload,
    // which needs a GPU: fetch the frame, and colour-map the basin.
    let shortest = (0..8)
        .map(|_| {
            let start = Instant::now();
            let frame = black_box(run.frame(SCENARIO_GRID_FRAMES - 1)).expect("the last frame");
            let heatmap = Heatmap::of_frame(header.grid, &frame, scale).expect("it fits the grid");
            black_box(heatmap);
            start.elapsed()
        })
        .min()
        .expect("at least one sample");

    // The shortest sample, for the same reason as above: what the bound is
    // about is the work, and the bound sits far enough above it — the work is
    // one decode and one colour-map of 32 000 cells — that scheduling noise
    // does not decide the test.
    assert!(
        shortest <= PERCEPTIBLE_LAG,
        "a frame of the scenario basin took {shortest:?}, past the {PERCEPTIBLE_LAG:?} \
         under which a drag feels instantaneous",
    );
}
