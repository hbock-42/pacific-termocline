//! Acceptance tests for the engine half of T-08.6 — a run taken a step at a
//! time.
//!
//! ADR-0012 puts the loop in the browser: a tab cannot block for the minutes a
//! full run takes, so it steps in chunks and yields between them. What that
//! costs is stated here — a run stepped in chunks has to be the same run as
//! one stepped straight through, or the browser and the CLI are two engines
//! wearing one name.
//!
//! The fs-free half of this file is what CI runs with `--no-default-features`,
//! so it names nothing behind the `fs` feature. The last test is native: it
//! holds the loop against the `run` command's own output, byte for byte.

use engine::{frame_of, RunLoop, Scenario};

/// A run small enough to step several times over in a test, held as text
/// rather than in a file — the shape the browser has (ADR-0012).
///
/// A 40° × 10° box at 1° under a steady easterly, 12 steps of an hour, a frame
/// every 4 steps.
const SCENARIO_TOML: &str = r#"
[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -5.0
northern_latitude_deg = 5.0
resolution_deg = 1.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 12
output_every_n_steps = 4

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0
"#;

// The values of `SCENARIO_TOML`, read off it by eye.
const DT_S: f64 = 3600.0;
const TOTAL_STEPS: u64 = 12;
const EVERY_N_STEPS: u64 = 4;
/// Frames at steps 0, 4, 8 and 12.
const EXPECTED_FRAME_COUNT: u64 = 4;

/// The scenario, built from text.
fn scenario() -> Scenario {
    Scenario::from_toml(SCENARIO_TOML).expect("the scenario text is valid")
}

/// Every frame of a run driven to the end of its schedule, `chunk` steps
/// between yields — which is what a browser's frame budget buys it.
///
/// The frames come back encoded, because that is the comparison that matters:
/// two runs that agree to the last bit of every field are the same run, and
/// nothing weaker distinguishes a chunked run from a broken one.
fn frames_of_run_in_chunks(chunk: u64) -> Vec<Vec<u8>> {
    let scenario = scenario();
    let mut run = RunLoop::of_scenario(&scenario, "chunked").expect("the scenario runs");
    let grid = run.header().grid;
    let mut frames = Vec::new();
    loop {
        // The yield a browser takes here: the loop is left, the tab draws, and
        // the next call resumes where this one stopped.
        let mut stepped = 0;
        while stepped < chunk {
            if let Some(saved) = run.take_frame() {
                let frame =
                    frame_of(saved.t_s, &grid, saved.state, saved.wind_stress).expect("a frame");
                frames.push(
                    bincode::serde::encode_to_vec(&frame, termocline_format::frame_encoding())
                        .expect("a frame encodes"),
                );
            }
            if !run.take_step() {
                return frames;
            }
            stepped += 1;
        }
    }
}

/// The loop hands out the frames the schedule promises, at the times it puts
/// them at.
///
/// Both expectations come from the schedule rather than from the run: frame
/// `k` is the state after `k · N` steps, at `k · N · dt` seconds.
#[test]
fn a_run_hands_out_the_frames_its_header_promises() {
    let scenario = scenario();
    let mut run = RunLoop::of_scenario(&scenario, "promised").expect("the scenario runs");
    assert_eq!(run.header().output.frame_count, EXPECTED_FRAME_COUNT);

    let mut times_s = Vec::new();
    loop {
        if let Some(saved) = run.take_frame() {
            times_s.push(saved.t_s);
        }
        if !run.take_step() {
            break;
        }
    }

    let expected_s: Vec<f64> = (0..EXPECTED_FRAME_COUNT)
        .map(|k| (k * EVERY_N_STEPS) as f64 * DT_S)
        .collect();
    assert_eq!(times_s, expected_s);
    assert_eq!(run.frames_taken(), EXPECTED_FRAME_COUNT);
    assert_eq!(run.steps_taken(), TOTAL_STEPS);
    assert!(run.is_finished());
}

/// A frame is handed out once per saved step, however often it is asked for.
///
/// A chunked driver asks at every chunk boundary, so a chunk that lands on a
/// saved step must not save it twice — a run of four frames with a fifth
/// wedged into it is not the run its header promises.
#[test]
fn the_frame_of_a_step_is_handed_out_once() {
    let scenario = scenario();
    let mut run = RunLoop::of_scenario(&scenario, "once").expect("the scenario runs");

    assert!(run.take_frame().is_some(), "step 0 is a saved step");
    assert!(
        run.take_frame().is_none(),
        "asking again at the same step gives nothing"
    );
    assert_eq!(run.frames_taken(), 1);
}

/// Stepping in chunks is stepping: the frames of a run yielded between every
/// step are bit-for-bit the frames of one taken straight through.
///
/// Byte equality rather than a tolerance, because there is nothing here to be
/// tolerant of: the chunk boundary changes when the caller is handed control,
/// not what the solver does with the state, and identical scenario in means
/// byte-identical output out (CODING_STANDARDS.md § *Correctness and
/// failure*). The chunk sizes straddle the output cadence of 4 — one below it,
/// one that divides it, one that does not, and one longer than the whole run —
/// so a boundary lands on a saved step, between two, and never at all.
#[test]
fn a_run_stepped_in_chunks_is_the_run_stepped_straight_through() {
    let straight_through = frames_of_run_in_chunks(TOTAL_STEPS + 1);
    assert_eq!(straight_through.len(), EXPECTED_FRAME_COUNT as usize);

    for chunk in [1, 3, 4, 5] {
        assert_eq!(
            frames_of_run_in_chunks(chunk),
            straight_through,
            "a run yielding every {chunk} steps is not the run stepped straight through"
        );
    }
}

/// The frames a browser holds are the frames the `run` command writes.
///
/// This is ADR-0012's reproducibility claim stated as a test: the web computes
/// a run instead of downloading one, and what makes those the same run is that
/// both are this loop. Compared as the bytes of `frames.bin`, which is the
/// archival contract of ADR-0004 — the strongest statement available, and the
/// one that would catch a browser path that quietly reordered a sum.
#[cfg(feature = "fs")]
#[test]
fn a_run_stepped_in_chunks_is_the_run_the_cli_writes() {
    let directory = std::env::temp_dir().join(format!(
        "termocline-run-loop-{}-{}",
        std::process::id(),
        line!()
    ));
    let scenario = scenario();
    engine::run_scenario(&scenario, "chunked", &directory).expect("the scenario runs");
    let written = std::fs::read(directory.join(engine::FRAME_FILE_NAME)).expect("frames.bin");
    std::fs::remove_dir_all(&directory).expect("the temporary run is removed");

    let stepped: Vec<u8> = frames_of_run_in_chunks(3).concat();
    assert_eq!(
        stepped, written,
        "the frames stepped in chunks are not the frames the run command wrote"
    );
}
