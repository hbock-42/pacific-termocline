//! Acceptance tests for T-08.6 — the browser computes the run it draws.
//!
//! ADR-0012's claim is that this costs the views nothing: a run computed in
//! the tab is a [`LoadedRun`] like any other, so the heatmap, the wind
//! overlay, the cross-section, the point time series and the comparison read
//! it without being told where it came from. Most of this file is that claim,
//! checked view by view.
//!
//! Nothing here is measured out of a run. The scenario values are read off the
//! TOML text by eye; the state after one step is the closed-form first step of
//! the momentum equation of ADR-0003, as `engine/tests/filesystem_free_api.rs`
//! states it; and the sizes are the arithmetic of the format — one `f64` per
//! point of each field — held against the 19.9 MB the browser scenario
//! actually writes natively.

use termocline_format::{RunHeader, Variable};
use visualizer::{
    BrowserScenario, Comparison, ComputeError, ComputedRun, CrossSection, FrameBudget, Heatmap,
    LoadedRun, PointSeries, WindOverlay,
};

/// A run small enough to compute inside a test: a 40° × 10° box at 1° under a
/// steady easterly, 12 steps of an hour, a frame every 4.
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
const EVERY_N_STEPS: u64 = 4;
const EXPECTED_NX: usize = 40;
const EXPECTED_NY: usize = 10;
/// Frames at steps 0, 4, 8 and 12.
const EXPECTED_FRAME_COUNT: u64 = 4;
/// τ₀ = −0.05 Pa, so the trades are easterly.
const EQUATORIAL_ZONAL_STRESS_PA: f64 = -0.05;
/// Ly, the meridional decay scale of the trade winds, in metres.
const MERIDIONAL_DECAY_SCALE_M: f64 = 361_000.0;

/// Steps enough to reach the end of `SCENARIO_TOML` however they are chunked.
const ENOUGH_STEPS: u64 = 100;

/// A run computed to the end of its schedule, in chunks of `chunk` steps.
fn computed_in_chunks(chunk: u64) -> ComputedRun {
    let mut run = ComputedRun::start(SCENARIO_TOML, "test scenario", FrameBudget::browser())
        .expect("the scenario is inside the browser's budget");
    while !run.is_finished() {
        run.advance_steps(chunk).expect("the run computes");
    }
    run
}

/// A run computed to the end of its schedule.
fn computed_run() -> ComputedRun {
    computed_in_chunks(ENOUGH_STEPS)
}

/// The frames a computed run holds come from the engine, and there are as many
/// as its header promises.
#[test]
fn the_engine_produces_the_frames_and_the_header_counts_them() {
    let computed = computed_run();
    let run = computed.run();

    assert_eq!(run.header().output.frame_count, EXPECTED_FRAME_COUNT);
    assert_eq!(run.frame_count(), EXPECTED_FRAME_COUNT);
    assert!(run.is_complete());
    assert_eq!(
        computed.progress(),
        (EXPECTED_FRAME_COUNT, EXPECTED_FRAME_COUNT)
    );

    // Frame `k` is the state after `k · N` steps, at `k · N · dt` seconds —
    // the output schedule's own arithmetic, not anything read out of a run.
    for k in 0..EXPECTED_FRAME_COUNT {
        let frame = run.frame(k).expect("the run holds every frame it counts");
        #[allow(clippy::cast_precision_loss)]
        let expected_s = (k * EVERY_N_STEPS) as f64 * DT_S;
        assert_eq!(frame.t_s(), expected_s);
        assert_eq!(frame.h().len(), EXPECTED_NX * EXPECTED_NY);
    }
}

/// The wind a computed frame records is the wind the scenario prescribes.
///
/// The strongest easterly anywhere in the frame is the one on the row nearest
/// the equator, and the trade winds are `τx(y) = τ₀·exp(−(y/Ly)²)`
/// (CONTEXT.md, *Wind stress*; `engine/src/forcing.rs`). A 1° grid over 10° of
/// latitude puts its nearest `u` row half a cell off the equator, so the
/// expected value is that formula evaluated there — arithmetic, not a number
/// read out of a run.
///
/// The check is that the browser's frames carry the forcing that drove them,
/// which is what the wind overlay draws.
#[test]
fn a_computed_frame_carries_the_stress_that_drove_it() {
    let scenario = engine::Scenario::from_toml(SCENARIO_TOML).expect("the scenario text is valid");
    let basin = scenario.basin();
    let scaled = basin.y_of_row_m(engine::U_STAGGERING, EXPECTED_NY / 2) / MERIDIONAL_DECAY_SCALE_M;
    let expected_pa = EQUATORIAL_ZONAL_STRESS_PA * (-scaled * scaled).exp();

    let computed = computed_run();
    let frame = computed.run().frame(0).expect("the first frame");
    let strongest = frame
        .tau_x()
        .iter()
        .copied()
        .fold(0.0_f64, |strongest, tau| strongest.min(tau));
    // Exactly: the frame records the very field the step read, so the only
    // arithmetic between the formula and the assertion is the one the engine
    // did (ADR-0009).
    assert_eq!(
        strongest, expected_pa,
        "the strongest easterly in the frame is {strongest} Pa, and the trade-wind profile puts \
         {expected_pa} Pa on the row nearest the equator"
    );
}

/// Chunking is a choice about when the tab gets to draw, not about what the
/// engine computes.
///
/// Byte equality of every frame, because there is nothing to be tolerant of:
/// identical scenario in, byte-identical output out (CODING_STANDARDS.md §
/// *Correctness and failure*). The chunk sizes straddle the output cadence of
/// 4, so a boundary lands on a saved step, between two, and never at all.
#[test]
fn a_run_computed_in_chunks_is_the_run_computed_in_one() {
    let whole = frames_of(computed_in_chunks(ENOUGH_STEPS).run());
    assert_eq!(whole.len(), EXPECTED_FRAME_COUNT as usize);

    for chunk in [1, 3, 4, 5] {
        assert_eq!(
            frames_of(computed_in_chunks(chunk).run()),
            whole,
            "a run advanced {chunk} steps at a time is not the run advanced in one"
        );
    }
}

/// Every frame of `run`, as the values a view would read.
fn frames_of(run: &LoadedRun) -> Vec<(f64, Vec<f64>, Vec<f64>)> {
    (0..run.frame_count())
        .map(|index| {
            let frame = run.frame(index).expect("a frame the run counts");
            (frame.t_s(), frame.h().to_vec(), frame.tau_x().to_vec())
        })
        .collect()
}

/// A partly computed run is a run: it holds the frames produced so far, and
/// nothing beyond them.
///
/// This is what makes progress visible rather than merely reported — the
/// scrubber, the map and the charts are over the frames that exist, and the
/// count they are bounded by is the run's own.
#[test]
fn a_run_being_computed_holds_only_the_frames_it_has_produced() {
    let mut computed = ComputedRun::start(SCENARIO_TOML, "partial", FrameBudget::browser())
        .expect("the scenario is inside the browser's budget");
    assert_eq!(computed.run().frame_count(), 0);
    assert!(!computed.run().is_complete());

    // Four steps is the cadence, so this is frames 0 and 1 and no more.
    computed
        .advance_steps(EVERY_N_STEPS)
        .expect("the run computes");
    assert_eq!(computed.run().frame_count(), 2);
    assert_eq!(computed.progress(), (2, EXPECTED_FRAME_COUNT));
    assert!(computed.run().frame(2).is_none());
    assert!(!computed.run().is_complete());
    assert!(!computed.is_finished());
}

/// The colour scale widens as frames arrive and never narrows, and the scale
/// of the finished run is the run-wide scale T-08.2 asks for.
///
/// The tilt grows through the run — the alizés pile water in the west from
/// rest — so a later frame is louder than an earlier one and the scale has to
/// follow it. What is asserted is the two properties a provisional scale must
/// have: monotonic while the run develops, and equal at the end to the scale
/// of every frame taken together.
#[test]
fn the_scale_widens_as_frames_arrive_and_ends_run_wide() {
    let mut computed = ComputedRun::start(SCENARIO_TOML, "widening", FrameBudget::browser())
        .expect("the scenario is inside the browser's budget");
    let mut widths_m = Vec::new();
    while !computed.is_finished() {
        computed
            .advance_steps(EVERY_N_STEPS)
            .expect("the run computes");
        widths_m.push(computed.run().anomaly_scale().half_range_m());
    }

    assert!(
        widths_m.windows(2).all(|pair| pair[1] >= pair[0]),
        "the scale narrowed as the run developed: {widths_m:?}"
    );
    let run = computed.run();
    let loudest_m = (0..run.frame_count())
        .flat_map(|index| run.frame(index).expect("a frame").h().to_vec())
        .fold(0.0_f64, |loudest, h_m| loudest.max(h_m.abs()));
    assert_eq!(
        run.anomaly_scale().half_range_m(),
        loudest_m,
        "the finished run's scale is not the largest anomaly anywhere in it"
    );
    // A run that has grown past the first frame has something to show, or the
    // widening is not being driven by the frames.
    assert!(loudest_m > 0.0);
}

/// Every view of the shell draws a computed run, unchanged.
///
/// The cheap path of ADR-0012 stated as a test: each of these takes a
/// `&LoadedRun` or a frame of one, and none of them was touched by this
/// ticket. If a computed run were not a run, this is where it would show.
#[test]
fn every_view_draws_a_computed_run() {
    let computed = computed_run();
    let run = computed.run();
    let grid = run.header().grid;
    let frame = run.frame(run.frame_count() - 1).expect("the last frame");

    Heatmap::of_frame(grid, &frame, run.anomaly_scale()).expect("the heatmap draws it");
    WindOverlay::of_frame(grid, &frame, run.wind_stress_scale()).expect("the overlay draws it");
    CrossSection::of_frame(grid, &frame, run.anomaly_scale()).expect("the section draws it");

    let point = visualizer::BasinPoint::at_map_fraction(grid, 0.5, 0.5)
        .expect("mid-basin is a cell of the grid");
    let series = PointSeries::at_point(run, point).expect("the time series walks it");
    assert_eq!(series.samples().len(), EXPECTED_FRAME_COUNT as usize);

    // And the comparison, which is the one view over two runs: two runs of the
    // same scenario are comparable, and share the frame index the panels are
    // drawn from.
    let other = computed_run();
    let comparison = Comparison::of(run, other.run()).expect("two runs of one scenario compare");
    assert_eq!(comparison.frame_count(), EXPECTED_FRAME_COUNT);

    // The metadata panel says what the run is, off the header the engine wrote
    // — the same rows a run read from disk fills in.
    let metadata = run.metadata();
    assert!(metadata
        .iter()
        .any(|row| row.label == "Grid"
            && row.value == format!("{EXPECTED_NX} × {EXPECTED_NY} cells")));
}

/// The frame budget is enforced before a step is taken, and the browser's
/// shipped scenarios are inside it.
///
/// The control run is the one ADR-0012 names as impossible: 731 frames of a
/// 320 × 100 basin. The refusal has to name the size, because a visitor who
/// cannot see why it was refused cannot pick a scenario that would work.
#[test]
fn a_run_too_big_for_a_tab_is_refused_before_it_starts() {
    let control = SCENARIO_TOML
        .replace("resolution_deg = 1.0", "resolution_deg = 0.5")
        .replace(
            "eastern_longitude_deg = 160.0",
            "eastern_longitude_deg = -80.0",
        )
        .replace(
            "southern_latitude_deg = -5.0",
            "southern_latitude_deg = -25.0",
        )
        .replace(
            "northern_latitude_deg = 5.0",
            "northern_latitude_deg = 25.0",
        )
        .replace("total_steps = 12", "total_steps = 17520")
        .replace("output_every_n_steps = 4", "output_every_n_steps = 24");

    let refused = ComputedRun::start(&control, "control", FrameBudget::browser())
        .err()
        .expect("the control run does not fit in a tab");
    let ComputeError::Budget(exceeded) = refused else {
        panic!("the control run was refused for something other than its size: {refused}");
    };
    assert_eq!(exceeded.frame_count, 731);
    assert!(
        exceeded.needed_bytes > exceeded.budget_bytes,
        "a run inside the budget was refused for exceeding it"
    );
    // ADR-0012 states the control run at 941 MB, to the three figures it is
    // written with. The estimate counts field values and not the format's
    // per-field length prefixes, so it lands just under; a thousandth is the
    // width of the last figure ADR-0012 gives, which is as close as two
    // numbers written that way can be held to each other.
    let stated_bytes = 941_000_000.0;
    #[allow(clippy::cast_precision_loss)]
    let estimated = exceeded.needed_bytes as f64;
    assert!(
        (estimated - stated_bytes).abs() / stated_bytes < 0.001,
        "the control run is estimated at {estimated} bytes against the 941 MB ADR-0012 states"
    );

    let message = exceeded.to_string();
    assert!(
        message.contains("940.6 MB") && message.contains("33.6 MB"),
        "the refusal does not say what the run costs and what it is allowed: {message}"
    );
}

/// Every scenario the browser ships is one it can actually compute.
///
/// A scenario that did not parse, or that did not fit the budget, would be a
/// button that fails when pressed — and the budget is the whole reason these
/// files are coarser than the engine's.
#[test]
fn every_shipped_scenario_fits_the_browser_budget() {
    for scenario in BrowserScenario::ALL {
        let computed = ComputedRun::start(scenario.toml, scenario.name, FrameBudget::browser())
            .unwrap_or_else(|error| panic!("{} does not start: {error}", scenario.name));
        let header = computed.run().header();
        assert_eq!(
            (header.grid.nx(), header.grid.ny()),
            (80, 25),
            "{} is not on the browser grid",
            scenario.name
        );
        assert!(
            FrameBudget::bytes_of(header) <= FrameBudget::browser().max_bytes(),
            "{} does not fit the browser's frame budget",
            scenario.name
        );
    }
}

/// The size a run is admitted on accounts for the size it turns out to be.
///
/// `termocline run` writes 19 935 776 bytes of `frames.bin` for this scenario
/// — an independent artifact, produced by the CLI rather than by anything
/// under test here. What the budget estimates is the field values alone, and
/// the difference between the two is not a fudge factor but the rest of what
/// the format puts in a frame, which is enumerable — so the two are asserted
/// to add up exactly rather than to agree to within something.
#[test]
fn the_estimated_size_of_a_run_accounts_for_the_size_it_writes() {
    let computed = ComputedRun::start(
        BrowserScenario::default_scenario().toml,
        "sizing",
        FrameBudget::browser(),
    )
    .expect("the browser scenario fits");
    let header = computed.run().header();
    let estimated = FrameBudget::bytes_of(header);

    /// Bytes an encoded frame carries beyond its field values: the frame's
    /// model time (one `f64`), `bincode`'s tag for the absent `T'` of an
    /// uncoupled run (one byte), and a length prefix per field — three bytes
    /// each at these field lengths, which is what `bincode`'s standard
    /// configuration spends on a length of 251..=65 535.
    const FORMAT_BYTES_PER_FRAME: u64 = 8 + 1 + 5 * 3;
    /// `frames.bin` of `scenarios/browser-steady-trades.toml`, as written by
    /// `termocline run` on the native engine.
    const MEASURED_BYTES: u64 = 19_935_776;

    assert_eq!(
        estimated + FORMAT_BYTES_PER_FRAME * header.output.frame_count,
        MEASURED_BYTES,
        "the field values plus the format's own bytes are not the run the engine writes"
    );
}

/// A header's size is the arithmetic of the format: one `f64` per point of
/// each field it declares, per frame.
///
/// Independent of any run: the staggered field lengths come from the grid, and
/// the frame count from the header.
#[test]
fn a_runs_size_is_one_f64_per_value_it_holds() {
    let computed = computed_run();
    let header: &RunHeader = computed.run().header();
    let per_frame: u64 = Variable::LINEAR_CORE
        .iter()
        .map(|variable| header.grid.field_len(*variable) as u64 * 8)
        .sum();
    assert_eq!(
        FrameBudget::bytes_of(header),
        per_frame * EXPECTED_FRAME_COUNT
    );
}

/// A run computed in the tab and the same run read back from bytes are the
/// same run.
///
/// The claim ADR-0012's cheap path rests on: [`LoadedRun`] has two origins and
/// one behaviour, so a view cannot tell which it was handed. The bytes here
/// are assembled the way a written run's are — the header as JSON, the frames
/// encoded one after another with nothing between them (ADR-0004) — and the
/// run read back has to agree with the computed one frame for frame and scale
/// for scale.
#[test]
fn a_computed_run_and_a_run_read_from_bytes_are_indistinguishable() {
    let computed = computed_run();
    let run = computed.run();

    let mut frames = Vec::new();
    for index in 0..run.frame_count() {
        let frame = run.frame(index).expect("a frame the run counts");
        frames.extend_from_slice(&termocline_format::encode_frame(&frame).expect("it encodes"));
    }
    let bytes = visualizer::RunBytes {
        header: serde_json::to_vec(run.header()).expect("the header serialises"),
        frames,
    };
    let read = LoadedRun::from_bytes("read back", bytes).expect("the bytes are a run");

    assert_eq!(read.frame_count(), run.frame_count());
    assert_eq!(read.anomaly_scale(), run.anomaly_scale());
    assert_eq!(read.wind_stress_scale(), run.wind_stress_scale());
    assert!(read.is_complete() && run.is_complete());
    assert_eq!(frames_of(&read), frames_of(run));
    assert_eq!(read.metadata(), run.metadata());
}
