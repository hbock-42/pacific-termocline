//! T-13.2 acceptance criteria: the teaching panel says *where* it is, *how
//! deep*, *when*, and it says all three in words a visitor already has.
//!
//! T-13.1 drew the ocean; this is what makes the picture readable by someone
//! who has never met the word thermocline. Four things are asserted here, one
//! per acceptance criterion: longitudes read as °E/°W and are derived from the
//! basin the *header* declares (the basin is a scenario parameter, so a
//! hardcoded 120°E–80°W would be a lie the moment a scenario moves it); both
//! coasts are named; depths are metres below the sea surface; and time is days
//! and months taken from the run's own output cadence rather than from a frame
//! index.
//!
//! Like every other view in this crate none of it knows what a GPU is: the
//! ticks, the coast names, the depth labels and the time phrases are values
//! ([ADR-0006]), so they are asserted here rather than looked at.
//!
//! # Where the expected values come from
//!
//! The geography is `CONTEXT.md`, *Basin* — the model domain is roughly
//! 120°E–80°W — and issue #131, which names the two coasts that basin ends at:
//! Indonesia / New Guinea in the west, South America in the east. The tick
//! longitudes and their positions are arithmetic on those bounds, done here
//! rather than read back out of the crate.
//!
//! The times are `engine/scenarios/steady-trades.toml`'s cadence, stated in
//! `tests/common`: 731 daily frames, so the run spans 730 days. Issue #131
//! calls that "2 years", and 730 / 365.2422 = 1.9987 is 2 years to the
//! precision a caption states. A month is a twelfth of the same tropical year,
//! 30.4368 days, so day 396 is 13.01 months in.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod common;

use common::{
    encoded_frames_with_h, header_on, steady_trades_header, EASTERN_WALL_H_M, FRAME_INTERVAL_S, NX,
    NY, PACIFIC, STEADY_TRADES_FRAMES, STEADY_TRADES_PARAMS, WESTERN_WALL_H_M,
};
use termocline_format::{BasinExtent, GridSpec, OutputTiming, RunHeader};
use visualizer::{
    BrowserScenario, CrossSection, LoadedRun, RunBytes, RunClock, SideView, Wall, PLAIN_WORDS,
    TRADE_WINDS,
};

/// How far a tick's computed position may sit from the fraction of the basin
/// its longitude is at. Dimensionless, and round-off on quantities of order
/// one.
const POSITION_TOLERANCE: f64 = 1e-12;

/// How far a stated day may sit from the day the cadence gives, in days.
/// Round-off on a division of two exactly-representable numbers.
const DAY_TOLERANCE: f64 = 1e-9;

/// The mean thermocline depth `H` of `steady-trades.toml`, in metres.
const MEAN_DEPTH_M: f64 = STEADY_TRADES_PARAMS.mean_depth_m;

/// A side view of a one-frame run over `extent`, its `h` uniformly the
/// measured western-wall anomaly so the panel has a depth to it.
fn side_view_over(extent: BasinExtent) -> SideView {
    let grid = GridSpec::new(NX, NY, extent).expect("a basin of the scenario's shape");
    let header = header_on(grid, "steady-trades", 1);
    let field = vec![WESTERN_WALL_H_M; NX * NY];
    let run = run_of(&header, |_| field.clone());
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    SideView::of_section(&section, MEAN_DEPTH_M)
}

/// A run of `header`'s shape whose frames carry the fields `h_m(index)` gives.
fn run_of(header: &RunHeader, h_m: impl Fn(u64) -> Vec<f64>) -> LoadedRun {
    let bytes = RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames_with_h(header, header.output.frame_count, h_m),
    };
    LoadedRun::from_bytes("run-steady-trades", bytes).expect("the run loads")
}

/// The side view of the scenario basin, 120°E–80°W.
fn pacific_side_view() -> SideView {
    side_view_over(PACIFIC)
}

#[test]
fn longitudes_read_as_degrees_east_and_west() {
    // The first acceptance criterion. The axis under the ocean is longitude,
    // written the way an atlas writes it — 140°E, 180°, 120°W — and never as a
    // cell index or a signed degrees-east number with a minus sign in it.
    let view = pacific_side_view();
    let labels: Vec<String> = view
        .longitude_ticks()
        .iter()
        .map(|tick| tick.label())
        .collect();

    // 160° of basin from 120°E: a tick every 20° puts nine of them on the
    // axis, and the antimeridian is 180° rather than either hemisphere's.
    assert_eq!(
        labels,
        ["120°E", "140°E", "160°E", "180°", "160°W", "140°W", "120°W", "100°W", "80°W",],
        "the axis reads as an atlas writes it"
    );
    assert!(
        labels.iter().all(|label| !label.contains('-')),
        "a hemisphere, never a minus sign"
    );
}

#[test]
fn the_longitude_axis_is_derived_from_the_basin_the_header_declares() {
    // The rest of the first criterion, and the reason it is worth a test: the
    // basin is a scenario parameter (`CONTEXT.md`, *Basin*), so an axis that
    // named 120°E–80°W whatever the header said would be wrong for every
    // scenario that truncates the basin differently.
    let moved = BasinExtent::new(150.0, -120.0, -25.0, 25.0);
    let view = side_view_over(moved);
    let labels: Vec<String> = view
        .longitude_ticks()
        .iter()
        .map(|tick| tick.label())
        .collect();

    // 90° of basin from 150°E, so the same 20° step gives five ticks and none
    // of them is anywhere near the scenario basin's western wall.
    assert_eq!(
        labels,
        ["160°E", "180°", "160°W", "140°W", "120°W"],
        "the axis follows the basin the header declares"
    );

    // And a tick sits where its longitude sits across the panel: the basin
    // spans 90° eastward from 150°E, so 180° is a third of the way across.
    let ticks = view.longitude_ticks();
    let antimeridian = ticks
        .iter()
        .find(|tick| tick.label() == "180°")
        .expect("the basin crosses the antimeridian");
    assert!(
        (antimeridian.x_fraction() - 30.0 / 90.0).abs() < POSITION_TOLERANCE,
        "180° sits a third of the way east across a basin from 150°E to 120°W"
    );
    assert!(
        ticks
            .windows(2)
            .all(|pair| pair[1].x_fraction() > pair[0].x_fraction()),
        "the ticks run west to east, as the ocean above them does"
    );
    assert!(ticks
        .iter()
        .all(|tick| (0.0..=1.0).contains(&tick.x_fraction())));
}

#[test]
fn both_coasts_are_named() {
    // The second criterion, in its own words. The basin is closed on all four
    // boundaries (`CONTEXT.md`, *Basin*) and those boundaries are land: the
    // maritime continent in the west and South America in the east. A reader
    // who is told which coast is which knows which way round the picture is.
    let view = pacific_side_view();
    let [west, east] = view.coasts();

    assert_eq!(west.wall(), Wall::West);
    assert_eq!(east.wall(), Wall::East);
    assert_eq!(west.land(), Some("Indonesia / New Guinea"));
    assert_eq!(east.land(), Some("South America"));
    // And each is labelled with the longitude the header puts that wall at, so
    // the name and the axis under it cannot disagree.
    assert!(
        west.label().contains("Indonesia / New Guinea") && west.label().contains("120°E"),
        "the western coast is named and placed: {}",
        west.label()
    );
    assert!(
        east.label().contains("South America") && east.label().contains("80°W"),
        "the eastern coast is named and placed: {}",
        east.label()
    );
}

#[test]
fn a_wall_in_open_ocean_is_not_given_a_coast_it_does_not_touch() {
    // The coasts are derived, like the longitudes: a scenario that truncates
    // the basin mid-Pacific has walls that are not the shore of anything, and
    // naming South America there would be a claim about the world rather than
    // about the run. Such a wall is labelled by where it is and nothing more.
    let mid_ocean = BasinExtent::new(170.0, -150.0, -25.0, 25.0);
    let [west, east] = side_view_over(mid_ocean).coasts();

    assert_eq!(west.land(), None);
    assert_eq!(east.land(), None);
    assert!(
        west.label().contains("170°E") && east.label().contains("150°W"),
        "an unnamed wall still says where it is: {} / {}",
        west.label(),
        east.label()
    );
    for label in [west.label(), east.label()] {
        assert!(
            !label.contains("South America") && !label.contains("Indonesia"),
            "no coast is claimed where the basin ends in open water: {label}"
        );
    }
}

#[test]
fn depths_are_stated_in_metres_below_the_sea_surface() {
    // The third criterion's first half. The vertical axis is depth, labelled
    // in metres, counted down from the sea surface — the depth the interface
    // is drawn at is `H + h` (`CONTEXT.md`), and a reader has to be able to
    // read it off the side of the panel.
    let view = pacific_side_view();
    let ticks = view.depth_ticks();

    assert!(ticks.len() >= 3, "an axis needs ticks to be an axis");
    assert_eq!(
        ticks[0].depth_m(),
        0.0,
        "the axis starts at the sea surface"
    );
    assert_eq!(ticks[0].label(), "0 m");
    assert!(
        ticks
            .windows(2)
            .all(|pair| pair[1].depth_m() > pair[0].depth_m()),
        "the ticks run downward from the surface"
    );
    for tick in &ticks {
        assert!(
            tick.label().ends_with(" m"),
            "every depth carries its unit: {}",
            tick.label()
        );
        assert!(
            tick.depth_m() <= view.deepest_drawn_depth_m(),
            "no tick is drawn past the foot of the panel"
        );
        // A tick sits at its own share of the panel's depth, measured down from
        // the surface: the same placement the interface itself gets.
        let expected = tick.depth_m() / view.deepest_drawn_depth_m();
        assert!((tick.fraction() - expected).abs() < POSITION_TOLERANCE);
    }
    // The scenario's mean thermocline depth is 150 m and the panel reaches
    // 282 m — `H` plus the run's anomaly range, half again — so an axis in
    // round hundreds of metres passes through the depth the thermocline
    // actually sits at.
    assert!(
        ticks.iter().any(|tick| tick.label() == "150 m"),
        "the axis is labelled through the depths the ocean is drawn at"
    );
    assert!(
        view.depth_axis_note().contains("below the sea surface"),
        "the axis says what its metres are measured from: {}",
        view.depth_axis_note()
    );
}

#[test]
fn time_is_days_and_months_taken_from_the_runs_own_cadence() {
    // The third criterion's second half: a visitor is told they are on day 396
    // of 730, thirteen months into a two-year run — not on frame 396 of 731,
    // which is a fact about the file and not about the ocean.
    let clock = RunClock::of_run(OutputTiming {
        frame_count: STEADY_TRADES_FRAMES,
        interval_s: FRAME_INTERVAL_S,
    });

    // 731 daily frames span 730 days: the span of `n` frames is `n - 1`
    // intervals.
    assert!((clock.last_day() - 730.0).abs() < DAY_TOLERANCE);
    assert!((clock.day_of_frame(396) - 396.0).abs() < DAY_TOLERANCE);

    let span = clock.span_phrase();
    assert!(
        span.contains("730 days") && span.contains("2 years"),
        "730 days is two years, and the run says both: {span}"
    );

    let moment = clock.moment(396);
    assert!(
        moment.contains("Day 396 of 730"),
        "the day, out of the run's own length: {moment}"
    );
    // 396 days is 13.01 months of a tropical year's twelfth, 30.4368 days.
    assert!(
        moment.contains("13 months"),
        "and the months a reader counts in: {moment}"
    );
    assert!(
        !moment.to_lowercase().contains("frame"),
        "a frame index is not a time: {moment}"
    );
}

#[test]
fn a_run_written_at_another_cadence_is_dated_by_that_cadence() {
    // "From the run's own cadence" is the whole of the criterion: the days are
    // computed from the interval the header declares, so a run written twice a
    // day is not silently dated as though its frames were days apart.
    let twice_daily = RunClock::of_run(OutputTiming {
        frame_count: 61,
        interval_s: FRAME_INTERVAL_S / 2.0,
    });

    assert!((twice_daily.last_day() - 30.0).abs() < DAY_TOLERANCE);
    assert!((twice_daily.day_of_frame(20) - 10.0).abs() < DAY_TOLERANCE);
    assert!(
        twice_daily.moment(20).contains("Day 10 of 30"),
        "{}",
        twice_daily.moment(20)
    );
    // A month of run is a month, however many frames it took to write it.
    let span = twice_daily.span_phrase();
    assert!(
        span.contains("30 days") && !span.contains("years"),
        "a thirty-day run is thirty days: {span}"
    );
}

#[test]
fn the_teaching_panel_glosses_every_word_it_uses_that_a_visitor_may_not_have() {
    // The fourth criterion. "Thermocline", "anomaly" and "trade winds" are the
    // three words the ticket names; wherever the teaching panel's own prose
    // uses one, the plain-language meaning is in the same sentence, so nothing
    // on the panel sends a reader elsewhere to find out what it means.
    let view = pacific_side_view();
    let mut prose = vec![view.caption(), view.depth_axis_note()];
    prose.extend(view.tilt_note());
    prose.extend(view.coasts().iter().map(|coast| coast.label()));

    for text in &prose {
        for term in PLAIN_WORDS {
            let lowered = text.to_lowercase();
            if lowered.contains(term.term()) {
                assert!(
                    lowered.contains(&term.plain().to_lowercase()),
                    "\"{}\" is used without its plain meaning in: {text}",
                    term.term()
                );
            }
        }
    }
    // And the panel does say what the thing in the picture is: the view exists
    // to show the thermocline, so the word and its meaning are both on it.
    let caption = view.caption().to_lowercase();
    assert!(caption.contains("thermocline"), "{caption}");
}

/// A side view of the measured equilibrium tilt: `h` falling linearly from the
/// western wall's +38.2 m to the eastern wall's −28.2 m (T-07.4, restated in
/// `tests/common`), uniform in latitude.
fn tilted_side_view() -> SideView {
    let mut field = Vec::with_capacity(NX * NY);
    for _ in 0..NY {
        for i in 0..NX {
            #[allow(clippy::cast_precision_loss)]
            let fraction = (i as f64 + 0.5) / NX as f64;
            field.push(WESTERN_WALL_H_M + (EASTERN_WALL_H_M - WESTERN_WALL_H_M) * fraction);
        }
    }
    let header = steady_trades_header(1);
    let run = run_of(&header, |_| field.clone());
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    SideView::of_section(&section, MEAN_DEPTH_M)
}

#[test]
fn the_wind_is_named_in_plain_words_and_only_where_the_ocean_shows_it() {
    // "Trade winds" replaces the project's own "alizés" (`CONTEXT.md`), and
    // like every word on the panel it arrives with its meaning attached. What
    // the sentence *claims* is measured off the frame it is drawn beside: the
    // equilibrium tilt of the control run is 66.4 m between the walls, which
    // is the warm water the trades pile in the west.
    let tilted = tilted_side_view().tilt_note().expect("both ends are drawn");
    assert!(tilted.contains("trade winds"), "{tilted}");
    assert!(
        tilted.contains(TRADE_WINDS.plain()),
        "the wind is explained where it is named: {tilted}"
    );
    // 66.4 m across the walls, a little less between the outermost cell
    // centres, and the note states it to the metre.
    assert!(
        tilted.contains("66 m deeper") || tilted.contains("65 m deeper"),
        "the tilt is the one the run has: {tilted}"
    );

    // And an ocean with no tilt in it is not told it has a warm pool piled in
    // the west: the field here is uniform, so the boundary is level.
    let level = pacific_side_view()
        .tilt_note()
        .expect("both ends are drawn");
    assert!(
        level.contains("level from coast to coast"),
        "a level ocean is described as level: {level}"
    );
    assert!(
        !level.contains("piled up"),
        "nothing is claimed that the picture does not show: {level}"
    );
}

#[test]
fn the_plain_words_are_themselves_plain() {
    // A gloss written in the vocabulary it is glossing explains nothing. The
    // three the ticket names are the three the key defines, and no definition
    // leans on another.
    let terms: Vec<&str> = PLAIN_WORDS.iter().map(|term| term.term()).collect();
    assert_eq!(terms, ["thermocline", "anomaly", "trade winds"]);

    for term in PLAIN_WORDS {
        let plain = term.plain().to_lowercase();
        for other in PLAIN_WORDS {
            assert!(
                !plain.contains(other.term()),
                "\"{}\" is explained with the word \"{}\"",
                term.term(),
                other.term()
            );
        }
        assert!(
            term.glossed().contains(term.term()) && term.glossed().contains(term.plain()),
            "a glossed term carries both halves: {}",
            term.glossed()
        );
    }
}

#[test]
fn the_scenarios_are_offered_in_words_the_key_explains() {
    // The scenario picker is on screen beside the teaching panel, and it named
    // the trade winds in French ("alizés", `CONTEXT.md`) — a term the project's
    // own glossary defines and a visitor does not. Nothing the picker says now
    // uses a word the plain-language key does not cover.
    for scenario in BrowserScenario::ALL {
        let text = format!("{} — {}", scenario.name, scenario.summary).to_lowercase();
        assert!(
            !text.contains("aliz"),
            "the picker speaks the reader's language: {text}"
        );
    }
}

#[test]
fn the_side_view_of_the_scenario_basin_names_its_own_geography() {
    // The four criteria together, on the run a visitor actually opens: the
    // control scenario's own header, drawn through the same path the shell
    // draws it through.
    let header = steady_trades_header(STEADY_TRADES_FRAMES);
    let field = vec![WESTERN_WALL_H_M; NX * NY];
    let run = run_of(&header, |_| field.clone());
    let frame = run.frame(0).expect("a run of frames has a frame 0");
    let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    let view = SideView::of_section(&section, run.header().physical_params.mean_depth_m);

    assert_eq!(view.coasts()[0].land(), Some("Indonesia / New Guinea"));
    assert_eq!(view.coasts()[1].land(), Some("South America"));
    assert_eq!(
        view.longitude_ticks().first().map(|t| t.label()),
        Some("120°E".to_owned())
    );
    assert!(!view.depth_ticks().is_empty());

    let clock = RunClock::of_run(run.header().output);
    assert!(clock.span_phrase().contains("2 years"));
}
