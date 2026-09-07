//! T-13.4 acceptance criteria: the captions under the ocean are read off the
//! run, they say only what the run supports, and their numbers are the run's.
//!
//! The ticket's one hard rule is that no caption fires on a schedule. "At
//! frame 40, say this" is a lie waiting to happen: the scenario, the grid, the
//! parameters and the length of the run are all a reader's to change (T-13.3
//! made the scenario a button), and a caption keyed to an index survives none
//! of it. So every caption here is a function of an [`EquatorialReading`] — a
//! value holding what was *measured* off the run — and each one is asserted
//! against fabricated readings that should and should not trigger it.
//!
//! The last three tests then take the scenarios a visitor actually presses
//! (T-13.3), run them, and hold the sentences to what each run turns out to
//! be: the control holds a tilt under a wind that does not change, the
//! relaxation lets it go, and the burst leaves ocean in the east that has not
//! moved. Those are the ticket's three register examples, on the runs they are
//! examples of.
//!
//! # Where the expected values come from
//!
//! The states are fabricated here, so the expected caption of each is read off
//! the state by hand rather than out of the crate. The one reading taken from
//! a real run is checked against T-07.4's measured equilibrium profile as
//! `tests/common` restates it — +38.2 m at the western wall and −28.2 m at the
//! eastern, so a 66.4 m difference between the two ends — and against the
//! −0.05 Pa of easterly stress every scenario in `engine/scenarios/` is
//! written with (`docs/planning/01-scientific-model.md`, *Wind stress*).
//!
//! # The claim no caption may make
//!
//! `docs/enso-oscillation-report.md` records the coupled model's measured
//! period as 1.03 years and its own acceptance criterion — a period in the
//! observed 2–7 year band — as **not met**. A caption that spoke of a cycle,
//! of an ENSO period, or of what the ocean is about to do would contradict the
//! project's own validation ledger, so the vocabulary test below sweeps every
//! caption the module can produce for that language.

mod common;

use common::{
    encoded_frames_with_fields, steady_trades_header, FrameFields, EASTERN_WALL_H_M, MID_BASIN_H_M,
    NX, NY, STEADY_TRADES_PARAMS, WESTERN_WALL_H_M,
};
use termocline_format::{RunHeader, Variable};
use visualizer::{
    Caption, CaptionTopic, ComputedRun, Departure, EquatorialReading, FrameBudget, LoadedRun,
    RunBytes, ScenarioPreset,
};

/// The mean thermocline depth `H` of the fixture scenario, in metres.
const MEAN_DEPTH_M: f64 = STEADY_TRADES_PARAMS.mean_depth_m;

/// Half the anomaly range of the fabricated readings below, in metres: the
/// equilibrium profile's own extreme, so the flat-tilt threshold the captions
/// use is the one the real run's scale would give.
const HALF_RANGE_M: f64 = 38.2;

/// The easterly stress every scenario in `engine/scenarios/` is written with,
/// in Pa, negative because the alizés blow east to west
/// (`docs/planning/01-scientific-model.md`, *Wind stress*).
const TRADE_STRESS_PA: f64 = -0.05;

/// How far a measured metre may sit from the metre the fixture wrote.
///
/// The fixture's fields are written and read back through `f64`, and the
/// equatorial reading averages the two rows either side of the equator, which
/// carry the same value: round-off on quantities of order ten.
const METRE_TOLERANCE: f64 = 1e-9;

/// How far a measured day may sit from the day the cadence gives, in days.
/// Round-off on a division of two exactly-representable numbers.
const DAY_TOLERANCE: f64 = 1e-9;

/// How far a measured pascal may sit from the pascal the fixture wrote, for
/// the same reason.
const PASCAL_TOLERANCE: f64 = 1e-12;

/// A reading of an ocean at rest under the full trades: no tilt, no history.
///
/// Every fabricated state below is this one with the fields the caption under
/// test reads changed, so what makes a caption fire is visible in the test
/// that fires it.
fn at_rest() -> EquatorialReading {
    EquatorialReading {
        west_depth_m: Some(MEAN_DEPTH_M),
        east_depth_m: Some(MEAN_DEPTH_M),
        earlier_tilt_m: None,
        lookback_days: 30.0,
        equatorial_stress_pa: TRADE_STRESS_PA,
        strongest_stress_pa: TRADE_STRESS_PA.abs(),
        anomaly_half_range_m: HALF_RANGE_M,
        departure: None,
    }
}

/// The same reading with the thermocline `west_m` down in the west and
/// `east_m` down in the east.
fn tilted(west_m: f64, east_m: f64) -> EquatorialReading {
    EquatorialReading {
        west_depth_m: Some(west_m),
        east_depth_m: Some(east_m),
        ..at_rest()
    }
}

/// The text of the caption on `topic`, or `None` where the reading did not
/// produce one.
fn caption_on(reading: &EquatorialReading, topic: CaptionTopic) -> Option<String> {
    reading
        .captions()
        .into_iter()
        .find(|caption| caption.topic() == topic)
        .map(|caption| caption.text().to_owned())
}

/// Every caption of `reading`, run together, lower-cased.
fn all_text(reading: &EquatorialReading) -> String {
    reading
        .captions()
        .iter()
        .map(Caption::text)
        .collect::<Vec<_>>()
        .join(" ")
        .to_lowercase()
}

#[test]
fn the_wind_caption_follows_the_stress_that_was_measured() {
    // The first criterion, on the caption the ticket's own register opens
    // with. It fires on the sign and the size of the measured equatorial
    // stress and on nothing else: full trades pile water west, a stress that
    // has fallen well below the run's strongest is slack, and a positive
    // stress is a wind blowing the other way.
    let blowing = caption_on(&at_rest(), CaptionTopic::Wind).expect("the trades are blowing");
    assert!(
        blowing.contains("piling warm water into the west"),
        "the trades are described by what they do: {blowing}"
    );
    assert!(
        blowing.contains("0.05 Pa"),
        "and by the stress that was measured: {blowing}"
    );

    // A fifth of the strongest stress in the run: the relaxation preset's own
    // floor, and below anything the seasonal cycle's ±20 % reaches.
    let slack = EquatorialReading {
        equatorial_stress_pa: TRADE_STRESS_PA * 0.2,
        ..at_rest()
    };
    let text = caption_on(&slack, CaptionTopic::Wind).expect("a slack wind is still a wind");
    assert!(
        text.contains("slack") && text.contains("0.01 Pa") && text.contains("0.05 Pa"),
        "the slack wind is stated against the strongest the run reached: {text}"
    );
    assert!(
        !text.contains("piling warm water into the west"),
        "a slack wind is not piling anything: {text}"
    );

    let reversed = EquatorialReading {
        equatorial_stress_pa: 0.04,
        ..at_rest()
    };
    let text = caption_on(&reversed, CaptionTopic::Wind).expect("a reversed wind is a wind");
    assert!(
        text.contains("west to east") && text.contains("0.04 Pa"),
        "a wind blowing the other way says so: {text}"
    );

    let calm = EquatorialReading {
        equatorial_stress_pa: 0.0,
        strongest_stress_pa: 0.0,
        ..at_rest()
    };
    let text = caption_on(&calm, CaptionTopic::Wind).expect("no wind is also a state");
    assert!(text.contains("No wind"), "{text}");
}

#[test]
fn the_tilt_caption_states_the_difference_between_the_two_ends_it_measured() {
    // The ticket's second example, in metres taken off the two ends of the
    // run's own equator. 188 m down in the west against 122 m in the east is
    // the 66 m difference T-07.4's equilibrium holds.
    let reading = tilted(
        MEAN_DEPTH_M + WESTERN_WALL_H_M,
        MEAN_DEPTH_M + EASTERN_WALL_H_M,
    );
    let text = caption_on(&reading, CaptionTopic::Tilt).expect("a tilted ocean has a tilt");
    assert!(
        text.contains("66 m") && text.contains("deeper in the west than in the east"),
        "the measured difference, and which way round it is: {text}"
    );

    // The same ocean turned over: the caption follows the sign rather than
    // assuming the west is always the deep end.
    let flipped = tilted(
        MEAN_DEPTH_M + EASTERN_WALL_H_M,
        MEAN_DEPTH_M + WESTERN_WALL_H_M,
    );
    let text = caption_on(&flipped, CaptionTopic::Tilt).expect("a tilt the other way is a tilt");
    assert!(
        text.contains("66 m") && text.contains("deeper in the east than in the west"),
        "{text}"
    );

    // And a different ocean gives a different number: the metres are the
    // run's, not a constant.
    let smaller = tilted(
        MEAN_DEPTH_M + MID_BASIN_H_M,
        MEAN_DEPTH_M + EASTERN_WALL_H_M,
    );
    let text = caption_on(&smaller, CaptionTopic::Tilt).expect("a tilt of 37 m is a tilt");
    assert!(
        text.contains("37 m"),
        "8.4 m over −28.2 m is 36.6 m: {text}"
    );
}

#[test]
fn a_level_thermocline_is_not_reported_as_a_tilt() {
    // The should-not half of the first criterion. An ocean at rest has two
    // ends at the same depth, and a caption that called that a tilt would be
    // describing a slope nobody can see. The threshold is a tenth of the run's
    // own anomaly range, so it moves with the run rather than being a metre
    // count fixed here.
    let level = caption_on(&at_rest(), CaptionTopic::Tilt).expect("a level ocean is a state");
    assert!(
        level.contains("nearly the same depth"),
        "a level ocean is level: {level}"
    );
    assert!(
        !level.contains("deeper in the west"),
        "and is not a tilt: {level}"
    );

    // A twentieth of the range either way is still level; a fifth is a tilt.
    let barely = tilted(MEAN_DEPTH_M + HALF_RANGE_M * 0.05, MEAN_DEPTH_M);
    assert!(caption_on(&barely, CaptionTopic::Tilt)
        .expect("a state")
        .contains("nearly the same depth"));
    let plainly = tilted(MEAN_DEPTH_M + HALF_RANGE_M * 0.4, MEAN_DEPTH_M);
    assert!(caption_on(&plainly, CaptionTopic::Tilt)
        .expect("a state")
        .contains("deeper in the west"));
}

#[test]
fn the_trend_caption_follows_the_change_that_was_measured() {
    // "Whether it is growing or collapsing" — read as the difference between
    // the tilt now and the tilt the run held a measured number of days ago,
    // never as a place in the run.
    let equilibrium = MEAN_DEPTH_M + WESTERN_WALL_H_M;
    let building = EquatorialReading {
        earlier_tilt_m: Some(20.0),
        ..tilted(equilibrium, MEAN_DEPTH_M + EASTERN_WALL_H_M)
    };
    let text = caption_on(&building, CaptionTopic::TiltTrend).expect("a trend was measured");
    assert!(
        text.contains("growing") && text.contains("46 m") && text.contains("30 days"),
        "66.4 m now against 20 m thirty days ago: {text}"
    );

    let collapsing = EquatorialReading {
        earlier_tilt_m: Some(66.4),
        ..tilted(MEAN_DEPTH_M + 10.0, MEAN_DEPTH_M)
    };
    let text = caption_on(&collapsing, CaptionTopic::TiltTrend).expect("a trend was measured");
    assert!(
        text.contains("shrinking") && text.contains("56 m"),
        "10 m now against 66.4 m thirty days ago: {text}"
    );

    // A tilt that has barely moved is holding, not growing.
    let holding = EquatorialReading {
        earlier_tilt_m: Some(66.0),
        ..tilted(equilibrium, MEAN_DEPTH_M + EASTERN_WALL_H_M)
    };
    let text = caption_on(&holding, CaptionTopic::TiltTrend).expect("a trend was measured");
    assert!(text.contains("holding"), "{text}");

    // And a run with nothing behind the frame on screen gets no trend at all,
    // rather than a trend measured against itself.
    assert_eq!(
        caption_on(&tilted(equilibrium, MEAN_DEPTH_M), CaptionTopic::TiltTrend),
        None,
        "no earlier tilt was measured, so nothing is said about one"
    );
}

#[test]
fn the_warm_water_slides_back_east_only_when_both_halves_of_that_are_measured() {
    // The ticket's third example, and the one compound state in the set: it
    // takes a wind that has gone slack *and* a tilt that is collapsing *and* a
    // west that is still the deep end. Any one of the three missing and the
    // caption is silent, because the sentence would then be describing
    // something the run is not doing.
    let sliding = EquatorialReading {
        equatorial_stress_pa: TRADE_STRESS_PA * 0.2,
        earlier_tilt_m: Some(66.4),
        ..tilted(MEAN_DEPTH_M + 10.0, MEAN_DEPTH_M)
    };
    let text = caption_on(&sliding, CaptionTopic::WarmWaterEast).expect("both halves hold");
    assert!(
        text.contains("sliding back east") && text.contains("slack"),
        "{text}"
    );

    // How the wind let go is measured too: a wind blowing the other way is not
    // a slack one, and the sentence says which was read.
    let reversed = EquatorialReading {
        equatorial_stress_pa: 0.04,
        ..sliding
    };
    let text = caption_on(&reversed, CaptionTopic::WarmWaterEast).expect("both halves hold");
    assert!(
        text.contains("sliding back east") && text.contains("pushing the other way"),
        "{text}"
    );

    let winds_still_blowing = EquatorialReading {
        equatorial_stress_pa: TRADE_STRESS_PA,
        ..sliding
    };
    assert_eq!(
        caption_on(&winds_still_blowing, CaptionTopic::WarmWaterEast),
        None,
        "the trades are at full strength: nothing is sliding anywhere"
    );

    let tilt_still_building = EquatorialReading {
        earlier_tilt_m: Some(5.0),
        ..sliding
    };
    assert_eq!(
        caption_on(&tilt_still_building, CaptionTopic::WarmWaterEast),
        None,
        "the tilt is growing, so the water is not going back"
    );

    let already_gone = EquatorialReading {
        earlier_tilt_m: Some(-66.4),
        ..tilted(MEAN_DEPTH_M, MEAN_DEPTH_M + 10.0)
    };
    assert_eq!(
        caption_on(&already_gone, CaptionTopic::WarmWaterEast),
        None,
        "the west is no longer the deep end: there is no pile left to slide"
    );
}

#[test]
fn the_front_caption_says_how_far_east_the_change_has_reached() {
    // "Where a wave front has reached", measured: the largest change over the
    // lookback window, and the easternmost longitude at which the change is
    // still half that peak. It is silent once the change has reached the
    // eastern wall, because then there is no leading edge left to point at.
    let travelling = EquatorialReading {
        departure: Some(Departure {
            peak_change_m: 20.0,
            half_peak_longitude_deg_east: -170.0,
            reaches_eastern_wall: false,
        }),
        ..at_rest()
    };
    let text = caption_on(&travelling, CaptionTopic::Front).expect("a front was measured");
    assert!(
        text.contains("170°W") && text.contains("30 days"),
        "the longitude the change reached, over the window it was measured in: {text}"
    );

    let arrived = EquatorialReading {
        departure: Some(Departure {
            reaches_eastern_wall: true,
            ..travelling.departure.expect("a departure was set")
        }),
        ..at_rest()
    };
    assert_eq!(
        caption_on(&arrived, CaptionTopic::Front),
        None,
        "the change is everywhere: there is no front to name"
    );

    // A change too small to see is not a front either — the same tenth of the
    // run's range the tilt caption is level below.
    let quiet = EquatorialReading {
        departure: Some(Departure {
            peak_change_m: 0.5,
            ..travelling.departure.expect("a departure was set")
        }),
        ..at_rest()
    };
    assert_eq!(
        caption_on(&quiet, CaptionTopic::Front),
        None,
        "0.5 m is noise"
    );

    assert_eq!(
        caption_on(&at_rest(), CaptionTopic::Front),
        None,
        "nothing was measured, so nothing is claimed"
    );
}

#[test]
fn a_column_the_model_broke_down_in_gets_no_depth_caption() {
    // A caption states a depth the run produced or it states nothing. Where
    // the upper layer has no thickness left the side view draws no water
    // (T-13.1), and a caption is owed the same refusal.
    let broken = EquatorialReading {
        west_depth_m: None,
        earlier_tilt_m: Some(20.0),
        ..tilted(MEAN_DEPTH_M + WESTERN_WALL_H_M, MEAN_DEPTH_M)
    };
    assert_eq!(caption_on(&broken, CaptionTopic::Tilt), None);
    assert_eq!(caption_on(&broken, CaptionTopic::TiltTrend), None);
    assert_eq!(caption_on(&broken, CaptionTopic::WarmWaterEast), None);
    // The wind was still measured, so it is still described.
    assert!(caption_on(&broken, CaptionTopic::Wind).is_some());
}

#[test]
fn no_caption_claims_a_cycle_a_period_or_a_prediction() {
    // The second acceptance criterion, and the one that protects the project's
    // own validation ledger: `docs/enso-oscillation-report.md` records the
    // measured period as 1.03 years and the 2–7 year band as unmet. A friendly
    // caption asserting the cycle anyway would quietly contradict it.
    //
    // Swept over every state the module can be in rather than over one, so a
    // caption added later cannot slip the vocabulary in through a branch this
    // test does not visit.
    let forbidden = [
        "el niño",
        "el nino",
        "la niña",
        "la nina",
        "enso",
        "cycle",
        "cyclic",
        "period",
        "oscillat",
        "2-7",
        "2–7",
        "will ",
        "going to",
        "expect",
        "predict",
        "forecast",
        "next year",
    ];
    for reading in every_state() {
        let text = all_text(&reading);
        for word in forbidden {
            assert!(
                !text.contains(word),
                "a caption claims more than the model supports (\"{word}\"): {text}"
            );
        }
    }
}

#[test]
fn every_caption_is_a_sentence_that_reads_without_colour() {
    // The fourth criterion. A caption carries text and nothing else — there is
    // no colour on it to lose — and it says what it means in words, so it does
    // not lean on the picture beside it: "the red patch" is unreadable to a
    // reader who cannot see red and meaningless to one whose panel is off.
    let colour_words = ["red", "blue", "orange", "colour", "color", "shaded"];
    for reading in every_state() {
        let captions = reading.captions();
        assert!(
            captions.len() <= visualizer::MAX_CAPTIONS,
            "the panel reserves room for {} captions and was given {}",
            visualizer::MAX_CAPTIONS,
            captions.len()
        );
        for caption in captions {
            let text = caption.text();
            assert!(
                text.ends_with('.') && text.chars().next().is_some_and(char::is_uppercase),
                "a caption is a sentence: {text}"
            );
            let lowered = text.to_lowercase();
            for word in colour_words {
                assert!(!lowered.contains(word), "a caption names a colour: {text}");
            }
        }
    }
}

#[test]
fn a_reading_is_taken_from_the_run_rather_than_from_the_frame_index() {
    // The first criterion where it actually bites: the reading is measured off
    // the frames. The fixture writes T-07.4's equilibrium profile — +38.2 m at
    // the western wall falling to −28.2 m at the eastern — under the −0.05 Pa
    // of easterly stress the scenarios carry, so what is asserted here is the
    // profile that went in.
    let run = run_of(&steady_trades_header(20), |_| equilibrium_fields());
    let reading = EquatorialReading::of_run(&run, 19).expect("the run's own frames fit its grid");

    // The outermost cell centres sit half a cell inside each wall, so the
    // linear profile is read at 38.2 − 66.4/(2·320) and −28.2 + the same.
    let half_cell_m = (WESTERN_WALL_H_M - EASTERN_WALL_H_M) / (2.0 * NX as f64);
    assert!(
        (reading.west_depth_m.expect("the west is drawn")
            - (MEAN_DEPTH_M + WESTERN_WALL_H_M - half_cell_m))
            .abs()
            < METRE_TOLERANCE
    );
    assert!(
        (reading.east_depth_m.expect("the east is drawn")
            - (MEAN_DEPTH_M + EASTERN_WALL_H_M + half_cell_m))
            .abs()
            < METRE_TOLERANCE
    );
    assert!(
        (reading.equatorial_stress_pa - TRADE_STRESS_PA).abs() < PASCAL_TOLERANCE,
        "the stress the frames carry, averaged along the equator"
    );
    // The window is days of model time rather than a count of frames: this run
    // is daily and only nineteen days long, so nineteen days is as far back as
    // it can be read, and nineteen is what the captions quote.
    assert!((reading.lookback_days - 19.0).abs() < DAY_TOLERANCE);

    // And the captions say the profile back: 66 m of difference, held, under
    // trades at full strength.
    let text = all_text(&reading);
    assert!(text.contains("66 m"), "{text}");
    assert!(text.contains("holding"), "an equilibrium holds: {text}");
    assert!(text.contains("piling warm water into the west"), "{text}");
}

#[test]
fn the_same_index_of_two_different_runs_is_captioned_differently() {
    // The rule the ticket says matters most, stated as a test: nothing here is
    // keyed to a frame number. Two runs of the same length and cadence, one
    // level and one tilted, are read at the same index and captioned by what
    // is in them.
    let level = run_of(&steady_trades_header(20), |_| {
        let header = steady_trades_header(20);
        FrameFields::calm(&header)
    });
    let tilted_run = run_of(&steady_trades_header(20), |_| equilibrium_fields());

    let level_text = all_text(&EquatorialReading::of_run(&level, 19).expect("a reading"));
    let tilted_text = all_text(&EquatorialReading::of_run(&tilted_run, 19).expect("a reading"));

    assert!(level_text.contains("nearly the same depth"), "{level_text}");
    assert!(tilted_text.contains("66 m"), "{tilted_text}");
    assert_ne!(level_text, tilted_text);
}

#[test]
fn a_growing_tilt_is_read_as_growing_wherever_in_the_run_it_happens() {
    // The other half of the same rule. The tilt builds over this run, so the
    // trend caption says "growing" at every frame far enough in for a lookback
    // to exist — and that is a fact about the frames, not about where in the
    // run they sit.
    let frames = 60;
    let run = run_of(&steady_trades_header(frames), |index| {
        // A tilt that grows linearly from nothing to the equilibrium profile
        // over the run: what T-07.4's spin-up does, written by hand.
        #[allow(clippy::cast_precision_loss)]
        let share = index as f64 / (frames - 1) as f64;
        scaled_equilibrium_fields(share)
    });
    for index in 40..frames {
        let reading = EquatorialReading::of_run(&run, index).expect("a reading");
        let text = all_text(&reading);
        assert!(
            text.contains("growing"),
            "frame {index} is part of a tilt that is still building: {text}"
        );
        assert!(
            caption_on(&reading, CaptionTopic::WarmWaterEast).is_none(),
            "nothing is sliding back east while the tilt builds"
        );
    }
    // And near the start there is no window to look back over, so nothing is
    // claimed about a trend at all.
    assert!(caption_on(
        &EquatorialReading::of_run(&run, 0).expect("a reading"),
        CaptionTopic::TiltTrend
    )
    .is_none());
}

/// Every state the caption module can be asked about, for the sweeps above.
///
/// The winds at full strength, slack, reversed and absent, crossed with an
/// ocean level, tilted west-deep, tilted east-deep and broken down, crossed
/// with a tilt growing, shrinking, holding and unmeasured, plus a front
/// travelling and a front arrived.
fn every_state() -> Vec<EquatorialReading> {
    let winds = [TRADE_STRESS_PA, TRADE_STRESS_PA * 0.2, 0.04, 0.0, f64::NAN];
    let oceans = [
        (Some(MEAN_DEPTH_M), Some(MEAN_DEPTH_M)),
        (
            Some(MEAN_DEPTH_M + WESTERN_WALL_H_M),
            Some(MEAN_DEPTH_M + EASTERN_WALL_H_M),
        ),
        (
            Some(MEAN_DEPTH_M + EASTERN_WALL_H_M),
            Some(MEAN_DEPTH_M + WESTERN_WALL_H_M),
        ),
        (None, Some(MEAN_DEPTH_M)),
        (Some(MEAN_DEPTH_M), None),
    ];
    let histories = [None, Some(0.0), Some(66.4), Some(-66.4), Some(20.0)];
    let departures = [
        None,
        Some(Departure {
            peak_change_m: 20.0,
            half_peak_longitude_deg_east: -170.0,
            reaches_eastern_wall: false,
        }),
        Some(Departure {
            peak_change_m: 20.0,
            half_peak_longitude_deg_east: 179.5,
            reaches_eastern_wall: true,
        }),
        Some(Departure {
            peak_change_m: 0.1,
            half_peak_longitude_deg_east: 150.0,
            reaches_eastern_wall: false,
        }),
    ];
    let mut states = Vec::new();
    for stress in winds {
        for (west, east) in oceans {
            for earlier in histories {
                for departure in departures {
                    states.push(EquatorialReading {
                        west_depth_m: west,
                        east_depth_m: east,
                        earlier_tilt_m: earlier,
                        equatorial_stress_pa: stress,
                        departure,
                        ..at_rest()
                    });
                }
            }
        }
    }
    states
}

/// T-07.4's equilibrium: `h` falling linearly from +38.2 m at the western wall
/// to −28.2 m at the eastern, uniform in latitude, under the steady easterly
/// stress the scenarios carry.
fn equilibrium_fields() -> FrameFields {
    scaled_equilibrium_fields(1.0)
}

/// The same profile at `share` of its full amplitude, for a spin-up.
fn scaled_equilibrium_fields(share: f64) -> FrameFields {
    let header = steady_trades_header(1);
    let mut h_m = Vec::with_capacity(NX * NY);
    for _ in 0..NY {
        for i in 0..NX {
            #[allow(clippy::cast_precision_loss)]
            let fraction = (i as f64 + 0.5) / NX as f64;
            h_m.push(share * (WESTERN_WALL_H_M + (EASTERN_WALL_H_M - WESTERN_WALL_H_M) * fraction));
        }
    }
    FrameFields {
        h_m,
        tau_x_pa: vec![TRADE_STRESS_PA; header.grid.field_len(Variable::ZonalWindStress)],
        ..FrameFields::calm(&header)
    }
}

/// A run of `header`'s shape whose frames carry the fields `fields(index)`
/// gives.
fn run_of(header: &RunHeader, fields: impl Fn(u64) -> FrameFields) -> LoadedRun {
    let bytes = RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames_with_fields(header, header.output.frame_count, fields),
    };
    LoadedRun::from_bytes("run-captions", bytes).expect("the run loads")
}

// ---------------------------------------------------------------------------
// The three shipped scenarios, captioned by their own runs.
//
// Everything above fabricates a state and asks what is said about it. These
// take the presets a visitor actually presses (T-13.3), run them, and hold the
// sentences to what the scenario is: the control holds a tilt under a wind
// that does not change, the relaxation lets it go, and the burst sends
// something east that has not arrived everywhere at once. Each is the ticket's
// register example on the run it is an example of.
// ---------------------------------------------------------------------------

/// Steps taken per call while running a preset to the end — the chunk
/// `tests/scenario_presets.rs` uses, and it changes nothing about the run.
const CHUNK_STEPS: u64 = 4_096;

/// The preset named `name`, computed to the end of its schedule.
fn computed(name: &str) -> ComputedRun {
    let preset = ScenarioPreset::ALL
        .into_iter()
        .find(|preset| preset.name() == name)
        .unwrap_or_else(|| panic!("no preset is called {name}"));
    let mut run = ComputedRun::of_preset(preset, FrameBudget::browser())
        .unwrap_or_else(|error| panic!("{name} does not start: {error}"));
    while !run.is_finished() {
        run.advance_steps(CHUNK_STEPS)
            .unwrap_or_else(|error| panic!("{name} does not compute: {error}"));
    }
    run
}

/// Seconds in a day.
const SECONDS_PER_DAY: f64 = 86_400.0;

/// The frame of `run` nearest `day` of model time.
fn frame_at_day(run: &LoadedRun, day: f64) -> u64 {
    let interval_days = run.header().output.interval_s / SECONDS_PER_DAY;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let index = (day / interval_days).round() as u64;
    index.min(run.frame_count() - 1)
}

#[test]
fn the_control_preset_is_captioned_as_a_tilt_the_wind_is_holding() {
    // The first of the ticket's three examples, on the run it describes. The
    // control's wind never changes, so a year in the ocean is on the
    // equilibrium T-07.4 measured: the trades are piling water west, there is
    // a difference of tens of metres between the two ends, and it is not
    // moving any more.
    let computed = computed("Normal trade winds");
    let run = computed.run();
    let index = frame_at_day(run, 500.0);
    let reading = EquatorialReading::of_run(run, index).expect("a frame of the run it came from");
    let text = all_text(&reading);

    assert!(text.contains("piling warm water into the west"), "{text}");
    assert!(
        text.contains("deeper in the west than in the east"),
        "{text}"
    );
    assert!(
        text.contains("holding"),
        "a wind that does not change holds the tilt it made: {text}"
    );
    // And the metres are the run's own. The 2° browser grid is not a validated
    // resolution (`crate::presets`), so what is asserted is the order the
    // linear equilibrium sets — tens of metres of difference — rather than
    // T-07.4's half-degree figure.
    let tilt_m = reading.tilt_m().expect("both ends are drawn");
    assert!(
        (10.0..200.0).contains(&tilt_m),
        "the equilibrium tilt of this ocean is tens of metres, and it measured {tilt_m} m"
    );
}

#[test]
fn the_relaxed_winds_preset_is_captioned_as_the_water_sliding_back_east() {
    // The ticket's third example, and the compound state, on the run that
    // produces it. The wind fades over the middle of the Pacific for a season
    // from day 365; somewhere in that season the caption has to appear, and
    // before the wind ever changes it must not.
    let computed = computed("The winds relax");
    let run = computed.run();

    let quiet = EquatorialReading::of_run(run, frame_at_day(run, 200.0)).expect("a frame");
    assert!(
        caption_on(&quiet, CaptionTopic::WarmWaterEast).is_none(),
        "the trades are still at full strength on day 200: {}",
        all_text(&quiet)
    );

    // The relaxation peaks at day 365 with a 90-day `e`-folding time, so the
    // season the wind is at its weakest in is the two months either side of
    // it.
    let mut sliding = false;
    let mut slack = false;
    for day in [330.0, 345.0, 360.0, 375.0, 390.0] {
        let reading = EquatorialReading::of_run(run, frame_at_day(run, day)).expect("a frame");
        let text = all_text(&reading);
        sliding |= text.contains("sliding back east");
        slack |= text.contains("slack");
    }
    assert!(slack, "the winds do go slack over the season");
    assert!(sliding, "and the water that was piled up runs back east");
}

#[test]
fn the_burst_preset_names_a_longitude_the_change_has_not_reached_past() {
    // "Where a wave front has reached", on the run that has one — and held to
    // exactly what the caption claims, which is that the ocean east of the
    // longitude it names has barely moved. That is checked against the run's
    // own frames rather than against the departure the reading computed, so
    // the sentence is verified rather than restated.
    //
    // Five days after the burst peaks over the warm pool on day 180, the pulse
    // it launched is well short of the 69 days a Kelvin wave needs to cross
    // this basin, so there is ocean in the east that has not felt it.
    let computed = computed("A westerly wind burst");
    let run = computed.run();
    let index = frame_at_day(run, 185.0);
    let reading = EquatorialReading::of_run(run, index).expect("a frame");
    let departure = reading.departure.expect("the burst changed the ocean");
    let text = caption_on(&reading, CaptionTopic::Front).expect("a front was measured");

    assert!(
        text.contains("reaches as far east as") && text.contains("has barely moved"),
        "{text}"
    );
    assert!(
        !departure.reaches_eastern_wall,
        "the change has not run out of ocean five days after the burst"
    );

    // The claim itself: the eastern end of the equator moved by less than half
    // the largest change anywhere, over the same window the caption quotes.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let window_frames =
        (reading.lookback_days / (run.header().output.interval_s / SECONDS_PER_DAY)).round() as u64;
    let eastern_h_m = |index: u64| {
        let grid = run.header().grid;
        let (nx, ny) = (grid.nx(), grid.ny());
        let frame = run.frame(index).expect("a frame the run holds");
        // The browser basin is 25 rows of two degrees, so the middle row is
        // centred on the equator, and the last cell of it is the eastern end.
        frame.h()[(ny / 2) * nx + nx - 1]
    };
    let moved_m = (eastern_h_m(index) - eastern_h_m(index - window_frames)).abs();
    assert!(
        moved_m < 0.5 * departure.peak_change_m,
        "the caption says the east has barely moved, and it moved {moved_m} m against a peak of \
         {} m",
        departure.peak_change_m
    );
    // And the longitude it names is inside the basin rather than at its edge.
    let named = departure.half_peak_longitude_deg_east;
    assert!(
        (named - 120.0).rem_euclid(360.0) < 160.0,
        "{named}\u{b0}E is not a meridian of this basin"
    );
}
