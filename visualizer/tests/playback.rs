//! T-09.2: the playback clock, without a window.
//!
//! Playback is the [`Scrubber`] driven by real time instead of by a hand: how
//! many frames a second of wall clock buys, what happens at the end of the
//! run, and what a pause does to the position. All three are arithmetic over
//! the same `u64` the scrubber already keeps, so — as with the scrubber and
//! the overlays (ADR-0006) — they are asserted here rather than watched.
//!
//! # Where the expected values come from
//!
//! The speed is a *definition*, not a measurement: at `f` frames a second, `t`
//! seconds of real time is `f · t` frames, and every count below is that
//! product worked out by hand. The run length is
//! `engine/scenarios/steady-trades.toml`'s 731 daily frames, as
//! `tests/scrubber.rs` derives it. Nothing here was read back off this code.

mod common;

use common::STEADY_TRADES_FRAMES;
use visualizer::{Playback, Scrubber, MAX_STALL_S, PLAYBACK_SPEEDS_FPS};

/// A scrubber over the control run, sitting where a freshly loaded run leaves
/// it.
fn over(frame_count: u64) -> Scrubber {
    let mut scrubber = Scrubber::new();
    scrubber.fit_to(frame_count);
    scrubber
}

/// One display frame of the repaint loop playback is stepped by.
///
/// A sixty-fourth of a second rather than the sixtieth a monitor actually
/// runs at, because 1/64 is exact in binary and 1/60 is not: every count
/// asserted below is then the product `f · t` worked out by hand, not that
/// product plus whatever a thousand inexact additions drifted by. The clock
/// under test does not care which of the two it is handed.
const DISPLAY_FRAME_S: f64 = 1.0 / 64.0;

/// Step `playback` through `seconds` of real time, one display frame at a
/// time, as the panel's repaint loop steps it.
fn run_for(playback: &mut Playback, scrubber: &mut Scrubber, seconds: f64) {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let steps = (seconds / DISPLAY_FRAME_S).round() as u64;
    for _ in 0..steps {
        playback.advance(scrubber, DISPLAY_FRAME_S);
    }
}

#[test]
fn a_freshly_loaded_run_is_not_playing() {
    assert!(!Playback::new().is_playing());
}

#[test]
fn a_paused_run_does_not_move_however_much_time_passes() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    run_for(&mut playback, &mut scrubber, 5.0);
    assert_eq!(scrubber.index(), 0);
}

#[test]
fn playing_advances_the_selected_number_of_frames_a_second() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    playback.set_frames_per_second(30.0);
    playback.play(&mut scrubber);
    run_for(&mut playback, &mut scrubber, 2.0);
    // 30 frames a second for two seconds is 60 frames on from frame 0.
    assert_eq!(scrubber.index(), 60);
}

#[test]
fn every_offered_speed_advances_at_the_speed_it_names() {
    // A second at f frames a second is f frames. Both halves of each pair are
    // written out here rather than one derived from the other, so that a speed
    // quietly changed in the menu fails the last assertion instead of moving
    // the expectation along with it.
    let named: [(f64, u64); 4] = [(5.0, 5), (15.0, 15), (30.0, 30), (60.0, 60)];
    for (fps, expected) in named {
        let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
        playback.set_frames_per_second(fps);
        playback.play(&mut scrubber);
        run_for(&mut playback, &mut scrubber, 1.0);
        assert_eq!(scrubber.index(), expected, "at {fps} frames a second");
    }
    assert_eq!(
        named.map(|(fps, _)| fps),
        PLAYBACK_SPEEDS_FPS,
        "the menu offers exactly the speeds asserted here"
    );
}

#[test]
fn a_speed_slower_than_the_repaint_rate_still_advances() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    // A frame every 200 ms: no single display frame is a whole run
    // frame, so a clock that dropped the remainder would never move at all.
    playback.set_frames_per_second(5.0);
    playback.play(&mut scrubber);
    run_for(&mut playback, &mut scrubber, 1.0);
    assert_eq!(scrubber.index(), 5);
}

#[test]
fn playback_stops_at_the_last_frame() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    playback.set_frames_per_second(60.0);
    playback.play(&mut scrubber);
    // 60 frames a second needs 731 / 60 s to reach the end of the run; twenty
    // seconds is well past it.
    run_for(&mut playback, &mut scrubber, 20.0);
    assert_eq!(scrubber.index(), 730);
    assert!(
        !playback.is_playing(),
        "the end of the run is where playback stops, not where it idles"
    );
}

#[test]
fn pausing_and_resuming_preserves_the_position() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    playback.set_frames_per_second(30.0);
    playback.play(&mut scrubber);
    run_for(&mut playback, &mut scrubber, 1.0);
    playback.pause();
    assert_eq!(scrubber.index(), 30);

    // Time passing while paused is time the run does not move through.
    run_for(&mut playback, &mut scrubber, 4.0);
    assert_eq!(scrubber.index(), 30);

    playback.play(&mut scrubber);
    assert_eq!(
        scrubber.index(),
        30,
        "resuming shows the frame it paused on"
    );
    run_for(&mut playback, &mut scrubber, 1.0);
    assert_eq!(scrubber.index(), 60, "and carries on from there");
}

#[test]
fn playing_from_the_last_frame_starts_the_run_again() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    scrubber.to_last();
    playback.play(&mut scrubber);
    assert_eq!(scrubber.index(), 0);
    assert!(playback.is_playing());
}

#[test]
fn a_stalled_repaint_loop_skips_no_further_than_one_stall() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    playback.set_frames_per_second(60.0);
    playback.play(&mut scrubber);
    // A backgrounded browser tab hands back one enormous gap when it wakes
    // (ADR-0006). It is one gap, not a fast-forward: at 60 frames a second the
    // documented quarter-second cap is 60 x 0.25 = 15 frames.
    playback.advance(&mut scrubber, 30.0);
    assert_eq!(scrubber.index(), 15);
    assert!(
        (MAX_STALL_S - 0.25).abs() < f64::EPSILON,
        "the 15 above is that cap times the speed; a different cap is a \
         different expectation, not a different assertion"
    );
}

#[test]
fn a_run_with_no_frames_never_starts_playing() {
    let (mut playback, mut scrubber) = (Playback::new(), over(0));
    playback.play(&mut scrubber);
    assert!(!playback.is_playing());
    run_for(&mut playback, &mut scrubber, 1.0);
    assert_eq!(scrubber.index(), 0);
}

#[test]
fn toggling_starts_a_paused_run_and_stops_a_playing_one() {
    let (mut playback, mut scrubber) = (Playback::new(), over(STEADY_TRADES_FRAMES));
    playback.toggle(&mut scrubber);
    assert!(playback.is_playing());
    playback.toggle(&mut scrubber);
    assert!(!playback.is_playing());
}

#[test]
fn the_offered_speeds_are_all_usable_and_ordered() {
    assert!(
        PLAYBACK_SPEEDS_FPS.windows(2).all(|pair| pair[0] < pair[1]),
        "the speed menu reads slowest first"
    );
    assert!(PLAYBACK_SPEEDS_FPS.iter().all(|fps| *fps > 0.0));
    assert!(
        PLAYBACK_SPEEDS_FPS.contains(&Playback::new().frames_per_second()),
        "the speed a run opens at is one the menu offers"
    );
}
