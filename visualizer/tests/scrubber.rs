//! T-08.3: the frame chooser's state, without a window.
//!
//! The scrubber is the part of stepping through a run that has a value in it —
//! which frame is chosen, and what each way of asking for another one does at
//! the ends of the run. Keeping it here rather than inside the panel is what
//! lets the ends of a run be asserted instead of dragged at (ADR-0006 asks the
//! same of everything else the shell draws).
//!
//! # Where the expected values come from
//!
//! The run length is `engine/scenarios/steady-trades.toml`'s: 17 520 steps of
//! an hour with a frame every 24 of them is 730 daily frames, and the frame at
//! t = 0 makes 731. What is asserted about it is the arithmetic of an index
//! into that many frames, not anything read back off this code.

use visualizer::Scrubber;

/// Frames `steady-trades.toml` writes: 17 520 / 24 = 730 daily frames, plus
/// the one at t = 0.
const STEADY_TRADES_FRAMES: u64 = 731;

/// A scrubber over a run of `frame_count` frames, sitting where a freshly
/// loaded run leaves it.
fn over(frame_count: u64) -> Scrubber {
    let mut scrubber = Scrubber::new();
    scrubber.fit_to(frame_count);
    scrubber
}

#[test]
fn a_freshly_loaded_run_starts_at_its_first_frame() {
    assert_eq!(over(STEADY_TRADES_FRAMES).index(), 0);
}

#[test]
fn the_last_frame_of_a_run_is_one_before_its_frame_count() {
    assert_eq!(over(STEADY_TRADES_FRAMES).last(), Some(730));
}

#[test]
fn a_run_with_no_frames_has_no_last_frame() {
    assert_eq!(over(0).last(), None);
}

#[test]
fn choosing_a_frame_past_the_end_of_the_run_lands_on_the_last_one() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.set_index(STEADY_TRADES_FRAMES);
    assert_eq!(scrubber.index(), 730);
}

#[test]
fn stepping_forward_moves_on_by_that_many_frames() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.step(30);
    assert_eq!(scrubber.index(), 30);
    scrubber.step(-7);
    assert_eq!(scrubber.index(), 23);
}

#[test]
fn stepping_past_the_end_of_the_run_stops_at_the_last_frame() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.to_last();
    scrubber.step(1);
    assert_eq!(scrubber.index(), 730);
}

#[test]
fn stepping_before_the_start_of_the_run_stops_at_the_first_frame() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.step(-1);
    assert_eq!(scrubber.index(), 0);
}

#[test]
fn the_ends_of_the_run_are_reachable_in_one_move() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.to_last();
    assert_eq!(scrubber.index(), 730);
    scrubber.to_first();
    assert_eq!(scrubber.index(), 0);
}

#[test]
fn a_shorter_run_pulls_the_chosen_frame_back_to_its_last() {
    let mut scrubber = over(STEADY_TRADES_FRAMES);
    scrubber.to_last();
    scrubber.fit_to(10);
    assert_eq!(scrubber.index(), 9);
}

#[test]
fn a_run_with_no_frames_leaves_nothing_to_choose() {
    let mut scrubber = over(0);
    scrubber.step(5);
    scrubber.to_last();
    assert_eq!(scrubber.index(), 0);
    assert_eq!(scrubber.last(), None);
}
