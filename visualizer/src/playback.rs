//! The playback clock: the [`Scrubber`] driven by real time rather than by a
//! hand on the slider.
//!
//! A run is 731 daily frames, and what those frames hold — a Kelvin wave
//! crossing the basin, a tilt collapsing — is *motion*. Stepped one frame at a
//! time it reads as 731 stills that happen to differ; played at thirty frames
//! a second it reads as the thing the model is of. That is the whole of what
//! this module adds, and it adds it without adding a second way to say which
//! frame is on screen: playback owns no index. It writes the scrubber's, the
//! same `u64` a drag writes, so every affordance downstream of the chooser —
//! the frame cache, the O(1) lookup, the run-wide colour scale — is on the
//! same path whether a frame arrived by hand or by clock.
//!
//! # Why the clock is a parameter
//!
//! [`Playback::advance`] is handed the seconds that have passed rather than
//! reading them. There is no instant stored here and nothing to subtract, so a
//! pause cannot bank time to spend on the resume: while paused the shell's
//! measurement is simply not spent. It also makes the one thing worth
//! asserting about playback — that `f` frames a second for `t` seconds is
//! `f · t` frames, that the end of the run stops it, that a pause holds the
//! position — arithmetic, and so assertable without a GPU or a wall clock, as
//! [ADR-0006] asks of everything the shell draws. `tests/playback.rs` is that
//! claim.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use crate::Scrubber;

/// Speeds the playback menu offers, in run frames per second of real time.
///
/// The scenario writes a frame a day (`steady-trades.toml`), so these are also
/// days a second: five is slow enough to follow one wave crest across the
/// basin, sixty puts the whole two-year run through in twelve seconds, and the
/// two in between are the ones a reader actually watches at. Named rather than
/// a free slider because a speed is chosen from across the room, not tuned.
pub const PLAYBACK_SPEEDS_FPS: [f64; 4] = [5.0, 15.0, 30.0, 60.0];

/// The speed a run opens at, in run frames per second of real time.
///
/// Thirty: the control run is 731 frames, which at this speed is a little
/// under twenty-five seconds — long enough to watch a tilt build, short enough
/// that a reader sees the end of the run without deciding to wait for it.
const DEFAULT_SPEED_FPS: f64 = 30.0;

/// The longest gap between repaints playback will spend, in seconds.
///
/// A browser tab that has been in the background hands back one enormous gap
/// when it wakes, and a minimised window does the same (ADR-0006). Spending it
/// would fast-forward the run by however long the reader was elsewhere, which
/// is not what they left it doing. A quarter of a second is longer than any
/// repaint interval a live window produces, so it caps only the gaps that are
/// stalls. Where a *repaint itself* takes longer than this — a large grid on a
/// slow machine — playback runs slower than the speed it names, which is the
/// honest outcome: the shell cannot show more frames a second than it can draw,
/// and skipping them to keep the clock would drop the frames the reader asked
/// to watch.
pub const MAX_STALL_S: f64 = 0.25;

/// Whether the run is playing, how fast, and how much of the next frame the
/// clock has already bought.
///
/// Starts paused: a run that began moving the moment it loaded would have
/// moved off frame zero before the reader had read the header, and it is also
/// what keeps a repaint that changes nothing changing nothing (`crate::app`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Playback {
    /// Whether real time is currently being spent on frames.
    playing: bool,
    /// Run frames per second of real time.
    frames_per_second: f64,
    /// The part of the next frame the clock has bought, in frames: always in
    /// `[0, 1)`.
    ///
    /// Without it, a speed slower than the repaint rate would never move at
    /// all — at five frames a second, no single display frame is a whole run
    /// frame, and a clock that dropped the remainder would drop all of it.
    carry_frames: f64,
}

impl Default for Playback {
    fn default() -> Self {
        Self::new()
    }
}

impl Playback {
    /// A paused clock at the default speed.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            playing: false,
            frames_per_second: DEFAULT_SPEED_FPS,
            carry_frames: 0.0,
        }
    }

    /// Whether the run is playing.
    #[must_use]
    pub const fn is_playing(&self) -> bool {
        self.playing
    }

    /// The speed, in run frames per second of real time.
    #[must_use]
    pub const fn frames_per_second(&self) -> f64 {
        self.frames_per_second
    }

    /// Play at `frames_per_second` run frames per second of real time.
    ///
    /// Takes effect from the next [`Playback::advance`]: the part-frame
    /// already bought was bought at the old speed, and re-pricing it would
    /// make a speed change jump the run.
    pub const fn set_frames_per_second(&mut self, frames_per_second: f64) {
        self.frames_per_second = frames_per_second;
    }

    /// Start playing `scrubber`'s run.
    ///
    /// From the last frame this starts the run again from its first: playback
    /// stops at the end rather than looping, so the only thing left for the
    /// button to mean there is "again". A run with no frames stays paused —
    /// there is nothing to play.
    pub const fn play(&mut self, scrubber: &mut Scrubber) {
        let Some(last) = scrubber.last() else {
            return;
        };
        if scrubber.index() == last {
            scrubber.to_first();
        }
        self.playing = true;
        self.carry_frames = 0.0;
    }

    /// Stop playing, holding the frame on screen.
    ///
    /// The part-frame the clock had bought is dropped rather than kept, so
    /// resuming starts a whole frame's interval — a resume that immediately
    /// jumped a frame would make a pause look like it had moved the run.
    pub const fn pause(&mut self) {
        self.playing = false;
        self.carry_frames = 0.0;
    }

    /// Play if paused, pause if playing.
    pub const fn toggle(&mut self, scrubber: &mut Scrubber) {
        if self.playing {
            self.pause();
        } else {
            self.play(scrubber);
        }
    }

    /// Spend `elapsed_s` seconds of real time on `scrubber`, stopping at the
    /// last frame of the run.
    ///
    /// Does nothing while paused, so time that passes then is time the run
    /// does not move through. Reaching the last frame both lands on it and
    /// stops: the run has an end, and running past it would either wrap the
    /// reader round to frame zero or sit there pretending to play.
    pub fn advance(&mut self, scrubber: &mut Scrubber, elapsed_s: f64) {
        if !self.playing {
            return;
        }
        let Some(last) = scrubber.last() else {
            self.playing = false;
            return;
        };
        let spent_s = elapsed_s.clamp(0.0, MAX_STALL_S);
        self.carry_frames += spent_s * self.frames_per_second;
        let whole = self.carry_frames.floor();
        self.carry_frames -= whole;
        // A speed and a gap are both finite and non-negative here, so the
        // count is too; `as` on it is the truncation of a value already whole.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let frames = whole as u64;
        if last - scrubber.index() <= frames {
            scrubber.to_last();
            self.pause();
            return;
        }
        scrubber.set_index(scrubber.index() + frames);
    }

    /// Draw the playback controls: the play/pause button, the speed menu, and
    /// the key that does the first of those without the mouse.
    ///
    /// `keyboard_free` says whether the keys are this control's to take, as it
    /// does for [`Scrubber::draw`] — a run URL is typed into a text field, and
    /// a space in it belongs to the caret.
    ///
    /// Draws nothing for a run with no frames: there is nothing to play. While
    /// playing it asks for the next repaint, because nothing else in the shell
    /// will — an idle panel repaints when the reader moves, and a run playing
    /// itself is the one thing on screen that moves when they do not.
    pub fn draw(&mut self, ui: &mut egui::Ui, scrubber: &mut Scrubber, keyboard_free: bool) {
        let Some(last) = scrubber.last() else {
            return;
        };
        let mut space_taken = false;
        ui.horizontal(|ui| {
            // Three labels, not two: paused on the last frame, the only thing
            // left for the button to do is start the run again, and it says so
            // rather than looking like a play button that has stopped working.
            let (label, hint) = match (self.playing, scrubber.index() == last) {
                (true, _) => ("⏸ Pause", "Pause (space)"),
                (false, false) => ("▶ Play", "Play (space)"),
                (false, true) => (
                    "↻ Replay",
                    "Play the run again from its first frame (space)",
                ),
            };
            let button = ui.button(label).on_hover_text(hint);
            if button.clicked() {
                self.toggle(scrubber);
            }
            // A focused button takes the space bar itself. Suppressing the key
            // here for *any* focused widget would kill it the moment the
            // reader touched the scrubber's slider, which does not want it —
            // so only this button's own focus counts, as `Scrubber::draw`
            // counts only the slider's for the arrow keys.
            space_taken = button.has_focus();

            ui.label("Speed:");
            let mut chosen = self.frames_per_second;
            egui::ComboBox::from_id_salt("playback-speed")
                .selected_text(speed_label(chosen))
                .show_ui(ui, |ui| {
                    for fps in PLAYBACK_SPEEDS_FPS {
                        ui.selectable_value(&mut chosen, fps, speed_label(fps));
                    }
                });
            // Through the setter, so the menu and the tests write the speed by
            // the same path and the part-frame rule lives on that path.
            self.set_frames_per_second(chosen);
        });
        if keyboard_free && !space_taken && ui.input(|input| input.key_pressed(egui::Key::Space)) {
            self.toggle(scrubber);
        }
        if self.playing {
            ui.ctx().request_repaint();
        }
    }
}

/// A speed as the menu says it.
///
/// In frames a second rather than in days a second: a frame is a day only
/// because `steady-trades.toml` writes one that often, and the run header is
/// what says how much model time a frame is worth.
fn speed_label(frames_per_second: f64) -> String {
    format!("{frames_per_second:.0} frames/s")
}
