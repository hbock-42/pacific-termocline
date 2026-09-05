//! The frame chooser: which frame of the run is on screen, and every way of
//! choosing another one.
//!
//! It is a `u64` and a bound, kept apart from the panel that draws it for the
//! same reason [`crate::heatmap`] is: what a reader would check by dragging —
//! that the ends of the run hold, that a shorter run cannot leave the chooser
//! pointing past its last frame — is arithmetic, and arithmetic is assertable
//! without a GPU. `tests/scrubber.rs` is that claim.
//!
//! Nothing here decodes a frame. The chooser names one; [`crate::LoadedRun`]
//! is what fetches it, and it is indexed so that fetching frame 730 costs what
//! fetching frame 0 costs.
//!
//! [`Scrubber::draw`] is the control itself — the slider, the steps either
//! side of it, and the keys that do the same without the mouse. It is here
//! rather than in [`crate::app`] because every one of those affordances is a
//! way of writing the same `u64`, and the panel that draws a frame should not
//! have to know how many of them there are.

/// The frame of a run the reader has chosen.
///
/// Clamped rather than fallible: every way of moving it saturates at the ends
/// of the run, because a scrubber dragged past the end of a run means "the
/// last frame", not an error. A chooser over a run of no frames sits at index
/// zero and has no [`Scrubber::last`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Scrubber {
    /// The chosen frame, by index into the run. Never past `frame_count`.
    index: u64,
    /// Frames in the run being scrubbed, as its header counts them.
    frame_count: u64,
}

impl Scrubber {
    /// A chooser over no run at all, sitting on the frame a run's first is.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            index: 0,
            frame_count: 0,
        }
    }

    /// Point this chooser at a run of `frame_count` frames, pulling the chosen
    /// frame back if that run is too short to hold it.
    ///
    /// Called every time the panel draws rather than only when a run is
    /// loaded: the frame count is the run's, and the chooser is the panel's,
    /// so this is where the two are made to agree.
    pub const fn fit_to(&mut self, frame_count: u64) {
        self.frame_count = frame_count;
        if let Some(last) = self.last() {
            if self.index > last {
                self.index = last;
            }
        } else {
            self.index = 0;
        }
    }

    /// The chosen frame, by index into the run.
    #[must_use]
    pub const fn index(&self) -> u64 {
        self.index
    }

    /// The last frame of the run, or `None` if it holds none.
    #[must_use]
    pub const fn last(&self) -> Option<u64> {
        self.frame_count.checked_sub(1)
    }

    /// Choose frame `index`, or the last frame if the run is shorter than that.
    pub const fn set_index(&mut self, index: u64) {
        self.index = match self.last() {
            Some(last) if index > last => last,
            Some(_) => index,
            None => 0,
        };
    }

    /// Move on by `frames`, stopping at whichever end of the run is reached.
    ///
    /// Negative steps back. Saturating, not wrapping: a reader stepping off
    /// the end of a run means to be at the end of it, and wrapping round to
    /// the other end would show them a frame two years away from the one they
    /// were looking at.
    pub const fn step(&mut self, frames: i64) {
        let moved = if frames < 0 {
            self.index.saturating_sub(frames.unsigned_abs())
        } else {
            self.index.saturating_add(frames.unsigned_abs())
        };
        self.set_index(moved);
    }

    /// Choose the run's first frame.
    pub const fn to_first(&mut self) {
        self.index = 0;
    }

    /// Choose the run's last frame, if it has one.
    pub const fn to_last(&mut self) {
        if let Some(last) = self.last() {
            self.index = last;
        }
    }

    /// Draw the scrubber: the slider, the steps either side of it, and the
    /// keys that do the same thing without the mouse.
    ///
    /// `keyboard_free` says whether the keys are this control's to take — a
    /// run URL is typed into a text field, and an arrow key inside it belongs
    /// to the caret. Every control here changes one number and nothing else;
    /// what puts a frame on screen is the panel reading that number back, so a
    /// drag, an arrow key and a jump to the end of the run all cost the same.
    ///
    /// Draws nothing for a run with no frames: there is nothing to choose.
    pub fn draw(&mut self, ui: &mut egui::Ui, keyboard_free: bool) {
        let Some(last) = self.last() else {
            return;
        };
        let mut arrows_taken = false;
        ui.horizontal(|ui| {
            let index = self.index;
            if step_button(ui, "⏮", "First frame (Home)", index > 0) {
                self.to_first();
            }
            if step_button(ui, "◀", "Back one frame (left arrow)", index > 0) {
                self.step(-1);
            }
            if step_button(ui, "▶", "On one frame (right arrow)", index < last) {
                self.step(1);
            }
            if step_button(ui, "⏭", "Last frame (End)", index < last) {
                self.to_last();
            }
            // The slider takes the whole rest of the row: it is dragged, and a
            // frame of a long run is worth more pixels than a short slider
            // gives it — at 731 frames a narrow one puts several frames under
            // every pixel and none of them within reach.
            let mut chosen = self.index;
            let slider = ui.add_sized(
                egui::vec2(ui.available_width(), ui.spacing().interact_size.y),
                egui::Slider::new(&mut chosen, 0..=last)
                    .integer()
                    .show_value(false),
            );
            if slider.changed() {
                self.set_index(chosen);
            }
            // A focused slider steps itself on the arrow keys. Taking them
            // here as well would move two frames for one press.
            arrows_taken = slider.has_focus();
        });
        if keyboard_free {
            self.take_keys(ui.ctx(), arrows_taken);
        }
    }

    /// Move the chooser by whatever the keyboard asked for this frame.
    fn take_keys(&mut self, ctx: &egui::Context, arrows_taken: bool) {
        let pressed = |key| ctx.input(|input| input.key_pressed(key));
        if pressed(egui::Key::Home) {
            self.to_first();
        }
        if pressed(egui::Key::End) {
            self.to_last();
        }
        let steps = [
            (egui::Key::PageUp, -FRAMES_PER_PAGE),
            (egui::Key::PageDown, FRAMES_PER_PAGE),
        ];
        let arrows = [(egui::Key::ArrowLeft, -1), (egui::Key::ArrowRight, 1)];
        for (key, frames) in steps
            .into_iter()
            .chain(arrows.into_iter().filter(|_| !arrows_taken))
        {
            if pressed(key) {
                self.step(frames);
            }
        }
    }
}

/// Frames a page key moves the chooser by.
///
/// Ten, against the arrow keys' one: the scenario writes a frame a day
/// (`steady-trades.toml`), so a page is a week and a half of model time — far
/// enough to see a change, near enough to still be reading the same event.
const FRAMES_PER_PAGE: i64 = 10;

/// One of the scrubber's step buttons: what it says, what it does when hovered
/// over, and whether there is anywhere for it to go.
fn step_button(ui: &mut egui::Ui, label: &str, hint: &str, enabled: bool) -> bool {
    ui.add_enabled(enabled, egui::Button::new(label))
        .on_hover_text(hint)
        .clicked()
}
