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
}
