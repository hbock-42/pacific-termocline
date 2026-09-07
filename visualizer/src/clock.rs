//! When the run is, in the units a reader lives in: days, months and years,
//! taken from the cadence the run's own header declares.
//!
//! A frame index is a fact about the file. "Frame 396 of 731" tells a visitor
//! nothing about the ocean, and it is not even a fixed amount of time: the
//! interval between frames is a scenario's to choose
//! ([`termocline_format::OutputTiming`]), so the same index is a different
//! moment in two runs. What a reader needs is the model time the frame sits
//! at, and the model time is `index · interval`.
//!
//! Months and years rather than only days because the phenomena are that long:
//! the seasonal cycle is a year, and the run a visitor opens is two of them
//! (`CONTEXT.md`, *Seasonal cycle*). A caption saying "day 396" alone leaves
//! the reader to divide.
//!
//! Device-free like everything beside it, so `tests/teaching_labels.rs`
//! asserts the phrases rather than looking at them ([ADR-0006]).
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use engine::TROPICAL_YEAR_S;
use termocline_format::OutputTiming;

use crate::run::SECONDS_PER_DAY;

/// The tropical year, in days: 365.2422, the equinox-to-equinox year the
/// seasons follow. Taken from the engine's own [`TROPICAL_YEAR_S`] rather than
/// restated, so a caption and the seasonal forcing it describes cannot come to
/// mean different years.
const DAYS_PER_YEAR: f64 = TROPICAL_YEAR_S / SECONDS_PER_DAY;

/// A twelfth of that year, in days: 30.4368.
///
/// The month a reader counts in is a twelfth of a year and not any particular
/// calendar month, because the model has no calendar — only elapsed time from
/// the start of the run.
const DAYS_PER_MONTH: f64 = DAYS_PER_YEAR / 12.0;

/// Below this many days elapsed, a duration is said in days.
///
/// Two months of a tropical year, to the day: past that the day count is a
/// number a reader has to divide, and below it "2 months" throws away the
/// difference between a run at day 31 and one at day 59.
const DAYS_SAID_IN_DAYS: f64 = 2.0 * DAYS_PER_MONTH;

/// And below this many, in months rather than in years.
///
/// Eighteen months. A year and a half is where a count of months stops being a
/// number a reader holds — "24 months" is arithmetic, "2 years" is a fact — and
/// it is above the twelve months at which a run is not yet a year old, so a
/// run part-way through its first year is still said in the months it has run.
const DAYS_SAID_IN_MONTHS: f64 = 18.0 * DAYS_PER_MONTH;

/// The clock of one run: how far apart its frames are in model time, and how
/// many of them there are.
///
/// A copy of the run's [`OutputTiming`] and nothing else, so a panel of a
/// comparison dates its frames by *its own* run's cadence — two runs written
/// at different intervals sit at different days at the same index, which is
/// exactly what [`crate::Mismatch::Cadence`] is about.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunClock {
    /// Model time between consecutive frames, in seconds.
    interval_s: f64,
    /// How many frames the run declares.
    frame_count: u64,
}

impl RunClock {
    /// The clock of a run whose header declares `output`.
    #[must_use]
    pub const fn of_run(output: OutputTiming) -> Self {
        Self {
            interval_s: output.interval_s,
            frame_count: output.frame_count,
        }
    }

    /// The model day frame `index` sits at, counted from the start of the run.
    #[must_use]
    pub fn day_of_frame(&self, index: u64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let intervals = index as f64;
        intervals * self.interval_s / SECONDS_PER_DAY
    }

    /// The model day of the last frame: the length of the run, in days.
    ///
    /// The span of `n` frames is `n - 1` intervals, as
    /// [`crate::LoadedRun::metadata`] counts it: a run of one frame spans no
    /// time at all.
    #[must_use]
    pub fn last_day(&self) -> f64 {
        self.day_of_frame(self.frame_count.saturating_sub(1))
    }

    /// Where in the run frame `index` is, in words: `Day 396 of 730 — 13
    /// months in`.
    ///
    /// The months are dropped for a run still in its first weeks, where they
    /// would only repeat the day count in a coarser unit.
    #[must_use]
    pub fn moment(&self, index: u64) -> String {
        self.day_phrase(self.day_of_frame(index))
    }

    /// The same for a day of the run named directly rather than by a frame:
    /// `Day 396 of 730 — 13 months in`.
    #[must_use]
    pub fn day_phrase(&self, day: f64) -> String {
        let head = format!("Day {} of {}", days_text(day), days_text(self.last_day()));
        if day < DAYS_SAID_IN_DAYS {
            return head;
        }
        format!("{head} — {} in", plain_duration(day))
    }

    /// How long the run is, in words: `730 days (2 years)`.
    ///
    /// Both units, because the days are what the frames are counted in and the
    /// years are what a reader has intuition about. A run too short for the
    /// second to say anything new says only the first.
    #[must_use]
    pub fn span_phrase(&self) -> String {
        let days = self.last_day();
        let head = format!("{} days", days_text(days));
        if days < DAYS_SAID_IN_DAYS {
            return head;
        }
        format!("{head} ({})", plain_duration(days))
    }
}

/// `days` of model time in the coarsest unit that still says something: days
/// for the first two months, then months, then years.
#[must_use]
pub fn plain_duration(days: f64) -> String {
    if days < DAYS_SAID_IN_DAYS {
        return format!("{} days", days_text(days));
    }
    if days < DAYS_SAID_IN_MONTHS {
        return format!("{:.0} months", days / DAYS_PER_MONTH);
    }
    let years = days / DAYS_PER_YEAR;
    let text = trimmed(format!("{years:.1}"));
    if text == "1" {
        return "1 year".to_owned();
    }
    format!("{text} years")
}

/// A number of days, at the precision the number itself carries: whole days
/// for a run written daily or slower, and a tenth or a hundredth for one
/// written faster, so a sub-daily cadence is not reported as day zero
/// throughout.
fn days_text(days: f64) -> String {
    if days >= 10.0 {
        return trimmed(format!("{days:.0}"));
    }
    if days >= 1.0 {
        return trimmed(format!("{days:.1}"));
    }
    trimmed(format!("{days:.2}"))
}

/// `text` with a trailing `.0` or `.00` dropped: a whole number is written
/// whole.
fn trimmed(text: String) -> String {
    match text.strip_suffix(".00").or_else(|| text.strip_suffix(".0")) {
        Some(whole) => whole.to_owned(),
        None => text,
    }
}
