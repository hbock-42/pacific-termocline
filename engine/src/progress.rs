//! What a run says about itself while it runs.
//!
//! A full example scenario is 17 520 steps and some 40 seconds of wall time, so
//! a run that printed nothing until it finished would be a black box for its
//! whole life. This module is the other side of that: a [`RunObserver`] the
//! runner calls at the few points a run has something to say, and one
//! implementation of it, [`RunProgress`], that turns those calls into a
//! progress line and a handful of structured log lines.
//!
//! # Why an observer rather than printing from the loop
//!
//! [`run_scenario`](crate::run::run_scenario) is a library function, and a
//! library that writes to a terminal has decided something that is not its to
//! decide. The runner therefore reports *events* — the run started, a step was
//! taken, a frame was written, the run finished — and the CLI chooses what to
//! do with them. `()` implements [`RunObserver`] as the observer that does
//! nothing, which is what a caller who wants a silent run passes.
//!
//! # The two styles, and why the choice is not the caller's taste
//!
//! A progress bar that redraws in place is the right thing on a terminal and
//! the wrong thing everywhere else: a run piped to a file, or into a CI log,
//! wants whole lines and no cursor control, because a carriage return in a
//! transcript is noise a reader has to undo. So [`ProgressStyle::of_terminal`]
//! picks the style from whether the stream is a terminal, and
//! [`ProgressStyle::Plain`] emits nothing but complete lines — no carriage
//! return, and no ANSI escape anywhere in this module.
//!
//! The two styles also differ in *cadence*, for the same reason. In place, the
//! bar redraws on a wall clock ([`PROGRESS_REDRAW_INTERVAL`]) because that is
//! what makes it look alive. In a log, the line is emitted once per
//! [`PLAIN_PROGRESS_PERCENT_STEP`] percent of the run, because a transcript
//! wants a bounded number of lines that mean something rather than one every
//! fifth of a second.
//!
//! # Where the clock is read
//!
//! Only in the [`RunObserver`] implementation. Every method that draws takes a
//! [`ProgressReport`] whose `elapsed` is an argument, so the cadence and the
//! ETA are functions of their inputs and a test can state both exactly.
//!
//! # Failure
//!
//! Progress is a side channel: it is written for a human watching, and it is
//! not part of what the run produces. A failed write — a closed pipe, a full
//! terminal buffer — is therefore discarded rather than propagated, because a
//! run that completed correctly did not fail on account of the progress bar
//! that described it.

use std::fmt::Write as _;
use std::io::{self, IsTerminal, Write};
use std::time::{Duration, Instant};

use crate::run::RunReport;
use crate::run_writer::OutputSchedule;

/// How often an in-place progress bar redraws, at most.
///
/// Five times a second: fast enough that the line reads as moving, slow enough
/// that a run at thousands of steps a second spends no measurable time
/// formatting it.
pub const PROGRESS_REDRAW_INTERVAL: Duration = Duration::from_millis(200);

/// How much of a run must complete between two plain progress lines, in
/// percent.
///
/// Ten, so a run of any length is described by ten lines and not by seventeen
/// thousand. A log is read after the fact, where what matters is that the run
/// was progressing and how fast, not the state of the bar at one moment.
pub const PLAIN_PROGRESS_PERCENT_STEP: u32 = 10;

/// How much a run says while it runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Verbosity {
    /// Nothing at all: no progress, no log lines. What `--quiet` asks for, so
    /// that a scripted or CI run's stderr carries only real problems.
    Quiet,
    /// Progress, and the events that bracket the run.
    #[default]
    Normal,
    /// The same, plus the run's per-frame detail. What `--verbose` asks for.
    Verbose,
}

impl Verbosity {
    /// Whether a message at `level` is printed at this verbosity.
    #[must_use]
    pub const fn prints(self, level: LogLevel) -> bool {
        match self {
            Self::Quiet => false,
            Self::Normal => matches!(level, LogLevel::Info),
            Self::Verbose => true,
        }
    }

    /// Whether progress is drawn at all at this verbosity.
    #[must_use]
    pub const fn draws_progress(self) -> bool {
        !matches!(self, Self::Quiet)
    }
}

/// The level of one log line.
///
/// Two levels, because the run has two kinds of thing to say: what a user
/// watching wants either way, and the per-frame detail they want when they are
/// debugging a run rather than waiting for one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    /// The run started; the run finished. Always worth a line.
    Info,
    /// Per-frame detail, behind `--verbose`.
    Debug,
}

impl LogLevel {
    /// The level as it appears in a log line's `level=` field.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Info => "info",
            Self::Debug => "debug",
        }
    }
}

/// How a progress line reaches its stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgressStyle {
    /// One line, redrawn in place with a carriage return, at
    /// [`PROGRESS_REDRAW_INTERVAL`]. For a terminal.
    InPlace,
    /// Whole lines, one per [`PLAIN_PROGRESS_PERCENT_STEP`] percent of the run,
    /// with no control characters. For a pipe, a file or a CI log.
    Plain,
}

impl ProgressStyle {
    /// The style for a stream, given whether it is a terminal.
    ///
    /// The whole rule: a terminal gets the bar, everything else gets lines.
    #[must_use]
    pub const fn of_terminal(is_terminal: bool) -> Self {
        if is_terminal {
            Self::InPlace
        } else {
            Self::Plain
        }
    }
}

/// How far along a run is, and how long it has taken to get there.
///
/// The one input a progress line is rendered from. `elapsed` is passed in
/// rather than read from a clock so that both the ETA and the cadence it feeds
/// are functions of their arguments.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProgressReport {
    /// Steps completed so far.
    steps_done: u64,
    /// Steps the run will take in total.
    total_steps: u64,
    /// Model time the run has reached, in seconds.
    model_time_s: f64,
    /// Wall time the run has spent so far.
    elapsed: Duration,
}

impl ProgressReport {
    /// A report of `steps_done` of `total_steps` completed, having reached
    /// `model_time_s` of model time in `elapsed` of wall time.
    #[must_use]
    pub const fn new(
        steps_done: u64,
        total_steps: u64,
        model_time_s: f64,
        elapsed: Duration,
    ) -> Self {
        Self {
            steps_done,
            total_steps,
            model_time_s,
            elapsed,
        }
    }

    /// Steps completed so far.
    #[must_use]
    pub const fn steps_done(self) -> u64 {
        self.steps_done
    }

    /// Steps the run will take in total.
    #[must_use]
    pub const fn total_steps(self) -> u64 {
        self.total_steps
    }

    /// Model time the run has reached, in seconds.
    #[must_use]
    pub const fn model_time_s(self) -> f64 {
        self.model_time_s
    }

    /// Wall time the run has spent so far.
    #[must_use]
    pub const fn elapsed(self) -> Duration {
        self.elapsed
    }

    /// The fraction of the run that is done, in `[0, 1]`.
    ///
    /// A run of no steps is finished before it starts, and reports 1 rather
    /// than dividing by zero.
    #[must_use]
    pub fn fraction_complete(self) -> f64 {
        if self.total_steps == 0 {
            return 1.0;
        }
        self.steps_done as f64 / self.total_steps as f64
    }

    /// The fraction of the run that is done, in percent, rounded down.
    ///
    /// Rounded down so that a run one step short of the end reports 99 and not
    /// 100: the last percent belongs to the run that has actually finished.
    #[must_use]
    pub fn percent_complete(self) -> u32 {
        // `fraction_complete` is in [0, 1], so the product is in [0, 100] and
        // the cast cannot saturate.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        {
            (self.fraction_complete() * 100.0).floor() as u32
        }
    }

    /// Steps a second, measured over the whole run so far, or `None` before the
    /// first step has finished.
    ///
    /// Averaged rather than instantaneous: every step of this scheme costs the
    /// same right-hand side, so the average is the better estimator of the
    /// next one and does not jitter with the machine's other work.
    #[must_use]
    pub fn steps_per_s(self) -> Option<f64> {
        let elapsed_s = self.elapsed.as_secs_f64();
        if self.steps_done == 0 || elapsed_s <= 0.0 {
            return None;
        }
        Some(self.steps_done as f64 / elapsed_s)
    }

    /// How much longer the run has to go at the rate measured so far, or
    /// `None` before the first step has finished.
    ///
    /// `elapsed · (total − done) / done`: the remaining work priced at the rate
    /// the run has actually managed. `None` before any step is done, because
    /// there is no rate yet and an ETA invented from nothing is worse than
    /// none.
    #[must_use]
    pub fn eta(self) -> Option<Duration> {
        if self.steps_done == 0 {
            return None;
        }
        let remaining = self.total_steps.saturating_sub(self.steps_done);
        Some(
            self.elapsed
                .mul_f64(remaining as f64 / self.steps_done as f64),
        )
    }

    /// The report as the line a progress reporter draws.
    ///
    /// Percent complete, steps, model time in days, wall time, and — while
    /// there is any run left — the ETA and the rate it was derived from, which
    /// is the set of numbers the ticket asks for. Free of control characters,
    /// so the same rendering serves a terminal and a log.
    ///
    /// Model time is shown in days rather than the seconds the run is
    /// integrated in: this is a line for a human watching a two-year run go
    /// past, not a record of the run, and the record — the header and the frame
    /// times — is in SI throughout.
    #[must_use]
    pub fn render(self) -> String {
        let mut line = String::new();
        // Writing into a `String` cannot fail.
        let _ = write!(
            line,
            "{:3}% | step {}/{} | model {:.1} d | elapsed {}",
            self.percent_complete(),
            self.steps_done,
            self.total_steps,
            self.model_time_s / SECONDS_PER_DAY,
            human_duration(self.elapsed),
        );
        // A finished run has an elapsed time, not an estimate; a run that has
        // taken no step yet has neither.
        if let Some(eta) = self.eta().filter(|eta| !eta.is_zero()) {
            let _ = write!(line, " | eta {}", human_duration(eta));
        }
        if let Some(rate) = self.steps_per_s() {
            let _ = write!(line, " | {rate:.0} steps/s");
        }
        line
    }
}

/// Seconds in a day, for rendering model time.
const SECONDS_PER_DAY: f64 = 86_400.0;

/// `elapsed` as a short human string: seconds under a minute, minutes and
/// seconds above it.
///
/// Rounded to tenths of a second *before* it is split into minutes and
/// seconds, because rounding afterwards renders 119.97 s as "1 m 60.0 s": the
/// seconds field rounds up past the minute the floor had already taken out of
/// it.
fn human_duration(elapsed: Duration) -> String {
    let tenths = (elapsed.as_secs_f64() * 10.0).round();
    let total_s = tenths / 10.0;
    if total_s < 60.0 {
        return format!("{total_s:.1} s");
    }
    let minutes = (total_s / 60.0).floor();
    let seconds = total_s - minutes * 60.0;
    format!("{minutes:.0} m {seconds:04.1} s")
}

/// What a run tells the outside world as it runs.
///
/// Every method has a do-nothing default, so an observer implements the events
/// it cares about; `()` implements the trait and observes nothing, which is the
/// silent run. Nothing here returns a `Result`: reporting is a side channel,
/// and a run does not fail because its narration did.
pub trait RunObserver {
    /// The run is about to take its first step.
    fn run_started(&mut self, _description: &str, _schedule: OutputSchedule) {}

    /// `steps_done` steps have been taken, reaching `model_time_s` of model
    /// time.
    fn step_taken(&mut self, _steps_done: u64, _model_time_s: f64) {}

    /// Frame `index` has been appended, holding the state at `t_s`.
    fn frame_written(&mut self, _index: u64, _t_s: f64) {}

    /// The run is over and its files are closed.
    fn run_finished(&mut self, _report: &RunReport) {}
}

/// The observer of a run nobody is watching.
impl RunObserver for () {}

/// Reports a run's progress and its events to a stream.
///
/// Holds the whole policy: what [`Verbosity`] prints, what [`ProgressStyle`]
/// draws and how often, and how a log line and a redrawn bar share one stream
/// without mangling each other.
#[derive(Debug)]
pub struct RunProgress<W: Write> {
    /// Where progress and log lines go. Usually stderr, so that the run's own
    /// output on stdout stays a machine-readable summary.
    writer: W,
    /// Whether progress redraws in place or appends lines.
    style: ProgressStyle,
    /// How much is printed.
    verbosity: Verbosity,
    /// Width of the in-place line currently on screen, so the next redraw can
    /// erase all of it. Zero when no bar is showing.
    drawn_width: usize,
    /// Wall time at the last in-place redraw, for [`PROGRESS_REDRAW_INTERVAL`].
    last_redraw: Option<Duration>,
    /// The last percent bucket a plain line was drawn for, for
    /// [`PLAIN_PROGRESS_PERCENT_STEP`].
    last_percent_bucket: u32,
    /// When the run started, for the [`RunObserver`] implementation — the one
    /// place in this module that reads a clock.
    started: Instant,
    /// The schedule the run is following, learned from
    /// [`RunObserver::run_started`], which is what turns a step count into a
    /// percentage.
    schedule: Option<OutputSchedule>,
}

impl<W: Write> RunProgress<W> {
    /// A reporter writing to `writer` in `style` at `verbosity`.
    pub fn new(writer: W, style: ProgressStyle, verbosity: Verbosity) -> Self {
        Self {
            writer,
            style,
            verbosity,
            drawn_width: 0,
            last_redraw: None,
            last_percent_bucket: 0,
            started: Instant::now(),
            schedule: None,
        }
    }

    /// Offer `report` to the reporter, which draws it if its cadence says so.
    ///
    /// A report of a *finished* run is never drawn here: [`Self::finish`] draws
    /// that one, and drawing it twice would put two 100 % lines in a log.
    pub fn observe(&mut self, report: &ProgressReport) {
        if !self.verbosity.draws_progress() || report.steps_done() >= report.total_steps() {
            return;
        }
        if self.is_due(report) {
            self.draw(report, false);
        }
    }

    /// Draw `report` as the run's last progress line, whatever the cadence
    /// says, and end the line.
    pub fn finish(&mut self, report: &ProgressReport) {
        if !self.verbosity.draws_progress() {
            return;
        }
        self.draw(report, true);
    }

    /// Write one structured log line, if `level` is printed at this verbosity.
    ///
    /// `fields` is the line's `key=value` body; the level is prepended. Any
    /// in-place bar on screen is erased first and redraws on the next tick, so
    /// a log line never lands in the middle of one.
    pub fn log(&mut self, level: LogLevel, fields: &str) {
        if !self.verbosity.prints(level) {
            return;
        }
        self.erase_bar();
        let _ = writeln!(self.writer, "level={} {fields}", level.as_str());
        let _ = self.writer.flush();
    }

    /// The stream, and whatever a `Vec<u8>` reporter has collected in it.
    #[must_use]
    pub fn into_writer(self) -> W {
        self.writer
    }

    /// Whether this style's cadence has come round for `report`.
    fn is_due(&self, report: &ProgressReport) -> bool {
        match self.style {
            ProgressStyle::InPlace => self.last_redraw.is_none_or(|last| {
                report.elapsed().saturating_sub(last) >= PROGRESS_REDRAW_INTERVAL
            }),
            ProgressStyle::Plain => {
                report.percent_complete() / PLAIN_PROGRESS_PERCENT_STEP > self.last_percent_bucket
            }
        }
    }

    /// Draw `report`, ending the line if the run is over.
    fn draw(&mut self, report: &ProgressReport, final_line: bool) {
        let line = report.render();
        self.last_redraw = Some(report.elapsed());
        self.last_percent_bucket = report.percent_complete() / PLAIN_PROGRESS_PERCENT_STEP;

        match self.style {
            ProgressStyle::InPlace => {
                // Erase by overwriting: the previous line's tail is blanked
                // with spaces, which needs no ANSI escape and so cannot leave
                // a sequence in a stream that turned out not to be a terminal
                // after all.
                let padding = self.drawn_width.saturating_sub(line.len());
                let _ = write!(self.writer, "\r{line}{:padding$}", "");
                self.drawn_width = if final_line { 0 } else { line.len() };
                if final_line {
                    let _ = writeln!(self.writer);
                }
            }
            ProgressStyle::Plain => {
                let _ = writeln!(self.writer, "{line}");
            }
        }
        let _ = self.writer.flush();
    }

    /// Blank the in-place bar on screen, if there is one, so something else may
    /// use the line.
    fn erase_bar(&mut self) {
        if self.drawn_width == 0 {
            return;
        }
        let _ = write!(self.writer, "\r{:width$}\r", "", width = self.drawn_width);
        self.drawn_width = 0;
    }
}

impl RunProgress<io::Stderr> {
    /// A reporter on stderr, in the style stderr deserves.
    ///
    /// Stderr rather than stdout because stdout is the run's result — T-06.1's
    /// one-line summary — and a script reading it should not have to filter a
    /// progress bar out of it. The style follows the same stream the progress
    /// goes to, so `2>log` gets plain lines while the terminal gets the bar.
    #[must_use]
    pub fn to_stderr(verbosity: Verbosity) -> Self {
        let stderr = io::stderr();
        let style = ProgressStyle::of_terminal(stderr.is_terminal());
        Self::new(stderr, style, verbosity)
    }
}

impl<W: Write> RunObserver for RunProgress<W> {
    fn run_started(&mut self, description: &str, schedule: OutputSchedule) {
        self.schedule = Some(schedule);
        self.started = Instant::now();
        self.log(
            LogLevel::Info,
            &format!(
                "event=run_started scenario={description} total_steps={} dt_s={} frames={}",
                schedule.total_steps(),
                schedule.dt_s(),
                schedule.frame_count(),
            ),
        );
    }

    fn step_taken(&mut self, steps_done: u64, model_time_s: f64) {
        let Some(schedule) = self.schedule else {
            return;
        };
        // The only allocation this makes per step is the one `render` makes
        // when the cadence is actually due; `observe` returns before that on
        // every other step (CODING_STANDARDS.md § *Performance*).
        self.observe(&ProgressReport::new(
            steps_done,
            schedule.total_steps(),
            model_time_s,
            self.started.elapsed(),
        ));
    }

    fn frame_written(&mut self, index: u64, t_s: f64) {
        // The verbosity is checked before the message is built, so a run that
        // is not logging frames allocates nothing per frame
        // (CODING_STANDARDS.md § *Performance*).
        if !self.verbosity.prints(LogLevel::Debug) {
            return;
        }
        self.log(
            LogLevel::Debug,
            &format!("event=frame_written index={index} t_s={t_s}"),
        );
    }

    fn run_finished(&mut self, report: &RunReport) {
        let elapsed = self.started.elapsed();
        let steps = report.steps_taken();
        // Without a schedule there is no model time to report, and a progress
        // line that stated one anyway would state a wrong one; the log line
        // below still says the run finished
        // (CODING_STANDARDS.md § *No silent clamping*).
        if let Some(schedule) = self.schedule {
            self.finish(&ProgressReport::new(
                steps,
                steps,
                schedule.model_time_at_step(steps),
                elapsed,
            ));
        }
        self.log(
            LogLevel::Info,
            &format!(
                "event=run_finished steps={steps} frames={} elapsed_s={:.3}",
                report.frames_written(),
                elapsed.as_secs_f64(),
            ),
        );
    }
}
