//! Acceptance tests for T-06.2 — progress reporting and logging.
//!
//! The criteria are that *a run's progress output updates at a reasonable
//! cadence and reports a sane ETA*, and that *`--quiet` suppresses it for
//! scripted/CI use*. Both are checked twice over: against
//! [`engine::progress`]'s own types, where a synthetic clock makes the cadence
//! and the ETA exact, and against the real `termocline run` binary, where the
//! progress a user actually sees is the thing under test.
//!
//! # Why the clock is passed in rather than read
//!
//! Wall-clock time is the one input a test cannot control, so
//! [`RunProgress`](engine::progress::RunProgress) is offered a
//! [`ProgressReport`] whose elapsed duration is an argument, and only its
//! `RunObserver` implementation reads a real clock. Every assertion below about *when* a line
//! is drawn therefore holds exactly, rather than flakily: the durations fed in
//! are the ones the reporter reasons about.
//!
//! # Tolerances
//!
//! The ETA assertions have none. An ETA is arithmetic, not physics: with `k` of
//! `N` steps done in wall time `T`, the remaining work at the rate so far is
//! `T · (N − k) / k`, and the cases below are chosen so that product is an
//! exact number of seconds (10 s at a quarter done ⇒ 30 s remaining). The
//! expected values come from that formula, not from running the reporter.
//!
//! The line-count bounds derive from the reporter's two published cadence
//! constants — [`PLAIN_PROGRESS_PERCENT_STEP`] and
//! [`PROGRESS_REDRAW_INTERVAL`] — rather than from an observed count, so a
//! reporter that redrew every step or once a run fails here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{self, Command};
use std::time::Duration;

use engine::progress::{
    ProgressReport, ProgressStyle, RunProgress, Verbosity, PLAIN_PROGRESS_PERCENT_STEP,
    PROGRESS_REDRAW_INTERVAL,
};
use engine::{ScenarioConfig, FRAME_FILE_NAME, HEADER_FILE_NAME};

mod common;

use common::ScratchDir;

/// This file's ticket, which labels the directories it leaves in the system
/// temp directory.
const TICKET: &str = "t062";

/// The full length of a shipped example: two years at a one-hour step. The
/// number a progress bar exists for, and the length the synthetic runs below
/// report on.
const FULL_RUN_STEPS: u64 = 17_520;
/// The one-hour timestep of the shipped examples, in seconds.
const EXAMPLE_DT_S: f64 = 3600.0;

/// How long a full run takes in release, in seconds — the figure from the
/// ticket. Used only to make the synthetic wall clock below resemble a real
/// one; nothing asserted depends on its value.
const FULL_RUN_WALL_S: f64 = 40.0;

/// Steps a run driven through the real binary takes, and its output cadence:
/// one day of model time saved four times, as in `tests/engine_run.rs`. Short
/// enough to be a test, long enough to have progress to report.
const TEST_TOTAL_STEPS: u64 = 24;
const TEST_OUTPUT_EVERY_N_STEPS: u64 = 8;

/// Control characters no output may contain when it is not going to a
/// terminal: a carriage return redraws in place, and an escape introduces an
/// ANSI sequence. Either one turns a piped log or a CI transcript into noise.
const CONTROL_CHARACTERS: [char; 2] = ['\r', '\u{1b}'];

// ---------------------------------------------------------------------------
// The report: percent complete, and the ETA derived from it
// ---------------------------------------------------------------------------

#[test]
fn the_eta_is_the_remaining_work_at_the_rate_measured_so_far() {
    // A quarter of the run done in 10 s means three quarters left at the same
    // rate: 10 s · (17520 − 4380) / 4380 = 30 s. Exact, from the formula, not
    // from the reporter.
    let quarter = FULL_RUN_STEPS / 4;
    let report = ProgressReport::new(
        quarter,
        FULL_RUN_STEPS,
        quarter as f64 * EXAMPLE_DT_S,
        Duration::from_secs(10),
    );

    assert_eq!(report.eta(), Some(Duration::from_secs(30)));
    assert_eq!(report.percent_complete(), 25);
    // 4380 steps in 10 s.
    assert_eq!(report.steps_per_s(), Some(438.0));
}

#[test]
fn a_run_that_has_taken_no_step_yet_has_no_eta_to_report() {
    // Nothing has been measured, so any ETA would be invented. The report says
    // it does not know rather than printing a zero or an infinity.
    let report = ProgressReport::new(0, FULL_RUN_STEPS, 0.0, Duration::from_secs(3));

    assert_eq!(report.eta(), None);
    assert_eq!(report.steps_per_s(), None);
    assert_eq!(report.percent_complete(), 0);
    let line = report.render();
    assert!(
        !line.contains("eta"),
        "a report with no ETA should not print one; it said: {line}"
    );
}

#[test]
fn a_finished_run_is_a_hundred_percent_done_with_nothing_remaining() {
    let report = ProgressReport::new(
        FULL_RUN_STEPS,
        FULL_RUN_STEPS,
        FULL_RUN_STEPS as f64 * EXAMPLE_DT_S,
        Duration::from_secs_f64(FULL_RUN_WALL_S),
    );

    assert_eq!(report.percent_complete(), 100);
    assert_eq!(report.fraction_complete(), 1.0);
    assert_eq!(report.eta(), Some(Duration::ZERO));
}

#[test]
fn a_progress_line_names_the_percent_the_model_time_the_wall_time_and_the_eta() {
    // What the ticket asks the line to carry: percent complete, simulated
    // time, wall time, and an ETA. 4380 steps of one hour is 182.5 days of
    // model time.
    let quarter = FULL_RUN_STEPS / 4;
    let report = ProgressReport::new(
        quarter,
        FULL_RUN_STEPS,
        quarter as f64 * EXAMPLE_DT_S,
        Duration::from_secs(10),
    );

    let line = report.render();

    for expected in ["25%", "4380", "17520", "182.5 d", "eta"] {
        assert!(
            line.contains(expected),
            "a progress line should name {expected}; it said: {line}"
        );
    }
    assert_no_control_characters(&line, "a rendered progress line");
}

// ---------------------------------------------------------------------------
// The cadence
// ---------------------------------------------------------------------------

#[test]
fn plain_progress_reports_once_per_percent_step_and_not_once_per_step() {
    // The whole point of the criterion: a 17 520-step run must not produce
    // 17 520 lines, and must not produce one. The plain reporter draws when
    // the percent complete crosses a multiple of PLAIN_PROGRESS_PERCENT_STEP,
    // so a full run draws at most 100 / step of those plus the final line, and
    // — since every one of those thresholds is crossed — at least that many.
    let drawn = plain_progress_lines(&full_run_reports());

    // The thresholds a run actually crosses while it is still running are 10 %
    // through 90 %: the hundredth percent belongs to the finished run, whose
    // line is the final one. Nine plus one, written as the constant it comes
    // from so that changing the cadence changes this number with it.
    let expected = u64::from(100 / PLAIN_PROGRESS_PERCENT_STEP - 1) + 1;
    assert_eq!(
        drawn.len() as u64,
        expected,
        "a full run's plain progress should be one line per percent step plus a final one; \
         it drew:\n{}",
        drawn.join("\n")
    );

    // Each line is further along than the one before it, and the last one is
    // the finished run.
    let percents: Vec<u32> = drawn.iter().map(|line| percent_of(line)).collect();
    assert!(
        percents.windows(2).all(|pair| pair[0] < pair[1]),
        "progress should only ever move forwards; it drew {percents:?}"
    );
    assert_eq!(
        percents.last().copied(),
        Some(100),
        "the last line of a finished run should say it is finished; it drew {percents:?}"
    );
}

#[test]
fn plain_progress_carries_an_eta_while_the_run_is_still_going() {
    let drawn = plain_progress_lines(&full_run_reports());
    let (finished, running) = drawn
        .split_last()
        .expect("a full run draws at least a final line");

    for line in running {
        assert!(
            line.contains("eta"),
            "a line drawn mid-run should carry an ETA; it said: {line}"
        );
    }
    assert!(
        !finished.contains("eta"),
        "a finished run has nothing remaining to estimate; it said: {finished}"
    );
}

#[test]
fn plain_progress_emits_no_control_characters() {
    // The style a run gets when its output is a pipe or a CI log rather than a
    // terminal: whole lines, no cursor games.
    let written = String::from_utf8(plain_progress_bytes(&full_run_reports()))
        .expect("progress output is UTF-8");

    assert_no_control_characters(&written, "plain progress output");
    assert!(
        written.ends_with('\n'),
        "plain progress should end its last line; it said: {written:?}"
    );
}

#[test]
fn in_place_progress_redraws_one_line_no_faster_than_the_redraw_interval() {
    // The interactive style: a run of FULL_RUN_WALL_S seconds redraws at most
    // once per PROGRESS_REDRAW_INTERVAL, however many steps it takes in
    // between, plus the final draw. The observations below arrive far faster
    // than that interval, so a reporter that drew on every one of them would
    // draw 17 520 times.
    let mut progress = RunProgress::new(Vec::new(), ProgressStyle::InPlace, Verbosity::Normal);
    for report in full_run_reports() {
        progress.observe(&report);
    }
    progress.finish(&finished_run_report());
    let written = String::from_utf8(progress.into_writer()).expect("progress output is UTF-8");

    let draws = written.matches('\r').count() as u64;
    let interval_budget =
        (FULL_RUN_WALL_S / PROGRESS_REDRAW_INTERVAL.as_secs_f64()).floor() as u64 + 1;
    assert!(
        draws <= interval_budget,
        "a {FULL_RUN_WALL_S} s run should redraw at most {interval_budget} times at a \
         {PROGRESS_REDRAW_INTERVAL:?} interval; it drew {draws} times"
    );
    // It must still be a progress bar rather than a single line at the end.
    assert!(
        draws > 1,
        "an interactive run should redraw as it goes; it drew {draws} times"
    );
    assert_eq!(
        written.matches('\n').count(),
        1,
        "an in-place redraw occupies one line, ended once the run is over; it said: {written:?}"
    );
    assert!(
        !written.contains('\u{1b}'),
        "redrawing in place needs a carriage return, not an ANSI escape; it said: {written:?}"
    );
}

#[test]
fn a_quiet_reporter_writes_nothing_at_all() {
    let mut progress = RunProgress::new(Vec::new(), ProgressStyle::Plain, Verbosity::Quiet);
    for report in full_run_reports() {
        progress.observe(&report);
    }
    progress.finish(&finished_run_report());

    assert!(
        progress.into_writer().is_empty(),
        "--quiet exists so that a scripted run says nothing"
    );
}

// ---------------------------------------------------------------------------
// The command
// ---------------------------------------------------------------------------

#[test]
fn a_run_reports_its_progress_to_stderr_as_it_goes() {
    let scratch = ScratchDir::new(TICKET, "progress");
    let output = run_example(scratch.path(), &[]);

    assert!(
        output.status.success(),
        "`run` failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();

    assert!(
        stderr.contains('%'),
        "a run should report how far along it is; stderr said: {stderr}"
    );
    assert!(
        stderr.contains("100%"),
        "a finished run should report that it finished; stderr said: {stderr}"
    );
    assert!(
        stderr.contains("event=run_started") && stderr.contains("event=run_finished"),
        "a run should log its start and its end; stderr said: {stderr}"
    );
    // The summary line of T-06.1 is what a script reads; progress must not
    // have moved onto stdout beside it.
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(
        stdout.lines().count(),
        1,
        "stdout is the run's summary and nothing else; it said: {stdout}"
    );
    assert!(
        stdout.contains("frames written to"),
        "stdout should still carry T-06.1's summary line; it said: {stdout}"
    );
}

#[test]
fn a_piped_run_emits_no_control_characters() {
    // A run whose stderr is a pipe — this test, a CI job, a shell
    // redirection — gets whole lines. `Command::output` gives the child pipes
    // rather than a terminal, so this is that case exactly.
    let scratch = ScratchDir::new(TICKET, "piped");
    let output = run_example(scratch.path(), &[]);

    assert_no_control_characters(
        &String::from_utf8_lossy(&output.stderr),
        "the stderr of a piped run",
    );
    assert_no_control_characters(
        &String::from_utf8_lossy(&output.stdout),
        "the stdout of a piped run",
    );
}

#[test]
fn quiet_suppresses_the_progress_output_but_not_the_summary() {
    let scratch = ScratchDir::new(TICKET, "quiet");
    let output = run_example(scratch.path(), &["--quiet"]);

    assert!(output.status.success(), "a quiet run should still run");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.is_empty(),
        "--quiet should leave stderr empty; it said: {stderr}"
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("frames written to"),
        "--quiet suppresses progress, not the run's result; stdout said: {stdout}"
    );
}

#[test]
fn verbose_logs_every_frame_and_normal_does_not() {
    // The levels of the ticket: a frame written is the run's per-frame detail,
    // which belongs behind --verbose, while start and end are always worth a
    // line.
    let verbose_dir = ScratchDir::new(TICKET, "verbose");
    let verbose = run_example(verbose_dir.path(), &["--verbose"]);
    let verbose_stderr = String::from_utf8_lossy(&verbose.stderr).into_owned();

    let frames = TEST_TOTAL_STEPS / TEST_OUTPUT_EVERY_N_STEPS + 1;
    assert_eq!(
        verbose_stderr.matches("event=frame_written").count() as u64,
        frames,
        "--verbose should log each of the {frames} frames; stderr said: {verbose_stderr}"
    );
    assert!(
        verbose_stderr.contains("level=debug"),
        "per-frame detail is debug-level; stderr said: {verbose_stderr}"
    );

    let normal_dir = ScratchDir::new(TICKET, "normal");
    let normal = run_example(normal_dir.path(), &[]);
    let normal_stderr = String::from_utf8_lossy(&normal.stderr).into_owned();
    assert!(
        !normal_stderr.contains("event=frame_written"),
        "a default run should not log every frame; stderr said: {normal_stderr}"
    );
}

#[test]
fn quiet_and_verbose_cannot_be_asked_for_together() {
    // Silence and detail are contradictory instructions, so the CLI says so
    // rather than picking one (CODING_STANDARDS.md § *No silent clamping*).
    let scratch = ScratchDir::new(TICKET, "conflict");
    let output = run_example(scratch.path(), &["--quiet", "--verbose"]);

    assert!(
        !output.status.success(),
        "--quiet with --verbose should be refused"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--quiet") && stderr.contains("--verbose"),
        "the refusal should name both flags; it said: {stderr}"
    );
}

#[test]
fn how_loudly_a_run_reports_does_not_change_what_it_writes() {
    // CODING_STANDARDS.md § *Correctness and failure*: identical scenario in,
    // byte-identical output. Progress is a side channel, so a quiet run and a
    // loud one must leave the same two files behind.
    let loud_dir = ScratchDir::new(TICKET, "loud-files");
    let quiet_dir = ScratchDir::new(TICKET, "quiet-files");

    let loud = run_example(loud_dir.path(), &["--verbose"]);
    let quiet = run_example(quiet_dir.path(), &["--quiet"]);
    assert!(
        loud.status.success() && quiet.status.success(),
        "both runs should succeed"
    );

    for name in [HEADER_FILE_NAME, FRAME_FILE_NAME] {
        assert_eq!(
            fs::read(loud_dir.path().join("run").join(name)).expect("the loud run wrote its files"),
            fs::read(quiet_dir.path().join("run").join(name))
                .expect("the quiet run wrote its files"),
            "a quiet run and a verbose one disagree about {name}"
        );
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Every step of a synthetic full-length run, as the reporter sees it.
///
/// The wall clock is linear in the step count — a real run's is close to it,
/// since every step costs the same right-hand side — which is what makes the
/// ETA in the assertions above an exact number.
fn full_run_reports() -> Vec<ProgressReport> {
    (1..=FULL_RUN_STEPS)
        .map(|step| {
            ProgressReport::new(
                step,
                FULL_RUN_STEPS,
                step as f64 * EXAMPLE_DT_S,
                elapsed_at_step(step),
            )
        })
        .collect()
}

/// The report of that run once it has finished.
fn finished_run_report() -> ProgressReport {
    ProgressReport::new(
        FULL_RUN_STEPS,
        FULL_RUN_STEPS,
        FULL_RUN_STEPS as f64 * EXAMPLE_DT_S,
        elapsed_at_step(FULL_RUN_STEPS),
    )
}

/// Wall time a synthetic run has spent by `step`, at a constant rate.
fn elapsed_at_step(step: u64) -> Duration {
    Duration::from_secs_f64(FULL_RUN_WALL_S * step as f64 / FULL_RUN_STEPS as f64)
}

/// The bytes a plain reporter writes for `observations`, finishing the run.
fn plain_progress_bytes(reports: &[ProgressReport]) -> Vec<u8> {
    let mut progress = RunProgress::new(Vec::new(), ProgressStyle::Plain, Verbosity::Normal);
    for report in reports {
        progress.observe(report);
    }
    progress.finish(&finished_run_report());
    progress.into_writer()
}

/// The progress lines a plain reporter drew for `observations`.
fn plain_progress_lines(reports: &[ProgressReport]) -> Vec<String> {
    String::from_utf8(plain_progress_bytes(reports))
        .expect("progress output is UTF-8")
        .lines()
        .map(str::to_owned)
        .collect()
}

/// The percent a progress line reports.
///
/// # Panics
/// If the line has no percent in it, since every progress line must.
fn percent_of(line: &str) -> u32 {
    let (before, _) = line
        .split_once('%')
        .unwrap_or_else(|| panic!("a progress line names a percent; it said: {line}"));
    before
        .rsplit(|c: char| !c.is_ascii_digit())
        .next()
        .filter(|digits| !digits.is_empty())
        .unwrap_or_else(|| panic!("a progress line names a percent; it said: {line}"))
        .parse()
        .expect("a percent is a number")
}

/// Assert `text` is free of the control characters a piped log must not carry.
fn assert_no_control_characters(text: &str, what: &str) {
    for control in CONTROL_CHARACTERS {
        assert!(
            !text.contains(control),
            "{what} should carry no {control:?}; it said: {text:?}"
        );
    }
}

/// `termocline run` on a shortened `steady-trades`, with `flags`, writing into
/// `directory/run`.
fn run_example(directory: &Path, flags: &[&str]) -> process::Output {
    let config = directory.join("steady-trades.toml");
    fs::write(&config, shortened_example("steady-trades"))
        .expect("the scratch directory is writable");

    Command::new(env!("CARGO_BIN_EXE_termocline"))
        .arg("run")
        .arg("--config")
        .arg(&config)
        .arg("--out")
        .arg(directory.join("run"))
        .args(flags)
        .output()
        .expect("the engine binary is built before its integration tests run")
}

/// `stem`'s shipped example with its run shortened to [`TEST_TOTAL_STEPS`], as
/// the text of a scenario file.
fn shortened_example(stem: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("scenarios/{stem}.toml"));
    let source = fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{} is a shipped example: {error}", path.display()));
    let mut config = ScenarioConfig::from_toml(&source).expect("a shipped example is a scenario");
    config.run.total_steps = TEST_TOTAL_STEPS;
    config.run.output_every_n_steps = TEST_OUTPUT_EVERY_N_STEPS;
    config
        .to_toml()
        .expect("a scenario read from TOML can be written back to it")
}
