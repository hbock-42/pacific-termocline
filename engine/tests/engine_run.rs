//! Acceptance tests for T-06.1 — the `run` command.
//!
//! The criterion is that *running each of the 3 example scenario configs from
//! T-03.4 end-to-end produces a valid, readable run directory*. So each of
//! `engine/scenarios/*.toml` is driven through the real binary, as a
//! subprocess, and the directory it leaves behind is opened with the
//! `termocline-format` reader — the same reader the visualizer will use — and
//! checked field by field against the scenario file it came from.
//!
//! # Why the runs here are shortened
//!
//! The shipped examples are two-year runs on a 320 × 100 basin: 17 520 steps
//! of the full right-hand side, three times over, which is not a test suite.
//! Each config is therefore read, its `[run]` section shortened to
//! [`TEST_TOTAL_STEPS`] steps, and written back out through
//! `ScenarioConfig::to_toml` — so the basin, the physics and every `[[wind]]`
//! entry the example ships are the ones that run, and only the run's *length*
//! is different. The shipped length is still exercised: `[run]`'s own numbers
//! are validated as written, and the frame count the full config promises is
//! checked against the cadence arithmetic by hand.
//!
//! # Tolerances
//!
//! Almost nothing here has one. A run directory is a transcription — the
//! header's parameters are the scenario's numbers, the frame times are
//! multiples of the output interval, the field lengths are the grid's — so
//! those are compared exactly, as IEEE-754 values, because a header that was
//! merely "close" to its scenario would misdescribe the run for good.
//!
//! The one physical assertion, on the sign of the thermocline tilt, has no
//! tolerance either: it is a statement about a sign, and the value it is
//! compared against is exactly zero.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{self, Command};

use engine::{Scenario, ScenarioConfig, Staggering, FRAME_FILE_NAME, HEADER_FILE_NAME};
use termocline_format::{RunReader, FORMAT_VERSION};

mod common;

use common::ScratchDir;

/// This file's ticket, which labels the directories it leaves in the system
/// temp directory.
const TICKET: &str = "t061";

/// The three example scenarios of T-03.4, by file stem. The list is written
/// out rather than read from the directory so that an example silently
/// deleted, or renamed, fails the build.
const EXAMPLE_STEMS: [&str; 3] = ["steady-trades", "seasonal-cycle", "wind-burst"];

/// Steps a shortened example takes here. One day of model time at the
/// examples' one-hour timestep: long enough that the trades have accelerated
/// the surface layer and the thermocline has begun to tilt, and far shorter
/// than the ≈ 69 days a Kelvin wave needs to cross the 17 800 km basin at
/// 3 m/s, so nothing has reflected off the eastern wall yet.
const TEST_TOTAL_STEPS: u64 = 24;
/// Output cadence of a shortened example, in steps: four frames from the run,
/// which is a decimated series rather than every step.
const TEST_OUTPUT_EVERY_N_STEPS: u64 = 8;

/// Frames a shortened example writes: the initial state plus one per interval
/// that fits. Written out by hand rather than taken from the schedule so that
/// the schedule and the test cannot drift together.
const TEST_FRAME_COUNT: u64 = TEST_TOTAL_STEPS / TEST_OUTPUT_EVERY_N_STEPS + 1;

/// The `[run]` numbers every shipped example carries: a one-hour step, two
/// years of it, saved once a day. Written here so that an example that
/// changes its run length announces it.
const EXAMPLE_DT_S: f64 = 3600.0;
const EXAMPLE_TOTAL_STEPS: u64 = 17_520;
const EXAMPLE_OUTPUT_EVERY_N_STEPS: u64 = 24;

/// The text of one shipped example.
fn example_source(stem: &str) -> String {
    let path = example_path(stem);
    fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{} is a shipped example: {error}", path.display()))
}

/// Where a shipped example lives, as an absolute path, so a test may hand it
/// to a subprocess whose working directory is not this crate's.
fn example_path(stem: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("scenarios/{stem}.toml"))
}

/// `stem`'s example with its run shortened to [`TEST_TOTAL_STEPS`], as the
/// text of a scenario file.
///
/// Everything else — the basin, the physics, the ordered `[[wind]]` list — is
/// the example's own, so what runs below is the shipped scenario at a
/// different length and nothing else.
fn shortened_example(stem: &str) -> String {
    let mut config =
        ScenarioConfig::from_toml(&example_source(stem)).expect("a shipped example is a scenario");
    config.run.total_steps = TEST_TOTAL_STEPS;
    config.run.output_every_n_steps = TEST_OUTPUT_EVERY_N_STEPS;
    config
        .to_toml()
        .expect("a scenario read from TOML can be written back to it")
}

/// The `run` subcommand of the engine binary.
fn run_cli(config: &Path, out: &Path) -> process::Output {
    Command::new(env!("CARGO_BIN_EXE_termocline"))
        .arg("run")
        .arg("--config")
        .arg(config)
        .arg("--out")
        .arg(out)
        .output()
        .expect("the engine binary is built before its integration tests run")
}

/// Write `source` as `<name>.toml` inside `directory` and run it into
/// `directory/run`, returning the run directory.
///
/// The file is named rather than fixed because the header records the
/// scenario's name, so the name is part of what the run is checked against.
///
/// # Panics
/// If the command failed, with both of its streams, since every use below
/// expects a run that succeeded.
fn run_scenario_text(source: &str, directory: &Path, name: &str) -> PathBuf {
    let config = directory.join(format!("{name}.toml"));
    fs::write(&config, source).expect("the scratch directory is writable");
    let out = directory.join("run");

    let output = run_cli(&config, &out);
    assert!(
        output.status.success(),
        "`run` failed on {}:\nstdout: {}\nstderr: {}",
        config.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    out
}

#[test]
fn every_shipped_example_is_a_scenario_the_engine_will_run() {
    // The examples as they ship, at their full two-year length: the criterion
    // is about these files, so the shortening below must not be what makes
    // them valid.
    for stem in EXAMPLE_STEMS {
        let scenario = Scenario::from_toml(&example_source(stem))
            .unwrap_or_else(|error| panic!("{stem}.toml is not a runnable scenario: {error}"));

        let schedule = scenario.output_schedule();
        assert_eq!(schedule.dt_s(), EXAMPLE_DT_S, "{stem}.toml: dt");
        assert_eq!(
            schedule.total_steps(),
            EXAMPLE_TOTAL_STEPS,
            "{stem}.toml: run length"
        );
        // 17 520 / 24 + 1: a frame a day for 730 days, plus the initial
        // state. Written out rather than read off the schedule.
        assert_eq!(
            schedule.frame_count(),
            EXAMPLE_TOTAL_STEPS / EXAMPLE_OUTPUT_EVERY_N_STEPS + 1,
            "{stem}.toml: frame count"
        );
        assert_eq!(
            schedule.interval_s(),
            EXAMPLE_OUTPUT_EVERY_N_STEPS as f64 * EXAMPLE_DT_S,
            "{stem}.toml: output interval"
        );
    }
}

#[test]
fn each_example_scenario_produces_a_readable_run_directory() {
    for stem in EXAMPLE_STEMS {
        let scratch = ScratchDir::new(TICKET, stem);
        let source = shortened_example(stem);
        let out = run_scenario_text(&source, scratch.path(), stem);

        assert!(
            out.join(HEADER_FILE_NAME).is_file(),
            "{stem}: the run has no header"
        );
        assert!(
            out.join(FRAME_FILE_NAME).is_file(),
            "{stem}: the run has no frames"
        );

        let scenario = Scenario::from_toml(&source).expect("the shortened example is a scenario");
        let bounds = scenario.bounds();
        let params = scenario.physical_params();

        let reader = RunReader::open(&out).unwrap_or_else(|error| {
            panic!("{stem}: the run directory does not read back: {error}")
        });
        let header = reader.header().clone();

        assert_eq!(header.format_version, FORMAT_VERSION, "{stem}: version");
        assert_eq!(header.grid.nx(), bounds.nx(), "{stem}: cells east-west");
        assert_eq!(header.grid.ny(), bounds.ny(), "{stem}: cells north-south");

        let extent = header.grid.extent();
        assert_eq!(
            extent.west_deg_east,
            bounds.western_longitude_deg(),
            "{stem}: western boundary"
        );
        assert_eq!(
            extent.east_deg_east,
            bounds.eastern_longitude_deg(),
            "{stem}: eastern boundary"
        );
        assert_eq!(
            extent.south_deg_north,
            bounds.southern_latitude_deg(),
            "{stem}: southern boundary"
        );
        assert_eq!(
            extent.north_deg_north,
            bounds.northern_latitude_deg(),
            "{stem}: northern boundary"
        );

        // The header states the ocean the run was integrated in; a run whose
        // header disagreed with its scenario is unusable however good its
        // frames are.
        assert_eq!(
            header.physical_params.mean_depth_m,
            params.mean_thermocline_depth_m(),
            "{stem}: H"
        );
        assert_eq!(
            header.physical_params.reduced_gravity_m_per_s2,
            params.reduced_gravity_m_per_s2(),
            "{stem}: g'"
        );
        assert_eq!(
            header.physical_params.beta_per_m_per_s,
            params.beta_per_m_per_s(),
            "{stem}: beta"
        );
        assert_eq!(
            header.physical_params.rayleigh_damping_per_s,
            params.rayleigh_damping_per_s(),
            "{stem}: r"
        );
        assert_eq!(
            header.physical_params.reference_density_kg_per_m3,
            params.reference_density_kg_per_m3(),
            "{stem}: rho_0"
        );

        assert_eq!(
            header.output.frame_count, TEST_FRAME_COUNT,
            "{stem}: frame count"
        );
        assert_eq!(
            header.output.interval_s,
            TEST_OUTPUT_EVERY_N_STEPS as f64 * EXAMPLE_DT_S,
            "{stem}: output interval"
        );
        // The header names the scenario the run came from — the config
        // file's own name — so `inspect` on the run says which of the examples
        // produced it.
        assert_eq!(header.scenario_description, stem, "{stem}: description");

        let frames = reader
            .collect::<Result<Vec<_>, _>>()
            .unwrap_or_else(|error| panic!("{stem}: a frame did not read back: {error}"));
        assert_eq!(
            frames.len() as u64,
            TEST_FRAME_COUNT,
            "{stem}: frames on disk"
        );

        let grid = header.grid.grid();
        for (index, frame) in frames.iter().enumerate() {
            // Saved steps are 0, N, 2N, …, so frame `k` is at `k` output
            // intervals of model time. Exact, not approximate: the times are
            // computed from the same step count either way.
            assert_eq!(
                frame.t_s(),
                index as f64 * header.output.interval_s,
                "{stem}: time of frame {index}"
            );
            // The header's own list, not every variable the format knows: a
            // frame carries what its run wrote, and these scenarios are the
            // linear core's.
            for spec in &header.variables {
                let variable = spec.variable;
                let (nx, ny) = grid.field_shape(variable.staggering());
                let field = frame
                    .field(variable)
                    .unwrap_or_else(|| panic!("{stem}: frame {index} is missing {variable:?}"));
                assert_eq!(
                    field.len(),
                    nx * ny,
                    "{stem}: shape of {variable:?} in frame {index}"
                );
            }
            assert_eq!(
                frame.sst_anomaly_k(),
                None,
                "{stem}: frame {index} of an uncoupled scenario has no SST anomaly to carry"
            );
        }
    }
}

#[test]
fn the_steady_trades_run_tilts_the_thermocline_up_to_the_east() {
    // The control scenario's alizés blow westward (τx < 0), so the surface
    // layer is driven west, water piles up against the western boundary and
    // the thermocline deepens there while it shoals in the east — the
    // east–west tilt of the equatorial Pacific (CONTEXT.md, *Warm pool*;
    // docs/planning/01-scientific-model.md). This is the sign of the response,
    // which is fixed by the sign of the stress; the magnitude after one day is
    // not what is being checked.
    let scratch = ScratchDir::new(TICKET, "tilt");
    let out = run_scenario_text(
        &shortened_example("steady-trades"),
        scratch.path(),
        "steady-trades",
    );

    let mut reader = RunReader::open(&out).expect("the run directory reads back");
    let grid = reader.header().grid.grid();
    let (nx, ny) = grid.field_shape(Staggering::CellCenter);
    let last = reader
        .by_ref()
        .last()
        .expect("the run wrote frames")
        .expect("the last frame reads back");

    let mut west_m = 0.0;
    let mut east_m = 0.0;
    for j in 0..ny {
        for i in 0..nx / 2 {
            west_m += last.h()[j * nx + i];
        }
        for i in nx / 2..nx {
            east_m += last.h()[j * nx + i];
        }
    }

    assert!(
        west_m > 0.0,
        "the trades should have deepened the thermocline in the western half; \
         the anomaly there sums to {west_m} m"
    );
    assert!(
        east_m < 0.0,
        "the trades should have shoaled the thermocline in the eastern half; \
         the anomaly there sums to {east_m} m"
    );
}

#[test]
fn the_same_scenario_twice_writes_the_same_bytes() {
    // CODING_STANDARDS.md § *Correctness and failure*: identical scenario in,
    // byte-identical output. Two runs of one config into two directories, then
    // both files compared byte for byte.
    let scratch = ScratchDir::new(TICKET, "deterministic");
    let source = shortened_example("wind-burst");

    let first_dir = scratch.path().join("first");
    let second_dir = scratch.path().join("second");
    fs::create_dir_all(&first_dir).expect("the scratch directory is writable");
    fs::create_dir_all(&second_dir).expect("the scratch directory is writable");

    let first = run_scenario_text(&source, &first_dir, "wind-burst");
    let second = run_scenario_text(&source, &second_dir, "wind-burst");

    for name in [HEADER_FILE_NAME, FRAME_FILE_NAME] {
        assert_eq!(
            fs::read(first.join(name)).expect("the first run wrote its files"),
            fs::read(second.join(name)).expect("the second run wrote its files"),
            "two runs of one scenario disagree about {name}"
        );
    }
}

#[test]
fn a_config_that_is_not_there_is_reported_rather_than_panicked() {
    // CODING_STANDARDS.md § *Correctness and failure*: invalid user input is a
    // `Result` all the way up, so what reaches the terminal names the file
    // rather than unwinding.
    let scratch = ScratchDir::new(TICKET, "missing");
    let missing = scratch.path().join("not-a-scenario.toml");

    let output = run_cli(&missing, &scratch.path().join("run"));

    assert!(!output.status.success(), "a missing config should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(&missing.display().to_string()),
        "the message should name the config that is not there; it said: {stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "a missing config is input, not a bug; it said: {stderr}"
    );
}

#[test]
fn an_unstable_timestep_is_refused_before_anything_is_written() {
    // The CFL bound for the example basin is 0.8·dx/c ≈ 1.48e4 s at
    // dx ≈ 55.6 km and c = 3 m/s (T-01.3), so a day-long step is far outside
    // it. The run must be refused rather than shortened
    // (CODING_STANDARDS.md § *No silent clamping*), and refused before the
    // header is written, since a run directory that exists should hold a run.
    let scratch = ScratchDir::new(TICKET, "unstable");
    let mut config = ScenarioConfig::from_toml(&example_source("steady-trades"))
        .expect("a shipped example is a scenario");
    config.run.dt_s = 86_400.0;
    config.run.total_steps = TEST_TOTAL_STEPS;
    config.run.output_every_n_steps = TEST_OUTPUT_EVERY_N_STEPS;
    let path = scratch.path().join("unstable.toml");
    fs::write(
        &path,
        config.to_toml().expect("the config writes back to TOML"),
    )
    .expect("the scratch directory is writable");
    let out = scratch.path().join("run");

    let output = run_cli(&path, &out);

    assert!(!output.status.success(), "an unstable timestep should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("panicked"),
        "an unstable timestep is input, not a bug; it said: {stderr}"
    );
    assert!(
        !out.join(HEADER_FILE_NAME).exists(),
        "a refused run should not have left a header behind"
    );
}
