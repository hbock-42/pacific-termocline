//! Acceptance tests for T-06.3 — scenario validation and actionable errors.
//!
//! The criterion is that *each known-bad scenario, from a small table of
//! deliberately broken examples, fails immediately with a message that says
//! what's wrong and how to fix it, not a panic/stack trace*. So this file is
//! that table: [`BROKEN_SCENARIOS`] holds one deliberately broken example per
//! way a scenario can be wrong, each as a one-line mutation of a template that
//! is itself checked to be valid, and every one of them is driven through the
//! same four assertions —
//!
//! 1. it is refused as a `Result`, never a panic (CODING_STANDARDS.md
//!    § *Correctness and failure*);
//! 2. the message names what is wrong — the offending value, and the field or
//!    section it came from;
//! 3. the message says how to fix it — the bound, the substitute value, or the
//!    knob to turn;
//! 4. and none of it reads like a stack trace.
//!
//! *Immediately* is the other half of the ticket, and it is checked twice: at
//! the library boundary, where `Scenario::from_toml` refuses the file rather
//! than letting `Solver::new` refuse it later; and through the real binary,
//! where a refused scenario has to leave the output directory untouched rather
//! than a half-written run behind.
//!
//! # Where the expected numbers come from
//!
//! Nothing here is read back from the engine. Every bound an example sits the
//! wrong side of is computed in the comment above it from the published
//! formula — the CFL bound of `docs/planning/adr/0003-numerical-scheme.md`,
//! the rotation bound of ADR-0007, the projection of
//! `docs/scenario-config-reference.md` — so an example that stops being a
//! violation fails here rather than passing for a new reason.
//!
//! # Tolerances
//!
//! None. Every assertion is either a substring of a message or a boolean about
//! a file existing; the arithmetic in the comments only has to be right to the
//! order of magnitude that puts an example on the wrong side of a bound, and each
//! is over a factor of two clear of it.

use std::fs;
use std::path::Path;
use std::process::{self, Command};

use engine::{Scenario, ScenarioError};

mod common;

use common::ScratchDir;

/// This file's ticket, which labels the directories it leaves in the system
/// temp directory.
const TICKET: &str = "t063";

/// A valid scenario, as the template every broken example mutates one line of.
///
/// The default Pacific at half a degree — 320 × 100 cells of 55 597.54 m — a
/// one-hour step well inside both stability bounds, and ten days of it saved
/// once a day.
const VALID_TOML: &str = r#"
[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = -80.0
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 240
output_every_n_steps = 24

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
"#;

/// One deliberately broken example, and what its refusal has to say.
struct BrokenScenario {
    /// What this example gets wrong, as the failure message reports it.
    what: &'static str,
    /// The line of [`VALID_TOML`] to replace, and what to replace it with. A
    /// `from` that is not in the template is a broken example rather than a
    /// broken engine, and is caught by
    /// [`every_broken_example_really_is_a_mutation_of_the_template`].
    from: &'static str,
    to: &'static str,
    /// Substrings that say *what is wrong*: the offending value, and the field
    /// or section it came from.
    names: &'static [&'static str],
    /// Substrings that say *how to fix it*: the bound, the substitute value,
    /// or the knob to turn.
    remedy: &'static [&'static str],
}

/// The table the acceptance criterion asks for: one deliberately broken
/// example per way a scenario can be wrong, in the order the loader checks
/// them.
const BROKEN_SCENARIOS: &[BrokenScenario] = &[
    // --- The file itself. ---
    BrokenScenario {
        what: "a boundary that is not a number",
        from: "western_longitude_deg = 120.0",
        to: r#"western_longitude_deg = "120E""#,
        // TOML reports the line it could not read, which is what a reader
        // needs to find it in a file of forty.
        names: &["western_longitude_deg = \"120E\""],
        remedy: &["expected"],
    },
    BrokenScenario {
        what: "a misspelled key, which would otherwise run a scenario nobody asked for",
        from: "equatorial_zonal_stress_pa",
        to: "equatorial_zonal_stres_pa",
        names: &["equatorial_zonal_stres_pa"],
        // `deny_unknown_fields` makes serde list the keys that do exist, which
        // is the fix.
        remedy: &["expected", "equatorial_zonal_stress_pa"],
    },
    BrokenScenario {
        what: "a forcing that does not exist",
        from: r#"type = "steady_trade_winds""#,
        to: r#"type = "hurricane""#,
        names: &["hurricane"],
        remedy: &["expected", "steady_trade_winds"],
    },
    // --- `[basin]`. ---
    BrokenScenario {
        what: "a cell size that is not a size",
        from: "resolution_deg = 0.5",
        to: "resolution_deg = -0.5",
        names: &["resolution_deg is -0.5"],
        remedy: &["greater than 0"],
    },
    BrokenScenario {
        what: "a northern boundary south of the southern one",
        from: "northern_latitude_deg = 25.0",
        to: "northern_latitude_deg = -25.0",
        names: &["northern_latitude_deg is -25", "southern_latitude_deg -25"],
        remedy: &["swap"],
    },
    BrokenScenario {
        // 120°E to 120°E is a basin of zero width, not one wrapped around the
        // planet (docs/scenario-config-reference.md, `[basin]`).
        what: "a basin of no width at all",
        from: "eastern_longitude_deg = -80.0",
        to: "eastern_longitude_deg = 120.0",
        names: &["longitude", "0 degrees"],
        remedy: &["widen"],
    },
    BrokenScenario {
        // The zonal span is 160°, and 160 / 0.3 = 533.33… cells: refused
        // rather than rounded, because rounding it would run a basin nobody
        // asked for.
        what: "a span that is not a whole number of cells",
        from: "resolution_deg = 0.5",
        to: "resolution_deg = 0.3",
        names: &["longitude", "resolution_deg", "0.3"],
        remedy: &["divides"],
    },
    BrokenScenario {
        // 160° by 50° at 0.01° is 16 000 × 5 000 = 8 × 10⁷ cells. A run holds
        // the state, RK4's five stage buffers and two wind-stress fields
        // resident, which is 22 `f64` per cell — call it 24, so 192 bytes —
        // and 8 × 10⁷ × 192 B ≈ 14 GiB, several times any budget a laptop can
        // meet. The engine has to say so before it starts allocating, not
        // when the allocator refuses.
        what: "a resolution so fine the run cannot be held in memory",
        from: "resolution_deg = 0.5",
        to: "resolution_deg = 0.01",
        names: &["16000", "5000"],
        remedy: &["coarsen resolution_deg"],
    },
    // --- `[physics]`. ---
    BrokenScenario {
        what: "a reduced gravity of zero, which would collapse the wave speed",
        from: "reduced_gravity_m_per_s2 = 0.06",
        to: "reduced_gravity_m_per_s2 = 0.0",
        names: &["reduced_gravity_m_per_s2 is 0"],
        remedy: &["greater than 0"],
    },
    BrokenScenario {
        what: "a negative damping coefficient, which would amplify rather than damp",
        from: "rayleigh_damping_per_s = 1.0e-7",
        to: "rayleigh_damping_per_s = -1.0e-7",
        names: &["rayleigh_damping_per_s is -0.0000001"],
        remedy: &["at least 0"],
    },
    // --- `[run]`. ---
    BrokenScenario {
        what: "a timestep of zero",
        from: "dt_s = 3600.0",
        to: "dt_s = 0.0",
        names: &["dt_s is 0"],
        remedy: &["greater than 0"],
    },
    BrokenScenario {
        what: "an output cadence of zero, which is not a cadence",
        from: "output_every_n_steps = 24",
        to: "output_every_n_steps = 0",
        names: &["every_n_steps is 0"],
        remedy: &["at least 1"],
    },
    BrokenScenario {
        // The run is 240 steps long, so a cadence of 480 saves the initial
        // state and nothing else: the run would take every one of its steps
        // and write none of them. That is the "output interval sane relative
        // to run length" of the ticket.
        what: "an output cadence longer than the run",
        from: "output_every_n_steps = 24",
        to: "output_every_n_steps = 480",
        names: &["every_n_steps is 480", "240 steps long"],
        remedy: &["at most 240"],
    },
    BrokenScenario {
        // dx = dy = 0.5 · 111 195.08 = 55 597.54 m and c = √(0.06 · 150) =
        // 3.0 m/s, so κ_max = 2√2/dx and the gravity-wave bound is
        // 0.8 · 2√2 / (c · κ_max) = 0.8 · dx / c ≈ 14 826 s. A day-long step
        // is nearly six times past it.
        what: "a timestep past the gravity-wave CFL bound",
        from: "dt_s = 3600.0",
        to: "dt_s = 86400.0",
        names: &["86400", "CFL"],
        remedy: &["at most", "coarsen the grid"],
    },
    // --- `[[wind]]`. ---
    BrokenScenario {
        what: "westerly trade winds",
        from: "equatorial_zonal_stress_pa = -0.05",
        to: "equatorial_zonal_stress_pa = 0.05",
        names: &["is 0.05 Pa"],
        remedy: &["negative"],
    },
];

/// A scenario that clears every bound above and is refused only by the
/// rotation bound of ADR-0007.
///
/// The template cannot show this one: at half a degree the gravity-wave bound
/// (≈ 14 826 s) is the tighter of the two, so any step past the rotation bound
/// is past the CFL bound first and the CFL message is the one that arrives.
/// Reaching to 60° of latitude in cells of 2° swaps them:
///
/// ```text
/// dx = 2 · 111 195.08 = 222 390.16 m,  so the CFL bound is 0.8·dx/c ≈ 59 304 s
/// |f|max = β · 60 · 111 195.08 = 2.3e-11 · 6.6717e6 = 1.5345e-4 s⁻¹
/// rotation bound = 0.8 · 2√2 / |f|max ≈ 14 746 s
/// ```
///
/// so `dt_s = 30 000 s` is comfortably inside the CFL bound and twice the
/// rotation bound. The basin still divides: 160° and 120° are 80 and 60 cells
/// of 2°.
const ROTATION_BOUND_TOML: &str = r#"
[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = -80.0
southern_latitude_deg = -60.0
northern_latitude_deg = 60.0
resolution_deg = 2.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 30000.0
total_steps = 240
output_every_n_steps = 24
"#;

// ---------------------------------------------------------------------------
// Criterion: every known-bad scenario is refused, by name, with a remedy.
// ---------------------------------------------------------------------------

#[test]
fn the_template_the_broken_examples_mutate_is_itself_valid() {
    // Otherwise every example below could pass for the wrong reason.
    Scenario::from_toml(VALID_TOML).expect("the template is a scenario the engine runs");
}

#[test]
fn every_broken_example_really_is_a_mutation_of_the_template() {
    for broken in BROKEN_SCENARIOS {
        assert!(
            VALID_TOML.contains(broken.from),
            "the example for {} replaces `{}`, which the template does not carry, so it would \
             be tested against an unmodified — and valid — scenario",
            broken.what,
            broken.from
        );
    }
}

#[test]
fn every_broken_example_is_refused_rather_than_panicking() {
    for broken in BROKEN_SCENARIOS {
        let error = refusal_of(broken);
        // The refusal is a value, so the caller chooses what to do with it;
        // reaching this line at all is the assertion.
        assert!(
            !error.to_string().is_empty(),
            "the example for {} was refused without saying anything",
            broken.what
        );
    }
}

#[test]
fn every_broken_example_names_what_is_wrong() {
    assert_every_refusal_carries(
        |broken| broken.names,
        "name the value or field that is wrong",
    );
}

#[test]
fn every_broken_example_says_how_to_fix_it() {
    assert_every_refusal_carries(|broken| broken.remedy, "say how to fix it");
}

#[test]
fn no_refusal_reads_like_a_stack_trace() {
    for broken in BROKEN_SCENARIOS {
        let message = refusal_of(broken).to_string();
        for leaked in ["panicked", "unwrap", "RUST_BACKTRACE", "src/"] {
            assert!(
                !message.contains(leaked),
                "the refusal of {} leaks `{leaked}` into a message a scenario author reads, \
                 got: {message}",
                broken.what
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Criterion: the refusal arrives up front — before the run, not partway
// through it.
// ---------------------------------------------------------------------------

#[test]
fn a_timestep_past_the_rotation_bound_is_refused_when_the_scenario_is_read() {
    // Both bounds on `dt_s` are part of reading the scenario, so one the
    // solver would refuse is one the loader has already refused: validation
    // is complete when `Scenario::from_toml` returns, and nothing downstream
    // gets to discover a new objection.
    let error = Scenario::from_toml(ROTATION_BOUND_TOML)
        .expect_err("a step past the rotation bound must be refused");
    let message = error.to_string();
    for expected in ["30000", "Coriolis"] {
        assert!(
            message.contains(expected),
            "the refusal should name `{expected}`, got: {message}"
        );
    }
    assert!(
        message.contains("closer to the equator"),
        "the refusal should say how to fix it, got: {message}"
    );
}

#[test]
fn the_rotation_bound_example_is_refused_only_by_rotation() {
    // Otherwise the test above would pass on the CFL bound instead, and the
    // rotation check could be deleted without anything noticing.
    let cfl_safe = ROTATION_BOUND_TOML.replace("dt_s = 30000.0", "dt_s = 14000.0");
    Scenario::from_toml(&cfl_safe)
        .expect("14 000 s is inside both the CFL bound (≈ 59 304 s) and the rotation one");
}

#[test]
fn a_basin_the_engine_can_hold_is_accepted() {
    // The memory budget is a bound, not a ban on refinement: 0.05° over the
    // Pacific is 3200 × 1000 = 3.2 × 10⁶ cells, which at the same 192 bytes a
    // cell is ≈ 586 MiB — a tenth of the cells the refused 0.01° example asks
    // for, times a hundredth.
    // At 0.05° the cells are ten times shorter, so the gravity-wave bound is
    // ten times shorter too — 0.8 · 5559.75 / 3.0 ≈ 1483 s — and the step has
    // to come down with the grid or the CFL bound would be what refused this.
    let refined = VALID_TOML
        .replace("resolution_deg = 0.5", "resolution_deg = 0.05")
        .replace("dt_s = 3600.0", "dt_s = 1000.0");
    Scenario::from_toml(&refined).expect("3.2 million cells is a run this build will start");
}

#[test]
fn a_refused_scenario_leaves_no_half_written_run_behind() {
    // "Fail fast rather than partway through a long run": the output
    // directory is untouched, so a second attempt is not competing with the
    // wreckage of the first.
    let scratch = ScratchDir::new(TICKET, "no-half-written-run");
    let config = scratch.path().join("broken.toml");
    let out = scratch.path().join("run");
    fs::write(
        &config,
        VALID_TOML.replace("dt_s = 3600.0", "dt_s = 86400.0"),
    )
    .expect("the scratch directory is writable");

    let output = run_cli(&config, &out);

    assert!(
        !output.status.success(),
        "a scenario the engine refuses must fail the command"
    );
    assert!(
        !out.exists(),
        "the run directory was created for a scenario the engine refused"
    );
}

#[test]
fn the_binary_reports_a_refused_scenario_on_stderr_rather_than_panicking() {
    let scratch = ScratchDir::new(TICKET, "no-stack-trace");
    let config = scratch.path().join("broken.toml");
    let out = scratch.path().join("run");
    fs::write(
        &config,
        VALID_TOML.replace("output_every_n_steps = 24", "output_every_n_steps = 480"),
    )
    .expect("the scratch directory is writable");

    let output = run_cli(&config, &out);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        !output.status.success(),
        "a scenario the engine refuses must fail the command, got stderr: {stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "the command panicked instead of reporting its input: {stderr}"
    );
    assert!(
        stderr.contains("every_n_steps") && stderr.contains("at most 240"),
        "stderr should carry the same actionable message the library returns, got: {stderr}"
    );
}

/// Assert that every refusal in the table carries the substrings `wanted`
/// picks out of its example, reporting `duty` when one is missing.
///
/// The two criteria — naming what is wrong and saying how to fix it — are the
/// same walk over the same table, differing only in which column they read, so
/// they are one helper rather than two loops that could drift apart.
fn assert_every_refusal_carries(
    wanted: fn(&BrokenScenario) -> &'static [&'static str],
    duty: &str,
) {
    for broken in BROKEN_SCENARIOS {
        let message = refusal_of(broken).to_string();
        for expected in wanted(broken) {
            assert!(
                message.contains(expected),
                "the refusal of {} should {duty}, by carrying `{expected}`; got: {message}",
                broken.what
            );
        }
    }
}

/// The refusal of one broken example, as the loader reports it.
fn refusal_of(broken: &BrokenScenario) -> ScenarioError {
    let toml = VALID_TOML.replace(broken.from, broken.to);
    Scenario::from_toml(&toml)
        .err()
        .unwrap_or_else(|| panic!("the engine accepted {}:\n{toml}", broken.what))
}

/// `termocline run` over `config`, writing into `out`.
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
