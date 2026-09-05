//! T-11.1 — the getting-started guide walks a reader through commands and a
//! scenario that still work.
//!
//! `docs/getting-started.md` was verified by following it, which is what its
//! acceptance criterion asks for; these tests are what keeps it verified. A
//! guide is prose, and prose does not fail to compile when a scenario field is
//! renamed, a shipped scenario is retuned or a document it links to is moved.
//!
//! So the three things in it that can go stale silently are checked here: the
//! scenario it tells the reader to write is parsed and built by the engine
//! itself, the run figures it quotes are recomputed from the scenario files
//! rather than trusted, and every repository path it links to is opened.
//!
//! What is *not* checked here is the shape of a terminal session — the
//! progress lines, the timings, the byte counts on disk. Those come from
//! running the guide, and `engine/tests/run_progress.rs` and
//! `engine/tests/inspect.rs` are where the output formats themselves are
//! pinned.

use std::path::{Path, PathBuf};

use engine::Scenario;

/// The document under test.
const GUIDE: &str = include_str!("../../docs/getting-started.md");

/// The repository root, from this test's own location.
fn repository_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the engine crate sits inside the repository")
        .to_path_buf()
}

/// The fenced TOML blocks of the guide, in the order they appear — the
/// scenarios it asks the reader to write.
fn scenario_blocks() -> Vec<String> {
    let fence = "```toml";
    let mut blocks = Vec::new();
    let mut rest = GUIDE;
    while let Some(start) = rest.find(fence) {
        let body = &rest[start + fence.len()..];
        let body = body.strip_prefix('\n').unwrap_or(body);
        let (block, after) = body.split_once("```").expect("a fence that is closed");
        blocks.push(block.to_owned());
        rest = after;
    }
    blocks
}

#[test]
fn the_scenario_the_guide_tells_the_reader_to_write_is_one_the_engine_accepts() {
    let blocks = scenario_blocks();
    assert_eq!(
        blocks.len(),
        1,
        "docs/getting-started.md is expected to hold exactly one scenario, the one section 5 \
         tells the reader to write; it holds {}. Check every TOML block below against the \
         engine before changing this count.",
        blocks.len()
    );

    let scenario = Scenario::from_toml(&blocks[0]).unwrap_or_else(|error| {
        panic!("the scenario in docs/getting-started.md § 5 is not one the engine accepts: {error}")
    });

    // The figures the guide quotes for it, each derived from the file rather
    // than from a previous run: the basin of `CONTEXT.md` is 160 degrees of
    // longitude by 50 of latitude, so one-degree cells make 160 x 50; and a
    // frame every 24 steps of a 2160-step run is 2160/24 = 90 frames after the
    // one at step 0.
    let grid = scenario.basin().grid();
    assert_eq!(
        (grid.nx(), grid.ny()),
        (160, 50),
        "docs/getting-started.md § 5 says the scenario is a 160 x 50 grid"
    );
    assert_eq!(
        scenario.output_schedule().frame_count(),
        91,
        "docs/getting-started.md § 5 says the scenario writes 91 frames"
    );
}

#[test]
fn the_run_figures_the_guide_quotes_for_the_control_scenario_are_the_scenario_s_own() {
    let path = repository_root().join("engine/scenarios/steady-trades.toml");
    let scenario = Scenario::load(&path)
        .expect("the control scenario the guide runs in section 2 is a scenario");

    // The figures sections 2 and 3 quote, each derived from the scenario file
    // rather than from the session they were taken in: the file asks for
    // 17 520 steps and a frame every 24 of them, which is 17520/24 = 730
    // frames after the one at step 0; and its half-degree cells cut the
    // 160 degrees of longitude and 50 of latitude of the basin (CONTEXT.md)
    // into 320 x 100. Retuning the scenario without retaking the session is
    // what this catches.
    let grid = scenario.basin().grid();
    assert_eq!(
        (grid.nx(), grid.ny()),
        (320, 100),
        "docs/getting-started.md quotes `grid: 320 x 100 cells` for {}",
        path.display()
    );
    assert_eq!(
        scenario.output_schedule().total_steps(),
        17_520,
        "docs/getting-started.md quotes 17520 steps for {}",
        path.display()
    );
    assert_eq!(
        scenario.output_schedule().frame_count(),
        731,
        "docs/getting-started.md quotes 731 frames for {}",
        path.display()
    );
}

#[test]
fn every_repository_path_the_guide_links_to_exists() {
    let root = repository_root();
    // Links are written relative to `docs/`, which is where the guide lives.
    let docs = root.join("docs");

    let mut checked = 0;
    let mut rest = GUIDE;
    while let Some(open) = rest.find("](") {
        let after = &rest[open + 2..];
        let (target, remainder) = after.split_once(')').expect("a link that is closed");
        rest = remainder;

        // External links and in-page anchors are somebody else's to keep
        // alive; only paths into this repository are checked.
        if target.starts_with("http") || target.starts_with('#') {
            continue;
        }
        let target = target.split('#').next().unwrap_or(target);
        assert!(
            docs.join(target).exists(),
            "docs/getting-started.md links to `{target}`, which is not in the repository"
        );
        checked += 1;
    }

    assert!(
        checked >= 5,
        "expected the guide to link to the reference documents around it; found {checked} \
         repository links, which suggests the extraction above stopped matching"
    );
}
