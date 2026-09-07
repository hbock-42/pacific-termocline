//! Acceptance tests for T-08.4 — the visualizer is published to GitHub Pages.
//!
//! Everything here is about one failure, the one the ticket names: GitHub
//! Pages serves a project site from a *subpath*, so a build made for `/`
//! produces a page whose wasm is fetched from the wrong place. It builds, the
//! workflow goes green, and the site is blank. Nothing in `cargo test` would
//! notice, because the mistake is in a YAML file and a document.
//!
//! So the deployment's three moving parts are pinned here:
//!
//! - the base path the site is built for, derived from the repository URL in
//!   `Cargo.toml` rather than written out again, so that the test disagrees
//!   with the workflow the day either one changes;
//! - the permissions and the actions `actions/deploy-pages` requires, which
//!   are the difference between a deploy and a 403;
//! - the links to the live site, which are what the ticket asks a reader to
//!   arrive at.
//!
//! What these cannot check is that the deployed page renders — that is a
//! browser looking at the site, and the pull request reports it.

use std::path::{Path, PathBuf};

/// The deploy workflow under test.
const WORKFLOW: &str = include_str!("../../.github/workflows/pages.yml");

/// The README, whose top the ticket asks the live site to appear at.
const README: &str = include_str!("../../README.md");

/// The guide, which the ticket asks to link the live site.
const GUIDE: &str = include_str!("../../docs/getting-started.md");

/// Where `trunk build` leaves the site, relative to the repository root — the
/// directory the workflow uploads, holding the wasm bundle and nothing else.
const DIST_DIR: &str = "visualizer/dist";

/// The workspace manifest, the independent source for where the site lives.
const WORKSPACE_MANIFEST: &str = include_str!("../../Cargo.toml");

/// The repository root, from this test's own location.
fn repository_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the visualizer crate sits inside the repository")
        .to_path_buf()
}

/// The `owner/name` of the repository, read off `repository = "…"` in the
/// workspace manifest.
///
/// The point of reading it rather than writing it down: every expectation
/// below is derived from this pair, so the base path, the site URL and the
/// workflow cannot drift apart from the repository they are for.
fn repository_slug() -> (String, String) {
    const PREFIX: &str = "repository = \"https://github.com/";
    let start = WORKSPACE_MANIFEST
        .find(PREFIX)
        .expect("the workspace manifest names the repository");
    let rest = &WORKSPACE_MANIFEST[start + PREFIX.len()..];
    let slug = rest
        .split('"')
        .next()
        .expect("the repository URL is a closed string");
    let (owner, name) = slug
        .split_once('/')
        .expect("a GitHub repository URL is owner/name");
    (owner.to_owned(), name.trim_end_matches('/').to_owned())
}

/// Where GitHub serves this repository's project site from: a subpath named
/// after the repository, which is the whole reason this file exists.
fn public_url() -> String {
    let (_, name) = repository_slug();
    format!("/{name}/")
}

/// The address of the live site.
fn site_url() -> String {
    let (owner, name) = repository_slug();
    format!("https://{owner}.github.io/{name}/")
}

/// The lines of the workflow that are inside the block introduced by `key` at
/// the given indentation — a mapping's body, by indentation, which is what
/// YAML nesting is.
fn block_under(key: &str, indent: usize) -> Vec<&'static str> {
    let opener = format!("{}{key}:", " ".repeat(indent));
    let mut lines = WORKFLOW.lines().skip_while(|line| **line != opener);
    lines.next().expect("the block's own line");
    lines
        .take_while(|line| {
            line.trim().is_empty()
                || line.starts_with(&" ".repeat(indent + 1))
                || line.trim_start().starts_with('#')
        })
        .collect()
}

#[test]
fn the_site_is_built_for_the_subpath_pages_serves_it_from() {
    let expected = format!("--public-url {}", public_url());
    assert!(
        WORKFLOW.contains(&expected),
        "the deploy workflow must build with `{expected}`: GitHub Pages serves this project \
         site from {}, and a build made for `/` asks for its wasm at the wrong path — a green \
         workflow and a blank page",
        public_url(),
    );
    assert!(
        WORKFLOW.contains("trunk build --release"),
        "the site is a release build; a debug wasm bundle is several times the size"
    );
}

#[test]
fn the_workflow_carries_the_permissions_deploy_pages_requires() {
    let permissions = block_under("permissions", 0).join("\n");
    for required in ["pages: write", "id-token: write"] {
        assert!(
            permissions.contains(required),
            "the workflow's permissions must include `{required}` — `actions/deploy-pages` \
             mints an OIDC token for the deployment and fails with a 403 without it; they are \
             currently:\n{permissions}"
        );
    }
}

#[test]
fn the_workflow_uploads_what_trunk_built_and_deploys_it() {
    let dist = format!("path: {DIST_DIR}");
    assert!(
        WORKFLOW.contains("actions/upload-pages-artifact@"),
        "the built site reaches the deployment as a Pages artifact"
    );
    assert!(
        WORKFLOW.contains(&dist),
        "the artifact is trunk's output directory and nothing else: per ADR-0012 the site is \
         the wasm bundle, with no run files to serve alongside it"
    );
    assert!(
        WORKFLOW.contains("actions/deploy-pages@"),
        "`actions/deploy-pages` is what publishes the artifact — the shape Pages' \
         `build_type: workflow` expects"
    );
    assert!(
        WORKFLOW.contains("name: github-pages"),
        "the deploying job runs in the `github-pages` environment, which is where the \
         deployment's URL and its branch policy live"
    );
}

#[test]
fn the_workflow_republishes_the_site_when_main_moves() {
    let triggers = block_under("on", 0).join("\n");
    assert!(
        triggers.contains("push:") && triggers.contains("main"),
        "the deliverable is a workflow that keeps the site current, so a push to `main` \
         publishes; the triggers are currently:\n{triggers}"
    );
    assert!(
        triggers.contains("workflow_dispatch:"),
        "a site can also need republishing without a commit — a failed deploy, a rerun — so \
         the workflow is manually dispatchable"
    );
}

#[test]
fn the_readme_offers_the_live_site_before_it_asks_for_a_toolchain() {
    let site = site_url();
    let link = README
        .find(&site)
        .unwrap_or_else(|| panic!("the README must link the live site at {site}"));
    let build = README
        .find("cargo build")
        .expect("the README tells a reader how to build the project");
    assert!(
        link < build,
        "the live site belongs above the fold — a visitor should meet `{site}` before the \
         first build instruction, because since ADR-0012 there is nothing to install and \
         nothing to download to see a run"
    );
    let opening = &README[..link];
    assert!(
        opening.lines().count() < 30,
        "the live-site link is {} lines into the README, which is past the opening: it was \
         asked for right after the description of the project",
        opening.lines().count()
    );
}

#[test]
fn the_getting_started_guide_links_the_live_site() {
    let site = site_url();
    assert!(
        GUIDE.contains(&site),
        "docs/getting-started.md must link the live site at {site}: the guide's first section \
         asks for a Rust toolchain, and the browser build asks for nothing at all"
    );
}

#[test]
fn the_workflow_the_tests_read_is_the_one_github_runs() {
    let path = repository_root().join(".github/workflows/pages.yml");
    assert!(
        path.is_file(),
        "the deploy workflow lives at {} — GitHub runs what is in .github/workflows, and a \
         test that read a file from anywhere else would pass against a site that never deploys",
        path.display()
    );
}
