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

/// The `owner/name` of the repository, read off `workspace.package.repository`
/// in the workspace manifest.
///
/// The point of reading it rather than writing it down: every expectation
/// below is derived from this pair, so the base path, the site URL and the
/// workflow cannot drift apart from the repository they are for. Parsed
/// rather than matched by substring, so a manifest that is reformatted, or
/// that grows a second `repository` key, does not quietly change what the
/// tests are asserting against.
fn repository_slug() -> (String, String) {
    let manifest: toml::Table =
        toml::from_str(WORKSPACE_MANIFEST).expect("the workspace manifest is TOML");
    let url = manifest["workspace"]["package"]["repository"]
        .as_str()
        .expect("workspace.package.repository is a URL");
    let slug = url
        .strip_prefix("https://github.com/")
        .expect("the repository is on GitHub, which is where Pages serves from")
        .trim_end_matches('/');
    let (owner, name) = slug
        .split_once('/')
        .expect("a GitHub repository URL is owner/name");
    (owner.to_owned(), name.to_owned())
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

/// The body of the workflow's top-level `key:` mapping — its indented lines,
/// which is what YAML nesting is.
fn top_level_block(key: &str) -> String {
    let opener = format!("{key}:");
    let mut lines = WORKFLOW.lines().skip_while(|line| **line != opener);
    lines.next().expect("the block's own line");
    lines
        .take_while(|line| line.starts_with(' ') || line.trim().is_empty())
        .collect::<Vec<_>>()
        .join("\n")
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
    let permissions = top_level_block("permissions");
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
    let triggers = top_level_block("on");
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
    let first_section = README.find("\n## ").expect("the README has sections");
    let second_section = README[first_section + 1..]
        .find("\n## ")
        .map_or(README.len(), |offset| first_section + 1 + offset);
    assert!(
        link < second_section,
        "the live site belongs in the README's first section, right after what the project \
         is; it is currently further down, where a visitor meets the build instructions first"
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
fn the_workflow_checks_the_bundle_it_built_before_publishing_it() {
    assert!(
        WORKFLOW.contains("dist/index.html"),
        "the flag is checked against the artifact, not only against this file: the workflow \
         greps the built `index.html` for the base path, so a `trunk` that stopped honouring \
         `--public-url` fails the build rather than publishing a blank page"
    );
    assert!(
        WORKFLOW.contains("GITHUB_STEP_SUMMARY"),
        "the bundle's size is reported by the run that built it — a number the ticket asks \
         for, and one that goes stale the first time it is written down by hand"
    );
    assert!(
        WORKFLOW.contains("--fail") && WORKFLOW.contains("page_url"),
        "the deployed site is fetched after the deploy: `it built` is not `it loads`, and a \
         404 on the wasm is what a wrong base path looks like from outside"
    );
}
