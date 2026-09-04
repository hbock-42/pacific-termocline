# Contributing

This repository is built by autonomous agents working one ticket at a time.
The workflow below is what they follow, and what a human contributor should
follow too.

## The unit of work is a ticket

Work is broken into `T-<epic>.<n>` tickets — merge-request-sized units, each
with a description, a deliverable, acceptance criteria, and dependencies.
**One ticket = one GitHub issue = one pull request = one squashed commit on
`main`.**

Tickets live as [GitHub issues](https://github.com/hbock-42/pacific-termocline/issues),
grouped into epics by milestone. The issue is authoritative. The epic files
under `docs/planning/epics/` are a frozen record of how the backlog was
originally specified and are not updated as work proceeds; where the two
differ, the issue wins.

Ordering comes from GitHub's native issue dependencies. A ticket whose
blockers are still open is not ready — pick one from the frontier instead:

```sh
gh api "repos/hbock-42/pacific-termocline/issues?state=open&per_page=100" \
  --jq '.[] | select(.issue_dependencies_summary.blocked_by == 0) | "\(.number)  \(.title)"'
```

## The gate

Before opening a pull request, run what CI runs:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

`main` is protected: pull request required, the `ci` check must pass, branches
must be up to date, and there are **no bypass actors** — nobody merges past a
red build, including the repository owner. See
[ADR-0005](docs/planning/adr/0005-autonomous-implementation-pipeline.md).

## Tests are not negotiable

Each acceptance-criteria bullet becomes at least one test, written failing
before the implementation exists. This is a scientific simulation: a suite
that is green because an assertion was weakened is worse than a red one,
because the error becomes invisible.

So when something will not pass, **the code is what changes** — never the
tolerance, never the test, never an ADR. If you cannot get green, stop and
open an issue labelled `needs-human` describing what you expected, what you
measured, and where they diverged. The full rule is in
[AGENTS.md](AGENTS.md#never-move-the-goalposts).

Expected values come from an independent source — an analytic solution, a
published result, a worked example — never from running the code and pasting
its output. Every tolerance carries a comment saying what it derives from.

## Before you write code

- [`CONTEXT.md`](CONTEXT.md) — the domain glossary. Use its terms and its
  symbols. `h` is a depth *anomaly*, not a total depth.
- [`CODING_STANDARDS.md`](CODING_STANDARDS.md) — units in names, `Result`
  versus panic, determinism, scope guards.
- [`docs/planning/adr/`](docs/planning/adr/) — the decisions already made.
  Contradicting one is an escalation, not a judgement call.

## Licence

Dual MIT / Apache-2.0. Contributions are accepted under the same terms.
