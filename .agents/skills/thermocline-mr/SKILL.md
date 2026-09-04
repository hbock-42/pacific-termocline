---
name: thermocline-mr
description: "Implement one ticket of the Pacific Thermocline backlog end to end: read the GitHub issue, write its acceptance criteria as failing tests, implement until the cargo gate is green, and open a pull request. Use when implementing a T-XX.Y ticket, when told to pick up an issue in this repo, or when dispatched into a worktree as a worker."
---

Implement exactly one ticket (`T-<epic>.<n>`), from a GitHub issue to an open
pull request. You are a **worker**: you finish at an open PR and never merge —
the coordinator owns a serial merge lane, so only one PR races `main` at a
time.

Vocabulary is in `CONTEXT.md`, rules in `CODING_STANDARDS.md`, and the standing
conventions in `AGENTS.md`. Read them before writing code, not after.

## 1. Resolve the ticket

`gh issue view <n> --comments` for the body: Description, Deliverable,
Acceptance criteria, Depends on.

Check the blockers are closed before starting:

```
gh api repos/hbock-42/pacific-termocline/issues/<n> --jq .issue_dependencies_summary
```

An open blocker means this ticket is not ready. Stop and report it rather
than working around the dependency — the graph is the schedule.

The issue is the spec. The epic files under `docs/planning/` are frozen
history and may contradict it; where they differ, the issue wins.

## 2. Turn acceptance criteria into failing tests

Each acceptance-criteria bullet becomes at least one test, written and
**failing** before the implementation exists. Use `/tdd` for the loop's rules —
seams, red before green, one vertical slice at a time.

Two rules from `CODING_STANDARDS.md` bind hardest here, because a scientific
suite that passes for the wrong reason is worse than no suite:

- Expected values come from an independent source — an analytic solution, a
  published result, a worked example. Never from running the code.
- Where a scheme has a known order of accuracy, assert the error *shrinks at
  that order* across at least two resolutions, rather than sitting under one
  fixed threshold.

Comment every tolerance with what it derives from.

## 3. Implement to green

Implement until the tests pass, then run the **gate** — the same three
commands CI runs, in this order:

```
cargo fmt --check
cargo clippy -- -D warnings
cargo test --workspace
```

All three green is the completion criterion. A PR opened before the gate is
green wastes a CI run and one of your two attempts.

## 4. Keep the goalposts where they are

When a test will not go green, the code is what changes. Never the test, never
the tolerance, never an ADR. The full rule is in `AGENTS.md` § *Never move the
goalposts* — read it if you feel the pull.

You get **two** attempts at a red build. If the second fails: stop, leave the
PR open, `gh issue edit <n> --add-label needs-human` with a comment stating
what you expected, what you measured, and where they diverged. Report failure
and finish. A stuck ticket surfaces as an open PR; it does not consume more
attempts.

Escalate the same way — without spending an attempt — the moment you find work
that contradicts an existing ADR, or a ticket whose acceptance criteria you
believe are wrong. Those are human decisions.

## 5. Review against the spec

Run `/code-review` with the fixed point `main`, and pass it the issue number so
its Spec axis reviews against the acceptance criteria. Fix what it finds, or
answer it in a PR comment saying why the finding does not apply. Its Standards
axis reads `CODING_STANDARDS.md`; treat a finding there as a defect, not a
suggestion.

## 6. Open the pull request

Branch `t-<epic>.<n>-<slug>`, one squashed-intent commit, then:

```
gh pr create --title "T-<epic>.<n> — <ticket title>" --body "Closes #<n>
<what changed, and any tolerance or design decision worth a reviewer's eye>"
```

Then **stop**. Do not merge, and do not enable auto-merge — handing a PR to
the merge lane while another is mid-rebase is exactly the race the lane
exists to prevent.

Report completion with the PR number and whether the gate was green. When
dispatched as an orca worker:

```
orca orchestration send --type worker_done --subject "T-<epic>.<n> PR #<pr>" \
  --body "<summary>" --task-id <t> --dispatch-id <d> \
  --outcome succeeded|failed --json
```
