# ADR-0005: Autonomous implementation pipeline (agents merge their own work)

## Status
Accepted

## Context
The backlog is ~55 fully-specified, merge-request-sized tickets with an
explicit dependency graph. The goal is throughput: implement them with
parallel agents and without a human reviewing each pull request. That
inverts the usual safety model — normally human review is the last gate
before `main`, and here there is none.

A future reader looking at a `main` branch where every commit was written,
reviewed and merged by an agent will reasonably ask what stopped bad work
from landing. This ADR is the answer.

## Decision

**Agents implement, open, and merge their own pull requests. CI and a
GitHub ruleset are the only gates, and they are enforced by the platform
rather than by agent judgement.**

Four mechanisms make that survivable:

1. **The gate is server-side.** A ruleset on `main` requires a pull request
   (0 approvals, so agents aren't deadlocked), requires the `ci` status check
   (`cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test
   --workspace`), requires branches to be up to date before merging
   (`strict`), and blocks force pushes and deletions. Crucially it has **no
   bypass actors** — the repo owner is an admin and the agents authenticate
   as that owner, so leaving admin bypass on would reduce every gate to a
   suggestion an agent could skip with `gh pr merge --admin`.

2. **Goalposts are immovable.** Agents may not loosen a documented tolerance,
   delete or ignore a failing test, narrow an assertion, or edit an ADR or the
   scientific model in order to reach green. Failure escalates to a
   `needs-human` issue instead. See AGENTS.md. This is the rule that matters
   most: for a scientific simulation, a green suite obtained by weakening
   assertions is strictly worse than a red one, because it hides the error.

3. **One PR races `main` at a time.** GitHub's merge queue is unavailable
   here (it requires an organization-owned repository; this is a personal
   account), so concurrent green PRs could each pass individually and still
   break `main` together. Instead workers stop at "PR open, CI green" and a
   coordinator agent owns a **serial merge lane**, rebasing and merging one
   PR at a time.

4. **Failure is bounded.** A worker gets two attempts at a red build, then
   reports failure, leaves the PR open, and labels the issue `needs-human`.
   The coordinator does not re-dispatch. A genuinely stuck ticket surfaces as
   an open PR rather than consuming agent runs against the same wall.

Parallelism is provided by orca worktrees (isolated checkouts, one agent
each), capped at 4 concurrent workers — a limit set by Rust build cost
(no shared `target/` between worktrees) and by the fact that the dependency
graph rarely offers more than four independent tickets at once.

## Considered options

- **Human review of every PR.** Rejected as the explicit project goal is
  speed without manual validation.
- **GitHub merge queue.** Unavailable without transferring the repository to
  an organization; the serial merge lane is the substitute.
- **Agents self-merging with no server-side ruleset**, relying on their own
  discipline to check CI. Rejected: agent discipline is exactly what fails
  under pressure to finish, and it leaves no audit trail.

## Consequences
- The repo owner is gated too: no direct pushes to `main`, no merging a red
  PR without deliberately editing the ruleset. This friction is the point.
- Bootstrapping has an ordering constraint. Required status checks match by
  check *name*, which does not exist until a workflow has run once, so a
  ruleset applied before the first CI run would block every PR including the
  agents'. Therefore `T-00.1` (workspace skeleton) and `T-00.2` (CI) run as a
  **supervised** smoke test, the ruleset is applied after CI has run once, and
  autonomous operation begins at `T-00.3`.
- If the repository is ever made private on a Free plan, rulesets and branch
  protection disappear entirely and every gate above degrades to agent
  self-discipline. Do not make this repository private without revisiting
  this ADR.
- Issue dependency edges (`blocked_by`) order the *work*, but GitHub does not
  enforce them against merges. Nothing at the platform level stops a PR
  closing a blocked issue; that ordering is the coordinator's responsibility.
