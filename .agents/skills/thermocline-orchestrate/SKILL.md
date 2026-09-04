---
name: thermocline-orchestrate
disable-model-invocation: true
description: "Run a wave of parallel ticket workers via orca and own the serial merge lane."
---

Drive the Pacific Thermocline backlog: dispatch up to **4** workers into orca
worktrees, supervise them, and merge their pull requests one at a time.

Two jobs, and the second is the one that keeps `main` honest:

- **Dispatch** — keep the ready frontier busy, never exceeding 4 in flight.
- **The merge lane** — merge exactly one PR at a time. GitHub's merge queue is
  unavailable on this repo (it needs an organization owner), so concurrent
  green PRs can each pass alone and still break `main` together. The lane is
  the substitute. See [ADR-0005](../../../docs/planning/adr/0005-autonomous-implementation-pipeline.md).

You are the only dispatcher: orca's nested worker depth is 1, so a worker
cannot fan out sub-workers even if it wants to.

## Before the first wave

Autonomous operation starts at `T-00.3`. `T-00.1`, `T-00.2` and `T-00.5` are
the supervised bootstrap — the ruleset cannot exist until CI has run once, and
until the ruleset exists nothing is gated. If `T-00.5` is open, stop and tell
the user; do not start a wave.

Confirm the gates are live:

```
gh api repos/hbock-42/pacific-termocline/rulesets
gh repo view hbock-42/pacific-termocline --json autoMergeAllowed,deleteBranchOnMerge
```

## Setting up the run

Orca tracks the DAG; it never schedules. You do the scheduling.

```
orca orchestration run-create --objective "Pacific Thermocline backlog" --json
orca orchestration task-create --spec "T-01.2 — <title> (issue #<n>)" --deps '<json>' --json
```

Build tasks from the open GitHub issues, taking each ticket's dependency edges
from its issue:

```
gh api repos/hbock-42/pacific-termocline/issues/<n>/dependencies/blocked_by
```

Keep dependency chains shallow — orca's own guidance is to avoid going deeper
than 3–4 steps in one run. Prefer one run per epic over one run for all 55
tickets.

## The wave loop

1. `orca orchestration task-list --ready --json` — the frontier. This is
   external memory: it survives your context being compacted, so trust it over
   your own recollection of what is in flight.
2. Dispatch up to 4 concurrent, preferring tickets whose files do not overlap:

```
orca orchestration worker-start --task <t> --worktree new-top-level \
  --name t-<epic>.<n> --agent claude --setup run --json
```

   The worker prompt names the ticket, its issue number, and the
   `thermocline-mr` skill, and tells it **not to send heartbeats** — send only
   `worker_done`, `escalation` or `question`. Every inbox message pokes the
   coordinator's terminal, and a heartbeat carries nothing `worker-list`
   cannot already tell you; a wave of heartbeating workers is pure noise in
   the human's input. `--name` becomes the worktree's branch, prefixed by
   the owner (`hbock-42/t-01-2-rk4-integrator`), so pass the ticket id in the
   name and let that branch stand — squash-merge makes it ephemeral anyway.

3. Wait on the mailbox — never a sleep/poll loop:

```
orca orchestration check --wait --types worker_done,escalation,question \
  --timeout-ms 900000 --json
```

   Acknowledge each delivery (`--ack <delivery_id>`), or the same batch
   replays. A timeout or `{count:0}` is a checkpoint, not a failure: tickets
   routinely run 15–60 minutes. Heartbeats and terminal activity mean alive,
   not done.
4. Answer a `question` with `orca orchestration reply`. Handle an `escalation`
   yourself or surface it to the user; an escalation is a decision a worker
   correctly refused to make alone.
5. On `worker_done`, take the merge lane below, then
   `orca orchestration worker-release --dispatch <d>` and refill to 4.

**Order matters in cleanup.** Wait for `worker_done`, then release, then remove
the worktree. Removing a worktree while its agent is still live kills the
worker before it can report, and the dispatch is recorded `failed` even though
the work merged fine — a lie in the record, and three consecutive failures
circuit-break a task that was never broken.

## The merge lane

One PR at a time, start to finish, before the next enters:

1. Rebase the PR branch onto current `main` and push.
2. Wait for CI: `gh pr checks <pr> --watch`.
3. Green ⇒ `gh pr merge <pr> --squash --delete-branch`.
4. Red after the rebase ⇒ the ticket's own gate passed but it conflicts with
   what landed since. Hand it back: re-dispatch that ticket once with the
   failure as context. That is a rebase failure, not one of the worker's two
   attempts.

Rebase in **your own clone**, never in the worker's worktree — its agent may
still be live, and two writers in one checkout is a race you will not enjoy
debugging:

```
git fetch origin && git checkout -B lane-rebase origin/<branch>
git rebase main
git push --force-with-lease origin lane-rebase:<branch>
```

`Cargo.lock` is the conflict you will actually hit, because every wave that
adds a dependency rewrites it. It is generated, so never hand-merge it: take
`main`'s copy (`git checkout --ours Cargo.lock`), finish the rebase, then
`cargo check --workspace` to regenerate it against the rebased tree and amend
it into the commit. Run the full gate locally before force-pushing.

The ruleset requires branches be up to date, so a stale PR cannot merge even
if you try. Let that be the backstop, not the plan.

## When a worker fails

A `worker_done` with `--outcome failed` means the worker spent both attempts.
**Do not re-dispatch it.** Confirm the issue carries `needs-human`, leave the
PR open, and treat everything downstream of that ticket as blocked. Orca
circuit-breaks a task after 3 consecutive failures; one dispatch per ticket
keeps you clear of that and keeps a stuck ticket cheap.

Report at the end of each wave: merged, failed, still blocked, and what the
next frontier looks like.
