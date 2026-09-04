# Pacific Thermocline — Agent Instructions

## Agent skills

### Issue tracker

Issues live as GitHub issues in `hbock-42/pacific-termocline` (via the `gh` CLI). See `docs/agents/issue-tracker.md`.

**GitHub issues are authoritative.** The epic files under `docs/planning/epics/`
are a frozen historical record of how the backlog was specified; they are not
updated as work proceeds. Read tickets from GitHub, never from the epic files.

### Triage labels

Default five-role vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `CONTEXT.md` at the root, ADRs at `docs/planning/adr/` (not the default `docs/adr/`). See `docs/agents/domain.md`.

### Project skills

- **`thermocline-mr`** — implement one ticket end to end (tests first, cargo
  gate, PR). The worker.
- **`thermocline-orchestrate`** — dispatch and supervise parallel workers via
  orca, and own the serial merge lane. The coordinator.

## Engineering conventions

### The CI gate

A ticket is done when CI is green on its pull request — not when the code
looks right. CI runs `cargo fmt --check`, `cargo clippy -- -D warnings`, and
`cargo test --workspace`, and a ruleset on `main` requires it. Run all three
locally before opening a PR; do not open a PR expecting CI to find the
problems for you.

`main` is protected with no bypass actors. Nothing merges without a green PR,
including changes made by the repo owner.

### Never move the goalposts

An agent that cannot get CI green **stops and escalates**. It does not make
the test agree with the code. Specifically, never do any of the following in
order to reach green:

- loosen a numeric tolerance stated in a ticket's acceptance criteria
- delete a failing test, or mark it `#[ignore]` or `#[should_panic]`
- narrow an assertion's range, or drop a test case
- edit `docs/planning/01-scientific-model.md` or any ADR

Any of these ⇒ stop, open a GitHub issue labelled `needs-human` describing the
discrepancy and what you expected, and move on. **Tolerances change only by
human decision.** This codebase is a scientific simulation: a suite that is
green because the assertions were weakened is worse than a red one, because
the failure becomes invisible.

### Coding standards

See [CODING_STANDARDS.md](CODING_STANDARDS.md) — units, `Result`-vs-panic,
determinism, tolerance justification, scope guards.

### ADR discipline

Write an ADR (in `docs/planning/adr/`, next number in sequence) only when all
three hold:

1. **Hard to reverse** — changing your mind later has a real cost.
2. **Surprising without context** — a future reader will ask "why on earth did
   they do it this way?"
3. **A real trade-off** — there were genuine alternatives and one was chosen
   for specific reasons.

If any is missing, skip it. And if your work *contradicts* an existing ADR,
say so explicitly rather than quietly overriding it — that is an escalation,
not a judgement call to make alone.

### Vocabulary

Use the terms in [CONTEXT.md](CONTEXT.md) in issue titles, type and function
names, and test names. If you need a domain concept that isn't in the glossary,
that is a signal: either the project doesn't use that language (reconsider), or
there is a real gap worth recording.
