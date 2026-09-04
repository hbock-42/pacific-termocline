# Epic 06 — Engine CLI & Scenario Runner

## Goal
Turn the engine into an actual runnable tool: a CLI that takes a scenario
config, runs the simulation with progress feedback, and writes output —
tying together Epics 01–05 into the finished engine.

## Scope
CLI argument parsing, run orchestration, progress/logging, error handling
for bad configs.

## Out of scope
Anything visual (Epic 08/09).

---

### T-06.1: `engine run` CLI command
- **Description:** `engine run --config scenario.toml --out runs/my-run/`
  loads the config (T-03.4), builds the grid/basin (T-04.1), constructs
  the wind forcing (Epic 03), runs the RK4 loop (Epic 02) with the CFL-safe
  timestep (T-01.3), and writes output via `RunWriter` (T-05.2).
- **Deliverable:** Working `engine` binary with this command, using `clap`
  (standard, well-supported Rust CLI parsing crate).
- **Acceptance criteria:** Running each of the 3 example scenario configs
  from T-03.4 end-to-end produces a valid, readable run directory.
- **Depends on:** T-05.2, T-04.2, T-03.4.

### T-06.2: Progress reporting and logging
- **Description:** Print run progress (e.g. simulated time / wall time,
  percent complete) and structured logs (using `tracing` or similar) at
  appropriate levels, so a multi-minute run isn't a silent black box.
- **Deliverable:** Progress output during `engine run`.
- **Acceptance criteria:** A run's progress output updates at a reasonable
  cadence and reports a sane ETA; `--quiet` flag suppresses it for
  scripted/CI use.
- **Depends on:** T-06.1.

### T-06.3: Config validation and actionable errors
- **Description:** Validate the scenario config up front (grid size
  sane, CFL-derived `dt` achievable, output interval sane relative to run
  length) and fail fast with a clear message rather than partway through a
  long run.
- **Deliverable:** Pre-flight validation step in `engine run`.
- **Acceptance criteria:** Each known-bad config (from a small table of
  deliberately broken examples) fails immediately with a message that says
  what's wrong and how to fix it, not a panic/stack trace.
- **Depends on:** T-06.1, T-01.3.

### T-06.4: `engine inspect` command
- **Description:** A lightweight companion command that prints a run's
  header metadata (grid size, params, scenario, frame count) to the
  terminal — useful for sanity-checking a run without needing the
  visualizer.
- **Deliverable:** `engine inspect --run runs/my-run/`.
- **Acceptance criteria:** Output matches the header written by T-05.2 for
  a known test run.
- **Depends on:** T-05.3.
