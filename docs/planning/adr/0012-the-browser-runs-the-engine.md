# ADR-0012: The browser runs the engine

## Status
Accepted. **Supersedes [ADR-0001](0001-engine-visualizer-split.md) on the web
target**, and narrows it natively. ADR-0006 (the visualizer must run in a
browser) stands and is what makes this possible.

## Context
ADR-0001 was written before any code existed. It made the engine and the
visualizer separate programs communicating only through files, and said
plainly that the visualizer "never links against the engine's simulation code".
Its stated costs included "no live *watch the simulation run* experience".

Three of that decision's four justifications have since expired or inverted:

- **Language independence** — ADR-0002 chose Rust for the visualizer, and both
  now live in one Cargo workspace. There is no second toolchain to protect.
- **No FFI or shared-memory complexity** — there is no FFI between two Rust
  crates in the same workspace. The cost the ADR was avoiding does not exist.
- **The engine is portable already.** Its simulation core — `solver`,
  `shallow_water`, `coriolis`, `boundary`, `forcing`, `state`, `params` —
  contains no filesystem access at all. Only `scenario` (config loading),
  `run_writer` (output) and `benchmark` touch the disk.

What did not expire is **reproducibility**: a run written to disk is a static
artifact you can archive, diff and replay. That is worth keeping, and this ADR
keeps it.

What forced the question was deployment. Serving pre-computed runs from GitHub
Pages means shipping the runs: the control run is 941 MB, against a 100 MB file
cap and a ~1 GB site cap. Every fix within ADR-0001 — a smaller demo run, HTTP
range requests over the fixed-size frames — buys a page that plays back *one
canned run*. The interesting thing about this model is not any single run; it
is what happens when you change the wind, or the feedback strength, and watch
the ocean answer.

## Decision

**On the web, the visualizer links the engine and computes runs itself. The
file format is not served to the browser at all.**

- `engine` builds for `wasm32-unknown-unknown`: the filesystem-touching
  modules go behind a feature the browser build does not enable.
- The browser holds a `Scenario` and steps it with `Solver::step`, in chunks
  small enough to keep the UI responsive, rather than calling a
  run-to-completion entry point.
- Frames are produced into memory. Every existing view already consumes
  `LoadedRun` rather than bytes, so the heatmap, scrubber, playback, wind
  overlay, cross-section, time series and comparison views are unchanged —
  what changes is where a `LoadedRun` comes from.
- `?run=`, drag-and-drop and the directory picker become **native-only**. The
  browser has nothing to load because it has nothing to fetch.

Natively, nothing changes: `termocline run` still writes header and frames,
`termocline inspect` still reads them, and the format remains the archival
contract of [ADR-0004].

## Consequences

- **The download problem disappears.** The site ships wasm, not runs. First
  paint costs a compile and a few steps rather than a transfer.
- **Memory replaces bandwidth as the binding constraint.** A frame of the
  0.5° basin is 1.29 MB, so the 731-frame control run is 941 MB — far past what
  a tab should hold. The browser's scenario must therefore be coarser, shorter,
  or decimated, and the UI should say which. This is a scenario choice, not an
  architectural one, but it is a real limit and must not be discovered by a
  visitor's tab dying.
- **The main thread must not be blocked.** A full run is 17,520 steps; at
  browser speed that is minutes. Stepping happens in chunks, yielding between
  them, so the UI stays live and progress is visible — which is a better
  experience than a progress bar over a download, and is the "watch the
  simulation run" ADR-0001 traded away.
- **Reproducibility is unchanged where it mattered.** The validated runs, the
  validation report and the benchmarks all still come from the native engine
  writing files. Nothing scientific rests on the browser.
- **The engine's public surface grows.** Making it usable step-by-step from
  another crate is a wider contract than "run this scenario to completion",
  and it is now load-bearing for the visualizer as well as the CLI.
- **A run can no longer be shared by URL on the web.** A scenario can, which
  is smaller and arguably more useful, but it is a real loss: two people
  looking at "the same run" now means two people computing it, and they get
  the same answer only because the engine is deterministic. It is — that is a
  tested property (`CODING_STANDARDS.md`, and the bit-for-bit tests of T-10.5
  and T-12.2) — but it is worth stating that the guarantee now carries weight
  it did not carry before.

[ADR-0004]: 0004-data-interchange-format.md
