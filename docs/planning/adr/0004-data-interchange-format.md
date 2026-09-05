# ADR-0004: Engine output format — versioned custom binary + JSON metadata

## Status
Accepted, revisitable. Amended by
[ADR-0010](0010-reading-runs-from-older-format-versions.md), which decides what
a reader does with a run from an older `format_version` — the question this
ADR's versioning left open.

## Context
The engine and visualizer only communicate through files (ADR-0001). The
format needs to: store a time series of 2D fields (`h`, `u`, `v`, `τx`,
`τy`) on a known grid, be readable efficiently from Rust without a heavy
dependency, and be self-describing enough that the visualizer never needs to
guess grid dimensions or units.

## Options considered
1. **NetCDF** — the standard in the earth-science community, but pulls in a
   C library dependency (`netcdf-c`) and adds real build/toolchain
   complexity across platforms for comparatively little benefit at this
   project's scale.
2. **Custom binary format**: a small JSON/TOML header (grid dimensions,
   physical parameters, units, variable list, timestep count) followed by a
   flat binary blob of `f64` (or `f32`) arrays, one per output timestep.
3. **Off-the-shelf serialization** (`bincode`/`postcard` over `serde`) of a
   plain Rust struct per timestep, plus a header.

## Decision
Option 3: a `serde`-based format — a JSON header (human-readable, easy to
inspect/diff, written once per run) plus a sequence of `bincode`-encoded
frames (one per saved timestep) in a companion binary file. Defined as its
own small crate, `termocline-format`, shared by both engine and visualizer,
so there is exactly one place the format is defined and both sides use the
same struct definitions — no hand-written parsing on either end.

## Consequences
- `termocline-format` crate has zero dependency on simulation logic or UI
  code — just data structures + serde. Both `engine` and `visualizer` depend
  on it (see Epic 00 for workspace layout).
- Human-readable header (JSON) makes debugging/inspecting a run's metadata
  trivial without needing the visualizer.
- Format is versioned from day one (a `format_version` field in the header)
  so later changes don't silently break old runs or require a "does this
  file work with this build" guessing game.
- NetCDF interop (for interchange with real oceanographic tools/data) is
  explicitly deferred — an optional future ticket, not a v1 requirement.
