# ADR-0002: Visualizer language — Rust (not Flutter/Dart) for v1

## Status
Accepted, revisitable. **Extended by
[ADR-0006](0006-web-visualizer.md)**: the visualizer must also run in a
browser. The language choice below stands — `eframe` compiles to
`wasm32` — but the "desktop scientific tool" framing and the `plotters`
suggestion are superseded there.

## Context
The engine must be Rust (hard requirement, performance-critical numerical
core). The visualizer was explicitly left open: Rust, or Dart/Flutter if a
different language is preferred. Since the two communicate only through
files (ADR-0001), the choice is genuinely free of technical constraints
either way.

## Options considered
1. **Rust**, native, using something like `egui`/`eframe` (immediate-mode
   UI, simple plotting) or `wgpu` directly for custom map rendering, plus
   `plotters` for time-series charts.
2. **Dart/Flutter**, cross-platform app UI, reading the engine's output
   files from disk (or a bundled copy) and rendering with `fl_chart` /
   custom `CustomPainter` widgets for the basin map.

## Decision
Rust, for v1.

Reasoning:
- Single toolchain and single language for the whole repo — one `cargo
  build` builds everything, one CI pipeline, no context-switching cost
  while the project is still being figured out.
- The visualization needs (2D scalar field over a lat/lon-like grid, vector
  field overlay for wind, line charts for time series, basic
  playback/scrubbing) are all well within reach of `egui` + `plotters`/
  custom immediate-mode drawing — no need for Flutter's mobile/app-store
  polish for what is, at least initially, a desktop scientific tool.
- Keeps the option in ADR-0001 (file-based contract) fully exercised without
  extra FFI/interop risk while the project is establishing its scientific
  correctness, which is the harder and more important problem right now.

## Consequences
- Visualizer crate(s) live in the same Cargo workspace as the engine
  (`visualizer/` alongside `engine/`), sharing only the data-format crate
  (ADR-0004), not simulation code.
- If a genuine need emerges later (e.g. wanting a polished mobile/tablet app
  for presenting results), the file-based contract from ADR-0001 means a
  Flutter rewrite of just the visualizer is a contained, non-disruptive
  change — the engine and data format are untouched either way. This ADR
  should be revisited then, not before.
