# ADR-0006: The visualizer must run in a browser

## Status
Accepted. Extends [ADR-0002](0002-visualizer-language-choice.md) (Rust for the
visualizer) and narrows a filesystem assumption in
[ADR-0001](0001-engine-visualizer-split.md).

## Context
ADR-0002 chose Rust with `egui`/`eframe` for the visualizer, describing it as
"a desktop scientific tool." A new requirement lands after that decision: the
visualizer must be **runnable on the web**, not only as a native binary.

The language choice survives this — `eframe` targets `wasm32-unknown-unknown`
and renders through `wgpu` (WebGPU where available, WebGL2 as fallback), so
the same Rust codebase serves both targets. What does not survive is the set
of assumptions the plan made about *how a run reaches the visualizer*.

## Decision

**The visualizer targets the browser as a first-class platform, built on
`eframe` + `wgpu`, compiled to `wasm32-unknown-unknown` and also runnable
natively from the same source.**

Three consequences follow, and they are the reason this ADR exists rather
than a one-line note:

1. **`wgpu` is the rendering backend**, not the `glow`/OpenGL path. It is the
   backend that reaches WebGPU, and the one whose native and web behaviour
   agree most closely.

2. **A browser has no filesystem.** ADR-0001's file-based contract holds — the
   engine still writes a header plus binary frames, and the visualizer still
   never links the physics — but "the visualizer opens a run directory" is a
   native-only affordance. On the web a run arrives by **user file selection
   or drag-and-drop**, or is **fetched over HTTP** from a served location.

3. **The reader must be source-agnostic.** `RunReader` is defined over a byte
   source (`impl Read + Seek`, or an in-memory buffer), with any
   path-taking constructor as a native-only convenience behind a feature gate.
   A reader written against `std::fs::File` would have to be rewritten for the
   web, so T-05.3 is specified this way from the start.

## Considered options

- **Native-only, as ADR-0002 assumed.** Rejected: the requirement is explicit.
- **A separate web visualizer (TypeScript/WebGL) alongside the native one.**
  Rejected: two renderers to keep in agreement, and the file-format contract
  would need a second implementation — exactly what ADR-0004 exists to avoid.
- **Web-only, dropping the native build.** Rejected: the native build stays
  useful for large local runs that would be painful to hand to a browser, and
  costs nothing while the source is shared.

## Consequences
- Epic 08's app shell is a browser app that also runs natively, not a desktop
  app. Its run-loading affordance is file selection or fetch, not a directory
  picker.
- Frame decimation and run size matter more than they would natively: a run
  crossing the network or sitting in browser memory has tighter limits than
  one memory-mapped from local disk. The lazy-iteration requirement in T-05.3
  is now load-bearing rather than merely prudent.
- CI should build the `wasm32-unknown-unknown` target once the visualizer
  exists, so the web target cannot silently rot while native stays green.
- `plotters` (floated in ADR-0002 for time-series charts) is no longer the
  obvious choice; `egui_plot` draws through the same `wgpu` surface and avoids
  a second rendering path. Not decided here — it belongs to Epic 09, and this
  note exists so that ticket knows the constraint.
