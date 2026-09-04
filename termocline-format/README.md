# termocline-format

The on-disk contract between the engine and the visualizer: a JSON header plus
a sequence of binary frames, defined once and shared by both sides.

Depends on neither simulation logic nor UI code — that independence is the
point (see [ADR-0004](../docs/planning/adr/0004-data-interchange-format.md)).

- `RunHeader` — written once per run, as JSON: format version, grid and basin
  extent, physical parameters, scenario description, the variable list with
  units, and the output cadence. Self-describing on purpose, so a reader never
  guesses at the shape or meaning of the frames beside it.
- `Frame` — one saved timestep: model time plus `h`, `u`, `v`, `τx` and `τy`,
  each a flat row-major buffer sized by where that variable sits on the
  C-grid (`GridSpec::field_len`).

The writer lands in T-05.2 and the reader in T-05.3; this crate holds the
types and nothing else.
