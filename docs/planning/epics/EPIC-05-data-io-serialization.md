# Epic 05 — Data I/O & Serialization

## Goal
Implement the `termocline-format` crate and engine-side writing per
[ADR-0004](../adr/0004-data-interchange-format.md), so simulation output can
be persisted and later read back by the visualizer (Epic 08) or by
validation tooling (Epic 07).

## Scope
The format crate itself, engine-side writer, and a reader usable by both the
visualizer and test/validation code.

## Out of scope
Any rendering (Epic 08/09).

---

### T-05.1: `termocline-format` crate — header + frame types
- **Description:** Define the run header (`format_version`, grid dimensions
  and physical extent, `PhysicalParams`, scenario description, variable
  list, units, timestep count/interval) and the per-timestep frame struct
  (`t: f64`, `h`, `u`, `v`, `τx`, `τy` fields), all `serde`-derivable per
  ADR-0004.
- **Deliverable:** `termocline-format` crate with these types and round-trip
  serialization tests (`bincode` for frames, JSON for header).
- **Acceptance criteria:** A header and a frame each survive a
  serialize→deserialize round trip bit-for-bit (for header) / value-for-value
  (for frame, floats).
- **Depends on:** T-00.1.

### T-05.2: Engine-side run writer
- **Description:** A `RunWriter` that the engine opens at the start of a run
  (writes the header once) and appends frames to at a configurable output
  interval (not necessarily every timestep — long runs need decimated
  output).
- **Deliverable:** `RunWriter` in the engine crate, using
  `termocline-format`.
- **Acceptance criteria:** A short test run produces a header file + frame
  file readable back via `termocline-format`'s reader, with the expected
  number of frames for the configured output interval.
- **Depends on:** T-05.1, T-02.5.

### T-05.3: Reader API
- **Description:** A `RunReader` in `termocline-format` that lazily iterates
  frames (don't require loading an entire long run into memory at once —
  runs at realistic resolution/duration could be large).
- **Deliverable:** `RunReader::open(path) -> impl Iterator<Item = Frame>` (or
  equivalent), plus a `RunReader::header()` accessor.
- **Acceptance criteria:** Reading a run produced by T-05.2 yields frames
  identical to what was written, in order; memory usage doesn't scale with
  total run length (verified by a test with a deliberately long run and a
  memory-bound assertion, or at minimum a design review note if a hard
  memory test proves impractical).
- **Depends on:** T-05.2.
