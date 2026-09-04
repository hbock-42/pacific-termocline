# Epic 08 — Visualizer Foundation

## Goal
Stand up the visualizer as a real (if minimal) application that can open a
run produced by the engine and show *something*, establishing the app
skeleton before layering on the full feature set in Epic 09.

## Scope
App shell, run-loading, a single basic rendered view. Rust + `egui`/`eframe`
per [ADR-0002](../adr/0002-visualizer-language-choice.md).

## Out of scope
The full visualization feature set (Epic 09), engine↔visualizer live
streaming (explicitly deferred per ADR-0001).

---

### T-08.1: `visualizer` app shell
- **Description:** A minimal `eframe`/`egui` desktop app that opens (a
  window, a menu to pick a run directory) and can load a run via
  `termocline-format::RunReader` (T-05.3).
- **Deliverable:** `visualizer` binary that opens and loads a run, no
  rendering of the data yet beyond printing header info in the UI.
- **Acceptance criteria:** Opening one of the example runs from Epic 06
  shows correct header metadata (grid size, scenario name, frame count) in
  the UI.
- **Depends on:** T-05.3, T-06.1.

### T-08.2: Single-frame thermocline-depth heatmap
- **Description:** Render `h` for one selected frame as a 2D heatmap
  (color-mapped scalar field) over the basin, the most basic possible
  version of the core visualization.
- **Deliverable:** Heatmap rendering of `h` for a chosen frame index.
- **Acceptance criteria:** Visual smoke test — rendered heatmap for a known
  test run's equilibrium frame shows deeper thermocline (one color extreme)
  in the west vs. shallower in the east, matching T-07.4's known result.
- **Depends on:** T-08.1.

### T-08.3: Frame scrubber
- **Description:** A slider/timeline control to step through a run's frames
  and update the heatmap (T-08.2) accordingly — the minimum viable
  "playback."
- **Deliverable:** Scrubber UI control wired to the heatmap view.
- **Acceptance criteria:** Dragging the scrubber updates the displayed frame
  with no perceptible lag on a moderately-sized test run.
- **Depends on:** T-08.2.
