# Epic 09 — Visualization Features

## Goal
Build out the full visualization feature set on top of the Epic 08
foundation: wind overlay, playback, cross-sections, and time-series charts —
enough to actually explore and understand a run's scientific content.

## Scope
Feature-level additions to the visualizer. No new engine work.

---

### T-09.1: Wind stress vector overlay
- **Description:** Overlay arrows/glyphs for `τx, τy` on top of the `h`
  heatmap (T-08.2), so the alizés forcing driving the response is visible
  alongside its effect.
- **Deliverable:** Toggleable wind-vector overlay layer.
- **Acceptance criteria:** Overlay correctly shows easterly arrows along the
  equator for the steady trade-wind scenario; toggling it on/off doesn't
  affect the underlying heatmap.
- **Depends on:** T-08.2.

### T-09.2: Playback controls
- **Description:** Play/pause/speed controls driving the frame scrubber
  (T-08.3) automatically, so a run can be watched rather than only
  manually stepped through.
- **Deliverable:** Play/pause/speed UI, auto-advancing the frame index.
- **Acceptance criteria:** Playback advances at the selected speed and
  correctly stops at the last frame; pause/resume preserves position.
- **Depends on:** T-08.3.

### T-09.3: Equatorial cross-section view
- **Description:** A secondary view plotting `h(x)` along the equator
  (y=0) as a line chart, updating with the current frame — this is the
  clearest way to see the west-high/east-low tilt and Kelvin wave
  propagation.
- **Deliverable:** Cross-section chart view (e.g. via `plotters` or
  `egui_plot`), synced to the same frame index as the heatmap.
- **Acceptance criteria:** Cross-section for a known equilibrium frame
  visually matches the analytic tilt from T-07.4.
- **Depends on:** T-08.3.

### T-09.4: Point time-series view
- **Description:** Click a point on the basin map to plot `h(t)` at that
  location across the whole run — the classic "Niño-index-style" time
  series view.
- **Deliverable:** Click-to-select-point + time-series chart.
- **Acceptance criteria:** Selecting a point near the eastern boundary
  during a wind-burst-anomaly run (T-03.3) shows the expected delayed
  thermocline-shoaling/deepening signal arriving after the western
  perturbation.
- **Depends on:** T-09.3.

### T-09.5: Run comparison / side-by-side
- **Description:** Load two runs at once and view their heatmaps
  side-by-side with synced frame index — useful for comparing scenarios
  (e.g. steady winds vs. wind-burst anomaly).
- **Deliverable:** Two-run side-by-side view.
- **Acceptance criteria:** Both panels stay frame-synced when scrubbing or
  playing either one.
- **Depends on:** T-08.3, T-09.2.
