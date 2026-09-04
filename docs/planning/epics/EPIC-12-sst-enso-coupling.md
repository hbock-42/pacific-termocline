# Epic 12 — SST & ENSO Coupling (Stretch)

## Goal
Turn the model from "ocean responds to prescribed wind" into a coupled
ocean-atmosphere system that can produce *emergent* ENSO-like oscillations,
per the Phase 2 extension in
[01-scientific-model.md](../01-scientific-model.md#phase-2-epic-12-stretch-sst-and-the-enso-feedback-loop).

## Status
**Stretch goal — explicitly deferred until Epics 01–07 (the validated linear
ocean core) are complete.** Not scheduled as part of the initial delivery;
kept here so the eventual scope and shape are already thought through when
the time comes.

## Scope (when undertaken)
Mixed-layer SST anomaly equation, its coupling to `h` via
upwelling/entrainment, and a simple statistical or Gill-type atmospheric
wind-stress response to SST (the Bjerknes feedback loop).

---

### T-12.1: SST anomaly state and equation
- **Description:** Add `T'(x, y, t)` to the ocean state with an advection +
  entrainment/upwelling equation coupling it to `h` and the mean upwelling
  implied by the wind forcing.
- **Deliverable:** Extended state type, new RHS term, config option to
  enable/disable the coupling (so the validated linear-only mode from
  Epics 01–07 remains available and regression-tested).
- **Acceptance criteria:** With the coupling disabled, all Epic 07
  validation tests still pass unchanged (proves the extension is additive,
  not a rewrite of the validated core).
- **Depends on:** Epic 07 complete.

### T-12.2: Statistical/Gill-type atmospheric wind response
- **Description:** A `WindStress` implementation whose `τx` responds to the
  current `T'` field (rather than being purely prescribed), closing the
  feedback loop.
- **Deliverable:** New `WindStress` implementation, composable with existing
  ones per T-03.3's `CompositeWind` pattern.
- **Acceptance criteria:** Feedback strength is a config parameter; at zero
  feedback strength, behavior matches the uncoupled model exactly (same
  regression argument as T-12.1).
- **Depends on:** T-12.1.

### T-12.3: Emergent-oscillation validation
- **Description:** Run the fully coupled model over a multi-year period and
  check for a self-sustained oscillation in an equatorial-Pacific SST index,
  with a period in the observed ENSO range (roughly 2–7 years), across a
  documented parameter range.
- **Deliverable:** Integration test + a written note analogous to the Epic 07
  validation report, since "does it oscillate at the right period" is the
  headline scientific claim of this epic.
- **Acceptance criteria:** Oscillation appears and its period falls in the
  documented target range for the chosen default parameters; sensitivity to
  the feedback-strength parameter is characterized (e.g. oscillation
  disappears below some threshold, consistent with delayed-oscillator
  theory).
- **Depends on:** T-12.2.
