# Epic 11 — Documentation

## Goal
Make the project usable and understandable by someone other than its
author: how to run a scenario, how to read the output, and the science
behind it in approachable form.

## Scope
User-facing docs only. Design docs (`docs/planning/`) already exist and
stay as-is; this epic is about docs for *using* the finished tool.

---

### T-11.1: "Running your first simulation" guide
- **Description:** Step-by-step walkthrough: install/build, pick or write a
  scenario config, run `engine run`, open the result in the visualizer.
- **Deliverable:** `docs/getting-started.md`.
- **Acceptance criteria:** Following the doc from a clean checkout produces
  a viewable run, verified by actually following it.
- **Depends on:** Epic 06, Epic 08 complete.

### T-11.2: Scenario config reference
- **Description:** Full reference of every `ScenarioConfig` field, valid
  ranges, and what each wind-forcing scenario's parameters mean physically.
- **Deliverable:** `docs/scenario-config-reference.md`.
- **Acceptance criteria:** Every field in the actual `ScenarioConfig` struct
  is documented; a CI check (doc-comment presence, or a small script diffing
  struct fields against the doc) keeps it from silently going stale.
- **Depends on:** T-03.4.

### T-11.3: "The physics, explained" primer
- **Description:** An approachable (less equation-dense than
  [01-scientific-model.md](../01-scientific-model.md)) explanation of what
  the thermocline is, why trade winds tilt it, and how that produces ENSO —
  aimed at a reader who wants to understand what they're looking at in the
  visualizer, not derive the equations themselves.
- **Deliverable:** `docs/the-physics-explained.md`.
- **Acceptance criteria:** Reviewed for accuracy against
  `01-scientific-model.md` (no contradictions), written for a non-specialist
  reader.
- **Depends on:** none (can be written any time after Epic 01–02 exist to
  reference).
