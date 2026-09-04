# Epic 04 — Boundaries and the Real Basin

## Goal
Replace the idealized test grid from Epic 02/03 with the real Pacific basin
geometry and physically correct boundary conditions, including the
Kelvin/Rossby reflection behavior called out specifically in the scientific
model doc.

## Scope
Basin geometry/truncation, closed-boundary (no-normal-flow) conditions on
all four edges, and validated wave reflection at the western and eastern
boundaries.

## Out of scope
Realistic coastline shapes (islands, Indonesian archipelago geometry) — the
basin is a rectangle for v1, per the scientific model doc's stated
simplification.

---

### T-04.1: Basin geometry as a config parameter
- **Description:** Make basin extent (lon/lat bounds) and resolution part of
  `ScenarioConfig` (from T-03.4) rather than a compile-time constant, with
  sensible Pacific-basin defaults (~120°E–80°W, ~25°S–25°N).
- **Deliverable:** Grid construction takes basin bounds from config.
- **Acceptance criteria:** Changing the config's basin bounds changes the
  resulting grid's physical extent, verified by a unit test checking
  physical (not just index) coordinates of grid corners.
- **Depends on:** T-03.4.

### T-04.2: No-normal-flow boundary conditions
- **Description:** Enforce `u = 0` at east/west boundary faces and `v = 0`
  at north/south boundary faces (closed basin), correctly on the C-grid
  staggering.
- **Deliverable:** Boundary-condition application step run each timestep (or
  each RK4 stage, whichever is correct — this is a design decision to
  resolve during implementation and document in code).
- **Acceptance criteria:** No flow ever appears at a boundary face
  regardless of forcing; unit test checks this holds after many steps under
  active wind forcing, not just at t=0.
- **Depends on:** T-04.1.

### T-04.3: Western boundary Kelvin→Rossby reflection validation
- **Description:** Construct a targeted test: initialize a Kelvin
  wave-like pulse near the western boundary, run forward, and confirm it
  reflects as Rossby wave energy propagating back eastward at approximately
  the theoretical `c/3` group speed for the gravest meridional mode.
- **Deliverable:** A dedicated integration test (not just a smoke test),
  documented with the theoretical prediction it's checked against.
- **Acceptance criteria:** Reflected signal's propagation speed matches
  theory within a documented, justified tolerance.
- **Depends on:** T-04.2.

### T-04.4: Eastern boundary Rossby→Kelvin reflection validation
- **Description:** Mirror of T-04.3 for the eastern boundary: an incident
  Rossby-wave-like pulse should reflect as an eastward-then-boundary-bound
  Kelvin wave.
- **Deliverable:** Integration test analogous to T-04.3.
- **Acceptance criteria:** Same as T-04.3, eastern boundary.
- **Depends on:** T-04.2.
