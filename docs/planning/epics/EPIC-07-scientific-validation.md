# Epic 07 — Scientific Validation

## Goal
Formally verify the engine against the theoretical benchmarks listed in
[01-scientific-model.md](../01-scientific-model.md#validation-targets-epic-07).
This is what makes the project a *scientific* simulation rather than just
"code that produces plausible-looking numbers" — several of these checks
were already sanity-tested inline in earlier epics; this epic makes them
rigorous, documented, and permanently part of the test suite.

## Scope
Analytic-comparison test suite, a short written validation report.

## Out of scope
New physics/features — this epic tests what already exists.

---

### T-07.1: Kelvin wave speed and non-dispersion test
- **Description:** Initialize a Kelvin-wave-like disturbance on the
  equator, run forward, measure its propagation speed and confirm it
  matches `c = √(g'H)` and stays eastward-only and non-dispersive
  (unchanging shape) as theory predicts.
- **Deliverable:** Automated integration test with a documented tolerance
  and the theoretical derivation referenced in comments/docs.
- **Acceptance criteria:** Measured phase speed within the documented
  tolerance of `√(g'H)` across at least two different grid resolutions
  (demonstrating the error shrinks with resolution, not a fixed offset).
- **Depends on:** Epic 04 complete.

### T-07.2: Rossby wave dispersion and `c/3` propagation test
- **Description:** Initialize a Rossby-wave-like disturbance, confirm
  westward propagation at approximately `c/3` for the gravest meridional
  mode, matching the equatorial Rossby dispersion relation.
- **Deliverable:** Automated integration test.
- **Acceptance criteria:** Analogous to T-07.1 for the Rossby case.
- **Depends on:** Epic 04 complete.

### T-07.3: Equatorial deformation radius test
- **Description:** Verify the meridional decay scale of both wave types
  matches `Le = √(c/β)`.
- **Deliverable:** Automated test fitting the meridional profile of a
  Kelvin-wave solution and comparing the decay scale to `Le`.
- **Acceptance criteria:** Fitted decay scale within documented tolerance of
  the theoretical `Le`.
- **Depends on:** T-07.1.

### T-07.4: Steady-state wind-driven tilt test
- **Description:** Run the steady trade-wind scenario (T-03.1) to
  equilibrium and compare the resulting `h(x)` tilt across the basin to the
  analytic Sverdrup/Stommel-type balance for that forcing.
- **Deliverable:** Automated test comparing equilibrium `h` profile to the
  analytic prediction.
- **Acceptance criteria:** Match within documented tolerance; test fails
  loudly (not silently passes) if equilibrium isn't reached within the
  configured run length.
- **Depends on:** T-03.1, Epic 04 complete.

### T-07.5: Conservation test (undamped, unforced limit)
- **Description:** Formalize the ad-hoc energy-conservation check from
  T-02.5 into a permanent, documented test: with `r=0` and `τ=0`, total
  energy should be conserved to within numerical-diffusion-scale tolerance
  over a long run.
- **Deliverable:** Automated long-run conservation test.
- **Acceptance criteria:** Energy drift over the test run stays within a
  documented, justified bound; bound is derived/explained, not just picked
  to make the test pass.
- **Depends on:** T-02.5.

### T-07.6: Validation report
- **Description:** A short written report (`docs/validation-report.md`)
  summarizing each test above: what was checked, the theoretical
  prediction, the measured result, and the tolerance/rationale. This is the
  document that answers "how do we know this simulation is scientifically
  correct."
- **Deliverable:** `docs/validation-report.md`.
- **Acceptance criteria:** Every test in this epic is represented in the
  report with its actual measured numbers from a real CI run, not
  placeholder text.
- **Depends on:** T-07.1 through T-07.5.
