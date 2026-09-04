# Epic 02 — Ocean Dynamics (Shallow-Water Core)

## Goal
Implement the actual reduced-gravity shallow-water equations from
[01-scientific-model.md](../01-scientific-model.md) on top of the Epic 01
numerical core, on an idealized doubly-periodic or simple rectangular grid
with placeholder (non-physical) boundaries — real basin boundaries are Epic
04, so this epic can be validated in isolation against known wave solutions.

## Scope
The `u, v, h` state, the RHS of the three governing PDEs (Coriolis, pressure
gradient via `g'∇h`, linear friction `r`), wired through the RK4 integrator.

## Out of scope
Wind forcing beyond a constant/test stub (Epic 03 owns the real forcing
scenarios), realistic basin boundaries (Epic 04).

---

### T-02.1: Shallow-water state type
- **Description:** Define the `OceanState { h: Field2D<f64>, u: Field2D<f64>,
  v: Field2D<f64> }` type on the C-grid staggering from Epic 01, plus the
  fixed physical parameters (`g'`, `H`, `r`, `β`, `ρ₀`) as a `PhysicalParams`
  struct.
- **Deliverable:** Types + constructors (e.g. `OceanState::at_rest(grid)`).
- **Acceptance criteria:** Compiles, documented units on every field
  (SI throughout — this matters for later validation against analytic
  formulas).
- **Depends on:** T-01.1.

### T-02.2: Coriolis + beta-plane term
- **Description:** Implement `f = β·y` and the Coriolis terms (`-f·v` in the
  u-equation, `+f·u` in the v-equation), correctly interpolated between the
  staggered `u`/`v` grid points.
- **Deliverable:** `coriolis_term(state, params) -> (du, dv)` contribution.
- **Acceptance criteria:** Unit test checks `f(y=0) == 0` and correct sign
  change across the equator; interpolation verified against a hand-computed
  small-grid example.
- **Depends on:** T-02.1.

### T-02.3: Pressure gradient + continuity terms
- **Description:** Implement `-g'·∇h` in the momentum equations and
  `-H·(∂u/∂x + ∂v/∂y)` (the divergence term) in the `h` equation, using the
  Epic 01 derivative operators.
- **Deliverable:** RHS contributions wired into a single `shallow_water_rhs
  (state, params, wind_stress) -> OceanState` (tendency).
- **Acceptance criteria:** Unit test: a Gaussian bump in `h` with zero
  velocity produces the expected initial acceleration direction (outward
  pressure-gradient flow) on a small test grid.
- **Depends on:** T-02.1.

### T-02.4: Linear friction term
- **Description:** Add the `-r·u`, `-r·v`, `-r·h` damping terms.
- **Deliverable:** Folded into `shallow_water_rhs`.
- **Acceptance criteria:** With `τ = 0` and `r > 0`, an initial perturbation
  decays monotonically to the rest state; unit test checks decay rate
  matches `exp(-r·t)` analytically for a single-mode test case.
- **Depends on:** T-02.3.

### T-02.5: Wire the RHS through RK4 on a test grid
- **Description:** Integrate `shallow_water_rhs` with the Epic 01 RK4
  stepper into a `step(state, dt, params, wind_stress_fn) -> OceanState`,
  runnable on a small idealized grid with a stubbed (constant) wind stress
  for now.
- **Deliverable:** End-to-end single-timestep and multi-timestep runs on a
  test grid, no I/O yet (that's Epic 05/06).
- **Acceptance criteria:**
  - Undamped, unforced (`r=0, τ=0`) run conserves total energy to within
    numerical-diffusion-scale tolerance over many steps (early version of
    the Epic 07 conservation test, run here for fast dev feedback).
  - Multi-step run doesn't NaN/blow up on the CFL-safe `dt` from Epic 01.
- **Depends on:** T-02.2, T-02.4, T-01.3.
