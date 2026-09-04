# Epic 01 — Numerical Core

## Goal
Build the generic numerical machinery the physics will run on: the C-grid
field representation with real spatial-derivative operators, the RK4
integrator, and CFL-based timestep selection — all independent of the
specific shallow-water equations, which land in Epic 02.

## Scope
Grid geometry, finite-difference operators, time integrator, CFL logic.

## Out of scope
The actual momentum/continuity equations (Epic 02), wind forcing (Epic 03),
boundary physics (Epic 04).

---

### T-01.1: C-grid spatial derivative operators
- **Description:** Implement the finite-difference operators needed on an
  Arakawa C-grid: `d/dx`, `d/dy` mapping center↔face values, and
  face-to-center / center-to-face interpolation, per ADR-0003.
- **Deliverable:** Operators in `termocline-grid` (or a new
  `termocline-numerics` crate if that separation reads cleaner once written).
- **Acceptance criteria:**
  - Unit tests against known analytic functions (e.g. derivative of
    `sin(kx)` matches `k·cos(kx)` to expected finite-difference truncation
    error, converging at the correct order as grid spacing shrinks).
  - Operators are generic over grid size, not hardcoded to a basin
    dimension.
- **Depends on:** T-00.4.

### T-01.2: Generic RK4 integrator
- **Description:** A generic `rk4_step<S>(state: S, dt: f64, rhs: impl
  Fn(&S, f64) -> S) -> S` (or equivalent trait-based design) usable for any
  state type — deliberately not coupled to the shallow-water state yet, so
  it can be unit tested against a trivial ODE independent of the ocean
  physics.
- **Deliverable:** RK4 integrator with its own unit tests.
- **Acceptance criteria:**
  - Tested against a known ODE with an analytic solution (e.g. exponential
    decay or simple harmonic oscillator); measured error shrinks at 4th
    order as `dt` is halved.
- **Depends on:** T-00.1.

### T-01.3: CFL-based timestep selection
- **Description:** Given grid spacing and the fastest wave speed
  `c = √(g'H)`, compute a safe `dt` and expose it so the engine (Epic 06)
  can refuse/clamp unsafe user-supplied timesteps rather than silently going
  unstable.
- **Deliverable:** `max_stable_dt(grid_spacing, wave_speed) -> f64` with a
  documented safety factor, plus a runtime check the engine calls before
  starting a run.
- **Acceptance criteria:**
  - Unit tests pin the formula and safety factor.
  - A deliberately-too-large `dt` produces a clear, actionable error rather
    than running and producing garbage.
- **Depends on:** T-01.2.
