# ADR-0007: The timestep is bounded by rotation as well as by gravity waves

## Status
Accepted. Extends [ADR-0003](0003-numerical-scheme.md), which fixed RK4 and
stated the CFL bound in terms of the fastest wave speed alone.

## Context
ADR-0003 says the engine "computes and enforces a safe `dt` from the grid
spacing", and T-01.3 implemented exactly that: `max_stable_dt` keeps
`c·κ_max·dt` inside RK4's stability region, where `c = √(g'·H)` is the Kelvin
wave speed. `CFL_SAFETY_FACTOR` (0.8) was documented as absorbing "the Coriolis
term, Rayleigh damping and the wind forcing, which move the eigenvalues off the
pure imaginary axis" — treating rotation as a perturbation of a bound set by
the waves.

T-02.5 wired the Coriolis term and the wave terms through the integrator
together for the first time, and that assumption does not hold. Rotation is
not a perturbation of the wave oscillation; it is a second oscillation with its
own frequency. The momentum pair

```text
∂u/∂t = +f·v      ∂v/∂t = −f·u
```

has eigenvalues `±i·f` with `f = β·y`, so RK4 follows it only while
`|f|·dt ≤ 2√2`. That limit involves no wave speed and no cell spacing, so
refining the grid tightens the wave bound while leaving the rotation bound
untouched — and on a coarse basin reaching far from the equator the rotation
bound is the binding one. At 625 km cells over ±2500 km, the wave bound admits
a 32-hour step while the inertial period at the meridional walls is 30 hours,
and a run at that step amplifies the wall rows by a factor of 70 **per step**:
the output is not inaccurate, it is meaningless.

Three options were on the table:

1. **Shrink `CFL_SAFETY_FACTOR`** until the wave bound covers rotation too.
   Rejected: no fixed factor can, because the ratio of the two bounds depends
   on the basin's meridional extent and on `dx`, and a factor small enough for
   the worst case would slow every well-resolved run for nothing.
2. **Fold rotation into `max_stable_dt`.** Rejected: `termocline-numerics` is
   deliberately physics-free (it knows a "fastest signal speed", not a Coriolis
   parameter), and the basin's position on the beta-plane is physics.
3. **Check the second bound where both terms are visible.** Chosen.

## Decision

**A timestep must satisfy two independent bounds, and the engine's `Solver`
enforces both.** The gravity-wave CFL bound of T-01.3 stays exactly as it is,
in `termocline-numerics`. On top of it, `Solver::new` rejects any `dt` longer
than `CFL_SAFETY_FACTOR · 2√2 / max|f|`, where `max|f| = β·max|y|` is taken at
whichever meridional boundary of the basin lies further from the equator.

The same safety factor governs both: they are one stability region read against
two oscillations, and there is no reason to trust one closer to its boundary
than the other.

## Consequences
- `Solver::new` returns `Result`, and `SolverError` distinguishes the two
  refusals so a scenario's error message names the bound it violated. Neither
  bound ever shortens the timestep (CODING_STANDARDS.md § *No silent
  clamping*).
- The safe timestep is no longer a property of the grid and the wave speed
  alone: a basin that reaches further from the equator has a shorter one at
  the same resolution. Scenario authors trading basin width against timestep
  need to know this; `max_stable_dt` on its own no longer answers the
  question.
- `CFL_SAFETY_FACTOR`'s rationale in `termocline-numerics/src/cfl.rs` is
  corrected to point here rather than to claim it absorbs rotation.
- The bound lives in the engine, so `termocline-numerics` stays physics-free.
  If a future scheme needs it in more than one place, the shared piece to lift
  is "the longest step RK4 can take on an oscillation of a given rate", which
  both bounds already are.
