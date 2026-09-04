# ADR-0003: Arakawa C-grid + RK4 time stepping

## Status
Accepted, revisitable during Epic 01

## Context
The governing equations (see [01-scientific-model.md](../01-scientific-model.md))
are the linear shallow-water equations on a beta-plane. Two independent
choices need making: how variables are arranged spatially (grid staggering)
and how the system is advanced in time.

## Decision

**Spatial: Arakawa C-grid.** `h` at cell centers, `u` at cell east/west
faces, `v` at cell north/south faces. This is the standard choice for
shallow-water/equatorial-wave models because it correctly represents
geostrophic adjustment and gives accurate wave dispersion at the grid scale
— the alternative (A-grid, everything co-located) is known to badly
mis-represent exactly the wave physics this project cares about most.

**Temporal: RK4 (classic 4th-order Runge-Kutta), not leapfrog.**
Leapfrog is the traditional choice in this literature and is cheaper per
step, but it requires a Robert-Asselin filter to suppress its computational
mode, which adds implicit numerical damping that has to be tuned and
justified. RK4 has no computational mode, is straightforward to implement
and to reason about for correctness/validation (Epic 07), and the extra
per-step cost (4 evaluations vs. 1) is affordable at the grid resolutions
this project targets. This can be revisited if profiling in Epic 10 shows
time-stepping cost dominating and a cheaper scheme is worth the added
complexity.

## Consequences
- The engine's core solver step is `state_{n+1} = rk4_step(state_n, dt,
  forcing_fn)`, operating on the C-grid layout — no computational mode to
  filter, but every RK4 stage must respect the C-grid staggering when
  evaluating spatial derivatives.
- Timestep `dt` is bounded by the CFL condition for the fastest wave speed
  `c = √(g'H)`; the engine computes and enforces a safe `dt` from the grid
  spacing rather than trusting a user-supplied one blindly (see the Epic 01 tickets).
