# Scientific Model

This is the physical model the engine implements. It is the standard
"1.5-layer reduced-gravity shallow-water" formulation used throughout
tropical ocean/ENSO dynamics research (Cane & Sarachik; McCreary; Cane–Zebiak
1985; Battisti 1988). It is simple enough to be tractable and fast, and rich
enough to reproduce real phenomena: equatorial Kelvin waves, Rossby waves,
wind-driven thermocline tilt, and (with the Epic 12 SST extension) an
emergent ENSO oscillation.

## 1.5-layer reduced-gravity ocean

Model the upper ocean as a single active layer of warm water sitting on an
infinitely deep, motionless abyss. State variables, functions of horizontal
position (x, y) and time t:

- `h(x, y, t)` — thermocline depth **anomaly** from a mean depth `H`
  (so total thermocline depth is `H + h`). This is the primary quantity of
  interest.
- `u(x, y, t)`, `v(x, y, t)` — zonal and meridional current anomalies of the
  upper layer.

## Governing equations

On an equatorial beta-plane (`f = β·y`, `f = 0` at the equator, valid for the
Pacific basin's latitude range):

```
∂u/∂t − f·v = −g'·∂h/∂x + τx/(ρ₀·H) − r·u
∂v/∂t + f·u = −g'·∂h/∂y + τy/(ρ₀·H) − r·v
∂h/∂t + H·(∂u/∂x + ∂v/∂y) = −r·h
```

Where:
- `g'` — reduced gravity, `g·Δρ/ρ₀` (Δρ = density contrast across the
  thermocline). Sets the wave speeds.
- `τx, τy` — wind stress components (the trade-wind / alizés forcing; see
  below).
- `ρ₀` — reference seawater density.
- `r` — linear (Rayleigh) damping/friction coefficient, standing in for
  unresolved dissipation and mixing. Keeps the model numerically stable and
  gives waves a realistic decay.
- `β` — meridional gradient of the Coriolis parameter at the equator
  (`≈ 2.3×10⁻¹¹ m⁻¹s⁻¹`).

These are linear shallow-water equations forced by wind stress — the
textbook equatorial-wave / adjustment problem. Nonlinear advection terms are
explicitly **out of scope for v1** (Epic 02 ships the linear model; a
nonlinear-advection ticket is a clearly separated, optional follow-up once the
linear core is validated).

## Wind forcing (the alizés)

`τx(x, y, t)`, `τy(x, y, t)` is an externally specified field, not something
the engine derives from an atmosphere model in v1. Three forcing scenarios
must be supported from the start (Epic 03):

1. **Steady trade winds** — spatially uniform (or simple x/y profile)
   easterly stress along the equatorial band, `τx < 0`, decaying off the
   equator. This is the control/equilibrium case: it should produce a
   thermocline that's deep in the west and shallow in the east (the observed
   mean state).
2. **Seasonal cycle** — the steady field modulated by an annual harmonic.
3. **Wind anomaly event** — a superimposed idealized westerly wind burst
   (Gaussian in x, y, and t) representing the kind of perturbation known to
   trigger El Niño onset.

Forcing is a pluggable function of `(x, y, t)`, not hardcoded, so new
scenarios can be added without touching the solver (see Epic 03).

## Domain and boundaries

- Basin approximating the equatorial Pacific: roughly 120°E–80°W in x,
  25°S–25°N in y (exact truncation is a config parameter, not fixed in the
  physics).
- **Meridional boundaries** (north/south edges): closed, no-normal-flow —
  physically justified since the mid-latitude subtropics aren't part of the
  reduced-gravity, equatorial-wave story we're after.
- **Western boundary** (~120°E, near Indonesia/New Guinea): closed,
  reflects incident Rossby energy partly back as Kelvin waves (this
  reflection is a well-known, physically important detail — Epic 04 covers
  it explicitly with its own validation test).
- **Eastern boundary** (~80°W, South America): closed, reflects incident
  Kelvin energy back as Rossby waves.

## Numerical scheme (summary — full detail in Epic 01)

- Spatial discretization: **Arakawa C-grid** finite differences (standard for
  shallow-water models — correctly represents geostrophic adjustment and
  wave dispersion at the grid scale).
- Time stepping: **leapfrog** with a **Robert–Asselin filter** to control the
  leapfrog computational mode, or RK4 as a simpler alternative — the choice
  is made and justified in Epic 01, not here.
- Explicit scheme ⇒ CFL condition ties max stable timestep to grid spacing
  and the fastest wave speed (`c = √(g'·H)`, the Kelvin wave speed); the
  engine must compute and enforce this rather than trusting user-supplied
  timesteps blindly.

## Validation targets (Epic 07)

Because this is meant to be a *scientific* simulation, correctness is judged
against known analytic/theoretical results, not just "it runs":
- Equatorial Kelvin wave phase speed `c = √(g'H)` and non-dispersive,
  eastward-only propagation on the equator.
- Equatorial Rossby wave dispersion relation and westward propagation,
  including the fact that Rossby waves travel at `c/3` for the gravest mode.
- Equatorial deformation radius `Le = √(c/β)` setting the meridional decay
  scale of both wave types.
- Steady-state thermocline tilt under constant easterly wind stress matching
  the analytic Sverdrup/Stommel-type balance.
- Energy/volume conservation in the undamped, unforced (`r = 0`, `τ = 0`)
  limit, to machine precision modulo numerical diffusion.

## Phase 2 (Epic 12, stretch): SST and the ENSO feedback loop

To get *emergent* ENSO rather than only prescribed-wind response, add a
mixed-layer SST anomaly equation (`T'`) coupled to `h` via
upwelling/entrainment, and feed `T'` back into `τx` (a simple statistical or
Gill-type atmosphere response — the Bjerknes feedback). This turns the model
from "ocean responds to wind" into "ocean and wind co-evolve," which is what
actually produces an oscillation. This is intentionally deferred: it depends
on the linear ocean core (Epics 01–07) being solid and validated first.
