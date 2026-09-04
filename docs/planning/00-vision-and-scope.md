# Vision & Scope

## What we're building

A physically-grounded simulation of the equatorial Pacific Ocean's
thermocline: the sharp temperature boundary separating the warm, wind-mixed
surface layer from the cold deep ocean. The thermocline's depth and tilt
across the Pacific basin are the central object of the simulation because
they are the mechanism behind ENSO (El Niño / La Niña): trade winds (alizés)
pile up warm water in the west, deepening the thermocline there and shoaling
it in the east; when the winds relax, that tilt collapses eastward as
equatorial Kelvin waves, surfacing warm water off South America.

The simulation should let a user:
- Run the ocean model forward in time under a chosen wind-forcing scenario
  (steady trade winds, seasonal cycle, an imposed westerly wind burst, etc.).
- Observe the thermocline depth field evolve across the basin.
- See the wind stress field that's driving it.
- Watch equatorial waves (Kelvin, Rossby) propagate and reflect off the
  western/eastern boundaries.
- (Stretch) Get emergent ENSO-like oscillations once SST feedback is coupled
  in (Bjerknes feedback), rather than only prescribed wind forcing.

## What this is not (initially)

- Not a full 3D general circulation model (GCM). We use a **reduced-gravity,
  shallow-water** approximation (see [01-scientific-model.md](01-scientific-model.md)) — the
  standard, well-validated simplification used in classic ENSO theory
  (Cane–Zebiak, Battisti, delayed-oscillator literature).
- Not assimilating real observational data initially. Idealized/analytic wind
  forcing first; real reanalysis wind products are a later epic if desired.
- Not real-time interactive (engine runs a scenario to completion or to a
  checkpoint; the visualizer plays back / steps through the result). Live
  streaming between engine and visualizer is a stretch goal, not a
  requirement for v1.

## Two components, one contract

- **Engine**: pure computation, Rust, no rendering, no UI. Input = a
  scenario config (grid size, physical parameters, wind-forcing description,
  run length). Output = a self-describing time series of the ocean state
  written to disk.
- **Visualizer**: reads that output and renders it. Never touches the
  physics. Could in principle be replaced or reimplemented without touching
  the engine, and vice versa — this boundary is load-bearing, see
  [ADR-0001](adr/0001-engine-visualizer-split.md).

## Process for this repo

Per explicit instruction: **no simulation or app code is written until the
epics and their tickets are specified.** Each epic in `docs/planning/epics/`
lists the tickets (merge-request-sized units of work) needed to deliver it, with
a description, deliverable, and acceptance criteria for each. Development
starts at Epic 00 and proceeds roughly in order, though later epics may be
reprioritized once the engine is minimally running.
