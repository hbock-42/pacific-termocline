# Pacific Thermocline

A scientific simulation of the equatorial Pacific Ocean thermocline: wind-driven
(trade winds / alizés) upper-ocean dynamics, thermocline depth variability, and
the physical mechanisms behind El Niño–Southern Oscillation (ENSO).

The project is split into two independent components:

- **Engine** (`engine/`, Rust) — the physics simulation core. No visuals. Runs
  headless, produces time-series output of the ocean state.
- **Visualizer** (`visualizer/`) — consumes the engine's output and renders it
  (maps, cross-sections, time series, playback). Language TBD per
  [ADR-0002](docs/planning/adr/0002-visualizer-language-choice.md).

Neither directory exists yet. **This repository is currently in the planning
phase.** Before any code is written, the scientific model, architecture, and
full backlog of epics/tickets are being specified under `docs/planning/`.

## Where to start reading

1. [`docs/planning/00-vision-and-scope.md`](docs/planning/00-vision-and-scope.md) — what we're building and why.
2. [`docs/planning/01-scientific-model.md`](docs/planning/01-scientific-model.md) — the physics and equations being simulated.
3. [`docs/planning/adr/`](docs/planning/adr/) — key architecture decisions.
4. [`docs/planning/epics/`](docs/planning/epics/) — the full backlog, epic by epic, each broken into ticket-sized units of work.
