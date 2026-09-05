# Pacific Thermocline

A scientific simulation of the equatorial Pacific Ocean thermocline: wind-driven
(trade winds / alizés) upper-ocean dynamics, thermocline depth variability, and
the physical mechanisms behind El Niño–Southern Oscillation (ENSO).

The project is split into two independent components:

- **Engine** (`engine/`, Rust) — the physics simulation core. No visuals. Runs
  headless, produces time-series output of the ocean state.
- **Visualizer** (`visualizer/`, Rust) — consumes the engine's output and
  renders it (maps, cross-sections, time series, playback), per
  [ADR-0002](docs/planning/adr/0002-visualizer-language-choice.md). It is a
  browser app that also runs natively
  ([ADR-0006](docs/planning/adr/0006-web-visualizer.md)).
- **`termocline-format/`** — the file format the two share, and their only
  coupling ([ADR-0001](docs/planning/adr/0001-engine-visualizer-split.md)).
- **`termocline-grid/`** — the shared 2D field and Arakawa C-grid geometry
  types: indexing and staggering, no physics
  ([ADR-0003](docs/planning/adr/0003-numerical-scheme.md)).

The scientific model, architecture, and full backlog were specified up front
under `docs/planning/`, and the backlog now lives as
[GitHub issues](https://github.com/hbock-42/pacific-termocline/issues).
[`docs/getting-started.md`](docs/getting-started.md) walks from a clean
checkout to a simulation on screen.

```sh
cargo build --workspace
cargo test --workspace
cargo bench -p engine    # the performance suite; see docs/benchmarks.md
```

## Running the visualizer

Natively, optionally naming a run directory to open:

```sh
cargo run -p visualizer --bin termocline-viz -- /tmp/run-demo
```

In a browser, from `visualizer/` ([trunk](https://trunkrs.dev) builds the
`wasm32-unknown-unknown` target and serves it on `localhost:8080`):

```sh
cargo install trunk
cd visualizer && trunk serve
```

A browser has no filesystem, so a run reaches the web build by dragging its
`header.json` and `frames.bin` onto the page, or over HTTP: serve a run
directory and open `?run=<url>`, e.g. `http://localhost:8080/?run=run-demo/`
with the run copied into `visualizer/dist/run-demo/`.

## Where to start reading

1. [`docs/getting-started.md`](docs/getting-started.md) — running your first simulation: build, run a scenario, inspect it, open it in the visualizer.
2. [`docs/planning/00-vision-and-scope.md`](docs/planning/00-vision-and-scope.md) — what we're building and why.
3. [`docs/planning/01-scientific-model.md`](docs/planning/01-scientific-model.md) — the physics and equations being simulated.
4. [`docs/the-physics-explained.md`](docs/the-physics-explained.md) — the same physics in plain language, for readers who want to understand a run rather than derive the equations.
5. [`docs/validation-report.md`](docs/validation-report.md) — how we know the simulation is scientifically correct: each scientific test, its analytic prediction, the measured result and the derived tolerance.
6. [`docs/enso-oscillation-report.md`](docs/enso-oscillation-report.md) — whether the coupled model produces an ENSO-like oscillation, what delayed-oscillator theory predicted beforehand, and where the model met that prediction and where it did not.
7. [`docs/scenario-config-reference.md`](docs/scenario-config-reference.md) — every field of a scenario TOML file, its units and its valid range.
8. [`docs/benchmarks.md`](docs/benchmarks.md) — what the performance suite measures, how to run it, and how to read its figures.
9. [`docs/planning/adr/`](docs/planning/adr/) — key architecture decisions.
10. [`docs/planning/epics/`](docs/planning/epics/) — the full backlog, epic by epic, each broken into ticket-sized units of work. Frozen: the GitHub issues are authoritative.
10. [`CONTEXT.md`](CONTEXT.md) — the domain glossary, physics terms with their symbols.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Work proceeds one ticket at a time, and
`main` is gated on CI with no bypass — see
[ADR-0005](docs/planning/adr/0005-autonomous-implementation-pipeline.md).

## Licence

Dual-licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option.
