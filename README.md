# Pacific Thermocline

A scientific simulation of the equatorial Pacific Ocean thermocline: wind-driven
(trade winds / alizés) upper-ocean dynamics, thermocline depth variability, and
the physical mechanisms behind El Niño–Southern Oscillation (ENSO).

The project is split into two independent components:

- **Engine** (`engine/`, Rust) — the physics simulation core. No visuals. Runs
  headless, produces time-series output of the ocean state.
- **Visualizer** (`visualizer/`, Rust) — consumes the engine's output and
  renders it (maps, cross-sections, time series, playback), per
  [ADR-0002](docs/planning/adr/0002-visualizer-language-choice.md).
- **`termocline-format/`** — the file format the two share, and their only
  coupling ([ADR-0001](docs/planning/adr/0001-engine-visualizer-split.md)).
- **`termocline-grid/`** — the shared 2D field and Arakawa C-grid geometry
  types: indexing and staggering, no physics
  ([ADR-0003](docs/planning/adr/0003-numerical-scheme.md)).

The physics is not implemented yet: the crates are skeletons, and the equations
land in Epics 01–04. The scientific model, architecture, and full backlog were
specified up front under `docs/planning/`, and the backlog now lives as
[GitHub issues](https://github.com/hbock-42/pacific-termocline/issues).

```sh
cargo build --workspace
cargo test --workspace
```

## Where to start reading

1. [`docs/planning/00-vision-and-scope.md`](docs/planning/00-vision-and-scope.md) — what we're building and why.
2. [`docs/planning/01-scientific-model.md`](docs/planning/01-scientific-model.md) — the physics and equations being simulated.
3. [`docs/the-physics-explained.md`](docs/the-physics-explained.md) — the same physics in plain language, for readers who want to understand a run rather than derive the equations.
4. [`docs/scenario-config-reference.md`](docs/scenario-config-reference.md) — every field of a scenario TOML file, its units and its valid range.
5. [`docs/planning/adr/`](docs/planning/adr/) — key architecture decisions.
6. [`docs/planning/epics/`](docs/planning/epics/) — the full backlog, epic by epic, each broken into ticket-sized units of work. Frozen: the GitHub issues are authoritative.
7. [`CONTEXT.md`](CONTEXT.md) — the domain glossary, physics terms with their symbols.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Work proceeds one ticket at a time, and
`main` is gated on CI with no bypass — see
[ADR-0005](docs/planning/adr/0005-autonomous-implementation-pipeline.md).

## Licence

Dual-licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option.
