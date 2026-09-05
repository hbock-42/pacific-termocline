# engine

The simulation core: pure computation, no rendering, no UI.

Input is a scenario (grid, physical parameters, wind forcing, run length);
output is a time series of the ocean state written through `termocline-format`.
Binary: `termocline`.

A scenario is a TOML file; every field it may carry is documented in
[`docs/scenario-config-reference.md`](../docs/scenario-config-reference.md), and
three worked examples live in [`scenarios/`](scenarios/).

A run reports its progress on stderr — percent complete, model time, wall time
and an ETA — redrawing one line on a terminal and emitting whole lines, free of
control characters, anywhere else. `--quiet` silences it for a scripted or CI
run; `--verbose` adds the run's per-frame log lines. stdout carries only the
run's one-line summary.

Physics lands in Epics 01–04, the CLI in Epic 06. See
[`docs/planning/01-scientific-model.md`](../docs/planning/01-scientific-model.md).
