# engine

The simulation core: pure computation, no rendering, no UI.

Input is a scenario (grid, physical parameters, wind forcing, run length);
output is a time series of the ocean state written through `termocline-format`.
Binary: `termocline`.

Physics lands in Epics 01–04, the CLI in Epic 06. See
[`docs/planning/01-scientific-model.md`](../docs/planning/01-scientific-model.md).
