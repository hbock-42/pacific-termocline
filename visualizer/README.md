# visualizer

Reads the engine's output files and renders them: basin maps, cross-sections,
time series, playback. Never touches the physics and never links against the
engine's simulation code ([ADR-0001](../docs/planning/adr/0001-engine-visualizer-split.md)).

Binary: `termocline-viz`. Rendering lands in Epics 08–09.
