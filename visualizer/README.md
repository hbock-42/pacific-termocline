# visualizer

Reads the engine's output files and renders them: basin maps, cross-sections,
time series, playback. Never touches the physics and never links against the
engine's simulation code ([ADR-0001](../docs/planning/adr/0001-engine-visualizer-split.md)).

Binary: `termocline-viz`. It loads a run — dropped files, a `?run=` URL, or a
directory natively — and draws the thermocline depth anomaly `h` of one chosen
frame as a colour map over the basin. Time series, cross-sections and playback
land in Epic 09.

The colour scale is ColorBrewer's 11-class `RdBu`, reversed, diverging about
zero: `h` is a signed anomaly, so zero is pinned to the neutral middle class
whatever the frame holds, and the scale reaches equally far either side of it.
Red is a deeper-than-average thermocline and blue a shallower one, which puts
the warm pool warm and the cold tongue cool. `src/heatmap.rs` has the rest of
the reasoning.
