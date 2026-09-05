# visualizer

Reads the engine's output files and renders them: basin maps, cross-sections,
time series, playback. Never touches the physics and never links against the
engine's simulation code ([ADR-0001](../docs/planning/adr/0001-engine-visualizer-split.md)).

Binary: `termocline-viz`. It loads a run — dropped files, a `?run=` URL, or a
directory natively — and draws the thermocline depth anomaly `h` of one chosen
frame as a colour map over the basin. Time series, cross-sections and playback
land in Epic 09.

The frame on screen is chosen with the scrubber above the map: drag it, step a
frame at a time with the arrow keys or the buttons beside it, a page of ten
with Page Up and Page Down, or jump to either end of the run with Home and End.
Any frame costs one decode wherever it sits in the run — the offset of each is
noted when the run is loaded (`src/run.rs`) — so dragging across a 731-frame
run does not get slower the further right it goes.

The colour scale is ColorBrewer's 11-class `RdBu`, reversed, diverging about
zero: `h` is a signed anomaly, so zero is pinned to the neutral middle class
whatever the frame holds, and the scale reaches equally far either side of it.
Red is a deeper-than-average thermocline and blue a shallower one, which puts
the warm pool warm and the cold tongue cool. `src/heatmap.rs` has the rest of
the reasoning.
