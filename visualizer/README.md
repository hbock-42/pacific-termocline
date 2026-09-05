# visualizer

Reads the engine's output files and renders them: basin maps, cross-sections,
time series, playback. Never touches the physics and never links against the
engine's simulation code ([ADR-0001](../docs/planning/adr/0001-engine-visualizer-split.md)).

Binary: `termocline-viz`. It loads a run — dropped files, a `?run=` URL, or a
directory natively — and draws the thermocline depth anomaly `h` of one chosen
frame as a colour map over the basin, with the wind stress `τ` that forced it
drawn over the top as arrows. Time series, cross-sections and playback land in
Epic 09.

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

The wind overlay is a toggleable layer of arrows, one every twelfth cell, each
pointing the way the stress there pushes the ocean and as long as that stress
is against the strongest in the run. The lattice is anchored on the middle of
the basin rather than its corner, so a basin symmetric about the equator gets a
row of arrows along the equator — where the trades are strongest, and where the
response they drive lives. An arrow that would come out shorter than the cell
it sits on is not drawn: the trades fall off as a Gaussian, and a mark standing
for 10⁻¹⁷ Pa is a dot that reads as data.

Under the steady trades the arrows point west along the equator, which is why
the thermocline under them tilts. The layer is drawn over the map rather than
into it, so turning it off leaves the map untouched; `src/wind.rs` has the
rest.
