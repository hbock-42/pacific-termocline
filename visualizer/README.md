# visualizer

Renders a run: basin maps, cross-sections, time series, playback.

Where the run comes from is the one thing that differs between the two targets.
Natively it is read — a directory, a pair of dropped files, or a URL — and the
visualizer never touches the physics. In a browser it is **computed**: per
[ADR-0012](../docs/planning/adr/0012-the-browser-runs-the-engine.md) the web
build links the engine, holds one of the scenarios in `scenarios/` and steps it
in the tab, because the run format is not served to the web at all. Both
origins end in the same `LoadedRun`, so every view below is the same code
either way.

A computed run is stepped between repaints, never in one call: `src/compute.rs`
steps until a wall-clock deadline — half of a 60 Hz frame, checked every eight
steps — so the tab keeps drawing and the run is watched as it develops. What it
may retain is capped: a run holds at most 33.6 MB of frames, checked against
the scenario's header before the first step, because with nothing downloaded it
is memory rather than bandwidth that a tab dies of. The browser scenarios are
the engine's coarsened to fit — 80 × 25 cells, 244 frames, 19.9 MB — and the
colour scale of a run still being computed covers the frames so far and widens
as it develops, which the shell says on screen rather than leaving to be
guessed.

Binary: `termocline-viz`. With *Compare two runs* ticked a second run opens
beside the first — another scenario computed, or another directory on the
command line. It draws the thermocline depth anomaly `h` of one chosen
frame as a colour map over the basin, with the wind stress `τ` that forced it
drawn over the top as arrows — or plays the run through. Time series and
cross-sections land in Epic 09.

The frame on screen is chosen with the scrubber above the map: drag it, step a
frame at a time with the arrow keys or the buttons beside it, a page of ten
with Page Up and Page Down, or jump to either end of the run with Home and End.
Any frame costs one decode wherever it sits in the run — the offset of each is
noted when the run is loaded (`src/run.rs`) — so dragging across a 731-frame
run does not get slower the further right it goes.

Playback is that same chooser driven by the clock instead of by a hand. ▶ Play
(or the space bar) starts it, the speed menu beside it picks how many frames a
second of real time it spends — a frame is a day under `steady-trades.toml`, so
thirty frames a second is a month a second — and it stops on the last frame
rather than looping. Pausing holds the frame it stopped on; on the last frame
the button reads ↻ Replay and starts the run over. A run opens paused, and a paused run asks for no
repaint, so nothing here costs anything until it is used. `src/playback.rs` has
the rest.

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

Two runs are drawn side by side on **one** colour scale — the one that covers
both of them — because two runs each on its own scale would look identical
however far apart they are: each would reach the ends of the same ramp. The
quieter run therefore occupies the middle of the scale, which is the true
reading, and both runs' own ranges are stated under the panels. The frame
index, the playback clock and the colour bar are shared too, so there is no
state in which the two panels could come to show different frames.

Two runs are refused rather than drawn where the side-by-side would claim
something they do not support: a shared grid — the same place on screen has to
be the same place in the ocean — or a shared meaning for the frame index, which
runs written at different cadences do not have. Everything else they may differ
in is what a comparison is for: the panels cover the frames both runs reach,
both draw `h`, which every run carries whether or not it couples SST, and what
else differs is stated beneath. `src/comparison.rs` has the rest.
