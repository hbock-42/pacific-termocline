# termocline-format

The on-disk contract between the engine and the visualizer: a JSON header plus
a sequence of binary frames, defined once and shared by both sides.

Depends on neither simulation logic nor UI code — that independence is the
point (see [ADR-0004](../docs/planning/adr/0004-data-interchange-format.md)).
The format lands in Epic 05.
