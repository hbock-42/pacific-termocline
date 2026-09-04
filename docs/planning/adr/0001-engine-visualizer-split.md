# ADR-0001: Split the engine and the visualizer into independent programs

## Status
Accepted

## Context
The user's explicit requirement: an engine with no visuals, and something
else able to use that engine to visualize it, possibly in a different
language. We need a firm contract between them so each can be built, tested,
and iterated on independently.

## Decision
The engine and visualizer are **separate binaries/crates communicating only
through a file-based data format** (see
[ADR-0004](0004-data-interchange-format.md)), not a shared in-process
library and not (in v1) a live IPC/network connection.

- The engine is a CLI: given a scenario config file, it runs the simulation
  to completion (or to a requested checkpoint) and writes output files.
- The visualizer is a separate program: given the engine's output files, it
  renders them. It never links against the engine's simulation code and
  never re-implements the physics.

## Consequences
- Clean reproducibility: a run's output is a static artifact you can replay,
  diff, or archive independently of engine versions.
- Enables the visualizer to be written in a different language/toolchain
  entirely (see ADR-0002) without any FFI or shared-memory complexity.
- Costs: no live "watch the simulation run" experience in v1 — you run the
  engine, then open the result. This is acceptable per the vision doc; live
  streaming is an explicit future epic (Epic 10 includes an optional
  streaming MR flagged as stretch) if the file-based workflow proves too
  slow to iterate with.
