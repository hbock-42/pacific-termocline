# Epic 10 — Performance

## Goal
Make the engine fast enough to run realistic-resolution, multi-decade
scenarios in reasonable wall-clock time, since "must be optimized" was an
explicit requirement, and profile before assuming any specific optimization
is needed.

## Scope
Benchmarking, profiling, and targeted optimization of the engine's hot path
(the RK4/RHS evaluation loop).

## Out of scope
Rewriting the numerical scheme (that's an ADR-0003 revisit, not a
performance epic).

---

### T-10.1: Benchmark harness
- **Description:** `cargo bench` (via `criterion`) benchmarks for the core
  RHS evaluation (Epic 02) and a full short scenario run, at a couple of
  representative grid resolutions.
- **Deliverable:** Benchmark suite, run in CI as a non-blocking report (not
  a pass/fail gate initially).
- **Acceptance criteria:** Benchmarks run reproducibly and report
  timestep-per-second / grid-cells-per-second figures.
- **Depends on:** Epic 06 complete.

### T-10.2: Profile and identify hot paths
- **Description:** Profile a realistic-resolution run (e.g. with `perf` /
  `cargo flamegraph`) and document where time actually goes, before
  optimizing anything.
- **Deliverable:** A short findings note (`docs/performance-notes.md`)
  identifying the actual hot path(s).
- **Acceptance criteria:** Findings are backed by an attached/described
  flamegraph or profile output, not guesswork.
- **Depends on:** T-10.1.

### T-10.3: Parallelize the grid-update loop (`rayon`)
- **Description:** If T-10.2 confirms the per-cell RHS evaluation is the
  bottleneck (expected, given it's embarrassingly parallel across grid
  cells), parallelize it with `rayon`.
- **Deliverable:** Parallelized RHS evaluation.
- **Acceptance criteria:** Benchmark from T-10.1 shows meaningful speedup
  on a multi-core machine; existing validation tests (Epic 07) still pass
  bit-for-bit or within the same documented tolerance (parallelism must not
  change results beyond floating-point summation-order noise).
- **Depends on:** T-10.2.

### T-10.4: Revisit `f32` vs `f64` where profiling justifies it
- **Description:** Only if T-10.2 shows memory-bandwidth-bound behavior:
  evaluate switching field storage to `f32` for the bulk grid data while
  keeping `f64` accumulation where precision matters (e.g. long-run energy
  conservation, per T-07.5).
- **Deliverable:** Either the change (if justified) or a documented decision
  not to make it (if profiling doesn't support it).
- **Acceptance criteria:** If made, Epic 07 validation tests still pass
  within a re-justified tolerance; if not made, the reasoning is recorded.
- **Depends on:** T-10.2.
