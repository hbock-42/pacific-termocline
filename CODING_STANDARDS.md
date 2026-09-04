# Coding Standards

Project-specific rules for this codebase. Deliberately short: anything `cargo
fmt` or `cargo clippy -- -D warnings` already enforces is out of scope, and so
is generic Rust advice. These are the rules a reviewer — human or agent —
should actually check, because nothing else will.

Vocabulary comes from [CONTEXT.md](CONTEXT.md). Numerical-scheme decisions come
from [ADR-0003](docs/planning/adr/0003-numerical-scheme.md).

## Physical quantities

- **Units are part of the name or the doc comment.** `mean_depth_m`,
  `wave_speed_m_per_s`, or a `///` line stating the unit. A bare `depth` on a
  physical quantity is a defect.
- **`h` is an anomaly, never a total depth.** Any code or comment treating it
  as absolute depth is wrong (see CONTEXT.md). Total depth is `H + h`.
- **Physical constants are named `const`s with a source.** `β`, `ρ₀`, `g`, and
  any tuned coefficient get a named constant and a comment citing where the
  value comes from. Inline numeric literals for physical quantities are a
  defect.
- **`f64` throughout the solver.** `f32` appears only where Epic 10 profiling
  justifies it, and never mixed silently into a `f64` computation.

## Correctness and failure

- **Invalid user input returns `Result`; broken invariants panic.** A scenario
  with an unstable timestep, a malformed config, or an unreadable output file
  is a `Result` with an actionable message naming the offending value and the
  bound it violated. Panics are for conditions that mean the code itself is
  wrong.
- **No silent clamping.** If the engine adjusts something the user asked for
  (a timestep, a grid dimension), it says so explicitly rather than quietly
  substituting a safe value.
- **Runs are deterministic.** Identical scenario in, byte-identical output.
  No unseeded randomness, no iteration-order dependence in anything reaching
  the output file.
- **No `unsafe` without an ADR.**

## Scope guards

- **The v1 core is linear.** Nonlinear advection terms are out of scope until
  the linear model is validated (see `01-scientific-model.md`). Adding them
  opportunistically inside another ticket is a defect, however tempting.
- **The grid knows about staggering; the physics doesn't.** Use the named
  C-grid offset types from `termocline-grid` rather than raw index arithmetic
  with magic `+1`/`-1` offsets.
- **The format crate is the contract.** Every public item in
  `termocline-format` carries a doc comment, and no simulation or UI logic
  lives there.

## Tests

- **Every tolerance is justified in a comment.** State what the bound is
  derived from — truncation order, machine epsilon, an analytic result — not
  that it happens to make the test pass.
- **Expected values come from an independent source**: an analytic solution, a
  published result, or a worked example. Never from running the code and
  pasting the output.
- **Convergence over point checks.** Where a scheme has a known order of
  accuracy, assert the error *shrinks at that order* across at least two
  resolutions rather than sitting under a single fixed threshold.

## Performance

- **No allocation in the inner time-stepping loop.** Buffers are allocated
  once per run and reused across steps.
- **Optimize against a measurement.** Performance changes cite a profile or a
  benchmark, per Epic 10. Speculative micro-optimization is Speculative
  Generality with extra steps.
