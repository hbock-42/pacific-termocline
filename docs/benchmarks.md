# Benchmarks

The engine's performance suite: what it measures, how to run it, and how to
read what it reports. Epic 10's later tickets — profiling, parallelisation, an
`f32` field layout — all change the hot path, and
[CODING_STANDARDS.md](../CODING_STANDARDS.md) § *Performance* says a
performance change cites a measurement. This is that measurement.

## Running them

```sh
cargo bench -p engine              # the full suite, about a minute
cargo bench -p engine -- --quick   # a single-pass estimate, about ten seconds
```

Criterion writes its reports to `target/criterion/`, including an HTML index at
`target/criterion/report/index.html`.

Criterion's throughput line reads `elem/s` whatever the element is, so each
benchmark prints one legend line before its results saying what its element is.
Read the `thrpt` figures against that line.

To compare a change against `main`, save a baseline on the *same machine*,
before and after:

```sh
git switch main   && cargo bench -p engine -- --save-baseline before
git switch -      && cargo bench -p engine -- --baseline before
```

Criterion then prints the change and whether it is statistically significant.
Comparing a figure from one machine against a figure from another says nothing
at all; neither does comparing against a run from a laptop on battery, or one
taken while a build was running.

## What is measured

Two benchmarks, each at two grid resolutions — the control scenario's own
0.5° basin (320 × 100 cells) and a 1.0° basin at a quarter of its cells
(160 × 50). Both run the physics of
[`engine/scenarios/steady-trades.toml`](../engine/scenarios/steady-trades.toml)
with only the resolution and the run length changed, so the figures are about
the simulation this project actually runs rather than about a benchmark-only
configuration. The workload definitions are in `engine/src/benchmark.rs`, and
`engine/tests/benchmark_workloads.rs` is what holds them to all of that.

### `rhs_evaluation` — **grid cells per second**

One `ShallowWaterRhs::evaluate` over a whole basin: the pressure gradient, the
surface stress, the Rayleigh damping and the continuity divergence. RK4 calls
it four times per timestep, so it is the hot path at its least diluted, and it
is where a change to the inner loop shows up first.

Criterion is told that one iteration processes `nx · ny` elements, so its
`thrpt` line is grid cells per second directly, and the two resolutions are
comparable to each other despite differing in size by a factor of four.

### `scenario_run` — **timesteps per second**

A whole short run through `run_scenario`, the same entry point the `run`
command uses: scenario build, solver construction, the time loop with its wind
re-sampled at each of RK4's four stages, and the run directory written. It is
240 steps, ten days of model time at the control scenario's one-hour timestep.

Criterion is told that one iteration takes 240 steps, so its `thrpt` line is
timesteps per second. Grid cells per second for a run is that figure times the
workload's cell count.

The run writes exactly two frames — the initial state and the final one — so
the filesystem is a constant term of the measurement rather than one that grows
with the run length.

## What is *not* measured

- **The timestep does not vary with the resolution.** `dt` is one hour at both
  grids, even though the CFL bound would let the coarser one take a longer
  step. The two cases differ in the grid alone, so a difference between them is
  attributable to the grid alone.
- **Nothing is measured across machines or across CI runs.** The `bench` job in
  CI is a report attached to a pull request, never a gate: a GitHub runner
  shares hardware with neighbours nobody controls, and a pass/fail threshold on
  a number like that would fail for reasons that have nothing to do with the
  change. Read a CI report against the other figures in the same run.
- **Accuracy is not a performance question.** Whether an optimisation changed
  the answer is Epic 07's validation suite, not this one.

## Reproducibility

Both workloads are deterministic by construction: the scenario is compiled into
the binary, the right-hand side is evaluated over a closed-form analytic state,
and the engine's runs are byte-reproducible
([CODING_STANDARDS.md](../CODING_STANDARDS.md) § *Correctness and failure*). Two
runs of a workload write identical output, and two evaluations produce identical
tendencies — so a difference between two measurements is the code changing,
not the input.

What is *not* controlled is the machine. On a quiet laptop the suite's
confidence intervals sit within about ±1% of the median, which is the
resolution of the instrument: a change smaller than that is not something these
benchmarks can see, and claiming it would be reading noise.
