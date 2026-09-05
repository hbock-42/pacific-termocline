# Performance notes

Where a timestep's time actually goes, measured before anything is optimised
(T-10.2). [`docs/benchmarks.md`](benchmarks.md) says *how fast* the engine is;
this note says *what it is spending that on*, so that the two tickets behind
it — a rayon-parallel inner loop (T-10.3) and an `f32` field layout
(T-10.4) — can cite a profile rather than an intuition, which
[CODING_STANDARDS.md](../CODING_STANDARDS.md) § *Performance* requires of them.

Nothing in the hot path changed for this ticket. The solver, the right-hand
side, the Coriolis term and the integrator are exactly as T-07.4 left them.

## The finding

**Re-sampling the wind stress is 71% of a timestep. The shallow-water
right-hand side — the loop everyone assumes is the hot path, and the one the
benchmark suite isolates — is 13%.** The engine evaluates the prescribed wind
at every C-grid face at every one of RK4's four stages: 257 680 evaluations per
step on the control basin, each a pair of virtual calls ending in a libm `exp`.
The stress those evaluations compute is a function of `y` alone and, for the
control scenario, constant in time, so all 257 680 of them carry 101 distinct
numbers — one per row of `τx`, plus the `τy` that is identically zero.

The consequence for Epic 10 is arithmetic. Parallelising or narrowing the
right-hand side attacks 13% of a step; even made free it leaves 87%.

## What was measured, and on what

| | |
|---|---|
| Machine | Apple M1 Pro, 10 cores, 32 GB, macOS 26.1, on mains power and otherwise idle |
| Toolchain | `rustc 1.90.0`; `--release` for the tables, the `profiling` profile (release plus debug info) for the sampled run |
| Workload | `BENCHMARK_WORKLOADS` — the control scenario [`engine/scenarios/steady-trades.toml`](../engine/scenarios/steady-trades.toml) at 320 × 100 (0.5°) and at 160 × 50 (1.0°), the same two grids `docs/benchmarks.md` reports |

The wall-clock figures below are this machine's. The *shares* are what the note
concludes from, because a share is a property of the program and a duration is
a property of the laptop.

As a cross-check on the setup, `cargo bench -p engine --bench scenario_run` on
this machine reports **474 steps/s** at 320 × 100 (2.11 ms per step) and
**1 865 steps/s** at 160 × 50 (536 µs). The instrumented steps below cost
2.095 ms and 533 µs — within 1% of the benchmark at both grids, which is what
says the clock between the phases is not paying for itself in the total.

The 30-second spin the sampler is pointed at runs slower than either: 464
steps/s alone, and 445 steps/s with `sample` attached. The first gap is a
laptop under sustained load; the second is the sampler's own stack walks, which
is also why 15 s of "one sample per millisecond" produced 11 343 samples rather
than 15 000. Neither moves the shares, which is all the sampled profile is read
for, but both are reasons not to read its steps/s as the engine's speed.

## Two instruments, and why both

A profile nobody can cross-check is a guess with decimal places. Two
independent instruments were used, and the note only claims what both of them
say.

**A phase decomposition** — `engine::profiling::StepProfiler`, driven by
`cargo run --release --example profile`. It steps the real `ShallowWaterRhs`,
`CoriolisTerm`, `NoNormalFlow` and `Rk4` with a clock between them. It is
exhaustive by construction: four phases are timed and the fifth, the RK4 stage
algebra, is the residual of the step, so the five shares sum to one and no work
can hide between them. `engine/tests/step_profile.rs` asserts that the
instrumented step reaches bit-identically the same state as
`Solver::step_forced_by`, which is what makes it a profile of the engine rather
than of a lookalike.

Its weakness is that it perturbs what it measures: 24 clock reads per step
(1.2 µs per step at the worst clock cost the example measured, against a step
of 2.1 ms — under 0.06%), and a compiler barrier at each phase boundary.

**A sampling profiler** — macOS `/usr/bin/sample` against
`cargo run --profile profiling --example profile -- spin 30`, which takes
uninstrumented steps through the real `Solver` and does nothing else. It reads
no clock and inserts no barrier; it interrupts the process and records the
stack. 11 343 samples over 15 s, kept verbatim at
[`docs/profiles/2026-09-05-m1-pro-320x100.sample`](profiles/2026-09-05-m1-pro-320x100.sample)
so that the excerpt below can be checked against the artefact rather than
taken on trust.

They agree:

| Phase | timed | sampled |
|---|---|---|
| wind stress sampling | 71.0% | 71.7% |
| shallow-water terms | 13.4% | 12.7% |
| coriolis | 6.6% | 6.5% |
| rk4 stage algebra | 9.0% | 9.0% |
| boundary condition | 0.1% | < 0.1% |

Within a percentage point on every row, from two instruments with different
failure modes. That is why the finding is stated as a conclusion rather than as
an observation.

## The sampled profile

The attached artefact's own summary, by self time (samples out of 11 343, Rust
names demangled and truncated here; the file has them mangled and in full):

```
Sort by top of stack, same collapsed (when >= 5):
    exp  (in libsystem_m.dylib)                                   2734   24.1%
    <CompositeWind as WindStress>::stress                         2160   19.0%
    <ScenarioWind as WindStress>::stress                          1943   17.1%
    ShallowWaterRhs::evaluate                                      924    8.1%
    <OceanState as StateVector>::add_scaled                        824    7.3%
    Solver::step_forced_by::{{closure}}   (solver.rs:376)          789    7.0%
    DYLD-STUB$$exp                                                 512    4.5%
    CoriolisTerm::add_to_tendency                                  316    2.8%
    CGridOperators::face_x_to_face_y                               218    1.9%
    CGridOperators::face_y_to_face_x                               202    1.8%
    _platform_memmove  (in libsystem_platform.dylib)               198    1.7%
    CGridOperators::ddx_center_to_face                             170    1.5%
    CGridOperators::ddx_face_to_center                             138    1.2%
    CGridOperators::ddy_face_to_center                             108    1.0%
    CGridOperators::ddy_center_to_face                             104    0.9%
```

Two rows are worth spelling out.

`solver.rs:376` is `stage_stress.sample(basin, wind, stage_t_s)`. Every one of
those 789 self-samples is on that line — none on line 377 (`rhs.evaluate`) or
378 (`coriolis.add_to_tendency`), whose costs appear as their own frames — so
the sampling loop of `write_component` inlined into the closure, and its time
belongs to the wind. Adding it to `exp`, the two `stress` implementations and
the `exp` call stub gives the 71.7% above.

`DYLD-STUB$$exp` is the lazy-binding stub in front of the libm call: the
`stress` implementation reaches `exp` through the dynamic linker rather than
through an inlined intrinsic, which is another 4.5% of the run spent on the
*approach* to a function it calls a quarter of a million times per step.

## The phase table

`cargo run --release --example profile`, 24 steps per grid after a quarter of a
second of warm-up:

```
== 160x50 cells, 24 steps, 532.737µs per step ==
   phase                       share    per step
   wind stress sampling        69.9%    372.432µs
   shallow-water terms         13.7%     73.015µs
   coriolis                     7.0%     37.284µs
   boundary condition           0.1%      0.765µs
   rk4 stage algebra            9.2%     49.239µs

== 320x100 cells, 24 steps, 2.095241ms per step ==
   phase                       share    per step
   wind stress sampling        71.0%   1.487591ms
   shallow-water terms         13.4%    281.064µs
   coriolis                     6.6%    137.413µs
   boundary condition           0.1%      1.602µs
   rk4 stage algebra            9.0%    187.569µs
```

**Cost is linear in the cell count, and the two grids spend their time the same
way.** 66.6 ns per cell per step at 1.0° against 65.5 ns at 0.5° — under 2%
apart across a factor of four in size — and every share matches to within a
point. Neither grid falls off a cache cliff the other avoids, so a conclusion
drawn at one resolution holds at the other, and the 0.25° grid the benchmark
suite deliberately omits should cost four times the 0.5° one rather than
something worse.

The warm-up is not a detail. Measured from a cold process the coarse grid
appears to cost 2.2 times the fine grid's time per cell, which is the
processor's clock ramping rather than anything about the engine. A fixed *step*
count warms the fine grid and leaves the coarse one cold, so the example warms
for a fixed *duration*.

## Inside the two evaluators

The phase table says the right-hand side is 13% of a step; it does not say
which of its loops that is. `engine::profiling::TermProfiler` splits the two
evaluators into the fourteen array kernels they are built from, calling the
same kernels in the same order over the same buffers —
`engine/tests/step_profile.rs` asserts its tendency is bit-identical to
`ShallowWaterRhs::evaluate` followed by `CoriolisTerm::add_to_tendency`.

At 320 × 100, 96 evaluations, as a share of the kernels' own total:

```
   kernel                       phase    share    per evaluation
   d(h)/dx  centre -> u face       sw     7.4%      7.648µs
   d(h)/dy  centre -> v face       sw     5.8%      5.940µs
   -g'.dh/dx + taux/(rho.H)        sw     6.8%      7.001µs
   -g'.dh/dy + tauy/(rho.H)        sw     7.1%      7.329µs
   -r.u                            sw     6.9%      7.129µs
   -r.v                            sw     7.2%      7.462µs
   d(u)/dx  u face -> centre       sw     5.5%      5.706µs
   d(v)/dy  v face -> centre       sw     5.3%      5.496µs
   -H.(du/dx + dv/dy)              sw     7.9%      8.129µs
   -r.h                            sw     6.7%      6.930µs
   v -> u faces (4-point)         rot     9.6%      9.871µs
   u -> v faces (4-point)         rot     9.6%      9.948µs
   +f.v                           rot     7.1%      7.300µs
   -f.u                           rot     7.0%      7.270µs
```

**There is no hot kernel inside the evaluators.** Fourteen passes over
basin-sized arrays, each between 5% and 10% of the total, with the two
four-point C-grid interpolations of the Coriolis term the largest at ~9.6%
apiece because each reads four neighbours per output point rather than one or
two. The 1.0° table has the same shape.

That flatness is itself the finding, and it is the shape of a
**memory-bandwidth-bound** computation rather than an arithmetic-bound one:
every kernel does a handful of flops per element and its cost tracks how many
arrays it streams. Two corroborations. First, the kernels timed one at a time
total 103.2 µs per evaluation against 104.6 µs for the two evaluators left
whole to be inlined and fused — the two agree to about 1.5%, in either
direction from run to run, so keeping values in registers across kernel
boundaries buys nothing measurable, which is what one expects when the limit is
the memory system. Second, the RK4 stage algebra — 9.0% of a step in
`add_scaled` and `memmove`, doing one multiply-add per element and nothing
else — costs about what a right-hand-side kernel costs, which only makes sense
if both are paying for the traffic rather than for the arithmetic.

At 0.5° a step streams roughly 4.6 MB through the six basin-sized
`OceanState`s RK4 holds, which does not fit in this machine's L2.

## Why the wind sampler costs what it does

The mechanism, from `engine/src/forcing.rs` and the call graph in the attached
profile. `WindStressField::sample` walks every face of both components and
calls `wind.stress(x, y, t)` once per point:

- **Two indirect calls per point.** The scenario's wind is a `CompositeWind`,
  a `Vec<Box<dyn WindStress>>`, whose `stress` dispatches virtually to a
  `ScenarioWind`, which matches on its own enum to reach `SteadyTradeWinds`.
  Neither call can be inlined, so the per-point body cannot be vectorised and
  the compiler cannot hoist anything out of the loop.
- **A libm `exp` per point.** `SteadyTradeWinds` applies a Gaussian meridional
  decay, `exp(−(y/Ly)²)`, reached through the dynamic linker's stub. `exp` and
  its stub alone are 28.6% of the whole run.
- **Recomputed where nothing changed.** The stress is a function of `y` alone,
  so every one of the 320 points in a row recomputes the same `exp`. And
  `SteadyTradeWinds` is constant in time, so all four RK4 stages of every one
  of a run's 17 520 steps recompute the same field.

The counting, at 320 × 100: `τx` lives on 321 × 100 east/west faces and `τy` on
320 × 101 north/south faces, so a stage evaluates the stress 64 420 times and a
step 257 680 times, at 5.8 ns each.

The redundancy is not a defect of the forcing module's design. The
`WindStress` trait is a pure function of `(x, y, t)` on purpose
(`forcing.rs`): it is what lets a scenario stack a steady field, a season and a
burst without any of them knowing about the grid, and it is what the four-stage
re-sampling in `Solver::step_forced_by` needs in order to integrate a
time-varying wind correctly. Sampling it per point per stage is the *simple*
implementation of that contract, and until this profile there was no evidence
it was the expensive one.

## What this means for the rest of Epic 10

Stated as measurements and their consequences, not as decisions — the tickets
are the place for those.

- **T-10.3 (rayon) and T-10.4 (`f32`) are aimed at 13% of a step.** By
  Amdahl, parallelising the shallow-water evaluator perfectly across this
  machine's 10 cores takes a step from 2.095 ms to 1.84 ms — 1.14×. Extended to
  every phase *except* the wind, including the RK4 algebra, the ceiling is
  still only 1.35×. Whatever those tickets do, they should be measured against
  the whole step and not only against `rhs_evaluation`, or they will report a
  speed-up nobody feels.
- **The evaluators look bandwidth-bound, which cuts both ways.** Fourteen flat
  kernels with nothing to gain from fusion is the profile of a computation
  limited by memory traffic. That is a weak case for threading the arrays
  (more cores, same bus) and a strong one for `f32`: halving the width of every
  field halves the traffic, and the term table says traffic is what is being
  paid for. If T-10.3 runs first, its own before/after is the check on this.
- **The wind sampler is a ticket that does not exist yet.** Epic 10's remaining
  two are rayon and `f32`, and neither touches the forcing. Caching the sampled
  field for a wind that does not vary in time, or hoisting the `y`-only factor
  out of the zonal loop, are changes to `forcing.rs` under the existing
  `WindStress` contract, and either would be worth more than both remaining
  tickets together. Raising that is what this note is for; choosing to do it is
  a human's call, and it is not made here.

## What this note does not measure

- **The run around the loop.** These figures are the time loop. Scenario build,
  solver construction and the run directory are in `scenario_run` in
  `docs/benchmarks.md`, which writes exactly two frames so that the filesystem
  is a constant rather than a term that grows.
- **Any scenario but the control one.** A `SeasonalTradeWinds` or a
  `WindBurstAnomaly` sample a genuinely time-varying stress: the *per-point*
  cost is the same, but the redundancy identified above is smaller, and a
  cached field would be wrong for them.
- **Accuracy.** Whether a future optimisation changed the answer is Epic 07's
  validation suite, not this note.
- **Any machine but this one.** The shares held across two grids on one laptop.
  A machine with a different memory system could move them, which is why the
  instrument is committed and the commands to re-run it are below.

## Reproducing this

```sh
# the phase and term tables
cargo run --release --example profile

# a hot path for an external sampler: uninstrumented steps, nothing else
cargo run --profile profiling --example profile -- spin 30 &
sample $! 15 1 -f spin.sample                # macOS
perf record -F 999 -g -p $! -- sleep 15      # Linux
```

The `profiling` cargo profile is release code plus the debug info a sampler
needs to name a return address; it is separate from `release` so that what
ships and what is profiled are the same instructions built from the same flags.

The instrument is `engine/src/profiling.rs` — what it can and cannot see is in
its module documentation — driven by `engine/examples/profile.rs`, and
`engine/tests/step_profile.rs` is what holds it to profiling the engine's own
step.
