# Performance notes

Where a timestep's time actually goes, measured before anything is optimised
(T-10.2). [`docs/benchmarks.md`](benchmarks.md) says *how fast* the engine is;
this note says *what it is spending that on*, so that the two tickets behind
it — a rayon-parallel inner loop (T-10.3) and an `f32` field layout
(T-10.4) — can cite a profile rather than an intuition, which
[CODING_STANDARDS.md](../CODING_STANDARDS.md) § *Performance* requires of them.

It now carries four measurements rather than one. Everything before
[*After T-10.5: the sampled field is cached*](#after-t-105-the-sampled-field-is-cached)
is T-10.2's baseline, taken with the solver, the right-hand side, the Coriolis
term and the integrator exactly as T-07.4 left them and nothing optimised at
all. [*After T-10.5: the sampled field is cached*](#after-t-105-the-sampled-field-is-cached)
is the same two instruments run again after the ticket that baseline provoked,
and it is where the current decomposition of a timestep is. Both are kept,
because the point of a before is to be compared with an after.
[*After T-10.3: parallelising the sweeps does not pay*](#after-t-103-parallelising-the-sweeps-does-not-pay)
is the third and
[*After T-10.4: `f32` field storage is measured and rejected*](#after-t-104-f32-field-storage-is-measured-and-rejected)
the fourth, and both are negative results: they are kept so that nobody spends
either ticket again.

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

  It was made: the paragraph above became T-10.5, and the last section of this
  note is what it did. The two bullets before it are left as they were written,
  because they are the reasoning that ticket was created from — but the
  arithmetic in the first of them is superseded, and the section below redoes
  it against the step that now exists.

## What this note does not measure

- **The run around the loop.** These figures are the time loop. Scenario build,
  solver construction and the run directory are in `scenario_run` in
  `docs/benchmarks.md`, which writes exactly two frames so that the filesystem
  is a constant rather than a term that grows.
- **Any scenario but the control one.** A `SeasonalTradeWinds` or a
  `WindBurstAnomaly` sample a genuinely time-varying stress: the *per-point*
  cost is the same, but the redundancy identified above is smaller. (T-10.5
  measured how much smaller, and what a field cache may and may not do about
  it: see below.)
- **Accuracy.** Whether a future optimisation changed the answer is Epic 07's
  validation suite, not this note.
- **Any machine but this one.** The shares held across two grids on one laptop.
  A machine with a different memory system could move them, which is why the
  instrument is committed and the commands to re-run it are below.

## After T-10.5: the sampled field is cached

The change is one sentence of design. A [`WindStress`](../engine/src/forcing.rs)
now declares its `TimeDependence` — `Steady` or `Varying`, defaulting to
`Varying` so that a wind which says nothing is re-sampled exactly as before —
and a `WindForcing` owns a wind together with the one field it is sampled into,
re-sampling only when the field in hand is not already the field of the instant
asked for. A run holds one across its whole time loop
(`Solver::step_with_forcing`), so the control scenario samples the wind **once
for the run** instead of 4 times a step.

### The benchmark

`cargo bench -p engine --bench scenario_run`, the same machine and the same two
workloads as the tables above. The two runs are back to back on this branch and
on `main`, with criterion's `--save-baseline`, so its comparison is between two
measurements taken minutes apart rather than against a figure from another day:

| workload | before (`main`) | after | criterion's own comparison |
|---|---|---|---|
| 320 × 100 | 513.11 ms / run, **467.7 steps/s** | 149.60 ms / run, **1 604.3 steps/s** | −70.90% [−71.08%, −70.71%], p < 0.05 |
| 160 × 50 | 137.40 ms / run, **1 746.7 steps/s** | 38.60 ms / run, **6 218.0 steps/s** | −70.93% [−72.12%, −70.20%], p < 0.05 |

**3.4× at both grids**, end to end through the entry point the `run` command
uses. That is the deliverable, and it is what a user feels: the 0.5° control
run that took half a second takes 150 ms.

It is also what the baseline predicted. A phase measured at 71.0% of a step,
removed, leaves 29.0% — a 3.45× ceiling, and −70.9% is that ceiling to within
the measurement. The agreement is the check on the baseline rather than a
coincidence: had the wind's share been mis-measured, this is where it would have
shown.

The 1.0° baseline run was the noisy one of the four (129.92 ms to 147.49 ms
across its ten samples, against a spread of under 3 ms in the other three), so
its interval is the wide one. The 0.5° figures are the ones to quote; the
earlier tables in this note are at that grid for the same reason.

### The new phase table

`cargo run --release --example profile`, 24 steps per grid after a quarter of a
second of warm-up, exactly as before:

```
== 160x50 cells, 24 steps, 162.979µs per step ==
   phase                       share    per step
   wind stress sampling         2.6%      4.253µs
   shallow-water terms         44.0%     71.659µs
   coriolis                    22.8%     37.126µs
   boundary condition           0.4%        652ns
   rk4 stage algebra           30.2%     49.286µs

== 320x100 cells, 24 steps, 634.357µs per step ==
   phase                       share    per step
   wind stress sampling         2.6%     16.635µs
   shallow-water terms         44.4%    281.946µs
   coriolis                    22.2%    140.850µs
   boundary condition           0.3%      2.111µs
   rk4 stage algebra           30.4%    192.814µs
```

**Wind stress sampling: 71.0% → 2.6%. The shallow-water right-hand side is now
the hot path it was always assumed to be, at 44.4%.** The two rows are not
quite like for like, and the second reading below is why: 71.0% was a cost every
step paid, while 2.6% is one sampling spread across the 24 steps the profiler
timed. A longer window makes it smaller, which is the opposite of what a
per-step share does.

Two readings that say the cache did what it claims rather than moving work
somewhere the table cannot see.

**Every other phase costs what it cost.** At 0.5°: shallow-water terms
281.1 µs → 281.9 µs, coriolis 137.4 µs → 140.9 µs, RK4 stage algebra
187.6 µs → 192.8 µs — within 3% of the baseline, in a decomposition whose five
shares sum to one by construction. The step got smaller because one phase left,
not because the clock moved.

**The 2.6% that is left is one sampling, amortised.** The profiler is rebuilt
after warm-up, so its 24 timed steps contain exactly one full sampling of the
field. 16.635 µs × 24 = 399 µs, against 372 µs for one sampling measured at the
baseline (1.4876 ms per step ÷ 4 stages): the same number to within 7%, on a
phase whose whole total is a few hundred microseconds. In a real run — the
control scenario is 17 520 steps — that one sampling is 0.002% of the run, and
what is left of the phase is the branch that decides not to sample.

### The sampled profile says the same thing

The baseline was only stated as a conclusion because two instruments with
different failure modes agreed on it, so the after is held to the same
standard. `/usr/bin/sample` against
`cargo run --profile profiling --example profile -- spin 30`, 15 s, 11 003
in-loop samples, kept at
[`docs/profiles/2026-09-05-m1-pro-320x100-cached.sample`](profiles/2026-09-05-m1-pro-320x100-cached.sample)
beside the baseline's artefact. (`spin` now holds a `WindForcing` across its
loop, because a sampler pointed at a loop that samples the wind differently
from a run would be profiling a program nobody runs.)

Its summary by self time, demangled and truncated here:

```
Sort by top of stack, same collapsed (when >= 5):
    ShallowWaterRhs::evaluate                                     3241   29.5%
    <OceanState as StateVector>::add_scaled                       2802   25.5%
    CoriolisTerm::add_to_tendency                                 1135   10.3%
    CGridOperators::face_y_to_face_x                               738    6.7%
    CGridOperators::face_x_to_face_y                               646    5.9%
    _platform_memmove  (in libsystem_platform.dylib)               603    5.5%
    CGridOperators::ddx_center_to_face                             562    5.1%
    CGridOperators::ddy_center_to_face                             455    4.1%
    CGridOperators::ddy_face_to_center                             421    3.8%
    CGridOperators::ddx_face_to_center                             388    3.5%
    boundary::hold_walls_at_rest                                    10    0.1%
```

**The wind is not in it at all.** `exp` was 24.1% of the baseline's samples,
the two `stress` implementations 36.1% between them, the lazy-binding stub
`DYLD-STUB$$exp` 4.5%, and the sampling loop inlined into the step closure
7.0%. None of the five now reaches the summary's five-sample threshold, in a
run of 48 000 steps that samples the field once.

Folded into phases, and set beside the timed table — the timed shares
renormalised over its four non-wind phases, because the spin's single sampling
is spread over 48 000 steps where the profiler's is spread over 24, so the two
instruments cannot be compared on that row:

| Phase | timed | sampled |
|---|---|---|
| shallow-water terms | 45.6% | 46.1% |
| coriolis | 22.8% | 22.9% |
| rk4 stage algebra | 31.2% | 31.0% |
| boundary condition | 0.3% | 0.1% |
| wind stress sampling | (one sampling, amortised) | below the threshold |

Within half a point on every row, from two instruments with different failure
modes — which is the agreement the baseline was stated on.

The spin itself reports 1 599 steps/s, against 464 alone and 445 under the
sampler at the baseline. As before, that is not the engine's speed — it is a
laptop under a 30-second sustained load with a sampler walking its stacks for
half of it — and the benchmark above is the figure to quote.

### What a wind that genuinely varies gets

The ticket is explicit that a cache correct only for steady winds is a bug, so
the reuse rule is stated in terms of the `WindStress` contract rather than of a
scenario. A held field may be returned on two grounds: the instant asked for is
the instant held, bit for bit — a `WindStress` is a pure function of
`(x, y, t)`, so the same `t` gives the same field — or the wind declared itself
`Steady`. There is no third ground.

That first ground pays a time-varying wind on its own. RK4's four stages ask
about three instants (`t`, `t + dt/2` twice, `t + dt`), and where a schedule's
`t_{n+1}` is the bit-identical `t_n + dt` — which `step as f64 * dt_s` gives for
the timesteps a scenario actually uses — a step's last stage has already
sampled what the next step's first stage asks for. So a `SeasonalTradeWinds` or
a `WindBurstAnomaly` run samples **twice per step instead of four times**, with
no promise about time dependence involved at all.

That is stated as a sampling count rather than as a speed-up because it is
asserted rather than measured: `engine/tests/wind_stress_cache.rs` counts the
evaluations at the `WindStress` itself, deriving the expected count from the
RK4 tableau. No benchmark workload has a time-varying wind — `BENCHMARK_WORKLOADS`
is the control scenario at two grids — so this note does not put a wall-clock
figure on it, and a halving of the sampling is not a halving of a step.

The correctness of all of it is `engine/tests/wind_stress_cache.rs` and the
suite that was already there: a cached run is compared bit for bit against an
uncached stepper that re-samples the wind at every stage, for steady, seasonal,
burst and composite forcings, and Epic 07's validations and
`engine/tests/wind_burst.rs` pass unchanged, with tolerances untouched.

### What this leaves for the rest of Epic 10

The Amdahl arithmetic the baseline did for T-10.3 and T-10.4 has to be redone,
because it was computed against a step that no longer exists. Against the
634 µs step at 0.5°:

- **T-10.3 (rayon) now attacks 44.4%, not 13.4%.** Parallelising the
  shallow-water evaluator perfectly across this machine's 10 cores takes a step
  from 634 µs to 381 µs — **1.67×**, against the 1.14× the same calculation gave
  before. Both evaluators together give 2.5×. The ticket's own instruction not
  to measure itself only against `rhs_evaluation` stands, and matters more now
  that the whole-step figure has something to show.

  That is a ceiling, and the ticket went and measured what was under it:
  [*After T-10.3*](#after-t-103-parallelising-the-sweeps-does-not-pay) is the
  answer, and it is that per-kernel threading costs several times more than it
  saves at this basin size. The bullet is left as it was written, for the same
  reason the two above it are.
- **T-10.4 (`f32`) now attacks 97% of a step.** The two evaluators and the RK4
  stage algebra are 97.0% of what is left, and the baseline's term table said
  all of it looks memory-bandwidth-bound — fourteen flat kernels, and an
  `add_scaled` that costs what a kernel costs. Halving the width of every field
  halves the traffic that argument says is being paid for.
- **The `exp` per point was not hoisted, deliberately.** The baseline's other
  suggestion — that the trade-wind stress is a function of `y` alone, so 320
  points of a row recompute one `exp` — is untouched. For a steady wind it is
  now worth 0.002% of a run. It is still the cost of the two samplings a
  time-varying wind takes per step, and that is where the case for it would have
  to be made, on a workload this suite does not have. Optimising it now would be
  the speculative micro-optimisation
  [CODING_STANDARDS.md](../CODING_STANDARDS.md) § *Performance* rules out.

## Reproducing this

```sh
# the phase and term tables
cargo run --release --example profile

# a hot path for an external sampler: uninstrumented steps, nothing else
cargo run --profile profiling --example profile -- spin 30 &
sample $! 15 1 -f spin.sample                # macOS
perf record -F 999 -g -p $! -- sleep 15      # Linux

# what a phase costs against the thread count, and where a row-split sweep
# starts to win (T-10.3) — both on a build that threads its sweeps
for n in 1 2 4 10; do RAYON_NUM_THREADS=$n cargo run --release --example profile; done
for n in 1 2 4 10; do RAYON_NUM_THREADS=$n cargo run --release -p termocline-grid \
    --example sweep_scaling; done

# what halving the width of a field buys the sweep over it (T-10.4)
cargo run --release -p termocline-grid --example width_scaling

# and what it costs: Epic 07's suite, tolerances untouched, over an engine
# whose prognostic states are stored at f32. A `--cfg` and not a cargo
# feature, because CI runs `--all-features` and this is not a build to gate on
RUSTFLAGS="--cfg f32_storage_probe" cargo test -p engine --no-fail-fast
```

The `profiling` cargo profile is release code plus the debug info a sampler
needs to name a return address; it is separate from `release` so that what
ships and what is profiled are the same instructions built from the same flags.

The instrument is `engine/src/profiling.rs` — what it can and cannot see is in
its module documentation — driven by `engine/examples/profile.rs`, and
`engine/tests/step_profile.rs` is what holds it to profiling the engine's own
step.

## After T-10.3: parallelising the sweeps does not pay

T-10.3 is the rayon ticket the section above re-armed: with the wind cached,
the two evaluators are 66.6% of a step rather than 20.0%, so perfect
parallelisation of them was worth **2.5×** where it had been worth 1.2×. That
is a ceiling worth chasing, so it was chased. It was measured, and it is not
there. **Threading the array sweeps makes the engine 3.6× slower at 0.5° and
9.1× slower at 1.0°**, and the ticket is closed on that measurement rather than
on the estimate that motivated it.

The implementation the numbers below were taken on is `812b337`, the first
commit of the pull request that closes T-10.3, reverted by its second. `main`
takes one squashed commit per ticket (CONTRIBUTING.md), so that commit never
reaches `main` — it stays in the pull request's own history, which is where to
check it out from and re-measure rather than take this section on trust. That
is the same reason the sampled profiles above are committed as artefacts.

### What was built

One primitive, `termocline_grid::sweep::write_rows`, and every kernel through
it: the six C-grid writers, the two four-point interpolations, the three
pointwise shallow-water kernels and the Coriolis accumulation. A `Field2D` is
row-major, so rows are contiguous and disjoint, and a sweep splits across them
with `par_chunks_exact_mut`.

**Correctness was never the difficulty, and it is worth saying why.** The
right-hand side contains no reduction — no sum whose order a work-stealing
scheduler could permute. Every kernel is a map: each output point is written
once, from the same inputs, in the same order, on whatever thread runs the row.
So the change is bit-identical by construction rather than by tolerance, and
`engine/tests/parallel_determinism.rs` — added by that commit — holds it to
that on 1, 2, 3 and 10 workers, three being the count that divides neither
benchmark basin's row count. Epic 07's validations, which compare against
analytic results rather than against another run, passed unchanged beside it.
The determinism CODING_STANDARDS.md § *Correctness and failure* requires was
never at risk. The cost was.

### The benchmark

`cargo bench -p engine`, the same machine and the same two workloads as
everything above, criterion's `--save-baseline` between two runs taken minutes
apart:

| workload | serial | rayon | criterion's own comparison |
|---|---|---|---|
| `rhs_evaluation` 320 × 100 | 68.87 µs | 302.18 µs | **+335.20%** [+330.60%, +339.51%], p < 0.05 |
| `rhs_evaluation` 160 × 50 | 17.58 µs | 226.43 µs | **+1179.5%** [+1163.0%, +1194.3%], p < 0.05 |
| `scenario_run` 320 × 100 | 150.72 ms | 544.83 ms | **+261.25%** [+256.36%, +265.48%], p < 0.05 |
| `scenario_run` 160 × 50 | 39.00 ms | 353.92 ms | **+805.62%** [+795.70%, +816.34%], p < 0.05 |

The last column is criterion's, computed from the two runs' whole sample
distributions; the two middle columns are their point estimates. Dividing the
middle columns gives a percentage a little different from the last one — +339%
against +335% on the first row — for that reason and not for any other.

Not a disappointing speed-up. A regression of a factor of several, at both
grids, in the isolated kernel benchmark and end to end through the entry point
the `run` command uses.

### Why: the sweeps are too short

`cargo run --release --example profile`, the two evaluator phases at 320 × 100,
against `RAYON_NUM_THREADS`:

| | shallow-water + coriolis, per step | against serial |
|---|---|---|
| serial (no rayon at all) | **421.6 µs** | 1.00× |
| rayon, 1 thread | 840.7 µs | 0.50× |
| rayon, 2 threads | 732.2 µs | 0.58× |
| rayon, 4 threads | 1 023.2 µs | 0.41× |
| rayon, 10 threads | 1 921.1 µs | 0.22× |

**The one-thread row is the fork/join cost with no parallelism in it at all.**
Going through rayon's split-and-join on a single thread doubles the evaluators
before a single row moves anywhere: 419 µs a step. A step performs 56 sweeps —
fourteen kernels at each of RK4's four stages — so that is about 7.5 µs of
scaffolding around a sweep the term table above says is 5 to 10 µs of work.
T-10.2's finding that the right-hand side is *fourteen flat kernels with no hot
one* turns out to be the same finding as *no sweep long enough to pay for being
handed to a thread pool*.

Granularity was tuned, in case the split was merely too fine. The table above
gives each worker at least 1 024 points; the coarsest possible split, one task
per thread per sweep, gives 683 µs at 2 threads, 1 002 µs at 4 and 2 055 µs at
10 — the same shape, no better. The 160 × 50 grid behaves the same way and
worse, as a smaller basin must: 109.7 µs serial against 350.7 / 411.0 / 687.5 /
1 193.3 µs at 1 / 2 / 4 / 10 threads.

### And where the memory system takes over

The table above is consistent with two different stories — sweeps too short to
amortise a fork, or a bandwidth-bound computation choking on extra cores — and
they point at different next tickets, so they are worth separating.
`termocline-grid/examples/sweep_scaling.rs` separates them by holding one sweep
fixed and growing the field. It streams three arrays for one flop, `c = a + b`,
which is the shape of every kernel in the term table, and reports nanoseconds
per point for a plain loop and for `write_rows`. Parallel against serial, so
above 1.00× threading wins:

| field | 1 thread | 2 threads | 4 threads | 10 threads |
|---|---|---|---|---|
| 320 × 100 — *the engine's 0.5° basin* | 0.52× | 0.52× | 0.42× | 0.15× |
| 640 × 200 | 0.78× | 0.94× | 1.03× | 0.48× |
| 1 280 × 400 | 1.03× | **1.80×** | **2.16×** | 1.24× |
| 2 560 × 800 | 0.95× | 1.50× | 1.55× | 1.44× |
| 5 120 × 1 600 | 0.80× | 1.07× | 1.61× | 1.26× |

Three things are in that table, and the first two are the answer.

**The crossover is at about 250 000 points, and the engine's fields are 32 000.**
A row-split sweep starts to win somewhere between 640 × 200 and 1 280 × 400 —
eight to sixteen times the 0.5° basin. Below it the fork and the join cost more
than the sweep, which is the 320 × 100 row: threading loses at *every* thread
count, by a factor of two at best and nearly seven at worst. So the first story
is the one that explains the engine's numbers.

**Where threading does win, it wins about twice, and never more.** The best
figure anywhere in the table is 2.16×, on a field sixteen times the engine's,
with four workers. Ten workers never beat four at any size, and past
2 560 × 800 the speed-up falls away again as the fields leave cache. That is
the second story, and it is real — it is the *ceiling* rather than the
explanation, and it is a low one, which is what the baseline note predicted in
as many words: *"a weak case for threading the arrays (more cores, same bus)"*.

**Neither story is about the number of cores.** This machine has ten and the
table never rewards more than four. That is worth carrying into T-10.4: the
constraint the engine is up against is how fast values move, not how many
places there are to compute them.

### What this closes, and what it does not

**T-10.3 is closed as measured-and-rejected.** The acceptance criterion is a
meaningful speed-up on a multi-core machine, and a change that is 3.6× slower
does not meet it; shipping it because the ticket exists would be the opposite
of CODING_STANDARDS.md § *Optimize against a measurement*. The estimate that
re-armed it — 2.5× if the evaluators parallelised perfectly — was an Amdahl
ceiling, and a ceiling is what a change cannot exceed rather than what it gets.
Nothing about it was wrong except the assumption underneath every Amdahl
calculation: that the phase can be parallelised for free.

**It does not close threading as such, but it bounds it.** What is measured
here is per-kernel parallelism, and a parallel region spanning the whole step —
a strip-decomposed solver holding a thread team across the time loop, with halo
rows and a barrier between kernels — would still synchronise 56 times a step,
but on a barrier between workers that are already awake rather than on a fork
and a join that wake them. The scaling table is what that design would be
chasing, and it says two things about it: at the engine's field size a sweep is
below the crossover *whatever* the synchronisation costs, and even far above
the crossover the prize is about 2×, not the 2.5× Amdahl offered or the 10× the
core count suggests. It would become interesting at a basin of a few hundred
thousand cells — 0.125°, four times finer than anything the benchmark suite
runs, which `benchmark.rs` calls "a benchmark nobody runs on a laptop". It is
not proposed here.

**Nor does it close the machine question.** These figures are one laptop's, and
a machine with more memory channels per core could move the crossover and the
ceiling both. That is the caveat every section of this note carries; it is why
the reverted commit and its two instruments are on the record.

**What is left is T-10.4.** It attacks 97% of a step rather than 66%, and it
attacks the quantity both instruments above say is actually being paid for —
traffic — by halving the width of every field, rather than adding claimants to
a bus that has just been measured saturating at four. The case for it is
stronger than it was, and it is stronger *because* of this measurement.

## After T-10.4: `f32` field storage is measured and rejected

T-10.4 is the ticket every section above has been arming. Its case was as
strong as a performance case gets: the section before this one measured the
engine's kernels as *bandwidth-bound* rather than arithmetic-bound, and a
bandwidth-bound computation is exactly the one that pays for the width of its
data. Halving that width should have been worth about `2×` on the 97% of a
step the two evaluators and the RK4 stage algebra now hold.

It is rejected, and not on the performance. **The prize is real and it is
about what was predicted; the accuracy cost is a validated Epic 07 tolerance,
and seven more derived budgets besides.**

### What was measured, and with what

Two instruments, because the ticket is two questions.

**The prize** — `termocline-grid/examples/width_scaling.rs`. It holds the
kernel fixed and varies the width: `c = a + b` over three arrays for one flop,
the shape of every kernel in the term table above, at the engine's own
320 × 100 and up the same size ladder T-10.3's thread scaling used. It is a
*ceiling*: no boundary handling, no interpolation, no scratch, nothing but
streaming.

**The cost** — `engine/src/precision.rs`, behind the `f32_storage_probe`
`--cfg`. It stores every prognostic state at `f32` and leaves the arithmetic at
`f64`, and it is an *exact* emulation of that pair rather than an approximation
of it: an `f32` widens to an `f64` losslessly, the arithmetic between two
stores is `f64` in both worlds, and rounding at every store leaves the field
holding the bits an `f32` field would hold. Without the `--cfg` the rounding is
not cheap, it is absent — the code is not compiled — so the engine that ships
is the engine `cargo test` validates, and `engine/tests/f32_field_storage.rs`
guards that.

A `--cfg` rather than a cargo feature for one reason: CI runs
`cargo test --workspace --all-features`, which would switch a feature on. An
instrument that the project's own gate enables is not an instrument.

**Be precise about which configuration this is.** The ticket describes "`f32`
for the bulk grid data while keeping `f64` accumulation where precision
matters (e.g. long-run energy conservation, per T-07.5)". What is measured here
is the first half taken to its conclusion: *every* stored state narrow,
including the one RK4 accumulates the step into, because a state field stored as
`f32` has nowhere else to keep the result of `state += w·dt·k`. Keeping an
accumulator wide is a different layout — it needs a second, wide copy of the
state, which is where the *"compensated accumulation"* bullet at the end of this
section goes — and this measurement does not settle it. It settles the layout
the ticket's own title names.

The point of a probe rather than a rewrite is that it re-runs **Epic 07's own
suite**, at the narrower width, with **not one tolerance touched**. A rewrite
would have taken the solver, the two evaluators, the operators, the forcing and
every test with it, and would have had to be believed before it could be
measured. The probe is measured first.

### The prize: about `2×`, as predicted

`cargo run --release --example width_scaling`, on the machine of *What was
measured, and on what* above:

| field | `f64` | `f32` | narrowing |
|---|---|---|---|
| 320 × 100 — *the engine's 0.5° basin* | 0.1959 ns/pt | 0.1003 ns/pt | **1.95×** |
| 640 × 200 | 0.1954 ns/pt | 0.0975 ns/pt | 2.00× |
| 1 280 × 400 | 0.2643 ns/pt | 0.0984 ns/pt | 2.69× |
| 2 560 × 800 | 0.3171 ns/pt | 0.1722 ns/pt | 1.84× |
| 5 120 × 1 600 | 0.3295 ns/pt | 0.1703 ns/pt | 1.94× |

**Halving the width halves the time, at every size including the engine's.**
That is the signature of a computation paying for traffic — a sweep doing one
flop per point does not get twice as fast because its operands got narrower
unless the operands were what it was waiting for. It is also the cleanest
confirmation yet of the baseline note's *bandwidth-bound* reading, taken by an
instrument that has nothing to do with the one that reading came from.

Read it as a ceiling, not as a speed-up. The engine's kernels interpolate,
handle boundaries and touch scratch, and the wind sampling and the run around
the loop do not narrow at all. What a narrowed engine would actually get is at
most this, on at most 97% of a step.

### The cost: eight derived budgets, one of them Epic 07's

`RUSTFLAGS="--cfg f32_storage_probe" cargo test -p engine --no-fail-fast`.
Thirteen tests fail. Five of them are the probe's own reach rather than the engine's
accuracy — `tests/step_profile.rs` and `tests/wind_stress_cache.rs` compare the
solver against hand-built reference steppers (`StepProfiler`, `UncachedStepper`)
that sit outside `solver::integrate` and so are not narrowed with it, and a
bit-identity assertion between a narrowed path and a wide one can only fail.
They are named here so the count below is the honest one.

The other **eight are the measurement**, and every one of them is a bound
somebody derived:

| suite | what it asserts | budget | at `f32` |
|---|---|---|---|
| **T-07.4** `steady_wind_tilt` | the tilt error falls at the predicted second-order ratio `0.2523` | ±4.465×10⁻³ | **0.2372** — 3.4× the budget |
| **T-07.4** `steady_wind_tilt` | a steady closed basin holds no net thermocline anomaly | 4.159×10⁻⁹ of the tilt | **1.677×10⁻⁸** — 4.0× |
| T-04.2 `no_normal_flow` | a closed, undamped basin leaks no volume through its coasts | 467.0 m³ | **1 285 225 m³** — 2 752× |
| T-02.4 `rayleigh_damping` | `exp(−r·t)` is reached at RK4's fourth order | order 4 ± 0.15 | **order 2.06** |
| T-02.4 `rayleigh_damping` | `Ė = −2·r·E` holds the whole way down | 10⁻⁶ relative | **1.005×10⁻⁶** at step 236 |
| T-03.2 `seasonal_wind` | a linear core generates no harmonics of the annual line | 10⁻⁹ of the line | **1.1×10⁻⁸** — 11× |
| T-03.3 `wind_burst` | a burst superposes exactly on the trades | 10⁻¹⁰ of the tilt | **43.136696 m against 43.136726 m** — the seventh figure, where the budget allows the eleventh |
| T-01.2 `time_stepping` | one step of a uniform `h` is RK4's amplification polynomial | 10⁻¹⁴ relative | **8.0×10⁻⁸** — 8×10⁶ |

**T-07.5 passes.** The spec names long-run energy conservation as the
precision-critical case, so it is worth reporting that all four tests in
`engine/tests/conservation.rs` are green at `f32`, including the eight-crossing
undamped drift against its derived bound and its second-order refinement check.
That is the one budget in the suite with an 11× margin, and 11× is enough to
absorb this. It is also the reason the ticket could not be closed on T-07.5
alone: the budgets that fail are the *tight* ones, and they are tight because
somebody derived where the run should land rather than allowing for where it
might.

The two T-07.4 rows are what close the ticket, and how they fail is worth
spelling out, because that file's *point* check passed.

**T-07.4's headline comparison — the settled tilt against the analytic damped
closed form — passes at `f32`.** It lands at 2.489×10⁻⁵ of the tilt against
its 2.5092×10⁻⁵ budget, marginally *closer* to the continuous solution than the
`f64` run's 2.5076×10⁻⁵. That is not the change being harmless. It is the
`f32` run having drifted 1.9×10⁻⁷ off the discrete closed form it used to sit
on: its departure *moved* by 1.9×10⁻⁷ of the tilt, about 120× the 1.6×10⁻⁹ of
margin the budget leaves, and it moved in the direction that happened to be
toward the continuous one. A bound is passed by an error of
any size below it, and `docs/validation-report.md` says so in as many words:
*"the point checks are bounds, and are generous by design"*.

**What catches it is the rate.** The same file's convergence test halves the
cell width and requires the departure to fall by the ratio the two closed forms
predict — the check the validation report keeps precisely so that a passed
bound is not mistaken for a validated run. At `f64` the measured and predicted
ratios agree to ten significant figures. At `f32` the fine run's error is no
longer the discretisation's, it is the round-off floor, so 2.489×10⁻⁵ becomes
5.905×10⁻⁶ where 6.28×10⁻⁶ was predicted: a ratio of 0.2372 against 0.2523,
outside the budget by 3.4×. **CODING_STANDARDS.md § *Convergence over point
checks* is what fails here, and it is what it is there for.**

### Why the arithmetic that predicted a bigger failure was wrong

Worth recording, because it is the estimate anyone re-arguing this ticket will
reach for first. T-07.4's tolerance carries a round-off term written as one
rounding per cell per step, none cancelling: `steps·N·ε`, which is 2.0×10⁻¹¹
of the tilt at `f64` and would be 1.1×10⁻² at `f32` — four hundred times the
whole tolerance. The measured departure is 2.5×10⁻⁵.

The bound is not wrong, it is worst-case by construction, and the reason the
run is nowhere near it is the physics: the channel is *damped*, so a rounding
made at step `n` has decayed by `e^{−r·(N−n)·dt}` by the end. Round-off does
not accumulate over 1 500 steps in a system that relaxes; it reaches a floor
set by the damping time. So the point check survives — and the convergence
check, which is a statement about that floor rather than about the total, does
not. Any future argument from `steps·N·ε` should expect to be off by three
orders of magnitude in a damped run, and should not conclude anything from
that in either direction.

### What the probe does not narrow, and which way it points

The probe rounds where a *state* is stored: the prognostic state, RK4's stage
state, and its four stage tendencies. A real `f32` field layout would also
narrow the divergence scratch of `ShallowWaterRhs`, the two interpolation
buffers of `CoriolisTerm`, the sampled `WindStressField`, and the gradient
written into the tendency before being turned into an acceleration in place.

Every one of those is a rounding not performed here, so the round-off a
narrowed engine injects is at least the round-off injected above.

That argument is about *magnitudes*, and it should not be stretched past them.
Six of the eight rows are magnitudes — a leak, a residual, a discrepancy — and
for those the bound reads straight: more rounding does not make a basin leak
less volume or a superposition close better. The two rows that are *rates* —
the convergence ratio and the measured order — are not monotone in added
round-off, and could move either way under a fuller narrowing.

Which is why it matters that the ticket does not close on a rate. Both T-07.4
rows are named above; the second, *a steady closed basin holds no net
thermocline anomaly*, is a magnitude, and `f32` storage misses it by 4.0×. The
convergence row is the more *interesting* failure — it is the one that explains
what went wrong — but the deciding one would still be there without it.

### What this closes, and what it does not

**T-10.4 is closed as measured-and-rejected**, on its own acceptance criterion:
*"if made, Epic 07 validation tests still pass within a re-justified
tolerance"*. They do not pass, and re-justifying the tolerance is not on the
table — AGENTS.md § *Never move the goalposts*, and a suite that is green
because a budget was widened to fit a faster engine is worse than a red one.
The deliverable is therefore its other branch, *"a documented decision not to
make it"*, and this section is it.

Nothing about the shipped engine changes, and `cargo bench -p engine` on this
branch says so rather than asserting it: `rhs_evaluation` 17.930 µs at 160 × 50
and 69.671 µs at 320 × 100, `scenario_run` 38.401 ms and 150.73 ms — the same
figures as T-10.3's `main` column above (68.87 µs, 150.72 ms, 39.00 ms) to
within their intervals. The probe costs nothing because with the feature off it
is not there.

**It does not close narrowing as such; it prices it.** What the measurement
rejects is `f32` *everywhere the state lives*, which is what the ticket
proposed. Three narrower things are untouched by it and would each need their
own measurement:

- **Output, not state.** `termocline-format` already writes what a visualizer
  reads; nothing in Epic 07 is asserted about a frame's width. That is a format
  ticket, not a solver one.
- **A mixed layout.** The failures cluster on quantities that are *exact
  identities* of the discrete scheme — volume, superposition, the amplification
  polynomial — rather than on quantities that are already approximations.
  A layout narrowing only the fields no identity is written in would have to
  name which those are, and Epic 07 currently says: fewer than one might hope.
- **Compensated accumulation** — the ticket's own "`f64` accumulation where
  precision matters", spelled out as a layout. The RK4 stage algebra is 30% of
  a step and one multiply-add per element, and every row above that fails on an
  *exact identity* of the scheme fails on something that identity accumulates.
  Whether narrow fields with a wide or Kahan-summed accumulator recover those
  identities is a real question, and this measurement does not answer it — it
  measures the layout with no wide accumulator anywhere. It is also a
  considerably larger change than the one just rejected, since a wide
  accumulator is a second copy of the state and gives back some of the traffic
  the narrowing bought, so it should be argued from a profile of its own.

**And it does not close the machine question**, which every section of this
note carries: these are one laptop's figures, and the 1.95× ceiling in
particular is a property of one memory system. The instrument is committed and
the command is below.
