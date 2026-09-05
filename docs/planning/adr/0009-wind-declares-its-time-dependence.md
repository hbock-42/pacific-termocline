# ADR-0009: A wind declares whether it varies in time

## Status
Accepted. Records a public trait method whose default is a correctness
contract, added by T-10.5.

## Context
T-10.2's profile ([`docs/performance-notes.md`](../../performance-notes.md))
found re-sampling the wind stress at **71% of a timestep**. The engine
evaluates the prescribed wind at every C-grid face at every one of RK4's four
stages — 257 680 evaluations per step on the control basin, each a pair of
virtual calls ending in a libm `exp` — and for `SteadyTradeWinds` every one of
them returns what the previous stage's did. The redundancy is not a defect of
the forcing module's design: `WindStress` is a pure function of `(x, y, t)` on
purpose, and sampling it per point per stage is the simple implementation of
that contract.

Avoiding it means reusing a field. Reuse is only correct while the field is
still the field of the instant being asked about, and *that* is a property of
the individual wind: `SeasonalTradeWinds` breathes with the year and
`WindBurstAnomaly` is a Gaussian in `t`, so a cache that never invalidated
would be right for the control scenario and silently wrong for both — the
failure T-10.5 names as "a bug, not an optimisation".

Nothing in the engine can work the property out for itself. `WindStress` is a
trait; its implementations may live outside this crate; and there is no way to
ask a Rust function whether it reads one of its arguments.

## Decision

**The trait asks the wind, and trusts the answer.**
`WindStress::time_dependence(&self) -> TimeDependence` returns `Steady` or
`Varying`, **and its default implementation returns `Varying`**. A held field
is reused on exactly two grounds: the instant asked for is bit-for-bit the
instant held — true of any pure function of `(x, y, t)`, so it needs no
declaration at all — or the wind declared itself `Steady`.

Three properties of that default carry the decision.

- **Silence is safe.** An implementation written before this method existed,
  or by someone who never read this ADR, is re-sampled at every instant and
  behaves exactly as it did before. The failure mode of forgetting is a slower
  run, not a wrong one.
- **The declaration is checkable where it is made.** `SteadyTradeWinds::stress`
  ignores its `t_s` parameter, three lines above the override that says so.
- **Composition is derived, not declared.** `CompositeWind` is `Steady` only if
  every component is, and `ScenarioWind` delegates: neither has an opinion of
  its own, so stacking a burst on the trades cannot lose the burst.

What a wrong declaration costs is nevertheless real, which is why this is an
ADR: a `WindStress` that returns `Steady` while depending on `t` freezes the
forcing of a whole run, and the run finishes and writes a plausible-looking
output file.

## Considered options

- **Ask the wind for its stress at two instants and compare.** Self-checking,
  and wrong: agreeing at two instants is not being constant, and a Gaussian in
  `t` agrees with itself either side of its peak. It would also cost two
  samplings to save four.
- **Decide it at the scenario loader.** `ScenarioConfig::build` knows which
  variant it built, so it could label the scenario steady without touching the
  trait. But then the property belongs to the loader rather than to the wind:
  a `WindStress` built in a test, in a benchmark, or by a future crate gets no
  answer, and the one place that knows the physics — the implementation whose
  `stress` does or does not read `t_s` — is not the place that states it.
- **Cache on the instant alone, with no declaration.** Sound with no promise
  from anybody, and it is half of what was built: RK4's four stages ask about
  three instants, and a step's last stage asks about what the next step's first
  asks about, so a wind that varies as fast as it likes samples twice a step
  instead of four times. It is not enough on its own — the control scenario
  would still re-sample every step, and the profile says that is 71% of it.
- **Separate every wind into a spatial shape times a scalar in time.** Every
  wind this engine ships is separable that way, and it would make seasonal and
  burst runs as cheap as steady ones. It is also a much larger contract: a
  composite of two differently-modulated components needs a cached field per
  component, and the trait grows from "what is the stress here" to "what is
  your shape, and what is your scale". Deferred, and reachable from here: a
  third variant of `TimeDependence` is a smaller change than a second trait.

## Consequences
- `TimeDependence` is public API. A wind that means to be cached must say so,
  and the compiler will not remind it — `engine/tests/wind_stress_cache.rs`
  counts the evaluations instead, and asserts the trait's default is `Varying`
  for a wind that declares nothing.
- Every wind that declares `Steady` needs a test that its run matches an
  uncached one bit for bit. That file has one per shipped wind, against a
  stepper that re-samples at every stage.
- The `Steady` claim is a *statement about `stress`*, so any change to a
  `stress` implementation has to be read against its `time_dependence`. Adding
  a `t`-dependent term to `SteadyTradeWinds` without changing its declaration
  would be a silent, run-wide error.
- A forcing owns its wind (`WindForcing<W>`), so a cached field cannot be
  paired with a wind other than the one it was sampled from. That is why the
  solver's older `step_forced_by`, which takes the wind per call, keeps its
  cache only for the length of one step.
