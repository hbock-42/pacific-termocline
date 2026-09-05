# ADR-0010: The atmospheric wind response is diagnosed per RK4 stage, not prescribed as a `[[wind]]` entry

## Status
Accepted. Records where the Epic 12 wind feedback sits in the engine's
structure, decided by T-12.2.

## Context
T-12.2 closes the Bjerknes loop: `τx` has to respond to `T'`
(`CONTEXT.md`, *Bjerknes feedback*). Everything about the forcing built so far
assumes the opposite. `WindStress` is documented as a **pure function of
`(x, y, t)`**, and both of the things built on that contract depend on it:

- `CompositeWind` sums implementations that each depend on nothing but
  position and time, and a scenario names them in `[[wind]]` entries that are
  read from a file before the run starts;
- the T-10.5 cache ([ADR-0009](0009-wind-declares-its-time-dependence.md))
  reuses a sampled field whenever the *instant* repeats — and RK4 asks about
  `t + dt/2` twice per step, with two different states.

A wind that reads `T'` is not a function of `(x, y, t)`. Left inside either
mechanism it breaks something: as a `[[wind]]` entry it would be a forcing the
config file cannot fully describe, and inside a `WindForcing` it would silently
serve the second `t + dt/2` stage the wind of the first, which is a wrong
integration rather than a slow one.

The obvious alternatives:

1. **Make it a `[[wind]]` type** and give `WindStress::stress` access to the
   state. Every prescribed forcing would then take an argument it must not
   read, the trait would stop being the pure function its contract and its
   cache are built on, and `WindStress` implementations outside this crate
   would break.
2. **Make it another right-hand-side term**, like `SstTerm`, writing directly
   into `∂u/∂t`. It would work, but it would put the wind stress in two places:
   the momentum equation's `τx/(ρ₀·H)` would no longer be the whole surface
   stress, and a run's frames would record a stress the ocean did not feel.
3. **Keep a snapshot of `T'` inside a `[[wind]]` implementation** and refresh
   it from outside. This is the decision below, minus the honesty about where
   the refresh happens — the failure mode is a stale snapshot serving a stage
   silently, which is exactly what the cache would do.

## Decision

**The response is a `WindStress`, but it is not part of the prescribed wind.**
Three parts:

- `SstWindResponse` implements `WindStress`, so it composes with every existing
  forcing exactly as T-03.3's deliverable asks — including inside a
  `CompositeWind`. Its stress is `μ·⟨T'⟩·exp(−(y/L_a)²)`, a pure function of
  position given the one number `⟨T'⟩` it holds.
- That number is set by `SstWindResponse::observe`, called **once per RK4
  stage, with the state of that stage**, from `CoupledWind::at`. So each
  stage's right-hand side is a function of that stage's state and time, which
  is what an ODE right-hand side is; no stage sees another stage's wind, and
  the purity the trait requires holds everywhere between two `observe` calls.
- `CoupledWind` holds the prescribed forcing and the response **apart**, and
  sums their fields. The prescribed half keeps ADR-0009's cache untouched; the
  response half is re-sampled every stage and never enters it. The sum is the
  same superposition `CompositeWind` performs — the equations are linear in the
  stress — so what a step reads, and what a frame records, is one field
  carrying the whole surface stress.

**The feedback strength is a parameter of `[sst]`, not a `[[wind]]` entry.**
`wind_feedback_strength_pa_per_k` and `wind_response_meridional_scale_m` live
in the section that switches the SST coupling on, because the response reads
`T'` and only that section makes `T'` exist. A scenario cannot ask for one
without the other, and a `[[wind]]` list stays what it has always been: the
winds the file prescribes.

## Consequences

- **The validated core is untouched, twice over.** A scenario with no `[sst]`
  section takes the `WindForcing` path unchanged. A scenario with one at
  `μ = 0` adds a field of exact zeros, and T-12.2's regression asserts the
  result bit for bit against the prescribed run — the same argument T-12.1
  made for its own extension.
- **`StageForcing` grew a state argument** and became public. The two
  prescribed forcings ignore it, which is what makes them prescribed; the
  signature is where a reader learns that a forcing may be diagnosed.
- **A coupled run costs one more wind-stress field**, counted into the memory
  budget (`docs/scenario-config-reference.md`).
- **The response is not a Gill model.** Nothing here integrates the atmospheric
  equations: the pattern is fixed and the amplitude is a regression on one SST
  index, which is the *statistical* atmosphere
  `docs/planning/01-scientific-model.md` § *Phase 2* asks for. If a later
  ticket wants the zonal structure of a real Gill solution — the westerly
  anomaly west of the heating — it is a new `WindStress` behind the same
  `observe`-then-sample contract, not a change to this one.
- **`stress` on a `SstWindResponse` is only meaningful after `observe`.** That
  is a real sharp edge, and it is why the type is constructed by
  `CoupledWind::new` in every path a run takes: the object that owns the
  response is the object that refreshes it.
