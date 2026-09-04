# Epic 03 — Wind Forcing (the Alizés)

## Goal
Implement the trade-wind forcing scenarios from
[01-scientific-model.md](../01-scientific-model.md) as pluggable functions
of `(x, y, t)`, feeding `τx, τy` into the Epic 02 dynamics.

## Scope
The three forcing scenarios (steady trade winds, seasonal cycle, wind-burst
anomaly), plus the plumbing that makes forcing a swappable strategy rather
than hardcoded into the solver.

## Out of scope
Deriving wind stress from a coupled atmosphere model (that's the Bjerknes
feedback in Epic 12).

---

### T-03.1: `WindStress` trait and steady trade-wind scenario
- **Description:** Define a `WindStress` trait (`fn stress(&self, x, y, t) ->
  (f64, f64)`), and implement the steady easterly trade-wind field: uniform
  or simple analytic `y`-decay profile, `τx < 0` on/near the equator, per the
  scientific model doc.
- **Deliverable:** `WindStress` trait + `SteadyTradeWinds` implementation in
  a new `forcing` module, plumbed into `shallow_water_rhs` from T-02.3 in
  place of the stub.
- **Acceptance criteria:** Steady forcing run to equilibrium produces a
  thermocline that is deeper in the west than the east (sanity-checked
  against the analytic tilt formula, full rigor in Epic 07).
- **Depends on:** T-02.5.

### T-03.2: Seasonal cycle scenario
- **Description:** `SeasonalTradeWinds`, modulating the steady field by an
  annual harmonic (`1 + a·cos(2π·t/T_year)`), amplitude and phase
  configurable.
- **Deliverable:** New `WindStress` implementation.
- **Acceptance criteria:** Output field's time series at a fixed point shows
  the expected annual periodicity (unit test via FFT or peak-detection on a
  short run).
- **Depends on:** T-03.1.

### T-03.3: Idealized westerly wind-burst anomaly
- **Description:** `WindBurstAnomaly`, a Gaussian-in-x, Gaussian-in-y,
  Gaussian-in-t westerly stress superimposed on a base scenario (composable
  — should be addable on top of steady or seasonal winds, not exclusive with
  them).
- **Deliverable:** `WindBurstAnomaly` implementation + a `CompositeWind`
  combinator so scenarios stack.
- **Acceptance criteria:** Injecting a burst on top of steady trade winds and
  running forward shows a visible (in raw field data — visualization is
  Epic 08/09) eastward-propagating thermocline-depth signal consistent with
  Kelvin wave behavior — a qualitative smoke test here, rigorous check in
  Epic 07.
- **Depends on:** T-03.1.

### T-03.4: Scenario configuration format
- **Description:** Decide how a scenario (grid size, physical params, which
  wind forcing + its parameters, run length) is described in a config file
  (TOML, matching Rust ecosystem convention) and deserialized into the
  concrete `WindStress` + `PhysicalParams` + grid setup.
- **Deliverable:** `ScenarioConfig` struct with `serde`-based TOML
  (de)serialization, and 3 example config files (one per scenario type).
- **Acceptance criteria:** Each example config loads and produces the
  corresponding `WindStress` implementation with the right parameters;
  invalid configs (bad grid size, unknown forcing type) fail with a clear
  error, not a panic.
- **Depends on:** T-03.1, T-03.2, T-03.3.
