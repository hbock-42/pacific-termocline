# Scenario configuration reference

A **scenario** is the engine's unit of input: a complete, runnable
specification of one simulation — grid, physical parameters, wind forcing and
run length (`CONTEXT.md`, *Scenario*). It is written as a TOML file and read by
[`engine/src/scenario.rs`](../engine/src/scenario.rs), which is the authority
this page documents; the tables below are checked against that file by
`engine/tests/scenario_config_reference.rs`, so a field that appears on one
side and not the other fails CI rather than going quietly stale.

Three worked examples live in [`engine/scenarios/`](../engine/scenarios/), one
per scenario of `docs/planning/01-scientific-model.md`.

## The shape of the file

<!-- scenario -->
```toml
[basin]                          # required, exactly once
nx = 200
ny = 60
dx_m = 50000.0
dy_m = 50000.0

[physics]                        # required, exactly once
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]                            # required, exactly once
dt_s = 3600.0
total_steps = 17520
output_every_n_steps = 24

[[wind]]                         # zero or more, summed in the order written
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0
```

<!-- fields: ScenarioConfig -->

| Section | Required | Meaning |
| --- | --- | --- |
| `basin` | required | The shape of the basin, the size of a cell, and where its southwest corner sits. |
| `physics` | required | The constants of the scenario's ocean. |
| `run` | required | How long the run is, and how often it is saved. |
| `wind` | optional | The `[[wind]]` entries, in the order they are summed. Omitted or empty is a calm ocean — the undriven limit of the model, not a mistake. |

<!-- end fields -->

TOML sections may appear in any order; the engine reads by name. Every section
is declared `#[serde(deny_unknown_fields)]`, so a key this reference does not
list is an error naming the key, not a silently ignored typo — a misspelled
parameter would otherwise run a scenario nobody asked for.

Units are part of every name. The engine is SI throughout and never rescales a
value on the way in: metres, seconds, pascals, kilograms per cubic metre. A
field's suffix is its unit.

## `[basin]`

The grid the solver runs on and where it sits relative to the equator. `y` is
measured north from the equator, so a basin straddling it has a negative
southern edge.

<!-- fields: BasinSection -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `nx` | integer | required | cells | At least 1. `nx = 0` is refused with *"nx is 0; a grid needs at least 1 cell on each axis"*. A negative or fractional value is a TOML type error. |
| `ny` | integer | required | cells | At least 1, refused the same way as `nx`. |
| `dx_m` | float | required | m | Finite and strictly greater than 0. Cell width east–west. |
| `dy_m` | float | required | m | Finite and strictly greater than 0. Cell height north–south. |
| `western_edge_x_m` | float | optional | m | Any finite position. Omitted means `0.0`. Zonal position of the basin's western wall; `x` increases eastward. |
| `southern_edge_y_m` | float | optional | m | Any finite position. Omitted means `−ny·dy_m/2`, which straddles the equator symmetrically ([`Basin::centered_on_equator`](../engine/src/basin.rs)) — the configuration the idealized scenarios run in, with the equatorial waveguide centred so no wave is trapped against a wall. |

<!-- end fields -->

The basin's zonal extent is `nx·dx_m` and its meridional extent is `ny·dy_m`.
The example above is 200 × 60 cells of 50 km: 10 000 km of Pacific by 3 000 km
of latitude.

Both edges must be finite; an infinity or a NaN is refused with *"…is …; it
must be a finite position"*.

## `[physics]`

The fixed constants of the scenario's ocean — the parameters of the 1.5-layer
reduced-gravity model. They do not vary with position or time; the wind stress,
which does, is `[[wind]]`.

`β` and `ρ₀` are properties of the planet rather than of the experiment, so
they carry defaults and a scenario states them only when it deliberately varies
them.

<!-- fields: PhysicsSection -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `reduced_gravity_m_per_s2` | float | required | m/s² | Finite and strictly greater than 0. Reduced gravity `g' = g·Δρ/ρ₀`, the buoyancy contrast across the thermocline. A zero would collapse the wave speed. |
| `mean_thermocline_depth_m` | float | required | m | Finite and strictly greater than 0. The resting upper-layer thickness `H` — a *total* depth, unlike the anomaly `h` the model solves for; the total depth at a point is `H + h`. |
| `rayleigh_damping_per_s` | float | required | s⁻¹ | Finite and at least 0. Linear damping coefficient `r`; its inverse `1/r` is the damping timescale. `0` is allowed and is the undamped limit the wave tests run in. Negative is refused, because it would amplify rather than damp. |
| `beta_per_m_per_s` | float | optional | m⁻¹s⁻¹ | Finite and strictly greater than 0. Omitted means `2.3e-11`, the equatorial value of `β = 2Ω·cos φ / R` at `φ = 0`. A zero would remove the equatorial waveguide entirely. |
| `reference_density_kg_per_m3` | float | optional | kg/m³ | Finite and strictly greater than 0. Omitted means `1025.0`, the standard Boussinesq reference density for the upper tropical ocean (Gill, *Atmosphere–Ocean Dynamics*, appendix 3). It enters only through the wind-stress term `τ/(ρ₀·H)`, so a zero would divide by zero. |

<!-- end fields -->

Two derived quantities follow from this section and matter to the rest of the
file:

- **Kelvin wave speed** `c = √(g'·H)`, in m/s — the fastest signal in the
  model, and therefore what bounds `[run].dt_s`. The example's
  `g' = 0.06 m/s²` and `H = 150 m` give `c = 3.0 m/s`, the observed
  first-baroclinic Kelvin speed of the equatorial Pacific.
- **Equatorial deformation radius** `Le = √(c/β)`, in m — the width of the
  equatorial waveguide, ≈ 3.61 × 10⁵ m for that `c` and the default `β`. It is
  the physically motivated choice for the meridional scales of `[[wind]]`.

An out-of-range parameter is refused by name: *"reduced_gravity_m_per_s2 is 0;
it must be finite and greater than 0"*.

## `[run]`

How long the run is and how often it writes a frame. The run length is given in
steps rather than in seconds, because that is what makes the frame count exact.

<!-- fields: RunSection -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `dt_s` | float | required | s | Finite, strictly greater than 0, and no longer than *both* the CFL-stable maximum and the rotation limit below. Length of one solver step. |
| `total_steps` | integer | required | steps | Any non-negative integer (`u64`). Steps the run takes from its initial state. `0` is a run of the initial state alone. |
| `output_every_n_steps` | integer | required | steps | At least 1. Steps between saved frames; `0` is refused with *"every_n_steps is 0; a run writes a frame every N steps, and N must be at least 1"*. |

<!-- end fields -->

### Frames

Output is *decimated*: the run steps at `dt_s` and saves every
`output_every_n_steps`-th step, starting from the initial state at step 0. So

```text
frames        = total_steps / output_every_n_steps + 1   (integer division)
frame spacing = output_every_n_steps · dt_s              seconds of model time
```

A run whose length is not a whole number of intervals simply stops at the last
interval that fits; the schedule never rounds the run up and never moves a
sample. The example writes `17520 / 24 + 1 = 731` frames — one a day over two
years, plus the initial state.

### The two bounds on `dt_s`

A timestep has to clear two separate stability bounds, and they come from two
different oscillations. Both are the same RK4 stability region — `2√2`, where
the classic four-stage method's region meets the imaginary axis — read against
a different eigenvalue, and both hold back the same `0.8` safety margin
(`CFL_SAFETY_FACTOR`, chosen in T-01.3 as project policy, not a measured
constant).

#### The gravity-wave CFL bound

Checked by the scenario loader, *last*, once `[basin]` and `[physics]` are
known:

```text
κ_max         = 2·√(1/dx_m² + 1/dy_m²)
max_stable_dt = 0.8 · 2√2 / (c · κ_max)          with c = √(g'·H)
```

Both axes enter through `κ_max`, so on an anisotropic grid the bound is
stricter than the smaller spacing alone would suggest: the fastest mode is the
diagonal one.

For the example — `dx = dy = 50 km`, `c = 3.0 m/s` — the bound is ≈ 13 333 s,
so `dt_s = 3600.0` sits well inside it.

#### The rotation bound

Checked when the solver is built rather than when the file is parsed, so a
scenario can satisfy everything above and still be refused here. The rotation
pair `u̇ = +f·v`, `v̇ = −f·u` has eigenvalues `±i·f`, so the step has to resolve
the fastest inertial oscillation in the basin — the one at whichever meridional
wall lies further from the equator (see
[ADR-0007](planning/adr/0007-rotation-timestep-bound.md)):

```text
|f|_max              = β · max(|southern_edge_y_m|, |southern_edge_y_m + ny·dy_m|)
max_stable_dt_rotation = 0.8 · 2√2 / |f|_max
```

Neither the spacing nor the wave speed appears: this bound depends on how far
north and south the basin reaches, so a meridionally tall basin is limited by
rotation while a wide, shallow one is limited by CFL. For the example — an
equator-centred basin reaching `±1.5 × 10⁶ m`, so `|f|_max = 3.45 × 10⁻⁵ s⁻¹` —
the rotation bound is ≈ 65 585 s, five times looser than the CFL bound, which
is why the CFL one is the binding constraint in all three worked scenarios. A
basin with both walls exactly on the equator has `f ≡ 0` and no rotation limit
at all.

Neither bound is ever applied silently. A timestep past either is **refused,
never shortened**: the error names the value asked for and the largest one this
basin allows, and fixing it is the scenario's job.

## `[[wind]]`

Zero or more entries, each an array-of-tables block. The forcings are summed
pointwise, in the order they appear in the file — the equations are linear in
the stress, so "a burst superimposed on the trades" is addition. The order is
part of the scenario because it fixes the floating-point result and therefore
keeps runs byte-reproducible.

Every entry carries a `type` key naming the forcing. An unknown value is
refused with a list of the ones that exist. There is no nesting: a composite is
the *list* of entries, not an entry.

All three forcings produce a purely zonal stress; `τy` is identically zero. The
alizés are zonal to the accuracy this model cares about, and a meridional
stress would drive an Ekman response the linear core has nothing to say about
yet.

Times (`peak_time_s`) are measured in seconds **from the start of the run**,
not from any calendar epoch.

### `type = "steady_trade_winds"`

The control scenario: a steady easterly stress that does not vary with `x` or
with `t`, optionally decaying away from the equator as a Gaussian.

```text
τx(x, y, t) = τ₀ · exp(−(y / Ly)²)          τy = 0
```

<!-- fields: WindSection::SteadyTradeWinds -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `equatorial_zonal_stress_pa` | float | required | Pa | Finite and **strictly negative**. The stress `τ₀` on the equator. The alizés blow from the east, which is `τx < 0`; a positive or zero value is refused, because it describes some other wind. |
| `meridional_decay_scale_m` | float | optional | m | Finite and strictly greater than 0. The `y` at which the stress has fallen to `1/e` of its equatorial value. Omitted means **no meridional structure at all** — the `Ly → ∞` limit, not a large `Ly`; that uniform profile is the one case with a closed-form steady state, which is why the analytic tilt check runs in it. |

<!-- end fields -->

The physically motivated scale is the equatorial deformation radius
`Le = √(c/β)`, the width of the waveguide the stress is meant to drive; it is
scenario input rather than a constant because it is exactly the knob the
forcing-sensitivity work varies.

### `type = "seasonal_trade_winds"`

The same field breathing with the year: a `steady_trade_winds` profile scaled
by an annual harmonic, the same factor everywhere in the basin at a given
instant.

```text
τ(x, y, t) = τ_steady(x, y) · (1 + a·cos(2π·(t − t_peak) / T_year))
```

<!-- fields: WindSection::SeasonalTradeWinds -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `equatorial_zonal_stress_pa` | float | required | Pa | Finite and strictly negative, exactly as for `steady_trade_winds`. This is `τ₀` *before* modulation. |
| `meridional_decay_scale_m` | float | optional | m | Finite and strictly greater than 0; omitted means no meridional structure. Same meaning as for `steady_trade_winds`. |
| `relative_amplitude` | float | required | — | Dimensionless, in the closed interval `[0, 1]`. The amplitude `a` of the annual harmonic. Outside that range the harmonic turns negative somewhere in the year and the stress flips westerly, which is a wind burst wearing a season's name; at `a = 1` the basin goes momentarily calm once a year, the strongest season this forcing can describe. |
| `peak_time_s` | float | required | s | Any finite instant, in seconds into the run. The moment the alizés are strongest — the phase, written as a time rather than an angle so it carries a unit. |

<!-- end fields -->

The period `T_year` is **not** configurable: it is the mean tropical year,
365.2422 mean solar days = 31 556 926.08 s (*Astronomical Almanac*). A scenario
wanting some other period is asking for a different forcing, not for a
differently tuned season.

The modulation is a pure scaling, so it does not move the wind's structure in
`y`: the alizés of March and of September have the same shape and different
strength. Meridional migration of the wind belt is real and deliberately not
modelled.

### `type = "wind_burst_anomaly"`

An idealized westerly wind burst: a positive-`τx` anomaly, Gaussian in `x`, in
`y` about the equator, and in `t`. It is an anomaly *against* the alizés, so it
is normally listed after a trade-wind entry rather than instead of one.

```text
τx(x, y, t) = τ_burst · exp(−((x − x₀)/Lx)²) · exp(−(y/Ly)²) · exp(−((t − t₀)/Lt)²)
```

<!-- fields: WindSection::WindBurstAnomaly -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `peak_zonal_stress_pa` | float | required | Pa | Finite and **strictly positive**. The peak stress `τ_burst`. A westerly burst blows against the alizés, so a negative or zero value is refused: it would describe a strengthening of the trades, or no burst at all. |
| `center_x_m` | float | required | m | Any finite position. The zonal centre `x₀` of the burst, in the same coordinate as `western_edge_x_m`. |
| `zonal_scale_m` | float | required | m | Finite and strictly greater than 0. The `e`-folding scale `Lx`: the distance east or west of `x₀` at which the stress has fallen to `1/e` of its peak. |
| `meridional_scale_m` | float | required | m | Finite and strictly greater than 0. The `e`-folding scale `Ly` about the equator. The physically motivated choice is the deformation radius `Le = √(c/β)`, the waveguide the burst is meant to excite. |
| `peak_time_s` | float | required | s | Any finite instant, in seconds into the run. The moment `t₀` of the burst's peak. |
| `duration_s` | float | required | s | Finite and strictly greater than 0. The temporal `e`-folding scale `Lt`: the stress is `1/e` of its peak this long before and after `t₀`. Observed bursts last on the order of 10 days. |

<!-- end fields -->

Note that the burst is *unbounded in time* in the same sense as a Gaussian: it
is never exactly zero, only exponentially small away from `t₀`. Setting
`peak_time_s` a few `duration_s` into the run lets the trade-driven tilt spin
up first.

## Composing forcings

A `[[wind]]` list of two entries is the sum of the two fields at every point
and instant. The westerly-burst example is exactly that:

<!-- scenario -->
```toml
[basin]
nx = 200
ny = 60
dx_m = 50000.0
dy_m = 50000.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 17520
output_every_n_steps = 24

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0

[[wind]]
type = "wind_burst_anomaly"
peak_zonal_stress_pa = 0.04
center_x_m = 2000000.0
zonal_scale_m = 1000000.0
meridional_scale_m = 361000.0
peak_time_s = 31556926.08
duration_s = 864000.0
```

`−0.05 Pa` of trades plus at most `+0.04 Pa` of burst still leaves the equator
easterly, which is the observed regime: the burst weakens the alizés, it does
not reverse them.

A scenario may also carry a seasonal cycle and a burst together, or several
bursts, or none at all. An empty list is calm.

## Errors

Nothing in this format is silently corrected. An invalid scenario is a returned
error naming the section it came from, the value that was wrong and the bound
it violated — never a panic, and never a substituted "safe" value.

The sections are validated in a fixed order, and the first failure is the one
reported:

1. `[basin]` — cell counts, then cell spacing, then edge positions.
2. `[physics]` — each parameter in the order of the table above.
3. `[run]` — timestep positive, then output cadence non-zero.
4. `[[wind]]` — each entry in file order.
5. The gravity-wave CFL bound on `dt_s`, last, because it needs the wave speed
   `[physics]` implies and the spacing `[basin]` implies.

A file that is not valid TOML, is missing a section, names a forcing that does
not exist, or carries an unknown key fails before any of that, at parse time.

The rotation bound on `dt_s` is checked later still — when the run builds its
solver, not when it reads its scenario — so it is the one refusal that arrives
after a file has otherwise been accepted.

## See also

- [`engine/scenarios/`](../engine/scenarios/) — the three worked examples.
- [`engine/src/scenario.rs`](../engine/src/scenario.rs) — the definition this
  page documents.
- [`docs/the-physics-explained.md`](the-physics-explained.md) — what these
  parameters mean physically, in plain language.
- [`docs/planning/01-scientific-model.md`](planning/01-scientific-model.md) —
  the equations and the three validation scenarios.
