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
[basin]                          # optional; omitted, it is the Pacific below
western_longitude_deg = 120.0    # 120°E
eastern_longitude_deg = -80.0    # 80°W, counted eastward across the dateline
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5             # cell size, both axes

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
| `basin` | optional | Which stretch of ocean the run covers, in degrees, and how finely it is cut into cells. Omitted means the equatorial Pacific of `CONTEXT.md` — 120°E–80°W by 25°S–25°N at half a degree. |
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

Which stretch of ocean the scenario runs on, in degrees, and how finely it is
cut into cells. The section is optional, and so is every key in it: omitted,
each falls back to the equatorial Pacific of `CONTEXT.md` (*Basin*) —
120°E–80°W by 25°S–25°N in cells of half a degree. A scenario states a
boundary only when it means something other than the basin this project is
about.

Longitude is degrees east and it wraps: `-80.0` and `280.0` name the same
meridian, and the zonal span is always measured **eastward from the western
boundary**, so a basin may cross the dateline — the Pacific has to. Two equal
longitudes are therefore a basin of zero width, not one wrapped around the
planet. Latitude is degrees north, so a basin straddling the equator has a
negative southern boundary.

<!-- fields: BasinSection -->

| Field | Type | Required | Unit | Valid values |
| --- | --- | --- | --- | --- |
| `western_longitude_deg` | float | optional | degrees east | Any finite longitude. Omitted means `120.0` — the Maritime Continent edge of the Pacific. |
| `eastern_longitude_deg` | float | optional | degrees east | Any finite longitude, counted eastward from the western one. Omitted means `-80.0`, the South American coast; `280.0` is the same meridian written the other way. |
| `southern_latitude_deg` | float | optional | degrees north | Finite and on the planet, `−90 ≤ φ ≤ 90`. Omitted means `-25.0`. |
| `northern_latitude_deg` | float | optional | degrees north | Finite, on the planet, and strictly north of `southern_latitude_deg`. Omitted means `25.0`. Equal or inverted latitudes are refused with *"northern_latitude_deg is …, which is not north of southern_latitude_deg …"*. |
| `resolution_deg` | float | optional | degrees | Finite, strictly greater than 0, and dividing *both* spans into a whole number of cells. Omitted means `0.5`. |

<!-- end fields -->

`resolution_deg` is one number and not two because the cells are square: on
the equatorial beta-plane a degree of longitude and a degree of latitude are
the same degree of arc, so a basin stated in degrees has `dx = dy`. An
anisotropic grid is a numerical decision, and it would arrive with the ADR
that justifies it rather than as a second key.

### From degrees to metres

The solver works in metres, and the projection is one multiplication:

```text
metres per degree of arc = R·π/180 ≈ 111 195.08 m     with R = 6 371 008.8 m
```

`R` is the IUGG mean radius of WGS-84 — the same radius `β = 2Ω·cos φ / R` is
quoted from, so the geometry and the rotation describe one planet. The
`cos(φ)` convergence of the meridians is exactly the term the beta-plane
approximation drops, so it is deliberately not applied here: reintroducing it
would place the grid on a geometry the equations are not solved on.

That gives the derived quantities the rest of the file depends on:

```text
nx            = zonal span / resolution_deg            cells, east–west
ny            = meridional span / resolution_deg       cells, north–south
dx_m = dy_m   = resolution_deg · 111 195.08            metres
```

`x` is measured **east from the western boundary**, which is therefore
`x = 0`; `y` is measured **north from the equator**, so the southern boundary
of the default basin sits at `−25 · 111 195.08 ≈ −2.78 × 10⁶ m`. Every `_x_m`
and `_y_m` elsewhere in the file — a wind burst's centre, for instance — is in
that frame.

The default basin is 160° by 50° at 0.5°, which is **320 × 100 cells** of
55 597.54 m: about 17 791 km of Pacific by 5 560 km of latitude.

Because every key defaults, a scenario that wants that basin need not say so —
and one that wants a coarser version of it states the single key it changes:

<!-- scenario -->
```toml
# No [basin]: the equatorial Pacific, at half a degree.

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 17520
output_every_n_steps = 24

[[wind]]
type = "seasonal_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0
relative_amplitude = 0.2
peak_time_s = 18144000.0
```

The three files in [`engine/scenarios/`](../engine/scenarios/) write the five
keys out anyway, so that a run's file records which basin it was on rather
than inheriting one that a later default could change.

### What is refused

A span that is not a whole number of cells is **refused, never rounded** —
rounding it would silently run a basin nobody asked for:

*"the basin spans 160.3 degrees of longitude, which is not a whole number of
cells of resolution_deg 0.5"*.

Whole is judged to a relative tolerance of `1e-9` of the cell count, which is
seven orders of magnitude looser than the binary rounding of decimal degrees
and still far tighter than any mis-specification worth catching: `1e-9` of a
half-degree cell is 56 µm. So a basin written in round degrees is never
refused, at any resolution.

The other refusals name their value the same way: a non-finite boundary or
resolution (*"…it must be a finite number of degrees"*), a latitude off the
planet, a non-positive `resolution_deg`, an axis shorter than a single cell,
and a resolution so fine that the cell count does not fit in a machine index.

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

Both axes enter through `κ_max`. A basin stated in degrees has square cells,
so the two terms are equal and the bound collapses to `0.8·dx_m/c`; the
general form is kept because the check is written against a `Spacing`, which
does not have to be square.

For the example — `resolution_deg = 0.5`, so `dx = dy = 55 597.54 m`, and
`c = 3.0 m/s` — the bound is ≈ 14 826 s, so `dt_s = 3600.0` sits well inside
it.

#### The rotation bound

Checked when the solver is built rather than when the file is parsed, so a
scenario can satisfy everything above and still be refused here. The rotation
pair `u̇ = +f·v`, `v̇ = −f·u` has eigenvalues `±i·f`, so the step has to resolve
the fastest inertial oscillation in the basin — the one at whichever meridional
wall lies further from the equator (see
[ADR-0007](planning/adr/0007-rotation-timestep-bound.md)):

```text
|f|_max                = β · max(|southern_latitude_deg|, |northern_latitude_deg|) · 111 195.08
max_stable_dt_rotation = 0.8 · 2√2 / |f|_max
```

Neither the spacing nor the wave speed appears: this bound depends on how far
north and south the basin reaches, so a meridionally tall basin is limited by
rotation while a wide, shallow one is limited by CFL. For the example — walls
at 25°S and 25°N, so `|f|_max = 6.39 × 10⁻⁵ s⁻¹` — the rotation bound is
≈ 35 390 s, more than twice the CFL bound, which is why the CFL one is the
binding constraint in all three worked scenarios. A basin whose walls are both
on the equator has `f ≡ 0` and no rotation limit at all; widening one to 60°
of latitude would make rotation the binding bound instead.

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
| `center_x_m` | float | required | m | Any finite position. The zonal centre `x₀` of the burst, in metres east of the basin's western boundary, which is `x = 0`. |
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
western_longitude_deg = 120.0
eastern_longitude_deg = -80.0
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5

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

1. `[basin]` — boundaries finite and on the planet, then `resolution_deg`
   positive, then the latitudes ordered, then each span a whole number of
   cells.
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
