# Pacific Thermocline

A physically-grounded simulation of the equatorial Pacific Ocean's thermocline
and the wave dynamics behind ENSO. This glossary is the canonical vocabulary
for the project: issue titles, function and type names, test names, and
documentation all use these terms and these symbols.

Numerical-method vocabulary (Arakawa C-grid, RK4, CFL condition, staggering)
is deliberately **not** defined here — it is implementation, and it lives in
[ADR-0003](docs/planning/adr/0003-numerical-scheme.md).

## Ocean state

**Thermocline**:
The sharp temperature boundary separating the warm, wind-mixed surface layer
from the cold deep ocean. Its depth and tilt across the basin are the central
object of this simulation.

**Thermocline depth anomaly** (`h`):
Departure of the thermocline from its mean depth, in metres. **Not** the total
depth: total depth is `H + h`. A positive `h` means a deeper-than-average
thermocline.
_Avoid_: thermocline depth (ambiguous — say anomaly or total), `depth`, `eta`

**Mean thermocline depth** (`H`):
The constant reference depth of the resting upper layer. A scenario parameter,
never a field.

**Upper layer**:
The single active layer of warm water in the 1.5-layer model, of thickness
`H + h`, sitting on a motionless abyss of infinite depth.
_Avoid_: mixed layer (a different physical concept), surface layer

**Reduced gravity** (`g'`):
The effective gravity `g·Δρ/ρ₀` felt at the density contrast across the
thermocline. Sets the wave speeds.

**Current anomaly** (`u`, `v`):
Zonal (`u`, eastward) and meridional (`v`, northward) velocity anomalies of the
upper layer.

**Rayleigh damping** (`r`):
Linear damping coefficient standing in for unresolved dissipation and mixing.
Gives waves a realistic decay and keeps the model stable.
_Avoid_: friction, viscosity (this is neither, physically)

## Forcing

**Alizés**:
The trade winds — the persistent easterly surface winds over the equatorial
Pacific that drive the whole system. Used interchangeably with "trade winds"
in prose; `alizes` (unaccented) in code identifiers.

**Wind stress** (`τx`, `τy`):
The externally specified surface stress field forcing the ocean, a function of
`(x, y, t)`. Easterly trade-wind stress is `τx < 0`. In v1 it is prescribed,
never derived from an atmosphere model.

**Wind forcing**:
A wind stress together with the field it is sampled into over one basin — what
a run holds and asks for the stress at each instant it integrates through, as
opposed to the stress *function* (`WindStress`) or the sampled field at one
instant (`WindStressField`). Whether a run re-samples for a given instant is
the wind's own declared time dependence
([ADR-0009](docs/planning/adr/0009-wind-declares-its-time-dependence.md)).

**Scenario**:
A complete, runnable specification of one simulation: grid, physical
parameters, wind-forcing description, and run length. The engine's unit of
input.
_Avoid_: config, experiment, case

**Seasonal cycle**:
The alizés breathing with the year: the steady wind-stress field scaled by an
annual harmonic, `1 + a·cos(2π(t − t_peak)/T_year)`, with `T_year` the tropical
year. It changes the winds' strength, not their shape or their direction.
_Avoid_: annual cycle, monsoon (a different phenomenon)

**Westerly wind burst**:
An idealized positive-`τx` anomaly (Gaussian in x, y and t) superimposed on the
trade winds, representing the class of perturbation known to trigger El Niño
onset.

## Waves and scales

**Kelvin wave**:
An equatorially trapped wave travelling **eastward only**, non-dispersive, at
speed `c`. The fast branch of the basin's adjustment.

**Rossby wave**:
An equatorially trapped wave travelling **westward**, dispersive; the gravest
meridional mode travels at `c/3`.

**Kelvin wave speed** (`c`):
`c = √(g'·H)`. The fastest signal in the model, and therefore the speed that
bounds the stable timestep.

**Equatorial deformation radius** (`Le`):
`Le = √(c/β)`. The meridional scale over which equatorial waves decay away from
the equator.

**Beta-plane** (`β`):
The linearization of the Coriolis parameter about the equator, `f = β·y`, with
`f = 0` exactly at the equator. `β ≈ 2.3×10⁻¹¹ m⁻¹s⁻¹`.

**Basin**:
The model domain approximating the equatorial Pacific, roughly 120°E–80°W by
25°S–25°N, closed on all four boundaries. Exact truncation is a scenario
parameter.

**Thermocline tilt**:
The steady-state east–west slope of `h` produced by sustained easterly wind
stress: deep in the west, shallow in the east. The observed mean state, and the
control case the model must reproduce.

## Coupled phenomena

**ENSO**:
El Niño–Southern Oscillation, the coupled ocean–atmosphere oscillation this
project exists to reproduce. **El Niño** is the warm phase (tilt collapses,
warm water surfaces in the east); **La Niña** the cold phase (tilt intensifies).

**SST anomaly** (`T'`):
Mixed-layer sea-surface temperature anomaly. Introduced only in the Epic 12
coupling extension, not part of the linear ocean core.

**Bjerknes feedback**:
The positive feedback loop that makes ENSO oscillate rather than merely
respond: weaker trade winds → flatter thermocline → warmer eastern SST →
weaker trade winds. Turning this loop on is what makes the model produce
*emergent* ENSO rather than a prescribed-wind response.

## Project vocabulary

**Ticket** (`T-<epic>.<n>`):
A merge-request-sized unit of work — the planning unit of this project.
Exactly one ticket = one GitHub issue = one pull request = one squashed commit
on `main`. Formerly written `MR-<epic>.<n>`; the frozen planning documents
under `docs/planning/` still use that older form.
_Avoid_: MR, merge request (this repo is on GitHub; "MR" invites GitLab
semantics), task, story

**Epic**:
A themed group of tickets with a shared goal and a position in the dependency
order. Represented as a GitHub milestone.
