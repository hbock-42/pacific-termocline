# Emergent-oscillation report

Whether the coupled model oscillates, and whether it oscillates at the period
theory predicts.

`docs/validation-report.md` is the ledger for Epic 07: the linear ocean core
against the closed-form consequences of its own equations. This is the ledger
for the Epic 12 stretch goal, and the question is a different one. Turning on
the Bjerknes feedback (T-12.1's SST equation, T-12.2's wind response) makes the
basin a *coupled system* whose behaviour nobody wrote down in advance — so the
only honest way to ask "does it do ENSO?" is to state what the delayed-
oscillator theory predicts, measure it, and compare.

The answer, up front, is **three predictions confirmed and one criterion
missed**:

| Prediction | Result |
|---|---|
| There is a critical feedback strength: below it a perturbation decays, above it a coupled mode grows | **Confirmed.** `0.06 < μ_c ≤ 0.08 Pa/K` |
| The growth rate rises with the feedback strength | **Confirmed.** −0.236, −0.193, −0.034 yr⁻¹ at μ = 0.02, 0.04, 0.06 Pa/K |
| The period is set by the basin wave-crossing time `L/c`, not by a thermal or damping timescale | **Confirmed.** Across a factor of 2 in `c` and 1.33 in `L`, `P·c/L` stays within 5% of 4.8 |
| The period falls in the observed ENSO band of roughly 2–7 years (the ticket's acceptance criterion) | **Not met.** `P = 1.03 years = 5.0·L/c` |

The missed criterion is escalated on
[issue #52](https://github.com/hbock-42/pacific-termocline/issues/52), labelled
`needs-human`: a green test suite is not a met acceptance criterion, and
whether to accept the negative result or to revisit the design decision
[Why the period is short](#why-the-period-is-short) points at is a human's
call rather than this ticket's.

The suite that produced every number below is
[`engine/tests/enso_oscillation.rs`](../engine/tests/enso_oscillation.rs).
[Why the period is short](#why-the-period-is-short) is the diagnosis, and
[What this does not establish](#what-this-does-not-establish) is the boundary
of the claim — which for a coupled model run at amplitudes this large is a
larger boundary than Epic 07's.

## Contents

- [What the delayed oscillator predicts](#what-the-delayed-oscillator-predicts)
- [The configuration](#the-configuration)
- [How a measurement is made](#how-a-measurement-is-made)
- [The threshold](#the-threshold)
- [The oscillation](#the-oscillation)
- [The period follows the basin crossing time](#the-period-follows-the-basin-crossing-time)
- [The period is short of the ENSO band](#the-period-is-short-of-the-enso-band)
- [The same claim, from what a run wrote to disk](#the-same-claim-from-what-a-run-wrote-to-disk)
- [Why the period is short](#why-the-period-is-short)
- [What this does not establish](#what-this-does-not-establish)
- [How these figures were produced](#how-these-figures-were-produced)

## What the delayed oscillator predicts

Suarez & Schopf (*J. Atmos. Sci.* 45, 1988) and Battisti & Hirst
(*J. Atmos. Sci.* 46, 1989) reduce the coupled equatorial basin to one delay
equation for an eastern SST anomaly `T`:

```text
dT/dt = a·T − b·T(t − δ)
```

- `a` is the instantaneous Bjerknes feedback (`CONTEXT.md`): warm east →
  weaker alizés → flatter thermocline → warmer east.
- `−b·T(t − δ)` is the same wind anomaly returning with the opposite sign,
  after an off-equatorial Rossby wave has carried it to the western wall and
  reflected back as an equatorial Kelvin wave.
- `δ` is the time that round trip takes. The gravest Rossby mode travels at
  `c/3` and the Kelvin wave at `c` — the two speeds
  [`docs/validation-report.md`](validation-report.md) already validated to five
  significant figures — so for a Rossby leg `L_R` and a Kelvin leg `L_K`,

  ```text
  δ = (3·L_R + L_K) / c
  ```

Substituting `T = e^{(σ + iω)t}` gives `σ + iω = a − b·e^{−(σ + iω)δ}`, and at
the margin `σ = 0`:

```text
cos(ωδ) = a/b            ω = b·sin(ωδ)
```

Three consequences, and they are what this report tests. None of them is a
number read out of the engine.

1. **A threshold.** Both `a` and `b` are proportional to the feedback strength
   `μ` of the T-12.2 wind response, while the damping that opposes them — the
   ocean's Rayleigh damping `r`, the mixed layer's thermal relaxation `ε_T`,
   the entrainment — is not. There is therefore a critical `μ` below which
   every perturbation decays.
2. **A period built out of `L/c`.** Every leg of the round trip is some
   fraction of the basin, so `δ ∝ L/c`; the phase condition fixes `ωδ`;
   therefore `P = 2π/ω ∝ L/c`, with a coefficient of order one. This is the
   falsifiable prediction: **change the basin's width or its wave speed and the
   period must follow, in proportion.** A period set by the mixed layer's
   125-day relaxation, by the ocean's 2.5-year damping, or by an artefact of
   the grid would not move at all.
3. **`P ≥ 4δ`, with equality only when `a = 0`.** `cos(ωδ) = a/b` with `a > 0`
   forces `ωδ < π/2`. So a measured period gives an *upper* bound on the
   model's effective delay — which is how [Why the period is
   short](#why-the-period-is-short) diagnoses the miss.

## The configuration

Every experiment is the same ocean, mixed layer and atmosphere; only the three
quantities the tests vary — feedback strength, basin width, reduced gravity —
differ, and each varies alone.

### The ocean

| symbol | value | what | source |
|---|---|---|---|
| `g'` | 0.05 m/s² | reduced gravity | Gill, *Atmosphere–Ocean Dynamics*, ch. 11; the value every Epic 07 suite runs at |
| `H` | 150 m | mean thermocline depth | as above |
| `β` | 2.3×10⁻¹¹ m⁻¹s⁻¹ | equatorial beta-plane gradient | `CONTEXT.md`, *Beta-plane* |
| `ρ₀` | 1025 kg/m³ | reference seawater density | `CONTEXT.md` |
| `1/r` | 2.5 years | Rayleigh damping timescale | Zebiak & Cane, *Mon. Wea. Rev.* 115, 1987, § 2a; Battisti, *J. Atmos. Sci.* 45, 1988 |

from which

```text
c  = √(g'·H) = 2.7386127875258306 m/s
Le = √(c/β)  = 345065.386842516 m   (345.07 km)
```

The damping is the one parameter that is *not* the Epic 07 value, and the
choice matters enough to say why. The steady-state validations of Epic 07 damp
on 100 days, because an undamped closed basin never stops ringing and has no
equilibrium to measure. A Rossby wave needs `3L/c ≈ 226 days` to cross this
basin: damping it on 100 days would remove the delayed branch of the feedback
before it returned, and the loop this report is about would not close at any
`μ`. That is not a tuning choice made to produce an oscillation — it is the
value the delayed-oscillator literature runs on, and it was fixed before any
of the runs below.

### The mixed layer and the atmosphere

| symbol | value | what | source |
|---|---|---|---|
| `H_m` | 50 m | mixed-layer depth | Zebiak & Cane 1987, § 2b |
| `r_s` | (2 days)⁻¹ | surface-layer Rayleigh drag | Zebiak & Cane 1987, § 2b (the engine's `DEFAULT_SURFACE_DRAG_PER_S`) |
| `∂T̄/∂x` | −4×10⁻⁷ K/m | mean zonal SST gradient | ≈ 7 K of warm-pool-to-cold-tongue contrast across this basin |
| `γ` | 0.1 K/m | `∂T_sub/∂h` | Zebiak & Cane 1987, § 2c |
| `ε_T` | (125 days)⁻¹ | thermal damping of `T'` | Zebiak & Cane 1987, § 2b |
| `L_a` | 2.3×10⁶ m | meridional scale of the wind response | the atmospheric equatorial Rossby radius, T-12.2's `DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M` |
| `τx` | −0.05 Pa | steady alizés | the mean easterly stress of the equatorial Pacific |
| `μ` | swept | feedback strength | the parameter this report characterises |

The mixed-layer numbers are the same set the T-12.1 and T-12.2 suites already
run at; nothing here re-tuned them.

### The basin and the numerics

| | |
|---|---|
| Basin | 120°E–80°W by 25°S–25°N (`CONTEXT.md`, *Basin*), `L = 17 791 km` |
| Resolution | 2°, so 80 × 25 cells and `Δx = Δy = 222.4 km` |
| Timestep | 30 000 s |
| `L/c` | 0.20586 tropical years (75.2 days) |
| `3L/c` | 0.61759 tropical years (226 days) |

Two degrees is coarse beside the Epic 07 wave suites, which resolve `Le` seven
to fourteen times over: here `Δy/Le = 0.645`. That is deliberate. What is
measured is a basin-crossing *time* — an integral over the whole waveguide —
rather than the shape of one wave's meridional profile, and every comparison
is made between two runs on the same grid. The tolerance on the scaling test
is built from what that costs.

The timestep is bounded by rotation and not by the gravity waves: `|f| = β·y`
reaches 6.4×10⁻⁵ s⁻¹ at 25°, whose inertial oscillation RK4 can follow up to
35 400 s ([ADR-0007](planning/adr/0007-rotation-timestep-bound.md)), where the
CFL bound of a 222.4 km cell at `c = 2.74 m/s` is 65 000 s. Cutting the
timestep to 12 000 s moved the measured period from 1.031339 to 1.031315
years — two parts in 10⁵ — so nothing below is timestep-limited.

## How a measurement is made

One experiment is one run of the coupled solver:

1. **Spin up.** The basin is held under the steady alizés for 20 years — eight
   `e`-foldings of the 2.5-year damping — and the eastern SST index is read
   off the state it settles to.
2. **Perturb.** A uniform 1 K warm anomaly, the scale of an observed El Niño,
   is added to `T'`.
3. **Record.** The index is followed for 44 more years as its *departure* from
   the value the spin-up left it at, sampled 2000 times per 32 years.

The index is `T'` averaged over the eastern third of the basin within 5° of the
equator — for the reference basin, 227°E–280°E, which is most of the observed
Niño-3 region (150°W–90°W). It is stated as a fraction of the basin rather than
as two longitudes so that it follows the basin when the scaling test makes the
basin narrower.

Three numbers are then read off the record:

- **Growth rate**, from the whole record: the log ratio of the r.m.s. departure
  over its second half to its first, per year. For `e^{σt}·cos(ωt + φ)` sampled
  over many cycles this is `σ`, and for a saturated cycle it is zero.
- **Amplitude**, from the last 32 years: half the peak-to-trough range, once
  the 1 K step has settled out. The step is every mode of the basin at once;
  what survives 12 years of settling is the slowest.
- **Period**, from the same 32 years: the discrete Fourier power evaluated on a
  frequency grid four times finer than the record's own resolution, its largest
  peak refined by a parabola through it and its two neighbours. Only periods
  between 0.2 and 8 years are considered — 8 years is past the slow edge of the
  observed ENSO band, so a model that *did* oscillate there would be seen to.

## The threshold

*`the_growth_rate_rises_with_the_feedback_strength`,
`the_perturbation_decays_below_the_threshold_and_grows_above_it`,
`an_open_loop_run_returns_to_the_state_the_trades_hold_it_at`*

| `μ` (Pa/K) | growth rate (yr⁻¹) | settled amplitude (K) | period (yr) |
|---|---|---|---|
| 0.00 (open loop) | −0.2480 | 9.55×10⁻⁷ | — |
| 0.02 | −0.2356 | 2.18×10⁻³ | — |
| 0.04 | −0.1932 | 4.20×10⁻² | 0.9632 |
| 0.06 | −0.0337 | 7.14×10⁻¹ | 0.9845 |
| 0.08 | +0.0004 | 5.787 | 1.0313 |

The growth rate climbs monotonically with the feedback strength, is negative at
every `μ ≤ 0.06 Pa/K`, and by `μ = 0.08 Pa/K` the run has stopped growing
because it has *saturated* — see [the next section](#the-oscillation). So

```text
0.06 Pa/K < μ_c ≤ 0.08 Pa/K
```

At `μ = 0` — T-12.1's model, in which `T'` reads the ocean and nothing reads
`T'` — the 1 K anomaly is gone: 9.5×10⁻⁷ K survives into the window, which is
what a relaxation no slower than the ocean's 2.5-year damping leaves after
twelve years. There is nothing for the loop to oscillate with, because `a` and
`b` are both zero.

No tolerance appears in either assertion. The monotonicity test compares three
numbers for order; the threshold test compares an amplitude against the 1 K
perturbation the run was given, and a run that ends smaller than it started
decayed.

**The two weakest rows have no period.** Their records are 10⁻³ K and below —
four to seven orders under the perturbation — and what a spectrum finds there
is the residue of the spin-up rather than a coupled mode. The dash records
that, rather than a number that would look like a 4-year oscillation and be
nothing of the kind.

## The oscillation

*`the_supercritical_run_settles_into_a_self_sustained_oscillation`*

Above the threshold the linear mode grows, and then stops: the coupled model
has one nonlinearity — the `w⁺ = max(w, 0)` clamp of
[`engine/src/sst.rs`](../engine/src/sst.rs), which is quadratic once the wind
depends on the state — and it bounds the cycle. Over the 32-year settled window
of the `μ = 0.08 Pa/K` run:

| | |
|---|---|
| largest departure, first half | 11.2765 K |
| largest departure, second half | 11.2627 K |
| ratio | 0.9988 |
| peak-to-trough swing | 11.5747 K (half-range 5.7873 K) |

That is a **limit cycle**, not a runaway: an unsaturated mode at this strength
would have grown by orders of magnitude across the same window. The waveform is
strongly asymmetric — a long warm phase and a short sharp cold one — which is
the shape a relaxation oscillator has and the shape the observed ENSO index
has.

The amplitude itself is *not* a result; see
[What this does not establish](#what-this-does-not-establish).

## The period follows the basin crossing time

*`the_period_scales_with_the_basin_crossing_time`*

The delayed oscillator's central claim. Three oceans are compared against the
reference at the same `μ = 0.08 Pa/K`, each changing one thing:

| ocean | `L` (km) | `c` (m/s) | `L/c` (yr) | period (yr) | `P·c/L` |
|---|---|---|---|---|---|
| reference | 17 791 | 2.7386 | 0.20586 | 1.03134 | 5.010 |
| narrowed to 120° of longitude | 13 343 | 2.7386 | 0.15440 | 0.73639 | 4.769 |
| `g' = 0.10`, so `c·√2` | 17 791 | 3.8730 | 0.14557 | 0.66580 | 4.574 |
| `g' = 0.20`, so `2c` | 17 791 | 5.4772 | 0.10293 | 0.47443 | 4.609 |

and as ratios against the reference, which is the form the test asserts on:

| ocean | predicted `(L/c)` ratio | measured period ratio | relative departure |
|---|---|---|---|
| narrowed to 120° | 0.7500 | 0.7141 | 4.8% |
| `c·√2` | 0.7071 | 0.6456 | 8.7% |
| `2c` | 0.5000 | 0.4600 | 8.0% |

Every other parameter — the mixed layer's 125-day relaxation, the ocean's
2.5-year damping, the atmosphere's width, the trades, the grid — is held fixed,
so the three alternatives this discriminates against (a period set by a thermal
time, by a damping time, or by the grid) all predict a ratio of exactly one.
The measured ratios move by up to a factor of two, in proportion to `L/c`, and
`P·c/L` stays inside ±5% of 4.8 across the set.

**The tolerance is 25% on the ratio, and it is a discrimination margin rather
than a precision claim.** Two reasons it is not tighter:

- `P = 4δ` holds only in the purely-delayed limit. `cos(ωδ) = a/b` and `a/b` is
  *not* held fixed when the basin or the wave speed changes, so the theory
  gives a proportionality with an order-one coefficient rather than an
  equality. The 4.6–5.0 spread in the last column is that effect, and it is
  physics rather than error.
- The largest numerical term does not cancel between two oceans. The meridional
  truncation of the Epic 07 budget, `(2m+1)·(Δy/Le)²` with `m = 0`, is 0.42 on
  this grid at `c = 2.74 m/s` and 0.21 at `c = 5.48 m/s`, because `Le = √(c/β)`
  grows with the wave speed. The 0.21 difference between them is itself of the
  order of the band.

What the band has to do is separate the prediction from its alternatives. The
smallest change it is asked to see is a factor of 0.71 — nearly three times the
band away from the ratio of one the alternatives predict — and the largest
measured departure is 8.7%, a third of it.

## The period is short of the ENSO band

*`the_period_falls_short_of_the_observed_enso_band`*

```text
measured period = 1.0313 tropical years = 5.010 · L/c
observed ENSO   = roughly 2 to 7 years
```

The ticket's acceptance criterion asks for a period in the observed range. It
is not met, by a factor of two at the band's fast edge. The suite records the
miss as an assertion rather than as prose — `period < 2 years` — so that
physics which lengthens the period into the band fails the test loudly instead
of quietly agreeing with a document nobody re-read. Both bounds in that
assertion come from observation and not from this model: two years is the fast
edge of the observed band, and `L/c` is the shortest timescale a delayed
oscillator can build a period out of at all.

## The same claim, from what a run wrote to disk

*`the_written_run_carries_the_oscillating_sst_index`*

The six results above are read off states a test holds in memory. T-05.4 put
`T'` in the frame format so that an SST index could be read back off *disk*
instead, and one test does exactly that: a coupled scenario in the config
format of [`docs/scenario-config-reference.md`](scenario-config-reference.md)
with `wind_feedback_strength_pa_per_k = 0.08`, driven through `run_scenario`
into a run directory, then reopened with `RunReader` and reduced to the same
eastern index frame by frame.

It differs from the experiments above in one way, deliberately: **nothing
perturbs it.** It starts at rest and runs for 45 years, so the coupled mode has
to grow out of the switch-on of the alizés — which is the scenario a user would
actually write. Over the last third of that run the index crosses its own mean
at least the four times the test asks for, which is two whole cycles.

The assertion is a crossing count rather than a period, because the frames are
saved every 200 steps (0.19 years) to keep the run directory small, and five
samples a cycle counts crossings honestly but does not locate a spectral peak.
What this test establishes is that the oscillation is in the *file*, not only
in the solver.

## Why the period is short

The measurement bounds the model's effective delay. `P ≥ 4δ`, so

```text
δ_eff ≤ P/4 = 0.2578 years = 94 days = 1.25 · L/c
```

and therefore `3·L_R + L_K ≤ 1.25·L`. Compare the classical configuration: a
wind anomaly centred in the middle of the basin and an index in the eastern
third gives `L_R = L/2` and `L_K = 5L/6`, hence

```text
δ_classical = (3·½ + ⅚)·L/c = 2.33·L/c = 175 days      4δ = 1.92 years
```

— which is inside the observed band at its fast edge. **The model's effective
delay is about half the classical one**, and the reason is a design decision of
T-12.2 recorded in
[ADR-0010](planning/adr/0010-wind-response-is-diagnosed-per-stage.md): the
statistical wind response is *zonally uniform*. It answers a basin-wide SST
index with a stress anomaly at every longitude, so Rossby waves are launched
everywhere — including hard against the western wall, where they reflect almost
immediately. The reflected Kelvin signal the eastern box sees is a
superposition over every source longitude, with delays running from nearly zero
to `4L/c`, and the short paths are both more numerous in the early response and
less damped. A single delay of `2.33·L/c` is what a *localised* central-Pacific
wind patch would give; a uniform one gives the shorter ensemble the model
measures.

That is the first candidate and the one the arithmetic points at. Three more,
in the order they would be worth trying:

1. **No zonal structure in the wind response at all.** Gill's (1980) solution
   puts the westerly anomaly *west* of the heating and an easterly *east* of
   it. T-12.2 deliberately carries none of that, and ADR-0010 says why: a
   Gill pattern needs the atmospheric model this project does not solve.
   Giving the response a centre and a width — even a prescribed one — is the
   smallest change that would lengthen `L_R`.
2. **A basin-wide index, not an eastern one.** The wind answers `⟨T'⟩` over the
   whole basin (T-12.2), which blurs the source of the delayed signal in the
   same direction.
3. **`a/b` is far from one.** The observed 4-year ENSO period needs
   `δ ≈ 365 days ≈ 4.85·L/c`, which is *longer than the largest round trip the
   basin has* (`4L/c`, a Rossby wave the full width plus a Kelvin wave back).
   So the real system cannot be at `P = 4δ` either: it must sit close to the
   phase condition, with the instantaneous feedback nearly balancing the
   delayed one. Nothing in this model was tuned to put it there, and the period
   moves only from 4.68 to 5.01 `L/c` across the whole `μ` sweep, so it is not
   near that limit.

None of these is in scope for T-12.3, and none should be attempted inside it:
the first two contradict a decision recorded in an ADR, which is a human's to
revisit.

## What this does not establish

Stating the boundary is part of the claim, and here it is a wide one.

- **The amplitude is not a result.** The reference cycle swings 11.57 K peak to
  trough — a half-range of 5.79 K, which is the "amplitude" column of the
  tables above — and reaches 11.28 K *below* the value the spin-up settled the
  eastern index at. Set against the 2–4 K a strong observed El Niño reaches in
  Niño-3, that is several times too large, and it
  is bounded by the `w⁺` clamp — a *kinematic* switch on the sign of the
  upwelling — rather than by any heat budget. The model has no climatological
  mean state for an anomaly to be small against, so nothing in it limits `T'`
  physically, and the spun-up state it swings about is itself a multi-kelvin
  departure. Read the period, the threshold and the scaling; do not read the
  amplitude.
- **The runs are far outside the linear core's validated range.** Epic 07
  validated the linear equations, and `CODING_STANDARDS.md` § *Scope guards*
  keeps the v1 core linear. A ±11 K SST anomaly and the thermocline excursion
  that goes with it are amplitudes at which the neglected nonlinear advection
  would matter in a real ocean. What is validated here is that *these
  equations* oscillate as delayed-oscillator theory says they must — not that
  the equatorial Pacific does it for the same reason.
- **One damping value.** The oscillation exists at `1/r = 2.5 years` and does
  not exist at all at the 100-day damping of the Epic 07 suites. The threshold
  and the period were not mapped as functions of `r`.
- **One resolution, for the period's absolute value.** Coarsening the grid from
  2° to 2.5° moves the period from 1.0313 to 1.0569 years, 2.5%. The scaling
  result is a comparison at fixed resolution and does not inherit that; the
  absolute `5.0·L/c` carries it. The *amplitude* over the same change goes from
  5.8 K to 58 K — a further reason to read only the period, the threshold and
  the scaling out of these runs.
- **No comparison against observations, except the band itself.** Every other
  number here is checked against the model's own theory. The 2–7 year band is
  the one place an observed quantity enters, and it is the one place the model
  fails.

## How these figures were produced

| | |
|---|---|
| Suite | `engine/tests/enso_oscillation.rs` |
| Command | `cargo test --workspace` — the command CI's gate runs |
| Toolchain | rustc 1.90.0, aarch64-apple-darwin |

Every figure above came out of the runs the assertions read. A passing Rust
test is silent, so they were read out by adding a temporary reporting test to
the file, which recomputes nothing — it calls the same `Experiment::record`,
`dominant_period_years`, `amplitude_k` and `growth_rate_per_year` the
assertions call, on the same experiments — and prints them. That
instrumentation was reverted before this document was committed; it changed no
assertion, no tolerance and no run. What is committed is therefore the *suite*,
not the tables: re-deriving a figure below means adding that reporting test
back, which is the same convention
[`docs/validation-report.md`](validation-report.md) records for Epic 07.

The runs are deterministic (`CODING_STANDARDS.md` § *Correctness and failure*):
no randomness, no iteration-order dependence. As in
[`docs/validation-report.md`](validation-report.md), a CI run on x86-64 Linux
may differ from these figures in the last bit or two, since the two platforms'
`exp` and `sqrt` do not round identically. Nothing asserted is decided closer
than the third significant figure, and the tightest margin above — the 8.7%
departure against a 25% band — has a factor of nearly three in hand.
