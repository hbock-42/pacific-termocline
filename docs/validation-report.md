# Validation report

How we know this simulation is scientifically correct.

The engine solves the linear 1.5-layer shallow-water equations on an
equatorial beta-plane ([ADR-0003](planning/adr/0003-numerical-scheme.md)).
That system has a small number of *closed-form* consequences — wave speeds, a
dispersion relation, a trapping scale, a steady wind-driven balance, an energy
invariant — and Epic 07 is the suite that holds the code to them. This document
is the ledger: for each check, what was compared, what theory predicts, what
the run measured, and what tolerance the comparison was made against and why.

Every number below the line "Measured" was produced by the run described in
[How these figures were produced](#how-these-figures-were-produced). None of
them is copied from a pull-request description, and none of them is a
placeholder. Where the tolerance is a sum of derived terms, the terms are shown
and not just the total.

## Contents

- [How these figures were produced](#how-these-figures-were-produced)
- [The ocean every test runs in](#the-ocean-every-test-runs-in)
- [How a tolerance is built](#how-a-tolerance-is-built)
- [Summary](#summary)
- [T-07.1 — Kelvin wave speed and non-dispersion](#t-071--kelvin-wave-speed-and-non-dispersion)
- [T-07.2 — Rossby dispersion and `c/3` propagation](#t-072--rossby-dispersion-and-c3-propagation)
- [T-07.3 — Equatorial deformation radius](#t-073--equatorial-deformation-radius)
- [T-07.4 — Steady wind-driven thermocline tilt](#t-074--steady-wind-driven-thermocline-tilt)
- [T-07.5 — Conservation in the undamped, unforced limit](#t-075--conservation-in-the-undamped-unforced-limit)
- [What this suite does not validate](#what-this-suite-does-not-validate)
- [Reproducing these figures](#reproducing-these-figures)

## How these figures were produced

| | |
|---|---|
| Commit | `e9579bf` (`main` at T-07.4, the last of the five suites to merge) |
| Command | `cargo test --workspace --all-features` — the command CI's `ci` job runs |
| Toolchain | rustc 1.90.0 (1159e78c4 2025-09-14) |
| Machine | Apple M1 Pro, aarch64-apple-darwin |
| Result | 41 test binaries, 334 tests, **0 failures** |

The five suites do not print their measurements when they pass — a passing
Rust test is silent, and each figure below lives inside the assertion that
consumes it. They were read out by adding a temporary reporting test to each of
the five files, which recomputes nothing: it reads the same
`measured_speed_m_per_s`, `fit`, `departure_from` and `worst_energy_drift`
the assertions read, off the same shared runs, and prints them. That
instrumentation was reverted before this document was committed; it changed no
assertion, no tolerance and no run.

These are figures from a real run, but — to be exact about it — from a real
*local* run, not from the CI job the acceptance criterion names. No CI artifact
carries them, because a passing suite prints nothing and the instrumentation
that read them out was reverted. Two properties are what make them
reproducible rather than anecdotal:

- **Debug and release agree bit for bit.** The same fifty measurement lines
  were taken under `cargo test` and `cargo test --release` and compared: every
  digit matches. That rules out the optimisation profile — CI's is debug — as
  anything that could move a digit.
- **The runs are deterministic** (CODING_STANDARDS.md § *Correctness and
  failure*): no randomness, no iteration-order dependence, one thread per run.

What those two do not rule out is the platform, which is the next paragraph.

CI runs `ubuntu-latest` on x86-64; this run was aarch64-darwin, and the two
platforms' libm implementations round `exp`, `sinh` and `asinh` to within an
ulp of each other but not identically. A CI run may therefore differ from the
figures below in the last bit or two. That is far below anything asserted:
every comparison here is decided at five significant figures or fewer. The
one comparison whose margin is not orders of magnitude — T-07.4's channel
profile, at 1.0006× — is decided at the fifth figure of a quantity whose
prediction and measurement were computed on the same machine, so a change of
libm moves both together. Where a number is printed to full
`f64` width below — the derived scales `c` and `Le`, the raw measured speeds,
the T-07.4 ratio that agrees with its prediction to ten figures — its last two
or three digits are a statement about this machine and not a portable one.

## The ocean every test runs in

All five suites use the same published equatorial-Pacific parameters — Gill,
*Atmosphere–Ocean Dynamics*, ch. 11; Cane & Sarachik 1981. Four of them read
them from `engine/tests/support/mod.rs`, the shared equatorial-wave kit;
`conservation.rs` predates that module and declares the same five values
itself, which is worth knowing when reading the two side by side:

| symbol | value | what |
|---|---|---|
| `g'` | 0.05 m/s² | reduced gravity, first baroclinic mode |
| `H` | 150 m | mean thermocline depth |
| `β` | 2.3×10⁻¹¹ m⁻¹s⁻¹ | equatorial beta-plane gradient |
| `ρ₀` | 1025 kg/m³ | reference seawater density |

from which the two derived scales every prediction is written in terms of are

```text
c  = √(g'·H) = 2.7386127875258306 m/s
Le = √(c/β)  = 345065.386842516 m   (345.07 km)
```

Both are computed in the tests from the definitions in `CONTEXT.md`, never
asked of the engine: an expected value read out of the code under test agrees
with it by construction.

`deformation_radius.rs` additionally runs two *counterfactual* oceans —
`g'` quadrupled and `β` doubled — for the reason given in
[T-07.3](#t-073--equatorial-deformation-radius).

## How a tolerance is built

No tolerance in this suite is a round number chosen because it passes. Each is
assembled from named error terms of its own configuration, and each is
therefore a **function of the resolution**, so that a finer grid is held to a
tighter bound rather than to the same one. The recurring entries:

- **Meridional truncation, `(2m+1)·(Δy/Le)²`.** The C-grid operators are second
  order and the structure they differentiate is `ψₘ`, which oscillates on
  `Le/√(2m+1)`. `m = 0` for the Kelvin wave, `m = 2` for the gravest Rossby
  mode. This is the dominant term almost everywhere.
- **Zonal truncation, `(Δx/σ)²/4`.** The centred difference reads a group speed
  low by `(kΔx)²/2`, and a Gaussian packet of width `σ` has `⟨k²⟩ = 1/(2σ²)`.
- **Stray energy, `⟨k̂²⟩ = 1/(2σ̂²)`.** The Rossby initial condition is the
  *long-wave* mode, exact only as `k̂ → 0`. The Kelvin branch is an exact
  solution at every wavenumber, so its share of this term is zero rather than
  small.
- **Equilibrium, `(e^{−rT} + e^{−rT/2})·√N`.** What a finite spin-up leaves of
  the initial transient (T-07.4 only).
- **A factor of two** (`TRUNCATION_SAFETY`) on the assembled budget — the
  stray-energy term included, not only the truncation ones — as the standing
  allowance for the leading `O(1)` coefficients none of them evaluates. It
  multiplies coarse and fine budgets alike, so it widens what a point check
  admits without touching what a convergence test measures.

Each ticket's own module header is the **normative** derivation of its budget:
the terms are declared there, next to the constants they are evaluated from,
and the code computes the tolerance from those constants rather than from
anything written here. What this document adds is the measurement beside them.
Where the two ever disagree, the module header and the code are right and this
file is stale — and the fastest way to find that out is that the arithmetic
below stops adding up.

One thing about those budgets is the structure of the whole suite, and is worth
stating plainly: **the point checks are bounds, and are generous by design.** A
bound is passed by an error of any size below it, so on its own it says little.
What ties each measured error to the discretisation is its *rate*: every ticket
also asserts that halving the cells divides the error by about four.

In T-07.1 to T-07.3 that order is bounded on **both** sides, at `[1.5, 3]`. Too
small an order is a scheme that is not second order; too large an order means
the fine error is no longer the truncation, so the point check's budget has
stopped describing the run. T-07.5's check is deliberately one-sided (`≥ 1.8`):
there the two terms of the bound fall at different rates and the run length is
fixed rather than the resolution, so falling *faster* than second order is not
a symptom of anything.

That is CODING_STANDARDS.md § *Convergence over point checks*, and it is why
the margins below — often three or four orders of magnitude — are not evidence
that the tolerances are slack.

## Summary

| ticket | claim | predicted | measured | budget | margin |
|---|---|---|---|---|---|
| T-07.1 | Kelvin pulse travels east at `c` | 2.7386 m/s | 2.7356 m/s (coarse) / 2.7379 (fine) | 5.09% / 1.27% | 46× / 46× |
| T-07.1 | error is second order in the cell size | 2 | **1.997** | `[1.5, 3]` | — |
| T-07.2 | gravest Rossby packet travels west at `c/3` | −0.90676 m/s | −0.90263 m/s (coarse) / −0.90559 (fine) | 21.4% / 5.38% | 47× / 42× |
| T-07.2 | error is second order in the cell size | 2 | **1.822** | `[1.5, 3]` | — |
| T-07.3 | both waves decay meridionally on `Le` | 345.07 km | 345.06 km (Kelvin) / 345.04 km (Rossby) | 5.09% / 7.17% | 1800× / 1000× |
| T-07.4 | steady tilt matches the damped closed form | `h(x)` sinh profile | 2.5076×10⁻⁵ of the tilt | 2.5092×10⁻⁵ | 1.0006× |
| T-07.5 | undamped energy drift stays within the derived bound | ≤ 1.2741×10⁻² | 1.1527×10⁻³ (16 rows) | 1.2741×10⁻² | 11.1× |

"Margin" is the tolerance divided by the measured error. Read it with the
paragraph above: a large margin means the *bound* is loose, and the
second-order rates are what say the error is the discretisation's.

---

## T-07.1 — Kelvin wave speed and non-dispersion

**Suite:** `engine/tests/kelvin_wave_propagation.rs` · **Issue:** #29 ·
**Claim:** `CONTEXT.md`, *Kelvin wave* — an equatorially trapped first-mode
disturbance travels **eastward only**, at `c = √(g'·H)`, **without dispersion**.

### Configuration

A Gaussian Kelvin pulse of amplitude 10 m and zonal width `σ = 1500 km`
(4.3 `Le`), launched 4 `σ` from the western wall of a closed, unforced,
undamped 20 000 km × 4000 km basin, and sampled at 0.5 and 5.0 transits of
that reference width — the reference rather than each run's own, so that the
narrower pulse of the non-dispersion test below shares the same clock and the
same flight. Two resolutions: 100 × 80 cells (Δx = 200 km, Δy = 50 km) and
200 × 160 (Δx = 100 km, Δy = 25 km). The pulse ends the flight 4 `σ` short of
the eastern wall, so no boundary takes part in what is measured.

The position is the energy-weighted zonal centroid of the eastward invariant's
`ψ₀` projection, `P₀[u/c + h/H]`, which no Rossby mode contributes to; the
speed is the displacement between the two samples over the steps' own elapsed
time.

### The speed — `the_kelvin_pulse_travels_east_at_the_reduced_gravity_wave_speed`

The Kelvin branch obeys `∂r/∂t + c·∂r/∂x = 0` exactly, at every wavenumber, so
the theoretical centroid speed is `c` with no physical bias term at all.

| | coarse (100×80) | fine (200×160) |
|---|---|---|
| predicted `c` | 2.7386127875258306 m/s | 2.7386127875258306 m/s |
| **measured** | **2.7355778987932604 m/s** | **2.7378525591688403 m/s** |
| error | 0.11082% | 0.027760% |
| budget | 5.0881% | 1.2720% |

Both are positive, which is the "eastward only" half of the claim asserted
directly.

**Where the budget comes from.** Two terms, both second order in the cell size:

| term | coarse | fine |
|---|---|---|
| meridional truncation `(Δy/Le)²` | 2.0996% | 0.52490% |
| zonal group-speed truncation `(Δx/σ)²/4` | 0.44444% | 0.11111% |
| × `TRUNCATION_SAFETY` = 2 | **5.0881%** | **1.2720%** |

Two further terms were derived and dropped as negligible: the RK4 phase error
`(ω·Δt)⁴/120 ≈ 1×10⁻⁸`, and centroid contamination from shed modes and clipped
tails `≈ 4×10⁻⁴`. They are five and one orders below the terms carried.

### The error is second order — `the_kelvin_speed_error_shrinks_at_the_schemes_second_order`

Halving both cell dimensions took the error from 0.11082% to 0.027760%:

> **measured convergence order 1.9971**, against the `[1.5, 3]` a second-order
> scheme owes.

This is the acceptance criterion's "the error shrinks with resolution, not a
fixed offset", and it is what makes the 46× margins above evidence about the
scheme rather than about the bound.

### Eastward only — `the_kelvin_pulse_carries_no_westward_energy`

The Kelvin branch has `q = u/c − h/H ≡ 0` identically, so `P₀[q]` must stay
empty for the whole flight.

| | westward share of the energy |
|---|---|
| early sample (t = 280 548 s) | 1.3237×10⁻⁶ |
| late sample (t = 2 745 367 s) | 1.1665×10⁻⁶ |
| ceiling `2·(Δy/Le)⁴` | 8.8167×10⁻⁴ |

The share *falls* over the flight rather than accumulating: nothing turned
round and nothing split. 660× inside the ceiling.

### Non-dispersion, as shape — `the_kelvin_pulse_keeps_its_zonal_shape`

RMS zonal width of the projected profile, measured in a window travelling with
the packet:

| | |
|---|---|
| early | 1 061 613.67 m |
| late | 1 061 680.33 m |
| **growth** | **0.0062793%** |
| bound | 0.89227% |

The bound is two derived terms times two: numerical dispersion
`(c·t/σ)²·(Δx/σ)⁴/16 = 4.9×10⁻⁴`, and shed modes at the window's lever arm
`(Δy/Le)⁴·(W/σ)² = 4.0×10⁻³`. The measured growth is 142× inside it.

### Non-dispersion, as spectral independence — `the_kelvin_speed_does_not_depend_on_the_packets_zonal_width`

The sharper statement: every wavenumber travels at `c`, so doubling the
packet's spectral content must not move its speed.

| | |
|---|---|
| wide pulse (`σ` = 1500 km) | 2.7355778987932604 m/s |
| narrow pulse (`σ` = 750 km) | 2.7265691609825990 m/s |
| difference | 0.32895% of `c` |
| bound (the two runs' zonal budgets) | 4.4444% |

This assertion discriminates. Applied to the *dispersive* gravest Rossby mode,
whose long-wave speed carries the bias `(4/9)·(Le/σ)²`, the same halving would
move the measured speed by 7.1% — outside this budget — and T-07.2 measures
exactly that effect on that branch (1.98%, at a width where the budget is
1.05%).

### Equatorial trapping — `the_kelvin_pulse_stays_on_the_equatorial_waveguide`

A Kelvin wave's thermocline anomaly is `ψ₀` and nothing else. At the end of the
flight the `ψ₂`/`ψ₀` cross-ratio of `h` was **1.7816×10⁻⁴**, against a ceiling
of `2·(Δy/Le)² =` 4.1992×10⁻²: 236× inside.

---

## T-07.2 — Rossby dispersion and `c/3` propagation

**Suite:** `engine/tests/rossby_wave_dispersion.rs` · **Issue:** #30 ·
**Claim:** `CONTEXT.md`, *Rossby wave* — an equatorial Rossby disturbance
travels **westward**, and the gravest meridional mode's long-wave speed is
`c/3`.

### Configuration

A gravest-mode (`n = 1`) Rossby packet of amplitude 10 m and width
`σ = 2800 km` (8.114 `Le`), including the `O(k̂)` meridional velocity of the
exact mode, launched 4 `σ` from the eastern wall of a 32 000 km × 4000 km
basin and sampled at 0.5 and 3.0 transits. Resolutions 128 × 80
(Δx = 250 km, Δy = 50 km) and 256 × 160.

The position is the energy-weighted zonal centroid of `P₀[u/c − h/H]`, the
channel the Kelvin branch contributes nothing to.

### The prediction is not `c/3`, and that is the point

This branch is **dispersive**, so a packet of finite zonal width travels
measurably slower than the long-wave limit for a reason that has nothing to do
with the grid. The expected speed is the energy-weighted mean group velocity of
the packet's own Gaussian spectrum,

```text
⟨c_g⟩ = ∫ c_g(k̂)·e^{−k̂²σ̂²} dk̂ / ∫ e^{−k̂²σ̂²} dk̂,
```

with `c_g` obtained by implicit differentiation of the equatorial dispersion
relation `ω̂³ − (k̂² + 3)·ω̂ − k̂ = 0` (Rossby root by Newton, quadrature by
Simpson). No expansion in `k̂`, and no number from the engine.

| | |
|---|---|
| long-wave `c/3` | −0.9128709291752769 m/s |
| dispersion relation, this packet | −0.9067608188421247 m/s |
| **predicted slowdown below `c/3`** | **0.66933%** |

### Speed and dispersion relation — `..._travels_west_at_a_third_of_the_kelvin_wave_speed`, `the_rossby_speed_matches_the_equatorial_dispersion_relation`

| | coarse (128×80) | fine (256×160) |
|---|---|---|
| **measured speed** | **−0.9026323164108231 m/s** | **−0.9055931345942254 m/s** |
| error vs the dispersion relation | 0.45530% | 0.12878% |
| budget | 21.421% | 5.3750% |
| measured slowdown below `c/3` | 1.1216% | 0.79724% |
| predicted slowdown | 0.66933% | 0.66933% |

Both speeds are negative — westward, as claimed. Both are also *nearer* `c/3`
than either neighbouring analytic speed, which is what identifies the mode
rather than merely a westward signal: the fine run misses `c/3 = 0.91287 m/s`
by 0.00728 m/s, where the `n = 2` mode's `c/5 = 0.54772 m/s` is 0.358 m/s away
and the Kelvin speed `c` is 1.833 m/s away.

**Where the budget comes from:**

| term | coarse | fine |
|---|---|---|
| meridional truncation `5·(Δy/Le)²` | 10.498% | 2.6245% |
| zonal group-speed truncation `(Δx/σ)²/4` | 0.19930% | 0.049825% |
| stray energy `⟨k̂²⟩²` weighed over the flight | 0.013181% | 0.013181% |
| × `TRUNCATION_SAFETY` = 2 | **21.421%** | **5.3750%** |

Dropped as negligible: RK4 phase error ≈ 2×10⁻¹¹, wall clipping ≈ `e^{−16}`.

The meridional term dominates because `ψ₂` oscillates on `Le/√5`, and it is
what makes this budget five times T-07.1's at the same `Δy`. Note that the
stray-energy term does *not* shrink with the cells — it is the floor a
convergence rate is read against, and it is why the reference packet is 8.1
`Le` wide.

### The error is second order — `the_rossby_speed_error_shrinks_at_the_schemes_second_order`

Halving both cell dimensions took the error from 0.45530% to 0.12878%:

> **measured convergence order 1.822**, against `[1.5, 3]`.

The rate is read against the *dispersion relation*, not against `c/3`: the gap
to `c/3` is the packet's own physics and would survive any refinement, so a
rate measured against it would be reading a floor. The order sits below
T-07.1's 1.997 because the resolution-independent stray-energy term is a larger
share of the fine run's error here.

### The dispersion itself — `narrowing_the_packet_slows_it_by_what_the_dispersion_relation_says`

The mirror image of T-07.1's non-dispersion test, and the strongest statement
in this file: narrow the packet by two and it must slow by an amount the
dispersion relation names *in advance*.

| | |
|---|---|
| predicted slowdown, `σ` = 2800 km | 0.66933% of `c/3` |
| predicted slowdown, `σ` = 1400 km | 2.6106% of `c/3` |
| **predicted difference** | **1.9413%** |
| **measured difference** | **1.9769%** |
| residual | 0.035664% |
| bound | 1.0483% |

(Measured speeds: −0.9055931345942254 m/s wide, −0.8875463285988355 m/s
narrow, both at the fine resolution.) The bound is *half* the effect it
measures, so a non-dispersive branch fails this test rather than passing it by
default. The residual is 29× inside the bound.

The bound's terms are the residue of the meridional truncation
(`(Δy/Le)²·5 ×` the predicted difference, since it scales both speeds by a
common relative factor), the two runs' zonal truncations, and the two
stray-energy terms — all × 2.

### Westward only — `the_rossby_packet_carries_no_eastward_energy`

| | eastward share | ceiling |
|---|---|---|
| coarse, early / late | 5.9066×10⁻⁷ / 3.1821×10⁻⁷ | 2.2042×10⁻² |
| fine, early / late | 3.6427×10⁻⁸ / 2.0478×10⁻⁸ | 1.3776×10⁻³ |

### The mode's meridional shape — `the_rossby_packet_keeps_the_gravest_modes_meridional_shape`

The gravest Rossby mode's thermocline anomaly is `2·ψ₀ + ψ₂/2`, i.e. a
`ψ₂`/`ψ₀` ratio of exactly 0.25 — double-lobed *off* the equator, which is what
tells it from a Kelvin wave.

| | coarse | fine |
|---|---|---|
| predicted ratio | 0.25 | 0.25 |
| **measured** | **0.25226** | **0.25092** |
| error | 0.90231% | 0.36627% |
| budget | 22.515% | 6.7678% |

---

## T-07.3 — Equatorial deformation radius

**Suite:** `engine/tests/deformation_radius.rs` · **Issue:** #31 ·
**Claim:** `CONTEXT.md`, *Equatorial deformation radius* — `Le = √(c/β)` is
the meridional scale over which equatorial waves decay away from the equator.

### What is measured

T-07.1 and T-07.2 measured how fast the waves travel; this measures how wide
they are. The zonal sum of one invariant, row by row, is fitted against
`ψₘ(y/L)` with the amplitude eliminated analytically and the trapping scale `L`
left free — a coarse scan over a bracket spanning a factor of nine, then
golden-section refinement. Each wave is fitted on the invariant and structure
Matsuno 1966 / Gill § 11.6 put it on:

| wave | invariant | shape |
|---|---|---|
| Kelvin | eastward `u/c + h/H` | `ψ₀` |
| gravest Rossby | westward `u/c − h/H` | `ψ₀` |
| gravest Rossby | eastward | `ψ₂` |

The profile is read at the **end** of a flight of several packet widths, never
at `t = 0`: at `t = 0` the state *is* the analytic profile and a fit of it would
measure only the fitting code. `Le` is imposed nowhere in the solver — it
emerges from `β` and the pressure gradient — so what these tests assert is that
a run *stays* on the waveguide.

Because a best-fit scale exists for any profile whatever, every point check
also asserts that the fit is a *good* one: `ε = √(1 − ρ)`, the profile's
departure from `ψₘ` at its own best-fitting scale as a relative amplitude, is
held to the **same** budget — which is right, because the budget is a budget on
contamination amplitude in the first place.

### The fits — `the_kelvin_pulse_decays_meridionally_on_the_deformation_radius`, `the_gravest_rossby_packet_decays_meridionally_on_the_same_radius`, `the_rossby_packets_off_equatorial_lobes_sit_on_the_same_radius`

Predicted in every row: `Le = 345065.386842516 m`.

| fit | grid | **fitted scale** | error | shape error `ε` | budget |
|---|---|---|---|---|---|
| Kelvin `ψ₀`, coarse | 100×80 | **345055.87 m** | 0.0027568% | 0.0097687% | 5.0881% |
| Kelvin `ψ₀`, fine | 200×160 | **345063.04 m** | 0.00068154% | 0.0023981% | 1.2720% |
| Rossby `ψ₀`, coarse | 128×160 | **345041.34 m** | 0.0069696% | 0.022533% | 7.1664% |
| Rossby `ψ₀`, fine | 256×320 | **345059.17 m** | 0.0018008% | 0.0052410% | 2.9307% |
| Rossby `ψ₂`, coarse | 128×160 | **345106.61 m** | 0.011945% | 0.034282% | 7.1664% |
| Rossby `ψ₂`, fine | 256×320 | **345076.16 m** | 0.0031231% | 0.0076119% | 2.9307% |

The `ψ₂` row is the sharpest of the three and deliberately not folded into the
`ψ₀` one: `ψ₂` has two off-equatorial lobes and a node at `ŷ = ±1/√2`, so a run
that had trapped the wave on the wrong scale would put the node in the wrong
place, and no amplitude can absorb that.

Two of those six fits are read by more than one test.
`the_rossby_packets_off_equatorial_lobes_sit_on_the_same_radius` makes the `ψ₂`
point check at the coarse resolution only; the fine `ψ₂` row is what the two
convergence tests below read the rate off, and it is quoted here because it is
the same fit of the same flight.

**Where the budgets come from:**

| term | Kelvin coarse | Kelvin fine | Rossby coarse | Rossby fine |
|---|---|---|---|---|
| meridional truncation `(2m+1)·(Δy/Le)²` | 2.0996% | 0.52490% | 2.6245% | 0.65613% |
| zonal truncation `(Δx/σ)²/4` | 0.44444% | 0.11111% | 0.19930% | 0.049825% |
| stray energy `⟨k̂²⟩` | 0 | 0 | 0.75938% | 0.75938% |
| × `TRUNCATION_SAFETY` = 2 | **5.0881%** | **1.2720%** | **7.1664%** | **2.9307%** |

The Kelvin row's stray-energy term is *zero* rather than small: that branch is
an exact solution at every wavenumber. The Rossby budget falls by less than
four under refinement because its stray-energy term does not shrink with the
cells at all — the entry the module header names as the reason the fine budget
is 2.93% and not 1.79%.

Dropped as negligible: wall clipping `ψ₀(±4.1) = 2×10⁻⁴`, fit quadrature
`< 10⁻¹⁰⁰` (Euler–Maclaurin on the whole line, not a finite interval), and RK4
time truncation `10⁻⁸`.

### Both quantities are second order — `the_fitted_decay_scale_converges_at_the_schemes_second_order`, `the_meridional_shape_error_converges_at_the_schemes_second_order`

| fit | fitted-scale order | shape-error order |
|---|---|---|
| Kelvin `ψ₀` | **2.0161** | **2.0263** |
| Rossby `ψ₀` | **1.9524** | **2.1041** |
| Rossby `ψ₂` | **1.9354** | **2.1712** |

All six inside `[1.5, 3]`. The shape error is the sharper of the two: a fitted
scale responds only to the part of a departure that is not orthogonal to it,
whereas `ε` sees the whole of it — and `ε` is precisely the quantity the
budget's leading entry is a bound on, so its rate is what says the budget
*describes* the scheme rather than merely bounding it.

### `Le` is a prediction, not a length — `the_decay_scale_follows_le_across_oceans`

One ocean cannot tell `√(c/β)` apart from any other number that happens to be
345 km. The same pulse, basin, packet and measurement are therefore run in
three oceans whose radii sit a factor of `√2` apart, changing a *different* one
of `Le`'s two parameters each time:

| ocean | `c` | predicted `Le` | **fitted scale** | error | budget | shape error `ε` |
|---|---|---|---|---|---|---|
| equatorial Pacific | 2.7386 m/s | 345 065.39 m | **345 055.87 m** | 0.0027568% | 5.0881% | 0.0097687% |
| `g'` × 4 (`c` doubled) | 5.4772 m/s | 487 996.15 m | **487 781.06 m** | 0.044076% | 2.9885% | 0.030469% |
| `β` × 2 (`c` unchanged) | 2.7386 m/s | 243 998.07 m | **243 155.52 m** | 0.34531% | 9.2873% | 0.028972% |

Each fit is held to its own ocean's budget — the same three terms as everywhere
else in this ticket, evaluated at that ocean's `Le`, which is why the budget
moves between rows although the grid does not. But the *discriminating* claim
is the ordering rather than the budget: each fit must be nearer its own ocean's
`√(c/β)` than either neighbour's, and the neighbours are 29% and 41% away where
the largest error above is 0.35%. A fit that returned a fixed length would fail
in two of the three.

---

## T-07.4 — Steady wind-driven thermocline tilt

**Suite:** `engine/tests/steady_wind_tilt.rs` · **Issue:** #32 ·
**Claim:** `CONTEXT.md`, *Thermocline tilt* — a sustained easterly stress
leaves the thermocline deep in the west and shallow in the east, and that
steady slope is the control case the model must reproduce.

### The prediction is not the textbook balance

The textbook tilt `g'·∂h/∂x = τx/(ρ₀·H)` is **not** a steady solution of the
equations the engine integrates. Steady continuity is `r·h + H·∇·u = 0`, so a
tilted thermocline is damped at `r·h`, mass must be fed to it, and the current
cannot be exactly zero. Keeping that term — a wind stress closed by linear
Rayleigh damping rather than by advection, which is the Sverdrup/Stommel-type
balance the ticket names, read on the equator where `f` vanishes — gives a
closed form:

```text
h(x) = (A/k)·sinh(k·(x − L/2))/cosh(k·L/2),   k = r/c,   A = τx/(ρ₀·H·g')
```

In a basin **one cell tall** this is exact rather than asymptotic: both `v`
rows are coast, so `NoNormalFlow` pins `v ≡ 0` and the rotation terms drop out
of the system rather than being small in it.

Configuration: `L` = 15 000 km, `τ₀` = −0.05 Pa (the control scenario's alizés),
and `δ = r·L/c = 2` — a harder ocean than the shipped scenario's `δ ≈ 0.6`,
chosen deliberately: at `δ = 2` the damped tilt is `tanh(1) = 76%` of the
undamped one, so a solver that had lost the `r·h` term would miss by 24% of the
tilt rather than hide inside a 2.5×10⁻⁵ tolerance.

| | |
|---|---|
| `r` | 3.6514837167011077×10⁻⁷ s⁻¹ |
| undamped tilt `A·L` | −97.561 m |
| **damped tilt `A·L·tanh(δ/2)/(δ/2)`** | **74.302 m** |

### The channel profile — `the_steady_channel_tilt_matches_the_analytic_damped_balance`

60 columns (Δx = 250 km), integrated 40 damping times (1500 steps) from rest.

**Equilibrium is asserted before anything is read.** Damping every prognostic
variable at the same rate makes an unforced departure decay as `e^{−r·t}` in
the discrete energy norm, and `r` is the *only* rate in the system, so there is
no slow creep to miss:

| | |
|---|---|
| drift over the run's second half | 1.8956×10⁻⁹ of the tilt |
| equilibrium bound `(e^{−rT} + e^{−rT/2})·√N` | 1.5966×10⁻⁸ |

Then the profile:

| | |
|---|---|
| **departure from the analytic damped tilt** | **2.5075939×10⁻⁵ of the tilt** |
| tolerance | **2.5091925×10⁻⁵** |
| — zonal truncation (`κ` for `k`, evaluated) | 2.5075939×10⁻⁵ |
| — equilibrium | 1.5966×10⁻⁸ |
| — round-off (`steps·N·ε`) | 1.9984×10⁻¹¹ |

This is the tightest comparison in the suite, and the closest thing here to a
direct hit. The measured departure and the *derived* truncation term differ by
1.4×10⁻¹⁵ of the tilt — round-off — so the run is not merely inside the
continuous closed form's neighbourhood, it is sitting on the **discrete**
closed form the C-grid elimination predicts, to the last bits an `f64` has.
The margin against the tolerance is only 1.0006×, and that is not luck either:
the tolerance is a derivation of where the run should land, not an allowance
for where it might.

There is no `TRUNCATION_SAFETY` factor here because nothing is left
unevaluated. The truncation term is the largest difference between the two
closed forms — the continuous one at `k = r/c`, and the discrete one at
`κ = (2/Δx)·asinh(k·Δx/2)`, which the C-grid elimination admits *exactly* —
sampled on the channel's own columns. No term of the
`κ = k·(1 − (kΔx)²/24 + …)` series is dropped on the way in, which is why this
budget is quoted to four significant figures where the wave suites' are quoted
to two.

Direction, on the same run: the settled thermocline stands at **+36.342 m**
against the western wall and **−36.342 m** against the eastern one (analytic:
±36.343 m). Deep in the west, shallow in the east.

### Second order, with the ratio predicted rather than assumed — `the_channel_tilt_error_falls_at_second_order_with_resolution`

The predicted ratio here is not "a quarter": the exact `asinh` carries
`(k·Δx)⁴` terms too, so the prediction is *evaluated* from the two closed forms
at each cell width rather than assumed.

| | |
|---|---|
| coarse departure (60 columns) | 2.5075939×10⁻⁵ |
| fine departure (120 columns) | 6.3256554×10⁻⁶ |
| **measured ratio** | **0.2522599593218** |
| **predicted ratio** | **0.2522599593312** |
| discrepancy | 9.4×10⁻¹² |
| tolerance (the equilibrium floor under each error) | 4.2061×10⁻³ |

Measured and predicted agree to ten significant figures — eight orders of
magnitude inside the tolerance the test allows. The ratio exceeds 1/4 by 0.9%,
which is the `(k·Δx)⁴` correction showing up exactly where the derivation says
it should.

### The equilibrium guard bites — `a_run_too_short_to_equilibrate_is_reported_as_not_equilibrated`

A steady-state test that assumed equilibrium instead of checking it would have
an answer — a wrong one — on a transient. The same channel run for **2** damping
times instead of 40 (76 steps):

| | measured | what an equilibrated run is allowed |
|---|---|---|
| drift over the second half | **0.26447** | 1.5966×10⁻⁸ |
| departure from the analytic profile | **0.13225** | 2.5092×10⁻⁵ |

Both halves of the criterion hold: the guard fires (by seven orders of
magnitude), *and* it fires on something that matters — such a run really is
13.2% away from the steady tilt, which is `e^{−2} = 13.5%` of the initial
departure still in place, as the energy bound says it should be. Without the
guard, a suite could pass a transient.

### The rotating basin

The channel is where the tilt has a closed form; it is not where the equatorial
Pacific is. The second configuration is the basin the shipped control scenario
describes — trades decaying on `Le`, beta-plane, closed walls — at 80 × 64 cells
(Δx = 172.53 km, Δy = 43.13 km), `δ = 2`, 25 damping times (3646 steps).
Equilibrium first: drift **1.8180×10⁻⁶** against a bound of 2.6666×10⁻⁴.

**No net anomaly** — `the_equilibrated_basin_carries_no_net_thermocline_anomaly`.
Summing the discrete continuity equation over the basin telescopes the
divergence to the flux through four walls that carry none, leaving `r·Σh = 0`
**exactly**, whatever shape the wind has.

| | |
|---|---|
| mean `h` | 5.898×10⁻¹⁷ m |
| tilt scale | 20.770 m |
| **ratio** | **2.8397×10⁻¹⁸** |
| tolerance (`e^{−rT}` + round-off) | 4.1589×10⁻⁹ |

**The Kelvin invariant's steady balance** —
`the_kelvin_invariant_obeys_the_steady_damped_balance`. Projecting the summed
zonal-momentum and continuity equations on `ψ₀` annihilates the meridional term
by one integration by parts, so `(r/c)·q₀ + dq₀/dx = X₀/c²` holds at every `x`
with no long-wave approximation and no truncation of the Rossby set. For
`τx = τ₀·exp(−(y/Ly)²)` the forcing `X₀ = τ₀/(ρ₀·H·√((Le/Ly)² + 1/2))` is
analytic, so this is a prediction and not a tautology. Asserted in integrated
form over the basin's middle half, clear of both walls' reflection layers:

| | |
|---|---|
| analytic `X₀·span/c²` | −0.24433145457206942 |
| **measured `[q₀] + (r/c)·∫q₀ dx`** | **−0.24430770367260538** |
| **residual** | **9.7208×10⁻⁵ of the forcing** |
| tolerance | **1.6089×10⁻²** |
| — meridional quadrature `(Δy/Le)²` | 1.5625×10⁻² |
| — waveguide tail (`ψ₀` outside `±4 Le`) | 6.6915×10⁻⁵ |
| — zonal quadrature `(k·Δx)²·(1/12 + 1/8)` | 1.3021×10⁻⁴ |
| — equilibrium | 2.6666×10⁻⁴ |
| — round-off | 4.1450×10⁻⁹ |

**The tilt itself** —
`the_equatorial_thermocline_deepens_to_the_west_across_the_whole_basin`. The
`ψ₀` projection of `h` stands at **+18.314 m** against the western wall and
**−28.944 m** against the eastern one, and falls **monotonically** across every
column pair of the interior half. Deep in the west, shallow in the east, with
no interior reversal. This one carries no tolerance and needs none: what is
asserted is a sign and an ordering, not a magnitude, and the magnitude is what
the channel comparison above already pinned to 2.5×10⁻⁵.

---

## T-07.5 — Conservation in the undamped, unforced limit

**Suite:** `engine/tests/conservation.rs` · **Issue:** #33 ·
**Claim:** with `r = 0` and `τ = 0` the discrete wave energy
`E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)` is conserved up to the scheme's own
truncation.

### Configuration

An equatorial channel 10 000 km × 1000 km started in the gravest zonal standing
mode (`h = 20·cos(π x/Lx)`, at rest), integrated for **32 basin crossings** at
`c` — about **3.7 simulated years**, 4548–18 108 steps, and some 148 periods of
the equatorial inertia-gravity motion (`2π/√(βc)` = 7.92×10⁵ s). Eight times the
length of the T-02.5 check this ticket formalises. Three meridional
resolutions, `Δx` fixed: the pressure-gradient/continuity pair is *exactly*
skew at every `Δx`, so the whole energy error is meridional.

### The bound

Two terms, both derived from the scheme and neither measured:

1. **The C-grid Coriolis pair's `O(Δy²)` skewness defect.** `u`'s equation
   evaluates `f = β·y` on cell-centre rows and `v`'s on face rows, half a cell
   apart, leaving `dE/dt = −H·(β·Δy²/4)·Σ_u u·∂v/∂y`. With `|∂v/∂y| ≤ |v|/Le`
   and `E ≥ H·Σ|u||v|` that is a *rate* `(1/4)·√(βc)·(Δy/Le)²`. The defect
   **oscillates rather than accumulates** — `u` and `∂v/∂y` are two components
   of one wave field — and the integral of `a·cos(ωt)` has amplitude `a/ω`, so
   at `ω = √(βc)` the excursion is `(1/4)·(Δy/Le)²`, independent of run length.
   This is what the long run exists to expose: a defect that accumulated would
   overshoot by the ~148 periods in the run, two orders of magnitude.
2. **RK4's numerical dissipation**, `N·θ⁶/72` from
   `|R(iθ)|² = 1 − θ⁶/72 + θ⁸/576`, with `θ = √(βc)·dt`.

Term 2's *size* is an estimate and is labelled as one in the source: `√(βc)` is
the frequency of the trapped motion the run is about, not a supremum over the
discrete spectrum. Term 2's *sign* is rigorous — `|R(iθ)| ≤ 1` inside the CFL
bound — which is why the gain test below leans on term 1 alone.

### The drift — `energy_drift_over_a_long_undamped_unforced_run_stays_within_the_derived_bound`

| rows | `Δy` | `dt` | steps | **worst drift** | bound | skewness | RK4 | margin |
|---|---|---|---|---|---|---|---|---|
| 16 | 62 500 m | 25 691.75 s | 4 548 | **1.1527×10⁻³** | 1.2741×10⁻² | 8.2016×10⁻³ | 4.5397×10⁻³ | 11.1× |
| 32 | 31 250 m | 12 893.84 s | 9 062 | **2.7895×10⁻⁴** | 2.1949×10⁻³ | 2.0504×10⁻³ | 1.4453×10⁻⁴ | 7.9× |
| 64 | 15 625 m | 6 452.96 s | 18 108 | **6.8893×10⁻⁵** | 5.1714×10⁻⁴ | 5.1260×10⁻⁴ | 4.5380×10⁻⁶ | 7.5× |

### Energy is never gained past the rigorous bound — `an_undamped_unforced_run_never_gains_energy_past_the_skewness_bound`

RK4 can only *remove* energy from this system, so every joule the basin gains
must come from the skewness defect — and this assertion uses term 1 alone, with
no frequency estimate anywhere in it.

| rows | **worst gain** | skewness bound | margin |
|---|---|---|---|
| 16 | **4.5566×10⁻⁴** | 8.2016×10⁻³ | 18.0× |
| 32 | **1.4892×10⁻⁴** | 2.0504×10⁻³ | 13.8× |
| 64 | **3.8612×10⁻⁵** | 5.1260×10⁻⁴ | 13.3× |

### The drift is second order — `the_long_run_energy_drift_falls_at_the_schemes_second_order_under_refinement`

| refinement | **measured order** |
|---|---|
| 16 → 32 rows | **2.0469** |
| 32 → 64 rows | **2.0176** |

Held to `≥ 1.8` (second order less the 0.2 that admits a sub-leading `O(Δy⁴)`
term worth up to ~18% of the leading one). No bound can fake this: it is the
*measured* drift falling at the rate term 1 is written about, and term 2 falls
faster still (`dt⁵` at fixed run length).

### The run is not vacuous — `the_conservation_run_exchanges_energy_between_potential_and_kinetic_form`

A run in which nothing moved would conserve energy perfectly and prove nothing.
The initial state is pure potential energy at rest; the peak kinetic fraction
reached was **93.41%**, **93.36%** and **93.31%** at 16, 32 and 64 rows, against
a floor of 50%.

### A note on volume

`01-scientific-model.md` names "energy/volume conservation" as one target and
this suite asserts only the energy half: `Σh` is a *linear* invariant of the
same equations, so a volume test would be a much weaker statement made with the
same machinery, and both review axes on T-07.5 called it scope creep. It is not
unmeasured, though — it is measured where it is a non-trivial claim, in T-07.4's
`the_equilibrated_basin_carries_no_net_thermocline_anomaly`, which holds a
*forced, damped, rotating* basin's `Σh` to
[2.8397×10⁻¹⁸](#the-rotating-basin) of the tilt.

---

## What this suite does not validate

Stating the boundary is part of the claim.

- **The v1 core is linear** (CODING_STANDARDS.md § *Scope guards*). Every
  prediction above is a property of the linearised equations. Nonlinear
  advection is out of scope, so nothing here says what the model would do at
  amplitudes where it matters.
- **The boundaries are validated elsewhere.** Every wave test above runs in
  open water, four packet widths clear of both zonal walls, deliberately. The
  reflection physics is T-04.3 (`western_boundary_reflection.rs`) and T-04.4
  (`eastern_boundary_reflection.rs`); the no-normal-flow condition itself is
  T-04.2.
- **No comparison against observations.** These are checks against *analytic
  solutions of the model's own equations*, not against the Pacific. They
  establish that the code solves the equations it claims to solve; whether
  those equations describe the ocean is the 1.5-layer model's own literature
  (Gill ch. 11; Cane & Sarachik 1981; Matsuno 1966).
- **One-directional error control.** The convergence tests establish an order,
  not an absolute error at any given production resolution. A scenario run at a
  coarser grid than these tests use inherits the same second-order behaviour
  with a correspondingly larger coefficient.
- **T-07.4's rotating-basin statements are modal, not a closed form.** The
  exact tilt profile is available only in the channel; in the rotating basin
  what is asserted is `Σh = 0`, the `ψ₀` balance, and the direction and
  monotonicity of the slope.

## Reproducing these figures

```sh
cargo test --workspace --all-features   # the gate, including these five suites
cargo test -p engine --test kelvin_wave_propagation \
                     --test rossby_wave_dispersion \
                     --test deformation_radius \
                     --test steady_wind_tilt \
                     --test conservation
```

A passing run prints no numbers — every figure above lives inside the assertion
that consumes it, and Rust reports an assertion's message only when it fails.
Each suite's module header derives its own budget from first principles and
names the accessor each quantity comes from, so a reader who wants a figure
regenerated can print it from the same accessor without re-deriving anything.
The failure messages are written to carry the same numbers: a broken run
reports what it measured, what theory predicted, and which term of the budget
it exceeded.
