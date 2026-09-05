# The physics, explained

This is the plain-language companion to
[`docs/planning/01-scientific-model.md`](planning/01-scientific-model.md).
That document states the equations the engine solves; this one explains what
they mean, so that when you open a run in the visualizer you know what the
colours and the moving fronts actually are.

No derivations here, and only a handful of symbols — the ones printed on the
visualizer's axes and used in scenario files. The glossary in
[`CONTEXT.md`](../CONTEXT.md) defines each of them precisely.

## 1. The thermocline: a lid on the cold ocean

The tropical Pacific is not a well-mixed bathtub. Sunlight and wind stir the
top hundred-and-something metres into a warm, buoyant layer; below that the
ocean is cold, dense and — on the timescales we care about — nearly still. The
two do not blend gradually. They are separated by a thin boundary across which
the temperature falls sharply: the **thermocline**.

Because the warm layer floats on the cold one, the thermocline behaves less
like a temperature contour and more like the *surface of an upside-down ocean*.
Push warm water into one region and the boundary bulges downward there; pull
warm water away and it rises. Its shape is the memory of everything the wind
has done recently, which is why it, and not sea-surface temperature, is the
primary quantity this simulation tracks.

The model is deliberately built around that picture. It represents the ocean as
**one active layer** of warm water of mean thickness `H`, sitting on an
infinitely deep motionless abyss — the "1.5-layer reduced-gravity" model, a
standard tool in tropical ocean dynamics. There is a single active layer, plus
half a layer's worth of information about the passive deep ocean beneath it:
hence 1.5.

Two consequences of that choice are worth carrying with you:

- The simulation has **no vertical structure to plot**. There is no temperature
  profile with depth. There is a warm layer, a cold abyss, and the surface
  between them.
- The state variable is `h`, the thermocline depth **anomaly** in metres — the
  departure from the mean depth `H`, not the depth itself. Total thermocline
  depth is `H + h`. Positive `h` means a deeper-than-average thermocline, so
  a thicker slab of warm water; negative `h` means cold water sits closer to
  the surface. When the visualizer shows `h = 0` everywhere, that is a flat
  resting ocean, not an ocean with no thermocline.

### Why the layer moves so slowly

A disturbance of the *sea surface* travels at √(gravity × depth) — for a
four-kilometre-deep ocean, some 200 m/s. The thermocline is different.
Warm water floating on cold water is a much weaker restoring system than water
floating on air, because the density contrast across the thermocline is a few
parts per thousand rather than a factor of eight hundred.

That weakening is captured by the **reduced gravity** `g'`: instead of the full
`g ≈ 9.8 m/s²`, the layer feels an effective gravity of only a few centimetres
per second squared. The wave speed that follows,

```
c = √(g'·H)
```

comes out around 2–3 m/s — walking pace. This is the **Kelvin wave speed**, the
fastest signal in the model, and it is the reason the tropical Pacific has a
memory measured in months: a disturbance takes roughly two to three months to
cross the basin. Everything about ENSO's timescale traces back to this number.

## 2. Why the trade winds tilt it

Over the equatorial Pacific the wind blows persistently from east to west: the
**trade winds**, the *alizés*. In the model these enter as a prescribed surface
stress `τx` — a force per unit area dragged along the top of the warm layer.
Easterly winds mean `τx < 0`, since they push water in the negative-`x`
(westward) direction.

Blow steadily westward on the warm layer and it piles up in the west. It cannot
pile up forever: as the layer thickens in the west and thins in the east, the
resulting pressure difference pushes back eastward. The system settles where
the wind's westward pull is balanced by the pressure gradient of the tilted
layer it created. That equilibrium is the **thermocline tilt**:

```
   west (Indonesia)                                  east (Peru)

        <---  <---  <---  trade winds  <---  <---  <---
   ~~~~~~~~~~~~~~~~~~~~ sea surface ~~~~~~~~~~~~~~~~~~~~~~~
   \
    \___              warm layer (thick in the west)
        \____
             \_____                                          thermocline
                   \______
                          \________
                                   \______  warm layer (thin in the east)
                                          \______
   cold abyss                                     cold water near the surface
```

In the real Pacific this amounts to a thermocline a few tens of metres deep off
South America and well over a hundred metres deep near Indonesia, with the sea
surface itself standing tens of centimetres higher in the west. It is the
observed mean state of the basin, and reproducing it under steady easterlies
is one of the analytic results the engine is validated against.

The tilt is also why the eastern Pacific is cold. With the thermocline near the
surface there, the same winds that tilt it also drag surface water away from
the coast and the equator, and cold water is drawn up from just below to
replace it. That upwelling has only a short distance to travel. Push the
thermocline down in the east and the upwelling brings up water that is no
longer cold — which is exactly what happens during El Niño, and the hinge on
which section 4 turns.

### The equator is a special place

Away from the equator, the Earth's rotation deflects moving water — right in
the northern hemisphere, left in the southern. At the equator that deflection
vanishes and then reverses sign as you cross. The model captures this with the
**beta-plane** approximation, `f = β·y`: the deflection is exactly zero on the
equator and grows linearly with distance north or south of it.

The consequence is that the equator acts as a waveguide. Water that strays
north gets bent back south, water that strays south gets bent back north, and
disturbances end up **trapped in a band centred on the equator**. The width of
that band is the **equatorial deformation radius**, `Le = √(c/β)` — a few
hundred kilometres, a few degrees of latitude. This is why the visualizer's
signals hug the equator and fade towards the top and bottom of the basin
instead of spreading everywhere.

## 3. The two waves you will see moving

Disturb the tilt — with a wind burst, or by switching the trade winds on at the
start of a run — and the basin does not simply slide to a new equilibrium. It
gets there by sending waves. The waveguide admits several kinds, but two of
them carry the story — and they are the two the engine is validated against.

**Kelvin waves** travel **eastward only**, at the full speed `c ≈ 2–3 m/s`,
right along the equator. They are non-dispersive: a pulse keeps its shape as it
travels, so a bump leaving the western Pacific arrives in the east two or three
months later still recognisably a bump. In the visualizer they are the crisp,
fast-moving front that crosses the basin left to right.

**Rossby waves** travel **westward only**, and the fastest (gravest) one moves
at just `c/3` — call it eight or nine months to cross the basin. They are
dispersive, so a pulse spreads and smears as it goes, and they are broader in
latitude than Kelvin waves, appearing as lobes straddling the equator rather
than a single sharp equatorial front.

That asymmetry — fast east, slow west, one direction each — is not a modelling
convenience. It falls out of the rotating shallow-water equations, and it is
the clock that sets ENSO's rhythm.

The waves do not disappear at the edges of the basin, because the basin is
closed on all four sides:

- A **Kelvin wave reaching the eastern boundary** (South America) cannot keep
  going. Its energy turns around and heads back west as Rossby waves. (In the
  real Pacific some of it also escapes poleward along the coast — the model's
  closed boundary has nowhere to send it.)
- A **Rossby wave reaching the western boundary** (the Indonesian maritime
  continent) reflects, and *part* of its energy comes back as an eastward
  Kelvin wave. Only part: the rest goes into other, slower modes. The
  partial-ness matters, and Epic 04 tests this reflection on its own.

So a single disturbance is not a one-way trip. It becomes a circuit: east as a
Kelvin wave, back west as Rossby waves, east again as a reflected Kelvin wave,
each lap taking longer than the last leg's speed alone suggests, each lap
weaker than the one before because damping (`r`, the **Rayleigh damping**
coefficient) steadily bleeds energy away. Watch a run long enough and you are
watching that circuit — the delayed feedback that keeps the basin oscillating
rather than simply relaxing.

## 4. How this produces ENSO

Now put the pieces together. The mean state — trade winds, tilted thermocline,
cold upwelling in the east — is a balance in which each part sustains the
others. That mutual sustaining is the **Bjerknes feedback**:

```
   strong trade winds
        │
        ▼
   steeper thermocline tilt: deep west, shallow east
        │
        ▼
   cold upwelling reaches the surface in the east
        │
        ▼
   big east–west temperature contrast across the basin
        │
        ▼
   strong trade winds   ← (the loop closes)
```

Every arrow reinforces the next. A loop like that amplifies whatever nudges it,
in either direction.

**El Niño** is the loop running backwards. Something weakens the trade winds —
in the model, an idealized **westerly wind burst**, a patch of `τx > 0`
superimposed on the trades. Warm water that the winds had been holding in the
west sloshes east as a Kelvin wave. It arrives two or three months later,
pushes the eastern thermocline down, and cuts the upwelling off from cold
water. The east warms, the east–west contrast shrinks, and the trade winds
weaken further — which sends more warm water east. The tilt collapses; warm
water surfaces in the east. That is the warm phase.

**La Niña** is the same loop turned up: stronger trades, steeper tilt, colder
east, stronger trades still. The cold phase.

**Why it does not simply stay stuck** is the part the waves explain. A wind
burst does not only launch that eastward Kelvin wave. It launches westward
Rossby waves at the same time, and those carry the *opposite* sign: where the
Kelvin wave pushes the thermocline down, the Rossby waves lift it. They creep
west at `c/3`, reflect off the western boundary into a Kelvin wave that keeps
that opposite sign, and only then head east — arriving many months after the
event that launched them, long after the feedback has done its work, to push
the eastern thermocline back up. That delayed return is what shuts an El Niño
down, and can overshoot into a La Niña. A positive feedback plus a delay is an
oscillator, and the delay here is a basin width divided by a wave speed. It is
no accident that ENSO's period is a few years.

## 5. What the engine actually simulates — and what it does not

Read that story back and you will notice it has two halves: the ocean
responding to wind, and the wind responding to the ocean. Version 1 of this
engine implements the first half only.

**In the model:**

- The three prognostic fields `h`, `u` and `v` — thermocline depth anomaly and
  the zonal and meridional current anomalies of the upper layer.
- Wind stress as a **prescribed** function of position and time. Three
  scenarios: steady trade winds (the control case that should reproduce the
  observed tilt), a seasonal cycle modulating them, and a westerly wind burst
  superimposed on top.
- The full equatorial wave dynamics — Kelvin and Rossby propagation, equatorial
  trapping, and reflection off closed western, eastern, northern and southern
  boundaries.
- Rayleigh damping `r`, a simple stand-in for all the mixing and dissipation
  the model does not resolve.

**Not in the model (v1):**

- **Sea-surface temperature**, unless you ask for it. There is no `T'` field in
  the linear core, so by default the engine simulates the ocean *dynamics* of
  ENSO — where the warm water is — not the temperature anomalies that ENSO is
  defined by. Since the eastern thermocline depth and eastern SST move
  together, `h` in the east is a good proxy, and that is how to read a default
  run. A scenario that carries an `[sst]` section gets a real mixed-layer `T'`
  integrated alongside, warmed and cooled by the water the trade-driven
  upwelling draws up from just below the thermocline (T-12.1); the linear core
  it rides on is unchanged, down to the last bit. Such a run writes `T'` into
  its frames as a sixth variable, in kelvin, and `termocline inspect` lists it;
  a run without the section records the anomaly as *absent* rather than as a
  basin of zeros, so a plot can never show a temperature the model never
  computed.
- **The atmospheric half of the Bjerknes loop.** The winds do not respond to
  the ocean. Feed the engine a wind burst and you get the ocean's response to
  that burst; you do not get a self-sustaining oscillation, because the arrow
  from "temperature contrast" back to "trade winds" is the one thing not
  wired up. Closing that loop — adding an SST equation and a simple atmospheric
  response — is the Epic 12 extension, and it is what turns a prescribed-wind
  response into an *emergent* ENSO.
- **Nonlinear advection.** The equations are linear: the currents transport
  nothing, including themselves. This is a good approximation for the wave
  dynamics and a poor one for strong events, and it is a deliberate scope
  choice — a linear core that can be checked against analytic solutions first.

## 6. Reading the visualizer

With all that in hand, a run should decode as follows.

- **The `h` map** is a top-down view of the basin, coloured by thermocline
  depth anomaly in metres. One end of the colour scale is a shallower-than-mean
  thermocline (negative `h`, cold water closer to the surface); the other is a
  deeper one (positive `h`, a thicker warm layer). The midpoint is the mean
  depth `H` — a flat resting ocean, not an ocean without a thermocline.
- **The tilt from section 2, seen edge-on** along the equator. Under steady
  trades it should settle deep in the west and shallow in the east and then
  stop changing. That steady tilt is the model's baseline, and matching it
  against the analytic balance is one of the engine's validation targets.
- **A front crossing the basin west to east in two to three months** is a
  Kelvin wave. **Broader lobes drifting east to west, three times slower and
  smearing as they go** are Rossby waves. Both stay within a few degrees of the
  equator — that is the deformation radius `Le` at work.
- **Fronts arriving at a boundary and something departing in the opposite
  direction** is reflection, not a numerical artefact. It is the mechanism
  behind the delay in section 4.
- **Everything fading over a long run** with no wind change is Rayleigh
  damping. Set `r = 0` with no wind and nothing should fade at all: energy is
  conserved in that limit, and that conservation is one of the engine's
  Epic 07 validation targets.

## Where to go next

- [`docs/planning/01-scientific-model.md`](planning/01-scientific-model.md) —
  the same physics as equations, with the numerical scheme and the list of
  analytic results the engine is validated against.
- [`CONTEXT.md`](../CONTEXT.md) — precise definitions of every term used here.
- [ADR-0003](planning/adr/0003-numerical-scheme.md) — how the equations are
  discretized: the Arakawa C-grid, the time-stepping scheme, and the CFL
  condition that ties the timestep to `c`.
