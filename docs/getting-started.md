# Running your first simulation

From a clean checkout to a simulation you can look at: build the two binaries,
run a scenario, inspect the run it wrote, open it in the visualizer, then write
a scenario of your own. Every command below was run from the repository root,
and the outputs quoted are what it printed.

Nothing here explains the physics — [`docs/the-physics-explained.md`](the-physics-explained.md)
does that, and is worth reading once you have a picture on screen to attach it
to.

## 1. Build

You need a Rust toolchain at **1.90 or newer** (`rustup` is the easy way to
get one) and a C linker, which every platform's usual build tools provide.

```sh
git clone https://github.com/hbock-42/pacific-termocline
cd pacific-termocline
cargo build --release --workspace
```

`--release` matters: the engine is a numerical solver, and an unoptimized
build integrates close to two orders of magnitude slower — the twelve-second
run below becomes a twenty-minute one.

That produces the two binaries this guide uses, under `target/release/`:

- `termocline` — the engine. Runs a scenario, writes a run.
- `termocline-viz` — the visualizer. Draws a run.

## 2. Run a scenario

A **scenario** is the engine's unit of input: one TOML file holding the basin,
the physical constants, the wind forcing and how long to integrate for. Three
worked ones ship in [`engine/scenarios/`](../engine/scenarios/). Start with the
control scenario, two years of steady easterly alizés over the equatorial
Pacific:

```sh
./target/release/termocline run \
  --config engine/scenarios/steady-trades.toml \
  --out /tmp/run-demo
```

It reports progress on stderr as it goes and finishes with a one-line summary
on stdout:

```
level=info event=run_started scenario=steady-trades total_steps=17520 dt_s=3600 frames=731
 10% | step 1752/17520 | model 73.0 d | elapsed 1.2 s | eta 10.9 s | 1447 steps/s
 ...
100% | step 17520/17520 | model 730.0 d | elapsed 12.0 s | 1464 steps/s
level=info event=run_finished steps=17520 frames=731 elapsed_s=11.966
engine/scenarios/steady-trades.toml: 17520 steps, 731 frames written to /tmp/run-demo
```

Twelve seconds on the laptop this was written on; the progress line reports
the rate yours actually achieves. Add `--quiet` to suppress the progress and
keep only that last line, or `--verbose` for per-frame detail.

**Check your disk first.** A run directory holds two files — a small JSON
header and every saved frame as raw `f64` — so its size is (frames) × (cells)
× (5 variables) × 8 bytes. This one is **941 MB**: 731 frames of a 320 × 100
grid. The 90-day scenario of section 5 is 28 MB, and that section explains
which knobs move the number.

## 3. Inspect what it wrote

`inspect` reads the run's header and nothing else, so it answers "what is in
this directory?" instantly whatever the run's length:

```sh
./target/release/termocline inspect --run /tmp/run-demo
```

```
run: /tmp/run-demo
format version: 1
scenario: steady-trades
grid: 320 x 100 cells
basin extent: 120.0 to -80.0 degrees east, -25.0 to 25.0 degrees north
mean thermocline depth H = 150.0 m
reduced gravity g' = 0.06 m s^-2
beta = 2.3e-11 m^-1 s^-1
Rayleigh damping r = 1e-7 s^-1
reference density rho_0 = 1025.0 kg m^-3
frames: 731, one every 86400.0 s
variables: h [m], u [m s^-1], v [m s^-1], tau_x [N m^-2], tau_y [N m^-2]
```

Every number is in the SI unit the header records — the command never
rescales a value on the way out, so what you read is what the run was
integrated with. `frames: 731, one every 86400.0 s` is a frame a day across
730 days.

## 4. Look at it

Natively, name the run directory on the command line:

```sh
./target/release/termocline-viz /tmp/run-demo
```

A window opens on the basin: thermocline depth anomaly `h` as a colour map,
red where the thermocline is deeper than its mean and blue where it is
shallower, with the wind stress that forced it drawn over the top as arrows.
The scrubber above the map chooses the frame — drag it, step with the arrow
keys, a page of ten with Page Up and Page Down, Home and End for either end of
the run.

You can also start it with no argument and drag the run's `header.json` and
`frames.bin` onto the window.

Under the steady trades, the story to look for is the tilt: the arrows point
west along the equator, and over the first months the west end of the basin
goes red while the east goes blue and then stops changing. That is the trades
piling warm water up in the west. §6 of
[the physics explained](the-physics-explained.md) reads the rest of the
picture.

## 5. Write your own scenario

Everything about a run is in its scenario file, and the fastest way to see the
model respond is to change one number. Put this in `quick-look.toml` — the
same steady trades as the control scenario, but 90 days on a one-degree grid
instead of two years on a half-degree one:

```toml
# A short, small first run: 90 days of steady alizés over the Pacific basin
# at one degree, a frame a day.

[basin]
resolution_deg = 1.0

[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 2160
output_every_n_steps = 24

[[wind]]
type = "steady_trade_winds"
equatorial_zonal_stress_pa = -0.05
meridional_decay_scale_m = 361000.0
```

```sh
./target/release/termocline run --config quick-look.toml --out /tmp/run-quick
```

Under a second, and 28 MB on disk: 91 frames of a 160 × 50 grid — small enough
to hand to the browser build in the next section.

Three fields decide that cost, and they are worth understanding before you
reach for a bigger run:

- `resolution_deg` sets the cell count, so it is quadratic in both the run's
  time and its size — and it also tightens the CFL bound on `dt_s`.
- `total_steps` × `dt_s` is the model time integrated.
- `output_every_n_steps` sets the frames written, and therefore the whole size
  of the run. Saving every step of the control scenario would write 22 GB.

## 6. See it in a browser

The visualizer also runs as a web app, and there it does not open a run at all:
it *computes* one. A browser has no filesystem and 941 MB is not a download, so
per [ADR-0012](planning/adr/0012-the-browser-runs-the-engine.md) the page links
the engine, holds a scenario and steps it in the tab.
[trunk](https://trunkrs.dev) builds and serves it:

```sh
cargo install trunk
rustup target add wasm32-unknown-unknown
cd visualizer && trunk serve
```

The first build takes a couple of minutes; when it prints `server listening
at: http://127.0.0.1:8080/`, open <http://localhost:8080>. The control
scenario starts computing on load and the map fills in as it goes — a progress
bar counts the frames, and the frame chooser follows the newest one until you
scrub back into the run. Nothing is downloaded and nothing is dropped on the
page: the run on screen was produced by the same engine `termocline run` is.

The scenarios the page offers are in `visualizer/scenarios/`, and they are the
engine's own coarsened to fit a tab: 80 × 25 cells at 2°, a frame every three
days, 244 frames — 19.9 MB of frames against the control run's 941 MB. That
limit is enforced rather than hoped for: a scenario whose frames would not fit
is refused before the first step, with the size it would have needed. Two
degrees still resolves the equatorial waveguide (`Le` ≈ 361 km), but a
*validated* run is a native run of `engine/scenarios/` — nothing scientific
rests on the browser.

Tick **Compare two runs** and each panel gets its own scenario picker, so the
trades and the trades-plus-a-westerly-burst compute side by side on one frame
index and one colour scale. That is the same comparison the native build gives
two run directories: `termocline-viz /tmp/run-quick /tmp/run-burst`.

Natively nothing about loading a written run has changed: a directory on the
command line, the **Open run directory…** button, a `?run=`-style URL in the
Run URL box, or the run's two files dropped on the window.

From here, the two other shipped scenarios are the interesting ones —
`seasonal-cycle.toml` breathes the trades with the year, and `wind-burst.toml`
fires a westerly wind burst that launches a Kelvin wave you can watch cross
the basin. Every field of the format, its unit and its valid range, is in
[`docs/scenario-config-reference.md`](scenario-config-reference.md).

## When something goes wrong

The engine treats bad input as an error to report, never a value to quietly
substitute: it exits non-zero and names what you asked for and the bound it
broke.

- **`dt is 36000 s, past the CFL-stable maximum of 29652.021 s for this grid
  spacing and wave speed; the run would go unstable. Set dt to at most
  29652.021 s, or coarsen the grid`** — the timestep and the grid are tied
  together by the wave speed `c = √(g'H)`, so a finer `resolution_deg` needs a
  smaller `dt_s`. The message names the bound it wants.
- **`this is not a scenario: TOML parse error ... unknown field
  'resolution_degrees', expected one of ...`** — a key this build does not
  define. Every section rejects unknown fields rather than ignoring them, so a
  misspelling is an error naming the key and listing the alternatives, never a
  silently skipped setting.
- **`.../header.json could not be opened`** from `inspect` or the visualizer —
  the path is not a run directory. A run is exactly the two files `run --out`
  wrote, and both must be there.
- **A run that takes minutes instead of seconds** — a debug build. Rebuild
  with `--release`.

## Where to go next

- [`docs/the-physics-explained.md`](the-physics-explained.md) — what the
  picture on screen means, in plain language.
- [`docs/scenario-config-reference.md`](scenario-config-reference.md) — every
  scenario field, its unit and its valid range.
- [`docs/validation-report.md`](validation-report.md) — how we know the
  simulation is right: each scientific test, its analytic prediction, and the
  measured result.
- [`CONTRIBUTING.md`](../CONTRIBUTING.md) — the workflow, and the gate a
  change has to pass.
