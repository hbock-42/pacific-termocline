//! Acceptance tests for T-04.4 — Kelvin→Rossby reflection at the eastern
//! boundary.
//!
//! An equatorial Kelvin wave travels **eastward only** (`CONTEXT.md`), so the
//! eastern wall is the one it strikes. It cannot pass, and it cannot stand
//! still: `docs/planning/01-scientific-model.md` states the outcome — "eastern
//! boundary (~80°W, South America): closed, reflects incident Kelvin energy
//! back as Rossby waves" — and the gravest of those Rossby modes carries the
//! reflected signal **west** at `c/3`. This file initializes an incident
//! Kelvin pulse in the interior, lets the closed basin of T-04.2 do the
//! reflection, and measures both speeds against the analytic theory below.
//!
//! Nothing here is measured out of a run and pasted back in as an expectation:
//! every number asserted is `c`, `c/3`, a ratio derived from the equatorial
//! wave modes, or a convergence rate derived from their dispersion relation,
//! and every tolerance is a sum of named error terms (CODING_STANDARDS.md
//! § Tests).
//!
//! # The modes, and how a run is decomposed into them
//!
//! Write the linear equations of [ADR-0003] in the equatorial variables:
//! `η = y/Le` with `Le = √(c/β)` the deformation radius, velocities scaled by
//! `c = √(g'·H)` and `h` by `H`. The meridional structures are the parabolic
//! cylinder functions `ψ_n(η) = H_n(η)·e^{−η²/2}` with `H_n` the Hermite
//! polynomials — `ψ_0 = e^{−η²/2}`, `ψ_1 = 2η·e^{−η²/2}`,
//! `ψ_2 = (4η² − 2)·e^{−η²/2}` — which satisfy
//!
//! ```text
//! (∂_η + η)·ψ_n = 2n·ψ_{n−1},        (∂_η − η)·ψ_n = −ψ_{n+1}
//! ```
//!
//! and are mutually orthogonal on the line. In the two combinations
//! `r = u/c + h/H` and `q = u/c − h/H` the system splits into an eastward and
//! a westward part,
//!
//! ```text
//! ∂r/∂t + c·∂r/∂x = (η·v − ∂v/∂η)/Le,   ∂q/∂t − c·∂q/∂x = (η·v + ∂v/∂η)/Le,
//! ```
//!
//! — `r` and `q` dimensionless, `v` in m/s, `∂/∂η` the derivative in the
//! stretched latitude — and the modes are what makes each side of that a
//! closed statement:
//!
//! - **Kelvin.** `v ≡ 0`, `u/c = h/H ∝ ψ_0`, so `q ≡ 0` and `r = 2·h/H`
//!   translates east at exactly `c`, without dispersion and without changing
//!   shape. This is an exact solution for *any* zonal profile, which is what
//!   makes a Gaussian pulse a legitimate initial condition rather than an
//!   approximation to one.
//! - **Long Rossby, mode `n ≥ 1`.** `v ∝ ψ_n`, and dropping `∂v/∂t` (the
//!   long-wave limit — the pulse here is many `Le` wide) leaves
//!
//!   ```text
//!   h/H = R·(ψ_{n−1} + ψ_{n+1}/(2(n+1))),   u/c = −R·(ψ_{n−1} − ψ_{n+1}/(2(n+1)))
//!   ```
//!
//!   with `R(x, t)` propagating **west** at `c/(2n+1)`, so `c/3` for the
//!   gravest mode `n = 1` (`CONTEXT.md`, *Rossby wave*). Its `q` carries one
//!   single structure, `ψ_{n−1}`, and its `r` carries `ψ_{n+1}`.
//!
//! Two projections therefore separate the two waves *exactly*, and this is how
//! the tests below measure anything at all:
//!
//! ```text
//! eastward(x) = P₀[r]  — only the Kelvin wave has ψ_0 in r  (mode n puts ψ_{n+1} there)
//! westward(x) = P₀[q]  — only the n = 1 Rossby wave has ψ_0 in q (mode n puts ψ_{n−1} there)
//! ```
//!
//! where `P₀[·]` is the `ψ_0` coefficient of the meridional profile at a
//! column. The westward one is blind to the Kelvin wave (`q ≡ 0`), to the
//! `n = 3, 5, …` Rossby waves the same reflection also launches (their `q` is
//! `ψ_2`, `ψ_4`), and — four deformation radii from the coast, where it is
//! evaluated — to the poleward-running coastal Kelvin wave trapped against the
//! wall. It is the reflected gravest mode and nothing else.
//!
//! # What the eastern wall does, analytically
//!
//! No normal flow means `u = 0` on the wall for every `y`. With the incident
//! Kelvin wave `u/c = A·ψ_0` and the reflected long Rossby modes above, the
//! `ψ_0`, `ψ_2`, `ψ_4`, … coefficients of the total `u` there give
//!
//! ```text
//! A − R₁ = 0,   R₁/4 − R₃ = 0,   R₃/8 − R₅ = 0,   …
//! ```
//!
//! so `R₁ = A`, `R₃ = A/4`, `R₅ = A/32`: the gravest mode takes the incident
//! amplitude, and the rest fall away fast. (The residual meridional mass flux
//! that this pointwise cancellation leaves is carried by Moore's coastal
//! Kelvin wave, trapped within `Le` of the wall and running poleward at `c` —
//! the "boundary-bound" half of the reflection, and the reason the
//! measurements below hold a four-`Le` clearance from the coast.)
//!
//! The coastal wave matters to the run beyond that clearance, and this is
//! where a *closed* basin differs from a textbook half-plane. Its wave guide
//! is the whole perimeter: from the equator up the eastern wall, west along
//! the northern one, down the western wall and back to the equator is
//! 30 000 km, which at `c` takes 127 days. The boundary-bound half of the
//! reflection therefore re-enters the equatorial wave guide *in the west*
//! about 140 days into the wider run, and leaves as an eastward Kelvin wave.
//! It is real physics rather than a numerical artefact, and it has two
//! consequences here: the eastward mode's emptiness is asserted at the first
//! post-reflection sample, before it arrives, and its own reflection off the
//! western wall is an entry in the Rossby speed's error budget. It does not
//! otherwise touch that measurement — a Kelvin wave has `q ≡ 0`, so it is
//! invisible to the projection the reflected packet is tracked in.
//!
//! Two further consequences are asserted directly:
//!
//! - the reflected `h` has the gravest mode's meridional shape,
//!   `ψ_0 + ψ_2/4` = `(η² + 1/2)·e^{−η²/2}` — deepest off the equator, at
//!   `η = ±√1.5`, and only just over half as deep on it;
//! - reflection *compresses the packet by three*. The wall replays the
//!   incident time signature into a wave that leaves at `c/3`, so an incident
//!   Gaussian of zonal width `σ_K` reflects into one of width
//!   `σ_R = σ_K/3`. That is what sets the dispersive error term below, and it
//!   is why the incident pulse here is 11 `Le` wide: the reflected one has to
//!   still be long compared with `Le` for `c/3` to be the speed it travels at.
//!
//! # Where the tolerances come from
//!
//! Both speeds are measured as the displacement of the energy-weighted zonal
//! centroid of the relevant projection between two sample times. For a linear
//! wave that centroid moves at the energy-weighted mean group velocity, which
//! is the whole content of the two error budgets.
//!
//! ## The Kelvin speed: `c`, and only numerical error
//!
//! The Kelvin wave is non-dispersive at every wavenumber, so the theoretical
//! centroid speed is `c` exactly and the tolerance is the discretisation's:
//!
//! | term | size | why |
//! |---|---|---|
//! | meridional truncation | `(Δy/Le)² = 2.1%` | the C-grid operators of T-01.1 are second order, and `Le` is the scale of the structure they differentiate; the `O(1)` coefficient is not evaluated, so this is a bound and not an estimate |
//! | zonal truncation | `(Δx/σ_K)²/6 = 0.05%` | the centred difference's phase error at the pulse's dominant wavenumber |
//! | measurement | `0.5%` | the initial condition is an exact *continuous* mode, not an exact discrete one, and sheds a small transient; the pulse is truncated by the western wall, and by the late sample its leading tail has begun to cross the eastern one |
//!
//! Summed and rounded up: [`KELVIN_SPEED_TOLERANCE`] = 3%.
//!
//! ## The Rossby speed: `c/3`, and a dispersive bias that dominates
//!
//! `c/3` is the *long-wave* limit. The gravest mode's exact dispersion
//! relation, `ω̂³ − (k̂² + 3)·ω̂ − k̂ = 0` in the units above, expands as
//! `ω̂ = −k̂/3 + (8/81)·k̂³ + O(k̂⁵)`, so a packet of finite width travels
//! measurably slower than `c/3`:
//!
//! ```text
//! c_g = −c/3 + (8/27)·c·k̂²,   k̂ = k·Le.
//! ```
//!
//! For the reflected Gaussian of width `σ_R` the energy spectrum is
//! `e^{−k²σ_R²}`, whose mean is `⟨k̂²⟩ = Le²/(2σ_R²)`, leaving a relative bias
//! of `(4/9)·(Le/σ_R)² = 3.3%` — a real physical effect of the finite pulse,
//! not an error of the model, and the largest term here.
//!
//! | term | size | why |
//! |---|---|---|
//! | dispersive bias | `(4/9)·(Le/σ_R)² = 3.3%` | above; the packet is not infinitely long |
//! | next order in `k̂` | `0.33·⟨k̂⁴⟩ = 0.14%` | the `k̂⁵` term of the same expansion, `⟨k̂⁴⟩ = (3/4)(Le/σ_R)⁴` |
//! | meridional truncation | `(Δy/Le)² = 2.1%` | as above |
//! | zonal truncation | `(Δx/σ_R)²/6 = 0.42%` | as above, at the *reflected* packet's width |
//! | measurement | `1%` | the coastal mask clips the packet's trailing edge 3.7 `σ_R` behind its centre; the perimeter wave's own reflection off the western wall puts a little westward energy far behind the packet, worth `≈0.3%` of the centroid at the amplitude the run's first sample shows it reaching |
//!
//! Summed: 7.0%, rounded up to [`ROSSBY_SPEED_TOLERANCE`] = 8% so that the
//! bound does not sit a hair inside a term estimated to one significant
//! figure. Every entry is a property of the configuration below, and none was
//! obtained by running the model.
//!
//! # The convergence this file asserts, and the one it does not
//!
//! CODING_STANDARDS.md § Tests asks that an error be shown to *shrink at the
//! scheme's order across two resolutions* rather than to sit under a fixed
//! threshold. The error that dominates here is not the scheme's: it is the
//! `(4/9)·(Le/σ_R)²` dispersive bias, which is a property of the continuous
//! problem and would not move if the grid were refined to nothing. So the
//! parameter that is varied is the one that bias depends on — the packet's
//! width — and [`the_measured_rossby_speed_converges_on_c_over_three_with_packet_width`]
//! runs the whole experiment a second time with an incident pulse
//! [`WIDTH_REFINEMENT`] times narrower and holds the two errors to the
//! `σ_R^{−2}` rate the dispersion relation predicts. That is a sharper
//! statement than either point check: it tests the *shape* of the error, and
//! it fails if the model's Rossby speed happens to be right for the wrong
//! reason.
//!
//! Refining the grid instead would measure neither the scheme nor the physics.
//! The Kelvin speed's `(Δy/Le)²` entry is a bound with an unevaluated
//! coefficient, and it sits alongside a resolution-*independent* floor of the
//! same order: at the late sample the pulse's leading tail has crossed the
//! wall, removing `≈5×10⁻⁵` of the energy from a lever arm of `2.9·σ_K` and
//! displacing the centroid by `≈0.01%` of the measured speed. A refinement
//! study at this configuration would be reading that floor.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::sync::OnceLock;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing, WaveSpeed,
    WindStressField, H_STAGGERING, U_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;
/// Reference seawater density `ρ₀`, in kg/m³ — `CONTEXT.md` and Gill, appendix 3.
const REFERENCE_DENSITY_KG_PER_M3: f64 = engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3;
/// Rayleigh damping `r` of this experiment, in s⁻¹.
///
/// Zero. Damping decays a wave without moving it, but it decays the incident
/// and the reflected halves of the run by different factors and would leave
/// the reflected packet's amplitude a function of the damping rather than of
/// the reflection; the speeds this file measures are cleanest in the undamped
/// limit, and nothing here needs a steady state to exist.
const UNDAMPED_PER_S: f64 = 0.0;

/// Cells along x of the reflection basin: 26 000 km of it.
///
/// Half again the equatorial Pacific's width (`CONTEXT.md`, *Basin*), and
/// deliberately so: an incident pulse long enough to reflect into a long
/// Rossby wave has to sit clear of both walls, and the reflected wave needs
/// room to run west without meeting the western boundary and reflecting a
/// second time.
const BASIN_NX: usize = 130;
/// Cells along y. Different from [`BASIN_NX`] so an x/y swap cannot pass.
const BASIN_NY: usize = 80;
/// Cell width, in metres: six to a `σ_R`, the narrowest structure the wider
/// run carries, which is what the `(Δx/σ_R)²/6` entry of the Rossby budget
/// costs.
const CELL_WIDTH_M: f64 = 2.0e5;
/// Cell height, in metres: about seven to an equatorial deformation radius,
/// which is what makes the `(Δy/Le)²` term of the two error budgets small.
const CELL_HEIGHT_M: f64 = 5.0e4;
/// Meridional extent of the basin, in metres: ±2000 km, about ±5.8 `Le`.
///
/// Wide enough that the equatorial modes are exponentially negligible at the
/// north and south walls (`e^{−16.8}`), so the closed meridional boundaries
/// take no part in the reflection being measured, and narrow enough to leave
/// the gravity-wave CFL bound the binding one rather than ADR-0007's rotation
/// bound.
const BASIN_LY_M: f64 = BASIN_NY as f64 * CELL_HEIGHT_M;

/// Zonal width `σ_K` of the incident Kelvin pulse, in metres — 11 `Le`.
///
/// The reflection compresses the packet by three (module header), so this is
/// what buys a reflected wave 3.7 `Le` wide and holds the dispersive bias of
/// its speed down to 3.3%.
const INCIDENT_WIDTH_M: f64 = 3.8e6;
/// Factor by which the second run narrows the incident pulse.
///
/// The dispersive bias goes as `σ_R^{−2}`, so this is a factor of 2.25 in the
/// error — far enough apart to measure the rate, close enough that the
/// narrower run's own bias (7.4%) is still a correction to `c/3` rather than a
/// different regime, and small enough that the narrower pulse fits the same
/// basin with more clearance than the wider one.
const WIDTH_REFINEMENT: f64 = 1.5;
/// Distance from the incident pulse's centre to the eastern wall at `t = 0`,
/// in widths.
///
/// Four `σ_K`: the pulse starts with `e^{−8}` of its amplitude at the coast,
/// so the run begins with no reflection under way, and still leaves 2.8 `σ_K`
/// of clearance behind it at the western wall in the wider run.
const INCIDENT_OFFSET_IN_WIDTHS: f64 = 4.0;
/// Peak thermocline depth anomaly of the incident pulse, in metres.
///
/// A downwelling Kelvin pulse of the scale a westerly wind burst leaves behind
/// (`CONTEXT.md`, *Westerly wind burst*). The core is linear, so every speed
/// and ratio here is independent of this number; it only sets the scale the
/// diagnostics are reported in.
const INCIDENT_AMPLITUDE_M: f64 = 10.0;

/// Clearance held between the eastern wall and the westward-mode measurements,
/// in deformation radii.
///
/// The coastal Kelvin wave the reflection also launches is trapped against the
/// wall with an `e^{−(x_E−x)/Le}` profile, so at four radii it is 1.8% in
/// amplitude and 0.03% in energy — below every other term in the budget.
const COASTAL_CLEARANCE_IN_RADII: f64 = 4.0;

/// When the incident pulse's speed is first sampled, in transits of one `σ_K`
/// at `c`. A tenth of the way to the late sample: late enough for the initial
/// condition's transient to have left the pulse, early enough to make the
/// baseline long.
const KELVIN_SAMPLE_EARLY_IN_TRANSITS: f64 = 0.125;
/// When it is sampled again. The pulse's centre is 2.75 `σ_K` short of the
/// eastern wall at this point, so what is measured is the incident wave alone.
const KELVIN_SAMPLE_LATE_IN_TRANSITS: f64 = 1.25;
/// When the reflected packet's speed is first sampled, in the same transits.
///
/// The incident centre reaches the wall after four transits and its trailing
/// `4·σ_K` after eight, so this is three quarters of a transit after the whole
/// reflection has happened, with the reflected centre 4.75 `σ_R` west of the
/// wall and better than three clear of the coastal mask.
const ROSSBY_SAMPLE_EARLY_IN_TRANSITS: f64 = 8.75;
/// When it is sampled again: three and a half transits later, over which the
/// packet travels 3.5 `σ_R` west. Its leading edge is still 4 `σ_K` short of
/// the western wall, so no second reflection has begun.
const ROSSBY_SAMPLE_LATE_IN_TRANSITS: f64 = 12.25;

/// Tolerance on the measured incident Kelvin speed, as a fraction of `c`.
/// Derived in the module header: 2.1% + 0.05% + 0.5%, rounded up.
const KELVIN_SPEED_TOLERANCE: f64 = 0.03;
/// Tolerance on the measured reflected Rossby speed, as a fraction of `c/3`.
/// Derived in the module header: 3.3% + 0.14% + 2.1% + 0.42% + 1%, rounded up.
const ROSSBY_SPEED_TOLERANCE: f64 = 0.08;
/// Tolerance on the rate at which that speed converges on `c/3` as the packet
/// is widened, as a fraction of the predicted `σ_R^{−2}`.
///
/// The predicted rate is the ratio of two leading-order biases, so what it
/// omits is the next order in each: `0.33·⟨k̂⁴⟩` is 4% of the wider run's bias
/// and 9% of the narrower one's, and they do not cancel. 20% is that 13%
/// with the numerical terms of the two budgets — which are common to both runs
/// only to the extent that the two packets have the same width, which is
/// exactly what differs — allowed for on top.
const ROSSBY_CONVERGENCE_TOLERANCE: f64 = 0.2;

/// The `ψ_2`-to-`ψ_0` ratio of the gravest Rossby mode's thermocline anomaly,
/// from `h/H = R·(ψ_0 + ψ_2/4)` (module header).
const GRAVEST_ROSSBY_SHAPE_RATIO: f64 = 0.25;
/// Tolerance on that ratio, as a fraction of itself.
///
/// The two structures share one zonal profile, so the ratio is exact in the
/// long-wave limit and the tolerance is what the corrections to it leave. The
/// `n = 3` mode contributes nothing measurable: it travels at `c/7`, so by the
/// late sample it is 4.7 `σ_R` east of the gravest packet's centre and only
/// `σ_K/7` wide, which puts the window's eastern edge six of its own widths
/// away from it. What is left is the `O(k̂²)` correction to the modal
/// structure at the packet's own wavenumber, `⟨k̂²⟩ = (1/2)(Le/σ_R)² = 3.7%`
/// with an unevaluated coefficient, and the same `(Δy/Le)² = 2.1%` truncation
/// as the speeds: 10% is the first of those doubled, plus the second.
const ROSSBY_SHAPE_TOLERANCE: f64 = 0.1;

/// Zonal half-width of the window the reflected packet's shape is read in, in
/// packet widths. Two `σ_R` either side of its centre holds 95% of it.
const SHAPE_WINDOW_IN_WIDTHS: f64 = 2.0;

/// Largest share of the run's energy the westward mode may hold before the
/// incident pulse has reached the coast.
///
/// The initial condition is an exact continuous Kelvin wave, whose `q` is
/// identically zero. What the discrete grid leaves behind is the same
/// `(Δy/Le)² = 2.1%` truncation the speed budgets carry, in amplitude, and so
/// `4×10⁻⁴` in energy; this is that doubled, for the coefficient the
/// truncation term does not evaluate.
const WESTWARD_SHARE_BEFORE_REFLECTION: f64 = 1.0e-3;
/// Largest share of the run's energy the eastward mode may hold once the
/// reflection is over.
///
/// At the first post-reflection sample the basin holds no eastward wave: the
/// incident one has been absorbed into the wall, the reflected one is
/// travelling west, and the coastal wave that carries the rest of the
/// reflection is still coming round the perimeter (module header). What is
/// left is the same discretisation leakage the westward mode carries before
/// the reflection, so this is bounded the same way and at the same figure.
const EASTWARD_SHARE_AFTER_REFLECTION: f64 = WESTWARD_SHARE_BEFORE_REFLECTION;

/// One of the two meridional structures this file projects a run onto.
///
/// Named rather than indexed by `n`, because `n` is not the interesting thing
/// about either of them: [`MeridionalStructure::Gravest`] is the Kelvin wave's
/// shape and the gravest Rossby mode's leading one, and
/// [`MeridionalStructure::Second`] is the partner that makes the Rossby mode's
/// thermocline anomaly deepest off the equator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MeridionalStructure {
    /// `ψ_0(η) = e^{−η²/2}`.
    Gravest,
    /// `ψ_2(η) = (4η² − 2)·e^{−η²/2}`.
    Second,
}

impl MeridionalStructure {
    /// This structure at `y` metres from the equator, for a waveguide of the
    /// given deformation radius.
    fn at(self, y_m: f64, deformation_radius_m: f64) -> f64 {
        let eta = y_m / deformation_radius_m;
        let envelope = (-0.5 * eta * eta).exp();
        match self {
            Self::Gravest => envelope,
            Self::Second => (4.0 * eta * eta - 2.0) * envelope,
        }
    }
}

/// The two `ψ_0` invariants of a state, column by column: the eastward
/// `P₀[u/c + h/H]` and the westward `P₀[u/c − h/H]` of the module header.
#[derive(Debug, Clone)]
struct Invariants {
    /// Where the Kelvin wave is, and how much of it there is.
    eastward: Vec<f64>,
    /// The same for the gravest Rossby mode.
    westward: Vec<f64>,
}

/// One reflection experiment: a basin, an ocean, and an incident pulse of a
/// given width.
///
/// The width is what the two runs differ in, and everything that depends on it
/// — where the pulse starts, when each sample is taken, how wide the reflected
/// packet is — is derived here rather than written down twice.
#[derive(Debug, Clone, Copy)]
struct Experiment {
    /// Shape, spacing and position of the basin.
    basin: Basin,
    /// The ocean the equations are written in terms of.
    params: PhysicalParams,
    /// `Le = √(c/β)`, in metres — cached because every projection needs it.
    deformation_radius_m: f64,
    /// Zonal width `σ_K` of the incident Kelvin pulse, in metres.
    incident_width_m: f64,
}

impl Experiment {
    /// The experiment with an incident pulse `incident_width_m` wide.
    fn new(incident_width_m: f64) -> Self {
        let params = PhysicalParams::new(
            PACIFIC_REDUCED_GRAVITY_M_PER_S2,
            PACIFIC_MEAN_DEPTH_M,
            UNDAMPED_PER_S,
            BETA_PER_M_PER_S,
            REFERENCE_DENSITY_KG_PER_M3,
        )
        .expect("the published equatorial-Pacific parameters are physical");
        let grid = Grid::new(BASIN_NX, BASIN_NY).expect("extents are non-zero");
        let spacing = Spacing::new(CELL_WIDTH_M, CELL_HEIGHT_M).expect("cell sizes are positive");
        let basin = Basin::new(grid, spacing, 0.0, -0.5 * BASIN_LY_M)
            .expect("both edges are finite positions");
        Self {
            basin,
            params,
            deformation_radius_m: (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s())
                .sqrt(),
            incident_width_m,
        }
    }

    /// The Kelvin wave speed `c = √(g'·H)`, in m/s (`CONTEXT.md`).
    fn wave_speed_m_per_s(self) -> f64 {
        self.params.kelvin_wave_speed_m_per_s()
    }

    /// Zonal width `σ_R = σ_K/3` of the reflected packet, in metres: the
    /// compression the module header derives.
    fn reflected_width_m(self) -> f64 {
        self.incident_width_m / 3.0
    }

    /// How long the incident pulse takes to travel one of its own widths, in
    /// seconds — the clock every sample time is stated in, so that the two
    /// runs sample the same geometry rather than the same calendar.
    fn transit_s(self) -> f64 {
        self.incident_width_m / self.wave_speed_m_per_s()
    }

    /// Zonal position of the incident pulse's centre at `t = 0`, in metres.
    fn incident_centre_x_m(self) -> f64 {
        self.basin.zonal_extent_m() - INCIDENT_OFFSET_IN_WIDTHS * self.incident_width_m
    }

    /// The easternmost column position the westward measurements read, in
    /// metres: the coastal Kelvin wave's clearance.
    fn coastal_mask_x_m(self) -> f64 {
        self.basin.zonal_extent_m() - COASTAL_CLEARANCE_IN_RADII * self.deformation_radius_m
    }

    /// The four sample times of the run, in seconds.
    fn sample_times_s(self) -> [f64; 4] {
        [
            KELVIN_SAMPLE_EARLY_IN_TRANSITS,
            KELVIN_SAMPLE_LATE_IN_TRANSITS,
            ROSSBY_SAMPLE_EARLY_IN_TRANSITS,
            ROSSBY_SAMPLE_LATE_IN_TRANSITS,
        ]
        .map(|transits| transits * self.transit_s())
    }

    /// `ψ_n` sampled on the cell-center rows, which is where both `h` and `u`
    /// sit meridionally.
    fn row_structure(self, structure: MeridionalStructure) -> Vec<f64> {
        (0..self.basin.grid().ny())
            .map(|j| {
                structure.at(
                    self.basin.y_of_row_m(H_STAGGERING, j),
                    self.deformation_radius_m,
                )
            })
            .collect()
    }

    /// The `ψ_n` coefficient, column by column, of a cell-centered field given
    /// as `value(i, j)`.
    ///
    /// The `ψ_n` are orthogonal on the line and the basin reaches 5.8 `Le`, so
    /// this discrete inner product is the modal coefficient to the accuracy of
    /// the row quadrature — the row spacing cancels between the projection and
    /// its normalisation.
    fn project_columns(
        self,
        structure: MeridionalStructure,
        value: impl Fn(usize, usize) -> f64,
    ) -> Vec<f64> {
        let weights = self.row_structure(structure);
        let normalisation: f64 = weights.iter().map(|weight| weight * weight).sum();
        (0..self.basin.grid().nx())
            .map(|i| {
                let projection: f64 = weights
                    .iter()
                    .enumerate()
                    .map(|(j, weight)| value(i, j) * weight)
                    .sum();
                projection / normalisation
            })
            .collect()
    }

    /// The thermocline depth anomaly's `ψ_n` coefficient, column by column, in
    /// units of the mean depth `H`.
    fn depth_projection(self, state: &OceanState, structure: MeridionalStructure) -> Vec<f64> {
        self.project_columns(structure, |i, j| {
            state.h().get(i, j).expect("a cell center") / self.params.mean_thermocline_depth_m()
        })
    }

    /// The two `ψ_0` invariants of `state` (module header).
    ///
    /// `u` is averaged from its two faces onto the cell center so that both
    /// invariants are read at one set of positions.
    fn invariants(self, state: &OceanState) -> Invariants {
        let zonal_current = self.project_columns(MeridionalStructure::Gravest, |i, j| {
            let west = state.u().get(i, j).expect("an east/west face");
            let east = state.u().get(i + 1, j).expect("an east/west face");
            0.5 * (west + east) / self.wave_speed_m_per_s()
        });
        let depth_anomaly = self.depth_projection(state, MeridionalStructure::Gravest);
        Invariants {
            eastward: zonal_current
                .iter()
                .zip(&depth_anomaly)
                .map(|(current, depth)| current + depth)
                .collect(),
            westward: zonal_current
                .iter()
                .zip(&depth_anomaly)
                .map(|(current, depth)| current - depth)
                .collect(),
        }
    }

    /// The columns of `profile` west of `eastern_limit_m`, as
    /// `(position in metres, amplitude)`.
    fn columns_west_of(
        self,
        profile: &[f64],
        eastern_limit_m: f64,
    ) -> impl Iterator<Item = (f64, f64)> + '_ {
        profile
            .iter()
            .enumerate()
            .map(move |(i, amplitude)| (self.basin.x_of_column_m(H_STAGGERING, i), *amplitude))
            .filter(move |(x_m, _)| *x_m <= eastern_limit_m)
    }

    /// The energy-weighted zonal centroid of `profile`, in metres, over the
    /// columns west of the coastal mask.
    ///
    /// The centroid of a linear wave packet moves at its energy-weighted mean
    /// group velocity, which is what both speed measurements read.
    fn energy_centroid_m(self, profile: &[f64]) -> f64 {
        let (weighted_position, weight) = self
            .columns_west_of(profile, self.coastal_mask_x_m())
            .fold((0.0, 0.0), |(moment, total), (x_m, amplitude)| {
                let energy = amplitude * amplitude;
                (moment + x_m * energy, total + energy)
            });
        assert!(
            weight > 0.0,
            "the profile carries no energy west of the coastal mask, so it has no centroid"
        );
        weighted_position / weight
    }

    /// `Σ profile²` over the same columns — the quantity that centroid weights
    /// by, and the one the before-and-after energy shares compare.
    fn energy(self, profile: &[f64]) -> f64 {
        self.columns_west_of(profile, self.coastal_mask_x_m())
            .map(|(_, amplitude)| amplitude * amplitude)
            .sum()
    }

    /// The incident Kelvin pulse: a Gaussian in `x` on the `ψ_0` waveguide,
    /// with `u = (c/H)·h` and `v = 0`.
    ///
    /// An exact solution of the continuous equations for any zonal profile
    /// (module header), so the run starts with one wave in it and no Rossby
    /// energy at all.
    fn initial_state(self) -> OceanState {
        let mut state = OceanState::at_rest(self.basin.grid());
        let centre_x_m = self.incident_centre_x_m();
        let current_amplitude_m_per_s = INCIDENT_AMPLITUDE_M * self.wave_speed_m_per_s()
            / self.params.mean_thermocline_depth_m();
        let profile = |x_m: f64| {
            let offset = (x_m - centre_x_m) / self.incident_width_m;
            (-0.5 * offset * offset).exp()
        };

        for j in 0..state.h().ny() {
            let waveguide = MeridionalStructure::Gravest.at(
                self.basin.y_of_row_m(H_STAGGERING, j),
                self.deformation_radius_m,
            );
            for i in 0..state.h().nx() {
                let x_m = self.basin.x_of_column_m(H_STAGGERING, i);
                *state.h_mut().get_mut(i, j).expect("a cell center") =
                    INCIDENT_AMPLITUDE_M * profile(x_m) * waveguide;
            }
            for i in 0..state.u().nx() {
                let x_m = self.basin.x_of_column_m(U_STAGGERING, i);
                *state.u_mut().get_mut(i, j).expect("an east/west face") =
                    current_amplitude_m_per_s * profile(x_m) * waveguide;
            }
        }
        state
    }

    /// Run the experiment: the incident pulse in a closed, unforced, undamped
    /// basin, sampled twice before it strikes the eastern wall and twice after
    /// it has reflected.
    fn run(self) -> Reflection {
        let wave_speed = WaveSpeed::new(self.wave_speed_m_per_s()).expect("a positive wave speed");
        let dt_s = max_stable_dt(self.basin.spacing(), wave_speed);
        let plane =
            BetaPlane::centered_on_equator(self.params, self.basin.spacing(), self.basin.grid());
        let mut solver = Solver::new(
            self.basin.grid(),
            self.basin.spacing(),
            self.params,
            plane,
            dt_s,
        )
        .unwrap_or_else(|error| {
            panic!("the experiment's own timestep must be admissible: {error}")
        });
        let calm = WindStressField::calm(self.basin.grid());
        let sample_steps = self
            .sample_times_s()
            .map(|t_s| (t_s / dt_s).round() as usize);

        let mut state = self.initial_state();
        let mut samples: Vec<Sample> = Vec::with_capacity(sample_steps.len());
        let last_step = sample_steps[sample_steps.len() - 1];
        for step in 0..=last_step {
            if sample_steps.contains(&step) {
                let invariants = self.invariants(&state);
                samples.push(Sample {
                    t_s: step as f64 * dt_s,
                    eastward_centroid_m: self.energy_centroid_m(&invariants.eastward),
                    westward_centroid_m: self.energy_centroid_m(&invariants.westward),
                    eastward_energy: self.energy(&invariants.eastward),
                    westward_energy: self.energy(&invariants.westward),
                    depth_in_gravest: self.depth_projection(&state, MeridionalStructure::Gravest),
                    depth_in_second: self.depth_projection(&state, MeridionalStructure::Second),
                });
            }
            if step < last_step {
                solver.step(&mut state, step as f64 * dt_s, |_t| &calm);
            }
        }

        let mut samples = samples.into_iter();
        let mut next = || samples.next().expect("every requested sample was taken");
        Reflection {
            incident: [next(), next()],
            reflected: [next(), next()],
        }
    }
}

/// One sample of a run: the two invariants' centroids and energies, and the
/// thermocline anomaly's two meridional coefficients, at a given time.
#[derive(Debug, Clone)]
struct Sample {
    /// When it was taken, in seconds — the step's time, not the requested one.
    t_s: f64,
    /// Centroid of `P₀[r]`, in metres: where the Kelvin wave is.
    eastward_centroid_m: f64,
    /// Centroid of `P₀[q]`, in metres: where the gravest Rossby wave is.
    westward_centroid_m: f64,
    /// `Σ P₀[r]²` west of the coastal mask.
    eastward_energy: f64,
    /// `Σ P₀[q]²` west of the coastal mask.
    westward_energy: f64,
    /// The `ψ_0` coefficient of `h/H`, column by column.
    depth_in_gravest: Vec<f64>,
    /// The `ψ_2` coefficient of `h/H`, column by column.
    depth_in_second: Vec<f64>,
}

/// What one run of the experiment leaves for the tests to read.
#[derive(Debug, Clone)]
struct Reflection {
    /// The incident pulse, sampled twice before it reaches the coast.
    incident: [Sample; 2],
    /// The reflected packet, sampled twice after the reflection is over.
    reflected: [Sample; 2],
}

/// The wider experiment, run once and shared by every test in this file.
///
/// Four tests read four different things out of one trajectory, so running it
/// once is both faster and stronger than running it four times: every
/// assertion is made about the same reflection.
fn wide_reflection() -> &'static Reflection {
    static REFLECTION: OnceLock<Reflection> = OnceLock::new();
    REFLECTION.get_or_init(|| Experiment::new(INCIDENT_WIDTH_M).run())
}

/// The same experiment with an incident pulse [`WIDTH_REFINEMENT`] times
/// narrower — the second point of the convergence test.
fn narrow_reflection() -> &'static Reflection {
    static REFLECTION: OnceLock<Reflection> = OnceLock::new();
    REFLECTION.get_or_init(|| Experiment::new(INCIDENT_WIDTH_M / WIDTH_REFINEMENT).run())
}

/// The speed a centroid moved between two samples, in m/s — positive eastward.
fn centroid_speed_m_per_s(from: &Sample, to: &Sample, centroid: impl Fn(&Sample) -> f64) -> f64 {
    (centroid(to) - centroid(from)) / (to.t_s - from.t_s)
}

/// How far the reflected packet's measured speed sits from `c/3`, as a
/// fraction of it.
fn rossby_speed_error(experiment: Experiment, reflection: &Reflection) -> f64 {
    let expected_m_per_s = -experiment.wave_speed_m_per_s() / 3.0;
    let measured_m_per_s = centroid_speed_m_per_s(
        &reflection.reflected[0],
        &reflection.reflected[1],
        |sample| sample.westward_centroid_m,
    );
    (measured_m_per_s - expected_m_per_s).abs() / expected_m_per_s.abs()
}

#[test]
fn the_incident_pulse_crosses_the_basin_eastward_at_the_kelvin_wave_speed() {
    // The first half of the ticket: the wave that strikes the eastern wall is
    // a Kelvin wave, travelling east at `c = √(g'·H)` (`CONTEXT.md`). Measured
    // on `P₀[u/c + h/H]`, which no Rossby mode contributes to, between two
    // samples taken before the pulse reaches the coast.
    let reflection = wide_reflection();
    let expected_m_per_s = Experiment::new(INCIDENT_WIDTH_M).wave_speed_m_per_s();

    let measured_m_per_s =
        centroid_speed_m_per_s(&reflection.incident[0], &reflection.incident[1], |sample| {
            sample.eastward_centroid_m
        });

    let error = (measured_m_per_s - expected_m_per_s).abs() / expected_m_per_s;
    assert!(
        error <= KELVIN_SPEED_TOLERANCE,
        "the incident pulse travelled at {measured_m_per_s} m/s, {:.1}% from the Kelvin wave \
         speed {expected_m_per_s} m/s; the budget allows {:.1}%",
        100.0 * error,
        100.0 * KELVIN_SPEED_TOLERANCE
    );
}

#[test]
fn the_eastern_wall_turns_the_eastward_wave_into_a_westward_one() {
    // What "reflects incident Kelvin energy back as Rossby waves"
    // (`01-scientific-model.md`) means for the two invariants: before the
    // pulse reaches the coast the westward one is empty, and after the
    // reflection the eastward one is.
    let reflection = wide_reflection();

    let before = &reflection.incident[1];
    let westward_share = before.westward_energy / (before.eastward_energy + before.westward_energy);
    assert!(
        westward_share < WESTWARD_SHARE_BEFORE_REFLECTION,
        "before the pulse reached the coast the westward mode already held {:.4}% of the energy, \
         past the {:.2}% the discretisation can leak",
        100.0 * westward_share,
        100.0 * WESTWARD_SHARE_BEFORE_REFLECTION
    );

    // The first post-reflection sample, not the second: by the second the
    // boundary-bound half of the reflection has come round the closed basin's
    // perimeter and re-entered the wave guide in the west (module header),
    // which puts an eastward wave back in the basin without any of it having
    // failed to reflect.
    let after = &reflection.reflected[0];
    let eastward_share = after.eastward_energy / (after.eastward_energy + after.westward_energy);
    assert!(
        eastward_share < EASTWARD_SHARE_AFTER_REFLECTION,
        "after the reflection the eastward mode still held {:.4}% of the energy; the incident \
         wave should have been turned into the westward one",
        100.0 * eastward_share
    );
}

#[test]
fn the_reflected_packet_travels_west_at_a_third_of_the_kelvin_wave_speed() {
    // The ticket's acceptance criterion: the reflected signal's propagation
    // speed matches theory — `c/3` westward for the gravest meridional Rossby
    // mode (`CONTEXT.md`, *Rossby wave*) — within the tolerance the module
    // header derives.
    let experiment = Experiment::new(INCIDENT_WIDTH_M);
    let reflection = wide_reflection();
    let expected_m_per_s = -experiment.wave_speed_m_per_s() / 3.0;

    let measured_m_per_s = centroid_speed_m_per_s(
        &reflection.reflected[0],
        &reflection.reflected[1],
        |sample| sample.westward_centroid_m,
    );

    assert!(
        measured_m_per_s < 0.0,
        "the reflected packet travelled east at {measured_m_per_s} m/s; a Rossby wave goes west"
    );
    let error = rossby_speed_error(experiment, reflection);
    assert!(
        error <= ROSSBY_SPEED_TOLERANCE,
        "the reflected packet travelled at {measured_m_per_s} m/s, {:.1}% from the gravest \
         Rossby mode's {expected_m_per_s} m/s; the budget allows {:.1}%",
        100.0 * error,
        100.0 * ROSSBY_SPEED_TOLERANCE
    );
}

#[test]
fn the_measured_rossby_speed_converges_on_c_over_three_with_packet_width() {
    // `c/3` is the long-wave limit, and what separates the measurement from it
    // is the dispersive bias `(4/9)·(Le/σ_R)²` of the module header. That is a
    // statement about the *shape* of the error, not merely its size, so it is
    // checked the way CODING_STANDARDS.md § Tests asks an order of accuracy to
    // be checked: at two widths, against the rate rather than a threshold.
    let wide = Experiment::new(INCIDENT_WIDTH_M);
    let narrow = Experiment::new(INCIDENT_WIDTH_M / WIDTH_REFINEMENT);

    let wide_error = rossby_speed_error(wide, wide_reflection());
    let narrow_error = rossby_speed_error(narrow, narrow_reflection());

    // `σ_R^{−2}`: narrowing the packet by a factor raises the bias by its
    // square.
    let expected_ratio = WIDTH_REFINEMENT * WIDTH_REFINEMENT;
    let measured_ratio = narrow_error / wide_error;

    let error = (measured_ratio - expected_ratio).abs() / expected_ratio;
    assert!(
        error <= ROSSBY_CONVERGENCE_TOLERANCE,
        "narrowing the incident pulse by {WIDTH_REFINEMENT} multiplied the speed error by \
         {measured_ratio:.2} ({:.1}% of it at the wider width, {:.1}% at the narrower), but the \
         dispersion relation predicts {expected_ratio:.2}; the budget allows {:.0}%",
        100.0 * wide_error,
        100.0 * narrow_error,
        100.0 * ROSSBY_CONVERGENCE_TOLERANCE
    );
}

#[test]
fn the_reflected_packet_carries_the_gravest_rossby_modes_meridional_shape() {
    // Speed alone does not identify a mode, so this is the other half of the
    // claim: the reflected thermocline anomaly has the `ψ_0 + ψ_2/4` structure
    // the eastern wall's `u = 0` condition produces (module header) — deepest
    // off the equator rather than on it.
    let experiment = Experiment::new(INCIDENT_WIDTH_M);
    let late = &wide_reflection().reflected[1];

    // Read the shape where the packet is, so the slower `n = 3` mode — still
    // 2.3 times closer to the coast — stays out of the window.
    let window_m = SHAPE_WINDOW_IN_WIDTHS * experiment.reflected_width_m();
    let centre_m = late.westward_centroid_m;
    let (cross_term, gravest_energy) = (0..experiment.basin.grid().nx())
        .filter(|i| (experiment.basin.x_of_column_m(H_STAGGERING, *i) - centre_m).abs() <= window_m)
        .fold((0.0, 0.0), |(cross, gravest), i| {
            (
                cross + late.depth_in_gravest[i] * late.depth_in_second[i],
                gravest + late.depth_in_gravest[i] * late.depth_in_gravest[i],
            )
        });
    let ratio = cross_term / gravest_energy;

    let error = (ratio - GRAVEST_ROSSBY_SHAPE_RATIO).abs() / GRAVEST_ROSSBY_SHAPE_RATIO;
    assert!(
        error <= ROSSBY_SHAPE_TOLERANCE,
        "the reflected packet's ψ₂/ψ₀ ratio is {ratio}, {:.1}% from the \
         {GRAVEST_ROSSBY_SHAPE_RATIO} of the gravest Rossby mode; the budget allows {:.0}%",
        100.0 * error,
        100.0 * ROSSBY_SHAPE_TOLERANCE
    );
}
