//! Acceptance test for T-04.3 — reflection at the closed western boundary.
//!
//! `docs/planning/01-scientific-model.md` § *Domain and boundaries* states the
//! process this file validates: the **western boundary** "reflects incident
//! Rossby energy partly back as Kelvin waves". The mirror image — incident
//! Kelvin energy returned as Rossby waves — happens at the *eastern* boundary
//! and belongs to T-04.4.
//!
//! # Why this is not the reflection issue #20 first asked for
//!
//! Issue #20 originally described the eastern process at the western wall: a
//! Kelvin pulse reflecting into Rossby energy "propagating back eastward at
//! c/3". No arrangement of the physics makes that statement true. A Kelvin
//! wave travels eastward only (`CONTEXT.md`, *Kelvin wave*), so it never
//! reaches the western wall; Rossby waves travel westward (`CONTEXT.md`,
//! *Rossby wave*), and `c/3` is the westward long-wave group speed of the
//! gravest meridional mode, so Rossby energy launched from a western wall
//! would leave the basin rather than propagate into it. T-04.3 and T-04.4 had
//! their wave pairs swapped relative to their boundary labels.
//!
//! That contradiction was **escalated before any of this was written**, not
//! resolved here: AGENTS.md § *Never move the goalposts* makes a ticket whose
//! acceptance criteria look wrong a human decision, not a judgement call for
//! the agent implementing it. The ruling was to keep this ticket's boundary
//! and correct its wave pair — an incident gravest-mode Rossby packet
//! reflecting as an eastward Kelvin wave — and to move the Kelvin→Rossby case
//! to T-04.4, where the eastern boundary already is. Issue #20 is retitled
//! *Western boundary Rossby→Kelvin reflection validation* to match. Nothing
//! below weakens the acceptance criterion the issue states: "reflected
//! signal's propagation speed matches theory within a documented, justified
//! tolerance" is exactly what is asserted, with the reflected signal being a
//! Kelvin wave and its theoretical speed `c`.
//!
//! # The analytic prediction
//!
//! Everything asserted below comes from equatorial long-wave theory (Matsuno
//! 1966; Gill, *Atmosphere–Ocean Dynamics*, § 11.6; Cane & Sarachik 1976–81),
//! not from running this engine.
//!
//! Write `ŷ = y/Le` for the meridional coordinate scaled by the equatorial
//! deformation radius `Le = √(c/β)`, and use the two dimensionless Riemann
//! variables
//!
//! ```text
//! q = u/c + h/H          r = u/c − h/H
//! ```
//!
//! In the long-wave limit the meridional momentum equation reduces to
//! geostrophy, `ŷ·u/c = −∂(h/H)/∂ŷ` in these units. Substituting
//! `u/c = (q + r)/2` and `h/H = (q − r)/2` turns that into
//!
//! ```text
//! (∂q/∂ŷ + ŷ·q) − (∂r/∂ŷ − ŷ·r) = 0
//! ```
//!
//! The Hermite functions `ψₘ(ŷ) = Hₘ(ŷ)·exp(−ŷ²/2)` satisfy the ladder
//! relations `(∂/∂ŷ + ŷ)ψₘ = 2m·ψₘ₋₁` and `(∂/∂ŷ − ŷ)ψₘ = −ψₘ₊₁`, so a
//! solution with `q = q̂·ψₙ₊₁` forces `r = −2(n+1)·q̂·ψₙ₋₁`. That single family
//! contains every wave this test uses:
//!
//! - **Kelvin wave** (`n = −1`): `q ∝ ψ₀`, `r = 0`, so `u/c = h/H`, and the
//!   phase speed is `+c` — eastward and non-dispersive.
//! - **Gravest Rossby mode** (`n = 1`): `q = A·ψ₂`, `r = −4A·ψ₀`, and the long
//!   wave speed is `−c/(2n+1) = −c/3` — westward, the `c/3` of `CONTEXT.md`.
//!
//! Reading the two back into the prognostic variables gives the initial
//! condition [`incident_rossby_packet`] builds, with a Gaussian zonal
//! envelope `E(x)`:
//!
//! ```text
//! h/H = A·E(x)·(2ŷ² + 1)·exp(−ŷ²/2)
//! u/c = A·E(x)·(2ŷ² − 3)·exp(−ŷ²/2)
//! v   = 0
//! ```
//!
//! The `h` of that mode is double-lobed off the equator — the familiar
//! off-equatorial Rossby signature — while the reflected Kelvin wave is
//! single-lobed *on* the equator, which is one of the checks below.
//!
//! ## What the wall does
//!
//! The Kelvin wave is the only eastward-propagating long wave, so it is the
//! only long wave a western boundary can radiate. Its amplitude follows from
//! requiring no net zonal mass flux through the wall, `∫u dy = 0`. With
//! `∫ψ₀ dŷ = √(2π)` and `∫ψ₂ dŷ = 2√(2π)`, the incident mode carries
//! `∫u_R dŷ = c·(2A − 4A)√(2π)/2 = −A√(2π)·c` and a Kelvin wave of
//! `q`-coefficient `2a_K` carries `∫u_K dŷ = a_K√(2π)·c`, so `a_K = A`: the
//! reflected Kelvin wave's coefficient on `ψ₀` in `q` is `2A`. Comparing
//! energy fluxes (`∫(u² + h²)dy` times the group speed) the Kelvin wave
//! carries `2A²√π·c` out of the `12A²√π·(c/3) = 4A²√π·c` that came in — half
//! of it. The rest goes into short Rossby waves, which are dispersive, have
//! eastward group velocity, and are absent from long-wave theory; that is
//! exactly the "**partly** back as Kelvin waves" of the scientific model doc,
//! and it is why the amplitude below is bounded loosely while the speeds are
//! not.
//!
//! # How a wave is measured
//!
//! `q` is projected onto the Hermite functions column by column. The
//! projections are orthogonal (`∫ψₗψₘ dŷ = 0` for `l ≠ m`), so the `ψ₀`
//! coefficient of `q` is the Kelvin amplitude with the incident Rossby mode
//! exactly removed, and the `ψ₂` coefficient is the gravest Rossby amplitude
//! with the Kelvin wave exactly removed. Each is recorded at fixed zonal
//! stations, and a speed is a time of flight between two of them: the peak of
//! the packet passes station `a` at `t_a` and station `b` at `t_b`, and the
//! speed is `(x_b − x_a)/(t_b − t_a)`. Peak times are refined by fitting a
//! parabola through the three samples around the largest one.
//!
//! Time of flight is the right measurement for a *non-dispersive* wave: the
//! reflected Kelvin waveform translates unchanged at `c` whatever its shape,
//! so the peak separation is `Δx/c` exactly, and neither the shape of the
//! incident packet nor the details of the reflection can bias it. Only
//! numerical dispersion can, which is what makes the convergence test below
//! meaningful.
//!
//! # Where the tolerances come from
//!
//! No tolerance here is a number that happened to make a test pass; each is
//! computed from the run's own resolution and packet width.
//!
//! - **Reflected Kelvin speed**: `(Δy/Le)²`, the leading truncation error of
//!   the second-order meridional discretisation of the waveguide. `ψ₀` varies
//!   on the scale `Le` itself, so its `O(1)` constant is taken as one; the
//!   general form is [`MeridionalMode::truncation_richness`]. The zonal
//!   contribution is `(k·Δx)²/24 ≈ 5×10⁻⁵` at this packet's dominant
//!   wavenumber and the RK4 time error is fourth order, so both are far below
//!   it. The point check is therefore generous on purpose; the *order* is
//!   pinned by [`reflected_kelvin_speed_converges_at_second_order`], per
//!   CODING_STANDARDS.md § *Convergence over point checks*.
//! - **Incident Rossby speed**: the same truncation, but `ψ₂` oscillates on
//!   the scale `Le/√5` rather than `Le`, so its bound is `5·(Δy/Le)²` — see
//!   [`MeridionalMode::truncation_richness`]. On top of it, and *one-sided*,
//!   the bias a finite packet carries: the gravest mode's group speed is
//!   `c_g/c = −(3 − κ²)/(κ² + 3)²` with `κ = k·Le`, which is `−(1/3)(1 − κ²)`
//!   to leading order, and a Gaussian envelope of e-folding half-width `σ` has
//!   `⟨κ²⟩ = Le²/(2σ²)`. Theory says which way that goes — a finite packet is
//!   *slower* than the long-wave limit, never faster — so it widens only the
//!   slow side of the band.
//!
//!   That band is wide, so the test does not lean on it alone: it also asserts
//!   the measured speed is nearer `c/3` than either neighbouring analytic
//!   speed, which is what makes it the *gravest* mode rather than the `c/5` of
//!   `n = 2` or the Kelvin branch.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::sync::OnceLock;

use engine::{
    max_stable_dt, Basin, BetaPlane, CGridOperators, Grid, OceanState, PhysicalParams, Solver,
    Spacing, WaveSpeed, WindStressField, H_STAGGERING, U_STAGGERING,
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
/// Rayleigh damping `r` of this validation, in s⁻¹. The reflection is a
/// statement about *speeds*, and the analytic wave speeds above are the
/// undamped ones; damping would decay the packet without moving it, so it is
/// switched off rather than corrected for.
const UNDAMPED_PER_S: f64 = 0.0;

/// Cell size of the coarse run, in metres — square, so that neither axis is
/// resolved at the other's expense. At `Le ≈ 345 km` this is 3.45 cells per
/// equatorial deformation radius: coarse enough that the meridional
/// truncation error dominates and the convergence test has something to
/// measure, fine enough that the waveguide is represented at all.
const COARSE_CELL_SIZE_M: f64 = 1.0e5;
/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Long enough that the reflected
/// Kelvin wave clears both stations before it reaches the eastern wall.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: ±1400 km, just over four
/// equatorial deformation radii either side of the equator.
///
/// That is where the walls have to be for them not to matter. `|ψ₂|` — the
/// broadest structure in play — peaks at `ŷ = √2.5` with the value
/// `8·exp(−1.25) = 2.29`, and at the wall `ŷ = 1400/345 = 4.06` it is down to
/// `(4ŷ² − 2)·exp(−ŷ²/2) = 0.017`, under a hundredth of its peak. The
/// incident mode is effectively unaware of the meridional boundaries, which is
/// what lets the analytic speeds above be compared against at all.
const BASIN_LY_M: f64 = 2.8e6;

/// Zonal position of the incident packet's centre at `t = 0`, in metres east
/// of the western wall. Three e-folding widths clear of it, so the wall starts
/// undisturbed.
const PULSE_CENTER_X_M: f64 = 4.0e6;
/// Zonal e-folding half-width `σ` of the packet's Gaussian envelope, in
/// metres. Wide compared with `Le` so the long-wave approximation the `c/3`
/// prediction rests on holds; narrow enough that the packet, and the reflected
/// Kelvin wave three times longer than it, fit inside the basin.
const PULSE_ZONAL_WIDTH_M: f64 = 1.0e6;
/// Thermocline depth anomaly at the packet's centre, on the equator, in
/// metres. The model is linear, so this only sets the scale everything else is
/// measured against.
const PULSE_EQUATORIAL_AMPLITUDE_M: f64 = 10.0;

/// Length of the run, in seconds — about 82 days.
///
/// The schedule it has to cover, at `c = 2.74 m/s` and `c/3 = 0.913 m/s`: the
/// incident packet's centre reaches the wall at `4.0×10⁶ / (c/3) = 4.4×10⁶ s`,
/// and the Kelvin wave it launches passes the two reflection stations at
/// `4.9×10⁶ s` and `6.0×10⁶ s`. The Kelvin wave's *leading edge* first touches
/// the eastern wall at about `5.8×10⁶ s`, and the Rossby waves that reflection
/// sends back are still 4000 km east of the stations when the run ends.
const RUN_DURATION_S: f64 = 7.1e6;
/// Fraction of the CFL-stable maximum this run's timestep takes. Below one so
/// that halving the cell size halves the timestep exactly, which is what makes
/// the convergence test a statement about the spatial discretisation.
const CFL_FRACTION: f64 = 0.9;

/// Zonal position of the western station where the reflected Kelvin wave is
/// timed, in metres. Clear of the wall, where the reflection is still forming.
const KELVIN_STATION_WEST_X_M: f64 = 1.0e6;
/// Zonal position of the eastern station where the reflected Kelvin wave is
/// timed, in metres. 3500 km of flight is about 15 timesteps' worth of packet
/// width, so the time of flight is well resolved.
const KELVIN_STATION_EAST_X_M: f64 = 4.5e6;
/// Zonal position of the eastern station where the incident Rossby packet is
/// timed, in metres — passed early, while the packet is still intact.
const ROSSBY_STATION_EAST_X_M: f64 = 3.4e6;
/// Zonal position of the western station where the incident Rossby packet is
/// timed, in metres. Kept well east of the wall so that the short Rossby waves
/// the reflection radiates, which stay trapped near it, do not contaminate the
/// timing.
const ROSSBY_STATION_WEST_X_M: f64 = 1.6e6;

/// Index of the western reflection station in a [`ReflectionRun`]'s records.
const KELVIN_WEST: usize = 0;
/// Index of the eastern reflection station.
const KELVIN_EAST: usize = 1;
/// Index of the eastern station the incident packet passes first.
const ROSSBY_EAST: usize = 2;
/// Index of the western station the incident packet passes second.
const ROSSBY_WEST: usize = 3;

/// Smallest convergence order the reflected Kelvin speed's error must show
/// when the cell size is halved.
///
/// The spatial discretisation is second order (ADR-0003), so the measured
/// order `log₂(coarse error / fine error)` should be 2. Requiring 1.5 leaves
/// margin for the sub-dominant terms the error also carries — fourth-order
/// time error, a zonal truncation four orders of magnitude smaller — while
/// still failing a first-order scheme, which is the point of asserting an
/// order at all rather than a bare shrinkage.
const MIN_CONVERGENCE_ORDER: f64 = 1.5;

/// Smallest share of the analytic Kelvin amplitude the reflection must return,
/// and the largest it may.
///
/// Long-wave theory pins the reflected `ψ₀` coefficient at `2A` (see the
/// module header), but it is a *long-wave* prediction: the engine enforces
/// `u = 0` pointwise rather than `∫u dy = 0`, so it also radiates short Rossby
/// waves that the theory omits, and the incident packet disperses on its way
/// to the wall. Both take amplitude away from the reflected Kelvin wave
/// without touching its speed. The band is therefore wide on purpose — its job
/// is to say that the reflection returns a *substantial part* of the incident
/// signal and not all of it, which is the "partly" of the scientific model
/// doc, not to pin a conversion ratio theory does not pin here.
const KELVIN_AMPLITUDE_BAND: (f64, f64) = (0.3, 1.5);

/// Factor by which the `ψ₀` content of the reflected signal must exceed its
/// `ψ₂` content for it to be called a Kelvin wave.
///
/// Theory says the reflected long wave is *purely* `ψ₀`, so the honest
/// assertion is qualitative: whatever residual Rossby structure the discrete
/// wall leaves behind, the reflected signal is dominated by the Kelvin mode.
const KELVIN_DOMINANCE: f64 = 2.0;

/// Largest `ψ₀` content the *initial condition* may carry at the western
/// reflection station, as a fraction of the reflected wave's peak there.
///
/// The packet is built on `ψ₂`, which is orthogonal to `ψ₀`, so the run starts
/// with exactly no Kelvin content in exact arithmetic; what survives is the
/// midpoint rule's quadrature error on a Gaussian-decaying integrand, which
/// falls faster than any power of `Δy`. One percent is a generous ceiling on
/// that, and asserting it is what makes the reflected wave *the wall's doing*:
/// a bug that seeded Kelvin energy at `t = 0` would otherwise reproduce the
/// measured speed for entirely the wrong reason.
const SEEDED_KELVIN_CEILING: f64 = 0.01;

/// The equatorial-Pacific parameter set this validation runs in.
fn pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        UNDAMPED_PER_S,
        BETA_PER_M_PER_S,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// One of the two meridional structures this validation separates `q` into.
///
/// Each names a wave by the Hermite function its `q` sits on, which is what
/// makes the two measurements independent: the Hermite functions are
/// orthogonal, so projecting onto one removes the other exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MeridionalMode {
    /// The equatorially trapped Kelvin wave, `q ∝ ψ₀` — the wave the western
    /// boundary radiates.
    Kelvin,
    /// The gravest meridional Rossby mode, `n = 1`. The ladder relations put
    /// its `q` on `ψₙ₊₁ = ψ₂`; it is the wave that arrives.
    GravestRossby,
}

impl MeridionalMode {
    /// Index `m` of the Hermite function `ψₘ` this mode's `q` sits on.
    const fn hermite_order(self) -> usize {
        match self {
            Self::Kelvin => 0,
            Self::GravestRossby => 2,
        }
    }

    /// The physicists' Hermite polynomial `Hₘ(ŷ)` of that order — the
    /// polynomial the equatorial wave problem separates in.
    fn hermite_polynomial(self, y_over_le: f64) -> f64 {
        match self {
            Self::Kelvin => 1.0,
            Self::GravestRossby => 4.0 * y_over_le * y_over_le - 2.0,
        }
    }

    /// The Hermite function `ψₘ(ŷ) = Hₘ(ŷ)·exp(−ŷ²/2)`.
    fn hermite_function(self, y_over_le: f64) -> f64 {
        self.hermite_polynomial(y_over_le) * (-0.5 * y_over_le * y_over_le).exp()
    }

    /// `∫ψₘ² dŷ = 2ᵐ·m!·√π`, the normalisation a projection divides by.
    ///
    /// Written as the product `∏ₖ₌₁..ₘ 2k`, which is that value, so the two
    /// orders share one definition instead of carrying a tabulated constant
    /// each.
    fn hermite_norm(self) -> f64 {
        let weight: f64 = (1..=self.hermite_order()).map(|k| 2.0 * k as f64).product();
        weight * std::f64::consts::PI.sqrt()
    }

    /// How much finer this mode's meridional structure is than the
    /// deformation radius, as the factor `2m + 1` multiplying `(Δy/Le)²` in
    /// its truncation bound.
    ///
    /// `ψₘ` oscillates on the scale `Le/√(2m+1)` — its classical turning point
    /// is at `ŷ = √(2m+1)` — so a second-order scheme's truncation error on it
    /// grows with the square of that local wavenumber. `ψ₀` gives one, which
    /// is why the Kelvin bound is `(Δy/Le)²` with no factor; `ψ₂` gives five,
    /// and taking one there would be optimistic rather than conservative.
    fn truncation_richness(self) -> f64 {
        (2 * self.hermite_order() + 1) as f64
    }
}

/// The equatorial waveguide of one run, as far as a meridional projection
/// needs to know it: the scale the Hermite functions are stretched by, and the
/// rows a column of `q` is sampled on.
struct Waveguide {
    /// Equatorial deformation radius `Le = √(c/β)`, in metres (`CONTEXT.md`).
    le_m: f64,
    /// Cell height, in metres — the quadrature weight of one row.
    dy_m: f64,
    /// Meridional positions of the cell-centre rows, in metres north of the
    /// equator.
    row_y_m: Vec<f64>,
}

impl Waveguide {
    /// The waveguide of `basin`, sampled on its cell-centre rows.
    fn new(basin: Basin, params: PhysicalParams) -> Self {
        Self {
            le_m: (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt(),
            dy_m: basin.spacing().dy_m(),
            row_y_m: (0..basin.grid().ny())
                .map(|j| basin.y_of_row_m(H_STAGGERING, j))
                .collect(),
        }
    }

    /// The coefficient of `ψₘ` in a meridional column of `q`, by discrete
    /// quadrature of `∫q·ψₘ dŷ / ∫ψₘ² dŷ`.
    ///
    /// The midpoint rule on the cell-centre rows, which is second order and
    /// symmetric about the equator — the same order as the scheme whose output
    /// it reads.
    fn coefficient(&self, column: &[f64], mode: MeridionalMode) -> f64 {
        let integral: f64 = column
            .iter()
            .zip(&self.row_y_m)
            .map(|(value, y_m)| value * mode.hermite_function(y_m / self.le_m))
            .sum();
        integral * (self.dy_m / self.le_m) / mode.hermite_norm()
    }

    /// The second-order meridional truncation bound for `mode`, as a fraction
    /// of a wave speed: `(2m + 1)·(Δy/Le)²`, with the remaining `O(1)`
    /// constant taken as one.
    fn truncation_bound(&self, mode: MeridionalMode) -> f64 {
        mode.truncation_richness() * (self.dy_m / self.le_m).powi(2)
    }

    /// The fraction by which a Gaussian packet of e-folding half-width `σ`
    /// lags the long-wave Rossby speed: `⟨κ²⟩ = Le²/(2σ²)`, from the
    /// leading-order group speed `c_g = −(c/3)(1 − κ²)`.
    fn packet_width_bias(&self) -> f64 {
        let width_in_radii = PULSE_ZONAL_WIDTH_M / self.le_m;
        0.5 / (width_in_radii * width_in_radii)
    }
}

/// Index of the cell-centre column whose position is closest to `x_m`.
fn column_nearest(x_m: f64, cell_m: f64) -> usize {
    (x_m / cell_m - 0.5).round().max(0.0) as usize
}

/// The time, in seconds, at which `series` peaks, refined below the sampling
/// interval by fitting a parabola through the largest sample and its two
/// neighbours.
///
/// # Panics
/// If the peak sits on either end of the record, which would mean the run is
/// too short or a station is misplaced rather than that the wave is slow.
fn peak_time_s(series: &[f64], dt_s: f64) -> f64 {
    let peak = peak_index(series);
    assert!(
        peak > 0 && peak + 1 < series.len(),
        "the packet peaks at sample {peak} of {}, on the edge of the record: \
         the run is too short, or the station is in the wrong place",
        series.len()
    );
    let (before, at, after) = (series[peak - 1], series[peak], series[peak + 1]);
    let curvature = before - 2.0 * at + after;
    assert!(
        curvature < 0.0,
        "the three samples around the largest one are not concave, so the record has no \
         resolved peak to time"
    );
    (peak as f64 + 0.5 * (before - after) / curvature) * dt_s
}

/// Index of the largest value in `values`.
fn peak_index(values: &[f64]) -> usize {
    values
        .iter()
        .enumerate()
        .max_by(|(_, left), (_, right)| {
            left.partial_cmp(right)
                .expect("an undamped linear run produces no NaN")
        })
        .map(|(index, _)| index)
        .expect("a recorded series is never empty")
}

/// What one zonal station saw over a run.
///
/// Both modes are projected at every station rather than only the one that
/// station exists to time: it costs one more dot product per step, and it
/// means the tests that ask what the reflected signal is *made of* read the
/// same record as the tests that time it.
struct StationRecord {
    /// Zonal position of the sampled column, in metres east of the western
    /// wall — the column's own centre, not the position asked for.
    x_m: f64,
    /// The `ψ₀` coefficient of `q` there, one value per recorded step.
    kelvin: Vec<f64>,
    /// The `ψ₂` coefficient of `q` there, one value per recorded step.
    gravest_rossby: Vec<f64>,
    /// `q` on this station's meridional column, one column of `ny` values per
    /// recorded step.
    ///
    /// Every step is kept because the step a wave peaks at is only known once
    /// the record is complete, and the meridional structure has to be read at
    /// that step rather than at an arbitrary one.
    columns: Vec<Vec<f64>>,
}

/// Everything one integration of the basin recorded.
struct ReflectionRun {
    /// Length of one step, in seconds.
    dt_s: f64,
    /// The waveguide the projections were taken against.
    waveguide: Waveguide,
    /// The four stations, indexed by [`KELVIN_WEST`], [`KELVIN_EAST`],
    /// [`ROSSBY_EAST`] and [`ROSSBY_WEST`].
    stations: [StationRecord; 4],
}

impl ReflectionRun {
    /// Zonal speed of the reflected Kelvin wave, in m/s, as a time of flight
    /// between the two reflection stations. Positive is eastward.
    fn reflected_kelvin_speed_m_per_s(&self) -> f64 {
        let (west, east) = (&self.stations[KELVIN_WEST], &self.stations[KELVIN_EAST]);
        (east.x_m - west.x_m)
            / (peak_time_s(&east.kelvin, self.dt_s) - peak_time_s(&west.kelvin, self.dt_s))
    }

    /// Zonal speed of the incident Rossby packet, in m/s. Negative is
    /// westward.
    fn incident_rossby_speed_m_per_s(&self) -> f64 {
        let (east, west) = (&self.stations[ROSSBY_EAST], &self.stations[ROSSBY_WEST]);
        (west.x_m - east.x_m)
            / (peak_time_s(&west.gravest_rossby, self.dt_s)
                - peak_time_s(&east.gravest_rossby, self.dt_s))
    }

    /// The station whose record the reflected wave's *structure* is read from,
    /// and the step at which it peaks there.
    ///
    /// The eastern reflection station, and not the western one, because a
    /// structure is a statement about the whole meridional column rather than
    /// about one projection of it. The reflected Kelvin wave passes the
    /// western station while the tail of the incident packet is still
    /// arriving — the two are orthogonal in projection, so a speed measured
    /// there is clean, but the column itself is their sum, and `ψ₂` stands
    /// more than twice as tall as `ψ₀` at its off-equatorial lobes. By the
    /// time the Kelvin wave reaches the eastern station 3500 km further on,
    /// the incident packet is `exp(−18)` of its amplitude there and the column
    /// is the reflected wave alone. Nothing is lost by moving: the Kelvin wave
    /// is non-dispersive and the run is undamped, so it arrives with the shape
    /// and the amplitude it left with.
    fn reflected_signal(&self) -> (&StationRecord, usize) {
        let station = &self.stations[KELVIN_EAST];
        (station, peak_index(&station.kelvin))
    }
}

/// The initial condition: a gravest-mode (`n = 1`) equatorial Rossby packet,
/// centred at [`PULSE_CENTER_X_M`] and travelling west.
///
/// The three fields are the long-wave solution written out in the module
/// header, sampled at each variable's own C-grid position. `v` is left at
/// rest: it is `O(Le/σ)` smaller than `u` in the long-wave limit, which is the
/// same approximation the `c/3` prediction itself rests on.
///
/// `u` is not zero at the western wall — the packet's tail reaches it with
/// `exp(−8) ≈ 3×10⁻⁴` of its amplitude — and the solver zeroes it on the way
/// into the first step, exactly as [`engine::NoNormalFlow`] promises. That is
/// the analytic solution being brought onto the closed basin's boundary
/// condition, and it is three orders of magnitude below the signal being
/// measured.
fn incident_rossby_packet(
    basin: Basin,
    params: PhysicalParams,
    waveguide: &Waveguide,
) -> OceanState {
    let grid = basin.grid();
    let mut state = OceanState::at_rest(grid);
    let le_m = waveguide.le_m;
    let mean_depth_m = params.mean_thermocline_depth_m();
    let wave_speed_m_per_s = params.kelvin_wave_speed_m_per_s();
    let amplitude = PULSE_EQUATORIAL_AMPLITUDE_M / mean_depth_m;

    let envelope = |x_m: f64| {
        let offset = (x_m - PULSE_CENTER_X_M) / PULSE_ZONAL_WIDTH_M;
        (-0.5 * offset * offset).exp()
    };

    let (h_nx, h_ny) = grid.field_shape(H_STAGGERING);
    for j in 0..h_ny {
        let y_hat = basin.y_of_row_m(H_STAGGERING, j) / le_m;
        let trapping = (-0.5 * y_hat * y_hat).exp();
        for i in 0..h_nx {
            let x_m = basin.x_of_column_m(H_STAGGERING, i);
            *state
                .h_mut()
                .get_mut(i, j)
                .expect("the loop bounds are the field's own shape") =
                mean_depth_m * amplitude * envelope(x_m) * (2.0 * y_hat * y_hat + 1.0) * trapping;
        }
    }

    let (u_nx, u_ny) = grid.field_shape(U_STAGGERING);
    for j in 0..u_ny {
        let y_hat = basin.y_of_row_m(U_STAGGERING, j) / le_m;
        let trapping = (-0.5 * y_hat * y_hat).exp();
        for i in 0..u_nx {
            let x_m = basin.x_of_column_m(U_STAGGERING, i);
            *state
                .u_mut()
                .get_mut(i, j)
                .expect("the loop bounds are the field's own shape") = wave_speed_m_per_s
                * amplitude
                * envelope(x_m)
                * (2.0 * y_hat * y_hat - 3.0)
                * trapping;
        }
    }

    state
}

/// Integrate the basin at `COARSE_CELL_SIZE_M / refinement` and record what
/// the four stations saw.
///
/// The physical configuration — basin, packet, stations, run length — is the
/// same at every refinement; only the discretisation changes, which is what
/// makes two runs comparable.
fn run_reflection(refinement: usize) -> ReflectionRun {
    let cell_m = COARSE_CELL_SIZE_M / refinement as f64;
    let nx = (BASIN_LX_M / cell_m).round() as usize;
    let ny = (BASIN_LY_M / cell_m).round() as usize;
    let grid = Grid::new(nx, ny).expect("the basin has cells on both axes");
    let spacing = Spacing::new(cell_m, cell_m).expect("the cell size is finite and positive");
    let params = pacific_params();
    let basin = Basin::centered_on_equator(grid, spacing);
    let plane = BetaPlane::centered_on_equator(params, spacing, grid);
    let waveguide = Waveguide::new(basin, params);

    let wave_speed = WaveSpeed::new(params.kelvin_wave_speed_m_per_s())
        .expect("the Kelvin wave speed of a physical ocean is positive");
    let dt_s = CFL_FRACTION * max_stable_dt(spacing, wave_speed);
    let steps = (RUN_DURATION_S / dt_s).ceil() as usize;

    let mut solver = Solver::new(grid, spacing, params, plane, dt_s)
        .expect("a timestep inside both the CFL and the rotation bound is accepted");
    let mut state = incident_rossby_packet(basin, params, &waveguide);
    let calm = WindStressField::calm(grid);
    let operators = CGridOperators::new(grid, spacing);
    let mut u_at_centers = grid.allocate(H_STAGGERING, 0.0);

    let columns = [
        column_nearest(KELVIN_STATION_WEST_X_M, cell_m),
        column_nearest(KELVIN_STATION_EAST_X_M, cell_m),
        column_nearest(ROSSBY_STATION_EAST_X_M, cell_m),
        column_nearest(ROSSBY_STATION_WEST_X_M, cell_m),
    ];
    let mut stations = columns.map(|i| StationRecord {
        x_m: basin.x_of_column_m(H_STAGGERING, i),
        kelvin: Vec::with_capacity(steps + 1),
        gravest_rossby: Vec::with_capacity(steps + 1),
        columns: Vec::with_capacity(steps + 1),
    });
    let mut column = vec![0.0; ny];

    for step in 0..=steps {
        operators.face_to_center_x(state.u(), &mut u_at_centers);
        for (record, i) in stations.iter_mut().zip(columns) {
            // `q = u/c + h/H` down one meridional column. Both parts are read
            // at the cell centres, which is what makes the projection below a
            // projection of one field.
            for (j, value) in column.iter_mut().enumerate() {
                let u_m_per_s = *u_at_centers
                    .get(i, j)
                    .expect("a station column is inside the basin");
                let h_m = *state
                    .h()
                    .get(i, j)
                    .expect("a station column is inside the basin");
                *value = u_m_per_s / params.kelvin_wave_speed_m_per_s()
                    + h_m / params.mean_thermocline_depth_m();
            }
            record
                .kelvin
                .push(waveguide.coefficient(&column, MeridionalMode::Kelvin));
            record
                .gravest_rossby
                .push(waveguide.coefficient(&column, MeridionalMode::GravestRossby));
            record.columns.push(column.clone());
        }

        if step < steps {
            solver.step(&mut state, step as f64 * dt_s, |_| &calm);
        }
    }

    ReflectionRun {
        dt_s,
        waveguide,
        stations,
    }
}

/// The coarse run, integrated once and shared by every test that reads it.
fn coarse_run() -> &'static ReflectionRun {
    static RUN: OnceLock<ReflectionRun> = OnceLock::new();
    RUN.get_or_init(|| run_reflection(1))
}

#[test]
fn the_western_boundary_returns_the_incident_rossby_packet_as_an_eastward_kelvin_wave() {
    // The acceptance criterion of T-04.3, at the boundary the scientific model
    // doc places this process at: the reflected signal's propagation speed is
    // the Kelvin speed `c = √(g'·H)`, eastward.
    let run = coarse_run();
    let params = pacific_params();
    let expected_m_per_s = params.kelvin_wave_speed_m_per_s();

    // The wall is what produced it. The packet is built on `ψ₂`, orthogonal to
    // the Kelvin wave's `ψ₀`, so the run starts with no Kelvin content to
    // propagate; without this the speed below would be reproduced just as well
    // by a bug that seeded a Kelvin wave at `t = 0`.
    let west = &run.stations[KELVIN_WEST];
    let seeded = west.kelvin[0].abs();
    let reflected_peak = west.kelvin[peak_index(&west.kelvin)];
    assert!(
        seeded <= SEEDED_KELVIN_CEILING * reflected_peak,
        "the initial condition already carries a ψ₀ coefficient of {seeded} at the western \
         station, against a reflected peak of {reflected_peak}: the Kelvin wave measured below \
         is not the wall's doing"
    );

    let measured_m_per_s = run.reflected_kelvin_speed_m_per_s();
    assert!(
        measured_m_per_s > 0.0,
        "the reflected signal travels at {measured_m_per_s} m/s, which is westward: a western \
         boundary can only radiate the eastward branch"
    );

    let tolerance = run.waveguide.truncation_bound(MeridionalMode::Kelvin);
    let error = (measured_m_per_s - expected_m_per_s).abs() / expected_m_per_s;
    assert!(
        error <= tolerance,
        "the reflected Kelvin wave travels at {measured_m_per_s} m/s against the analytic \
         c = {expected_m_per_s} m/s: a relative error of {error}, past the {tolerance} the \
         meridional resolution allows"
    );
}

#[test]
fn reflected_kelvin_speed_converges_at_second_order() {
    // The point check above is bounded by `(Δy/Le)²` with its constant taken
    // as one, which is generous. This is the assertion that the error really
    // is the second-order truncation of the scheme and not a fixed offset that
    // happens to fit under it: halving the cell size must shrink it by the
    // order the scheme claims (CODING_STANDARDS.md § Convergence over point
    // checks).
    let params = pacific_params();
    let expected_m_per_s = params.kelvin_wave_speed_m_per_s();
    let relative_error = |run: &ReflectionRun| {
        (run.reflected_kelvin_speed_m_per_s() - expected_m_per_s).abs() / expected_m_per_s
    };

    let coarse = relative_error(coarse_run());
    let fine = relative_error(&run_reflection(2));
    let order = (coarse / fine).log2();

    assert!(
        order >= MIN_CONVERGENCE_ORDER,
        "halving the cell size took the reflected Kelvin speed's error from {coarse} to {fine}, \
         a convergence order of {order}: short of the {MIN_CONVERGENCE_ORDER} a second-order \
         scheme owes"
    );
}

#[test]
fn the_incident_rossby_packet_travels_west_at_a_third_of_the_kelvin_speed() {
    // The other half of the reflection: the wave that arrives is the gravest
    // meridional Rossby mode, whose long-wave speed is `−c/3` (`CONTEXT.md`,
    // *Rossby wave*). Measured on the same run, through the `ψ₂` channel the
    // Kelvin wave is orthogonal to.
    let run = coarse_run();
    let wave_speed_m_per_s = pacific_params().kelvin_wave_speed_m_per_s();
    let expected_m_per_s = wave_speed_m_per_s / 3.0;
    let measured_m_per_s = run.incident_rossby_speed_m_per_s();

    assert!(
        measured_m_per_s < 0.0,
        "the incident packet travels at {measured_m_per_s} m/s, which is eastward: the gravest \
         Rossby mode goes west"
    );
    let measured_speed_m_per_s = -measured_m_per_s;

    // Two terms, both derived in the module header, and only one of them
    // two-sided. The resolution's truncation could go either way; the lag a
    // packet of finite zonal width carries has a sign — a finite packet is
    // slower than the long-wave limit, never faster — so it widens the slow
    // side alone.
    let truncation = run
        .waveguide
        .truncation_bound(MeridionalMode::GravestRossby);
    let slowest_m_per_s = expected_m_per_s * (1.0 - truncation - run.waveguide.packet_width_bias());
    let fastest_m_per_s = expected_m_per_s * (1.0 + truncation);
    assert!(
        measured_speed_m_per_s >= slowest_m_per_s && measured_speed_m_per_s <= fastest_m_per_s,
        "the incident packet travels at {measured_speed_m_per_s} m/s against the analytic \
         c/3 = {expected_m_per_s} m/s, outside the [{slowest_m_per_s}, {fastest_m_per_s}] m/s \
         the resolution and the packet's width allow"
    );

    // That band is wide enough to be worth pinning down further. The analytic
    // speeds either side of it are the `n = 2` Rossby mode at `c/5` and the
    // Kelvin branch at `c`, and the measured speed must be nearer `c/3` than
    // to either — which is what identifies it as the *gravest* mode rather
    // than merely a westward signal of about the right speed.
    let miss = (measured_speed_m_per_s - expected_m_per_s).abs();
    for (name, other_m_per_s) in [
        ("the n = 2 mode's c/5", wave_speed_m_per_s / 5.0),
        ("the Kelvin speed c", wave_speed_m_per_s),
    ] {
        assert!(
            miss < (measured_speed_m_per_s - other_m_per_s).abs(),
            "the incident packet's {measured_speed_m_per_s} m/s is nearer {name} \
             ({other_m_per_s} m/s) than the gravest mode's c/3 = {expected_m_per_s} m/s"
        );
    }
}

#[test]
fn the_reflected_signal_is_an_equatorially_trapped_kelvin_wave_carrying_part_of_the_incident_energy(
) {
    // Speed alone does not make a signal a Kelvin wave. Long-wave theory also
    // says what it looks like — single-lobed on the equator, on `ψ₀` and not
    // on `ψ₂` — and roughly how much of the incident wave comes back.
    let run = coarse_run();
    let (station, peak) = run.reflected_signal();

    let kelvin = station.kelvin[peak];
    let rossby = station.gravest_rossby[peak];
    assert!(
        kelvin.abs() > KELVIN_DOMINANCE * rossby.abs(),
        "at its peak the reflected signal carries a ψ₀ coefficient of {kelvin} against a ψ₂ \
         coefficient of {rossby}: theory says the reflected long wave is purely ψ₀"
    );

    // `ψ₀` peaks on the equator; the incident `ψ₂` has a *minimum* there
    // between two off-equatorial lobes, so this is the check that tells the
    // two structures apart. The basin has an even number of rows, so the
    // equator falls between the two innermost of them.
    let profile = &station.columns[peak];
    let equator_row = profile.len() / 2;
    let peak_row = peak_index(profile);
    assert!(
        peak_row == equator_row || peak_row + 1 == equator_row,
        "the reflected signal peaks at row {peak_row} of {}, whose centre is {} m from the \
         equator: an equatorially trapped Kelvin wave peaks on the equator",
        profile.len(),
        run.waveguide.row_y_m[peak_row]
    );

    // The long-wave prediction: a `ψ₀` coefficient of `2A` in `q`, from the
    // zero-net-mass-flux condition at the wall.
    let incident_amplitude = PULSE_EQUATORIAL_AMPLITUDE_M / PACIFIC_MEAN_DEPTH_M;
    let predicted = 2.0 * incident_amplitude;
    let (floor, ceiling) = KELVIN_AMPLITUDE_BAND;
    assert!(
        kelvin > floor * predicted && kelvin < ceiling * predicted,
        "the reflection returned a ψ₀ coefficient of {kelvin} against the long-wave prediction \
         of {predicted}: outside the [{floor}, {ceiling}] band the short Rossby waves and the \
         packet's dispersion justify"
    );
}
