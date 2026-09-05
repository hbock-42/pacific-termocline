//! Acceptance tests for T-07.5 — conservation in the undamped, unforced limit.
//!
//! T-02.5 checked, over four basin crossings, that a run with `r = 0` and
//! `τ = 0` keeps its energy. This file is that check made permanent and made
//! *long*: the same limit, run for eight times as far, against a bound written
//! down from the scheme rather than read off a run.
//!
//! ```text
//! ∂u/∂t = +f·v − g'·∂h/∂x
//! ∂v/∂t = −f·u − g'·∂h/∂y
//! ∂h/∂t =       −H·(∂u/∂x + ∂v/∂y)
//! ```
//!
//! The conserved quantity is the discrete wave energy
//! `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)`, in which the pressure-gradient and
//! continuity pair is exactly skew (see below), so `E` is what this
//! discretisation conserves rather than merely what the continuous equations
//! do.
//!
//! # Where the bound comes from
//!
//! With `r = 0` and `τ = 0` the semi-discrete system is `ẋ = (W + C)·x`, `W`
//! the pressure-gradient/continuity pair and `C` the beta-plane rotation.
//! Exactly two things stop `E` from being constant, and both are properties of
//! the scheme, derived here and never measured:
//!
//! ## 1. The Coriolis skewness defect — `O(Δy²)`, and bounded in time
//!
//! `W` is **exactly** skew in `E`. Summation by parts on the C-grid leaves a
//! boundary term `h·u` at the walls, and the operators of T-01.1 leave every
//! wall face at zero, so a basin that starts with quiescent walls and is never
//! forced keeps them quiescent and the boundary term vanishes for all time.
//! `W` therefore contributes nothing to `dE/dt`, at any `Δx` and any `Δy`.
//!
//! `C` is skew only to `O(Δy²)`. The u equation evaluates `f = β·y` on the
//! cell-center rows and the v equation on the north/south-face rows, half a
//! cell apart, so the two four-point averages fail to cancel by a residual
//! weight `±β·Δy/2` and leave (T-02.5, `time_stepping.rs`)
//!
//! ```text
//! dE/dt = −H·(β·Δy²/4)·Σ_u u·∂v/∂y + O(Δy⁴).
//! ```
//!
//! Turning that into a bound on the *relative* drift takes two inequalities.
//! The solution is trapped in the equatorial waveguide, whose meridional scale
//! is the deformation radius `Le = √(c/β)` (`CONTEXT.md`), so
//! `|∂v/∂y| ≤ |v|/Le`; and `E ≥ (H/2)·Σ(u² + v²) ≥ H·Σ|u|·|v|` by the
//! arithmetic–geometric mean inequality. Together,
//!
//! ```text
//! |dE/dt| / E  ≤  β·Δy²/(4·Le)  =  (1/4)·√(β·c)·(Δy/Le)²,
//! ```
//!
//! using `β·Le = √(β·c)`. That is a *rate*, and the run is long, so the rate
//! alone would allow an enormous drift. It does not, because the defect is
//! **oscillatory rather than one-signed**: `u` and `∂v/∂y` are two components
//! of the same wave field, and their product changes sign every quarter period
//! of it. Integrating a rate `a·cos(ω·t)` gives `(a/ω)·sin(ω·t)`, an excursion
//! of `a/ω` and not one that grows — so at the waveguide's own frequency
//! `ω = c/Le = √(β·c)` the `√(β·c)` cancels and what is left is
//!
//! ```text
//! |ΔE/E|  ≤  (1/4)·(Δy/Le)²
//! ```
//!
//! that does **not** grow with the length of the run. This is the term the
//! long run exists to expose: a defect that accumulated instead of oscillating
//! would multiply this by the number of wave periods in the run — some 150 of
//! them here, `2π/√(β·c) ≈ 7.9×10⁵ s` apiece — and so would break the bound by
//! more than two orders of magnitude.
//!
//! ## 2. RK4's numerical dissipation — `O(dt⁶)`, and linear in the step count
//!
//! This is the "numerical diffusion" the ticket names. On the purely imaginary
//! spectrum of a skew system, one RK4 step multiplies a mode of frequency `ω`
//! by `R(i·θ)` with `θ = ω·dt` and
//! `R(z) = 1 + z + z²/2 + z³/6 + z⁴/24`, whence
//!
//! ```text
//! |R(i·θ)|² = 1 − θ⁶/72 + θ⁸/576.
//! ```
//!
//! Two things follow, and they are of very different strengths.
//!
//! The **sign** is rigorous. `|R(i·θ)| ≤ 1` for every `θ` inside the CFL bound
//! the solver enforces, so RK4 can only ever *remove* energy from this system,
//! whatever frequencies it happens to contain. The energy therefore cannot be
//! pushed *up* past the skewness excursion of §1 by anything in the time
//! discretisation, and
//! [`an_undamped_unforced_run_never_gains_energy_past_the_skewness_bound`]
//! asserts exactly that, against the `(1/4)·(Δy/Le)²` of §1 alone.
//!
//! The **size** of the loss is an estimate, and is labelled as one. Over `N`
//! steps it is `N·θ⁶/72` with `θ = ω·dt`, and the `ω` that belongs there is
//! the energy-weighted frequency of the modes the run actually excites.
//! [`rk4_dissipation_bound`] uses `√(β·c)`, the equatorial inertia-gravity
//! frequency of the trapped motion the beta-plane makes out of the gravest
//! zonal mode. That is the scale of the motion the run is *about*; what it is
//! not is a supremum over the discrete spectrum, whose grid-scale end sits at
//! `θ = O(1)` where this expression is meaningless. The linearity of the v1
//! core (CODING_STANDARDS.md § Scope guards) is what makes the substitution
//! reasonable — there is no cascade to carry energy to the grid scale, and the
//! only coupling between modes is `C`, whose meridional scale is `Le` — but
//! "reasonable" is the honest word for it, not "proved".
//!
//! ## What each assertion is worth
//!
//! [`derived_drift_bound`] is the sum of the two terms, and
//! [`energy_drift_over_a_long_undamped_unforced_run_stays_within_the_derived_bound`]
//! holds the run to it. No coefficient in either term was chosen by looking at
//! a run — the check runs the other way, and the margin that remains is the
//! `O(1)` shape factors the inequalities above gave away.
//!
//! But a sum whose second term rests on an estimate is not the strongest thing
//! here, and the file does not lean on it alone:
//!
//! - [`an_undamped_unforced_run_never_gains_energy_past_the_skewness_bound`]
//!   uses only §1 and the sign of §2, both rigorous.
//! - [`the_long_run_energy_drift_falls_at_the_schemes_second_order_under_refinement`]
//!   holds the *measured* drift to the `Δy²` rate of §1 across three
//!   resolutions (CODING_STANDARDS.md § Tests), which no bound can fake:
//!   §2 falls as `dt⁵` at a fixed run length, faster still.
//! - [`the_conservation_run_exchanges_energy_between_potential_and_kinetic_form`]
//!   rules out the vacuous pass of a run in which nothing moves.
//!
//! Refining in `y` alone is the right study for all of them, because `W` is
//! exactly skew at every `Δx`: the whole energy error is meridional, and a
//! zonal refinement would only make the runs more expensive.

use std::sync::OnceLock;

use engine::{
    max_stable_dt, BetaPlane, Field2D, Grid, OceanState, PhysicalParams, Solver, Spacing,
    Staggering, WaveSpeed, WindStressField, H_STAGGERING,
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
/// Rayleigh damping `r` of the limit under test, in s⁻¹. The ticket's `r = 0`:
/// with any damping at all the energy budget is T-02.4's decay law, not a
/// conservation statement.
const UNDAMPED_PER_S: f64 = 0.0;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Deliberately different from
/// [`BASIN_LY_M`] so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: an equatorial channel
/// reaching ±500 km, or about 1.4 equatorial deformation radii either side of
/// the equator.
///
/// The same channel `time_stepping.rs` runs on, and for the same reason: at
/// the ±2500 km of `CONTEXT.md` the rotation bound of ADR-0007, not the CFL
/// bound, would set the timestep, which is a different experiment.
const BASIN_LY_M: f64 = 1.0e6;

/// Cells across the basin in `x`, at every resolution of the study.
///
/// Fixed rather than refined: the pressure-gradient/continuity pair is exactly
/// skew at every `Δx`, so `Δx` contributes nothing to the energy error and
/// refining it would only lengthen the runs. See the module comment.
const BASIN_CELLS_X: usize = 16;
/// Cells across the basin in `y`: the same channel resolved at three
/// resolutions, each a halving of `Δy` — the spacing the drift bound depends
/// on.
const MERIDIONAL_RESOLUTIONS: [usize; 3] = [16, 32, 64];

/// Amplitude of the test thermocline depth anomaly, in metres. A 20 m
/// departure is the scale of an observed equatorial Pacific anomaly.
const H_AMPLITUDE_M: f64 = 20.0;

/// Length of the conservation run, in crossings of the basin at the Kelvin
/// wave speed `c = √(g'·H)`.
///
/// Eight times the four crossings of the T-02.5 check this ticket formalizes:
/// about 3.7 years of simulated time, between 4.5×10³ and 1.8×10⁴ steps
/// depending on resolution, and some 150 periods of the equatorial
/// inertia-gravity motion whose skewness defect the bound is written around.
/// That last number is the one that matters: the bound of the module comment
/// is one period's excursion, so a defect that accumulated rather than
/// oscillated would overshoot it by two orders of magnitude here.
const LONG_RUN_CROSSINGS: f64 = 32.0;

/// Coefficient of the Coriolis skewness excursion `(Δy/Le)²`.
///
/// The `1/4` of `dE/dt = −H·(β·Δy²/4)·Σ_u u·∂v/∂y`, carried through the two
/// inequalities of the module comment — `|∂v/∂y| ≤ |v|/Le` and
/// `E ≥ H·Σ|u|·|v|` — and divided by the frequency `√(β·c)` the defect
/// oscillates at. Every factor is the discretisation's; none is fitted.
const CORIOLIS_SKEWNESS_COEFFICIENT: f64 = 0.25;

/// Denominator of RK4's per-step energy loss `θ⁶/72` on an imaginary
/// eigenvalue.
///
/// From `|R(i·θ)|² = 1 − θ⁶/72 + θ⁸/576`, the modulus of the RK4 amplification
/// polynomial quoted in `termocline-numerics`' CFL derivation (Hairer &
/// Wanner, *Solving Ordinary Differential Equations I*, § II.2). The `θ⁸` term
/// is positive, so dropping it keeps this a bound.
const RK4_ENERGY_LOSS_DENOMINATOR: f64 = 72.0;

/// How far below second order the measured convergence of the energy drift may
/// sit.
///
/// The same 0.2 T-02.5 arrived at for the same quantity, and for the same
/// reason: the slack is the sub-leading `O(Δy⁴)` correction, which is what
/// keeps a measurement at finite resolution off the asymptote. Writing the
/// drift as `A·Δy²·(1 + c·Δy²)`, a halving of `Δy` moves the measured order by
/// about `(3/4)·ε/ln 2` for a correction worth a fraction `ε` of the leading
/// term, so 0.2 admits a sub-leading term up to about 18% of the leading one.
/// It is reproduced here rather than shared because the two files pin
/// different runs: `time_stepping.rs` refines both axes over four crossings,
/// this one refines `Δy` alone over thirty-two, and a single constant would
/// tie the two studies' slack together for no reason.
///
/// The check is one-sided: falling faster than second order is not a defect,
/// falling slower is — and first-order or resolution-independent behaviour,
/// which is what a mis-wired step shows, is far outside it.
const DRIFT_ORDER_TOLERANCE: f64 = 0.2;

/// Fraction of the initial energy that must reach kinetic form at some point
/// in a run for the run to count as having exercised the equations.
///
/// A guard against a vacuous pass, not a physical prediction: the initial
/// state is pure potential energy and at rest, so a run in which nothing ever
/// moves would conserve energy perfectly and prove nothing. The gravest
/// standing mode puts *all* of its energy into motion twice per period in the
/// non-rotating limit, so any threshold below one is met by a working step;
/// a half is well clear of both a working run and a dead one.
const KINETIC_EXCHANGE_FLOOR: f64 = 0.5;

/// The equatorial-Pacific parameter set, undamped.
fn undamped_pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        UNDAMPED_PER_S,
        BETA_PER_M_PER_S,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// Kelvin wave speed `c = √(g'·H)`, in m/s, written out from the definition in
/// `CONTEXT.md` rather than asked of the code under test.
fn kelvin_wave_speed_m_per_s() -> f64 {
    (PACIFIC_REDUCED_GRAVITY_M_PER_S2 * PACIFIC_MEAN_DEPTH_M).sqrt()
}

/// Equatorial deformation radius `Le = √(c/β)`, in metres — the meridional
/// scale of the waveguide (`CONTEXT.md`), and the length the drift bound
/// measures `Δy` against.
fn equatorial_deformation_radius_m() -> f64 {
    (kelvin_wave_speed_m_per_s() / BETA_PER_M_PER_S).sqrt()
}

/// The equatorial inertia-gravity frequency `√(β·c) = c/Le`, in s⁻¹: the
/// fastest frequency a smooth, equatorially trapped solution of the linear
/// core carries, and the `ω` of both terms of the bound.
fn equatorial_wave_frequency_per_s() -> f64 {
    (BETA_PER_M_PER_S * kelvin_wave_speed_m_per_s()).sqrt()
}

/// A basin of [`BASIN_CELLS_X`] by `cells_y` cells spanning [`BASIN_LX_M`] by
/// [`BASIN_LY_M`]. `dx ≠ dy`, so an x/y swap cannot pass.
fn basin(cells_y: usize) -> (Grid, Spacing) {
    let grid = Grid::new(BASIN_CELLS_X, cells_y).expect("extents are non-zero");
    let spacing = Spacing::new(
        BASIN_LX_M / BASIN_CELLS_X as f64,
        BASIN_LY_M / cells_y as f64,
    )
    .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

/// Position of a field point, in metres from the basin's southwest corner.
fn position_m(spacing: Spacing, staggering: Staggering, i: usize, j: usize) -> (f64, f64) {
    let (offset_x, offset_y) = staggering.offset_in_cells();
    (
        (i as f64 + offset_x) * spacing.dx_m(),
        (j as f64 + offset_y) * spacing.dy_m(),
    )
}

/// A field at `staggering` holding `g` sampled at each of its points.
fn sample(
    grid: Grid,
    spacing: Spacing,
    staggering: Staggering,
    g: impl Fn(f64, f64) -> f64,
) -> Field2D<f64> {
    let mut field = grid.allocate(staggering, 0.0_f64);
    for j in 0..field.ny() {
        for i in 0..field.nx() {
            let (x_m, y_m) = position_m(spacing, staggering, i, j);
            *field.get_mut(i, j).expect("in-bounds point") = g(x_m, y_m);
        }
    }
    field
}

/// The gravest zonal standing mode of the closed basin: `h = A·cos(π·x/Lx)`,
/// at rest, uniform in `y`.
///
/// It is an exact eigenvector of the discrete wave operator — the C-grid
/// difference of `cos(k·x)` at cell centers is `sin(k·x)` at faces, which
/// vanishes on both walls — so its velocities start and stay exactly zero on
/// the four walls, which is the condition under which the discrete energy's
/// boundary term vanishes and `W` is exactly skew.
fn gravest_zonal_mode(grid: Grid, spacing: Spacing) -> OceanState {
    let wavenumber_per_m = std::f64::consts::PI / BASIN_LX_M;
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, _y_m| {
        H_AMPLITUDE_M * (wavenumber_per_m * x_m).cos()
    });
    state
}

/// Sum of the squares of a field's points.
fn sum_of_squares(field: &Field2D<f64>) -> f64 {
    field.as_slice().iter().map(|value| value * value).sum()
}

/// The kinetic half of the discrete energy, `(H/2)·Σ(u² + v²)`, in m³/s².
fn kinetic_energy(state: &OceanState, params: PhysicalParams) -> f64 {
    0.5 * params.mean_thermocline_depth_m()
        * (sum_of_squares(state.u()) + sum_of_squares(state.v()))
}

/// The discrete energy `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)`, summed over grid
/// points, in m³/s² — energy per unit reference density and per unit cell
/// area, which is constant across a run and so drops out of every ratio taken
/// here.
///
/// These are the weights that make the pressure-gradient and continuity pair
/// skew under summation by parts on the C-grid; see the module comment.
fn wave_energy(state: &OceanState, params: PhysicalParams) -> f64 {
    0.5 * params.reduced_gravity_m_per_s2() * sum_of_squares(state.h())
        + kinetic_energy(state, params)
}

/// The largest relative energy drift the Coriolis pair's skewness defect can
/// produce on a grid of meridional spacing `dy_m`.
///
/// `(1/4)·(Δy/Le)²`, derived in the module comment. Independent of the length
/// of the run, because the defect oscillates rather than accumulates.
fn coriolis_skewness_bound(dy_m: f64) -> f64 {
    CORIOLIS_SKEWNESS_COEFFICIENT * (dy_m / equatorial_deformation_radius_m()).powi(2)
}

/// The largest relative energy loss `steps` RK4 steps of `dt_s` seconds can
/// produce on this solution.
///
/// `N·θ⁶/72` with `θ = √(β·c)·dt`, derived in the module comment.
fn rk4_dissipation_bound(dt_s: f64, steps: usize) -> f64 {
    let theta = equatorial_wave_frequency_per_s() * dt_s;
    steps as f64 * theta.powi(6) / RK4_ENERGY_LOSS_DENOMINATOR
}

/// The bound on the relative energy drift of an undamped, unforced run of
/// `steps` steps of `dt_s` seconds on a grid of meridional spacing `dy_m`.
///
/// The sum of the two terms of the module comment: a spatial one that does not
/// grow with the run, and a temporal one that grows linearly in the step count
/// and falls as `dt⁵` for a fixed run length.
fn derived_drift_bound(dy_m: f64, dt_s: f64, steps: usize) -> f64 {
    coriolis_skewness_bound(dy_m) + rk4_dissipation_bound(dt_s, steps)
}

/// What one resolution of the conservation study measured.
#[derive(Debug, Clone, Copy)]
struct ConservationRun {
    /// Cells across the basin in `y`.
    cells_y: usize,
    /// Meridional cell spacing `Δy`, in metres.
    dy_m: f64,
    /// Length of one step, in seconds — the CFL-safe maximum for this grid.
    dt_s: f64,
    /// Steps taken, covering [`LONG_RUN_CROSSINGS`] basin crossings.
    steps: usize,
    /// Largest `|E(t)/E(0) − 1|` over the run, in either direction.
    worst_energy_drift: f64,
    /// Largest `E(t)/E(0) − 1` over the run, counting only the excursions in
    /// which the energy *rose*. Zero for a run that never gained any.
    worst_energy_gain: f64,
    /// Largest kinetic fraction `KE(t)/E(0)` the run reached.
    peak_kinetic_fraction: f64,
}

/// One undamped, unforced run of [`LONG_RUN_CROSSINGS`] basin crossings at
/// `cells_y` meridional cells, at the CFL-safe timestep.
fn measure_conservation_run(cells_y: usize) -> ConservationRun {
    let (grid, spacing) = basin(cells_y);
    let params = undamped_pacific_params();
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    let dt_s = max_stable_dt(spacing, wave_speed);
    let run_length_s = LONG_RUN_CROSSINGS * BASIN_LX_M / kelvin_wave_speed_m_per_s();
    let steps = (run_length_s / dt_s).round() as usize;

    let mut solver = Solver::new(
        grid,
        spacing,
        params,
        BetaPlane::centered_on_equator(params, spacing, grid),
        dt_s,
    )
    .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"));

    let mut state = gravest_zonal_mode(grid, spacing);
    let initial_energy = wave_energy(&state, params);
    let calm = WindStressField::calm(grid);

    let mut worst_energy_drift = 0.0_f64;
    let mut worst_energy_gain = 0.0_f64;
    let mut peak_kinetic_fraction = 0.0_f64;
    for n in 0..steps {
        solver.step(&mut state, n as f64 * dt_s, |_t_s| &calm);
        let relative_drift = wave_energy(&state, params) / initial_energy - 1.0;
        worst_energy_drift = worst_energy_drift.max(relative_drift.abs());
        worst_energy_gain = worst_energy_gain.max(relative_drift);
        peak_kinetic_fraction =
            peak_kinetic_fraction.max(kinetic_energy(&state, params) / initial_energy);
    }

    ConservationRun {
        cells_y,
        dy_m: spacing.dy_m(),
        dt_s,
        steps,
        worst_energy_drift,
        worst_energy_gain,
        peak_kinetic_fraction,
    }
}

/// The study, run once and read by every test below.
///
/// The runs are the expensive part of this file — tens of thousands of steps
/// at the finest resolution — and every assertion here is about the same three
/// of them, so they are integrated once and shared. The runs are deterministic
/// (CODING_STANDARDS.md § Correctness and failure), so neither sharing them nor
/// integrating them on one thread each can make one test's result depend on
/// another's; the finest resolution costs more than the other two together, so
/// running them side by side is most of the file's wall clock.
fn conservation_study() -> &'static [ConservationRun] {
    static STUDY: OnceLock<Vec<ConservationRun>> = OnceLock::new();
    STUDY.get_or_init(|| {
        std::thread::scope(|scope| {
            let runs: Vec<_> = MERIDIONAL_RESOLUTIONS
                .iter()
                .map(|&cells_y| scope.spawn(move || measure_conservation_run(cells_y)))
                .collect();
            runs.into_iter()
                .map(|run| run.join().expect("a conservation run must not panic"))
                .collect()
        })
    })
}

#[test]
fn energy_drift_over_a_long_undamped_unforced_run_stays_within_the_derived_bound() {
    // The acceptance criterion. `r = 0` and `τ = 0` leave a system whose only
    // departures from exact energy conservation are the C-grid Coriolis pair's
    // `O(Δy²)` skewness defect and RK4's `O(dt⁶)` dissipation; the module
    // comment derives both, and their sum is what the run has to stay inside.
    for run in conservation_study() {
        let bound = derived_drift_bound(run.dy_m, run.dt_s, run.steps);
        assert!(
            run.worst_energy_drift <= bound,
            "at {}x{} cells the energy drifted by a relative {} over {} steps \
             ({LONG_RUN_CROSSINGS} basin crossings), past the derived bound {bound} \
             (skewness {}, RK4 dissipation {})",
            BASIN_CELLS_X,
            run.cells_y,
            run.worst_energy_drift,
            run.steps,
            coriolis_skewness_bound(run.dy_m),
            rk4_dissipation_bound(run.dt_s, run.steps),
        );
    }
}

#[test]
fn the_long_run_energy_drift_falls_at_the_schemes_second_order_under_refinement() {
    // A bound can be met for the wrong reason — by a drift that is small but
    // does not come from the truncation error the bound is written about. This
    // is the check that it does: the leading term is `O(Δy²)`, so halving `Δy`
    // must quarter the drift (CODING_STANDARDS.md § Tests).
    let drifts: Vec<f64> = conservation_study()
        .iter()
        .map(|run| run.worst_energy_drift)
        .collect();

    for (coarse, pair) in drifts.windows(2).enumerate() {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order >= 2.0 - DRIFT_ORDER_TOLERANCE,
            "the energy drift must fall at least at second order under meridional \
             refinement, but {} to {} meridional cells gives order {order} (drifts {drifts:?})",
            MERIDIONAL_RESOLUTIONS[coarse],
            MERIDIONAL_RESOLUTIONS[coarse + 1],
        );
    }
}

#[test]
fn the_conservation_run_exchanges_energy_between_potential_and_kinetic_form() {
    // The guard against a vacuous pass: a run that never moved would conserve
    // energy exactly and mean nothing. The initial state is pure potential
    // energy at rest, so a working step has to turn most of it into motion.
    for run in conservation_study() {
        assert!(
            run.peak_kinetic_fraction >= KINETIC_EXCHANGE_FLOOR,
            "at {}x{} cells only a fraction {} of the initial energy ever reached \
             kinetic form; the run cannot be said to have exercised the equations",
            BASIN_CELLS_X,
            run.cells_y,
            run.peak_kinetic_fraction,
        );
    }
}

#[test]
fn an_undamped_unforced_run_never_gains_energy_past_the_skewness_bound() {
    // The half of the acceptance criterion that rests on nothing estimated.
    // RK4 has `|R(i·θ)| ≤ 1` everywhere inside the CFL bound the solver
    // enforces, so the time discretisation can only *remove* energy from this
    // system; every joule the basin gains has to come from the Coriolis pair's
    // skewness defect, whose excursion is the `(1/4)·(Δy/Le)²` of §1 alone.
    // No frequency estimate enters, so this bound holds whatever the run's
    // spectrum turns out to be.
    for run in conservation_study() {
        let bound = coriolis_skewness_bound(run.dy_m);
        assert!(
            run.worst_energy_gain <= bound,
            "at {}x{} cells the energy rose by a relative {} over {} steps \
             ({LONG_RUN_CROSSINGS} basin crossings), past the skewness bound {bound}; \
             RK4 cannot add energy, so nothing else could have put it there",
            BASIN_CELLS_X,
            run.cells_y,
            run.worst_energy_gain,
            run.steps,
        );
    }
}
