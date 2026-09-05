//! Acceptance tests for T-07.4 — the steady wind-driven thermocline tilt.
//!
//! `CONTEXT.md` (*Thermocline tilt*) states the claim this file validates: a
//! sustained easterly stress leaves the thermocline **deep in the west and
//! shallow in the east**, and that steady slope is "the observed mean state,
//! and the control case the model must reproduce". `docs/planning/01-scientific-model.md`
//! names the balance it comes from, and ADR-0003 fixes the discretisation
//! whose truncation error is what stands between the two.
//!
//! Nothing below is measured out of a run and pasted back in as an
//! expectation. Every asserted number is a closed-form solution of the
//! equations of ADR-0003, an exact conservation statement about them, or a
//! tolerance assembled from named error terms of the configuration
//! (CODING_STANDARDS.md § *Tests*).
//!
//! # The balance
//!
//! Written for a steady state, the linear 1.5-layer equations of ADR-0003 are
//!
//! ```text
//! r·u − f·v + g'·∂h/∂x = X,     X = τx/(ρ₀·H)
//! r·v + f·u + g'·∂h/∂y = 0
//! r·h + H·(∂u/∂x + ∂v/∂y) = 0.
//! ```
//!
//! The textbook tilt is what the first of those says when the current has
//! stopped: `g'·∂h/∂x = X`, a thermocline sloping down to the west under an
//! easterly stress. It is *not* an exact solution of the three, though, and
//! the reason is the third: a tilted thermocline is damped at `r·h`, so mass
//! has to be fed to it, so the current cannot be exactly zero after all. This
//! file's headline test is against the closed-form solution that keeps that
//! term rather than against the balance that drops it, and the single
//! parameter separating the two is
//!
//! ```text
//! δ = r·L/c,
//! ```
//!
//! the basin width in units of the distance `c/r` a signal covers in one
//! damping time.
//!
//! That pair — a wind stress closed by linear friction rather than by
//! advection — is the Sverdrup/Stommel-type balance the ticket asks for, read
//! on the equator where `f` vanishes and the meridional half of the Sverdrup
//! relation has nothing to say. The rotating basin, where it does, is the
//! second configuration below.
//!
//! # Where the closed form is exact: the zonal channel
//!
//! A basin one cell tall is all coast in the meridional direction, so
//! [`NoNormalFlow`](engine::NoNormalFlow) holds `v` at exactly zero on both of
//! its rows for all time. The rotation terms then contribute nothing — `f·v`
//! is zero in the zonal equation, and the `−f·u` of the meridional one is
//! discarded at the wall — and the system above collapses, with no
//! approximation at all, to the one-dimensional pair
//!
//! ```text
//! r·u + g'·dh/dx = X,        r·h + H·du/dx = 0,        u(0) = u(L) = 0.
//! ```
//!
//! Eliminating `u` gives `d²h/dx² = k²·h` with `k = r/c` and `c = √(g'·H)`,
//! and the two no-normal-flow conditions pick the odd solution and fix its
//! amplitude:
//!
//! ```text
//! h(x) = (A/k)·sinh(k·(x − L/2)) / cosh(k·L/2),      A = X/g' = τx/(ρ₀·H·g'),
//! u(x) = (X/r)·[1 − cosh(k·(x − L/2))/cosh(k·L/2)].
//! ```
//!
//! `A` is the slope of the undamped balance, and `A·L` the tilt it predicts.
//! As `k → 0` the solution above tends to `h = A·(x − L/2)` — the straight
//! thermocline of the textbook balance, with the zero basin mean that steady
//! state forces on it. At finite `δ` it is a sinh: the tilt is
//! `A·L·tanh(δ/2)/(δ/2)`, less than `A·L`, and the profile is measurably
//! curved. [`DAMPING_IN_BASIN_WIDTHS`] is chosen at `δ = 2` for exactly that
//! reason: a solver that had lost the `r·h` term and settled on the straight
//! line instead would sit `24%` of the tilt away at the walls — ten thousand
//! times the tolerance below — so it would fail loudly rather than pass.
//!
//! ## The discrete solution, and therefore the truncation term
//!
//! The same elimination runs through on the C-grid. With `h` at cell centres
//! and `u` on the faces between them, `h_i = C·sinh(κ·(x_i − L/2))` satisfies
//! the discrete pair exactly, for the wavenumber `κ` that solves
//!
//! ```text
//! sinh(κ·Δx/2) = k·Δx/2,     i.e.   κ = (2/Δx)·asinh(k·Δx/2),
//! ```
//!
//! and both wall conditions then hold at once, again exactly, with
//! `C = A/(k·cosh(κ·L/2))`. So the grid's entire error on this solution is
//! `κ` in place of `k`, which is `κ = k·(1 − (k·Δx)²/24 + …)`: second order,
//! as ADR-0003's scheme promises. [`Channel::zonal_truncation`] does not
//! expand that — it evaluates both closed forms on the run's own columns and
//! takes the largest difference between them, so the leading term is a
//! consequence of the tolerance rather than an assumption in it.
//!
//! Nothing else of the scheme contributes. A linear system's RK4 fixed point
//! is the fixed point of its right-hand side, so time stepping adds no error
//! *at* equilibrium; and `f·v` and the wind stress at the walls are exactly
//! zero rather than nearly so.
//!
//! # Reaching equilibrium, and saying so when it is not reached
//!
//! Damping every prognostic variable at the same rate `r` makes the unforced
//! system `ẋ = (L − r)·x` with `L` skew in the discrete energy
//! `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)` (`shallow_water.rs`), so `Ė = −2·r·E`
//! and any departure from the steady state decays as `e^{−r·t}` in that norm
//! however the waves shuffle it around. A run started from rest carries a
//! departure of exactly the steady state itself, so after `T` seconds
//!
//! ```text
//! ‖h(T) − h_steady‖₂ ≤ e^{−r·T}·‖h_steady‖₂ ≤ e^{−r·T}·√N·max|h_steady|,
//! ```
//!
//! and `r` is not merely the leading rate but the *only* one: every eigenvalue
//! of `L − r` has real part exactly `−r`, so there is no slow creep for a
//! drift measured over the run's second half to miss.
//!
//! which is [`Channel::equilibrium_bound`] and the twin for the basin. That
//! bound is the *whole* justification for calling a finite run "steady", and
//! the acceptance criterion asks for it to be checked rather than assumed. So
//! every run here is sampled twice — half way and at the end — and the drift
//! between the two samples is asserted against that bound before any profile
//! is compared.
//!
//! The guard is not decoration:
//! [`a_run_too_short_to_equilibrate_is_reported_as_not_equilibrated`] runs the
//! same channel for a twentieth of the spin-up and asserts both halves of the
//! criterion — that the guard says so, and that the profile such a run
//! produces really is outside the tolerance a silent pass would have accepted.
//!
//! # The rotating basin
//!
//! The channel is where the tilt has a closed form; it is not where the
//! equatorial Pacific is. The second configuration is the basin the shipped
//! control scenario describes — trade winds decaying away from the equator on
//! the deformation radius `Le = √(c/β)`, over a rotating beta-plane — and two
//! statements about it are exact rather than approximate, which is what makes
//! them worth asserting.
//!
//! **The basin carries no net anomaly.** Summing the continuity equation over
//! every cell telescopes the divergence to the flux through the four walls,
//! which no-normal-flow makes zero, leaving `r·Σh = 0`. A steady closed basin
//! has `Σh = 0` whatever shape the wind has, and
//! [`the_equilibrated_basin_carries_no_net_thermocline_anomaly`] holds the run
//! to that with a tolerance of the equilibrium bound and accumulated
//! round-off, nothing else.
//!
//! **The Kelvin invariant obeys a first-order balance.** Scale `y` by `Le`,
//! `x` by `Le`, `u` by `c` and `h` by `H`, and add the zonal-momentum and
//! continuity equations of the steady system:
//!
//! ```text
//! ε·q + ∂q/∂x̂ + (∂v̂/∂ŷ − ŷ·v̂) = F,      q = u/c + h/H,  ε = r·Le/c.
//! ```
//!
//! Projecting on `ψ₀ = e^{−ŷ²/2}` annihilates the meridional term exactly —
//! one integration by parts turns `∫(∂v̂/∂ŷ)ψ₀ dŷ` into `∫ŷ·v̂·ψ₀ dŷ` — so
//! with `q₀ = P₀[u]/c + P₀[h]/H` and `X₀ = P₀[τx]/(ρ₀·H)`,
//!
//! ```text
//! (r/c)·q₀ + dq₀/dx = X₀/c²
//! ```
//!
//! holds at every `x`, with no long-wave approximation and no truncation of
//! the Rossby set: the modes this relation does not mention are `ψ₂`, `ψ₄`, …
//! of `q` and `ψ₀`, `ψ₂`, … of `u/c − h/H`, and none of them has a `ψ₀`
//! component of `q` to contribute (Matsuno 1966; Gill, *Atmosphere–Ocean
//! Dynamics*, § 11.6 — the decomposition is `tests/support/mod.rs`'s, and the
//! same one the wave tests of T-07.1 and T-07.2 read). For a stress
//! `τx = τ₀·exp(−(y/Ly)²)` the right-hand side is analytic,
//! `X₀ = τ₀/(ρ₀·H·√((Le/Ly)² + 1/2))`, so the relation is a prediction and not
//! a tautology. [`the_kelvin_invariant_obeys_the_steady_damped_balance`]
//! asserts it in integrated form over the middle half of the basin, where a
//! quadrature replaces the derivative and no wall boundary layer is inside the
//! interval.
//!
//! # Where the tolerances come from
//!
//! Every entry below is a property of the configuration and none was obtained
//! by running anything. Collected, they are:
//!
//! | term | what it is | channel | basin |
//! |---|---|---|---|
//! | zonal truncation | `κ` for `k`, evaluated exactly | `≈ 2.6e−5` | — |
//! | meridional quadrature | `(Δy/Le)²`, `ψ₀`'s [`Waveguide::truncation_bound`] | — | `1.6e−2` |
//! | waveguide tail | `ψ₀` outside the basin's `±Ŷ` | — | `≈ 7e−5` |
//! | zonal quadrature | trapezoid and face averaging, `(k·Δx)²·(1/12 + 1/8)` | — | `≈ 1.3e−4` |
//! | equilibrium | `(e^{−r·T} + e^{−r·T/2})·√N` | `≈ 2e−8` | `≈ 3e−4` |
//! | round-off | `steps·N·ε` | negligible | negligible |
//!
//! and the second-order claim itself is asserted rather than assumed:
//! [`the_channel_tilt_error_falls_at_second_order_with_resolution`] halves the
//! cell width and requires the error to fall by the ratio the two closed forms
//! predict.

mod support;

use std::sync::OnceLock;

use engine::{
    Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, SteadyTradeWinds, WaveSpeed,
    H_STAGGERING,
};
use termocline_numerics::{max_stable_dt, Spacing};

use support::{
    equatorial_deformation_radius_m, kelvin_wave_speed_m_per_s, pacific_damped_params,
    zonal_current_at_cell_centre_in_c, MeridionalStructure, Waveguide,
};

/// Zonal stress `τ₀` of the alizés on the equator, in Pa.
///
/// The mean easterly stress over the equatorial Pacific, and the value
/// `engine/scenarios/steady-trades.toml` runs the control scenario at.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// `δ = r·L/c`, the basin width in units of the distance a signal covers in
/// one damping time — the one parameter that separates the damped closed form
/// from the textbook balance.
///
/// Two, which is a hard ocean rather than a realistic one: the shipped control
/// scenario sits at `δ ≈ 0.6`. It is deliberately hard. At `δ = 2` the damped
/// tilt is `tanh(1) = 76%` of the undamped one and the profile in between is
/// visibly a sinh, so a solver that had lost the `r·h` term of the continuity
/// equation — the term that makes the two differ at all — would miss by `24%`
/// of the tilt rather than hide inside a tolerance of `2.5e−5`. It also
/// shortens the spin-up, which is proportional to `1/r`.
const DAMPING_IN_BASIN_WIDTHS: f64 = 2.0;

/// How many damping times `1/r` a run is integrated for before its state is
/// called steady.
///
/// The transient decays as `e^{−r·T}` in the energy norm (module header), so
/// this is what sets the equilibrium term of every tolerance below. Forty puts
/// that term three orders of magnitude under the truncation term it is added
/// to, which is what keeps the convergence test measuring the grid rather than
/// the clock.
const SPIN_UP_IN_DAMPING_TIMES: f64 = 40.0;

/// The spin-up of the run that is deliberately too short to have equilibrated.
///
/// Two damping times leaves `e^{−2} ≈ 14%` of the initial departure in place —
/// a transient, not a steady state, and the thing the acceptance criterion
/// asks this file to refuse rather than average over.
const TOO_SHORT_SPIN_UP_IN_DAMPING_TIMES: f64 = 2.0;

/// Zonal extent `L` of the channel, in metres: 15 000 km, the width of the
/// equatorial Pacific of `CONTEXT.md` to two figures.
const CHANNEL_LENGTH_M: f64 = 1.5e7;

/// Columns of the coarse channel. The fine run has twice as many.
const COARSE_CHANNEL_COLUMNS: usize = 60;

/// How much of the coarse cell width the refined run uses.
const CHANNEL_REFINEMENT: usize = 2;

// ---------------------------------------------------------------------------
// The zonal channel: where the damped tilt has a closed form
// ---------------------------------------------------------------------------

/// One tilt experiment: a zonal channel at some refinement, forced by steady
/// alizés with no meridional structure at all.
///
/// One cell tall, which is what makes it a channel rather than a basin: both
/// of its `v` rows are coast, so `v` is held at exactly zero and the rotation
/// terms drop out of the system rather than being small in it (module header).
#[derive(Debug, Clone, Copy)]
struct Channel {
    /// Shape, spacing and position of the channel.
    basin: Basin,
    /// The ocean the equations are written in terms of, damped at `r`.
    params: PhysicalParams,
}

impl Channel {
    /// The channel at `1/refinement` of the coarse cell width.
    fn new(refinement: usize) -> Self {
        let columns = COARSE_CHANNEL_COLUMNS * refinement;
        let cell_m = CHANNEL_LENGTH_M / columns as f64;
        let grid = Grid::new(columns, 1).expect("the channel has cells on both axes");
        let spacing = Spacing::new(cell_m, cell_m).expect("the cell size is finite and positive");
        Self {
            basin: Basin::centered_on_equator(grid, spacing),
            params: pacific_damped_params(
                DAMPING_IN_BASIN_WIDTHS * kelvin_wave_speed_m_per_s() / CHANNEL_LENGTH_M,
            ),
        }
    }

    /// Number of columns.
    fn columns(self) -> usize {
        self.basin.grid().nx()
    }

    /// Cell width `Δx`, in metres.
    fn cell_m(self) -> f64 {
        self.basin.spacing().dx_m()
    }

    /// The slope `A = τ₀/(ρ₀·H·g')` of the undamped balance, in m/m.
    ///
    /// Negative under easterly stress: the thermocline deepens westward.
    fn undamped_slope_per_m(self) -> f64 {
        TRADE_WIND_STRESS_PA
            / (self.params.reference_density_kg_per_m3()
                * self.params.mean_thermocline_depth_m()
                * self.params.reduced_gravity_m_per_s2())
    }

    /// The damping wavenumber `k = r/c`, in m⁻¹ — the inverse of the distance
    /// a signal covers in one damping time.
    fn damping_wavenumber_per_m(self) -> f64 {
        self.params.rayleigh_damping_per_s() / kelvin_wave_speed_m_per_s()
    }

    /// The wavenumber `κ = (2/Δx)·asinh(k·Δx/2)` the *discrete* pair admits in
    /// place of `k`, in m⁻¹ (module header).
    fn discrete_wavenumber_per_m(self) -> f64 {
        let half_cell_in_wavenumbers = 0.5 * self.damping_wavenumber_per_m() * self.cell_m();
        2.0 * half_cell_in_wavenumbers.asinh() / self.cell_m()
    }

    /// `h(x) = (A/k)·sinh(κ·(x − L/2))/cosh(κ·L/2)`, in metres: the closed-form
    /// steady thermocline anomaly at zonal wavenumber `wavenumber_per_m`.
    ///
    /// The amplitude carries `k` and the shape carries the wavenumber asked
    /// for, which is what the discrete elimination gives: passing `k` returns
    /// the continuous solution the assertions are made against, and passing
    /// `κ` the one the grid can represent.
    fn analytic_depth_anomaly_m(self, x_m: f64, wavenumber_per_m: f64) -> f64 {
        let half_width = 0.5 * self.basin.zonal_extent_m();
        self.undamped_slope_per_m() / self.damping_wavenumber_per_m()
            * (wavenumber_per_m * (x_m - half_width)).sinh()
            / (wavenumber_per_m * half_width).cosh()
    }

    /// The closed-form profile at `wavenumber_per_m`, sampled on the channel's
    /// own columns, in metres.
    fn analytic_profile_m(self, wavenumber_per_m: f64) -> Vec<f64> {
        (0..self.columns())
            .map(|i| {
                let x_m = self.basin.x_of_column_m(H_STAGGERING, i);
                self.analytic_depth_anomaly_m(x_m, wavenumber_per_m)
            })
            .collect()
    }

    /// The zonal truncation term, as a fraction of the tilt: the largest gap
    /// between the discrete closed form and the continuous one, over the
    /// channel's columns.
    ///
    /// Second order in `Δx` by construction — `κ = k·(1 − (k·Δx)²/24 + …)` —
    /// but evaluated rather than expanded, so no term of that series is
    /// dropped on the way into a tolerance.
    fn zonal_truncation(self) -> f64 {
        relative_deviation(
            &self.analytic_profile_m(self.discrete_wavenumber_per_m()),
            &self.analytic_profile_m(self.damping_wavenumber_per_m()),
        )
    }

    /// The equilibrium term, as a fraction of the tilt: what the module
    /// header's `e^{−r·T}·√N` energy bound leaves of the initial transient
    /// after `spin_up_in_damping_times`, at the half-way sample and the final
    /// one together.
    fn equilibrium_bound(self, spin_up_in_damping_times: f64) -> f64 {
        equilibrium_bound(spin_up_in_damping_times, self.columns())
    }

    /// The round-off term, as a fraction of the tilt.
    fn round_off(self, steps: usize) -> f64 {
        round_off(steps, self.columns())
    }

    /// The tolerance on the whole profile, as a fraction of the tilt: the
    /// three terms of the module header's table.
    ///
    /// `spin_up_in_damping_times` is the run the tolerance is *derived for*,
    /// which is not always the run it is applied to: what makes a transient a
    /// failure is that it misses the tolerance an equilibrated run was
    /// promised, not that it misses one relaxed to fit its own youth.
    fn profile_tolerance(self, spin_up_in_damping_times: f64, steps: usize) -> f64 {
        self.zonal_truncation()
            + self.equilibrium_bound(spin_up_in_damping_times)
            + self.round_off(steps)
    }

    /// Integrate the channel from rest under steady, meridionally uniform
    /// alizés for `spin_up_in_damping_times` damping times, sampling the
    /// thermocline half way and at the end.
    fn run(self, spin_up_in_damping_times: f64) -> Integration<Vec<f64>> {
        let wind = SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA)
            .expect("the alizés of the control scenario are easterly");
        integrate(
            self.basin,
            self.params,
            &wind,
            spin_up_in_damping_times,
            |state| equator_row(self.basin, state),
        )
    }
}

// ---------------------------------------------------------------------------
// The rotating equatorial basin: the shipped control scenario's forcing
// ---------------------------------------------------------------------------

/// Zonal extent of the basin, in deformation radii. Forty `Le` is 13 800 km,
/// the equatorial Pacific to the accuracy this configuration cares about.
const BASIN_LENGTH_IN_RADII: f64 = 40.0;

/// Meridional half-extent of the basin, in deformation radii.
///
/// Four, which puts the walls where `ψ₀` has fallen to `e^{−8}` of its
/// equatorial value: far enough that the waveguide is complete inside the
/// basin, and the [`EquatorialBasin::waveguide_tail`] term of the tolerance is
/// what quantifies the remainder.
const BASIN_HALF_WIDTH_IN_RADII: f64 = 4.0;

/// Columns per deformation radius.
const BASIN_COLUMNS_PER_RADIUS: f64 = 2.0;

/// Rows per deformation radius.
///
/// Eight, because the meridional quadrature of the `ψ₀` projection is the
/// leading term of this configuration's tolerance and it is second order:
/// `(Δy/Le)² = 1/64`.
const BASIN_ROWS_PER_RADIUS: f64 = 8.0;

/// Meridional decay scale `Ly` of the alizés, in deformation radii.
///
/// One, as in `engine/scenarios/steady-trades.toml`: the stress is confined to
/// the waveguide it is meant to drive.
const TRADE_WIND_DECAY_IN_RADII: f64 = 1.0;

/// How many damping times the basin is integrated for.
///
/// Shorter than [`SPIN_UP_IN_DAMPING_TIMES`] because it can be: this
/// configuration's tolerance is dominated by the meridional quadrature, and
/// twenty-five damping times already puts the equilibrium term two orders of
/// magnitude below it.
const BASIN_SPIN_UP_IN_DAMPING_TIMES: f64 = 25.0;

/// The rotating equatorial basin the control scenario describes.
#[derive(Debug, Clone, Copy)]
struct EquatorialBasin {
    /// Shape, spacing and position of the basin.
    basin: Basin,
    /// The ocean the equations are written in terms of, damped at `r`.
    params: PhysicalParams,
    /// `Le = √(c/β)`, in metres — cached because every projection needs it.
    deformation_radius_m: f64,
}

impl EquatorialBasin {
    /// The basin of the constants above.
    fn new() -> Self {
        let radius_m = equatorial_deformation_radius_m();
        let cell_width_m = radius_m / BASIN_COLUMNS_PER_RADIUS;
        let cell_height_m = radius_m / BASIN_ROWS_PER_RADIUS;
        let grid = Grid::new(
            (BASIN_LENGTH_IN_RADII * BASIN_COLUMNS_PER_RADIUS).round() as usize,
            (2.0 * BASIN_HALF_WIDTH_IN_RADII * BASIN_ROWS_PER_RADIUS).round() as usize,
        )
        .expect("the basin has cells on both axes");
        let spacing = Spacing::new(cell_width_m, cell_height_m)
            .expect("the cell sizes are finite and positive");
        let basin = Basin::centered_on_equator(grid, spacing);
        Self {
            basin,
            params: pacific_damped_params(
                DAMPING_IN_BASIN_WIDTHS * kelvin_wave_speed_m_per_s() / basin.zonal_extent_m(),
            ),
            deformation_radius_m: radius_m,
        }
    }

    /// The alizés of the control scenario: easterly, decaying away from the
    /// equator on [`TRADE_WIND_DECAY_IN_RADII`] deformation radii.
    fn wind(self) -> SteadyTradeWinds {
        SteadyTradeWinds::with_meridional_decay(
            TRADE_WIND_STRESS_PA,
            TRADE_WIND_DECAY_IN_RADII * self.deformation_radius_m,
        )
        .expect("the alizés of the control scenario are easterly and decay on a positive scale")
    }

    /// `X₀ = P₀[τx]/(ρ₀·H)`, in m/s²: the `ψ₀` coefficient of the surface
    /// stress the Kelvin balance is forced by.
    ///
    /// Analytic. For `τx = τ₀·exp(−(y/Ly)²)` the projection
    /// `∫τx·ψ₀ dŷ / ∫ψ₀² dŷ` is two Gaussian integrals,
    /// `τ₀·√(π/(a + 1/2))/√π` with `a = (Le/Ly)²`, which is the
    /// `1/√(a + 1/2)` below.
    fn gravest_stress_forcing_m_per_s2(self) -> f64 {
        let decay_in_radii_squared = 1.0 / (TRADE_WIND_DECAY_IN_RADII * TRADE_WIND_DECAY_IN_RADII);
        TRADE_WIND_STRESS_PA
            / (self.params.reference_density_kg_per_m3() * self.params.mean_thermocline_depth_m())
            / (decay_in_radii_squared + 0.5).sqrt()
    }

    /// The waveguide the `ψ₀` projections are taken on.
    fn waveguide(self) -> Waveguide {
        Waveguide::new(self.basin, self.params)
    }

    /// `q₀ = P₀[u]/c + P₀[h]/H`, column by column: the Kelvin invariant of a
    /// state, dimensionless.
    ///
    /// `u` is read at the cell centre by
    /// [`zonal_current_at_cell_centre_in_c`], so that it and the depth anomaly
    /// sit at one set of positions — the same average
    /// `support::gravest_current_projection` takes. What differs is the
    /// normalisation: [`Waveguide::column_coefficients`] divides by the
    /// analytic `∫ψ₀² dŷ`, because what this is compared against is a stress
    /// coefficient written down from theory.
    fn kelvin_invariant(self, state: &OceanState) -> Vec<f64> {
        let waveguide = self.waveguide();
        let columns = self.basin.grid().nx();
        let wave_speed_m_per_s = kelvin_wave_speed_m_per_s();
        let mean_depth_m = self.params.mean_thermocline_depth_m();
        let current =
            waveguide.column_coefficients(columns, MeridionalStructure::Gravest, |i, j| {
                zonal_current_at_cell_centre_in_c(state, i, j, wave_speed_m_per_s)
            });
        let depth = waveguide.column_coefficients(columns, MeridionalStructure::Gravest, |i, j| {
            state.h().get(i, j).expect("a cell centre") / mean_depth_m
        });
        current
            .iter()
            .zip(&depth)
            .map(|(current, depth)| current + depth)
            .collect()
    }

    /// `P₀[h]/H`, column by column: the waveguide's own thermocline anomaly,
    /// dimensionless.
    ///
    /// The equatorial profile of a rotating basin, read the way the waveguide
    /// defines it. A basin with an even number of rows has no cell row *on*
    /// the equator — which is what makes its meridional quadrature symmetric —
    /// so the `ψ₀` coefficient is the equatorial anomaly, not a substitute for
    /// one.
    fn gravest_depth_projection(self, state: &OceanState) -> Vec<f64> {
        let mean_depth_m = self.params.mean_thermocline_depth_m();
        self.waveguide().column_coefficients(
            self.basin.grid().nx(),
            MeridionalStructure::Gravest,
            |i, j| state.h().get(i, j).expect("a cell centre") / mean_depth_m,
        )
    }

    /// The meridional quadrature term, as a fraction: `ψ₀`'s second-order
    /// [`Waveguide::truncation_bound`], `(Δy/Le)²`.
    fn meridional_truncation(self) -> f64 {
        self.waveguide()
            .truncation_bound(MeridionalStructure::Gravest)
    }

    /// The waveguide-tail term, as a fraction: the share of `∫ψ₀ dŷ` that lies
    /// outside the basin's walls and is therefore missing from every
    /// projection.
    ///
    /// `∫_Ŷ^∞ e^{−ŷ²/2} dŷ ≤ e^{−Ŷ²/2}/Ŷ` on both flanks, over the whole
    /// `∫ψ₀ dŷ = √(2π)`.
    fn waveguide_tail(self) -> f64 {
        let half_width = BASIN_HALF_WIDTH_IN_RADII;
        2.0 * (-0.5 * half_width * half_width).exp()
            / (half_width * (2.0 * std::f64::consts::PI).sqrt())
    }

    /// The zonal quadrature term, as a fraction: `(k·Δx)²·(1/12 + 1/8)`.
    ///
    /// Two second-order errors on the same cell width, both taken against the
    /// only zonal scale the Kelvin invariant has — `c/r`, since `q₀` obeys a
    /// first-order equation with an `x`-independent right-hand side and is
    /// therefore an exponential of that scale and nothing shorter. `1/12` is
    /// the trapezoid rule's on `∫q₀ dx`, and `1/8` the face-to-centre average
    /// of `u`.
    fn zonal_quadrature(self) -> f64 {
        let cell_in_wavenumbers =
            self.params.rayleigh_damping_per_s() / kelvin_wave_speed_m_per_s() * self.cell_m();
        cell_in_wavenumbers * cell_in_wavenumbers * (1.0 / 12.0 + 1.0 / 8.0)
    }

    /// The equilibrium term, as a fraction.
    fn equilibrium_bound(self) -> f64 {
        equilibrium_bound(BASIN_SPIN_UP_IN_DAMPING_TIMES, self.cells())
    }

    /// The round-off term, as a fraction.
    fn round_off(self, steps: usize) -> f64 {
        round_off(steps, self.cells())
    }

    /// The columns of the basin's middle half, west to east inclusive.
    ///
    /// The interval every zonal statement of this file is made over. Both
    /// walls carry a reflection layer that long-wave theory makes
    /// *infinitely thin* — the Rossby modes a wall generates decay eastward on
    /// `c/((2n + 1)·r)`, which shortens without limit in `n` — so no grid
    /// resolves one, and the few columns nearest a wall carry the grid's
    /// rendering of it rather than the interior solution. A quarter of the
    /// basin either side is a wide berth: at `δ = 2` even the `n = 1` mode has
    /// decayed by `e^{−3/2}` across it, and every higher one by more.
    fn interior_columns(self) -> (usize, usize) {
        let columns = self.basin.grid().nx();
        (columns / 4, 3 * columns / 4)
    }

    /// How many cells the basin has.
    fn cells(self) -> usize {
        self.basin.grid().nx() * self.basin.grid().ny()
    }

    /// Cell width `Δx`, in metres.
    fn cell_m(self) -> f64 {
        self.basin.spacing().dx_m()
    }

    /// Zonal position of column `i`'s centre, in metres east of the western
    /// wall.
    fn column_x_m(self, i: usize) -> f64 {
        self.basin.x_of_column_m(H_STAGGERING, i)
    }

    /// The tolerance on the Kelvin balance's residual, as a fraction of the
    /// forcing it is balanced against: the four terms of the module header's
    /// table that bear on a projected quantity.
    fn kelvin_tolerance(self, steps: usize) -> f64 {
        self.meridional_truncation()
            + self.waveguide_tail()
            + self.zonal_quadrature()
            + self.equilibrium_bound()
            + self.round_off(steps)
    }

    /// Integrate the basin from rest under the control scenario's alizés,
    /// sampling the whole state half way and at the end.
    fn run(self) -> Integration<OceanState> {
        let wind = self.wind();
        integrate(
            self.basin,
            self.params,
            &wind,
            BASIN_SPIN_UP_IN_DAMPING_TIMES,
            Clone::clone,
        )
    }
}

/// Refuse to read `run` as a steady state unless its thermocline stopped
/// moving over the second half of the integration.
///
/// The acceptance criterion of T-07.4 in the rotating basin: three tests ask
/// three questions of one equilibrium, and each of them would have an answer —
/// a wrong one — on a transient. The drift is read on `h` itself rather than
/// on the projection each test happens to want, because `h` is the prognostic
/// variable the criterion is about and the energy bound is stated for it.
fn assert_basin_has_settled(basin: EquatorialBasin, run: &Integration<OceanState>) {
    let drift = relative_deviation(run.half_way.h().as_slice(), run.settled.h().as_slice());
    let bound = basin.equilibrium_bound();
    assert!(
        drift <= bound,
        "the basin is still adjusting after {BASIN_SPIN_UP_IN_DAMPING_TIMES} damping times: its \
         thermocline moved by {drift:.3e} of the tilt over the second half of the run, against \
         the {bound:.3e} the energy bound allows. This run has not reached equilibrium, so \
         nothing read from it is a steady state"
    );
}

/// The basin run, integrated once and shared by the tests that read it.
///
/// Three tests ask three different questions of one equilibrium, and the
/// integration is the expensive part of all three.
fn settled_basin() -> &'static Integration<OceanState> {
    static RUN: OnceLock<Integration<OceanState>> = OnceLock::new();
    RUN.get_or_init(|| EquatorialBasin::new().run())
}

// ---------------------------------------------------------------------------
// The integration, and the terms every configuration's tolerance shares
// ---------------------------------------------------------------------------

/// One run's two samples, in whatever form the caller asked for them.
struct Integration<T> {
    /// How many steps the run took.
    steps: usize,
    /// The sample half way through.
    half_way: T,
    /// The sample at the end.
    settled: T,
}

impl Integration<Vec<f64>> {
    /// How far the profile still moved over the second half of the run, as a
    /// fraction of the tilt — what the equilibrium guard reads.
    fn drift(&self) -> f64 {
        relative_deviation(&self.half_way, &self.settled)
    }

    /// How far the settled profile is from `analytic_m`, as a fraction of the
    /// tilt.
    fn departure_from(&self, analytic_m: &[f64]) -> f64 {
        relative_deviation(&self.settled, analytic_m)
    }
}

/// Integrate `basin` from rest under `wind` for `spin_up_in_damping_times`
/// damping times, reading `sample` from the state half way through and at the
/// end.
///
/// The step count is made even so that the half-way sample lands on a step
/// rather than between two, which is what lets the drift the equilibrium guard
/// reads be a difference of two states the solver actually produced.
///
/// # Panics
/// If the timestep the CFL bound hands back is one the solver refuses, which
/// would mean the two bounds of `solver.rs` disagree rather than that the run
/// is misconfigured.
fn integrate<T>(
    basin: Basin,
    params: PhysicalParams,
    wind: &SteadyTradeWinds,
    spin_up_in_damping_times: f64,
    sample: impl Fn(&OceanState) -> T,
) -> Integration<T> {
    let wave_speed =
        WaveSpeed::new(kelvin_wave_speed_m_per_s()).expect("the Kelvin wave speed is positive");
    let dt_s = max_stable_dt(basin.spacing(), wave_speed);
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());
    let mut solver = Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
        .unwrap_or_else(|error| {
            panic!("the experiment's own timestep must be admissible: {error}")
        });

    let spin_up_s = spin_up_in_damping_times / params.rayleigh_damping_per_s();
    let half_steps = (0.5 * spin_up_s / dt_s).ceil() as usize;
    let steps = 2 * half_steps;

    let mut state = OceanState::at_rest(basin.grid());
    let mut half_way = None;
    for step in 0..steps {
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, wind);
        if step + 1 == half_steps {
            half_way = Some(sample(&state));
        }
    }
    Integration {
        steps,
        half_way: half_way.expect("the half-way step is inside the run"),
        settled: sample(&state),
    }
}

/// The equilibrium term of a tolerance, as a fraction: what the module
/// header's energy bound leaves of the initial transient at the two sample
/// times together.
///
/// `√N` is the price of turning an energy-norm bound into a pointwise one —
/// `max|h| ≤ ‖h‖₂` and `‖h_steady‖₂ ≤ √N·max|h_steady|` — and it is kept
/// rather than dropped because it is what makes the bound true of every cell
/// rather than of the basin on average.
fn equilibrium_bound(spin_up_in_damping_times: f64, cells: usize) -> f64 {
    ((-spin_up_in_damping_times).exp() + (-0.5 * spin_up_in_damping_times).exp())
        * (cells as f64).sqrt()
}

/// What the energy bound leaves of the initial transient at the end of a run
/// of `spin_up_in_damping_times` damping times, as a fraction of the steady
/// state — in the *mean*, where the `√N` of [`equilibrium_bound`] is exactly
/// cancelled by the averaging: `|mean h'| ≤ ‖h'‖₂/√N ≤ e^{−r·T}·max|h_steady|`.
fn settled_transient(spin_up_in_damping_times: f64) -> f64 {
    (-spin_up_in_damping_times).exp()
}

/// The round-off term of a tolerance, as a fraction: one rounding per cell per
/// step, none of them cancelling.
///
/// The most an `f64` run can accumulate, and three orders of magnitude looser
/// than a random walk over the same arithmetic would give.
fn round_off(steps: usize, cells: usize) -> f64 {
    steps as f64 * cells as f64 * f64::EPSILON
}

/// The thermocline anomaly along the equatorial row of a state, in metres.
///
/// # Panics
/// If the basin has more than one row, which would mean the caller took a
/// channel's measurement of something that is not a channel.
fn equator_row(basin: Basin, state: &OceanState) -> Vec<f64> {
    assert_eq!(
        basin.grid().ny(),
        1,
        "the equatorial row of a channel is its only row"
    );
    (0..basin.grid().nx())
        .map(|i| *state.h().get(i, 0).expect("a cell centre"))
        .collect()
}

/// The largest gap between two profiles, as a fraction of the largest value the
/// `reference` one reaches.
///
/// # Panics
/// If the profiles are of different lengths, or if the reference profile is
/// everywhere zero and so has no scale to measure a gap against.
fn relative_deviation(measured: &[f64], reference: &[f64]) -> f64 {
    assert_eq!(
        measured.len(),
        reference.len(),
        "two profiles of one basin have one length"
    );
    let scale = reference
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    assert!(
        scale > 0.0,
        "the reference profile is everywhere zero, so it sets no scale"
    );
    measured
        .iter()
        .zip(reference)
        .fold(0.0_f64, |largest, (measured, reference)| {
            largest.max((measured - reference).abs())
        })
        / scale
}

/// `∫ values dx` by the trapezoid rule at spacing `step_m`.
///
/// # Panics
/// If fewer than two samples are handed in, which would mean the interval the
/// caller asked to integrate over is not an interval.
fn trapezoid(values: &[f64], step_m: f64) -> f64 {
    assert!(
        values.len() >= 2,
        "an integration interval needs at least two samples"
    );
    let ends = 0.5 * (values[0] + values[values.len() - 1]);
    let interior: f64 = values[1..values.len() - 1].iter().sum();
    (ends + interior) * step_m
}

// ---------------------------------------------------------------------------
// The tests
// ---------------------------------------------------------------------------

#[test]
fn the_steady_channel_tilt_matches_the_analytic_damped_balance() {
    let channel = Channel::new(1);
    let run = channel.run(SPIN_UP_IN_DAMPING_TIMES);

    let equilibrium = channel.equilibrium_bound(SPIN_UP_IN_DAMPING_TIMES);
    assert!(
        run.drift() <= equilibrium,
        "the channel is still moving after {SPIN_UP_IN_DAMPING_TIMES} damping times: the profile \
         changed by {:.3e} of the tilt over the second half of the run, against the {equilibrium:.3e} \
         the energy bound allows. This run has not reached equilibrium, so nothing below it is a \
         steady state",
        run.drift()
    );

    let analytic_m = channel.analytic_profile_m(channel.damping_wavenumber_per_m());
    let tolerance = channel.profile_tolerance(SPIN_UP_IN_DAMPING_TIMES, run.steps);
    let departure = run.departure_from(&analytic_m);
    assert!(
        departure <= tolerance,
        "the settled thermocline departs from the analytic damped tilt by {departure:.3e} of the \
         tilt, against a tolerance of {tolerance:.3e} = {:.3e} truncation + {:.3e} equilibrium + \
         {:.3e} round-off",
        channel.zonal_truncation(),
        equilibrium,
        channel.round_off(run.steps),
    );

    // The tilt is the headline number of `CONTEXT.md`: deep in the west,
    // shallow in the east. Its analytic value is `A·L·tanh(δ/2)/(δ/2)`, and
    // the profile has just been held to the closed form it comes from, so this
    // asserts the *direction* the comparison above cannot — that the sign
    // convention has not been reproduced backwards on both sides at once.
    let west_m = run.settled[0];
    let east_m = run.settled[run.settled.len() - 1];
    assert!(
        west_m > 0.0 && east_m < 0.0,
        "easterly alizés must leave the thermocline deep in the west and shallow in the east, \
         but the run ends at {west_m:.1} m against the western wall and {east_m:.1} m against \
         the eastern one"
    );
}

#[test]
fn the_channel_tilt_error_falls_at_second_order_with_resolution() {
    let coarse = Channel::new(1);
    let fine = Channel::new(CHANNEL_REFINEMENT);

    let departures = [coarse, fine].map(|channel| {
        let run = channel.run(SPIN_UP_IN_DAMPING_TIMES);
        let equilibrium = channel.equilibrium_bound(SPIN_UP_IN_DAMPING_TIMES);
        assert!(
            run.drift() <= equilibrium,
            "the {}-column channel has not reached equilibrium: it drifted {:.3e} of the tilt \
             over the second half of the run, against {equilibrium:.3e}",
            channel.columns(),
            run.drift()
        );
        run.departure_from(&channel.analytic_profile_m(channel.damping_wavenumber_per_m()))
    });

    // The predicted ratio is the two closed forms' own — `κ` against `k` at
    // each cell width — and it is `1/4` only to leading order: the exact
    // `asinh` carries `(k·Δx)⁴` terms too, so the prediction is evaluated
    // rather than assumed to be a quarter.
    let predicted = fine.zonal_truncation() / coarse.zonal_truncation();
    let measured = departures[1] / departures[0];
    // What the ratio of two truncation terms cannot account for: the
    // equilibrium floor under each of them, which is the larger share of the
    // finer error.
    let tolerance = fine.equilibrium_bound(SPIN_UP_IN_DAMPING_TIMES) / departures[1]
        + coarse.equilibrium_bound(SPIN_UP_IN_DAMPING_TIMES) / departures[0];
    assert!(
        (measured - predicted).abs() <= tolerance,
        "halving the cell width should divide the departure from the analytic tilt by \
         {predicted:.4} — the second-order ratio the discrete and continuous closed forms \
         predict — but {:.3e} became {:.3e}, a ratio of {measured:.4}, outside the {tolerance:.3e} \
         the equilibrium floor accounts for",
        departures[0],
        departures[1],
    );
}

#[test]
fn a_run_too_short_to_equilibrate_is_reported_as_not_equilibrated() {
    let channel = Channel::new(1);
    let run = channel.run(TOO_SHORT_SPIN_UP_IN_DAMPING_TIMES);

    // Both halves of the acceptance criterion. First: the guard fires. It is
    // the *settled* run's allowance the guard is read against — the question a
    // steady-state test asks is whether this run may be treated as the one the
    // tolerance was derived for, not whether it is self-consistently short.
    let equilibrium = channel.equilibrium_bound(SPIN_UP_IN_DAMPING_TIMES);
    assert!(
        run.drift() > equilibrium,
        "a run of {TOO_SHORT_SPIN_UP_IN_DAMPING_TIMES} damping times is a transient, but its \
         second half drifted only {:.3e} of the tilt — under the {equilibrium:.3e} an \
         equilibrated run is allowed, so the guard that is supposed to catch it would wave it \
         through",
        run.drift()
    );

    // Second: it fires on something that matters. A transient this young is
    // genuinely outside the tolerance, so a suite without the guard would not
    // merely be lucky — it would be wrong.
    let analytic_m = channel.analytic_profile_m(channel.damping_wavenumber_per_m());
    let tolerance = channel.profile_tolerance(SPIN_UP_IN_DAMPING_TIMES, run.steps);
    let departure = run.departure_from(&analytic_m);
    assert!(
        departure > tolerance,
        "a run of {TOO_SHORT_SPIN_UP_IN_DAMPING_TIMES} damping times should still be \
         {:.0}% short of the steady tilt, but it departs from the analytic profile by only \
         {departure:.3e}, inside the {tolerance:.3e} tolerance — which would make the \
         equilibrium guard the only thing standing between a transient and a green suite",
        100.0 * (-TOO_SHORT_SPIN_UP_IN_DAMPING_TIMES).exp()
    );
}

#[test]
fn the_equilibrated_basin_carries_no_net_thermocline_anomaly() {
    let basin = EquatorialBasin::new();
    let run = settled_basin();
    assert_basin_has_settled(basin, run);

    let depths_m = run.settled.h().as_slice();
    let cells = depths_m.len() as f64;
    let scale_m = depths_m
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    let mean_m = depths_m.iter().sum::<f64>() / cells;

    // Exact, not asymptotic: summing the discrete continuity equation over the
    // basin telescopes the divergence to the flux through four walls that
    // carry none, leaving `r·Σh = 0`. So the only things between the run and
    // zero are the transient it has not finished shedding and the arithmetic
    // it did on the way.
    let tolerance = settled_transient(BASIN_SPIN_UP_IN_DAMPING_TIMES) + basin.round_off(run.steps);
    assert!(
        (mean_m / scale_m).abs() <= tolerance,
        "a steady closed basin holds no net thermocline anomaly, but the settled state averages \
         {mean_m:.3e} m against a {scale_m:.1} m tilt — {:.3e} of it, against a tolerance of \
         {tolerance:.3e}",
        (mean_m / scale_m).abs()
    );
}

#[test]
fn the_kelvin_invariant_obeys_the_steady_damped_balance() {
    let basin = EquatorialBasin::new();
    let run = settled_basin();
    assert_basin_has_settled(basin, run);

    let settled = basin.kelvin_invariant(&run.settled);
    let equilibrium = basin.equilibrium_bound();

    // Over the basin's middle half, which leaves `q₀`'s own `c/r` as the only
    // zonal scale the quadrature terms are taken against.
    let (west, east) = basin.interior_columns();
    let cell_m = basin.cell_m();
    let span_m = basin.column_x_m(east) - basin.column_x_m(west);

    // `(r/c)·∫q₀ dx + [q₀] = X₀·span/c²`, the module header's balance
    // integrated over that interval.
    let damping_wavenumber_per_m =
        basin.params.rayleigh_damping_per_s() / kelvin_wave_speed_m_per_s();
    let wave_speed_m_per_s = kelvin_wave_speed_m_per_s();
    let forced = basin.gravest_stress_forcing_m_per_s2() * span_m
        / (wave_speed_m_per_s * wave_speed_m_per_s);
    let balanced = settled[east] - settled[west]
        + damping_wavenumber_per_m * trapezoid(&settled[west..=east], cell_m);

    let residual = (balanced - forced).abs() / forced.abs();
    let tolerance = basin.kelvin_tolerance(run.steps);
    assert!(
        residual <= tolerance,
        "the settled basin's Kelvin invariant misses its steady balance by {residual:.3e} of the \
         stress forcing it: {balanced:.6e} against the analytic {forced:.6e}. The tolerance is \
         {tolerance:.3e} = {:.3e} meridional quadrature + {:.3e} waveguide tail + {:.3e} zonal \
         quadrature + {:.3e} equilibrium + {:.3e} round-off",
        basin.meridional_truncation(),
        basin.waveguide_tail(),
        basin.zonal_quadrature(),
        equilibrium,
        basin.round_off(run.steps),
    );
}

#[test]
fn the_equatorial_thermocline_deepens_to_the_west_across_the_whole_basin() {
    let basin = EquatorialBasin::new();
    let run = settled_basin();
    assert_basin_has_settled(basin, run);

    // `CONTEXT.md`, *Thermocline tilt*: "deep in the west, shallow in the
    // east". Read on `ψ₀`, which is what the waveguide's own anomaly is.
    let profile = basin.gravest_depth_projection(&run.settled);
    let mean_depth_m = basin.params.mean_thermocline_depth_m();
    let (west, east) = (profile[0], profile[profile.len() - 1]);
    assert!(
        west > 0.0 && east < 0.0,
        "steady easterly alizés must leave the equatorial thermocline deep in the west and \
         shallow in the east, but the settled basin ends at {:.1} m against the western wall and \
         {:.1} m against the eastern one",
        west * mean_depth_m,
        east * mean_depth_m
    );

    // And monotonically so across the interior: the steady response to a
    // zonally uniform stress has no interior structure to reverse the slope,
    // so a column deeper than the one to its west would be a reflection
    // artefact rather than a tilt. The walls' own reflection layers are
    // excluded — [`EquatorialBasin::interior_columns`] says why they are not a
    // slope the continuum solution has.
    let (interior_west, interior_east) = basin.interior_columns();
    for (offset, pair) in profile[interior_west..=interior_east]
        .windows(2)
        .enumerate()
    {
        let column = interior_west + offset;
        assert!(
            pair[1] <= pair[0],
            "the equatorial thermocline deepens eastward between columns {column} and {}: \
             {:.3} m then {:.3} m, which is not a tilt",
            column + 1,
            pair[0] * mean_depth_m,
            pair[1] * mean_depth_m
        );
    }
}
