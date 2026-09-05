//! Acceptance tests for T-07.1 — the equatorial Kelvin wave's propagation
//! speed and its non-dispersion.
//!
//! `CONTEXT.md` (*Kelvin wave*) states the three claims this file validates: an
//! equatorially trapped disturbance of the first baroclinic mode travels
//! **eastward only**, at `c = √(g'·H)`, and **without dispersion** — it keeps
//! its shape as it goes. `docs/planning/01-scientific-model.md` makes the same
//! statement of the 1.5-layer core, and ADR-0003 fixes the discretisation whose
//! truncation error is the only thing standing between the two.
//!
//! Nothing below is measured out of a run and pasted back in as an
//! expectation. Every asserted number is `c` itself, zero, or a tolerance
//! assembled from named error terms of the configuration (CODING_STANDARDS.md
//! § *Tests*), and the acceptance criterion — the measured phase speed within
//! the documented tolerance of `√(g'H)` at more than one resolution, with the
//! error *shrinking* rather than sitting at a fixed offset — is asserted twice
//! over: as a point check against a per-resolution budget in
//! [`the_kelvin_pulse_travels_east_at_the_reduced_gravity_wave_speed`], and as
//! a convergence rate in
//! [`the_kelvin_speed_error_shrinks_at_the_schemes_second_order`].
//!
//! # The mode, and how a run is decomposed onto it
//!
//! Write the linear equations of [ADR-0003] in the equatorial variables:
//! `η = y/Le` with `Le = √(c/β)` the equatorial deformation radius, velocities
//! scaled by `c` and `h` by the mean thermocline depth `H`. The meridional
//! structures are the parabolic cylinder functions
//! `ψ_n(η) = H_n(η)·e^{−η²/2}` built on the Hermite polynomials —
//! `ψ_0 = e^{−η²/2}`, `ψ_2 = (4η² − 2)·e^{−η²/2}` — which are mutually
//! orthogonal on the line. In the two combinations
//!
//! ```text
//! r = u/c + h/H          q = u/c − h/H
//! ```
//!
//! the system splits into an eastward and a westward part, and the Kelvin wave
//! is the branch with `v ≡ 0` and `u/c = h/H ∝ ψ_0`:
//!
//! ```text
//! q ≡ 0,     r = 2·h/H,     ∂r/∂t + c·∂r/∂x = 0.
//! ```
//!
//! That last equation is the whole of this ticket. It is an *exact* solution of
//! the continuous equations for **any** zonal profile `r(x − c·t)`, which is
//! what makes the two claims below theorems rather than approximations: the
//! speed is `c` at every wavenumber, so the packet translates unchanged and a
//! Gaussian pulse is a legitimate initial condition rather than a
//! near-solution. Every Rossby mode `n ≥ 1` puts `ψ_{n−1}` in `q` and
//! `ψ_{n+1}` in `r` (Matsuno 1966; Gill, *Atmosphere–Ocean Dynamics*, § 11.6),
//! so the two projections
//!
//! ```text
//! eastward(x) = P₀[r]      westward(x) = P₀[q]
//! ```
//!
//! separate the Kelvin wave from everything else exactly: no Rossby mode has
//! `ψ_0` in `q` except the gravest, and none has `ψ_0` in `r` at all. This is
//! the same decomposition `western_boundary_reflection.rs` and
//! `eastern_boundary_reflection.rs` measure their reflections with; here it is
//! turned on the incident wave alone, in a basin the pulse never reaches the
//! end of.
//!
//! # What is measured
//!
//! The pulse's position is the energy-weighted zonal centroid of `P₀[r]`,
//! which for a linear wave moves at the energy-weighted mean group velocity,
//! and the speed is its displacement between two sample times divided by the
//! elapsed time. The two times are the *steps'* own times, not the requested
//! ones, so the denominator carries no rounding.
//!
//! A centroid reads a group velocity, and the criterion names a *phase* speed.
//! On this branch they are one number: `ω = c·k` is linear, so
//! `ω/k = ∂ω/∂k = c` at every wavenumber, and it is precisely the Kelvin wave's
//! non-dispersion that makes the two coincide. What separates them in a
//! *discrete* run is the truncation, and the budget below is stated at the
//! group speed's `(kΔx)²/2` rather than the phase speed's `(kΔx)²/6` — the
//! larger of the two, so the bound covers the quantity actually measured.
//!
//! The pulse's shape is the RMS zonal width of the same profile about that
//! centroid, over a window [`SHAPE_WINDOW_IN_WIDTHS`] wide either side of it.
//! A non-dispersive wave's width does not change; a dispersive one's grows.
//!
//! # Where the tolerances come from
//!
//! ## The speed: `c`, and nothing but numerical error
//!
//! Because the Kelvin branch is non-dispersive *at every wavenumber*, the
//! theoretical centroid speed is `c` exactly and there is no physical bias term
//! of the kind the reflected Rossby packet of T-04.4 carries. The budget is
//! purely the discretisation's, and both of its terms are second order in the
//! cell size, which is why the whole error converges at second order:
//!
//! | term | coarse size | why |
//! |---|---|---|
//! | meridional truncation | `(Δy/Le)² = 2.1%` | the C-grid operators of T-01.1 are second order, and `Le` is the scale of the `ψ_0` structure they differentiate; the discrete geostrophic balance that holds `v ≡ 0` is only satisfied to this order, so this much of the wave is shed into other modes |
//! | zonal group-speed truncation | `(Δx/σ)²/4 = 0.44%` | the centred difference has `ω = c·k·(1 − (kΔx)²/6)`, hence a group velocity `c·(1 − (kΔx)²/2)`; the energy spectrum of a Gaussian of width `σ` has `⟨k²⟩ = 1/(2σ²)`, giving a mean group-speed deficit of `(Δx/σ)²/4` |
//! | RK4 phase error | `1×10⁻⁸` | fourth order: `(ω·Δt)⁴/120` at the pulse's dominant frequency `ω = c/σ` and this run's CFL-stable timestep |
//! | centroid contamination | `4×10⁻⁴` | the modes shed by the first term hold `(Δy/Le)⁴ = 4×10⁻⁴` of the energy at a lever arm of a few `σ`, and the pulse's tails are clipped by the walls at 4 `σ` (`e^{−16}` in energy) |
//!
//! The last two are five and one orders of magnitude below the first two and
//! are not carried into the number. What is carried is the sum of the first
//! two, times [`TRUNCATION_SAFETY`]: neither leading coefficient is evaluated
//! (the meridional one depends on how the discrete Coriolis stencil acts on
//! `ψ_0`, the zonal one on the fourth moment of the packet's spectrum), and a
//! factor of two is the standing allowance for an unevaluated `O(1)`
//! coefficient. It multiplies both resolutions equally, so it changes what the
//! point checks allow and not what the convergence test measures.
//!
//! At the coarse resolution that is 5.1%, at the fine one 1.3% — the budget
//! itself shrinks by four when the cells are halved, so passing it at both
//! resolutions is already the "not a fixed offset" the acceptance criterion
//! asks for. Both are bounds and not estimates, so the point check is
//! generous by design; what pins the error's *size* to the discretisation
//! rather than to a coincidence is
//! [`the_kelvin_speed_error_shrinks_at_the_schemes_second_order`], which
//! measures the rate.
//!
//! ## The shape: unchanged, to two terms
//!
//! The continuous wave's RMS width is constant. Two things move the measured
//! one, and both are computed in [`Experiment::width_growth_bound`]:
//!
//! - **Numerical dispersion.** With `ω = c·k − a·k³` and `a = c·Δx²/6`, a
//!   Gaussian whose spectrum is `e^{−k²σ²/2}` acquires the variance
//!   `Var(x) = σ² + 9a²t²·Var(k²)`, and `Var(k²) = 2⟨k²⟩² = 1/(2σ⁴)` for that
//!   spectrum, so the relative width growth after a flight of `c·t` is
//!   `(c·t/σ)²·(Δx/σ)⁴/16` — `4.9×10⁻⁴` over this run.
//! - **Shed modes.** The `(Δy/Le)²` of the speed budget again, now as an
//!   amplitude: content of that size at a lever arm of the window half-width
//!   `W` moves the variance by `(Δy/Le)⁴·(W/σ)²` relative — `4.0×10⁻³` here,
//!   the larger of the two.
//!
//! Summed and multiplied by the same [`TRUNCATION_SAFETY`]: 0.9% at the coarse
//! resolution. That bound catches gross dispersion — a scheme that smeared the
//! packet would fail it — but it is not sharp enough to separate the Kelvin
//! branch from a weakly dispersive one, so it is not the whole of the
//! non-dispersion claim. The next section is.
//!
//! ## Non-dispersion, stated as an independence rather than a bound
//!
//! Width growth is one face of dispersion. The sharper one is that a
//! non-dispersive wave's speed does not depend on which wavenumbers it is made
//! of, so [`the_kelvin_speed_does_not_depend_on_the_packets_zonal_width`] runs
//! the same experiment with a pulse [`NARROW_WIDTH_FACTOR`] times narrower —
//! doubling the spectral content — and holds the two speeds to the sum of
//! their zonal budgets, 4.4%. For contrast, the same halving applied to the
//! gravest Rossby mode, whose long-wave speed carries the dispersive bias
//! `(4/9)·(Le/σ)²`, would move its measured speed by 7.1%: the assertion
//! discriminates between a non-dispersive branch and a dispersive one at this
//! configuration, which is what makes it worth making.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

mod support;

use std::sync::OnceLock;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing, WaveSpeed,
    WindStressField, H_STAGGERING, U_STAGGERING,
};

use support::{
    equatorial_deformation_radius_m, gaussian_envelope, kelvin_wave_speed_m_per_s, pacific_params,
    MeridionalStructure,
};

/// Zonal extent of the test basin, in metres — 20 000 km, the order of the
/// equatorial Pacific's width (`CONTEXT.md`, *Basin*).
///
/// The pulse starts four of its own widths from the western wall and ends the
/// run four from the eastern one, so neither boundary takes part in what is
/// measured: this is a validation of the *interior* wave, and the boundaries'
/// own physics is the subject of T-04.3 and T-04.4.
const BASIN_LX_M: f64 = 2.0e7;
/// Meridional extent of the test basin, in metres: ±2000 km, about ±5.8 `Le`.
///
/// Far enough out that `ψ_0` is `e^{−16.8}` at the northern and southern walls,
/// so the equatorial waveguide does not feel them, and near enough that the
/// gravity-wave CFL bound stays the binding one rather than ADR-0007's rotation
/// bound.
const BASIN_LY_M: f64 = 4.0e6;
/// Cell width of the coarse run, in metres.
///
/// About seven to a pulse width, which is what the `(Δx/σ)²/4` entry of the
/// speed budget costs.
const COARSE_CELL_WIDTH_M: f64 = 2.0e5;
/// Cell height of the coarse run, in metres.
///
/// About seven to an equatorial deformation radius, which is what the
/// `(Δy/Le)²` entry costs — the dominant term of every budget in this file.
const COARSE_CELL_HEIGHT_M: f64 = 5.0e4;
/// Factor by which the coarse run refines the two cell sizes above: one, it
/// being the run they state. Named so that neither resolution is a bare
/// literal where a run is asked for.
const COARSE_REFINEMENT: usize = 1;
/// Factor by which the fine run refines both cell dimensions.
///
/// Both axes together, so that both second-order terms of the speed budget
/// shrink by the same four and the measured convergence rate is the scheme's
/// rather than a mixture of two rates.
const FINE_REFINEMENT: usize = 2;

/// Zonal width `σ` of the reference pulse's Gaussian envelope, in metres —
/// 4.3 `Le`.
///
/// The Kelvin branch is an exact solution at every wavenumber, so no long-wave
/// approximation constrains this from below; what does is the zonal truncation
/// term, which grows as `σ⁻²`. What constrains it from above is the basin: the
/// pulse needs four widths of clearance at each wall and five widths of flight
/// between them.
const PULSE_WIDTH_M: f64 = 1.5e6;
/// Factor by which the second, narrower pulse of the non-dispersion test is
/// compressed. Two, so its spectral content is doubled.
const NARROW_WIDTH_FACTOR: f64 = 2.0;
/// Zonal position of every pulse's centre at `t = 0`, in metres east of the
/// western wall.
///
/// Four reference widths, so the wall sees `e^{−8} = 3×10⁻⁴` of the pulse's
/// amplitude and the run starts with an undisturbed boundary. The narrower
/// pulse starts at the same place, so both runs measure the same flight.
const PULSE_CENTRE_X_M: f64 = 6.0e6;
/// Peak thermocline depth anomaly of the pulse, in metres — a downwelling
/// Kelvin pulse of the scale a westerly wind burst leaves behind (`CONTEXT.md`,
/// *Westerly wind burst*).
///
/// The core is linear, so every speed, width and ratio below is independent of
/// this number; it only sets the scale the diagnostics are reported in.
const PULSE_AMPLITUDE_M: f64 = 10.0;

/// When the pulse is first sampled, in transits of one reference width at `c`.
///
/// Half a width in: far enough that the transient the continuous initial
/// condition sheds on the discrete grid has separated from the packet, near
/// enough to leave a long baseline for the speed.
const SAMPLE_EARLY_IN_TRANSITS: f64 = 0.5;
/// When it is sampled again, in the same transits. The centre has then
/// travelled five widths and is four short of the eastern wall.
const SAMPLE_LATE_IN_TRANSITS: f64 = 5.0;

/// Zonal half-width of the window the pulse's shape is read in, in reference
/// widths. Three `σ` either side of the centroid holds all but `e^{−9}` of it.
const SHAPE_WINDOW_IN_WIDTHS: f64 = 3.0;

/// Factor applied to every truncation-derived bound in this file, for the
/// leading `O(1)` coefficients the truncation terms do not evaluate.
///
/// Two. It multiplies the coarse and the fine budgets alike, so it widens what
/// the point checks admit without touching the rate
/// [`the_kelvin_speed_error_shrinks_at_the_schemes_second_order`] measures.
const TRUNCATION_SAFETY: f64 = 2.0;

/// Smallest convergence order the Kelvin speed's error must show when both cell
/// dimensions are halved.
///
/// The spatial discretisation is second order (ADR-0003) and both terms of the
/// speed budget are second order in the cell size, so the measured order
/// `log₂(coarse error / fine error)` should be 2. Requiring 1.5 leaves margin
/// for the sub-dominant terms of the table above, which do not scale that way,
/// while still failing a first-order scheme — which is the point of asserting
/// an order rather than a bare shrinkage.
const MIN_CONVERGENCE_ORDER: f64 = 1.5;

/// Largest convergence order that same refinement may show.
///
/// A second-order scheme cannot converge faster than second order, so an order
/// well above 2 is not a better result: it means the fine run's error is not
/// the truncation at all — two terms of the budget cancelling by accident, or
/// the measurement's own floor — and the ratio the test reads is then a
/// coincidence rather than a rate. Three is 2 plus the same half-order of slack
/// [`MIN_CONVERGENCE_ORDER`] leaves on the other side.
const MAX_CONVERGENCE_ORDER: f64 = 3.0;

/// One propagation experiment: a basin at some refinement, carrying a Gaussian
/// Kelvin pulse of some zonal width.
///
/// Refinement and width are the two things the runs of this file differ in —
/// the first for the convergence of the speed, the second for its independence
/// of the spectrum — and everything derived from them lives here rather than
/// being written down once per run.
#[derive(Debug, Clone, Copy)]
struct Experiment {
    /// Shape, spacing and position of the basin.
    basin: Basin,
    /// The ocean the equations are written in terms of.
    params: PhysicalParams,
    /// `Le = √(c/β)`, in metres — cached because every projection needs it.
    deformation_radius_m: f64,
    /// Zonal width `σ` of this run's pulse, in metres.
    pulse_width_m: f64,
}

impl Experiment {
    /// The experiment at `1/refinement` of the coarse cell size, carrying a
    /// pulse `pulse_width_m` wide.
    fn new(refinement: usize, pulse_width_m: f64) -> Self {
        let cell_width_m = COARSE_CELL_WIDTH_M / refinement as f64;
        let cell_height_m = COARSE_CELL_HEIGHT_M / refinement as f64;
        let grid = Grid::new(
            (BASIN_LX_M / cell_width_m).round() as usize,
            (BASIN_LY_M / cell_height_m).round() as usize,
        )
        .expect("the basin has cells on both axes");
        let spacing = Spacing::new(cell_width_m, cell_height_m)
            .expect("the cell sizes are finite and positive");
        let params = pacific_params();
        Self {
            basin: Basin::centered_on_equator(grid, spacing),
            params,
            deformation_radius_m: equatorial_deformation_radius_m(),
            pulse_width_m,
        }
    }

    /// The Kelvin wave speed `c = √(g'·H)`, in m/s: the analytic
    /// [`kelvin_wave_speed_m_per_s`], never the engine's own.
    fn wave_speed_m_per_s(self) -> f64 {
        kelvin_wave_speed_m_per_s()
    }

    /// The two sample times of every run, in seconds.
    ///
    /// Stated in transits of the *reference* width so that the narrow run
    /// samples the same clock, and therefore the same flight, as the wide one:
    /// the two speeds are then comparable without a correction.
    fn sample_times_s(self) -> [f64; 2] {
        let transit_s = PULSE_WIDTH_M / self.wave_speed_m_per_s();
        [SAMPLE_EARLY_IN_TRANSITS, SAMPLE_LATE_IN_TRANSITS].map(|transits| transits * transit_s)
    }

    /// The tolerance on this run's measured speed, as a fraction of `c`.
    ///
    /// The two second-order truncation terms of the module header, times
    /// [`TRUNCATION_SAFETY`]. It is a function of the run rather than a
    /// constant because the point of the acceptance criterion is that a finer
    /// grid is held to a tighter bound.
    fn speed_tolerance(self) -> f64 {
        TRUNCATION_SAFETY * self.meridional_truncation() + self.zonal_speed_budget()
    }

    /// This run's share of a speed budget from the zonal group-speed
    /// truncation alone, as a fraction of `c`: `(Δx/σ)²/4`, times
    /// [`TRUNCATION_SAFETY`].
    ///
    /// Named rather than inlined because it is claimed twice — once inside
    /// [`Experiment::speed_tolerance`], and once as the whole budget of
    /// [`the_kelvin_speed_does_not_depend_on_the_packets_zonal_width`], where
    /// the meridional term is common to the two runs being compared.
    fn zonal_speed_budget(self) -> f64 {
        TRUNCATION_SAFETY * 0.25 * self.zonal_truncation_in_widths()
    }

    /// `(Δy/Le)²`: the second-order meridional truncation of the waveguide, as
    /// a fraction.
    fn meridional_truncation(self) -> f64 {
        let cell_in_radii = self.basin.spacing().dy_m() / self.deformation_radius_m;
        cell_in_radii * cell_in_radii
    }

    /// `(Δx/σ)²`: the square of the zonal cell size in pulse widths, which the
    /// group-speed and the shape budgets both build on.
    fn zonal_truncation_in_widths(self) -> f64 {
        let cell_in_widths = self.basin.spacing().dx_m() / self.pulse_width_m;
        cell_in_widths * cell_in_widths
    }

    /// The largest relative change in the pulse's RMS width the run's own
    /// numerics account for — the two terms derived in the module header,
    /// times [`TRUNCATION_SAFETY`].
    fn width_growth_bound(self) -> f64 {
        let flight_in_widths =
            self.wave_speed_m_per_s() * self.sample_times_s()[1] / self.pulse_width_m;
        let dispersion = flight_in_widths
            * flight_in_widths
            * self.zonal_truncation_in_widths()
            * self.zonal_truncation_in_widths()
            / 16.0;
        let shed_modes = self.meridional_truncation()
            * self.meridional_truncation()
            * SHAPE_WINDOW_IN_WIDTHS
            * SHAPE_WINDOW_IN_WIDTHS;
        TRUNCATION_SAFETY * (dispersion + shed_modes)
    }

    /// The thermocline depth anomaly's `ψ_n` coefficient, column by column, in
    /// units of the mean depth `H`.
    fn depth_projection(self, state: &OceanState, structure: MeridionalStructure) -> Vec<f64> {
        support::depth_projection(
            self.basin,
            self.deformation_radius_m,
            self.params,
            state,
            structure,
        )
    }

    /// The two `ψ_0` invariants of `state`, column by column, and the `ψ_0`
    /// depth projection they are built from.
    fn invariants(self, state: &OceanState) -> support::Invariants {
        support::invariants(
            self.basin,
            self.deformation_radius_m,
            self.params,
            state,
            self.wave_speed_m_per_s(),
        )
    }

    /// Zonal position of column `i`'s centre, in metres east of the western
    /// wall.
    fn column_x_m(self, i: usize) -> f64 {
        self.basin.x_of_column_m(H_STAGGERING, i)
    }

    /// `Σ profile²` over the whole basin — the weight the centroid below
    /// divides by, and the quantity the eastward/westward split is stated in.
    fn energy(self, profile: &[f64]) -> f64 {
        support::energy(profile.iter().copied())
    }

    /// The energy-weighted zonal centroid of `profile`, in metres.
    ///
    /// Over the whole basin, with no window: a window would need a position to
    /// be centred on, and taking that from theory is exactly the circularity a
    /// speed measurement must not have. The cost is that whatever the run sheds
    /// is weighed too, which the module header bounds at `4×10⁻⁴` of the
    /// centroid — three orders below the terms the tolerance is built from.
    ///
    /// # Panics
    /// If the profile carries no energy at all, which would mean the run never
    /// had a wave in it.
    fn energy_centroid_m(self, profile: &[f64]) -> f64 {
        support::energy_centroid_m(
            profile
                .iter()
                .enumerate()
                .map(|(i, amplitude)| (self.column_x_m(i), *amplitude)),
        )
    }

    /// The RMS zonal width of `profile` about its own centroid, in metres, over
    /// the columns within [`SHAPE_WINDOW_IN_WIDTHS`] of it.
    ///
    /// The window travels with the packet, so the same fraction of the same
    /// wave is weighed at both sample times and the comparison between them is
    /// a statement about the packet's shape rather than about where it is.
    ///
    /// It is an *energy*-weighted RMS, so for a Gaussian amplitude profile of
    /// width `σ` it reads `σ/√2` rather than `σ`; only the ratio of two of them
    /// is asserted, so that factor cancels. The window reaches
    /// [`SHAPE_WINDOW_IN_WIDTHS`]`·√2 = 4.2` of those RMS widths, which is what
    /// keeps it from clipping the growth it exists to measure: a packet would
    /// have to spread by many times the bound before the window, rather than
    /// the dispersion, set the number.
    fn rms_width_m(self, profile: &[f64]) -> f64 {
        support::rms_width_m(
            profile,
            |i| self.column_x_m(i),
            SHAPE_WINDOW_IN_WIDTHS * self.pulse_width_m,
        )
    }

    /// The initial condition: an equatorial Kelvin pulse, Gaussian in `x` on
    /// the `ψ_0` waveguide, with `u = (c/H)·h` and `v = 0`.
    ///
    /// An exact solution of the continuous equations for any zonal profile
    /// (module header), so the run starts with one wave in it, travelling east,
    /// and no Rossby energy at all — which is what the westward projection of
    /// [`the_kelvin_pulse_carries_no_westward_energy`] checks has stayed true.
    fn initial_state(self) -> OceanState {
        let mut state = OceanState::at_rest(self.basin.grid());
        let current_amplitude_m_per_s =
            PULSE_AMPLITUDE_M * self.wave_speed_m_per_s() / self.params.mean_thermocline_depth_m();
        let profile = |x_m: f64| gaussian_envelope(x_m, PULSE_CENTRE_X_M, self.pulse_width_m);

        for j in 0..state.h().ny() {
            let waveguide = MeridionalStructure::Gravest.at(
                self.basin.y_of_row_m(H_STAGGERING, j),
                self.deformation_radius_m,
            );
            for i in 0..state.h().nx() {
                let x_m = self.basin.x_of_column_m(H_STAGGERING, i);
                *state.h_mut().get_mut(i, j).expect("a cell centre") =
                    PULSE_AMPLITUDE_M * profile(x_m) * waveguide;
            }
            for i in 0..state.u().nx() {
                let x_m = self.basin.x_of_column_m(U_STAGGERING, i);
                *state.u_mut().get_mut(i, j).expect("an east/west face") =
                    current_amplitude_m_per_s * profile(x_m) * waveguide;
            }
        }
        state
    }

    /// Run the experiment: the pulse in a closed, unforced, undamped basin,
    /// sampled twice while it is in open water.
    fn run(self) -> Propagation {
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
                    eastward: invariants.eastward,
                    westward: invariants.westward,
                    depth_in_second: self.depth_projection(&state, MeridionalStructure::Second),
                    depth_in_gravest: invariants.depth_in_gravest,
                });
            }
            if step < last_step {
                solver.step(&mut state, step as f64 * dt_s, |_t| &calm);
            }
        }

        let mut samples = samples.into_iter();
        let mut next = || samples.next().expect("every requested sample was taken");
        Propagation {
            samples: [next(), next()],
        }
    }

    /// The measured zonal speed of the pulse, in m/s — positive eastward.
    fn measured_speed_m_per_s(self, propagation: &Propagation) -> f64 {
        let [early, late] = &propagation.samples;
        (self.energy_centroid_m(&late.eastward) - self.energy_centroid_m(&early.eastward))
            / (late.t_s - early.t_s)
    }

    /// How far that speed sits from `c`, as a fraction of `c`.
    fn speed_error(self, propagation: &Propagation) -> f64 {
        (self.measured_speed_m_per_s(propagation) - self.wave_speed_m_per_s()).abs()
            / self.wave_speed_m_per_s()
    }
}

/// One sample of a run: the two `ψ_0` invariants and the thermocline anomaly's
/// two meridional coefficients, column by column, at a given time.
#[derive(Debug, Clone)]
struct Sample {
    /// When it was taken, in seconds — the step's own time, not the requested
    /// one, so that a speed's denominator is exact.
    t_s: f64,
    /// `P₀[u/c + h/H]`, column by column: where the Kelvin wave is.
    eastward: Vec<f64>,
    /// `P₀[u/c − h/H]`, column by column: what, if anything, is going west.
    westward: Vec<f64>,
    /// The `ψ_0` coefficient of `h/H`, column by column.
    depth_in_gravest: Vec<f64>,
    /// The `ψ_2` coefficient of `h/H`, column by column.
    depth_in_second: Vec<f64>,
}

/// What one run leaves for the tests to read: the pulse early in its flight and
/// late in it.
#[derive(Debug, Clone)]
struct Propagation {
    /// The two samples, in time order.
    samples: [Sample; 2],
}

/// The coarse run at the reference width, integrated once and shared.
///
/// Several tests read several different things out of one trajectory, which is
/// both cheaper than a run each and stronger: every assertion is then made
/// about the same wave.
fn coarse_run() -> &'static Propagation {
    static RUN: OnceLock<Propagation> = OnceLock::new();
    RUN.get_or_init(|| Experiment::new(COARSE_REFINEMENT, PULSE_WIDTH_M).run())
}

/// The same experiment with both cell dimensions refined by
/// [`FINE_REFINEMENT`] — the second resolution the acceptance criterion asks
/// for.
fn fine_run() -> &'static Propagation {
    static RUN: OnceLock<Propagation> = OnceLock::new();
    RUN.get_or_init(|| Experiment::new(FINE_REFINEMENT, PULSE_WIDTH_M).run())
}

/// The coarse experiment with a pulse [`NARROW_WIDTH_FACTOR`] times narrower —
/// the second spectrum the non-dispersion test compares.
fn narrow_run() -> &'static Propagation {
    static RUN: OnceLock<Propagation> = OnceLock::new();
    RUN.get_or_init(|| {
        Experiment::new(COARSE_REFINEMENT, PULSE_WIDTH_M / NARROW_WIDTH_FACTOR).run()
    })
}

/// The two resolutions of the acceptance criterion, coarse first.
fn resolutions() -> [(Experiment, &'static Propagation); 2] {
    [
        (
            Experiment::new(COARSE_REFINEMENT, PULSE_WIDTH_M),
            coarse_run(),
        ),
        (Experiment::new(FINE_REFINEMENT, PULSE_WIDTH_M), fine_run()),
    ]
}

#[test]
fn the_kelvin_pulse_travels_east_at_the_reduced_gravity_wave_speed() {
    // The acceptance criterion of T-07.1: the measured phase speed is within
    // the documented tolerance of `c = √(g'·H)` (`CONTEXT.md`, *Kelvin wave*)
    // at more than one grid resolution — and each resolution is held to its own
    // budget, so the finer run has to be better rather than merely no worse.
    for (experiment, propagation) in resolutions() {
        let expected_m_per_s = experiment.wave_speed_m_per_s();
        let measured_m_per_s = experiment.measured_speed_m_per_s(propagation);
        let cells = experiment.basin.grid();

        assert!(
            measured_m_per_s > 0.0,
            "on the {}×{} grid the pulse travelled at {measured_m_per_s} m/s, which is westward: \
             an equatorial Kelvin wave goes east only",
            cells.nx(),
            cells.ny()
        );

        let error = experiment.speed_error(propagation);
        let tolerance = experiment.speed_tolerance();
        assert!(
            error <= tolerance,
            "on the {}×{} grid the pulse travelled at {measured_m_per_s} m/s, {:.2}% from the \
             Kelvin wave speed {expected_m_per_s} m/s; that grid's truncation budget allows {:.2}%",
            cells.nx(),
            cells.ny(),
            100.0 * error,
            100.0 * tolerance
        );
    }
}

#[test]
fn the_kelvin_speed_error_shrinks_at_the_schemes_second_order() {
    // The other half of the acceptance criterion: "demonstrating the error
    // shrinks with resolution, not a fixed offset". Both terms of the speed
    // budget are second order in the cell size (module header) and both cell
    // dimensions are halved, so the error must fall by about four
    // (CODING_STANDARDS.md § *Convergence over point checks*).
    let [(coarse, coarse_run), (fine, fine_run)] = resolutions();

    let coarse_error = coarse.speed_error(coarse_run);
    let fine_error = fine.speed_error(fine_run);
    assert!(
        fine_error > 0.0,
        "the fine run reproduced c to the last bit, so there is no error left to measure a rate \
         on: what this run reads is the measurement's floor and not the scheme's order"
    );
    let order = (coarse_error / fine_error).log2();

    // Bounded on both sides: too small an order is a scheme that is not second
    // order, and too large a one is a fine error that is no longer the
    // truncation — either way the point check's budget has stopped describing
    // what the run does.
    assert!(
        (MIN_CONVERGENCE_ORDER..=MAX_CONVERGENCE_ORDER).contains(&order),
        "refining the cells by {FINE_REFINEMENT} took the Kelvin speed's error from {:.3}% to \
         {:.3}%, a convergence order of {order:.2}: outside the \
         [{MIN_CONVERGENCE_ORDER}, {MAX_CONVERGENCE_ORDER}] a second-order scheme owes",
        100.0 * coarse_error,
        100.0 * fine_error
    );
}

#[test]
fn the_kelvin_pulse_carries_no_westward_energy() {
    // "Eastward-only" (`CONTEXT.md`, *Kelvin wave*), read through the
    // decomposition that can tell the two directions apart: the Kelvin branch
    // has `q ≡ 0` exactly, so `P₀[q]` is empty for the whole flight and
    // anything in it is either a Rossby wave the run should not contain or the
    // discretisation's leakage.
    //
    // The initial condition sets `u/c = h/H`, so the early sample starts from
    // `q ≡ 0` by construction and it is the *late* one that carries the claim:
    // a wave that had turned round, or split, would have put energy there over
    // the flight. That the wave moved at all is what the speed test asserts,
    // and this one leans on that rather than repeating it.
    let (experiment, propagation) = resolutions()[0];

    // The `(Δy/Le)²` of the speed budget is an amplitude, so it is that
    // squared in energy, times `TRUNCATION_SAFETY` for the coefficient it does
    // not evaluate.
    let leakage = experiment.meridional_truncation();
    let ceiling = TRUNCATION_SAFETY * leakage * leakage;

    for sample in &propagation.samples {
        let westward = experiment.energy(&sample.westward);
        let eastward = experiment.energy(&sample.eastward);
        let share = westward / (eastward + westward);
        assert!(
            share <= ceiling,
            "{:.0} s into the run the westward mode held {:.4}% of the energy, past the {:.4}% \
             the meridional discretisation can leak: the pulse is not travelling east only",
            sample.t_s,
            100.0 * share,
            100.0 * ceiling
        );
    }
}

#[test]
fn the_kelvin_pulse_keeps_its_zonal_shape() {
    // Non-dispersion, in the form the deliverable states it: the packet's shape
    // does not change as it travels. `∂r/∂t + c·∂r/∂x = 0` translates any
    // profile unchanged, so the RMS zonal width of `P₀[r]` at the late sample
    // is the width at the early one, up to the two numerical terms
    // [`Experiment::width_growth_bound`] adds up.
    let (experiment, propagation) = resolutions()[0];
    let [early, late] = &propagation.samples;

    let early_width_m = experiment.rms_width_m(&early.eastward);
    let late_width_m = experiment.rms_width_m(&late.eastward);
    let growth = (late_width_m - early_width_m).abs() / early_width_m;
    let bound = experiment.width_growth_bound();

    assert!(
        growth <= bound,
        "over its flight the pulse's RMS zonal width went from {early_width_m} m to \
         {late_width_m} m, a change of {:.2}%; a non-dispersive wave keeps its shape, and this \
         run's numerics account for {:.2}%",
        100.0 * growth,
        100.0 * bound
    );
}

#[test]
fn the_kelvin_speed_does_not_depend_on_the_packets_zonal_width() {
    // The sharp statement of non-dispersion: every wavenumber travels at `c`,
    // so a packet with twice the spectral content travels at the same speed.
    // The two runs share a grid, a start position and a clock, so what differs
    // between their measured speeds is the spectrum alone — and the only budget
    // that is not common to them is the zonal group-speed truncation, which is
    // four times larger for the narrower pulse.
    let wide = Experiment::new(COARSE_REFINEMENT, PULSE_WIDTH_M);
    let narrow = Experiment::new(COARSE_REFINEMENT, PULSE_WIDTH_M / NARROW_WIDTH_FACTOR);

    let wide_m_per_s = wide.measured_speed_m_per_s(coarse_run());
    let narrow_m_per_s = narrow.measured_speed_m_per_s(narrow_run());

    let difference = (narrow_m_per_s - wide_m_per_s).abs() / wide.wave_speed_m_per_s();
    let bound = wide.zonal_speed_budget() + narrow.zonal_speed_budget();

    assert!(
        difference <= bound,
        "narrowing the pulse by {NARROW_WIDTH_FACTOR} moved the measured speed from \
         {wide_m_per_s} m/s to {narrow_m_per_s} m/s, {:.2}% of c; a non-dispersive wave's speed \
         does not depend on its spectrum, and the two zonal truncations allow {:.2}%",
        100.0 * difference,
        100.0 * bound
    );
}

#[test]
fn the_kelvin_pulse_stays_on_the_equatorial_waveguide() {
    // The meridional half of "unchanging shape": the Kelvin wave's thermocline
    // anomaly is `ψ_0` and nothing else — a single lobe centred on the equator
    // — so its `ψ_2` content stays at the zero it started from, to the order
    // the meridional discretisation can produce.
    let (experiment, propagation) = resolutions()[0];
    let late = &propagation.samples[1];

    let (cross_term, gravest_energy) = late
        .depth_in_gravest
        .iter()
        .zip(&late.depth_in_second)
        .fold((0.0, 0.0), |(cross, gravest), (in_gravest, in_second)| {
            (
                cross + in_gravest * in_second,
                gravest + in_gravest * in_gravest,
            )
        });
    let ratio = (cross_term / gravest_energy).abs();
    let ceiling = TRUNCATION_SAFETY * experiment.meridional_truncation();

    assert!(
        ratio <= ceiling,
        "at the end of its flight the pulse's thermocline anomaly carried a ψ₂/ψ₀ ratio of \
         {ratio}, past the {ceiling} the meridional truncation allows: a Kelvin wave is ψ₀ alone",
    );
}
