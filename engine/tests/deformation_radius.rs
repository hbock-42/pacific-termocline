//! Acceptance tests for T-07.3 — the meridional decay scale of the equatorial
//! waves is the deformation radius `Le = √(c/β)`.
//!
//! `CONTEXT.md` (*Equatorial deformation radius*) states the claim this file
//! validates: `Le = √(c/β)` is "the meridional scale over which equatorial
//! waves decay away from the equator". T-07.1 and T-07.2 measured how fast the
//! two waves travel; this one measures how wide they are, for **both** wave
//! types, by fitting the meridional profile a run actually carries and
//! comparing the fitted scale against the analytic `Le`.
//!
//! Nothing below is measured out of a run and pasted back in as an
//! expectation. Every asserted number is `Le` itself or a tolerance assembled
//! from named error terms of the configuration (CODING_STANDARDS.md § *Tests*).
//!
//! # What is fitted, and to what
//!
//! Write the linear equatorial equations in the equatorial variables of
//! [ADR-0003] — `ŷ = y/Le`, velocities scaled by `c`, `h` by the mean
//! thermocline depth `H` — and they separate on the parabolic cylinder
//! functions `ψₘ(ŷ) = Hₘ(ŷ)·e^{−ŷ²/2}`. In the two invariants
//!
//! ```text
//! eastward = u/c + h/H          westward = u/c − h/H
//! ```
//!
//! each wave is one `ψₘ` of one invariant (Matsuno 1966; Gill,
//! *Atmosphere–Ocean Dynamics*, § 11.6), and that is what makes a *shape* fit
//! possible at all — a raw field mixes the two branches, and a mixture of two
//! Gaussians of the same scale is not a Gaussian of that scale:
//!
//! | wave | invariant | shape | why |
//! |---|---|---|---|
//! | Kelvin | eastward | `ψ₀(ŷ)` | `u/c = h/H ∝ ψ₀`, so `r = 2h/H` and `q ≡ 0` |
//! | gravest Rossby | westward | `ψ₀(ŷ)` | `h/H ∝ (2ŷ²+1)e^{−ŷ²/2}`, `u/c ∝ (2ŷ²−3)e^{−ŷ²/2}`, so `q = −4·ψ₀` |
//! | gravest Rossby | eastward | `ψ₂(ŷ)` | the same two combine to `r = (4ŷ²−2)e^{−ŷ²/2}` |
//!
//! Every one of those shapes is stretched by `Le` and by nothing else, which is
//! the statement being tested. The measurement is
//! [`support::fitted_trapping_scale_m`]: the zonal sum of the invariant, row by
//! row, fitted against `ψₘ(y/L)` with the amplitude eliminated analytically and
//! the scale `L` left free. `L` is searched for over a bracket spanning
//! [`SCALE_BRACKET_IN_RADII`] — a factor of nine, wide enough that the answer
//! is found rather than assumed, and stated in radii only because the search
//! has to start somewhere.
//!
//! The third row is worth its own test rather than being folded into the
//! second. `ψ₀` and `ψ₂` decay on the same `Le` but do not look alike — `ψ₂`
//! has two off-equatorial lobes and a node at `ŷ = ±1/√2` — so a fit that
//! recovered `Le` from both is a statement about the waveguide's scale and not
//! about one convenient Gaussian.
//!
//! The profile is read once, at the end of a flight of several packet widths,
//! not at `t = 0`. At `t = 0` the state *is* the analytic profile and a fit of
//! it would measure only the fitting code. What the tests below measure is that
//! the run **stays** on the waveguide: `Le` is not imposed anywhere in the
//! solver — it emerges from `β` and the pressure gradient — so a wave that the
//! discrete equations trapped on the wrong scale, or failed to trap, would have
//! spread or narrowed by the time it is read.
//!
//! # Three oceans, so that `Le` is a prediction rather than a length
//!
//! One ocean cannot tell `√(c/β)` apart from any other number that happens to
//! be 345 km. [`the_decay_scale_follows_le_across_oceans`] therefore fits the
//! same Kelvin pulse in three: the equatorial Pacific, one with `g'` quadrupled
//! (`c` doubled, `Le` up by `√2`) and one with `β` doubled (`Le` down by `√2`).
//! Nothing else changes — not the basin, not the packet, not the measurement.
//! The three predicted radii are 244 km, 345 km and 488 km, each 41% from its
//! neighbour, and every budget below is under 10%, so the test asserts each fit
//! is nearer its own ocean's `Le` than either neighbour's: a fit that returned a
//! fixed length would fail two of the three.
//!
//! # Where the tolerances come from
//!
//! Every entry is a property of the configuration below and none was obtained by
//! running the model. What moves a *fitted scale* is anything that perturbs the
//! profile's shape: if the profile is `ψₘ(ŷ) + ε·ψₘ₊₂(ŷ)`, then because
//! `∂/∂L ψₘ(y/L)` is itself a combination of `ψₘ₋₂` and `ψₘ₊₂`, a contamination
//! of relative amplitude `ε` moves the best-fit scale by `O(ε)`. So the budget
//! is a budget on contamination amplitudes, entry for entry:
//!
//! | term | coarse size (Kelvin / Rossby) | why |
//! |---|---|---|
//! | meridional truncation | `(Δy/Le)² = 2.1%` / `5·(Δy/Le)² = 2.6%` | the C-grid operators of T-01.1 are second order, and the scale they differentiate the wave on is `Le/√(2m+1)` for the finest `ψₘ` it carries — [`MeridionalStructure::truncation_richness`], `ψ₀` for the Kelvin wave and `ψ₂` for the Rossby mode. This much of the wave is shed into neighbouring structures over the flight |
//! | zonal truncation | `(Δx/σ)²/4 = 0.44%` / `0.20%` | the centred difference sees `k_eff = sin(kΔx)/Δx`, so wavenumbers are misread by `(kΔx)²/2`; with `⟨k²⟩ = 1/(2σ²)` for a Gaussian of width `σ` that is `(Δx/σ)²/4`. It reaches the meridional shape because a wavenumber-dependent error makes the packet's shape vary with `x`, and the zonal sum then reads a mixture rather than one profile |
//! | stray energy | — / `⟨k̂²⟩ = 0.76%` | the Rossby initial condition is the **long-wave** mode, exact only as `k̂ → 0` and correct to `O(k̂²)` in amplitude; the Kelvin branch has no such term, being an exact solution at every wavenumber |
//! | wall clipping | `ψ₀(±4.1) = 2×10⁻⁴` | the meridional walls stand at `±5.8 Le`, where `ψ₀` is `5×10⁻⁸`, and at `±4.1 Le` in the widest-waveguide ocean of the three; the profile the fit reads is the whole profile to those digits |
//! | fit quadrature | `< 10⁻¹⁰⁰` | the fit's sums are the midpoint rule on an analytic function decaying on the whole line, where Euler–Maclaurin gives an error `~e^{−2π²(Le/Δy)²}` rather than the `O(Δy²)` of a finite interval |
//! | RK4 time truncation | `10⁻⁸` | fourth order: `(ω·Δt)⁴/120` at the packet's dominant frequency and this run's CFL-stable timestep |
//!
//! The last three are two, six and many orders of magnitude below the first two
//! and are not carried into the number. What is carried is the sum of the terms
//! that are, times [`TRUNCATION_SAFETY`]: none of the leading coefficients is
//! evaluated, and a factor of two is the standing allowance for an unevaluated
//! `O(1)` coefficient. It multiplies both resolutions equally, so it changes
//! what the point checks admit and not what the convergence tests measure.
//!
//! At the coarse resolution that is 5.1% for the Kelvin fit and 7.2% for the
//! Rossby one, and 1.3% and 2.9% at the fine one — the two truncation terms
//! shrink by four when the cells are halved, the stray-energy term being the
//! one that does not, which is why the Rossby budget falls by less than four.
//! Passing a budget that shrinks with the cells at both resolutions is already
//! the "not a fixed offset" of CODING_STANDARDS.md § *Convergence over point
//! checks*. Both are bounds and not estimates, so the point checks
//! are generous by design; what pins the *physics* is
//! [`the_decay_scale_follows_le_across_oceans`], whose separation is six times
//! the budget it is asserted against, and what pins the error's *size* to the
//! discretisation is the pair of rate tests below.
//!
//! # Whether the profile is that shape at all, and how fast that improves
//!
//! A best-fit scale exists for any profile whatever, so every point check also
//! asserts the fit is a *good* one. The quantity it asserts on is `ε` itself:
//! the profile's departure from `ψₘ` at its own best-fitting scale, as a
//! relative amplitude, which is `√(1 − ρ)` for the correlation `ρ` of
//! [`support::shape_correlation`]. That is the very quantity the budget above is
//! a budget on, so it is held to the same number rather than to a second
//! tolerance.
//!
//! It is also the quantity the sharper of the two rate tests is read on.
//! [`the_meridional_shape_error_converges_at_the_schemes_second_order`] asserts
//! `ε` falls by four when both cell dimensions are halved — which is the
//! assertion that the budget's leading entry describes the scheme rather than
//! merely bounding it — and
//! [`the_fitted_decay_scale_converges_at_the_schemes_second_order`] asserts the
//! same of the fitted scale itself, the quantity the acceptance criterion names.
//! Both are made for both waves and for both of the structures the Rossby mode
//! is fitted on. The shape error is the sharper of the two because a scale
//! responds only to the part of a departure that is not orthogonal to it, while
//! the shape error sees the whole of it.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

mod support;

use std::sync::OnceLock;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing, WaveSpeed,
    WindStressField,
};

use support::{
    deformation_radius_of_m, params_with_reduced_gravity_and_beta, wave_speed_of_m_per_s,
    Invariant, MeridionalStructure, Packet, BETA_PER_M_PER_S, PACIFIC_REDUCED_GRAVITY_M_PER_S2,
};

/// Meridional extent of every test basin, in metres: ±2000 km.
///
/// 5.8 `Le` in the equatorial Pacific and 4.1 `Le` in the fastest of the three
/// oceans, so the profile the fit reads is complete to `ψ₀(±4.1) = 2×10⁻⁴` at
/// worst — the wall-clipping entry of the budget. Near enough that the
/// gravity-wave CFL bound stays the binding one rather than ADR-0007's rotation
/// bound.
const BASIN_LY_M: f64 = 4.0e6;
/// Factor by which the coarse run refines the cell sizes below: one, it being
/// the run they state. Named so that neither resolution is a bare literal where
/// a run is asked for.
const COARSE_REFINEMENT: usize = 1;
/// Factor by which the fine run refines both cell dimensions.
///
/// Both axes together, so that both truncation terms of the budget shrink by
/// the same four and the measured convergence rate is the scheme's rather than
/// a mixture of two rates.
const FINE_REFINEMENT: usize = 2;

/// Peak thermocline depth anomaly of every packet, in metres — the scale a
/// westerly wind burst leaves behind (`CONTEXT.md`, *Westerly wind burst*).
///
/// The core is linear and a fitted *scale* is invariant under a change of
/// amplitude besides, so this number appears in no assertion; it only sets the
/// units the diagnostics are reported in.
const PACKET_AMPLITUDE_M: f64 = 10.0;

/// Factor applied to every truncation-derived bound in this file, for the
/// leading `O(1)` coefficients the truncation terms do not evaluate.
///
/// Two. It multiplies the coarse and the fine budgets alike, so it widens what
/// the point checks admit without touching the rate
/// [`the_fitted_decay_scale_converges_at_the_schemes_second_order`] measures.
const TRUNCATION_SAFETY: f64 = 2.0;

/// The range of trapping scales the fit searches, as multiples of the ocean's
/// own `Le`.
///
/// A factor of nine from end to end. Wide enough that the fit finds the scale
/// rather than being told it: the three oceans' radii span a factor of two
/// between them and all three fall well inside a single bracket of this width,
/// so a run whose wave was trapped on a neighbouring ocean's radius — or on
/// half or twice its own — would be reported as such rather than clipped to the
/// edge. Narrow enough that the coarse scan of the fit resolves the objective's
/// humps.
const SCALE_BRACKET_IN_RADII: (f64, f64) = (1.0 / 3.0, 3.0);

/// Smallest convergence order the fitted scale's error must show when both cell
/// dimensions are halved.
///
/// The spatial discretisation is second order (ADR-0003) and both terms of the
/// budget are second order in the cell size, so the measured order
/// `log₂(coarse error / fine error)` should be 2. Requiring 1.5 leaves margin
/// for the resolution-independent stray-energy term, which does not scale that
/// way, while still failing a first-order scheme — which is the point of
/// asserting an order rather than a bare shrinkage.
const MIN_CONVERGENCE_ORDER: f64 = 1.5;
/// Largest convergence order that same refinement may show.
///
/// A second-order scheme cannot converge faster than second order, so an order
/// well above 2 is not a better result: it means the fine run's error is not the
/// truncation at all — two terms of the budget cancelling by accident, or the
/// measurement's own floor — and the ratio the test reads is then a coincidence
/// rather than a rate. Three is 2 plus the same half-order of slack
/// [`MIN_CONVERGENCE_ORDER`] leaves on the other side.
const MAX_CONVERGENCE_ORDER: f64 = 3.0;

/// Reduced gravity `g'` of the fast ocean, in m/s².
///
/// Four times the Pacific's, which doubles `c = √(g'H)` and so multiplies
/// `Le = √(c/β)` by `√2`. Four rather than some other factor because a `√2` in
/// `Le` is six times the largest budget in this file and still leaves the walls
/// at 4.1 `Le`.
const FAST_OCEAN_REDUCED_GRAVITY_M_PER_S2: f64 = 4.0 * PACIFIC_REDUCED_GRAVITY_M_PER_S2;
/// Beta-plane gradient of the strongly rotating ocean, in m⁻¹s⁻¹.
///
/// Twice the Pacific's, which divides `Le` by `√2` while leaving `c` — and so
/// the packet, the timestep and the flight — exactly as they were. Between them
/// the two altered oceans move `Le` by one factor of `√2` each way, and by
/// changing a different one of its two parameters each time.
const STRONGLY_ROTATING_BETA_PER_M_PER_S: f64 = 2.0 * BETA_PER_M_PER_S;

/// Which of the two equatorial waves a run carries.
///
/// The waves differ in the basin they need, the packet that suits them and the
/// invariant they are read in, and every one of those differences lives here
/// rather than being written down once per run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Wave {
    /// The eastward-travelling Kelvin wave of `CONTEXT.md`.
    Kelvin,
    /// The gravest meridional Rossby mode, travelling west at `c/3`.
    GravestRossby,
}

/// The configuration one wave's runs are built from.
///
/// Basin, packet, flight and speed differ between the two waves and nothing
/// else does, so they are one table per wave rather than one `match` per
/// number: the seven values that have to be read together to see whether a run
/// fits in its basin then sit together.
struct WaveSetup {
    /// Zonal extent of the basin, in metres.
    basin_lx_m: f64,
    /// Cell width of the coarse run, in metres.
    coarse_cell_width_m: f64,
    /// Cell height of the coarse run, in metres.
    coarse_cell_height_m: f64,
    /// Zonal e-folding half-width `σ` of the packet's Gaussian envelope, in
    /// metres.
    packet_width_m: f64,
    /// Zonal position of the packet's centre at `t = 0`, in metres east of the
    /// western wall.
    packet_centre_x_m: f64,
    /// The wave's zonal speed, as a signed multiple of `c` — eastward positive.
    speed_in_c: f64,
    /// How far the packet travels before its profile is read, in its own
    /// widths.
    flight_in_widths: f64,
    /// The wave's name, for the message an assertion fails with.
    name: &'static str,
}

/// The Kelvin pulse's runs.
const KELVIN_SETUP: WaveSetup = WaveSetup {
    // 20 000 km, the order of the equatorial Pacific's width (`CONTEXT.md`,
    // *Basin*). The pulse starts four of its own widths from the western wall
    // and ends the flight four short of the eastern one, so no boundary takes
    // part in what is measured: this is a validation of the *interior* wave,
    // and the boundaries' own physics is the subject of T-04.3 and T-04.4.
    basin_lx_m: 2.0e7,
    // Seven cells to a pulse width, which is what the `(Δx/σ)²/4` entry of the
    // budget costs.
    coarse_cell_width_m: 2.0e5,
    // Seven cells to a deformation radius, which is what the `(Δy/Le)²` entry
    // costs — the dominant term of this wave's budget.
    coarse_cell_height_m: 5.0e4,
    // 1500 km, 4.3 `Le`. The Kelvin branch is an exact solution at every
    // wavenumber, so nothing constrains this from below but the zonal
    // truncation term, which grows as `σ⁻²`; what constrains it from above is
    // the basin's room for the clearances and the flight.
    packet_width_m: 1.5e6,
    // Four widths from the western wall, which the eastward pulse travels away
    // from, so that wall sees `e^{−8} = 3×10⁻⁴` of the pulse and the run starts
    // with an undisturbed boundary.
    packet_centre_x_m: 6.0e6,
    // Eastward at `c` (`CONTEXT.md`, *Kelvin wave*), the speed T-07.1 measured.
    speed_in_c: 1.0,
    // Five widths, which leaves the pulse four short of the eastern wall. Long
    // enough that the meridional error the budget bounds has had the whole
    // flight to accumulate — reading the profile at `t = 0` would measure the
    // fit and not the solver.
    flight_in_widths: 5.0,
    name: "Kelvin pulse",
};

/// The gravest Rossby packet's runs.
const GRAVEST_ROSSBY_SETUP: WaveSetup = WaveSetup {
    // 32 000 km — half again the Kelvin basin, for a packet nearly twice as
    // wide and the same four widths of clearance at each wall.
    basin_lx_m: 3.2e7,
    // Eleven cells to a packet width.
    coarse_cell_width_m: 2.5e5,
    // Fourteen cells to a deformation radius, twice the Kelvin run's. The
    // dominant term of every budget here is `(2m+1)·(Δy/Le)²` and this mode's
    // `2m+1` is five where the Kelvin wave's is one, so resolving the waveguide
    // twice as finely is what holds the two waves to budgets of the same order
    // rather than admitting this fit five times the error for the same physics.
    coarse_cell_height_m: 2.5e4,
    // 2800 km, 8.1 `Le`. Wider than the Kelvin pulse because this mode's
    // initial condition is the long-wave one, whose stray-energy term also goes
    // as `σ⁻²` but from a much larger coefficient.
    packet_width_m: 2.8e6,
    // Four widths from the eastern wall, which the westward packet travels away
    // from.
    packet_centre_x_m: 2.08e7,
    // Westward at `c/3` (`CONTEXT.md`, *Rossby wave*), the speed T-07.2
    // measured: the long-wave root of `ω̂³ − (k̂² + 3)·ω̂ − k̂ = 0` is
    // `ω̂ = −k̂/(2n+1)`, which is `−k̂/3` for the gravest meridional mode.
    speed_in_c: -1.0 / 3.0,
    // Two widths, which leaves the packet 5.4 short of the western wall.
    flight_in_widths: 2.0,
    name: "gravest Rossby packet",
};

impl Wave {
    /// The table this wave's runs are built from.
    const fn setup(self) -> &'static WaveSetup {
        match self {
            Self::Kelvin => &KELVIN_SETUP,
            Self::GravestRossby => &GRAVEST_ROSSBY_SETUP,
        }
    }

    /// The invariant this wave's `ψ₀` content sits in, and the structure it sits
    /// on there — the first row of the module header's table for the Kelvin
    /// wave, the second for the Rossby mode.
    const fn gravest_signature(self) -> (Invariant, MeridionalStructure) {
        match self {
            Self::Kelvin => (Invariant::Eastward, MeridionalStructure::Gravest),
            Self::GravestRossby => (Invariant::Westward, MeridionalStructure::Gravest),
        }
    }

    /// Every (invariant, structure) pair this wave's profile is fitted on — one
    /// row of the module header's table each.
    ///
    /// The Kelvin wave has exactly one: it is `ψ₀` of the eastward invariant
    /// and identically nothing anywhere else. The gravest Rossby mode has two,
    /// and that difference is what
    /// [`the_fitted_decay_scale_converges_at_the_schemes_second_order`] turns
    /// on.
    fn signatures(self) -> Vec<(Invariant, MeridionalStructure)> {
        match self {
            Self::Kelvin => vec![self.gravest_signature()],
            Self::GravestRossby => vec![
                self.gravest_signature(),
                (Invariant::Eastward, MeridionalStructure::Second),
            ],
        }
    }

    /// The finest meridional structure this wave carries: `ψ₀` for the Kelvin
    /// wave, `ψ₂` for the gravest Rossby mode.
    ///
    /// What the second-order meridional truncation of a budget is scaled by, via
    /// [`MeridionalStructure::truncation_richness`].
    const fn finest_structure(self) -> MeridionalStructure {
        match self {
            Self::Kelvin => MeridionalStructure::Gravest,
            Self::GravestRossby => MeridionalStructure::Second,
        }
    }
}

/// One decay-scale experiment: a wave, in an ocean, in a basin at some
/// refinement.
///
/// Wave, ocean and refinement are the three things the runs of this file differ
/// in — the first for "both wave types", the second for the physics of
/// `Le = √(c/β)`, the third for the convergence — and everything derived from
/// them lives here rather than being written down once per run.
#[derive(Debug, Clone, Copy)]
struct Experiment {
    /// Which wave the packet is.
    wave: Wave,
    /// Shape, spacing and position of the basin.
    basin: Basin,
    /// The ocean the equations are written in terms of.
    params: PhysicalParams,
    /// `c = √(g'·H)` of that ocean, in m/s — from the analytic definition,
    /// never from the engine.
    wave_speed_m_per_s: f64,
    /// `Le = √(c/β)` of that ocean, in metres — likewise, and the value every
    /// assertion below is made against.
    deformation_radius_m: f64,
}

impl Experiment {
    /// The experiment carrying `wave` in `params`' ocean, at `1/refinement` of
    /// that wave's coarse cell size.
    fn new(wave: Wave, params: PhysicalParams, refinement: usize) -> Self {
        let cell_width_m = wave.setup().coarse_cell_width_m / refinement as f64;
        let cell_height_m = wave.setup().coarse_cell_height_m / refinement as f64;
        let grid = Grid::new(
            (wave.setup().basin_lx_m / cell_width_m).round() as usize,
            (BASIN_LY_M / cell_height_m).round() as usize,
        )
        .expect("the basin has cells on both axes");
        let spacing = Spacing::new(cell_width_m, cell_height_m)
            .expect("the cell sizes are finite and positive");
        let wave_speed_m_per_s = wave_speed_of_m_per_s(params.reduced_gravity_m_per_s2());
        Self {
            wave,
            basin: Basin::centered_on_equator(grid, spacing),
            params,
            wave_speed_m_per_s,
            deformation_radius_m: deformation_radius_of_m(
                wave_speed_m_per_s,
                params.beta_per_m_per_s(),
            ),
        }
    }

    /// The Gaussian packet this run starts from.
    fn packet(self) -> Packet {
        Packet {
            amplitude_m: PACKET_AMPLITUDE_M,
            centre_x_m: self.wave.setup().packet_centre_x_m,
            width_m: self.wave.setup().packet_width_m,
        }
    }

    /// How long the run integrates for, in seconds: the flight of
    /// [`Wave::flight_in_widths`] packet widths at this wave's own zonal speed.
    fn flight_time_s(self) -> f64 {
        self.wave.setup().flight_in_widths * self.wave.setup().packet_width_m
            / (self.wave.setup().speed_in_c * self.wave_speed_m_per_s).abs()
    }

    /// `(2m+1)·(Δy/Le)²`: the second-order meridional truncation of the
    /// waveguide, as a fraction.
    ///
    /// The factor is [`MeridionalStructure::truncation_richness`] of the finest
    /// structure this wave carries, which is the scale a second-order scheme's
    /// error on it is set by.
    fn meridional_truncation(self) -> f64 {
        let cell_in_radii = self.basin.spacing().dy_m() / self.deformation_radius_m;
        self.wave.finest_structure().truncation_richness() * cell_in_radii * cell_in_radii
    }

    /// `(Δx/σ)²/4`: the zonal truncation, as a fraction.
    fn zonal_truncation(self) -> f64 {
        let cell_in_widths = self.basin.spacing().dx_m() / self.wave.setup().packet_width_m;
        0.25 * cell_in_widths * cell_in_widths
    }

    /// `⟨k̂²⟩ = 1/(2σ̂²)`: the amplitude the initial condition sheds into other
    /// branches, as a fraction.
    ///
    /// The Rossby initial condition is the long-wave mode, which is the exact
    /// one only as `k̂ → 0` and correct to `O(k̂²)` in amplitude; the mean square
    /// wavenumber of the packet's energy spectrum is what that `O(k̂²)` is worth
    /// here. The Kelvin branch is an exact solution at every wavenumber, so its
    /// share of this term is zero rather than small.
    fn stray_amplitude(self) -> f64 {
        match self.wave {
            Wave::Kelvin => 0.0,
            Wave::GravestRossby => {
                let width_in_radii = self.wave.setup().packet_width_m / self.deformation_radius_m;
                0.5 / (width_in_radii * width_in_radii)
            }
        }
    }

    /// The tolerance on this run's fitted decay scale, as a fraction of `Le`:
    /// the three terms of the module header, times [`TRUNCATION_SAFETY`].
    ///
    /// A function of the run rather than a constant, because the point of the
    /// acceptance criterion is that a finer grid is held to a tighter bound.
    fn scale_tolerance(self) -> f64 {
        TRUNCATION_SAFETY
            * (self.meridional_truncation() + self.zonal_truncation() + self.stray_amplitude())
    }

    /// The bracket the trapping-scale fit searches, in metres.
    fn scale_bracket_m(self) -> (f64, f64) {
        let (smallest, largest) = SCALE_BRACKET_IN_RADII;
        (
            smallest * self.deformation_radius_m,
            largest * self.deformation_radius_m,
        )
    }

    /// Run the experiment: the packet in a closed, unforced, undamped basin,
    /// integrated through its whole flight and read once at the end.
    fn run(self) -> Flight {
        let wave_speed = WaveSpeed::new(self.wave_speed_m_per_s).expect("a positive wave speed");
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
        let steps = (self.flight_time_s() / dt_s).round() as usize;

        let mut state = self.initial_state();
        for step in 0..steps {
            solver.step(&mut state, step as f64 * dt_s, |_t| &calm);
        }

        Flight {
            row_y_m: support::row_positions_m(self.basin),
            eastward: self.meridional_profile(&state, Invariant::Eastward),
            westward: self.meridional_profile(&state, Invariant::Westward),
        }
    }

    /// The initial condition: this wave's analytic packet, from the one
    /// definition of it the suite shares.
    fn initial_state(self) -> OceanState {
        let (basin, params, le_m, c) = (
            self.basin,
            self.params,
            self.deformation_radius_m,
            self.wave_speed_m_per_s,
        );
        match self.wave {
            Wave::Kelvin => support::kelvin_pulse_state(basin, params, le_m, c, self.packet()),
            Wave::GravestRossby => {
                support::gravest_rossby_packet_state(basin, params, le_m, c, self.packet())
            }
        }
    }

    /// The zonal sum of one invariant, row by row — the profile the fit reads.
    fn meridional_profile(self, state: &OceanState, invariant: Invariant) -> Vec<f64> {
        support::invariant_meridional_profile(
            self.basin,
            self.params,
            state,
            self.wave_speed_m_per_s,
            invariant,
        )
    }

    /// Fit the flight's `structure` profile: the measurement this ticket is
    /// about, and the goodness of it.
    ///
    /// One entry point rather than one per number the tests read, so that the
    /// golden-section search runs once per profile and the scale an assertion
    /// reports is the same scale the shape error was measured at.
    fn fit(
        self,
        flight: &Flight,
        invariant: Invariant,
        structure: MeridionalStructure,
    ) -> ScaleFit {
        let profile = flight.profile(invariant);
        let scale_m = support::fitted_trapping_scale_m(
            &flight.row_y_m,
            profile,
            structure,
            self.scale_bracket_m(),
        );
        let explained = support::shape_correlation(&flight.row_y_m, profile, structure, scale_m);
        ScaleFit {
            scale_m,
            scale_error: (scale_m - self.deformation_radius_m).abs() / self.deformation_radius_m,
            shape_error: (1.0 - explained).max(0.0).sqrt(),
        }
    }
}

/// What fitting one profile yields.
#[derive(Debug, Clone, Copy)]
struct ScaleFit {
    /// The trapping scale, in metres, the profile is best fitted by.
    scale_m: f64,
    /// How far that scale sits from `Le = √(c/β)`, as a fraction of `Le`.
    scale_error: f64,
    /// How far the profile departs from `ψₘ` at that scale, as a relative
    /// amplitude.
    ///
    /// `√(1 − ρ)`: the correlation `ρ` is the fraction of the profile's energy
    /// the fitted shape explains, so `1 − ρ` is the energy in the departure and
    /// its square root is that departure's amplitude. This is the `ε` the
    /// module header's budget is a budget on, measured directly rather than
    /// through the fitted scale's response to it.
    shape_error: f64,
}

/// What one run leaves for the tests to read: the meridional profile of each
/// invariant at the end of the flight, and the rows they are sampled on.
#[derive(Debug, Clone)]
struct Flight {
    /// Meridional positions of the cell-centre rows, in metres.
    row_y_m: Vec<f64>,
    /// `Σᵢ (u/c + h/H)`, row by row.
    eastward: Vec<f64>,
    /// `Σᵢ (u/c − h/H)`, row by row.
    westward: Vec<f64>,
}

impl Flight {
    /// The profile of one invariant.
    fn profile(&self, invariant: Invariant) -> &[f64] {
        match invariant {
            Invariant::Eastward => &self.eastward,
            Invariant::Westward => &self.westward,
        }
    }
}

/// Assert that `experiment`'s flight is fitted by `structure` on the ocean's own
/// deformation radius, and that the fit explains the profile.
///
/// The two halves of one point check, so that every use of it makes both: a
/// best-fit scale exists for any profile whatever, and a scale reported for a
/// profile that is not that shape would be meaningless rather than wrong.
fn assert_trapped_on_the_deformation_radius(
    experiment: Experiment,
    flight: &Flight,
    invariant: Invariant,
    structure: MeridionalStructure,
) {
    let fit = experiment.fit(flight, invariant, structure);
    let fitted_m = fit.scale_m;
    let expected_m = experiment.deformation_radius_m;
    let error = fit.scale_error;
    let tolerance = experiment.scale_tolerance();
    let cells = experiment.basin.grid();

    assert!(
        error <= tolerance,
        "on the {}×{} grid the {}'s meridional profile decayed on {fitted_m:.0} m, {:.2}% from \
         the deformation radius √(c/β) = {expected_m:.0} m; that grid's truncation budget allows \
         {:.2}%",
        cells.nx(),
        cells.ny(),
        experiment.wave.setup().name,
        100.0 * error,
        100.0 * tolerance
    );

    // A contamination of relative amplitude `ε` moves the fitted scale by
    // `O(ε)` and leaves `ε²` of the profile's energy unexplained (module
    // header), so the same budget, already assembled as an amplitude, bounds
    // the departure directly — not a second tolerance.
    assert!(
        fit.shape_error <= tolerance,
        "at its best-fitting scale the {} profile still departed from ψ{} by {:.3}% of its own \
         amplitude, past the {:.3}% the same truncation budget allows: the profile is not that \
         shape, so the scale fitted to it does not mean what it says",
        experiment.wave.setup().name,
        structure.hermite_order(),
        100.0 * fit.shape_error,
        100.0 * tolerance
    );
}

/// The equatorial-Pacific ocean of `CONTEXT.md`, `Le = 345 km`.
fn pacific() -> PhysicalParams {
    support::pacific_params()
}

/// The same ocean with `g'` quadrupled: `c` doubled, `Le = 488 km`.
fn fast_ocean() -> PhysicalParams {
    params_with_reduced_gravity_and_beta(FAST_OCEAN_REDUCED_GRAVITY_M_PER_S2, BETA_PER_M_PER_S)
}

/// The same ocean with `β` doubled: `c` unchanged, `Le = 244 km`.
fn strongly_rotating_ocean() -> PhysicalParams {
    params_with_reduced_gravity_and_beta(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        STRONGLY_ROTATING_BETA_PER_M_PER_S,
    )
}

/// One Pacific experiment, at the resolution named.
fn pacific_experiment(wave: Wave, refinement: usize) -> Experiment {
    Experiment::new(wave, pacific(), refinement)
}

/// The four Pacific flights — two waves at two resolutions — each integrated
/// once and shared.
///
/// Several tests read several different things out of one flight, which is both
/// cheaper than a run each and stronger: every assertion about a wave is then
/// made about the same wave.
fn pacific_flight(wave: Wave, refinement: usize) -> &'static Flight {
    static FLIGHTS: [[OnceLock<Flight>; 2]; 2] = [
        [OnceLock::new(), OnceLock::new()],
        [OnceLock::new(), OnceLock::new()],
    ];
    let wave_index = match wave {
        Wave::Kelvin => 0,
        Wave::GravestRossby => 1,
    };
    let resolution_index = usize::from(refinement == FINE_REFINEMENT);
    FLIGHTS[wave_index][resolution_index].get_or_init(|| pacific_experiment(wave, refinement).run())
}

/// The two resolutions of one wave, coarse first — the pair every point check
/// is asserted over and every rate is read from.
fn resolutions(wave: Wave) -> [(Experiment, &'static Flight); 2] {
    [COARSE_REFINEMENT, FINE_REFINEMENT].map(|refinement| {
        (
            pacific_experiment(wave, refinement),
            pacific_flight(wave, refinement),
        )
    })
}

#[test]
fn the_kelvin_pulse_decays_meridionally_on_the_deformation_radius() {
    // The acceptance criterion of T-07.3 for the first of the two wave types:
    // the decay scale fitted to the wave the run actually carries is within the
    // documented tolerance of `Le = √(c/β)` (`CONTEXT.md`, *Equatorial
    // deformation radius*) — at more than one grid resolution, each held to its
    // own budget, so the finer run has to be better rather than merely no worse.
    for (experiment, flight) in resolutions(Wave::Kelvin) {
        let (invariant, structure) = Wave::Kelvin.gravest_signature();
        assert_trapped_on_the_deformation_radius(experiment, flight, invariant, structure);
    }
}

#[test]
fn the_gravest_rossby_packet_decays_meridionally_on_the_same_radius() {
    // The same criterion for the second wave type — "both wave types" is the
    // ticket's own wording. The Rossby mode's `ψ₀` content lives in the
    // *westward* invariant, which the Kelvin branch has none of at all, so this
    // fit reads the Rossby packet and nothing else.
    for (experiment, flight) in resolutions(Wave::GravestRossby) {
        let (invariant, structure) = Wave::GravestRossby.gravest_signature();
        assert_trapped_on_the_deformation_radius(experiment, flight, invariant, structure);
    }
}

#[test]
fn the_rossby_packets_off_equatorial_lobes_sit_on_the_same_radius() {
    // `Le` is the scale of the *waveguide*, not of one convenient Gaussian. The
    // gravest Rossby mode's eastward invariant is `ψ₂` — two off-equatorial
    // lobes with a node at `ŷ = ±1/√2` and a sign change across it — and it is
    // stretched by the same `Le` as everything else. Fitting that shape is a
    // sharper statement than fitting `ψ₀`: a run that had trapped the wave on
    // the wrong scale would put the node in the wrong place, and no amplitude
    // can absorb that.
    let (experiment, flight) = resolutions(Wave::GravestRossby)[0];
    assert_trapped_on_the_deformation_radius(
        experiment,
        flight,
        Invariant::Eastward,
        MeridionalStructure::Second,
    );
}

#[test]
fn the_fitted_decay_scale_converges_at_the_schemes_second_order() {
    // CODING_STANDARDS.md § *Convergence over point checks*: the point checks
    // above are bounds, and a bound is passed by an error of any size below it.
    // What ties the error to the discretisation is its rate. Both truncation
    // terms of every budget here are second order in the cell size (module
    // header) and both cell dimensions are halved, so the error must fall by
    // about four — for both waves, and for both of the structures the Rossby
    // mode is fitted on.
    for wave in [Wave::Kelvin, Wave::GravestRossby] {
        for (invariant, structure) in wave.signatures() {
            let [(coarse, coarse_flight), (fine, fine_flight)] = resolutions(wave);
            let coarse_error = coarse.fit(coarse_flight, invariant, structure).scale_error;
            let fine_error = fine.fit(fine_flight, invariant, structure).scale_error;
            assert_second_order(
                wave,
                structure,
                "fitted decay scale",
                coarse_error,
                fine_error,
            );
        }
    }
}

#[test]
fn the_meridional_shape_error_converges_at_the_schemes_second_order() {
    // The same rate, read on the quantity every budget in this file is a budget
    // on: `ε`, the amplitude by which a run's meridional profile departs from
    // the `ψₘ` the theory says it is. The dominant entry bounds that departure
    // at `(2m+1)·(Δy/Le)²`, so if the budget describes what the scheme does and
    // not merely what it is allowed, halving both cell dimensions must quarter
    // it.
    //
    // It is the sharper of the two rate tests: a fitted scale responds only to
    // the part of a departure that is not orthogonal to it, whereas this sees
    // the whole departure.
    for wave in [Wave::Kelvin, Wave::GravestRossby] {
        for (invariant, structure) in wave.signatures() {
            let [(coarse, coarse_flight), (fine, fine_flight)] = resolutions(wave);
            let coarse_error = coarse.fit(coarse_flight, invariant, structure).shape_error;
            let fine_error = fine.fit(fine_flight, invariant, structure).shape_error;
            assert_second_order(
                wave,
                structure,
                "departure from ψₘ",
                coarse_error,
                fine_error,
            );
        }
    }
}

/// Assert that halving both cell dimensions took `quantity`'s error from
/// `coarse_error` to `fine_error` at the scheme's second order.
///
/// Bounded on both sides: too small an order is a scheme that is not second
/// order, and too large a one is a fine error that is no longer the truncation
/// — two terms of the budget cancelling by accident, or the measurement's own
/// floor — and the ratio is then a coincidence rather than a rate. Either way
/// the point checks' budget has stopped describing what the run does.
fn assert_second_order(
    wave: Wave,
    structure: MeridionalStructure,
    quantity: &str,
    coarse_error: f64,
    fine_error: f64,
) {
    assert!(
        fine_error > 0.0,
        "the fine {} run reproduced ψ{} exactly, so there is no {quantity} error left to measure \
         a rate on: what this run reads is the measurement's floor and not the scheme's order",
        wave.setup().name,
        structure.hermite_order()
    );
    let order = (coarse_error / fine_error).log2();
    assert!(
        (MIN_CONVERGENCE_ORDER..=MAX_CONVERGENCE_ORDER).contains(&order),
        "refining the cells by {FINE_REFINEMENT} took the {}'s {quantity}, read on ψ{}, from \
         {:.4}% to {:.4}%, a convergence order of {order:.2}: outside the \
         [{MIN_CONVERGENCE_ORDER}, {MAX_CONVERGENCE_ORDER}] a second-order scheme owes",
        wave.setup().name,
        structure.hermite_order(),
        100.0 * coarse_error,
        100.0 * fine_error
    );
}

#[test]
fn the_decay_scale_follows_le_across_oceans() {
    // The physics, as opposed to the numerics: `Le = √(c/β)` is a *prediction*,
    // and one ocean cannot tell it apart from a length that happens to fit. The
    // same pulse, basin, packet and measurement in three oceans whose radii sit
    // a factor of `√2` apart — 244 km, 345 km, 488 km — must therefore give
    // three different answers, each its own ocean's.
    //
    // Stated as "nearer its own than either neighbour's" rather than as a bare
    // budget, because that is the discriminating claim: the neighbours are 29%
    // and 41% away and the largest budget here is 9.3%, so a fit that returned
    // a fixed length would fail in two of the three oceans.
    let oceans = [
        ("the equatorial Pacific", pacific()),
        ("the fast ocean", fast_ocean()),
        ("the strongly rotating ocean", strongly_rotating_ocean()),
    ];
    let radii_m: Vec<f64> = oceans
        .iter()
        .map(|(_, params)| Experiment::new(Wave::Kelvin, *params, COARSE_REFINEMENT))
        .map(|experiment| experiment.deformation_radius_m)
        .collect();

    for (index, (name, params)) in oceans.iter().enumerate() {
        let experiment = Experiment::new(Wave::Kelvin, *params, COARSE_REFINEMENT);
        let flight = experiment.run();
        let (invariant, structure) = Wave::Kelvin.gravest_signature();
        assert_trapped_on_the_deformation_radius(experiment, &flight, invariant, structure);

        let fitted_m = experiment.fit(&flight, invariant, structure).scale_m;
        let own_m = experiment.deformation_radius_m;
        for (other, other_m) in radii_m.iter().enumerate() {
            if other == index {
                continue;
            }
            assert!(
                (fitted_m - own_m).abs() < (fitted_m - other_m).abs(),
                "in {name} the pulse decayed on {fitted_m:.0} m, which is nearer the \
                 {other_m:.0} m of {} than its own √(c/β) = {own_m:.0} m: the fitted scale is \
                 not following the ocean it was measured in",
                oceans[other].0
            );
        }
    }
}
