//! Acceptance tests for T-07.2 — the gravest equatorial Rossby mode's westward
//! propagation at `c/3`, and its dispersion.
//!
//! `CONTEXT.md` (*Rossby wave*) states the claim this file validates: an
//! equatorial Rossby disturbance of the first baroclinic mode travels
//! **westward**, and the gravest meridional mode's long-wave speed is `c/3`
//! with `c = √(g'·H)`. `docs/planning/01-scientific-model.md` makes the same
//! statement of the 1.5-layer core, and ADR-0003 fixes the discretisation whose
//! truncation error is what stands between the two.
//!
//! This is T-07.1's ticket for the other branch, and it is deliberately its
//! mirror: the same basin-interior geometry, the same modal projection, the
//! same centroid measurement, the same two-resolution convergence. What is not
//! a mirror is *why* the measured speed differs from the headline number. The
//! Kelvin branch is non-dispersive, so any departure from `c` is numerical. The
//! Rossby branch is not, so a packet of finite zonal width travels measurably
//! slower than `c/3` for a reason that has nothing to do with the grid — and
//! the deliverable's "matching the equatorial Rossby dispersion relation" is
//! precisely the claim that the departure is *that* one.
//!
//! Nothing below is measured out of a run and pasted back in as an expectation.
//! Every asserted number is `c/3`, a quantity computed from the analytic
//! dispersion relation alone ([`mean_group_speed_in_c`]), zero, or a tolerance
//! assembled from named error terms of the configuration (CODING_STANDARDS.md
//! § *Tests*).
//!
//! # The mode
//!
//! Write the linear equations of [ADR-0003] in the equatorial variables:
//! `ŷ = y/Le` with `Le = √(c/β)` the equatorial deformation radius, `x̂ = x/Le`,
//! velocities scaled by `c`, `h` by the mean thermocline depth `H`, and time by
//! `Le/c`. Separating on the parabolic cylinder functions
//! `ψₘ(ŷ) = Hₘ(ŷ)·e^{−ŷ²/2}` with the meridional velocity on `ψₙ` gives the
//! equatorial dispersion relation
//!
//! ```text
//! ω̂³ − (k̂² + 2n + 1)·ω̂ − k̂ = 0,     k̂ = k·Le,  ω̂ = ω·Le/c.
//! ```
//!
//! Its small-`|ω̂|` root is the Rossby branch. `n = 1` is the gravest
//! meridional Rossby mode, so [`gravest_rossby_frequency`] solves
//! `ω̂³ − (k̂² + 3)·ω̂ − k̂ = 0`; at `k̂ → 0` that root is `ω̂ = −k̂/3`, which is
//! the westward `c/3` of `CONTEXT.md`.
//!
//! With `v = V(ŷ)e^{i(k̂x̂−ω̂t̂)}` and `V = ψ₁`, the ladder relations
//! `(∂/∂ŷ + ŷ)ψₘ = 2m·ψₘ₋₁` and `(∂/∂ŷ − ŷ)ψₘ = −ψₘ₊₁` give the two invariants
//! of this file exactly, at every wavenumber and not only in the long-wave
//! limit:
//!
//! ```text
//! eastward = u/c + h/H = i·ψ₂/(ω̂ − k̂)      westward = u/c − h/H = 2i·ψ₀/(ω̂ + k̂)
//! ```
//!
//! So the gravest Rossby mode is `ψ₂` of the eastward invariant and `ψ₀` of the
//! westward one, purely; the Kelvin wave is `ψ₀` of the eastward one and
//! nothing else (Matsuno 1966; Gill, *Atmosphere–Ocean Dynamics*, § 11.6).
//! `P₀[u/c − h/H]` is therefore the channel that holds this wave and no Kelvin
//! wave at all, and it is the profile every measurement below reads. This is
//! the same decomposition `kelvin_wave_propagation.rs` and the two boundary
//! reflection tests use, turned here on a Rossby packet in open water; the
//! machinery is `tests/support/mod.rs`.
//!
//! # The initial condition
//!
//! Reading those relations back at `ω̂ = −k̂/3` gives the long-wave fields, with
//! a Gaussian zonal envelope `E(x)` of e-folding half-width `σ`:
//!
//! ```text
//! h/H = A·E(x)·(2ŷ² + 1)·e^{−ŷ²/2}
//! u/c = A·E(x)·(2ŷ² − 3)·e^{−ŷ²/2}
//! v/c = (8/3)·A·Le·(dE/dx)·ŷ·e^{−ŷ²/2}
//! ```
//!
//! The `h` of that mode is double-lobed off the equator — the familiar
//! off-equatorial Rossby signature — which is what
//! [`the_rossby_packet_keeps_the_gravest_modes_meridional_shape`] checks it
//! still is at the end of the run.
//!
//! The meridional velocity is written down rather than left at rest, which is
//! where this file departs from `western_boundary_reflection.rs`. For the exact
//! mode `V` is `O(k̂)` relative to `u`, so a run started from `v = 0` is the
//! Rossby mode plus an `O(k̂)` admixture of the `n = 1` inertia-gravity pair,
//! carrying `O(k̂²)` — better than a percent — of the energy into branches this
//! test is not measuring. A boundary reflection times a *peak* at a station and
//! is insensitive to that; a centroid weighs the whole basin and is not.
//! Including `v` costs three lines and takes the admixture to `O(k̂²)` in
//! amplitude and `O(k̂⁴)` in energy, which is the [`Experiment::stray_energy`]
//! term of the budget below.
//!
//! # What is measured
//!
//! The packet's position is the energy-weighted zonal centroid of
//! `P₀[u/c − h/H]` over the whole basin, and its speed is the displacement
//! between two sample times divided by the elapsed time. The two times are the
//! *steps'* own times, not the requested ones, so the denominator carries no
//! rounding.
//!
//! A centroid is the right measurement for a *dispersive* wave, and the reason
//! is the whole point of this ticket. The centroid of a linear packet moves at
//! the energy-weighted mean of the group velocity over its own spectrum,
//!
//! ```text
//! ⟨c_g⟩ = ∫ c_g(k)·|Â(k)|² dk / ∫ |Â(k)|² dk,
//! ```
//!
//! exactly and with no expansion in `k̂`. The Gaussian envelope has
//! `|Â(k)|² = e^{−k²σ²}`, and `c_g(k) = c·dω̂/dk̂` comes from implicit
//! differentiation of the cubic above, so [`mean_group_speed_in_c`] evaluates
//! that quotient by quadrature and hands back the speed this run *should*
//! show — from the dispersion relation alone, without the engine. Expanded, it
//! is the familiar
//!
//! ```text
//! ω̂ = −k̂/3 + (8/81)·k̂³,   c_g = −c/3 + (8/27)·c·k̂²,
//! ⟨c_g⟩ = −(c/3)·(1 − (4/9)·(Le/σ)²),
//! ```
//!
//! a `0.68%` slowdown at [`PACKET_WIDTH_M`] and `2.7%` at the narrower width —
//! but the quadrature is used in preference to the expansion, so that the `k̂⁵`
//! term is not a budget entry.
//!
//! One departure from T-07.1 is worth naming rather than leaving implicit.
//! Its acceptance criterion says *phase* speed, and what a centroid reads is a
//! *group* speed. On the Kelvin branch that distinction is empty — `ω = c·k` is
//! linear, so `ω/k = ∂ω/∂k = c` at every wavenumber — which is exactly why
//! T-07.1 could use the words interchangeably. On this branch they are
//! genuinely different numbers, and the group speed is the one that has to be
//! measured: it is what a packet's envelope travels at, so it is what a
//! position measured twice can read at all. `c/3` is the `k̂ → 0` limit of both,
//! which is what keeps the headline claim the same claim.
//!
//! Nothing else in the basin disturbs that centroid. The Kelvin branch has no
//! `ψ₀` in the westward invariant at all, so whatever Kelvin energy the run
//! carries is invisible to it by construction; the packet stays four widths
//! clear of both zonal walls; and what the initial condition sheds into other
//! branches is the [`Experiment::stray_energy`] term.
//!
//! # Where the tolerances come from
//!
//! Every entry is a property of the configuration below and none was obtained
//! by running the model. The two truncation terms are second order in the cell
//! size, which is why the whole error converges at second order.
//!
//! | term | coarse size | why |
//! |---|---|---|
//! | meridional truncation | `5·(Δy/Le)² = 10.5%` | the C-grid operators of T-01.1 are second order; the finest structure this wave carries is `ψ₂`, which oscillates on `Le/√5` rather than `Le`, so its bound carries that factor — [`MeridionalStructure::truncation_richness`] |
//! | zonal group-speed truncation | `(Δx/σ)²/4 = 0.20%` | the centred difference sees `k_eff = sin(kΔx)/Δx`, so a group speed is read low by `(kΔx)²/2`; the energy spectrum of a Gaussian of width `σ` has `⟨k²⟩ = 1/(2σ²)` |
//! | stray energy | `0.013%` | the initial condition is the long-wave mode, exact to `O(k̂²)` in amplitude and so `⟨k̂²⟩² = 5.8×10⁻⁵` in energy; weighed at a lever arm of half the basin over the distance the packet travels — [`Experiment::stray_energy`] |
//! | RK4 phase error | `2×10⁻¹¹` | fourth order: `(ω·Δt)⁴/120` at the packet's dominant frequency `ω = c/(3σ)` and this run's CFL-stable timestep |
//! | wall clipping | `e^{−16}` | the packet starts the run 4 `σ` from the eastern wall and ends it 4.9 `σ` from the western one, so neither ever sees more than `e^{−8}` of its amplitude |
//!
//! The last two are six and many orders below the rest and are not carried into
//! the number. What is carried is the sum of the first three, times
//! [`TRUNCATION_SAFETY`]: none of the leading coefficients is evaluated, and a
//! factor of two is the standing allowance for an unevaluated `O(1)`
//! coefficient. It multiplies both resolutions equally, so it changes what the
//! point checks allow and not what the convergence test measures.
//!
//! At the coarse resolution that is 21%, at the fine one 5.4% — the budget
//! itself shrinks by nearly four when the cells are halved, so passing it at
//! both resolutions is already the "not a fixed offset" the acceptance
//! criterion asks for. Both are bounds and not estimates, so the point checks
//! are generous by design; what pins the error's *size* to the discretisation
//! is [`the_rossby_speed_error_shrinks_at_the_schemes_second_order`], and what
//! pins the *physics* is
//! [`narrowing_the_packet_slows_it_by_what_the_dispersion_relation_says`],
//! which varies the one parameter no grid refinement touches and whose budget
//! is half the effect it measures, the meridional truncation having cancelled
//! out of it.
//!
//! That last test is why the reference packet is as wide as it is. The
//! stray-energy term is the only one that does not shrink with the cells, so it
//! is the floor a convergence rate is read against; it goes as `σ⁻⁴` while the
//! truncation terms do not depend on `σ` at all, and at [`PACKET_WIDTH_M`] it
//! sits two orders of magnitude below the fine run's own truncation budget,
//! which is what leaves a rate there to be read.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

mod support;

use std::sync::OnceLock;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing, WaveSpeed,
    WindStressField, H_STAGGERING, U_STAGGERING, V_STAGGERING,
};

use support::{
    equatorial_deformation_radius_m, gaussian_envelope, kelvin_wave_speed_m_per_s, pacific_params,
    MeridionalStructure,
};

/// Zonal extent of the test basin, in metres — 32 000 km.
///
/// Two thirds again the equatorial Pacific's width (`CONTEXT.md`, *Basin*), and
/// deliberately so: the packet starts four of its own widths from the eastern
/// wall and ends the run 4.9 from the western one, so neither boundary takes
/// part in what is measured. This is a validation of the *interior* wave; the
/// boundaries' own physics is the subject of T-04.3 and T-04.4.
const BASIN_LX_M: f64 = 3.2e7;
/// Meridional extent of the test basin, in metres: ±2000 km, about ±5.8 `Le`.
///
/// Far enough out that `ψ₂` is `3×10⁻⁶` of its peak at the northern and
/// southern walls, so the equatorial waveguide does not feel them, and near
/// enough that the gravity-wave CFL bound stays the binding one rather than
/// ADR-0007's rotation bound.
const BASIN_LY_M: f64 = 4.0e6;
/// Cell width of the coarse run, in metres — eleven to a packet width, which is
/// what the `(Δx/σ)²/4` entry of the speed budget costs.
const COARSE_CELL_WIDTH_M: f64 = 2.5e5;
/// Cell height of the coarse run, in metres — about seven to an equatorial
/// deformation radius, which is what the `5·(Δy/Le)²` entry costs. The dominant
/// term of every budget in this file.
const COARSE_CELL_HEIGHT_M: f64 = 5.0e4;
/// Factor by which the coarse run refines the two cell sizes above: one, it
/// being the run they state. Named so that neither resolution is a bare literal
/// where a run is asked for.
const COARSE_REFINEMENT: usize = 1;
/// Factor by which the fine run refines both cell dimensions.
///
/// Both axes together, so that both truncation terms of the speed budget shrink
/// by the same four and the measured convergence rate is the scheme's rather
/// than a mixture of two rates.
const FINE_REFINEMENT: usize = 2;

/// Zonal e-folding half-width `σ` of the reference packet's Gaussian envelope,
/// in metres — 8.1 `Le`.
///
/// Wide enough that the packet is a recognisable long wave, that its `c/3`
/// slowdown is a 0.68% correction rather than a different regime, and — the
/// binding constraint — that the `σ⁻⁴` stray-energy floor stays well under the
/// fine run's truncation, which is what leaves a convergence rate to measure.
/// Narrow enough that four widths of clearance at the eastern wall, 4.9 at the
/// western one and two and a half widths of flight between them fit in
/// [`BASIN_LX_M`].
const PACKET_WIDTH_M: f64 = 2.8e6;
/// Factor by which the second packet of the dispersion test is narrowed.
///
/// The slowdown goes as `σ⁻²`, so this is a factor of four in it: 0.68% at the
/// reference width against 2.7% here, a gap of 2 percentage points that the
/// budget of [`narrowing_the_packet_slows_it_by_what_the_dispersion_relation_says`]
/// covers twice over. Far enough apart to measure against, close enough that
/// the narrower packet's own slowdown is still a correction to `c/3` rather
/// than a different regime, and it fits the same basin with more clearance than
/// the wider one.
const NARROW_WIDTH_FACTOR: f64 = 2.0;
/// Zonal position of every packet's centre at `t = 0`, in metres east of the
/// western wall.
///
/// Four reference widths short of the eastern wall, so that wall sees `e^{−8}`
/// of the packet's amplitude and the run starts with an undisturbed boundary.
/// Both packets start at the same place, so both measure the same flight.
const PACKET_CENTRE_X_M: f64 = 2.08e7;
/// Peak thermocline depth anomaly of the packet, in metres — the scale a
/// westerly wind burst leaves behind (`CONTEXT.md`, *Westerly wind burst*).
///
/// The core is linear, so every speed and ratio below is independent of this
/// number; it only sets the scale the diagnostics are reported in.
const PACKET_AMPLITUDE_M: f64 = 10.0;

/// When the packet is first sampled, in transits of one reference width at
/// `c/3`.
///
/// Half a width in: far enough that the transient the continuous initial
/// condition sheds on the discrete grid has separated from the packet, near
/// enough to leave a long baseline for the speed.
const SAMPLE_EARLY_IN_TRANSITS: f64 = 0.5;
/// When it is sampled again, in the same transits. The centre has then
/// travelled two and a half widths and is 4.9 short of the western wall.
const SAMPLE_LATE_IN_TRANSITS: f64 = 3.0;

/// Zonal half-width of the window the packet's meridional shape is read in, in
/// packet widths. Two `σ` either side of the centroid holds 95% of it.
const SHAPE_WINDOW_IN_WIDTHS: f64 = 2.0;

/// Factor applied to every truncation-derived bound in this file, for the
/// leading `O(1)` coefficients the truncation terms do not evaluate.
///
/// Two. It multiplies the coarse and the fine budgets alike, so it widens what
/// the point checks admit without touching the rate
/// [`the_rossby_speed_error_shrinks_at_the_schemes_second_order`] measures.
const TRUNCATION_SAFETY: f64 = 2.0;

/// Smallest convergence order the Rossby speed's error must show when both cell
/// dimensions are halved.
///
/// The spatial discretisation is second order (ADR-0003) and both truncation
/// terms of the speed budget are second order in the cell size, so the measured
/// order `log₂(coarse error / fine error)` should be 2. Requiring 1.5 leaves
/// margin for the resolution-independent stray-energy term, which does not
/// scale that way, while still failing a first-order scheme — which is the
/// point of asserting an order rather than a bare shrinkage.
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

/// The `ψ₂`-to-`ψ₀` ratio of the gravest Rossby mode's thermocline anomaly.
///
/// `h/H = A·E(x)·(2ŷ² + 1)·e^{−ŷ²/2}` and `2ŷ² + 1 = 2·ψ₀ + ψ₂/2`, so the
/// ratio is `(1/2)/2 = 1/4`. It is what makes the mode's thermocline signature
/// double-lobed *off* the equator, where the Kelvin wave's is single-lobed on
/// it.
const GRAVEST_ROSSBY_SHAPE_RATIO: f64 = 0.25;

/// Meridional mode number `n` of the gravest Rossby wave, as the `2n + 1` it
/// enters the dispersion relation as.
///
/// One constant rather than two, because the `3` of `ω̂³ − (k̂² + 3)ω̂ − k̂ = 0`
/// and the `3` of `c/3` are the same number: the long-wave root of that cubic
/// is `ω̂ = −k̂/(2n + 1)`, so a mode's long-wave speed is `c/(2n + 1)`. Writing
/// the speed as `c` over this term rather than over a literal is what keeps the
/// two from drifting apart.
const GRAVEST_ROSSBY_MERIDIONAL_TERM: f64 = 3.0;
/// The same `2n + 1` for the *second* meridional Rossby mode, `n = 2`.
///
/// Its long-wave speed `c/5` is the nearest analytic speed on the slow side of
/// `c/3`, and
/// [`the_rossby_packet_travels_west_at_a_third_of_the_kelvin_wave_speed`] uses
/// it as one of the two neighbours the measurement has to be nearer `c/3` than.
const SECOND_ROSSBY_MERIDIONAL_TERM: f64 = 5.0;

/// How far out in `k̂` the spectral quadrature of [`mean_group_speed_in_c`]
/// integrates, in units of `1/σ̂`.
///
/// Eight. The energy spectrum is `e^{−k̂²σ̂²}`, so the tail beyond this holds
/// `e^{−64}` of the packet — twenty orders of magnitude and more below the
/// smallest term in any budget here.
const SPECTRUM_CUTOFF_IN_INVERSE_WIDTHS: f64 = 8.0;
/// Number of Simpson intervals that quadrature uses. Even, as Simpson's rule
/// requires. The integrand is analytic and slowly varying on the scale `1/σ̂`,
/// so its fourth-order error at this many points is below `10⁻¹²`.
const SPECTRUM_QUADRATURE_INTERVALS: usize = 4096;
/// How small `|ω̂³ − (k̂² + 3)ω̂ − k̂|` must be for a root to count as found.
///
/// The cubic's coefficients are `O(1)` and its Rossby root is simple and well
/// separated from the two inertia-gravity roots, so Newton's method reaches
/// machine precision here; this is that, with room for the residual's own
/// rounding.
const DISPERSION_ROOT_TOLERANCE: f64 = 1.0e-14;
/// How many Newton steps that is allowed to take. Newton doubles the correct
/// digits each step from a starting point this close, so this is far more than
/// enough and exists only to make a non-converging root an error rather than a
/// hang.
const DISPERSION_ROOT_MAX_STEPS: usize = 100;

/// The Rossby root `ω̂(k̂)` of the equatorial dispersion relation
/// `ω̂³ − (k̂² + 3)·ω̂ − k̂ = 0` for the gravest meridional mode.
///
/// Solved by Newton's method from the long-wave root `−k̂/(k̂² + 3)`, which is
/// what the cubic reduces to when its `ω̂³` term is dropped — the same branch,
/// and near enough to it that Newton stays on it.
///
/// # Panics
/// If the iteration does not converge, which would mean the starting point has
/// left the Rossby branch rather than that the run is wrong.
fn gravest_rossby_frequency(k_le: f64) -> f64 {
    let meridional = k_le * k_le + GRAVEST_ROSSBY_MERIDIONAL_TERM;
    let mut omega = -k_le / meridional;
    for _ in 0..DISPERSION_ROOT_MAX_STEPS {
        let residual = omega * omega * omega - meridional * omega - k_le;
        if residual.abs() <= DISPERSION_ROOT_TOLERANCE {
            return omega;
        }
        omega -= residual / (3.0 * omega * omega - meridional);
    }
    panic!("Newton's method did not reach the Rossby root of the dispersion relation at k̂ = {k_le}")
}

/// The group speed `dω̂/dk̂` of that branch, in units of `c`.
///
/// By implicit differentiation of the cubic: `3ω̂²ω̂' − 2k̂ω̂ − (k̂² + 3)ω̂' = 1`.
/// At `k̂ = 0` it is `−1/3`, the westward long-wave speed of `CONTEXT.md`.
fn gravest_rossby_group_speed_in_c(k_le: f64) -> f64 {
    let omega = gravest_rossby_frequency(k_le);
    (2.0 * k_le * omega + 1.0)
        / (3.0 * omega * omega - k_le * k_le - GRAVEST_ROSSBY_MERIDIONAL_TERM)
}

/// The energy-weighted mean group speed of a Gaussian packet `width_in_radii`
/// deformation radii wide, in units of `c`. Negative: this branch goes west.
///
/// `∫c_g(k̂)·e^{−k̂²σ̂²} dk̂ / ∫e^{−k̂²σ̂²} dk̂` by Simpson's rule over the half-line
/// — both `c_g` and the weight are even in `k̂`, so the two halves are equal and
/// cancel between numerator and denominator. This is the speed the run's
/// centroid should show, computed from the dispersion relation and the packet's
/// own spectrum with no reference to the engine.
fn mean_group_speed_in_c(width_in_radii: f64) -> f64 {
    let cutoff = SPECTRUM_CUTOFF_IN_INVERSE_WIDTHS / width_in_radii;
    let step = cutoff / SPECTRUM_QUADRATURE_INTERVALS as f64;
    let (mut moment, mut weight) = (0.0, 0.0);
    for interval in 0..=SPECTRUM_QUADRATURE_INTERVALS {
        let k_le = step * interval as f64;
        let simpson = if interval == 0 || interval == SPECTRUM_QUADRATURE_INTERVALS {
            1.0
        } else if interval % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let energy = simpson * (-k_le * k_le * width_in_radii * width_in_radii).exp();
        moment += energy * gravest_rossby_group_speed_in_c(k_le);
        weight += energy;
    }
    moment / weight
}

/// One propagation experiment: a basin at some refinement, carrying a gravest
/// Rossby packet of some zonal width.
///
/// Refinement and width are the two things the runs of this file differ in —
/// the first for the convergence of the numerics, the second for the physics of
/// the dispersion — and everything derived from them lives here rather than
/// being written down once per run.
#[derive(Debug, Clone, Copy)]
struct Experiment {
    /// Shape, spacing and position of the basin.
    basin: Basin,
    /// The ocean the equations are written in terms of.
    params: PhysicalParams,
    /// `Le = √(c/β)`, in metres — cached because every projection needs it.
    deformation_radius_m: f64,
    /// Zonal width `σ` of this run's packet, in metres.
    packet_width_m: f64,
}

impl Experiment {
    /// The experiment at `1/refinement` of the coarse cell size, carrying a
    /// packet `packet_width_m` wide.
    fn new(refinement: usize, packet_width_m: f64) -> Self {
        let cell_width_m = COARSE_CELL_WIDTH_M / refinement as f64;
        let cell_height_m = COARSE_CELL_HEIGHT_M / refinement as f64;
        let grid = Grid::new(
            (BASIN_LX_M / cell_width_m).round() as usize,
            (BASIN_LY_M / cell_height_m).round() as usize,
        )
        .expect("the basin has cells on both axes");
        let spacing = Spacing::new(cell_width_m, cell_height_m)
            .expect("the cell sizes are finite and positive");
        Self {
            basin: Basin::centered_on_equator(grid, spacing),
            params: pacific_params(),
            deformation_radius_m: equatorial_deformation_radius_m(),
            packet_width_m,
        }
    }

    /// The Kelvin wave speed `c = √(g'·H)`, in m/s: the analytic
    /// [`kelvin_wave_speed_m_per_s`], never the engine's own.
    fn wave_speed_m_per_s(self) -> f64 {
        kelvin_wave_speed_m_per_s()
    }

    /// The long-wave speed of the gravest Rossby mode, in m/s — `−c/3`, the
    /// headline number of `CONTEXT.md`. Negative: westward.
    fn long_wave_speed_m_per_s(self) -> f64 {
        -self.wave_speed_m_per_s() / GRAVEST_ROSSBY_MERIDIONAL_TERM
    }

    /// The speed this run's centroid should show, in m/s, from the dispersion
    /// relation and this packet's own spectrum. Negative: westward.
    fn predicted_speed_m_per_s(self) -> f64 {
        self.wave_speed_m_per_s() * mean_group_speed_in_c(self.width_in_radii())
    }

    /// How far short of `c/3` that prediction falls, as a fraction of `c/3`.
    ///
    /// A property of the packet and not of the grid: it is what the *continuous*
    /// equations do to a disturbance of finite zonal extent, and it would not
    /// move if the grid were refined to nothing.
    fn predicted_slowdown(self) -> f64 {
        1.0 - self.predicted_speed_m_per_s() / self.long_wave_speed_m_per_s()
    }

    /// The packet's width in deformation radii, `σ̂ = σ/Le` — the one number the
    /// dispersion relation needs about it.
    fn width_in_radii(self) -> f64 {
        self.packet_width_m / self.deformation_radius_m
    }

    /// `⟨k̂²⟩ = 1/(2σ̂²)`, the mean square wavenumber of the packet's energy
    /// spectrum in units of `1/Le`.
    fn mean_square_wavenumber(self) -> f64 {
        0.5 / (self.width_in_radii() * self.width_in_radii())
    }

    /// The two sample times of every run, in seconds.
    ///
    /// Stated in transits of the *reference* width at `c/3` so that the narrow
    /// run samples the same clock, and therefore very nearly the same flight, as
    /// the wide one: the two speeds are then comparable without a correction.
    fn sample_times_s(self) -> [f64; 2] {
        let transit_s = PACKET_WIDTH_M / self.long_wave_speed_m_per_s().abs();
        [SAMPLE_EARLY_IN_TRANSITS, SAMPLE_LATE_IN_TRANSITS].map(|transits| transits * transit_s)
    }

    /// `5·(Δy/Le)²`: the second-order meridional truncation of the waveguide,
    /// as a fraction.
    ///
    /// The factor is [`MeridionalStructure::truncation_richness`] of `ψ₂` — the
    /// finest meridional structure this wave carries, and the one a
    /// second-order scheme's error on this mode is set by.
    fn meridional_truncation(self) -> f64 {
        let cell_in_radii = self.basin.spacing().dy_m() / self.deformation_radius_m;
        MeridionalStructure::Second.truncation_richness() * cell_in_radii * cell_in_radii
    }

    /// `(Δx/σ)²/4`: the zonal group-speed truncation, as a fraction.
    fn zonal_truncation(self) -> f64 {
        let cell_in_widths = self.basin.spacing().dx_m() / self.packet_width_m;
        0.25 * cell_in_widths * cell_in_widths
    }

    /// The share of a speed the energy shed by the initial condition can move
    /// the centroid by.
    ///
    /// The long-wave fields are the exact mode to `O(k̂²)` in amplitude, so
    /// `⟨k̂²⟩²` of the run's energy is in branches this measurement does not
    /// follow. Bound its lever arm by the whole basin and its displacement by
    /// the whole measurement baseline, both of which are the most it could be.
    fn stray_energy(self) -> f64 {
        let shed_share = self.mean_square_wavenumber() * self.mean_square_wavenumber();
        let [early_s, late_s] = self.sample_times_s();
        let travel_m = self.long_wave_speed_m_per_s().abs() * (late_s - early_s);
        shed_share * (0.5 * BASIN_LX_M / travel_m)
    }

    /// The tolerance on this run's measured speed, as a fraction of it: the
    /// three terms of the module header, times [`TRUNCATION_SAFETY`].
    ///
    /// A function of the run rather than a constant, because the point of the
    /// acceptance criterion is that a finer grid is held to a tighter bound.
    fn speed_tolerance(self) -> f64 {
        TRUNCATION_SAFETY
            * (self.meridional_truncation() + self.zonal_truncation() + self.stray_energy())
    }

    /// The thermocline depth anomaly's `ψₘ` coefficient, column by column, in
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

    /// The two `ψ₀` invariants of `state`, column by column.
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

    /// `Σ profile²` over the whole basin.
    fn energy(self, profile: &[f64]) -> f64 {
        support::energy(profile.iter().copied())
    }

    /// The energy-weighted zonal centroid of `profile`, in metres.
    ///
    /// Over the whole basin, with no window: a window would need a position to
    /// be centred on, and taking that from theory is exactly the circularity a
    /// speed measurement must not have.
    fn energy_centroid_m(self, profile: &[f64]) -> f64 {
        support::energy_centroid_m(
            profile
                .iter()
                .enumerate()
                .map(|(i, amplitude)| (self.column_x_m(i), *amplitude)),
        )
    }

    /// The initial condition: the gravest-mode equatorial Rossby packet of the
    /// module header, Gaussian in `x`, centred at [`PACKET_CENTRE_X_M`] and
    /// travelling west.
    fn initial_state(self) -> OceanState {
        let mut state = OceanState::at_rest(self.basin.grid());
        let grid = self.basin.grid();
        let mean_depth_m = self.params.mean_thermocline_depth_m();
        let amplitude = PACKET_AMPLITUDE_M / mean_depth_m;
        let envelope = |x_m: f64| gaussian_envelope(x_m, PACKET_CENTRE_X_M, self.packet_width_m);
        // `dE/dx` of that envelope, in m⁻¹.
        let envelope_slope_per_m = |x_m: f64| {
            -(x_m - PACKET_CENTRE_X_M) / (self.packet_width_m * self.packet_width_m) * envelope(x_m)
        };

        let (h_nx, h_ny) = grid.field_shape(H_STAGGERING);
        for j in 0..h_ny {
            let y_hat = self.basin.y_of_row_m(H_STAGGERING, j) / self.deformation_radius_m;
            let trapping = (-0.5 * y_hat * y_hat).exp();
            for i in 0..h_nx {
                let x_m = self.basin.x_of_column_m(H_STAGGERING, i);
                *state
                    .h_mut()
                    .get_mut(i, j)
                    .expect("the loop bounds are the field's own shape") = mean_depth_m
                    * amplitude
                    * envelope(x_m)
                    * (2.0 * y_hat * y_hat + 1.0)
                    * trapping;
            }
        }

        let (u_nx, u_ny) = grid.field_shape(U_STAGGERING);
        for j in 0..u_ny {
            let y_hat = self.basin.y_of_row_m(U_STAGGERING, j) / self.deformation_radius_m;
            let trapping = (-0.5 * y_hat * y_hat).exp();
            for i in 0..u_nx {
                let x_m = self.basin.x_of_column_m(U_STAGGERING, i);
                *state
                    .u_mut()
                    .get_mut(i, j)
                    .expect("the loop bounds are the field's own shape") = self
                    .wave_speed_m_per_s()
                    * amplitude
                    * envelope(x_m)
                    * (2.0 * y_hat * y_hat - 3.0)
                    * trapping;
            }
        }

        // `(8/3)·ŷ·e^{−ŷ²/2}` is `(4/3)·ψ₁`, and writing it that way is the
        // point: the mode is defined by its meridional velocity sitting on
        // `ψ₁` (module header), and every other field of it follows from that.
        let (v_nx, v_ny) = grid.field_shape(V_STAGGERING);
        for j in 0..v_ny {
            let waveguide = MeridionalStructure::First.at(
                self.basin.y_of_row_m(V_STAGGERING, j),
                self.deformation_radius_m,
            );
            for i in 0..v_nx {
                let x_m = self.basin.x_of_column_m(V_STAGGERING, i);
                *state
                    .v_mut()
                    .get_mut(i, j)
                    .expect("the loop bounds are the field's own shape") = self
                    .wave_speed_m_per_s()
                    * (4.0 / 3.0)
                    * amplitude
                    * self.deformation_radius_m
                    * envelope_slope_per_m(x_m)
                    * waveguide;
            }
        }

        state
    }

    /// Run the experiment: the packet in a closed, unforced, undamped basin,
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

    /// The measured zonal speed of the packet, in m/s — negative westward.
    fn measured_speed_m_per_s(self, propagation: &Propagation) -> f64 {
        let [early, late] = &propagation.samples;
        (self.energy_centroid_m(&late.westward) - self.energy_centroid_m(&early.westward))
            / (late.t_s - early.t_s)
    }

    /// How far that speed sits from the dispersion relation's prediction, as a
    /// fraction of the prediction.
    fn speed_error(self, propagation: &Propagation) -> f64 {
        (self.measured_speed_m_per_s(propagation) - self.predicted_speed_m_per_s()).abs()
            / self.predicted_speed_m_per_s().abs()
    }

    /// How far short of `c/3` the *measured* speed falls, as a fraction of
    /// `c/3` — the quantity the dispersion relation predicts as
    /// [`Experiment::predicted_slowdown`].
    fn measured_slowdown(self, propagation: &Propagation) -> f64 {
        1.0 - self.measured_speed_m_per_s(propagation) / self.long_wave_speed_m_per_s()
    }
}

/// One sample of a run: the two `ψ₀` invariants and the thermocline anomaly's
/// two meridional coefficients, column by column, at a given time.
#[derive(Debug, Clone)]
struct Sample {
    /// When it was taken, in seconds — the step's own time, not the requested
    /// one, so that a speed's denominator is exact.
    t_s: f64,
    /// `P₀[u/c + h/H]`, column by column: what, if anything, is going east.
    eastward: Vec<f64>,
    /// `P₀[u/c − h/H]`, column by column: where the Rossby packet is.
    westward: Vec<f64>,
    /// The `ψ₀` coefficient of `h/H`, column by column.
    depth_in_gravest: Vec<f64>,
    /// The `ψ₂` coefficient of `h/H`, column by column.
    depth_in_second: Vec<f64>,
}

/// What one run leaves for the tests to read: the packet early in its flight
/// and late in it.
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
    RUN.get_or_init(|| Experiment::new(COARSE_REFINEMENT, PACKET_WIDTH_M).run())
}

/// The same experiment with both cell dimensions refined by
/// [`FINE_REFINEMENT`] — the second resolution the acceptance criterion asks
/// for.
fn fine_run() -> &'static Propagation {
    static RUN: OnceLock<Propagation> = OnceLock::new();
    RUN.get_or_init(|| Experiment::new(FINE_REFINEMENT, PACKET_WIDTH_M).run())
}

/// The fine experiment with a packet [`NARROW_WIDTH_FACTOR`] times narrower —
/// the second spectrum the dispersion test compares.
///
/// At the *fine* resolution, and not the coarse one, because that test's budget
/// is the two runs' zonal truncations and those are what the refinement shrinks;
/// at the coarse resolution the bound would exceed the effect it is measuring.
fn narrow_run() -> &'static Propagation {
    static RUN: OnceLock<Propagation> = OnceLock::new();
    RUN.get_or_init(|| Experiment::new(FINE_REFINEMENT, PACKET_WIDTH_M / NARROW_WIDTH_FACTOR).run())
}

/// The two resolutions of the acceptance criterion, coarse first.
fn resolutions() -> [(Experiment, &'static Propagation); 2] {
    [
        (
            Experiment::new(COARSE_REFINEMENT, PACKET_WIDTH_M),
            coarse_run(),
        ),
        (Experiment::new(FINE_REFINEMENT, PACKET_WIDTH_M), fine_run()),
    ]
}

#[test]
fn the_rossby_packet_travels_west_at_a_third_of_the_kelvin_wave_speed() {
    // The headline claim of T-07.2, at more than one grid resolution and with
    // each resolution held to its own budget, so the finer run has to be better
    // rather than merely no worse: the gravest meridional Rossby mode travels
    // *westward* at approximately `c/3` (`CONTEXT.md`, *Rossby wave*).
    //
    // "Approximately" has a size and a sign here, and both come from theory.
    // The band is `c/3` widened by the run's truncation budget on the fast
    // side, and by that budget *plus the packet's own dispersive slowdown* on
    // the slow side: a packet of finite zonal width is slower than the
    // long-wave limit and never faster, so that term is one-sided.
    for (experiment, propagation) in resolutions() {
        let long_wave_m_per_s = experiment.long_wave_speed_m_per_s();
        let measured_m_per_s = experiment.measured_speed_m_per_s(propagation);
        let cells = experiment.basin.grid();

        assert!(
            measured_m_per_s < 0.0,
            "on the {}x{} grid the packet travelled at {measured_m_per_s} m/s, which is eastward: \
             the gravest equatorial Rossby mode goes west",
            cells.nx(),
            cells.ny()
        );

        let tolerance = experiment.speed_tolerance();
        let slowest_m_per_s =
            long_wave_m_per_s.abs() * (1.0 - tolerance - experiment.predicted_slowdown());
        let fastest_m_per_s = long_wave_m_per_s.abs() * (1.0 + tolerance);
        let measured_speed_m_per_s = measured_m_per_s.abs();
        assert!(
            (slowest_m_per_s..=fastest_m_per_s).contains(&measured_speed_m_per_s),
            "on the {}x{} grid the packet travelled west at {measured_speed_m_per_s} m/s, outside \
             the [{slowest_m_per_s}, {fastest_m_per_s}] m/s that grid's truncation budget of \
             {:.2}% and the packet's own {:.2}% dispersive slowdown allow around \
             c/3 = {:.6} m/s",
            cells.nx(),
            cells.ny(),
            100.0 * tolerance,
            100.0 * experiment.predicted_slowdown(),
            long_wave_m_per_s.abs()
        );

        // That band is wide enough to be worth pinning down further. The
        // analytic speeds either side of it are the `n = 2` Rossby mode at
        // `c/5` and the Kelvin branch at `c`, and the measured speed must be
        // nearer `c/3` than to either -- which is what identifies it as the
        // *gravest* mode rather than merely a westward signal of about the
        // right speed.
        let miss = (measured_speed_m_per_s - long_wave_m_per_s.abs()).abs();
        for (name, other_m_per_s) in [
            (
                "the n = 2 mode's c/5",
                experiment.wave_speed_m_per_s() / SECOND_ROSSBY_MERIDIONAL_TERM,
            ),
            ("the Kelvin speed c", experiment.wave_speed_m_per_s()),
        ] {
            assert!(
                miss < (measured_speed_m_per_s - other_m_per_s).abs(),
                "the packet's {measured_speed_m_per_s} m/s is nearer {name} ({other_m_per_s} m/s) \
                 than the gravest mode's c/3 = {:.6} m/s",
                long_wave_m_per_s.abs()
            );
        }
    }
}

#[test]
fn the_rossby_speed_matches_the_equatorial_dispersion_relation() {
    // The deliverable's second half: not merely "about `c/3`", but the speed
    // the equatorial Rossby dispersion relation gives for *this* packet.
    // `mean_group_speed_in_c` computes it from the cubic and the packet's own
    // Gaussian spectrum, so the expected value here comes from theory and not
    // from the engine (CODING_STANDARDS.md, § *Tests*).
    for (experiment, propagation) in resolutions() {
        let predicted_m_per_s = experiment.predicted_speed_m_per_s();
        let measured_m_per_s = experiment.measured_speed_m_per_s(propagation);
        let error = experiment.speed_error(propagation);
        let tolerance = experiment.speed_tolerance();
        let cells = experiment.basin.grid();

        assert!(
            error <= tolerance,
            "on the {}x{} grid the packet travelled at {measured_m_per_s} m/s, {:.2}% from the \
             {predicted_m_per_s} m/s the dispersion relation gives a packet {:.2} deformation \
             radii wide; that grid's truncation budget allows {:.2}%",
            cells.nx(),
            cells.ny(),
            100.0 * error,
            experiment.width_in_radii(),
            100.0 * tolerance
        );
    }
}

#[test]
fn the_rossby_speed_error_shrinks_at_the_schemes_second_order() {
    // The other half of T-07.1's acceptance criterion, carried over: the error
    // shrinks with resolution rather than sitting at a fixed offset. Both
    // truncation terms of the speed budget are second order in the cell size
    // (module header) and both cell dimensions are halved, so the error must
    // fall by about four (CODING_STANDARDS.md, § *Convergence over point
    // checks*).
    //
    // The error measured is the one against the *dispersion relation*, not
    // against `c/3`: the gap to `c/3` is the packet's own physics and would
    // survive any refinement, so a rate measured against it would be reading a
    // floor rather than the scheme.
    let [(coarse, coarse_run), (fine, fine_run)] = resolutions();

    let coarse_error = coarse.speed_error(coarse_run);
    let fine_error = fine.speed_error(fine_run);
    assert!(
        fine_error > 0.0,
        "the fine run reproduced the predicted speed to the last bit, so there is no error left \
         to measure a rate on: what this run reads is the measurement's floor and not the \
         scheme's order"
    );
    let order = (coarse_error / fine_error).log2();

    // Bounded on both sides: too small an order is a scheme that is not second
    // order, and too large a one is a fine error that is no longer the
    // truncation -- either way the point check's budget has stopped describing
    // what the run does.
    assert!(
        (MIN_CONVERGENCE_ORDER..=MAX_CONVERGENCE_ORDER).contains(&order),
        "refining the cells by {FINE_REFINEMENT} took the Rossby speed's error from {:.3}% to \
         {:.3}%, a convergence order of {order:.2}: outside the \
         [{MIN_CONVERGENCE_ORDER}, {MAX_CONVERGENCE_ORDER}] a second-order scheme owes",
        100.0 * coarse_error,
        100.0 * fine_error
    );
}

#[test]
fn narrowing_the_packet_slows_it_by_what_the_dispersion_relation_says() {
    // The sharp statement of dispersion, and the mirror image of T-07.1's
    // `the_kelvin_speed_does_not_depend_on_the_packets_zonal_width`: on the
    // Kelvin branch every wavenumber travels at `c`, so a narrower packet
    // travels at the same speed; on this branch the group speed depends on the
    // wavenumber, so a narrower packet is *measurably slower*, by an amount the
    // dispersion relation names in advance.
    //
    // What is compared is the difference between the two runs' slowdowns, and
    // not each slowdown on its own. The two runs share a grid, a start
    // position and a clock, so the meridional truncation -- the dominant term
    // of either run's own budget, and by far the largest thing separating a
    // measured speed from a predicted one -- is a relative error common to
    // both and cancels from the difference to the extent that the two speeds
    // are equal. What is left is the terms below.
    let wide = Experiment::new(FINE_REFINEMENT, PACKET_WIDTH_M);
    let narrow = Experiment::new(FINE_REFINEMENT, PACKET_WIDTH_M / NARROW_WIDTH_FACTOR);

    let predicted = narrow.predicted_slowdown() - wide.predicted_slowdown();
    let measured = narrow.measured_slowdown(narrow_run()) - wide.measured_slowdown(fine_run());

    // The residue of the meridional truncation, which scales both speeds by
    // the same relative factor and so moves their difference by that factor
    // times the difference itself; the two zonal truncations, which differ
    // because they go as `σ⁻²` and the widths are what differ; and the two
    // stray-energy terms, for the same reason. Times `TRUNCATION_SAFETY`, as
    // everywhere else in this file.
    let bound = TRUNCATION_SAFETY
        * (wide.meridional_truncation() * predicted
            + wide.zonal_truncation()
            + narrow.zonal_truncation()
            + wide.stray_energy()
            + narrow.stray_energy());

    assert!(
        (measured - predicted).abs() <= bound,
        "narrowing the packet by {NARROW_WIDTH_FACTOR} slowed it by {:.3}% of c/3, where the \
         dispersion relation predicts {:.3}% ({:.3}% at the wider width, {:.3}% at the narrower); \
         the budget allows a {:.3}% discrepancy",
        100.0 * measured,
        100.0 * predicted,
        100.0 * wide.predicted_slowdown(),
        100.0 * narrow.predicted_slowdown(),
        100.0 * bound
    );
}

#[test]
fn the_rossby_packet_carries_no_eastward_energy() {
    // "Westward" (`CONTEXT.md`, *Rossby wave*), read through the decomposition
    // that can tell the two directions apart. The gravest Rossby mode's
    // eastward invariant is `ψ₂` and nothing else, at every wavenumber (module
    // header), so `P₀[u/c + h/H]` is empty for the whole flight and anything in
    // it is either a Kelvin wave the run should not contain or the
    // discretisation's leakage.
    for (experiment, propagation) in resolutions() {
        // The `5·(Δy/Le)²` of the speed budget is an amplitude, so it is that
        // squared in energy, times `TRUNCATION_SAFETY` for the coefficient it
        // does not evaluate.
        let leakage = experiment.meridional_truncation();
        let ceiling = TRUNCATION_SAFETY * leakage * leakage;

        for sample in &propagation.samples {
            let eastward = experiment.energy(&sample.eastward);
            let westward = experiment.energy(&sample.westward);
            let share = eastward / (eastward + westward);
            assert!(
                share <= ceiling,
                "{:.0} s into the run the eastward mode held {:.4}% of the energy, past the {:.4}% \
                 the meridional discretisation can leak: the packet is not travelling west only",
                sample.t_s,
                100.0 * share,
                100.0 * ceiling
            );
        }
    }
}

#[test]
fn the_rossby_packet_keeps_the_gravest_modes_meridional_shape() {
    // Speed alone does not identify a mode, so this is the other half of the
    // claim: the packet's thermocline anomaly still has the `2·ψ₀ + ψ₂/2`
    // structure of the gravest Rossby mode at the end of its flight -- deepest
    // off the equator rather than on it, which is what tells it from a Kelvin
    // wave.
    for (experiment, propagation) in resolutions() {
        let late = &propagation.samples[1];

        // Read the shape where the packet is, so that whatever the run has left
        // elsewhere in the basin stays out of the window.
        let centre_m = experiment.energy_centroid_m(&late.westward);
        let window_m = SHAPE_WINDOW_IN_WIDTHS * experiment.packet_width_m;
        let (cross_term, gravest_energy) = (0..experiment.basin.grid().nx())
            .filter(|i| (experiment.column_x_m(*i) - centre_m).abs() <= window_m)
            .fold((0.0, 0.0), |(cross, gravest), i| {
                (
                    cross + late.depth_in_gravest[i] * late.depth_in_second[i],
                    gravest + late.depth_in_gravest[i] * late.depth_in_gravest[i],
                )
            });
        let ratio = cross_term / gravest_energy;

        // Two terms. The meridional truncation, which is what a second-order
        // scheme can do to the `ψ₂` half of the structure; and the `O(k̂²)`
        // correction to the modal structure itself at the packet's own
        // wavenumber, the initial condition being the long-wave mode. Times
        // `TRUNCATION_SAFETY` for the coefficients neither evaluates.
        let tolerance = TRUNCATION_SAFETY
            * (experiment.meridional_truncation() + experiment.mean_square_wavenumber());
        let error = (ratio - GRAVEST_ROSSBY_SHAPE_RATIO).abs() / GRAVEST_ROSSBY_SHAPE_RATIO;
        assert!(
            error <= tolerance,
            "at the end of its flight the packet's thermocline anomaly carried a psi2/psi0 ratio \
             of {ratio}, {:.1}% from the {GRAVEST_ROSSBY_SHAPE_RATIO} of the gravest Rossby mode; \
             the budget allows {:.1}%",
            100.0 * error,
            100.0 * tolerance
        );
    }
}
