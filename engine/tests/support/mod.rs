//! The equatorial-wave analysis kit the scientific-validation tests share.
//!
//! `tests/*.rs` are separate crates, so anything more than one of them needs is
//! otherwise one copy per test binary, drifting apart. Four of them —
//! `western_boundary_reflection.rs`, `eastern_boundary_reflection.rs`,
//! `kelvin_wave_propagation.rs` and `rossby_wave_dispersion.rs` — measure the
//! same three things in the same way, and this module is the one definition of
//! each:
//!
//! - **The equatorial-Pacific ocean** the validations are run in, and the two
//!   speeds written out from theory rather than asked of the engine:
//!   `c = √(g'·H)` and `Le = √(c/β)`.
//! - **The meridional decomposition.** The linear equatorial equations separate
//!   on the parabolic cylinder functions `ψₘ(ŷ) = Hₘ(ŷ)·e^{−ŷ²/2}`, and in the
//!   two invariants
//!
//!   ```text
//!   eastward = u/c + h/H          westward = u/c − h/H
//!   ```
//!
//!   each wave sits on one `ψₘ` of one invariant: the Kelvin wave is `ψ₀` of
//!   the eastward one and nothing else, and the gravest Rossby mode is `ψ₂` of
//!   the eastward one and `ψ₀` of the westward one (Matsuno 1966; Gill,
//!   *Atmosphere–Ocean Dynamics*, § 11.6). Projecting a run onto a `ψₘ` is
//!   therefore how a test isolates one wave from another, and
//!   [`MeridionalStructure`], [`Waveguide`] and [`project_columns`] are how it
//!   is done.
//! - **How a wave's position and speed are read**: the energy-weighted zonal
//!   centroid of a projected profile ([`energy_centroid_m`]), which moves at the
//!   packet's energy-weighted mean group velocity, or the time of flight of its
//!   peak between two fixed zonal stations ([`peak_time_s`]).
//! - **The initial conditions themselves**: [`kelvin_pulse_state`] and
//!   [`gravest_rossby_packet_state`], the two analytic packets every run in
//!   this suite starts from.
//! - **How a wave's meridional decay scale is read**: the zonal sum of one
//!   invariant, row by row ([`invariant_meridional_profile`]), fitted against
//!   the `ψₘ` the theory says it is with the trapping scale left free
//!   ([`fitted_trapping_scale_m`]). The projections above hold the scale fixed
//!   at `Le` and read an amplitude; this reads the scale itself.
//!
//! Nothing here asserts anything or carries a tolerance: a tolerance is a
//! property of one run's configuration and belongs in the test that states that
//! configuration. What this module provides is the machinery a tolerance is
//! *computed with* — [`MeridionalStructure::truncation_richness`] and
//! [`Waveguide::truncation_bound`] — and the measurement itself.
//!
//! `common/mod.rs` is the other shared test module; it holds what is not
//! physics (scratch directories for tests that write real runs).

// Each test binary uses the part of the kit its own wave needs, so the
// remainder is dead code in that binary and live in the next one.
#![allow(dead_code)]

use engine::{Basin, OceanState, PhysicalParams, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
pub const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
pub const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
pub const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;
/// Reference seawater density `ρ₀`, in kg/m³ — `CONTEXT.md` and Gill, appendix 3.
pub const REFERENCE_DENSITY_KG_PER_M3: f64 = engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3;
/// Rayleigh damping `r` of a wave-speed validation, in s⁻¹.
///
/// Zero. Every analytic speed these tests compare against is the *undamped*
/// one, and damping decays a packet's amplitude without moving it; running
/// damped would only put a decay factor between the two samples of a
/// measurement that reads positions.
pub const UNDAMPED_PER_S: f64 = 0.0;

/// The equatorial-Pacific parameter set the validations run in.
///
/// # Panics
/// If the published parameters above are rejected as unphysical, which would
/// mean the engine's validation is wrong rather than the values.
pub fn pacific_params() -> PhysicalParams {
    params_with_reduced_gravity_and_beta(PACIFIC_REDUCED_GRAVITY_M_PER_S2, BETA_PER_M_PER_S)
}

/// The same ocean with `g'` and `β` replaced, and `H`, `ρ₀` and the undamped
/// `r` left at the Pacific values above.
///
/// `Le = √(c/β) = √(√(g'H)/β)` is a *prediction*, and a validation that only
/// ever ran in one ocean could not tell that prediction apart from a fixed
/// length that happened to fit. These are the only two parameters `Le` depends
/// on, so varying them — and nothing else — moves the predicted radius while
/// leaving the basin, the packet and the measurement alone.
///
/// # Panics
/// If the pair is rejected as unphysical, which is the caller asking for an
/// ocean that does not exist rather than a fault in the engine.
pub fn params_with_reduced_gravity_and_beta(
    reduced_gravity_m_per_s2: f64,
    beta_per_m_per_s: f64,
) -> PhysicalParams {
    PhysicalParams::new(
        reduced_gravity_m_per_s2,
        PACIFIC_MEAN_DEPTH_M,
        UNDAMPED_PER_S,
        beta_per_m_per_s,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the caller's reduced gravity and beta are physical")
}

/// Kelvin wave speed `c = √(g'·H)`, in m/s, written out from the definition in
/// `CONTEXT.md` rather than asked of the code under test.
///
/// It is the value the assertions are made against, so it must not come from
/// the engine: an expected value read out of the thing being measured agrees
/// with it by construction (CODING_STANDARDS.md § *Tests*). The engine computes
/// the same number from the same two parameters, which is what makes the
/// comparison a test rather than a tautology.
pub fn kelvin_wave_speed_m_per_s() -> f64 {
    wave_speed_of_m_per_s(PACIFIC_REDUCED_GRAVITY_M_PER_S2)
}

/// The same `c = √(g'·H)` for an ocean of reduced gravity
/// `reduced_gravity_m_per_s2`, in m/s.
pub fn wave_speed_of_m_per_s(reduced_gravity_m_per_s2: f64) -> f64 {
    (reduced_gravity_m_per_s2 * PACIFIC_MEAN_DEPTH_M).sqrt()
}

/// Equatorial deformation radius `Le = √(c/β)`, in metres — the meridional
/// scale of the waveguide (`CONTEXT.md`), written out from the same definition
/// and for the same reason.
pub fn equatorial_deformation_radius_m() -> f64 {
    deformation_radius_of_m(kelvin_wave_speed_m_per_s(), BETA_PER_M_PER_S)
}

/// The same `Le = √(c/β)` for an ocean of wave speed `wave_speed_m_per_s` and
/// beta-plane gradient `beta_per_m_per_s`, in metres.
pub fn deformation_radius_of_m(wave_speed_m_per_s: f64, beta_per_m_per_s: f64) -> f64 {
    (wave_speed_m_per_s / beta_per_m_per_s).sqrt()
}

/// A Gaussian zonal envelope of e-folding half-width `width_m`, centred on
/// `centre_m`.
///
/// The zonal profile every packet in these tests is built from: it is smooth,
/// its spectrum is a Gaussian of known moments — `⟨k²⟩ = 1/(2σ²)` in energy —
/// and those moments are what the zonal truncation and dispersive-bias terms of
/// the tolerances are computed from.
pub fn gaussian_envelope(x_m: f64, centre_m: f64, width_m: f64) -> f64 {
    let offset = (x_m - centre_m) / width_m;
    (-0.5 * offset * offset).exp()
}

/// One of the parabolic cylinder functions the equatorial waveguide separates
/// on.
///
/// Named by how gravely they oscillate rather than by an index, because the
/// index is not the interesting thing about any of them:
/// [`MeridionalStructure::Gravest`] is the Kelvin wave's whole shape and the
/// gravest Rossby mode's leading one, [`MeridionalStructure::First`] is the
/// meridional velocity of that Rossby mode, and
/// [`MeridionalStructure::Second`] is the partner that makes its thermocline
/// anomaly deepest off the equator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeridionalStructure {
    /// `ψ₀(ŷ) = e^{−ŷ²/2}`.
    Gravest,
    /// `ψ₁(ŷ) = 2ŷ·e^{−ŷ²/2}`.
    First,
    /// `ψ₂(ŷ) = (4ŷ² − 2)·e^{−ŷ²/2}`.
    Second,
}

impl MeridionalStructure {
    /// Index `m` of the Hermite function this structure is.
    pub const fn hermite_order(self) -> usize {
        match self {
            Self::Gravest => 0,
            Self::First => 1,
            Self::Second => 2,
        }
    }

    /// The physicists' Hermite polynomial `Hₘ(ŷ)` of that order — the
    /// polynomial the equatorial wave problem separates in.
    pub fn hermite_polynomial(self, y_over_le: f64) -> f64 {
        match self {
            Self::Gravest => 1.0,
            Self::First => 2.0 * y_over_le,
            Self::Second => 4.0 * y_over_le * y_over_le - 2.0,
        }
    }

    /// The Hermite function `ψₘ(ŷ) = Hₘ(ŷ)·exp(−ŷ²/2)`.
    pub fn hermite_function(self, y_over_le: f64) -> f64 {
        self.hermite_polynomial(y_over_le) * (-0.5 * y_over_le * y_over_le).exp()
    }

    /// The same function at `y_m` metres from the equator, for a waveguide of
    /// the given deformation radius.
    pub fn at(self, y_m: f64, deformation_radius_m: f64) -> f64 {
        self.hermite_function(y_m / deformation_radius_m)
    }

    /// `∫ψₘ² dŷ = 2ᵐ·m!·√π`, the normalisation an analytic projection divides
    /// by.
    ///
    /// Written as the product `∏ₖ₌₁..ₘ 2k`, which is that value, so the orders
    /// share one definition instead of carrying a tabulated constant each.
    pub fn hermite_norm(self) -> f64 {
        let weight: f64 = (1..=self.hermite_order()).map(|k| 2.0 * k as f64).product();
        weight * std::f64::consts::PI.sqrt()
    }

    /// How much finer this structure is than the deformation radius, as the
    /// factor `2m + 1` multiplying `(Δy/Le)²` in a truncation bound.
    ///
    /// `ψₘ` oscillates on the scale `Le/√(2m+1)` — its classical turning point
    /// is at `ŷ = √(2m+1)` — so a second-order scheme's truncation error on it
    /// grows with the square of that local wavenumber. `ψ₀` gives one, which is
    /// why a Kelvin bound is `(Δy/Le)²` with no factor; `ψ₂` gives five, and
    /// taking one there would be optimistic rather than conservative.
    pub fn truncation_richness(self) -> f64 {
        (2 * self.hermite_order() + 1) as f64
    }
}

/// The equatorial waveguide of one run, as far as a meridional projection needs
/// to know it: the scale the Hermite functions are stretched by, and the rows a
/// column is sampled on.
pub struct Waveguide {
    /// Equatorial deformation radius `Le = √(c/β)`, in metres (`CONTEXT.md`).
    pub le_m: f64,
    /// Cell height, in metres — the quadrature weight of one row.
    pub dy_m: f64,
    /// Meridional positions of the cell-centre rows, in metres north of the
    /// equator.
    pub row_y_m: Vec<f64>,
}

impl Waveguide {
    /// The waveguide of `basin`, sampled on its cell-centre rows.
    pub fn new(basin: Basin, params: PhysicalParams) -> Self {
        Self {
            le_m: (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt(),
            dy_m: basin.spacing().dy_m(),
            row_y_m: row_positions_m(basin),
        }
    }

    /// The coefficient of `ψₘ` in one meridional column, by discrete quadrature
    /// of `∫column·ψₘ dŷ / ∫ψₘ² dŷ`.
    ///
    /// The midpoint rule on the cell-centre rows, which is second order and
    /// symmetric about the equator — the same order as the scheme whose output
    /// it reads. Because it divides by the *analytic* `∫ψₘ² dŷ`, the number it
    /// returns is on the same scale as an amplitude written down from theory,
    /// which is what lets a test compare one against the other.
    pub fn coefficient(&self, column: &[f64], structure: MeridionalStructure) -> f64 {
        let integral: f64 = column
            .iter()
            .zip(&self.row_y_m)
            .map(|(value, y_m)| value * structure.hermite_function(y_m / self.le_m))
            .sum();
        integral * (self.dy_m / self.le_m) / structure.hermite_norm()
    }

    /// The second-order meridional truncation bound for `structure`, as a
    /// fraction of a wave speed: `(2m + 1)·(Δy/Le)²`, with the remaining `O(1)`
    /// constant taken as one.
    pub fn truncation_bound(&self, structure: MeridionalStructure) -> f64 {
        structure.truncation_richness() * (self.dy_m / self.le_m).powi(2)
    }
}

/// Meridional positions of `basin`'s cell-centre rows, in metres north of the
/// equator — where `h` sits, where `u` is averaged to, and the abscissae of
/// every meridional quadrature and fit below.
pub fn row_positions_m(basin: Basin) -> Vec<f64> {
    (0..basin.grid().ny())
        .map(|j| basin.y_of_row_m(H_STAGGERING, j))
        .collect()
}

/// `ψₘ` sampled on the cell-centre rows, which is where `h` sits and where `u`
/// is averaged to.
pub fn row_structure(
    basin: Basin,
    deformation_radius_m: f64,
    structure: MeridionalStructure,
) -> Vec<f64> {
    (0..basin.grid().ny())
        .map(|j| structure.at(basin.y_of_row_m(H_STAGGERING, j), deformation_radius_m))
        .collect()
}

/// The `ψₘ` coefficient, column by column, of a cell-centred field given as
/// `value(i, j)`.
///
/// The `ψₘ` are orthogonal on the line and the basins these tests use reach
/// several deformation radii, so this discrete inner product is the modal
/// coefficient to the accuracy of the row quadrature; the row spacing cancels
/// between the projection and its normalisation. Unlike
/// [`Waveguide::coefficient`] it divides by the *discrete* `Σψₘ²`, which makes
/// it the least-squares coefficient of the rows actually present — the right
/// normalisation when what is compared is one projection against another rather
/// than a projection against an analytic amplitude.
pub fn project_columns(
    basin: Basin,
    deformation_radius_m: f64,
    structure: MeridionalStructure,
    value: impl Fn(usize, usize) -> f64,
) -> Vec<f64> {
    let weights = row_structure(basin, deformation_radius_m, structure);
    let normalisation: f64 = weights.iter().map(|weight| weight * weight).sum();
    (0..basin.grid().nx())
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

/// The thermocline depth anomaly's `ψₘ` coefficient, column by column, in units
/// of the mean depth `H`.
pub fn depth_projection(
    basin: Basin,
    deformation_radius_m: f64,
    params: PhysicalParams,
    state: &OceanState,
    structure: MeridionalStructure,
) -> Vec<f64> {
    project_columns(basin, deformation_radius_m, structure, |i, j| {
        state.h().get(i, j).expect("a cell centre") / params.mean_thermocline_depth_m()
    })
}

/// The `ψ₀` coefficient of the zonal current, column by column, in units of
/// `c`.
///
/// `u` is averaged from its two faces onto the cell centre so that it and the
/// depth anomaly are read at one set of positions, which is what makes their
/// sum and difference the invariants below.
pub fn gravest_current_projection(
    basin: Basin,
    deformation_radius_m: f64,
    state: &OceanState,
    wave_speed_m_per_s: f64,
) -> Vec<f64> {
    project_columns(
        basin,
        deformation_radius_m,
        MeridionalStructure::Gravest,
        |i, j| {
            let west = state.u().get(i, j).expect("an east/west face");
            let east = state.u().get(i + 1, j).expect("an east/west face");
            0.5 * (west + east) / wave_speed_m_per_s
        },
    )
}

/// The two `ψ₀` invariants of a state, column by column.
#[derive(Debug, Clone)]
pub struct Invariants {
    /// `P₀[u/c + h/H]`: where the Kelvin wave is, and how much of it there is.
    pub eastward: Vec<f64>,
    /// `P₀[u/c − h/H]`: the same for the gravest Rossby mode.
    pub westward: Vec<f64>,
    /// The `ψ₀` coefficient of `h/H` alone, column by column — the half of the
    /// two invariants a meridional-shape test reads on its own.
    pub depth_in_gravest: Vec<f64>,
}

/// The two `ψ₀` invariants of `state`, and the depth projection they are built
/// from.
pub fn invariants(
    basin: Basin,
    deformation_radius_m: f64,
    params: PhysicalParams,
    state: &OceanState,
    wave_speed_m_per_s: f64,
) -> Invariants {
    let current =
        gravest_current_projection(basin, deformation_radius_m, state, wave_speed_m_per_s);
    let depth_in_gravest = depth_projection(
        basin,
        deformation_radius_m,
        params,
        state,
        MeridionalStructure::Gravest,
    );
    Invariants {
        eastward: current
            .iter()
            .zip(&depth_in_gravest)
            .map(|(current, depth)| current + depth)
            .collect(),
        westward: current
            .iter()
            .zip(&depth_in_gravest)
            .map(|(current, depth)| current - depth)
            .collect(),
        depth_in_gravest,
    }
}

/// `Σ amplitude²` over whatever columns are handed in — the weight a centroid
/// divides by, and the quantity an eastward/westward split is stated in.
pub fn energy(amplitudes: impl Iterator<Item = f64>) -> f64 {
    amplitudes.map(|amplitude| amplitude * amplitude).sum()
}

/// The energy-weighted zonal centroid of a profile given as
/// `(position in metres, amplitude)` pairs.
///
/// The centroid of a linear wave packet moves at its energy-weighted mean group
/// velocity, which is what a centroid speed measurement reads.
///
/// # Panics
/// If the profile carries no energy at all, which would mean the run never had
/// a wave in it.
pub fn energy_centroid_m(profile: impl Iterator<Item = (f64, f64)>) -> f64 {
    let (moment, weight) = profile.fold((0.0, 0.0), |(moment, weight), (x_m, amplitude)| {
        let energy = amplitude * amplitude;
        (moment + x_m * energy, weight + energy)
    });
    assert!(
        weight > 0.0,
        "the profile carries no energy, so it has no centroid"
    );
    moment / weight
}

/// The RMS zonal width of `profile` about its own centroid, in metres, over the
/// columns within `window_m` of it.
///
/// The window travels with the packet, so the same fraction of the same wave is
/// weighed at every sample time and a comparison between two of them is a
/// statement about the packet's shape rather than about where it is.
///
/// It is an *energy*-weighted RMS, so for a Gaussian amplitude profile of width
/// `σ` it reads `σ/√2` rather than `σ`; only the ratio of two of them is ever
/// asserted, so that factor cancels.
///
/// # Panics
/// If the window holds no energy, which would mean the packet is not where its
/// centroid says.
pub fn rms_width_m(profile: &[f64], column_x_m: impl Fn(usize) -> f64, window_m: f64) -> f64 {
    let centre_m = energy_centroid_m(
        profile
            .iter()
            .enumerate()
            .map(|(i, amplitude)| (column_x_m(i), *amplitude)),
    );
    let (moment, weight) = profile
        .iter()
        .enumerate()
        .map(|(i, amplitude)| (column_x_m(i) - centre_m, amplitude * amplitude))
        .filter(|(offset_m, _)| offset_m.abs() <= window_m)
        .fold((0.0, 0.0), |(moment, weight), (offset_m, energy)| {
            (moment + offset_m * offset_m * energy, weight + energy)
        });
    assert!(
        weight > 0.0,
        "the shape window holds no energy, so the packet is not where its centroid says"
    );
    (moment / weight).sqrt()
}

/// Index of the cell-centre column whose position is closest to `x_m`.
pub fn column_nearest(x_m: f64, cell_m: f64) -> usize {
    (x_m / cell_m - 0.5).round().max(0.0) as usize
}

/// Index of the largest value in `values`.
///
/// # Panics
/// If `values` is empty, or contains a NaN.
pub fn peak_index(values: &[f64]) -> usize {
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

/// The time, in seconds, at which `series` peaks, refined below the sampling
/// interval by fitting a parabola through the largest sample and its two
/// neighbours.
///
/// # Panics
/// If the peak sits on either end of the record, which would mean the run is
/// too short or a station is misplaced rather than that the wave is slow.
pub fn peak_time_s(series: &[f64], dt_s: f64) -> f64 {
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

/// Which of the two invariants `u/c ± h/H` a measurement reads.
///
/// The module header states which wave lives in which: the Kelvin wave is the
/// whole of [`Invariant::Eastward`], and the gravest Rossby mode puts its `ψ₀`
/// content in [`Invariant::Westward`] and its `ψ₂` content in the eastward one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Invariant {
    /// `u/c + h/H`.
    Eastward,
    /// `u/c − h/H`.
    Westward,
}

impl Invariant {
    /// The invariant's value from the two scaled fields at one cell centre.
    pub fn of(self, current_in_c: f64, depth_in_mean_depths: f64) -> f64 {
        match self {
            Self::Eastward => current_in_c + depth_in_mean_depths,
            Self::Westward => current_in_c - depth_in_mean_depths,
        }
    }
}

/// The meridional profile of one invariant: its zonal sum, row by row.
///
/// [`project_columns`] and the functions built on it collapse the *meridional*
/// axis and keep the zonal one, because what they read is where a wave is. A
/// decay-scale fit needs the opposite, so this collapses the zonal axis
/// instead. For a separable packet `A·E(x)·ψₘ(y/Le)` the sum is
/// `A·(ΣᵢE)·ψₘ(y/Le)`: the same meridional shape, with the envelope reduced to
/// one constant that the fit's free amplitude absorbs.
///
/// Summing over the whole basin rather than a window around the packet is
/// deliberate — a window needs a centre, and taking that from theory is the
/// circularity a measurement must not have. The cost is that whatever the run
/// has shed elsewhere is weighed too, which is the leading term of the fit's
/// error budget rather than a neglected one.
///
/// `u` is averaged from its two faces onto the cell centre, so it and the
/// depth anomaly are read at one set of positions — the same averaging
/// [`gravest_current_projection`] makes, and what makes their sum and
/// difference the invariants at all.
pub fn invariant_meridional_profile(
    basin: Basin,
    params: PhysicalParams,
    state: &OceanState,
    wave_speed_m_per_s: f64,
    invariant: Invariant,
) -> Vec<f64> {
    (0..basin.grid().ny())
        .map(|j| {
            (0..basin.grid().nx())
                .map(|i| {
                    let west = state.u().get(i, j).expect("an east/west face");
                    let east = state.u().get(i + 1, j).expect("an east/west face");
                    let current_in_c = 0.5 * (west + east) / wave_speed_m_per_s;
                    let depth_in_mean_depths = state.h().get(i, j).expect("a cell centre")
                        / params.mean_thermocline_depth_m();
                    invariant.of(current_in_c, depth_in_mean_depths)
                })
                .sum()
        })
        .collect()
}

/// How well `ψₘ(y/scale_m)` explains `profile`, as a number in `[0, 1]`.
///
/// The squared cosine between the profile and the model,
/// `(Σd·m)² / (Σd²·Σm²)`, which is the fraction of the profile's energy a
/// least-squares fit of the model to it accounts for. The model's amplitude
/// does not appear: for a one-parameter family scaled by a free amplitude, the
/// best amplitude is `Σd·m / Σm²` and substituting it leaves exactly this. So
/// maximising this over `scale_m` *is* the least-squares fit, with the
/// amplitude eliminated analytically rather than searched for.
///
/// # Panics
/// If the profile carries no energy, which would mean the run has no wave in
/// it.
pub fn shape_correlation(
    row_y_m: &[f64],
    profile: &[f64],
    structure: MeridionalStructure,
    scale_m: f64,
) -> f64 {
    let (cross, model_energy, profile_energy) = row_y_m
        .iter()
        .zip(profile)
        .map(|(y_m, value)| (structure.at(*y_m, scale_m), *value))
        .fold((0.0, 0.0, 0.0), |(cross, model, data), (shape, value)| {
            (
                cross + shape * value,
                model + shape * shape,
                data + value * value,
            )
        });
    assert!(
        profile_energy > 0.0,
        "the meridional profile carries no energy, so it has no shape to fit"
    );
    cross * cross / (model_energy * profile_energy)
}

/// How many points [`fitted_trapping_scale_m`] samples its bracket at before
/// refining.
///
/// The objective is unimodal for the `ψ₀` fit — for two Gaussians of scales
/// `L` and `Le` it is `4LLe/(L+Le)²`, whose only stationary point is `L = Le` —
/// but a `ψₘ` with nodes can put a second, lower hump in it where the model's
/// lobes fall between the profile's. A scan first makes the refinement start in
/// the right hump instead of trusting unimodality; 64 points over a bracket
/// spanning a factor of nine resolve the humps of every structure used here,
/// which are `O(Le)` wide, many times over.
const TRAPPING_SCALE_SCAN_POINTS: usize = 64;
/// Relative width the golden-section refinement narrows the bracket to.
///
/// `10⁻⁶`, which is four orders of magnitude below the smallest term of any
/// decay-scale budget these tests state, so the fit's own convergence never
/// appears in a comparison. It is not taken to machine precision because the
/// objective is quadratic at its maximum, where rounding flattens it below
/// about `√ε ≈ 10⁻⁸` of the scale.
const TRAPPING_SCALE_PRECISION: f64 = 1.0e-6;

/// The trapping scale, in metres, that best fits `profile` as `ψₘ(y/scale)`.
///
/// The measurement the deformation-radius validation is built on: the theory
/// says the profile is `ψₘ(y/Le)` with `Le = √(c/β)`, so leaving the scale free
/// and fitting it recovers `Le` from the run, to be compared against the
/// analytic one. The amplitude is eliminated analytically
/// ([`shape_correlation`]) and only the scale is searched for: a coarse scan of
/// `bracket_m` to pick the right hump of the objective, then golden-section
/// refinement within the two scan points either side of the best.
///
/// # Panics
/// If the best scan point is an endpoint of `bracket_m`, which means the
/// profile's scale is outside the bracket — the run is not carrying the wave
/// the caller thinks it is, rather than the fit having failed to converge.
pub fn fitted_trapping_scale_m(
    row_y_m: &[f64],
    profile: &[f64],
    structure: MeridionalStructure,
    bracket_m: (f64, f64),
) -> f64 {
    let (smallest_m, largest_m) = bracket_m;
    assert!(
        0.0 < smallest_m && smallest_m < largest_m,
        "a trapping-scale bracket runs from a positive scale up to a larger one, not \
         ({smallest_m}, {largest_m})"
    );
    let step_m = (largest_m - smallest_m) / (TRAPPING_SCALE_SCAN_POINTS - 1) as f64;
    let scan_at = |point: usize| smallest_m + point as f64 * step_m;
    let objective = |scale_m: f64| shape_correlation(row_y_m, profile, structure, scale_m);

    let best = peak_index(
        &(0..TRAPPING_SCALE_SCAN_POINTS)
            .map(|point| objective(scan_at(point)))
            .collect::<Vec<_>>(),
    );
    assert!(
        best > 0 && best + 1 < TRAPPING_SCALE_SCAN_POINTS,
        "the profile is best fitted at the {}edge of the bracket [{smallest_m} m, {largest_m} m], \
         so its meridional scale lies outside it",
        if best == 0 { "lower " } else { "upper " }
    );

    // Golden-section search for the maximum inside the bracketing triple the
    // scan found. The interval shrinks by the golden ratio each iteration, so
    // the loop terminates; it is written as a `while` on the width rather than
    // a fixed count so the stopping condition is the stated precision itself.
    let golden = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut low_m, mut high_m) = (scan_at(best - 1), scan_at(best + 1));
    while high_m - low_m > TRAPPING_SCALE_PRECISION * high_m {
        let span_m = high_m - low_m;
        let (left_m, right_m) = (high_m - golden * span_m, low_m + golden * span_m);
        if objective(left_m) > objective(right_m) {
            high_m = right_m;
        } else {
            low_m = left_m;
        }
    }
    0.5 * (low_m + high_m)
}

/// The Gaussian wave packet a validation run is started from.
///
/// The zonal half of an initial condition: how tall it is, where it is, and how
/// wide. The meridional half is fixed by which wave the packet is, and belongs
/// to the state constructor rather than to this.
#[derive(Debug, Clone, Copy)]
pub struct Packet {
    /// Scale of the packet's thermocline depth anomaly, in metres.
    ///
    /// The core is linear, so every speed, width and ratio a validation reads
    /// is independent of this; it only sets the units the diagnostics come out
    /// in.
    pub amplitude_m: f64,
    /// Zonal position of the packet's centre at `t = 0`, in metres east of the
    /// western wall.
    pub centre_x_m: f64,
    /// Zonal e-folding half-width `σ` of the Gaussian envelope, in metres.
    pub width_m: f64,
}

/// An equatorial Kelvin pulse: Gaussian in `x` on the `ψ₀` waveguide, with
/// `u = (c/H)·h` and `v = 0`.
///
/// An exact solution of the continuous linear equations for *any* zonal
/// profile — the eastward invariant obeys `∂r/∂t + c·∂r/∂x = 0` and the
/// westward one is identically zero — so a run started from it carries one
/// wave, travelling east, and no Rossby energy at all.
pub fn kelvin_pulse_state(
    basin: Basin,
    params: PhysicalParams,
    deformation_radius_m: f64,
    wave_speed_m_per_s: f64,
    packet: Packet,
) -> OceanState {
    let mut state = OceanState::at_rest(basin.grid());
    let current_amplitude_m_per_s =
        packet.amplitude_m * wave_speed_m_per_s / params.mean_thermocline_depth_m();
    let profile = |x_m: f64| gaussian_envelope(x_m, packet.centre_x_m, packet.width_m);

    for j in 0..state.h().ny() {
        let waveguide = MeridionalStructure::Gravest
            .at(basin.y_of_row_m(H_STAGGERING, j), deformation_radius_m);
        for i in 0..state.h().nx() {
            let x_m = basin.x_of_column_m(H_STAGGERING, i);
            *state.h_mut().get_mut(i, j).expect("a cell centre") =
                packet.amplitude_m * profile(x_m) * waveguide;
        }
        for i in 0..state.u().nx() {
            let x_m = basin.x_of_column_m(U_STAGGERING, i);
            *state.u_mut().get_mut(i, j).expect("an east/west face") =
                current_amplitude_m_per_s * profile(x_m) * waveguide;
        }
    }
    state
}

/// The gravest-mode equatorial Rossby packet, Gaussian in `x` and travelling
/// west.
///
/// The long-wave mode of Matsuno 1966: `v ∝ ψ₁` with the zonal slope of the
/// envelope, `h/H ∝ (2ŷ² + 1)·e^{−ŷ²/2}` and `u/c ∝ (2ŷ² − 3)·e^{−ŷ²/2}`. In
/// the invariants those combine to `ψ₂` eastward and `−4·ψ₀` westward, which
/// is the decomposition every measurement of this mode reads it through. It is
/// the exact mode only as `k̂ → 0`, so a run started from it sheds an amplitude
/// `O(⟨k̂²⟩)` into the other branches — the stray-energy term of a budget.
pub fn gravest_rossby_packet_state(
    basin: Basin,
    params: PhysicalParams,
    deformation_radius_m: f64,
    wave_speed_m_per_s: f64,
    packet: Packet,
) -> OceanState {
    let mut state = OceanState::at_rest(basin.grid());
    let grid = basin.grid();
    let mean_depth_m = params.mean_thermocline_depth_m();
    let amplitude = packet.amplitude_m / mean_depth_m;
    let envelope = |x_m: f64| gaussian_envelope(x_m, packet.centre_x_m, packet.width_m);
    // `dE/dx` of that envelope, in m⁻¹.
    let envelope_slope_per_m =
        |x_m: f64| -(x_m - packet.centre_x_m) / (packet.width_m * packet.width_m) * envelope(x_m);

    let (h_nx, h_ny) = grid.field_shape(H_STAGGERING);
    for j in 0..h_ny {
        let y_hat = basin.y_of_row_m(H_STAGGERING, j) / deformation_radius_m;
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
        let y_hat = basin.y_of_row_m(U_STAGGERING, j) / deformation_radius_m;
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

    // `(8/3)·ŷ·e^{−ŷ²/2}` is `(4/3)·ψ₁`, and writing it that way is the point:
    // the mode is defined by its meridional velocity sitting on `ψ₁`, and every
    // other field of it follows from that.
    let (v_nx, v_ny) = grid.field_shape(V_STAGGERING);
    for j in 0..v_ny {
        let waveguide =
            MeridionalStructure::First.at(basin.y_of_row_m(V_STAGGERING, j), deformation_radius_m);
        for i in 0..v_nx {
            let x_m = basin.x_of_column_m(V_STAGGERING, i);
            *state
                .v_mut()
                .get_mut(i, j)
                .expect("the loop bounds are the field's own shape") = wave_speed_m_per_s
                * (4.0 / 3.0)
                * amplitude
                * deformation_radius_m
                * envelope_slope_per_m(x_m)
                * waveguide;
        }
    }

    state
}
