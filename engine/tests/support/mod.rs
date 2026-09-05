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

use engine::{Basin, OceanState, PhysicalParams, H_STAGGERING};

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
///
/// It is the value the assertions are made against, so it must not come from
/// the engine: an expected value read out of the thing being measured agrees
/// with it by construction (CODING_STANDARDS.md § *Tests*). The engine computes
/// the same number from the same two parameters, which is what makes the
/// comparison a test rather than a tautology.
pub fn kelvin_wave_speed_m_per_s() -> f64 {
    (PACIFIC_REDUCED_GRAVITY_M_PER_S2 * PACIFIC_MEAN_DEPTH_M).sqrt()
}

/// Equatorial deformation radius `Le = √(c/β)`, in metres — the meridional
/// scale of the waveguide (`CONTEXT.md`), written out from the same definition
/// and for the same reason.
pub fn equatorial_deformation_radius_m() -> f64 {
    (kelvin_wave_speed_m_per_s() / BETA_PER_M_PER_S).sqrt()
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
            row_y_m: (0..basin.grid().ny())
                .map(|j| basin.y_of_row_m(H_STAGGERING, j))
                .collect(),
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
