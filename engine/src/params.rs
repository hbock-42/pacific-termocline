//! The fixed physical parameters of the 1.5-layer reduced-gravity model.
//!
//! Everything here is SI and stays SI: `g'` in m/s², `H` in m, `r` in s⁻¹, `β`
//! in m⁻¹s⁻¹, `ρ₀` in kg/m³. Epic 07 validates the solver against analytic
//! formulas — `c = √(g'·H)`, `Le = √(c/β)` — which hold only in a consistent
//! unit system, so a parameter never carries a unit its name does not state
//! and is never rescaled on the way in.
//!
//! These are the *constants* of a scenario, not fields: they do not vary with
//! position or time. The wind stress, which does, is Epic 03's.

use std::fmt;

/// Meridional gradient of the Coriolis parameter at the equator, in m⁻¹s⁻¹.
///
/// `β = 2Ω·cos(φ)/R` evaluated at `φ = 0`; the value quoted for the equatorial
/// beta-plane in `CONTEXT.md` and in `docs/planning/01-scientific-model.md`.
pub const EQUATORIAL_BETA_PER_M_PER_S: f64 = 2.3e-11;

/// Reference seawater density `ρ₀`, in kg/m³.
///
/// The standard Boussinesq reference density for the upper tropical ocean
/// (Gill, *Atmosphere–Ocean Dynamics*, appendix 3). It enters the momentum
/// equations only through the wind-stress term `τ/(ρ₀·H)`.
pub const SEAWATER_REFERENCE_DENSITY_KG_PER_M3: f64 = 1025.0;

/// Why a set of physical parameters was rejected.
///
/// These describe invalid *input* — a scenario asking for an unphysical
/// ocean — so they are returned rather than panicked, and each names the
/// offending parameter and the value it carried.
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicalParamsError {
    /// A parameter that must be strictly positive and finite was not.
    NotPositive {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
    /// A parameter that must be non-negative and finite was negative (or not a
    /// number).
    Negative {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
}

impl fmt::Display for PhysicalParamsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPositive { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and greater than 0"
            ),
            Self::Negative { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and at least 0"
            ),
        }
    }
}

impl std::error::Error for PhysicalParamsError {}

/// The fixed physical parameters of one scenario's ocean.
///
/// Constructed once per run and read from every right-hand-side evaluation, so
/// it is `Copy` and validated at the boundary rather than at each use.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalParams {
    /// Reduced gravity `g'`, in m/s².
    reduced_gravity_m_per_s2: f64,
    /// Mean thermocline depth `H`, in metres. The resting upper-layer
    /// thickness — a total depth, unlike the anomaly `h`.
    mean_depth_m: f64,
    /// Rayleigh damping coefficient `r`, in s⁻¹.
    rayleigh_damping_per_s: f64,
    /// Meridional gradient of the Coriolis parameter `β`, in m⁻¹s⁻¹.
    beta_per_m_per_s: f64,
    /// Reference seawater density `ρ₀`, in kg/m³.
    reference_density_kg_per_m3: f64,
}

impl PhysicalParams {
    /// The parameter set `(g', H, r, β, ρ₀)`, in SI units.
    ///
    /// `r` may be zero — the undamped limit is one of the validation targets
    /// in `docs/planning/01-scientific-model.md` — but never negative, which
    /// would amplify rather than damp. The other four must be strictly
    /// positive: a zero `g'` or `H` collapses the wave speed, a zero `β`
    /// removes the equatorial waveguide, and a zero `ρ₀` divides by zero in
    /// the wind-stress term.
    ///
    /// # Errors
    /// [`PhysicalParamsError::NotPositive`] or
    /// [`PhysicalParamsError::Negative`], naming the first parameter that
    /// failed and the value it carried.
    pub fn new(
        reduced_gravity_m_per_s2: f64,
        mean_depth_m: f64,
        rayleigh_damping_per_s: f64,
        beta_per_m_per_s: f64,
        reference_density_kg_per_m3: f64,
    ) -> Result<Self, PhysicalParamsError> {
        check_positive("reduced_gravity_m_per_s2", reduced_gravity_m_per_s2)?;
        check_positive("mean_depth_m", mean_depth_m)?;
        check_non_negative("rayleigh_damping_per_s", rayleigh_damping_per_s)?;
        check_positive("beta_per_m_per_s", beta_per_m_per_s)?;
        check_positive("reference_density_kg_per_m3", reference_density_kg_per_m3)?;
        Ok(Self {
            reduced_gravity_m_per_s2,
            mean_depth_m,
            rayleigh_damping_per_s,
            beta_per_m_per_s,
            reference_density_kg_per_m3,
        })
    }

    /// Reduced gravity `g'`, in m/s².
    #[must_use]
    pub const fn reduced_gravity_m_per_s2(self) -> f64 {
        self.reduced_gravity_m_per_s2
    }

    /// Mean thermocline depth `H`, in metres.
    #[must_use]
    pub const fn mean_depth_m(self) -> f64 {
        self.mean_depth_m
    }

    /// Rayleigh damping coefficient `r`, in s⁻¹.
    #[must_use]
    pub const fn rayleigh_damping_per_s(self) -> f64 {
        self.rayleigh_damping_per_s
    }

    /// Meridional gradient of the Coriolis parameter `β`, in m⁻¹s⁻¹.
    #[must_use]
    pub const fn beta_per_m_per_s(self) -> f64 {
        self.beta_per_m_per_s
    }

    /// Reference seawater density `ρ₀`, in kg/m³.
    #[must_use]
    pub const fn reference_density_kg_per_m3(self) -> f64 {
        self.reference_density_kg_per_m3
    }

    /// Kelvin wave speed `c = √(g'·H)`, in m/s.
    ///
    /// The fastest signal in the model, and therefore the speed that bounds
    /// the stable timestep (`CONTEXT.md`). Both factors are checked positive
    /// at construction, so the square root is real.
    #[must_use]
    pub fn kelvin_wave_speed_m_per_s(self) -> f64 {
        (self.reduced_gravity_m_per_s2 * self.mean_depth_m).sqrt()
    }

    /// Equatorial deformation radius `Le = √(c/β)`, in metres.
    ///
    /// The meridional scale over which equatorial waves decay away from the
    /// equator (`CONTEXT.md`).
    #[must_use]
    pub fn equatorial_deformation_radius_m(self) -> f64 {
        (self.kelvin_wave_speed_m_per_s() / self.beta_per_m_per_s).sqrt()
    }
}

fn check_positive(parameter: &'static str, value: f64) -> Result<(), PhysicalParamsError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(PhysicalParamsError::NotPositive { parameter, value });
    }
    Ok(())
}

fn check_non_negative(parameter: &'static str, value: f64) -> Result<(), PhysicalParamsError> {
    if !value.is_finite() || value < 0.0 {
        return Err(PhysicalParamsError::Negative { parameter, value });
    }
    Ok(())
}
