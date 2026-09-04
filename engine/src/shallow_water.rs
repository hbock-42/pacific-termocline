//! The right-hand side of the linear 1.5-layer shallow-water equations.
//!
//! The equations are the linear ones of `docs/planning/01-scientific-model.md`.
//! This module owns two of their terms — the pressure gradient `−g'·∇h` in the
//! momentum equations, the divergence `−H·(∂u/∂x + ∂v/∂y)` in the `h` equation
//! — plus the surface stress that makes the wind-stress argument of
//! [`shallow_water_rhs`] mean something:
//!
//! ```text
//! ∂u/∂t = −g'·∂h/∂x + τx/(ρ₀·H)
//! ∂v/∂t = −g'·∂h/∂y + τy/(ρ₀·H)
//! ∂h/∂t = −H·(∂u/∂x + ∂v/∂y)
//! ```
//!
//! The Coriolis terms (T-02.2) and the Rayleigh damping (T-02.4) fold into the
//! same evaluation later; nothing else belongs here. In particular the v1 core
//! is linear (CODING_STANDARDS.md § Scope guards): the divergence carries the
//! *mean* depth `H`, never the total `H + h`, and there is no advection.
//!
//! The spatial derivatives come from `termocline-numerics`, which is where the
//! C-grid neighbour arithmetic lives. Two consequences are worth stating,
//! because they are contracts rather than accidents:
//!
//! - A center→face difference is undefined on the basin's four boundary faces,
//!   which have a cell on one side only; the operators write zero there. So on
//!   a wall the acceleration this module produces is the wind stress alone,
//!   until Epic 04 gives the boundary a condition of its own.
//! - A face→center difference is defined at every cell, so `∂h/∂t` has no such
//!   gap.

use termocline_grid::{Field2D, Grid, H_STAGGERING};
use termocline_numerics::{CGridOperators, Spacing};

use crate::params::PhysicalParams;
use crate::state::OceanState;
use crate::wind::WindStress;

/// The right-hand side of the shallow-water equations for one basin, one cell
/// spacing and one parameter set, holding the scratch space it needs.
///
/// One evaluator is built per run and reused for every stage of every step, so
/// a whole simulation allocates its intermediate fields exactly once
/// (CODING_STANDARDS.md § Performance). [`shallow_water_rhs`] is the
/// allocating convenience wrapper for tests and one-off evaluation.
#[derive(Debug, Clone)]
pub struct ShallowWaterRhs {
    /// Physical parameters `(g', H, r, β, ρ₀)` of the ocean being integrated.
    params: PhysicalParams,
    /// The C-grid derivative operators at this basin's cell spacing.
    operators: CGridOperators,
    /// `∂u/∂x` at cell centers, in s⁻¹. Scratch, rewritten every evaluation.
    zonal_divergence_per_s: Field2D<f64>,
    /// `∂v/∂y` at cell centers, in s⁻¹. Scratch, rewritten every evaluation.
    meridional_divergence_per_s: Field2D<f64>,
}

impl ShallowWaterRhs {
    /// An evaluator for `grid` at `spacing`, for an ocean with `params`.
    #[must_use]
    pub fn new(grid: Grid, spacing: Spacing, params: PhysicalParams) -> Self {
        Self {
            params,
            operators: CGridOperators::new(grid, spacing),
            zonal_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            meridional_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
        }
    }

    /// Write the time derivative of `state` under `wind_stress` into
    /// `tendency`.
    ///
    /// `tendency` is an [`OceanState`] read as rates: `∂h/∂t` in m/s, `∂u/∂t`
    /// and `∂v/∂t` in m/s². Every point of it is written, so the same buffer
    /// can be reused across RK4 stages without carrying a stage's values into
    /// the next.
    ///
    /// # Panics
    /// If `state`, `wind_stress` or `tendency` covers a different basin from
    /// the one this evaluator was built for. A shape mismatch means the
    /// calling code is wrong, which is what panics are for
    /// (CODING_STANDARDS.md § Correctness and failure).
    pub fn evaluate(
        &mut self,
        state: &OceanState,
        wind_stress: &WindStress,
        tendency: &mut OceanState,
    ) {
        self.check_grid("state", state.grid());
        self.check_grid("wind stress", wind_stress.grid());
        self.check_grid("tendency", tendency.grid());

        // Momentum: the pressure gradient lands straight on the velocity
        // faces, so it is written into the tendency and then turned into an
        // acceleration in place.
        self.operators
            .ddx_center_to_face(state.h(), tendency.u_mut());
        self.operators
            .ddy_center_to_face(state.h(), tendency.v_mut());
        let minus_g_prime_m_per_s2 = -self.params.reduced_gravity_m_per_s2();
        // Mass of the upper layer per unit area, in kg/m²: the `ρ₀·H` a
        // surface stress in Pa is divided by to become an acceleration.
        let layer_mass_kg_per_m2 =
            self.params.reference_density_kg_per_m3() * self.params.mean_thermocline_depth_m();
        turn_gradient_into_acceleration(
            tendency.u_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_x_pa(),
            layer_mass_kg_per_m2,
        );
        turn_gradient_into_acceleration(
            tendency.v_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_y_pa(),
            layer_mass_kg_per_m2,
        );

        // Continuity: both halves of the divergence land on the cell centers,
        // where `h` lives.
        self.operators
            .ddx_face_to_center(state.u(), &mut self.zonal_divergence_per_s);
        self.operators
            .ddy_face_to_center(state.v(), &mut self.meridional_divergence_per_s);
        let minus_mean_depth_m = -self.params.mean_thermocline_depth_m();
        let thickness_rate = tendency.h_mut().as_mut_slice().iter_mut();
        let divergence = self
            .zonal_divergence_per_s
            .as_slice()
            .iter()
            .zip(self.meridional_divergence_per_s.as_slice());
        for (rate, (zonal, meridional)) in thickness_rate.zip(divergence) {
            *rate = minus_mean_depth_m * (zonal + meridional);
        }
    }

    /// Panic unless `grid` is the basin this evaluator was built for.
    fn check_grid(&self, role: &str, grid: Grid) {
        assert!(
            grid == self.operators.grid(),
            "{role} covers {grid:?}, but this right-hand side was built for {:?}",
            self.operators.grid()
        );
    }
}

/// Overwrite a gradient of `h`, in m/m, with the acceleration it produces:
/// `−g'·∂h/∂· + τ/(ρ₀·H)`, in m/s².
///
/// In-place because the gradient operators write straight into the tendency
/// buffer, which is where the acceleration has to end up; `field` therefore
/// arrives holding a gradient and leaves holding an acceleration.
fn turn_gradient_into_acceleration(
    field: &mut Field2D<f64>,
    minus_g_prime_m_per_s2: f64,
    stress_pa: &Field2D<f64>,
    layer_mass_kg_per_m2: f64,
) {
    for (value, stress) in field.as_mut_slice().iter_mut().zip(stress_pa.as_slice()) {
        *value = minus_g_prime_m_per_s2 * *value + stress / layer_mass_kg_per_m2;
    }
}

/// The time derivative of `state` under `wind_stress`, as a freshly allocated
/// tendency.
///
/// The named deliverable of T-02.3, and the convenient form for a test or a
/// one-off evaluation. It allocates a state and an evaluator per call, so the
/// time loop uses [`ShallowWaterRhs::evaluate`] instead
/// (CODING_STANDARDS.md § Performance).
///
/// # Panics
/// If `wind_stress` covers a different basin from `state`.
#[must_use]
pub fn shallow_water_rhs(
    state: &OceanState,
    params: PhysicalParams,
    spacing: Spacing,
    wind_stress: &WindStress,
) -> OceanState {
    let mut evaluator = ShallowWaterRhs::new(state.grid(), spacing, params);
    let mut tendency = OceanState::at_rest(state.grid());
    evaluator.evaluate(state, wind_stress, &mut tendency);
    tendency
}
