//! The closed basin's boundary condition: no flow through the coast.
//!
//! The basin of `CONTEXT.md` is "closed on all four boundaries", and on the
//! Arakawa C-grid of [ADR-0003] that statement is about four lines of points.
//! The zonal current anomaly `u` lives on east/west faces, indexed `0..=nx`,
//! so the western and eastern coasts *are* the `u` columns `i = 0` and
//! `i = nx`; the meridional anomaly `v` lives on north/south faces, so the
//! southern and northern coasts are the `v` rows `j = 0` and `j = ny`. The
//! normal velocity on those lines is not a degree of freedom of a closed
//! basin — it is the coast — and [`NoNormalFlow`] is what says so.
//!
//! Only the *normal* component is constrained. The `v` faces along a
//! north–south coast, and the `u` faces along an east–west one, are tangential
//! there; the linear v1 core has no lateral friction, so the closed basin is
//! free-slip and those velocities are left alone.
//!
//! # Each RK4 stage, not each completed step
//!
//! The deliverable of T-04.2 leaves open which of the two this is. It is each
//! stage, and the reason is the continuity equation. RK4 evaluates the whole
//! right-hand side at four stage states per step, each of them
//! `state + a·dt·k` for a previous stage's tendency `k`. If the wall
//! *acceleration* were left alone, every one of those intermediate states
//! would carry a wall velocity of order `dt·τ/(ρ₀·H)`, and `∂h/∂t` would pick
//! it up as a divergence: summed over the basin the discrete continuity
//! equation telescopes to the flux through the four walls, so a stage that
//! moves water through the coast changes the basin's total volume anomaly. A
//! condition applied only to the completed step would tidy the wall away
//! afterwards and leave that leak in `h`.
//!
//! So the condition is applied twice per step, and the two applications are
//! different jobs:
//!
//! - [`NoNormalFlow::apply_to_tendency`], at every stage, sets the wall
//!   acceleration to zero. Because RK4 combines states and tendencies
//!   linearly, a state whose walls are at rest and tendencies whose walls are
//!   at rest give stage states and a result whose walls are at rest — exactly,
//!   since `0.0 + a·dt·0.0` is `0.0` in IEEE-754, not merely nearly zero.
//! - [`NoNormalFlow::apply_to_state`], once at the start of a step, puts the
//!   incoming state on the boundary condition. That is what makes the
//!   guarantee unconditional: a state built by hand, or restored from a run,
//!   satisfies it from its first step rather than carrying a wall flow that
//!   the tendency rule would then faithfully preserve.
//!
//! Both are `O(nx + ny)` — the perimeter, not the area — and neither
//! allocates, so applying them per stage costs nothing measurable against the
//! right-hand side they follow.
//!
//! # Who applies it
//!
//! [`Solver`](crate::Solver) does, on both of its ways in. The condition is a
//! statement about the system being integrated rather than about any one term,
//! so it belongs where the terms are composed and not inside
//! [`ShallowWaterRhs`](crate::ShallowWaterRhs) or
//! [`CoriolisTerm`](crate::CoriolisTerm) — which is also what keeps those two
//! answerable, on their own, to the question of what the momentum equations do
//! with a stress at a wall. The consequence is that the guarantee is the
//! solver's and not the engine's: code that drives [`Rk4`](crate::Rk4) and a
//! right-hand side directly, as some of the Epic 01 and Epic 02 tests do, gets
//! no boundary condition unless it applies this one itself.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use termocline_grid::Field2D;

use crate::state::OceanState;

/// Value a held wall carries, in whatever unit the field it is written into is
/// read in: m/s in a state, m/s² in a tendency.
///
/// One constant for both, as in [`OceanState`]'s own `AT_REST`, because the
/// condition is the same statement read twice — the normal velocity at the
/// wall is zero, and therefore so is its rate of change.
const AT_REST: f64 = 0.0;

/// No normal flow through the four walls of a closed basin: `u = 0` on the
/// western and eastern boundary faces, `v = 0` on the southern and northern
/// ones.
///
/// Stateless — the walls of a basin are wherever its fields end, which every
/// [`OceanState`] already knows — so this is a namespace rather than a
/// configured object, and it can be applied to a state the solver did not
/// build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoNormalFlow;

impl NoNormalFlow {
    /// Bring `state` onto the boundary condition: set the normal velocity to
    /// exactly zero on all four walls.
    ///
    /// The tangential velocities and the thermocline depth anomaly are
    /// untouched — the basin is free-slip, and `h` is a cell-centered variable
    /// with no wall points at all.
    pub fn apply_to_state(state: &mut OceanState) {
        hold_walls_at_rest(state);
    }

    /// Hold the walls at rest through a right-hand-side evaluation: set the
    /// normal *acceleration* to exactly zero on all four walls of `tendency`.
    ///
    /// `tendency` is an [`OceanState`] read as rates, so this zeroes `∂u/∂t`
    /// on the zonal walls and `∂v/∂t` on the meridional ones — whatever the
    /// pressure gradient, the rotation, the damping and above all the surface
    /// stress wrote there.
    pub fn apply_to_tendency(tendency: &mut OceanState) {
        hold_walls_at_rest(tendency);
    }

    /// Bring a *diagnosed* horizontal flow onto the same condition: set
    /// `zonal_m_per_s` to zero on the western and eastern walls and
    /// `meridional_m_per_s` on the southern and northern ones.
    ///
    /// The two fields are at the `u` and `v` staggerings but are not the
    /// prognostic currents — the Epic 12 surface layer's wind-driven flow is
    /// the caller ([`crate::sst::SurfaceLayer`]). It gets the same condition
    /// for the same reason: a closed basin's coast is a coast to every flow in
    /// the model, not only to the one being integrated, and the surface
    /// layer's divergence is about to be read as an upwelling.
    ///
    /// It also removes the one place that divergence could have been read off
    /// an undefined number. Each mixed-layer component needs *both* stress
    /// components, so one of them arrives interpolated from the other's faces
    /// — and a C-grid interpolation is undefined on a wall face, which has
    /// cells on one side only. Holding the wall at rest says what the physics
    /// says there, instead of differencing whatever the interpolation left.
    pub fn apply_to_surface_flow(
        zonal_m_per_s: &mut Field2D<f64>,
        meridional_m_per_s: &mut Field2D<f64>,
    ) {
        zero_first_and_last_columns(zonal_m_per_s);
        zero_first_and_last_rows(meridional_m_per_s);
    }
}

/// Zero the normal component on the four walls of an [`OceanState`], read
/// either as a state or as a tendency.
///
/// A face field's wall lines are its own first and last lines along the axis
/// it is staggered on — a `u` field has one more column than the basin has
/// cells precisely so that both coasts are addressable (`termocline-grid`) —
/// so this is written in terms of the field's extents rather than in `nx` and
/// `ny` with an offset, per CODING_STANDARDS.md § Scope guards.
fn hold_walls_at_rest(state: &mut OceanState) {
    zero_first_and_last_columns(state.u_mut());
    zero_first_and_last_rows(state.v_mut());
}

/// Set the westernmost and easternmost columns of an east/west-face field to
/// rest. On a one-cell-wide basin they are the same column, written twice.
fn zero_first_and_last_columns(face: &mut Field2D<f64>) {
    let last = face.nx() - 1;
    for j in 0..face.ny() {
        for i in [0, last] {
            *face
                .get_mut(i, j)
                .expect("a field's own first and last columns are in bounds") = AT_REST;
        }
    }
}

/// Set the southernmost and northernmost rows of a north/south-face field to
/// rest. The meridional twin of [`zero_first_and_last_columns`].
fn zero_first_and_last_rows(face: &mut Field2D<f64>) {
    let last = face.ny() - 1;
    for j in [0, last] {
        for i in 0..face.nx() {
            *face
                .get_mut(i, j)
                .expect("a field's own first and last rows are in bounds") = AT_REST;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::NoNormalFlow;
    use crate::state::OceanState;
    use termocline_grid::Grid;

    /// A state whose every point carries `value`, in whatever unit its field
    /// is read in.
    fn filled(grid: Grid, value: f64) -> OceanState {
        let mut state = OceanState::at_rest(grid);
        state.u_mut().as_mut_slice().fill(value);
        state.v_mut().as_mut_slice().fill(value);
        state.h_mut().as_mut_slice().fill(value);
        state
    }

    #[test]
    fn only_the_normal_component_at_the_walls_is_touched() {
        // Free-slip: the tangential velocity along a coast, the interior, and
        // the cell-centered `h` all keep the value they arrived with.
        let grid = Grid::new(3, 2).expect("extents are non-zero");
        let mut state = filled(grid, 1.0);

        NoNormalFlow::apply_to_state(&mut state);

        for j in 0..state.u().ny() {
            assert_eq!(*state.u().get(0, j).expect("the western wall"), 0.0);
            assert_eq!(*state.u().get(3, j).expect("the eastern wall"), 0.0);
            for i in 1..3 {
                assert_eq!(*state.u().get(i, j).expect("an interior u face"), 1.0);
            }
        }
        for i in 0..state.v().nx() {
            assert_eq!(*state.v().get(i, 0).expect("the southern wall"), 0.0);
            assert_eq!(*state.v().get(i, 2).expect("the northern wall"), 0.0);
            assert_eq!(*state.v().get(i, 1).expect("an interior v face"), 1.0);
        }
        assert!(state.h().as_slice().iter().all(|h_m| *h_m == 1.0));
    }

    #[test]
    fn a_one_cell_basin_is_all_coast() {
        // Degenerate but reachable: with one cell on each axis the western and
        // eastern walls are the same field's two columns and there is no
        // interior at all, so the whole velocity field is held at rest.
        let grid = Grid::new(1, 1).expect("extents are non-zero");
        let mut tendency = filled(grid, 1.0);

        NoNormalFlow::apply_to_tendency(&mut tendency);

        assert!(tendency.u().as_slice().iter().all(|rate| *rate == 0.0));
        assert!(tendency.v().as_slice().iter().all(|rate| *rate == 0.0));
    }
}
