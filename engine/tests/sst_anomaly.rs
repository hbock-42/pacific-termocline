//! T-12.1 — the mixed-layer SST anomaly `T'`, its equation, and the proof that
//! turning it on changes nothing about the validated linear core.
//!
//! The equation under test is the linearized mixed-layer temperature equation
//! of the intermediate coupled models (Zebiak & Cane, *Mon. Wea. Rev.* 115,
//! 1987, § 2b; Battisti, *J. Atmos. Sci.* 45, 1988), written in the anomaly
//! variables `CONTEXT.md` already names:
//!
//! ```text
//! ∂T'/∂t = −u'·∂T̄/∂x + (w⁺/H_m)·(γ·h − T') − ε_T·T'
//! ```
//!
//! Three strands are checked here, and each has an independent source:
//!
//! - **The wind-implied upwelling.** The mixed layer's steady Rayleigh-drag
//!   momentum balance has a closed-form solution, and its `r_s → 0` limit is
//!   Ekman's (1905) transport `−τx/(ρ₀·f)`. Both are written out from theory
//!   below and never asked of the engine.
//! - **The equation's assembly.** Each term is isolated by a configuration in
//!   which the other two vanish identically, so the expected tendency is a
//!   product of numbers stated in the test.
//! - **Additivity.** The acceptance criterion of the ticket: with the coupling
//!   disabled the validated core is untouched. Asserted here in its strongest
//!   form — the same forced basin, stepped with and without the coupling,
//!   agrees on `h`, `u` and `v` bit for bit.

use engine::sst::{
    mixed_layer_velocity_m_per_s, SstParams, SstParamsError, SstTerm, SurfaceLayer,
    DEFAULT_SURFACE_DRAG_PER_S,
};
use engine::{
    Basin, BetaPlane, Grid, OceanState, PhysicalParams, Scenario, ScenarioConfig, ScenarioError,
    Solver, Spacing, StateVector, WindStressField, H_STAGGERING,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s² (Gill, *Atmosphere–Ocean Dynamics*, ch. 11).
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H`, in metres.
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: a 100-day decay.
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0 / (100.0 * SECONDS_PER_DAY);
/// Meridional gradient of the Coriolis parameter, in m⁻¹s⁻¹ (`CONTEXT.md`).
const BETA_PER_M_PER_S: f64 = 2.3e-11;
/// Seconds in a day.
const SECONDS_PER_DAY: f64 = 86_400.0;

/// Mixed-layer depth `H_m`, in metres. The wind-mixed surface layer of the
/// equatorial Pacific is 50 m deep (Zebiak & Cane 1987, § 2b).
const MIXED_LAYER_DEPTH_M: f64 = 50.0;
/// Zonal gradient of the mean SST, in K/m: the cold tongue is about 6 K colder
/// than the warm pool over the basin's 15 000 km, and the gradient is negative
/// because the ocean cools eastward.
const MEAN_ZONAL_SST_GRADIENT_K_PER_M: f64 = -4.0e-7;
/// Sensitivity `γ = ∂T_sub/∂h` of the entrained water's temperature to the
/// thermocline depth anomaly, in K/m (Zebiak & Cane 1987, § 2c).
const SUBSURFACE_SENSITIVITY_K_PER_M: f64 = 0.1;
/// Thermal damping `ε_T` of an SST anomaly, in s⁻¹: a 125-day relaxation to
/// the climatological surface heat flux (Zebiak & Cane 1987, § 2b).
const THERMAL_DAMPING_PER_S: f64 = 1.0 / (125.0 * SECONDS_PER_DAY);

/// Equatorial trade-wind stress, in Pa. Easterly, so negative
/// (`CONTEXT.md`, *Wind stress*).
const TRADE_WIND_STRESS_PA: f64 = -0.05;

fn physical_params() -> PhysicalParams {
    PhysicalParams::new(
        REDUCED_GRAVITY_M_PER_S2,
        MEAN_THERMOCLINE_DEPTH_M,
        RAYLEIGH_DAMPING_PER_S,
        BETA_PER_M_PER_S,
        SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("these are the standard equatorial-Pacific parameters")
}

fn sst_params() -> SstParams {
    SstParams::new(
        MIXED_LAYER_DEPTH_M,
        DEFAULT_SURFACE_DRAG_PER_S,
        MEAN_ZONAL_SST_GRADIENT_K_PER_M,
        SUBSURFACE_SENSITIVITY_K_PER_M,
        THERMAL_DAMPING_PER_S,
    )
    .expect("these are the standard Zebiak-Cane mixed-layer parameters")
}

/// A basin `ny` rows tall spanning `meridional_extent_m`, centred on the
/// equator, with `ny` odd so that one row of cell centers lies exactly on it.
fn equatorial_basin(nx: usize, ny: usize, dx_m: f64, meridional_extent_m: f64) -> Basin {
    assert!(
        ny % 2 == 1,
        "an odd row count puts a center row on the equator"
    );
    let grid = Grid::new(nx, ny).expect("a basin has cells");
    let spacing = Spacing::new(dx_m, meridional_extent_m / ny as f64).expect("cells have width");
    Basin::centered_on_equator(grid, spacing)
}

/// The index of the cell-center row that lies on the equator.
const fn equator_row(ny: usize) -> usize {
    (ny - 1) / 2
}

// ---------------------------------------------------------------------------
// The wind-implied upwelling
// ---------------------------------------------------------------------------

/// Meridional Ekman transport, in m²/s, from Ekman (1905): the wind-driven
/// transport of the surface layer is `−τx/(ρ₀·f)`, ninety degrees to the right
/// of the stress in the northern hemisphere. Written out from theory.
fn ekman_transport_m2_per_s(tau_x_pa: f64, coriolis_per_s: f64) -> f64 {
    -tau_x_pa / (SEAWATER_REFERENCE_DENSITY_KG_PER_M3 * coriolis_per_s)
}

#[test]
fn the_mixed_layer_transport_converges_on_the_ekman_transport_as_the_drag_weakens() {
    // The surface layer's steady balance `r_s·u − f·v = τx/(ρ₀·H_m)`,
    // `r_s·v + f·u = τy/(ρ₀·H_m)` reduces to Ekman's transport in the limit
    // `r_s/|f| → 0`, and its relative departure from it is exactly
    // `(r_s/f)²/(1 + (r_s/f)²)` — second order in the small parameter. So the
    // error must fall by four when `r_s/|f|` is halved, which is the
    // convergence CODING_STANDARDS.md asks for in place of a point check.
    let coriolis_per_s = BETA_PER_M_PER_S * 1.0e6; // f at 1000 km from the equator.
    let layer_mass_kg_per_m2 = SEAWATER_REFERENCE_DENSITY_KG_PER_M3 * MIXED_LAYER_DEPTH_M;
    let expected = ekman_transport_m2_per_s(TRADE_WIND_STRESS_PA, coriolis_per_s);

    let mut samples = Vec::new();
    for halvings in 0..3 {
        let surface_drag_per_s = DEFAULT_SURFACE_DRAG_PER_S / f64::from(1 << halvings);
        let (_, v_ml_m_per_s) = mixed_layer_velocity_m_per_s(
            TRADE_WIND_STRESS_PA,
            0.0,
            coriolis_per_s,
            layer_mass_kg_per_m2,
            surface_drag_per_s,
        );
        let transport_m2_per_s = v_ml_m_per_s * MIXED_LAYER_DEPTH_M;
        let smallness = (surface_drag_per_s / coriolis_per_s).powi(2);
        samples.push((
            smallness,
            ((transport_m2_per_s - expected) / expected).abs(),
        ));
    }

    for (smallness, error) in samples.iter().copied() {
        // The departure is not merely second order, it is known exactly: the
        // solution's transport is `1/(1 + x)` of Ekman's with
        // `x = (r_s/f)²`, so the relative error is `x/(1 + x)`. Asserting the
        // closed form rather than a slack band around a ratio of four leaves
        // nothing to a tolerance — the only inexactness left is the
        // floating-point arithmetic that produced both sides.
        let exact = smallness / (1.0 + smallness);
        // The budget is a few ulps of *one*, not of `exact`. Both sides are
        // formed by subtracting two transports that agree to better than a
        // part in fifty, so the cancellation carries the rounding of the
        // order-one operands into a difference two orders smaller; scaling the
        // budget by the small difference would be asking the arithmetic for
        // precision it never had.
        assert!(
            (error - exact).abs() <= 16.0 * f64::EPSILON,
            "at (r_s/f)² = {smallness:e} the transport should fall short of Ekman's by \
             exactly {exact:e} of itself, but it falls short by {error:e}"
        );
    }

    // ...and that closed form is second order in `r_s/f`, which is what the
    // halvings are there to make legible: each one very nearly quarters the
    // departure. "Very nearly" is itself exact — `x/(1 + x)` with `x` quartered
    // gives a ratio of `4·(1 + x)/(1 + 4x)`, approaching four from below — so
    // this is asserted against its own closed form too, and not against a band
    // around four.
    for pair in samples.windows(2) {
        let ((coarse_smallness, coarse), (fine_smallness, fine)) = (pair[0], pair[1]);
        assert!(
            coarse_smallness > fine_smallness,
            "the samples must run from the strongest drag to the weakest"
        );
        let observed = coarse / fine;
        let exact = 4.0 * (1.0 + fine_smallness) / (1.0 + 4.0 * fine_smallness);
        assert!(
            exact < 4.0,
            "the ratio approaches four from below, so a prediction at or above it is a \
             mis-derivation rather than a near miss"
        );
        // A part in a billion. The two errors being divided each carry the
        // cancellation noise bounded above — a few ulps of one against a
        // difference of order `x` — so the ratio inherits about `ε/x ≈ 10⁻¹⁴`
        // of relative noise; a part in a billion is far inside any physical
        // slack while leaving that arithmetic alone.
        assert!(
            (observed - exact).abs() <= 1.0e-9 * exact,
            "halving the surface drag should shrink the departure from the Ekman transport \
             by exactly {exact}, but it went from {coarse:e} to {fine:e}, a ratio of {observed}"
        );
    }
}

/// Equatorial upwelling implied by a uniform easterly stress, in m/s, from the
/// closed-form solution of the surface layer's balance:
///
/// ```text
/// w(y) = −(β·τx/ρ₀) · (r_s² − β²y²) / (r_s² + β²y²)²
/// ```
///
/// which at the equator is `−β·τx/(ρ₀·r_s²)`. Positive — upward — for the
/// easterly `τx < 0` of the alizés, which is the equatorial upwelling the cold
/// tongue is made of. Written out from theory, not from the engine.
fn analytic_equatorial_upwelling_m_per_s(tau_x_pa: f64, surface_drag_per_s: f64) -> f64 {
    -BETA_PER_M_PER_S * tau_x_pa
        / (SEAWATER_REFERENCE_DENSITY_KG_PER_M3 * surface_drag_per_s * surface_drag_per_s)
}

/// The relative error the C-grid difference makes on that upwelling at a cell
/// height of `dy_m`, to leading order — the tolerance every fixed-resolution
/// check of the equatorial upwelling is entitled to, derived rather than
/// observed.
///
/// The discrete `∂v_ml/∂y` at a center row sitting on the equator is the
/// centred difference `(v(dy/2) − v(−dy/2))/dy`, whose Taylor series is
/// `v'(0) + (dy²/24)·v'''(0) + O(dy⁴)`. Writing `v_ml = (a/r_s)·f(u)` with
/// `u = βy/r_s`, `a = −τx/(ρ₀·H_m)` and `f(u) = u/(1 + u²) = u − u³ + u⁵ − …`,
/// the two derivatives are `f'(0) = 1` and `f'''(0) = −6`, so
///
/// ```text
/// (dy²/24)·v'''(0) / v'(0) = −(1/4)·(β·dy/r_s)²
/// ```
///
/// Second order in `dy`, as [ADR-0003]'s C-grid difference promises, and set
/// by the one meridional scale the upwelling has: the half-width `r_s/β` over
/// which the Ekman singularity is smoothed.
///
/// [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md
fn upwelling_truncation_fraction(dy_m: f64, surface_drag_per_s: f64) -> f64 {
    let cells_per_equatorial_scale = BETA_PER_M_PER_S * dy_m / surface_drag_per_s;
    0.25 * cells_per_equatorial_scale * cells_per_equatorial_scale
}

#[test]
fn the_equatorial_upwelling_converges_on_the_analytic_ekman_divergence_at_second_order() {
    // The upwelling is the divergence of the mixed-layer flow, and the C-grid
    // difference that computes it is second-order accurate (ADR-0003). So
    // refining `dy` must shrink the departure from the closed form above at
    // second order — the order the scheme claims, asserted across three
    // resolutions rather than against one fixed threshold.
    let params = physical_params();
    let sst = sst_params();
    let expected =
        analytic_equatorial_upwelling_m_per_s(TRADE_WIND_STRESS_PA, sst.surface_drag_per_s());
    let meridional_extent_m = 4.0e6;

    let mut samples = Vec::new();
    for ny in [81_usize, 161, 321] {
        let basin = equatorial_basin(4, ny, 1.0e5, meridional_extent_m);
        let plane = BetaPlane::of_basin(params, basin);
        let mut layer = SurfaceLayer::new(basin.grid(), basin.spacing(), plane, params, sst);
        let wind = WindStressField::uniform(basin.grid(), TRADE_WIND_STRESS_PA, 0.0);
        layer.diagnose(&wind);
        let measured = *layer
            .upwelling_m_per_s()
            .get(2, equator_row(ny))
            .expect("the equator row is inside the basin");
        samples.push((meridional_extent_m / ny as f64, (measured - expected).abs()));
    }

    for pair in samples.windows(2) {
        let ((coarse_dy, coarse_error), (fine_dy, fine_error)) = (pair[0], pair[1]);
        let order = (coarse_error / fine_error).ln() / (coarse_dy / fine_dy).ln();
        // Second order is what the centred C-grid difference gives. The three
        // resolutions are chosen fine enough to be in the asymptotic regime:
        // the upwelling's meridional structure turns over at `|βy| = r_s`,
        // about 250 km from the equator, so cells of 50 km and below resolve
        // it and the ±0.1 band is the room the remaining higher-order terms
        // need.
        assert!(
            (order - 2.0).abs() < 0.1,
            "the upwelling should converge at second order, but the error went from \
             {coarse_error:e} at dy = {coarse_dy} m to {fine_error:e} at dy = {fine_dy} m, \
             an observed order of {order}"
        );
    }
}

#[test]
fn the_upwelling_is_upward_at_the_equator_under_the_alizes_and_downward_under_a_westerly() {
    // Easterly stress drives poleward Ekman flow on both flanks of the
    // equator, and mass conservation fills the divergence from below: this is
    // the sign that makes the cold tongue cold. A westerly reverses it.
    let params = physical_params();
    let sst = sst_params();
    let ny = 41;
    let basin = equatorial_basin(4, ny, 1.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let mut layer = SurfaceLayer::new(basin.grid(), basin.spacing(), plane, params, sst);

    layer.diagnose(&WindStressField::uniform(
        basin.grid(),
        TRADE_WIND_STRESS_PA,
        0.0,
    ));
    let easterly = *layer.upwelling_m_per_s().get(2, equator_row(ny)).unwrap();
    assert!(
        easterly > 0.0,
        "easterly trades must upwell at the equator, but w = {easterly} m/s"
    );

    layer.diagnose(&WindStressField::uniform(
        basin.grid(),
        -TRADE_WIND_STRESS_PA,
        0.0,
    ));
    let westerly = *layer.upwelling_m_per_s().get(2, equator_row(ny)).unwrap();
    assert!(
        westerly < 0.0,
        "a westerly must downwell at the equator, but w = {westerly} m/s"
    );
}

#[test]
fn the_wind_driven_flow_does_not_cross_the_coast() {
    // A closed basin's coast is a coast to the surface layer too, and the
    // upwelling is that layer's divergence — so a wall carrying flow would put
    // a fictitious `w` in the perimeter cells that entrainment would then feed
    // on. It matters most under a wind with a meridional component, which is
    // where each mixed-layer component needs the *other* stress interpolated
    // onto its faces, and where a C-grid interpolation is undefined on a wall.
    // Exactly zero, not nearly: the condition is an assignment.
    let params = physical_params();
    let sst = sst_params();
    let ny = 21;
    let basin = equatorial_basin(6, ny, 1.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let mut layer = SurfaceLayer::new(basin.grid(), basin.spacing(), plane, params, sst);

    // Both components non-zero, which no shipped forcing produces today — the
    // point is that the boundary condition does not depend on that staying
    // true.
    layer.diagnose(&WindStressField::uniform(
        basin.grid(),
        TRADE_WIND_STRESS_PA,
        0.02,
    ));

    let zonal = layer.zonal_flow_m_per_s();
    let last_column = zonal.nx() - 1;
    for j in 0..zonal.ny() {
        for i in [0, last_column] {
            assert_eq!(
                *zonal.get(i, j).expect("a wall column is in bounds"),
                0.0,
                "the mixed layer flows through the meridional wall at column {i}, row {j}"
            );
        }
    }
    let meridional = layer.meridional_flow_m_per_s();
    let last_row = meridional.ny() - 1;
    for i in 0..meridional.nx() {
        for j in [0, last_row] {
            assert_eq!(
                *meridional.get(i, j).expect("a wall row is in bounds"),
                0.0,
                "the mixed layer flows through the zonal wall at column {i}, row {j}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The equation, term by term
// ---------------------------------------------------------------------------

/// A coupled state over `grid` with every field set to a constant.
fn uniform_coupled_state(grid: Grid, h_m: f64, u_m_per_s: f64, sst_anomaly_k: f64) -> OceanState {
    let mut state = OceanState::at_rest_with_sst_anomaly(grid);
    state.h_mut().as_mut_slice().fill(h_m);
    state.u_mut().as_mut_slice().fill(u_m_per_s);
    state
        .sst_anomaly_k_mut()
        .expect("this state carries an SST anomaly")
        .as_mut_slice()
        .fill(sst_anomaly_k);
    state
}

#[test]
fn zonal_advection_of_the_mean_gradient_is_the_only_term_left_under_a_calm_ocean() {
    // With no wind there is no upwelling and so no entrainment, and with
    // `T' = 0` there is no thermal damping. What survives is `−u'·∂T̄/∂x`, and
    // a uniform `u'` interpolates onto the cell centers exactly, so the
    // expected tendency is a product of two numbers stated above — machine
    // precision, not a physical tolerance.
    let params = physical_params();
    let sst = sst_params();
    let basin = equatorial_basin(6, 21, 1.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let mut term = SstTerm::new(basin.grid(), basin.spacing(), plane, params, sst);

    let zonal_current_m_per_s = 0.2;
    let state = uniform_coupled_state(basin.grid(), 0.0, zonal_current_m_per_s, 0.0);
    let mut tendency = OceanState::at_rest_with_sst_anomaly(basin.grid());
    term.add_to_tendency(&state, &WindStressField::calm(basin.grid()), &mut tendency);

    let expected_k_per_s = -zonal_current_m_per_s * MEAN_ZONAL_SST_GRADIENT_K_PER_M;
    let measured = *tendency
        .sst_anomaly_k()
        .expect("the tendency carries an SST anomaly")
        .get(3, 10)
        .expect("an interior cell");
    assert!(
        (measured - expected_k_per_s).abs() <= 8.0 * f64::EPSILON * expected_k_per_s.abs(),
        "a uniform current over a uniform mean gradient should warm at {expected_k_per_s} K/s, \
         but the tendency is {measured} K/s"
    );
}

#[test]
fn thermal_damping_is_the_only_term_left_for_a_still_anomaly_under_a_westerly() {
    // A westerly downwells at the equator, and only upwelling entrains, so the
    // entrainment term is identically zero there; with `u' = 0` the advection
    // term is too. What is left is `−ε_T·T'`, again a product of two stated
    // numbers.
    let params = physical_params();
    let sst = sst_params();
    let ny = 21;
    let basin = equatorial_basin(6, ny, 1.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let mut term = SstTerm::new(basin.grid(), basin.spacing(), plane, params, sst);

    let anomaly_k = 1.5;
    // A deep thermocline anomaly is carried too, to show it cannot leak in
    // through an entrainment term that is switched off.
    let state = uniform_coupled_state(basin.grid(), 20.0, 0.0, anomaly_k);
    let mut tendency = OceanState::at_rest_with_sst_anomaly(basin.grid());
    let westerly = WindStressField::uniform(basin.grid(), -TRADE_WIND_STRESS_PA, 0.0);
    term.add_to_tendency(&state, &westerly, &mut tendency);

    let expected_k_per_s = -THERMAL_DAMPING_PER_S * anomaly_k;
    let measured = *tendency
        .sst_anomaly_k()
        .unwrap()
        .get(3, equator_row(ny))
        .expect("the equator row is inside the basin");
    assert!(
        (measured - expected_k_per_s).abs() <= 8.0 * f64::EPSILON * expected_k_per_s.abs(),
        "a still anomaly under a downwelling wind should decay at {expected_k_per_s} K/s, \
         but the tendency is {measured} K/s"
    );
}

#[test]
fn entrainment_carries_the_thermocline_anomaly_into_the_mixed_layer() {
    // The coupling the ticket is about: a deeper thermocline (`h > 0`) puts
    // warmer water under the mixed layer, and the upwelling the alizés imply
    // pumps it in at `w/H_m · γ·h`. With `u' = 0` and `T' = 0` that is the
    // whole tendency, so every factor of the expected value is written out
    // from theory — the closed-form equatorial upwelling above, and the two
    // mixed-layer constants stated at the top of this file. Nothing is read
    // back from the engine.
    let params = physical_params();
    let sst = sst_params();
    let ny = 161;
    let basin = equatorial_basin(6, ny, 1.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let mut term = SstTerm::new(basin.grid(), basin.spacing(), plane, params, sst);

    let thermocline_anomaly_m = 10.0;
    let state = uniform_coupled_state(basin.grid(), thermocline_anomaly_m, 0.0, 0.0);
    let mut tendency = OceanState::at_rest_with_sst_anomaly(basin.grid());
    let alizes = WindStressField::uniform(basin.grid(), TRADE_WIND_STRESS_PA, 0.0);
    term.add_to_tendency(&state, &alizes, &mut tendency);

    let expected_k_per_s =
        analytic_equatorial_upwelling_m_per_s(TRADE_WIND_STRESS_PA, sst.surface_drag_per_s())
            / MIXED_LAYER_DEPTH_M
            * SUBSURFACE_SENSITIVITY_K_PER_M
            * thermocline_anomaly_m;
    let measured = *tendency
        .sst_anomaly_k()
        .unwrap()
        .get(3, equator_row(ny))
        .unwrap();
    assert!(
        expected_k_per_s > 0.0,
        "a deeper thermocline under upwelling must warm the mixed layer, not cool it"
    );
    // The only inexact step between theory and tendency is the C-grid
    // difference that produced `w`, so the tolerance is that difference's own
    // leading truncation term, derived above — with a tenth of itself in hand
    // for the `O(dy⁴)` remainder, which at this resolution is another two
    // orders smaller again.
    let tolerance_k_per_s = 1.1
        * upwelling_truncation_fraction(basin.spacing().dy_m(), sst.surface_drag_per_s())
        * expected_k_per_s;
    assert!(
        (measured - expected_k_per_s).abs() <= tolerance_k_per_s,
        "entrainment should warm at {expected_k_per_s} K/s to within {tolerance_k_per_s} K/s, \
         but the tendency is {measured} K/s"
    );
}

// ---------------------------------------------------------------------------
// The state extension
// ---------------------------------------------------------------------------

#[test]
fn the_linear_core_state_carries_no_sst_anomaly() {
    // The extension is opt-in: the state the validated Epics 01-07 core is
    // integrated in is the three-variable one it has always been.
    let grid = Grid::new(4, 3).unwrap();
    assert!(OceanState::at_rest(grid).sst_anomaly_k().is_none());
    assert!(OceanState::at_rest_with_sst_anomaly(grid)
        .sst_anomaly_k()
        .is_some());
}

#[test]
fn the_state_vector_operations_carry_the_sst_anomaly() {
    // RK4 combines states through `assign` and `add_scaled`; a fourth
    // prognostic variable that those two did not touch would simply never be
    // integrated.
    let grid = Grid::new(2, 2).unwrap();
    let source = uniform_coupled_state(grid, 0.0, 0.0, 3.0);
    let mut target = OceanState::at_rest_with_sst_anomaly(grid);

    target.assign(&source);
    assert_eq!(target.sst_anomaly_k().unwrap().as_slice(), &[3.0; 4]);

    target.add_scaled(2.0, &source);
    assert_eq!(target.sst_anomaly_k().unwrap().as_slice(), &[9.0; 4]);
}

#[test]
#[should_panic(expected = "SST anomaly")]
fn combining_a_coupled_state_with_an_uncoupled_one_is_a_bug_and_panics() {
    // Shape is not in the type here, so a mismatch has to be caught at run
    // time rather than truncated (CODING_STANDARDS.md § Correctness and
    // failure).
    let grid = Grid::new(2, 2).unwrap();
    let mut coupled = OceanState::at_rest_with_sst_anomaly(grid);
    coupled.assign(&OceanState::at_rest(grid));
}

// ---------------------------------------------------------------------------
// The acceptance criterion: the extension is additive
// ---------------------------------------------------------------------------

/// Step a forced basin `steps` times and return its `h`, `u` and `v`.
fn stepped_core_fields(couple_sst: bool, steps: u64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let params = physical_params();
    let basin = equatorial_basin(24, 21, 2.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let dt_s = 1800.0;
    let (mut solver, mut state) = if couple_sst {
        (
            Solver::coupled_to_sst(
                basin.grid(),
                basin.spacing(),
                params,
                plane,
                dt_s,
                sst_params(),
            )
            .expect("half an hour is well inside both timestep bounds"),
            OceanState::at_rest_with_sst_anomaly(basin.grid()),
        )
    } else {
        (
            Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
                .expect("half an hour is well inside both timestep bounds"),
            OceanState::at_rest(basin.grid()),
        )
    };
    // A thermocline anomaly to give the coupling something to feed on: if the
    // SST equation could reach back into the core at all, this is the run in
    // which it would.
    state.h_mut().as_mut_slice().fill(5.0);
    let alizes = WindStressField::uniform(basin.grid(), TRADE_WIND_STRESS_PA, 0.0);

    for step in 0..steps {
        solver.step(&mut state, step as f64 * dt_s, |_t| &alizes);
    }
    (
        state.h().as_slice().to_vec(),
        state.u().as_slice().to_vec(),
        state.v().as_slice().to_vec(),
    )
}

#[test]
fn enabling_the_coupling_leaves_the_validated_core_bit_for_bit_unchanged() {
    // The ticket's acceptance criterion, in its strongest available form. The
    // Epic 07 validation suite passing with the coupling disabled shows the
    // *default* path is untouched; this shows the extension is additive even
    // when it is switched on, because `T'` appears in no term of the
    // shallow-water equations and RK4 combines the four variables
    // component-wise. Equality is exact — anything else would mean the core's
    // arithmetic had been re-associated.
    let uncoupled = stepped_core_fields(false, 200);
    let coupled = stepped_core_fields(true, 200);
    assert_eq!(
        uncoupled.0, coupled.0,
        "the thermocline anomaly `h` differs"
    );
    assert_eq!(uncoupled.1, coupled.1, "the zonal current `u` differs");
    assert_eq!(uncoupled.2, coupled.2, "the meridional current `v` differs");
}

#[test]
fn the_coupled_run_actually_evolves_the_sst_anomaly() {
    // The companion to the test above: "nothing changed" is only interesting
    // if the coupling was doing something. A basin under the alizés with a
    // deeper-than-average thermocline must warm.
    let params = physical_params();
    let basin = equatorial_basin(24, 21, 2.0e5, 4.0e6);
    let plane = BetaPlane::of_basin(params, basin);
    let dt_s = 1800.0;
    let mut solver = Solver::coupled_to_sst(
        basin.grid(),
        basin.spacing(),
        params,
        plane,
        dt_s,
        sst_params(),
    )
    .expect("half an hour is well inside both timestep bounds");
    let mut state = OceanState::at_rest_with_sst_anomaly(basin.grid());
    state.h_mut().as_mut_slice().fill(5.0);
    let alizes = WindStressField::uniform(basin.grid(), TRADE_WIND_STRESS_PA, 0.0);
    for step in 0..200 {
        solver.step(&mut state, f64::from(step) * dt_s, |_t| &alizes);
    }
    let warmest = state
        .sst_anomaly_k()
        .unwrap()
        .as_slice()
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        warmest > 0.0,
        "upwelling of a deeper-than-average thermocline must warm the mixed layer somewhere, \
         but the largest anomaly reached is {warmest} K"
    );
}

// ---------------------------------------------------------------------------
// The scenario switch
// ---------------------------------------------------------------------------

/// A scenario file with an optional `[sst]` section appended.
fn scenario_toml(sst_section: &str) -> String {
    format!(
        "[basin]\n\
         resolution_deg = 2.0\n\
         \n\
         [physics]\n\
         reduced_gravity_m_per_s2 = 0.05\n\
         mean_thermocline_depth_m = 150.0\n\
         rayleigh_damping_per_s = 1e-7\n\
         \n\
         [run]\n\
         dt_s = 1800.0\n\
         total_steps = 10\n\
         output_every_n_steps = 5\n\
         {sst_section}"
    )
}

#[test]
fn a_scenario_without_an_sst_section_is_the_uncoupled_linear_model() {
    let scenario = Scenario::from_toml(&scenario_toml("")).expect("a valid scenario");
    assert!(scenario.sst_params().is_none());
}

#[test]
fn the_sst_section_switches_the_coupling_on_and_round_trips_through_toml() {
    let source = scenario_toml(
        "\n[sst]\n\
         mixed_layer_depth_m = 50.0\n\
         mean_zonal_sst_gradient_k_per_m = -4e-7\n\
         subsurface_temperature_sensitivity_k_per_m = 0.1\n\
         thermal_damping_per_s = 9.26e-8\n",
    );
    let config = ScenarioConfig::from_toml(&source).expect("a valid scenario");
    let scenario = config.build().expect("a runnable scenario");
    let sst = scenario.sst_params().expect("the section switches it on");
    assert_eq!(sst.mixed_layer_depth_m(), 50.0);
    assert_eq!(sst.subsurface_temperature_sensitivity_k_per_m(), 0.1);
    // Omitted, so it takes the Zebiak-Cane two-day surface drag.
    assert_eq!(sst.surface_drag_per_s(), DEFAULT_SURFACE_DRAG_PER_S);

    let reparsed =
        ScenarioConfig::from_toml(&config.to_toml().expect("TOML can hold these numbers"))
            .expect("what the engine writes, the engine reads");
    assert_eq!(reparsed, config);
}

#[test]
fn an_unphysical_mixed_layer_is_refused_by_name() {
    let source = scenario_toml(
        "\n[sst]\n\
         mixed_layer_depth_m = 0.0\n\
         mean_zonal_sst_gradient_k_per_m = -4e-7\n\
         subsurface_temperature_sensitivity_k_per_m = 0.1\n\
         thermal_damping_per_s = 9.26e-8\n",
    );
    let error = Scenario::from_toml(&source).expect_err("a zero mixed layer is not an ocean");
    assert!(
        matches!(error, ScenarioError::Sst(SstParamsError::NotPositive { parameter, .. })
            if parameter == "mixed_layer_depth_m"),
        "the error should name the offending parameter, but it is {error}"
    );
    assert!(error.to_string().contains("[sst]"));
}

#[test]
fn a_state_at_rest_is_the_shape_its_grid_asks_for() {
    // The SST anomaly is a cell-centered field, beside `h`, because that is
    // where the entrainment term needs both of them (ADR-0003).
    let grid = Grid::new(5, 3).unwrap();
    let state = OceanState::at_rest_with_sst_anomaly(grid);
    let (nx, ny) = grid.field_shape(H_STAGGERING);
    let sst = state.sst_anomaly_k().unwrap();
    assert_eq!((sst.nx(), sst.ny()), (nx, ny));
    assert!(sst.as_slice().iter().all(|&value| value == 0.0));
}
