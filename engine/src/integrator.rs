//! Time integration: classic 4th-order Runge-Kutta, generic over the state.
//!
//! [ADR-0003] chose RK4 over leapfrog for the solver's time step, so the core
//! of a run is `state_{n+1} = rk4.step(state_n, t, dt, rhs)`. Nothing in this
//! module knows about the shallow-water equations, the C-grid or even the
//! number of prognostic variables: it advances anything that behaves like a
//! vector over the reals, which is what [`StateVector`] captures. That keeps
//! the scheme testable against ODEs with analytic solutions, independently of
//! the ocean physics.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

/// A state the integrator can advance: an element of a real vector space.
///
/// The two operations are deliberately in-place. The textbook signature
/// `Fn(&S, f64) -> S` would allocate a fresh state for each of the four stages
/// of every step, which for a basin-sized state is an allocation in the inner
/// time-stepping loop (CODING_STANDARDS.md § Performance). Implementors write
/// into `self` instead, and [`Rk4`] owns the scratch states.
///
/// Both shipped implementations are shape-safe by construction: `f64` is a
/// single component, and `[f64; N]` fixes its length in the type. An
/// implementation whose shape is *not* in the type — a basin-sized field,
/// arriving with the physics — must panic rather than truncate when asked to
/// combine mismatched shapes: that is a bug in the caller rather than bad user
/// input, and truncating would be exactly the silent clamping
/// CODING_STANDARDS.md § Correctness and failure forbids.
pub trait StateVector: Clone {
    /// Overwrite `self` with a copy of `source`.
    fn assign(&mut self, source: &Self);

    /// `self += factor * other`, component-wise.
    fn add_scaled(&mut self, factor: f64, other: &Self);
}

impl StateVector for f64 {
    fn assign(&mut self, source: &Self) {
        *self = *source;
    }

    fn add_scaled(&mut self, factor: f64, other: &Self) {
        *self += factor * other;
    }
}

impl<const N: usize> StateVector for [f64; N] {
    fn assign(&mut self, source: &Self) {
        self.copy_from_slice(source);
    }

    fn add_scaled(&mut self, factor: f64, other: &Self) {
        // Both operands are `[f64; N]` for the same `N`, so the type already
        // rules out the shape mismatch the trait promises to panic on — there
        // is nothing here for `zip` to silently truncate.
        for (value, term) in self.iter_mut().zip(other.iter()) {
            *value += factor * term;
        }
    }
}

/// Node of the two midpoint stages in the RK4 tableau, as a fraction of `dt`.
const MIDPOINT_NODE: f64 = 0.5;
/// Weight of the endpoint stages `k1` and `k4` in the RK4 tableau.
const ENDPOINT_WEIGHT: f64 = 1.0 / 6.0;
/// Weight of the midpoint stages `k2` and `k3` in the RK4 tableau.
const MIDPOINT_WEIGHT: f64 = 1.0 / 3.0;

/// The classic 4th-order Runge-Kutta integrator, holding its own scratch
/// space.
///
/// One `Rk4` is built per run from a prototype state and then reused for every
/// step, so a whole simulation allocates the four stage buffers and the stage
/// state exactly once. It carries no time and no state of its own: the caller
/// owns `state`, `t` and `dt`, which leaves timestep selection (T-01.3) and
/// output cadence outside the integrator.
#[derive(Debug, Clone)]
pub struct Rk4<S> {
    k1: S,
    k2: S,
    k3: S,
    k4: S,
    stage: S,
}

impl<S: StateVector> Rk4<S> {
    /// An integrator sized for states shaped like `prototype`.
    ///
    /// The prototype's *values* are irrelevant — only its shape is kept, as
    /// the shape of the scratch buffers.
    #[must_use]
    pub fn new(prototype: &S) -> Self {
        Self {
            k1: prototype.clone(),
            k2: prototype.clone(),
            k3: prototype.clone(),
            k4: prototype.clone(),
            stage: prototype.clone(),
        }
    }

    /// Advance `state` from time `t` to `t + dt` in place.
    ///
    /// `rhs(state, t, out)` writes the time derivative of `state` at time `t`
    /// into `out`; it is called exactly four times per step, at `t`,
    /// `t + dt/2`, `t + dt/2` and `t + dt`. `dt` is in seconds, as is `t`, and
    /// `rhs` therefore returns state units per second.
    pub fn step<F>(&mut self, state: &mut S, t: f64, dt: f64, rhs: &mut F)
    where
        F: FnMut(&S, f64, &mut S),
    {
        let half_dt = MIDPOINT_NODE * dt;

        // The four stages are written out rather than looped: each one has its
        // own node and its own predecessor, so a loop would need the tableau as
        // data, and the tableau is easier to audit against ADR-0003 in this
        // form.

        rhs(state, t, &mut self.k1);

        self.stage.assign(state);
        self.stage.add_scaled(half_dt, &self.k1);
        rhs(&self.stage, t + half_dt, &mut self.k2);

        self.stage.assign(state);
        self.stage.add_scaled(half_dt, &self.k2);
        rhs(&self.stage, t + half_dt, &mut self.k3);

        self.stage.assign(state);
        self.stage.add_scaled(dt, &self.k3);
        rhs(&self.stage, t + dt, &mut self.k4);

        state.add_scaled(ENDPOINT_WEIGHT * dt, &self.k1);
        state.add_scaled(MIDPOINT_WEIGHT * dt, &self.k2);
        state.add_scaled(MIDPOINT_WEIGHT * dt, &self.k3);
        state.add_scaled(ENDPOINT_WEIGHT * dt, &self.k4);
    }
}

#[cfg(test)]
mod tests {
    use super::{Rk4, StateVector};

    #[test]
    fn add_scaled_accumulates_component_wise() {
        // Written out by hand: [1, 2] + 3 * [10, 20] = [31, 62].
        let mut state = [1.0_f64, 2.0];
        state.add_scaled(3.0, &[10.0, 20.0]);
        assert_eq!(state, [31.0, 62.0]);
    }

    #[test]
    fn the_right_hand_side_is_evaluated_four_times_per_step() {
        // Four stages is the definition of classic RK4; a step that called the
        // right-hand side a different number of times would be a different
        // scheme.
        let mut evaluations = 0_u32;
        let mut integrator = Rk4::new(&0.0_f64);
        let mut state = 0.0_f64;
        integrator.step(&mut state, 0.0, 0.1, &mut |_y: &f64, _t, out: &mut f64| {
            evaluations += 1;
            *out = 1.0;
        });
        assert_eq!(evaluations, 4);
    }

    #[test]
    fn the_stage_times_are_the_butcher_tableau_nodes() {
        // Nodes c = [0, 1/2, 1/2, 1] of the classic RK4 tableau, offset from
        // the step's start time.
        let mut times = Vec::new();
        let mut integrator = Rk4::new(&0.0_f64);
        let mut state = 0.0_f64;
        integrator.step(&mut state, 10.0, 2.0, &mut |_y: &f64, t, out: &mut f64| {
            times.push(t);
            *out = 0.0;
        });
        assert_eq!(times, vec![10.0, 11.0, 11.0, 12.0]);
    }
}
