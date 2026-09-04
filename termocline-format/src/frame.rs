//! One saved timestep of a run.

use serde::{Deserialize, Serialize};

use crate::{FormatError, GridSpec, Variable};

/// The ocean state and the forcing at one instant of model time.
///
/// Each field is a flat, row-major buffer (`x` varies fastest) holding one
/// value per point of that variable's staggered position on the grid; the
/// lengths come from [`GridSpec::field_len`], so a reader sizes the fields
/// from the header rather than from the frame.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    t_s: f64,
    h: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    tau_x: Vec<f64>,
    tau_y: Vec<f64>,
}

impl Frame {
    /// A frame at model time `t_s` on `grid`.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if any field does not carry exactly one
    /// value per point of its staggered position on `grid`.
    pub fn new(
        t_s: f64,
        grid: &GridSpec,
        h: Vec<f64>,
        u: Vec<f64>,
        v: Vec<f64>,
        tau_x: Vec<f64>,
        tau_y: Vec<f64>,
    ) -> Result<Self, FormatError> {
        let frame = Self {
            t_s,
            h,
            u,
            v,
            tau_x,
            tau_y,
        };
        frame.validate(grid)?;
        Ok(frame)
    }

    /// Check every field against the grid the header describes.
    ///
    /// [`Frame::new`] calls this, but deserialization does not: a reader that
    /// has just decoded a frame from disk calls it to confirm the frame and
    /// the header agree.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] naming the first variable whose field does
    /// not carry one value per point of its staggered position on `grid`.
    pub fn validate(&self, grid: &GridSpec) -> Result<(), FormatError> {
        for variable in Variable::ALL {
            let actual = self.field(variable).len();
            let expected = grid.field_len(variable);
            if actual != expected {
                return Err(FormatError::FieldShape {
                    variable,
                    expected,
                    actual,
                });
            }
        }
        Ok(())
    }

    /// The field of `variable`, row-major.
    #[must_use]
    pub fn field(&self, variable: Variable) -> &[f64] {
        match variable {
            Variable::ThermoclineDepthAnomaly => &self.h,
            Variable::ZonalCurrentAnomaly => &self.u,
            Variable::MeridionalCurrentAnomaly => &self.v,
            Variable::ZonalWindStress => &self.tau_x,
            Variable::MeridionalWindStress => &self.tau_y,
        }
    }

    /// Model time of this frame, in seconds since the start of the run.
    #[must_use]
    pub const fn t_s(&self) -> f64 {
        self.t_s
    }

    /// Thermocline depth anomaly `h`, in metres, at cell centers.
    #[must_use]
    pub fn h(&self) -> &[f64] {
        &self.h
    }

    /// Zonal current anomaly `u`, in m s^-1, at cell east/west faces.
    #[must_use]
    pub fn u(&self) -> &[f64] {
        &self.u
    }

    /// Meridional current anomaly `v`, in m s^-1, at cell north/south faces.
    #[must_use]
    pub fn v(&self) -> &[f64] {
        &self.v
    }

    /// Zonal wind stress `τx`, in N m^-2, where `u` lives.
    #[must_use]
    pub fn tau_x(&self) -> &[f64] {
        &self.tau_x
    }

    /// Meridional wind stress `τy`, in N m^-2, where `v` lives.
    #[must_use]
    pub fn tau_y(&self) -> &[f64] {
        &self.tau_y
    }
}
