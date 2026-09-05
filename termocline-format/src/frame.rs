//! One saved timestep of a run.

use serde::{Deserialize, Serialize};

use crate::{FormatError, GridSpec, RunHeader, Variable};

/// The ocean state and the forcing at one instant of model time.
///
/// Each field is a flat, row-major buffer (`x` varies fastest) holding one
/// value per point of that variable's staggered position on the grid; the
/// lengths come from [`GridSpec::field_len`], so a reader sizes the fields
/// from the header rather than from the frame.
///
/// # The absent SST anomaly
///
/// The five fields of the linear core are always there. `T'` is an
/// [`Option`] because it is not: a scenario without an `[sst]` section never
/// integrates one, and `CONTEXT.md` is explicit that the anomaly is not part
/// of the linear ocean core. A buffer of zeros would round-trip perfectly
/// well and would claim, in a unit the reader is told to believe, that the
/// whole basin sat at exactly its climatological temperature — a physical
/// statement about a run that never made it. `None` says the run has no SST
/// to report, which is the true one, and it costs one byte a frame rather
/// than one `f64` a cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    /// Written as `t`, the name the format gives it; `_s` states the unit on
    /// the Rust side.
    #[serde(rename = "t")]
    t_s: f64,
    h: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    tau_x: Vec<f64>,
    tau_y: Vec<f64>,
    /// `T'` in kelvin, or `None` for a run that did not couple SST. Never
    /// skipped when it is `None`: the field is part of the frame's shape, and
    /// a `bincode` reader counts on every frame having the same one.
    sst: Option<Vec<f64>>,
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
            sst: None,
        };
        frame.validate(grid)?;
        Ok(frame)
    }

    /// The same frame carrying the mixed-layer SST anomaly `T'`, in kelvin, at
    /// cell centers.
    ///
    /// The only way a frame comes to hold `T'`, and the only difference
    /// between a coupled run's frames and an uncoupled run's: a writer builds
    /// the core frame either way and adds the sixth field when — and only
    /// when — the run has one to add.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `sst_anomaly_k` does not carry exactly
    /// one value per cell of `grid`.
    pub fn with_sst_anomaly(
        mut self,
        grid: &GridSpec,
        sst_anomaly_k: Vec<f64>,
    ) -> Result<Self, FormatError> {
        self.sst = Some(sst_anomaly_k);
        self.validate(grid)?;
        Ok(self)
    }

    /// Check every field against the grid the header describes.
    ///
    /// [`Frame::new`] calls this, but deserialization does not: a reader that
    /// has just decoded a frame from disk calls it to confirm the frame and
    /// the header agree.
    ///
    /// A variable this frame does not carry is not checked: the absence of
    /// `T'` is a fact about the run, not a field of the wrong length.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] naming the first variable whose field does
    /// not carry one value per point of its staggered position on `grid`.
    pub fn validate(&self, grid: &GridSpec) -> Result<(), FormatError> {
        for variable in Variable::ALL {
            let Some(field) = self.field(variable) else {
                continue;
            };
            let actual = field.len();
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

    /// Check this frame against the whole of `header`: the grid its fields
    /// must cover, and the variables its run declared.
    ///
    /// What a reader calls, because a frame is bytes until the header says
    /// what shape they should be *and* what they should mean. A frame carrying
    /// a variable its header does not list would be a field nothing announced;
    /// one missing a variable its header lists would be a run that promised
    /// data it does not have. Either way the run is refused rather than read
    /// under labels that do not fit it.
    ///
    /// # Errors
    /// The errors of [`Frame::validate`], and
    /// [`FormatError::UndeclaredVariable`] naming the first variable the frame
    /// and the header disagree about.
    pub fn validate_against(&self, header: &RunHeader) -> Result<(), FormatError> {
        self.validate(&header.grid)?;
        for variable in Variable::ALL {
            let declared = header.carries(variable);
            if self.field(variable).is_some() != declared {
                return Err(FormatError::UndeclaredVariable { variable, declared });
            }
        }
        Ok(())
    }

    /// The field of `variable`, row-major, or `None` if this run does not
    /// carry that variable. The per-variable accessors below are the same
    /// buffers under the names the model uses; this one is for code that walks
    /// a variable list rather than naming a field.
    ///
    /// Only [`Variable::SstAnomaly`] is ever `None`; the five of
    /// [`Variable::LINEAR_CORE`] are in every frame. The header's `variables`
    /// is what says which to expect, and it is written from the same fact.
    ///
    /// Units and staggering of each field come from [`Variable::unit`] and
    /// [`Variable::staggering`], which the header also writes out.
    #[must_use]
    pub fn field(&self, variable: Variable) -> Option<&[f64]> {
        match variable {
            Variable::ThermoclineDepthAnomaly => Some(&self.h),
            Variable::ZonalCurrentAnomaly => Some(&self.u),
            Variable::MeridionalCurrentAnomaly => Some(&self.v),
            Variable::ZonalWindStress => Some(&self.tau_x),
            Variable::MeridionalWindStress => Some(&self.tau_y),
            Variable::SstAnomaly => self.sst.as_deref(),
        }
    }

    /// Model time of this frame, in seconds since the start of the run.
    #[must_use]
    pub const fn t_s(&self) -> f64 {
        self.t_s
    }

    /// Thermocline depth anomaly `h`.
    #[must_use]
    pub fn h(&self) -> &[f64] {
        &self.h
    }

    /// Zonal current anomaly `u`.
    #[must_use]
    pub fn u(&self) -> &[f64] {
        &self.u
    }

    /// Meridional current anomaly `v`.
    #[must_use]
    pub fn v(&self) -> &[f64] {
        &self.v
    }

    /// Zonal wind stress `τx`.
    #[must_use]
    pub fn tau_x(&self) -> &[f64] {
        &self.tau_x
    }

    /// Meridional wind stress `τy`.
    #[must_use]
    pub fn tau_y(&self) -> &[f64] {
        &self.tau_y
    }

    /// Mixed-layer SST anomaly `T'`, in kelvin, or `None` if this run did not
    /// couple SST.
    ///
    /// `None` means "this run has no SST to report", and a caller has to say
    /// what it does about that rather than reading a zero it would have no way
    /// to tell from a real one.
    #[must_use]
    pub fn sst_anomaly_k(&self) -> Option<&[f64]> {
        self.sst.as_deref()
    }

    /// Whether this frame carries the Epic 12 SST anomaly.
    #[must_use]
    pub const fn carries_sst_anomaly(&self) -> bool {
        self.sst.is_some()
    }
}

/// A frame as [`crate::FORMAT_VERSION`] 1 wrote it: the five fields of the
/// linear core and nothing after them.
///
/// Version 1 predates the SST coupling, so it has no place in its bytes where
/// `T'` could have been — not even a `None` tag. Decoding a version 1 frame
/// with the current [`Frame`] would read the *next* frame's first byte as the
/// tag and desynchronize the rest of the file, which is why the old layout is
/// kept here as a type rather than approximated by the new one.
///
/// Deserialize only: nothing writes version 1 any more.
#[derive(Deserialize)]
pub(crate) struct FrameV1 {
    #[serde(rename = "t")]
    t_s: f64,
    h: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    tau_x: Vec<f64>,
    tau_y: Vec<f64>,
}

impl From<FrameV1> for Frame {
    /// A version 1 frame is a current frame without an SST anomaly — `None`
    /// rather than zeros, because a run written before the coupling existed
    /// has no `T'`, and saying it was everywhere zero would be inventing one.
    fn from(frame: FrameV1) -> Self {
        let FrameV1 {
            t_s,
            h,
            u,
            v,
            tau_x,
            tau_y,
        } = frame;
        Self {
            t_s,
            h,
            u,
            v,
            tau_x,
            tau_y,
            sst: None,
        }
    }
}
