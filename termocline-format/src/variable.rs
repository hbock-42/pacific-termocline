//! The variables a run writes, and what a reader needs to know about each.

use serde::{Deserialize, Serialize};
use termocline_grid::{Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// A field a run writes once per frame.
///
/// The set is closed: it is the state and forcing of the 1.5-layer model
/// (`docs/planning/01-scientific-model.md`) plus the SST anomaly of the Epic
/// 12 coupling extension, and adding a member is a format change, hence a
/// [`crate::FORMAT_VERSION`] bump.
///
/// Not every run carries every member. [`Variable::LINEAR_CORE`] is what every
/// run has; the SST anomaly is there only when the scenario asked for the
/// coupling, and a run's header says which of the two lists it wrote.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Variable {
    /// Thermocline depth anomaly `h` — a departure from the mean depth `H`,
    /// never a total depth.
    ThermoclineDepthAnomaly,
    /// Zonal (eastward) current anomaly `u`.
    ZonalCurrentAnomaly,
    /// Meridional (northward) current anomaly `v`.
    MeridionalCurrentAnomaly,
    /// Zonal component of the wind stress forcing, `τx`. Easterly trade-wind
    /// stress is negative.
    ZonalWindStress,
    /// Meridional component of the wind stress forcing, `τy`.
    MeridionalWindStress,
    /// Mixed-layer sea-surface-temperature anomaly `T'` of the Epic 12
    /// coupling — a departure from the climatological SST, never an absolute
    /// temperature. Present only in a run whose scenario asked for the
    /// coupling; `CONTEXT.md` is explicit that it is not part of the linear
    /// ocean core.
    SstAnomaly,
}

impl Variable {
    /// Every variable the format knows, in the order a header lists them.
    ///
    /// The *format's* list, not any one run's: a run carries
    /// [`Variable::LINEAR_CORE`] and, when its scenario asked for the
    /// coupling, [`Variable::SstAnomaly`] as well. What a given run carries is
    /// its header's `variables`.
    pub const ALL: [Self; 6] = [
        Self::ThermoclineDepthAnomaly,
        Self::ZonalCurrentAnomaly,
        Self::MeridionalCurrentAnomaly,
        Self::ZonalWindStress,
        Self::MeridionalWindStress,
        Self::SstAnomaly,
    ];

    /// The variables of the linear ocean core, which every run carries.
    ///
    /// The three prognostic fields of the 1.5-layer model and the two
    /// components of the wind stress that forced them — the whole of a run
    /// that does not couple SST, and the first five entries of one that does.
    pub const LINEAR_CORE: [Self; 5] = [
        Self::ThermoclineDepthAnomaly,
        Self::ZonalCurrentAnomaly,
        Self::MeridionalCurrentAnomaly,
        Self::ZonalWindStress,
        Self::MeridionalWindStress,
    ];

    /// The symbol this variable is written with in `CONTEXT.md`, unaccented so
    /// it doubles as the field name a non-Rust reader looks for.
    ///
    /// The SST anomaly is the one variable whose symbol is not a
    /// transliteration of the one `CONTEXT.md` uses: `T'` unaccented is
    /// `t_prime`, which would sit in a frame beside the `t` that is the
    /// frame's model time and read as a variant of it. `sst` is what
    /// `CONTEXT.md` names the quantity in words, and it cannot be mistaken for
    /// the clock.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::ThermoclineDepthAnomaly => "h",
            Self::ZonalCurrentAnomaly => "u",
            Self::MeridionalCurrentAnomaly => "v",
            Self::ZonalWindStress => "tau_x",
            Self::MeridionalWindStress => "tau_y",
            Self::SstAnomaly => "sst",
        }
    }

    /// The SI unit of the values, so a reader never has to assume one.
    #[must_use]
    pub const fn unit(self) -> &'static str {
        match self {
            Self::ThermoclineDepthAnomaly => "m",
            Self::ZonalCurrentAnomaly | Self::MeridionalCurrentAnomaly => "m s^-1",
            Self::ZonalWindStress | Self::MeridionalWindStress => "N m^-2",
            Self::SstAnomaly => "K",
        }
    }

    /// Where this variable sits on the Arakawa C-grid of [ADR-0003].
    ///
    /// Each wind-stress component sits where the current it forces sits: `τx`
    /// enters the `u` equation, `τy` the `v` equation. The SST anomaly sits at
    /// cell centers with `h`, which is the field the entrainment term
    /// multiplies it against point for point.
    ///
    /// [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md
    #[must_use]
    pub const fn staggering(self) -> Staggering {
        match self {
            Self::ThermoclineDepthAnomaly | Self::SstAnomaly => H_STAGGERING,
            Self::ZonalCurrentAnomaly | Self::ZonalWindStress => U_STAGGERING,
            Self::MeridionalCurrentAnomaly | Self::MeridionalWindStress => V_STAGGERING,
        }
    }
}

/// One entry of the header's variable list: what a frame's field means, and in
/// what unit, spelled out for a reader that does not share these Rust types.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariableSpec {
    /// The variable this entry describes.
    pub variable: Variable,
    /// Its symbol, as written in `CONTEXT.md`.
    pub symbol: String,
    /// Its SI unit.
    pub unit: String,
}

impl VariableSpec {
    /// The header entry describing `variable`.
    #[must_use]
    pub fn of(variable: Variable) -> Self {
        Self {
            variable,
            symbol: variable.symbol().to_owned(),
            unit: variable.unit().to_owned(),
        }
    }
}
