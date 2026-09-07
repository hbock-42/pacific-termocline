//! The words the visualizer puts on screen that a visitor may not have, and
//! what each one means in words they do.
//!
//! This project's vocabulary is fixed by `CONTEXT.md` and used everywhere:
//! *thermocline*, *anomaly*, *alizés*. That is the right vocabulary for the
//! code, the tests and the issues — and the wrong one for a stranger looking at
//! a picture of the ocean, who is owed either a plain word or an explanation of
//! the technical one.
//!
//! So there are exactly two moves here, and no third. A term is *replaced* by
//! plain words where plain words say the same thing — "alizés" becomes "trade
//! winds", which is the same wind in a language the reader already speaks — or
//! it is *glossed*, kept and immediately explained, where the term is itself
//! worth learning. "Thermocline" is worth learning: it is what the whole
//! picture is of.
//!
//! The glosses live here rather than inline in the panels so that one term has
//! one explanation wherever it appears, and so that
//! `tests/teaching_labels.rs` can hold the panel's prose to it.

/// A term the project's vocabulary needs and a visitor does not have, together
/// with what it means in plain words.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlainTerm {
    /// The term as the project writes it, in lower case.
    term: &'static str,
    /// What it means, in words that use no term from this list.
    plain: &'static str,
}

impl PlainTerm {
    /// The term as the project writes it.
    #[must_use]
    pub const fn term(&self) -> &'static str {
        self.term
    }

    /// What it means, in plain words.
    #[must_use]
    pub const fn plain(&self) -> &'static str {
        self.plain
    }

    /// The term with its meaning attached: `thermocline (the boundary
    /// between …)`.
    ///
    /// One sentence, not a footnote: a gloss a reader has to go and look up
    /// somewhere else on the screen is a gloss they will not read.
    #[must_use]
    pub fn glossed(&self) -> String {
        format!("{} ({})", self.term, self.plain)
    }
}

/// "Thermocline" — the object of the whole simulation, and the one word here
/// worth a visitor learning.
pub const THERMOCLINE: PlainTerm = PlainTerm {
    term: "thermocline",
    plain: "the boundary between the warm water on top and the cold water \
            beneath",
};

/// "Anomaly" — the model's fields are departures from a resting state
/// (`CONTEXT.md`, *Thermocline depth anomaly*), which is not what the word
/// means in ordinary speech.
pub const ANOMALY: PlainTerm = PlainTerm {
    term: "anomaly",
    plain: "how far something sits from its usual value",
};

/// "Trade winds" — the plain-English replacement for *alizés*
/// (`CONTEXT.md`, *Alizés*, which gives the two as interchangeable), glossed
/// in turn because the name says nothing about what the wind does.
pub const TRADE_WINDS: PlainTerm = PlainTerm {
    term: "trade winds",
    plain: "the steady winds that blow from east to west along the equator",
};

/// The plain-language key the teaching panel carries, in the order it shows
/// them.
///
/// Three, because these are the three the panels use. A term that has to be
/// added here is a term that has appeared on screen without an explanation,
/// which is what `tests/teaching_labels.rs` checks the panel's own prose for.
pub const PLAIN_WORDS: [PlainTerm; 3] = [THERMOCLINE, ANOMALY, TRADE_WINDS];
