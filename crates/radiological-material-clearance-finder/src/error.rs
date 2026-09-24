//! The one error type every fallible function in this crate returns.

use std::fmt;

/// Everything that can go wrong, split by what the caller can do about it.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// A string that cannot be read as a single nuclide.
    NuclideName(String),
    /// A nuclide that parses but is absent from the decay data tables, so
    /// "stable" would be a guess rather than a fact.
    UnknownNuclide {
        /// The canonical nuclide name.
        name: String,
        /// The full message.
        message: String,
    },
    /// A quantity that cannot be derived from what the material was given,
    /// such as a volumetric activity for a material with no density.
    InsufficientData(String),
    /// A limit set name that is not registered.
    UnknownLimitSet(String),
    /// A cooling series whose index climbs back above the threshold after
    /// clearing.
    Ingrowth(String),
    /// A file that could not be read.
    Io(String),
    /// Any other invalid input.
    Invalid(String),
}

impl Error {
    pub(crate) fn unknown_nuclide(name: &str, what: &str) -> Self {
        Error::UnknownNuclide {
            name: name.to_string(),
            message: format!(
                "no {what} for '{name}'. It parses as a nuclide but is not in the \
                 vendored tables, so it is neither known-stable nor known-radioactive. \
                 Check the spelling, or supply it through the decay_data argument."
            ),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::UnknownNuclide { message, .. } => f.write_str(message),
            Error::NuclideName(m)
            | Error::InsufficientData(m)
            | Error::UnknownLimitSet(m)
            | Error::Ingrowth(m)
            | Error::Io(m)
            | Error::Invalid(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for Error {}

/// `Result` with this crate's [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
