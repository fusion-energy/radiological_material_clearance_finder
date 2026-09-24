//! Nuclide name parsing, validation and normalisation.
//!
//! The canonical form used throughout this crate is the one OpenMC and GND use,
//! `Co60` and `Ag108_m1`, because the most common source of an inventory is an
//! OpenMC material and matching it avoids a translation layer at the boundary.
//!
//! Regulatory tables spell nuclides differently again (`Co-60`, `Ag-108m`,
//! `U-238sec`, `Sr-90+`), so [`parse_regulatory`] handles those and separates
//! the secular-equilibrium marker from the nuclide identity.

use crate::error::{Error, Result};

/// Element symbols indexed by atomic number, with index 0 unused.
pub const SYMBOLS: [&str; 119] = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// The secular equilibrium marker on a regulatory table row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Marker {
    /// `+`: the tabulated daughters are covered by the parent's value.
    Plus,
    /// `sec`: the whole decay chain is covered.
    Sec,
}

impl Marker {
    /// The marker as the tables write it, `"+"` or `"sec"`.
    pub fn as_str(self) -> &'static str {
        match self {
            Marker::Plus => "+",
            Marker::Sec => "sec",
        }
    }
}

/// Case-insensitive atomic number of an element symbol.
fn z_of(symbol: &str) -> Option<u32> {
    SYMBOLS
        .iter()
        .skip(1)
        .position(|s| s.eq_ignore_ascii_case(symbol))
        .map(|i| i as u32 + 1)
}

fn capitalise(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(first) => first
            .to_uppercase()
            .chain(chars.flat_map(char::to_lowercase))
            .collect(),
        None => String::new(),
    }
}

/// Match `symbol, optional separator, mass number, optional metastable marker`.
///
/// The marker covers "m", "m1", "m2", "_m1", "_m2", the "n" that the German
/// regulation uses for the second metastable state ("Hf-178n" is Hf178_m2),
/// and the space-separated form EUR-Lex uses ("Zn-69 m").
fn match_pattern(name: &str) -> Option<(&str, u32, u8)> {
    let text = name.trim_matches(char::is_whitespace);
    let letters = text.bytes().take_while(u8::is_ascii_alphabetic).count();
    if !(1..=2).contains(&letters) {
        return None;
    }
    let symbol = &text[..letters];
    let mut rest = &text[letters..];

    // One optional whitespace or hyphen separator.
    if let Some(c) = rest.chars().next() {
        if c == '-' || c.is_whitespace() {
            rest = &rest[c.len_utf8()..];
        }
    }

    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    if !(1..=3).contains(&digits) {
        return None;
    }
    let mass: u32 = rest[..digits].parse().ok()?;
    rest = &rest[digits..];

    if rest.is_empty() {
        return Some((symbol, mass, 0));
    }

    // One optional whitespace or underscore before the marker.
    if let Some(c) = rest.chars().next() {
        if c == '_' || c.is_whitespace() {
            rest = &rest[c.len_utf8()..];
        }
    }
    let mut chars = rest.chars();
    let meta = chars.next()?.to_ascii_lowercase();
    if meta != 'm' && meta != 'n' {
        return None;
    }
    let tail = chars.as_str();
    let level = match tail.chars().next() {
        None => {
            // "n" is the second metastable state, "m" the first.
            if meta == 'n' {
                2
            } else {
                1
            }
        }
        Some(d) if d.is_ascii_digit() && tail.len() == 1 => d as u8 - b'0',
        Some(_) => return None,
    };
    Some((symbol, mass, level))
}

/// Split a nuclide name into element symbol, mass number and metastable level.
///
/// Accepts any spelling such as `Co60`, `Co-60`, `Ag-108m` or `Ag108_m1`. The
/// level is 0 for a ground state and 1 or 2 for a metastable state. The symbol
/// is returned in its canonical capitalisation regardless of the input's.
///
/// A bare element such as `Fe` is an error, because a clearance limit applies
/// to a nuclide and silently guessing a mass number would be wrong.
pub fn parse(name: &str) -> Result<(&'static str, u32, u8)> {
    let Some((symbol_in, mass, level)) = match_pattern(name) else {
        let trimmed = name.trim_matches(char::is_whitespace);
        if (1..=2).contains(&trimmed.len()) && trimmed.bytes().all(|b| b.is_ascii_alphabetic()) {
            return Err(Error::NuclideName(format!(
                "'{name}' is an element, not a nuclide. Clearance limits are per \
                 nuclide, so give a mass number, for example {}56.",
                capitalise(trimmed)
            )));
        }
        return Err(Error::NuclideName(format!(
            "cannot read '{name}' as a nuclide name"
        )));
    };

    let z = z_of(symbol_in).ok_or_else(|| {
        Error::NuclideName(format!("unknown element symbol '{symbol_in}' in '{name}'"))
    })?;

    // The parser validates syntax, not the drip lines, but these two bounds
    // catch the common typo class without pretending to know nuclear
    // stability. The heaviest nuclide in AME2020 is Og295.
    if mass > 300 {
        return Err(Error::NuclideName(format!(
            "'{name}' has mass number {mass}, above the heaviest known nuclide"
        )));
    }
    if mass < z {
        return Err(Error::NuclideName(format!(
            "'{name}' has mass number {mass} below its proton number {z}, which is \
             not a real nuclide"
        )));
    }
    Ok((SYMBOLS[z as usize], mass, level))
}

fn canonical(symbol: &str, mass: u32, level: u8) -> String {
    if level == 0 {
        format!("{symbol}{mass}")
    } else {
        format!("{symbol}{mass}_m{level}")
    }
}

/// Return a nuclide name in canonical form.
///
/// `Co-60` and `co 60` both become `Co60`; `Ag-108m` becomes `Ag108_m1` and
/// `Hf-178n` becomes `Hf178_m2`. A trailing `+` is stripped, since it marks a
/// secular equilibrium value rather than a different nuclide, but `sec` is not
/// accepted here because a `sec` value is a separate table entry. Use
/// [`parse_regulatory`] when reading a regulatory table.
pub fn normalise(name: &str) -> Result<String> {
    let trimmed = name.trim_end_matches(char::is_whitespace);
    let stripped = match trimmed.strip_suffix('+') {
        Some(rest) => rest.trim_end_matches(char::is_whitespace),
        None => name,
    };
    let (symbol, mass, level) = parse(stripped)?;
    Ok(canonical(symbol, mass, level))
}

/// Read a nuclide name as spelled in a regulatory table.
///
/// Regulatory tables mark parent nuclides whose limit already accounts for
/// daughters in secular equilibrium, with `+` for the daughters tabulated
/// alongside and `sec` for the whole decay chain. The marker is not part of the
/// nuclide's identity but it does select a different limit value, so it is
/// returned separately rather than discarded.
pub fn parse_regulatory(label: &str) -> Result<(String, Option<Marker>)> {
    let trimmed = label.trim_end_matches(char::is_whitespace);
    let (text, marker) = if let Some(rest) = trimmed.strip_suffix('+') {
        (rest, Some(Marker::Plus))
    } else if trimmed.len() >= 3
        && trimmed.is_char_boundary(trimmed.len() - 3)
        && trimmed[trimmed.len() - 3..].eq_ignore_ascii_case("sec")
    {
        (&trimmed[..trimmed.len() - 3], Some(Marker::Sec))
    } else {
        (label, None)
    };
    let (symbol, mass, level) = parse(text)?;
    Ok((canonical(symbol, mass, level), marker))
}

/// The element symbol of a nuclide, for example `Co` for `Co60`.
pub fn element(name: &str) -> Result<&'static str> {
    Ok(parse(name)?.0)
}

/// The mass number of a nuclide, for example 60 for `Co60`.
pub fn mass_number(name: &str) -> Result<u32> {
    Ok(parse(name)?.1)
}

/// The metastable level of a nuclide, 0 for a ground state.
pub fn metastable_state(name: &str) -> Result<u8> {
    Ok(parse(name)?.2)
}

/// The proton number of a nuclide, for example 27 for `Co60`.
pub fn atomic_number(name: &str) -> Result<u32> {
    let symbol = parse(name)?.0;
    Ok(z_of(symbol).expect("parse returns a known symbol"))
}

/// Whether a string can be read as a single nuclide.
pub fn is_valid(name: &str) -> bool {
    parse(name).is_ok()
}

/// The ground state name of a canonical nuclide name.
pub(crate) fn ground_state(canonical: &str) -> &str {
    match canonical.find("_m") {
        Some(i) => &canonical[..i],
        None => canonical,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spellings_normalise() {
        for (given, expected) in [
            ("Co60", "Co60"),
            ("Co-60", "Co60"),
            ("co 60", "Co60"),
            ("CO-60", "Co60"),
            ("Ag-108m", "Ag108_m1"),
            ("Ag108_m1", "Ag108_m1"),
            ("Ag108m2", "Ag108_m2"),
            ("Hf-178n", "Hf178_m2"),
            ("Zn-69 m", "Zn69_m1"),
            ("Sr-90+", "Sr90"),
            ("Og295", "Og295"),
            ("Rn222", "Rn222"),
            ("Sn113", "Sn113"),
        ] {
            assert_eq!(normalise(given).unwrap(), expected, "{given}");
        }
    }

    #[test]
    fn regulatory_markers_are_separated() {
        assert_eq!(
            parse_regulatory("Sr-90+").unwrap(),
            ("Sr90".into(), Some(Marker::Plus))
        );
        assert_eq!(
            parse_regulatory("U-238sec").unwrap(),
            ("U238".into(), Some(Marker::Sec))
        );
        assert_eq!(
            parse_regulatory("Th-232 sec").unwrap(),
            ("Th232".into(), Some(Marker::Sec))
        );
        assert_eq!(parse_regulatory("Co-60").unwrap(), ("Co60".into(), None));
    }

    #[test]
    fn bad_names_are_rejected() {
        for bad in ["", "Xx60", "H-500", "not a nuclide", "U50", "Fe", "Co6000"] {
            assert!(parse(bad).is_err(), "{bad}");
        }
        assert!(parse("Fe")
            .unwrap_err()
            .to_string()
            .contains("element, not a nuclide"));
    }
}
