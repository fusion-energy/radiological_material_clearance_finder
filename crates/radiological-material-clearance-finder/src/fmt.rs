//! Number formatting for messages and reports.

/// Format like C's and Python's `%g`: `precision` significant figures, fixed
/// or scientific by magnitude, with trailing zeros removed.
pub(crate) fn g(value: f64, precision: usize) -> String {
    let p = precision.max(1);
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    // Rounding to p significant figures first gives the exponent the choice
    // between fixed and scientific is made on, as %g does.
    let scientific = format!("{:.*e}", p - 1, value);
    let (mantissa, exponent) = scientific.split_once('e').expect("{:e} always has an e");
    let exponent: i32 = exponent.parse().expect("{:e} exponent is an integer");
    if exponent < -4 || exponent >= p as i32 {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", strip_zeros(mantissa), sign, exponent.abs())
    } else {
        let decimals = (p as i32 - 1 - exponent).max(0) as usize;
        strip_zeros(&format!("{:.*}", decimals, value)).to_string()
    }
}

fn strip_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

#[cfg(test)]
mod tests {
    use super::g;

    #[test]
    fn matches_percent_g() {
        assert_eq!(g(4000.0, 6), "4000");
        assert_eq!(g(11.244691988476553, 4), "11.24");
        assert_eq!(g(1.0, 6), "1");
        assert_eq!(g(1e-5, 4), "1e-05");
        assert_eq!(g(123456789.0, 4), "1.235e+08");
        assert_eq!(g(0.0001234, 3), "0.000123");
        assert_eq!(g(9.9996, 4), "10");
        assert_eq!(g(-2.5, 3), "-2.5");
    }
}
