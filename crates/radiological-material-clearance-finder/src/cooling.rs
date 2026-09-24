//! Finding when a material becomes clearable.
//!
//! This crate computes no decay of its own. The caller supplies a series of
//! materials at increasing cooling times, which is what a depletion or
//! transmutation calculation already produces, and this finds where the
//! clearance index crosses the threshold.

use crate::error::{Error, Result};
use crate::fmt::g;
use crate::index::{clearance_index, ClearanceOptions};
use crate::limits::LimitSet;
use crate::material::Material;

/// Evaluate the clearance index at each cooling time.
///
/// `series` pairs a cooling time in seconds with the material at that time,
/// in any order. Returns `(time, index)` ordered by time. Two entries at the
/// same time are an error, since they could disagree.
pub fn index_series<'a, I>(
    series: I,
    limit_set: &LimitSet,
    options: ClearanceOptions,
) -> Result<Vec<(f64, f64)>>
where
    I: IntoIterator<Item = (f64, &'a Material)>,
{
    let mut entries: Vec<(f64, &Material)> = series.into_iter().collect();
    if let Some((t, _)) = entries.iter().find(|(t, _)| !t.is_finite()) {
        return Err(Error::Invalid(format!("cooling time {t:?} is not finite")));
    }
    entries.sort_by(|a, b| a.0.total_cmp(&b.0));
    if let Some(pair) = entries.windows(2).find(|w| w[0].0 == w[1].0) {
        return Err(Error::Invalid(format!(
            "cooling time {:?} appears more than once in the series",
            pair[0].0
        )));
    }
    entries
        .into_iter()
        .map(|(t, m)| Ok((t, clearance_index(m, limit_set, options)?.index)))
        .collect()
}

/// Find the cooling time at which a material first meets a set of limits.
///
/// The index falls roughly exponentially with cooling time, so this
/// interpolates logarithmically in the index between the two samples that
/// bracket the threshold. Interpolating linearly instead would place the
/// crossing systematically late, by a factor that grows with the spacing of
/// the samples.
///
/// Returns the cooling time in seconds at which the index first falls below the
/// threshold. If the earliest sample already meets it, that sample's time is
/// returned rather than zero, since the series says nothing about anything
/// earlier. `None` if the series never meets it.
///
/// An index that climbs back above the threshold at a later time means a
/// daughter is growing in faster than its parent decays, so the material does
/// not stay clear, and reporting only the first crossing would be misleading.
/// That is an [`Error::Ingrowth`] unless `allow_ingrowth` is set.
pub fn time_to_clear<'a, I>(
    series: I,
    limit_set: &LimitSet,
    options: ClearanceOptions,
    allow_ingrowth: bool,
) -> Result<Option<f64>>
where
    I: IntoIterator<Item = (f64, &'a Material)>,
{
    let threshold = limit_set.threshold();
    let indexes = index_series(series, limit_set, options)?;
    if indexes.len() < 2 {
        return Err(Error::Invalid(
            "a cooling series needs at least two times to interpolate between".into(),
        ));
    }

    let crossing = if indexes[0].1 < threshold {
        Some(indexes[0].0)
    } else {
        indexes.windows(2).find(|w| w[1].1 < threshold).map(|w| {
            let ((time_a, index_a), (time_b, index_b)) = (w[0], w[1]);
            interpolate(time_a, index_a, time_b, index_b, threshold)
        })
    };
    let Some(crossing) = crossing else {
        return Ok(None);
    };

    if !allow_ingrowth {
        if let Some(&(time, value)) = indexes
            .iter()
            .find(|&&(t, v)| t > crossing && v >= threshold)
        {
            return Err(Error::Ingrowth(format!(
                "the index first falls below {} at {} s but climbs back to {} at {} s, so \
                 the material does not stay clear. Pass allow_ingrowth to take the first \
                 crossing anyway.",
                g(threshold, 6),
                g(crossing, 4),
                g(value, 4),
                g(time, 4)
            )));
        }
    }
    Ok(Some(crossing))
}

/// Interpolate the crossing time, logarithmically in the index.
fn interpolate(time_a: f64, index_a: f64, time_b: f64, index_b: f64, threshold: f64) -> f64 {
    if index_a == index_b {
        return time_b;
    }
    let fraction = if index_a <= 0.0 || index_b <= 0.0 {
        // Linear when a log is not defined.
        (index_a - threshold) / (index_a - index_b)
    } else {
        (index_a.ln() - threshold.ln()) / (index_a.ln() - index_b.ln())
    };
    time_a + fraction.clamp(0.0, 1.0) * (time_b - time_a)
}
